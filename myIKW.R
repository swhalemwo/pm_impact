library(PanelMatch)
library(data.table)
library(fixest)

prepare_matched_data <- function(pm_object, panel_data) {
  df_raw <- as.data.table(panel_data)
  
  matched_data_list <- lapply(names(pm_object$att), function(set_name) {
    set <- pm_object$att[[set_name]]
    weights <- attr(set, "weights")
    
    parts <- strsplit(set_name, "\\.")[[1]]
    unit_id <- as.integer(parts[1])
    time_id <- as.integer(parts[2])
    
    lead_periods <- pm_object$lead
    
    set_dfs <- lapply(lead_periods, function(l) {
      combined_ids <- c(unit_id, as.integer(set))
      curr_time <- time_id + l
      
      matched_set_df <- df_raw[wbcode2 %in% combined_ids & year == curr_time]
      
      if(nrow(matched_set_df) == 0) return(NULL)
      
      w <- rep(1, nrow(matched_set_df))
      ctrl_idx <- which(matched_set_df$wbcode2 %in% as.integer(set))
      matched_ids <- matched_set_df$wbcode2[ctrl_idx]
      w[ctrl_idx] <- weights[match(matched_ids, as.integer(set))]
      
      matched_set_df[, matched_weights := w]
      matched_set_df[, set_id := set_name]
      matched_set_df[, lead := l]
      return(matched_set_df)
    })
    
    return(rbindlist(set_dfs))
  })
  
  return(rbindlist(matched_data_list))
}

data(dem)
dem <- as.data.frame(lapply(dem, function(x) if(is.list(x)) unlist(x) else x))

pm_data <- PanelData(panel.data = dem,
                     unit.id = "wbcode2",
                     time.id = "year",
                     treatment = "dem",
                     outcome = "y")

pm_results <- PanelMatch(panel.data = pm_data, lag = 4, 
                         refinement.method = "ps.match", 
                         match.missing = TRUE, 
                         covs.formula = ~tradewb + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                         size.match = 5, qoi = "att",
                         lead = 0:4, 
                         forbid.treatment.reversal = FALSE)

pe_results <- PanelEstimate(sets = pm_results, panel.data = pm_data, se.method = "bootstrap")

matched_df <- prepare_matched_data(pm_results, pm_data)

lm_est <- lm(y ~ dem, data = matched_df, weights = matched_df$matched_weights)

cat("PanelEstimate (pooled ATT):\n")
print(mean(pe_results$estimate))
cat("\nReplication (OLS):\n")
print(coef(lm_est)["dem"])
