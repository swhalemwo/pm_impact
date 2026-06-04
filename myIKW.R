library(PanelMatch)
library(data.table)
library(fixest)

# Correct extraction: 
# PanelMatch keeps original data in the 'panel.data' attribute/slot
# The PanelMatch object is a list, but we can't just call $panel.data on it
# The PanelData object is what we passed to PanelMatch

prepare_matched_data <- function(pm_object, panel_data) {
  # pm_data is the PanelData object
  # which IS a data.frame in this package
  df_raw <- as.data.frame(panel_data)
  
  matched_data_list <- lapply(names(pm_object$att), function(set_name) {
    set <- pm_object$att[[set_name]]
    weights <- attr(set, "weights")
    
    parts <- strsplit(set_name, "\\.")[[1]]
    unit_id <- as.integer(parts[1])
    time_id <- as.integer(parts[2])
    
    # Set contains the control units
    combined_ids <- c(unit_id, as.integer(set))
    
    # We get data for these units for time_id
    matched_set_df <- df_raw[df_raw$wbcode2 %in% combined_ids & df_raw$year == time_id, ]
    
    matched_set_df$matched_weights <- 1
    # Control units get the weights
    if (!is.null(weights)) {
      matched_set_df$matched_weights[matched_set_df$wbcode2 %in% as.integer(set)] <- weights
    }
    
    matched_set_df$set_id <- set_name
    return(matched_set_df)
  })
  
  return(rbindlist(matched_data_list))
}

data(dem)
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

# Replicate with OLS
lm_est <- lm(y ~ dem, data = matched_df, weights = matched_weights)

cat("PanelEstimate result (t+0):\n")
print(pe_results$estimate["t+0"])
cat("\nReplication (OLS) result:\n")
print(coef(lm_est)["dem"])
