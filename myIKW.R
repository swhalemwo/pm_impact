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

prep_pm_long <- function(pm_object, dt_panel) {
    # pm_object$att[1] %>% names
    
    ## process the sets into one df    
    dt_sets <- keep(pm_object$att, ~len(.x) > 0) %>%
        imap(~attr(.x, "weights") %>% vadt(name = "wbcode2", vlu = "weight") %>%
                 .[, type := "crol"] %>%
                 rbind(data.table(wbcode2 = strsplit(.y, "\\.")[[1]][1], weight = 1, type = "treat")) %>%
                 .[, `:=`(year = strsplit(.y, "\\.")[[1]][2], set = .y)]) %>% rbindlist %>%
        .[, `:=`(wbcode2 = as.integer(wbcode2), year = as.integer(year))]

    l_leads <- attr(pm_object,"lead")
    dv <- attr(pm_object, "outcome")
    l_y_lead <- paste0(dv, "_l" ,l_leads)

    ## construct leading DVs (long-term impact)
    dt_yl <- dt_panel %>% copy %>%
        .[, (l_y_lead) := map(l_leads, ~shift(get(dv), -.x)), wbcode2]
        
    ## dt_yl %>% na.omit
    ## dt_yl[, .N, .(wbcode2)]

    dt_sets_wcov <- merge(dt_sets, dt_yl, by= c("wbcode2", "year")) %>%
        .[, setlen := .N, set]

    return(dt_sets_wcov)

}

dt_pm_long <- prep_pm_long(pm_object, adt(dem))


## OLS: for contemporaneous effect essentially
feols(y ~ type | set, dt_pm_long, weights = ~weight) # -1.29


## hmm doesn't work
    
rx <- feols(c(y_l0, y_l1, y_l2, y_l3, y_l4) ~ type | set, dt_pm_long, weights = ~weight) # -1.29

dt_mnl <- map(rx, ~tidy(.x) %>% adt) %>% rbindlist(idcol = "spec")


## try manual difference
dt_diff <- dt_pm_long[weight != 0][order(set, -weight)] %>%
    .[type == "crol", .(y_crol = mean(y)), set] %>%
    merge(dt_pm_long[type == "treat", .(wbcode2, set, y_treat = y)], by = "set") %>%
    .[, y_diff := y_treat - y_crol]


dt_diff[, mean(y_diff)] # also -1.29

## try boostrap

r_boot <- boot(dt_diff$y_diff, \(data, indices) mean(data[indices]), R = 1000)
r_boot$t %>% mean # still -1.21

## compare estimates: 
rbind(data.table(vlu = pe_results$bootstrapped.estimates[,1], src = "PM"),
      data.table(vlu = r_boot$t[,1], src = "mnl")) %>%
    ggplot(aes(x=vlu, color = src)) + geom_density()
## mnl is much wider
    

  
  
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


apply(pe_results$bootstrapped.estimates, 2, mean)

matched_df <- prepare_matched_data(pm_results, pm_data)

lm_est <- lm(y ~ dem, data = matched_df, weights = matched_df$matched_weights)

cat("PanelEstimate (pooled ATT):\n")
print(mean(pe_results$estimate))
cat("\nReplication (OLS):\n")
print(coef(lm_est)["dem"])
