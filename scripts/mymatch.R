library(data.table)

#' Match on treatment history by year with summary or full return option
match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    history_col = "exn_pm_cum",
                                    unit_id = "PeopleID", time_id = "year", 
                                    return_summary = TRUE) {
  
  d <- copy(dt)
  setnames(d, c(unit_id, time_id), c("uid", "tid"))
  
  d[, treat_bin := as.integer(get(treatment_col) > 0)]
  d[, hist_val := get(history_col)]
  
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(hist_val, 1:L), by = uid]
  d <- d[complete.cases(d[, ..lag_cols])]
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  
  years <- sort(unique(d$tid))
  
  if (return_summary) {
    total_matches_count <- 0
    treated_set <- c()
    
    for (yr in years) {
      d_yr <- d[tid == yr]
      treated <- d_yr[treat_bin == 1]
      controls <- d_yr[treat_bin == 0]
      
      if (nrow(treated) > 0 && nrow(controls) > 0) {
        # Create an index for faster lookup of control counts
        ctrl_counts <- controls[, .N, by = hist_str]
        setkey(ctrl_counts, hist_str)
        
        # Merge this count into treated units
        treated_counts <- ctrl_counts[treated, on = .(hist_str), nomatch = NULL]
        
        # Filter where N > 0
        matched_treated <- treated_counts[N > 0]
        
        total_matches_count <- total_matches_count + sum(matched_treated$N)
        treated_set <- c(treated_set, matched_treated$uid)
      }
    }
    
    summary <- data.table(
      total_treated_with_treatment = uniqueN(d[treat_bin == 1, uid]),
      treated_matched = uniqueN(treated_set),
      avg_matches_per_treated = ifelse(uniqueN(treated_set) > 0, total_matches_count / uniqueN(treated_set), 0),
      total_matches = total_matches_count
    )
    return(summary)
    
  } else {
    results <- list()
    for (yr in years) {
      d_yr <- d[tid == yr]
      treated <- d_yr[treat_bin == 1]
      controls <- d_yr[treat_bin == 0]
      
      if (nrow(treated) > 0 && nrow(controls) > 0) {
        setkey(treated, hist_str)
        setkey(controls, hist_str)
        
        matched <- controls[treated, .(control_id = uid, treated_id = i.uid, year = yr), 
                            nomatch = NULL, 
                            on = .(hist_str), 
                            allow.cartesian = TRUE]
        results[[as.character(yr)]] <- matched
      }
    }
    return(rbindlist(results))
  }
}
