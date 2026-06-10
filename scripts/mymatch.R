library(data.table)

match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    unit_id = "PeopleID", time_id = "year", 
                                    binary_treatment = TRUE) {
  
  d <- copy(dt)
  # Standardize column names for join
  setnames(d, c(unit_id, time_id), c("uid", "tid"))
  
  if (binary_treatment) {
    d[, treat_bin := as.integer(get(treatment_col) > 0)]
  } else {
    d[, treat_bin := get(treatment_col)]
  }
  
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(treat_bin, 1:L), by = uid]
  d <- d[complete.cases(d[, ..lag_cols])]
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  
  years <- sort(unique(d$tid))
  results <- list()
  
  ## years <- 2010

  for (yr in years) {
    d_yr <- d[tid == yr]
    treated <- d_yr[treat_bin == 1]
    controls <- d_yr[treat_bin == 0]
    
    if (nrow(treated) > 0 && nrow(controls) > 0) {
      setkey(treated, hist_str)
      setkey(controls, hist_str)
      
      # Join without dynamic column resolution
      matched <- controls[treated, .(control_id = uid, treated_id = i.uid, year = yr, hist_str), 
                          nomatch = NULL, 
                          on = .(hist_str), 
                          allow.cartesian = TRUE]

      ## matched[, .N, hist_str]
      ## matched[, .(uniqueN(treated_id)), hist_str]
      
      results[[as.character(yr)]] <- matched
    }
  }
  
  return(rbindlist(results))
}
