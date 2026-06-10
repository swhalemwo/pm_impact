library(data.table)

match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    unit_id = "PeopleID", time_id = "year", 
                                    binary_treatment = TRUE) {
  
  d <- copy(dt)
  setorderv(d, c(unit_id, time_id))
  
  if (binary_treatment) {
    d[, treat_bin := as.integer(get(treatment_col) > 0)]
  } else {
    d[, treat_bin := get(treatment_col)]
  }
  
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(treat_bin, 1:L), by = unit_id]
  d <- d[complete.cases(d[, ..lag_cols])]
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  
  treated <- d[treat_bin == 1]
  controls <- d[treat_bin == 0]
  
  setkey(treated, year, hist_str)
  setkey(controls, year, hist_str)
  
  # Allow cartesian product because one treated unit can match multiple controls
  matched <- controls[treated, .(control_id = PeopleID, treated_id = i.PeopleID, year, hist_str), 
                      nomatch = NULL, 
                      on = .(year, hist_str), 
                      allow.cartesian = TRUE]
  
  return(matched)
}
