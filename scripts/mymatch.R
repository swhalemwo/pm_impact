library(data.table)

#' Match on treatment history
#' 
#' @param dt A data.table with unit.id, time.id, and treatment columns.
#' @param L Number of lags to match on.
#' @param treatment_col The name of the treatment column.
#' @param unit_id The name of the unit id column.
#' @param time_id The name of the time id column.
#' @param binary_treatment If TRUE, treats exn_pm > 0 as treatment.
#' @return A data.table of matched sets.
match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    unit_id = "PeopleID", time_id = "year", 
                                    binary_treatment = TRUE) {
  
  # Work on a copy
  d <- copy(dt)
  setorderv(d, c(unit_id, time_id))
  
  # 1. Binarize
  if (binary_treatment) {
    d[, treat_bin := as.integer(get(treatment_col) > 0)]
  } else {
    d[, treat_bin := get(treatment_col)]
  }
  
  # 2. Create L lag history columns
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(treat_bin, 1:L), by = unit_id]
  
  # Remove rows where history is incomplete
  d <- d[complete.cases(d[, ..lag_cols])]
  
  # 3. Create history string
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  
  # 4. Join
  treated <- d[treat_bin == 1]
  controls <- d[treat_bin == 0]
  
  setkey(treated, year, hist_str)
  setkey(controls, year, hist_str)
  
  # Using allow.cartesian=TRUE to allow multiple controls per treated unit
  matched <- controls[treated, .(control_id = PeopleID, treated_id = i.PeopleID, year, hist_str), 
                      nomatch = NULL, 
                      on = .(year, hist_str),
                      allow.cartesian = TRUE]
  
  return(matched)
}
