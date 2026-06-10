library(data.table)

#' Match on treatment history using a two-tier exact matching strategy.
#'
#' @description
#' This function performs a hierarchical matching process:
#' 1. Tier 1: Exact matching on treatment history, nationality, and gender.
#' 2. Tier 2: Exact matching on treatment history only for units not matched in Tier 1.
#'
#' @param dt A data.table containing the panel data.
#' @param L Integer. The number of lags of treatment history to match on.
#' @param treatment_col String. The column name of the treatment variable.
#' @param history_col String. The column name for the cumulative history variable.
#' @param unit_id String. The column name of the unit identifier.
#' @param time_id String. The column name of the time identifier.
#' @param return_summary Boolean. If TRUE, returns a diagnostic summary; if FALSE, returns a table of matches.
#' @return A data.table containing either summary statistics or all matched (control, treated) pairs.
match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    history_col = "exn_pm_cum",
                                    unit_id = "PeopleID", time_id = "year", 
                                    return_summary = TRUE) {
  
  d <- prepare_data(dt, L, treatment_col, history_col, unit_id, time_id)
  
  # Tier 1: Exact Match (History + Nationality + Gender)
  d[, tier1_key := paste(hist_str, Nationality, Gender, sep = "_")]
  matches_t1 <- execute_tier_match(d, "tier1_key")
  
  # Tier 2: Exact Match (History only for orphans)
  orphaned_treated <- d[treat_bin == 1 & !(uid %in% matches_t1$uid)]
  matches_t2 <- execute_tier_match(orphaned_treated, "hist_str", controls_pool = d[treat_bin == 0])
  
  all_matches <- rbind(matches_t1, matches_t2, fill=TRUE)
  
  if (return_summary) return(calculate_summary(d, all_matches))
  return(all_matches)
}

#' Prepare data for matching
#' @param dt Input data.table.
#' @param L Number of lags.
#' @param t_col Name of treatment column.
#' @param h_col Name of history column.
#' @param uid Name of unit id column.
#' @param tid Name of time id column.
#' @return A prepared data.table with history signatures.
prepare_data <- function(dt, L, t_col, h_col, uid, tid) {
  d <- copy(dt)
  setnames(d, c(uid, tid), c("uid", "tid"))
  d[, treat_bin := as.integer(get(t_col) > 0)]
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(get(h_col), 1:L), by = uid]
  d <- d[complete.cases(d[, ..lag_cols])]
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  return(d)
}

#' Execute tier-specific matching
#' @param d The data.table to match within.
#' @param key_col Column name representing the exact match grouping.
#' @param controls_pool Optional pool of controls (if NULL, derives from d).
#' @return A data.table containing match counts per treated unit.
execute_tier_match <- function(d, key_col, controls_pool = NULL) {
  results <- list()
  years <- sort(unique(d$tid))
  for (yr in years) {
    ctrls <- if(is.null(controls_pool)) d[tid == yr & treat_bin == 0] else controls_pool[tid == yr]
    trtd <- d[tid == yr & treat_bin == 1]
    
    if (nrow(trtd) == 0 || nrow(ctrls) == 0) next
    
    setkeyv(trtd, key_col); setkeyv(ctrls, key_col)
    
    # Identify matches for EACH treated row using by=.EACHI
    counts <- ctrls[trtd, .(uid = i.uid, N = .N), on = key_col, nomatch = NULL, by = .EACHI]
    results[[as.character(yr)]] <- counts[N > 0]
  }
  return(rbindlist(results, fill=TRUE))
}

#' Calculate summary statistics for matches
#' @param d The full prepared dataset.
#' @param matches Data.table of matches (uid and match count N).
#' @return A summary data.table.
calculate_summary <- function(d, matches) {
  data.table(
    total_treated = uniqueN(d[treat_bin == 1, uid]),
    treated_matched = uniqueN(matches$uid),
    avg_matches_per_treated = mean(matches$N),
    total_matches = sum(matches$N)
  )
}
