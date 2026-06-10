library(data.table)

# --- Core Matching Logic ---

match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    history_col = "exn_pm_cum",
                                    unit_id = "PeopleID", time_id = "year", 
                                    return_summary = TRUE) {
  
  d <- prepare_data(dt, L, treatment_col, history_col, unit_id, time_id)
  
  # Check if we have enough data - maybe L is too big for a single year?
  # L=3 requires 4 years of data (t, t-1, t-2, t-3). 
  # If tid=2010, we need 2007-2010.
  if(nrow(d) == 0) {
    # If empty, let's try a smaller L for the test or warn
    warning("Data is empty after processing. Ensure L lags are available.")
    return(data.table())
  }

  # Tier 1: Exact Match (History + Nationality + Gender)
  d[, tier1_key := paste(hist_str, Nationality, Gender, sep = "_")]
  matches_t1 <- execute_tier_match(d, "tier1_key", return_summary = return_summary)
  
  matched_uid <- if(return_summary) unique(matches_t1$uid) else unique(matches_t1$treated_id)
  
  orphaned_treated <- d[treat_bin == 1 & !(uid %in% matched_uid)]
  matches_t2 <- execute_tier_match(orphaned_treated, "hist_str", controls_pool = d[treat_bin == 0], return_summary = return_summary)
  
  all_matches <- rbind(matches_t1, matches_t2, fill=TRUE)
  
  if (return_summary) return(calculate_summary(d, all_matches))
  return(all_matches)
}

# --- Helper Functions ---

prepare_data <- function(dt, L, t_col, h_col, uid, tid) {
  d <- copy(dt)
  setnames(d, c(uid, tid), c("uid", "tid"))
  d[, treat_bin := as.integer(get(t_col) > 0)]
  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(get(h_col), 1:L), by = uid]
  # We return even if incomplete to see what's happening
  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  return(d)
}

execute_tier_match <- function(d, key_col, controls_pool = NULL, return_summary = TRUE) {
  results <- list()
  years <- sort(unique(d$tid))
  for (yr in years) {
    ctrls <- if(is.null(controls_pool)) d[tid == yr & treat_bin == 0] else controls_pool[tid == yr]
    trtd <- d[tid == yr & treat_bin == 1]
    
    if (nrow(trtd) == 0 || nrow(ctrls) == 0) next
    
    setkeyv(trtd, key_col); setkeyv(ctrls, key_col)
    
    if (return_summary) {
      counts <- ctrls[trtd, .(uid = i.uid, N = .N), on = key_col, nomatch = NULL, by = .EACHI]
      results[[as.character(yr)]] <- counts[N > 0]
    } else {
      matched <- ctrls[trtd, .(control_id = uid, treated_id = i.uid, year = yr), 
                       on = key_col, nomatch = NULL, allow.cartesian = TRUE]
      results[[as.character(yr)]] <- matched
    }
  }
  return(rbindlist(results, fill=TRUE))
}

calculate_summary <- function(d, matches) {
  data.table(
    total_treated = uniqueN(d[treat_bin == 1, uid]),
    treated_matched = uniqueN(matches$uid),
    avg_matches_per_treated = ifelse(uniqueN(matches$uid) > 0, mean(matches$N), 0),
    total_matches = sum(matches$N)
  )
}
