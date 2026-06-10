library(data.table)
library(MatchIt)

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
#' @param return_type Option: "summary", "matches", "matches-info".
#' @return A data.table containing either summary statistics or all matched pairs.
match_treatment_history <- function(dt, L, treatment_col = "exn_pm", 
                                    history_col = "exn_pm_cum",
                                    unit_id = "PeopleID", time_id = "year", 
                                    return_type = "summary") {
  
  d <- prepare_data(dt, L, treatment_col, history_col, unit_id, time_id)
  
  # Tier 1: Exact Match (History + Nationality + Gender)
  d[, tier1_key := paste(hist_str, Nationality, Gender, sep = "_")]
  matches_t1 <- execute_tier_match(d, "tier1_key", return_type = return_type, tier = 1)
  
  # Tier 2: Exact Match (History only for orphans)
  matched_uid <- if(return_type == "summary") unique(matches_t1$uid) else unique(matches_t1$treated_id)
  orphaned_treated <- d[treat_bin == 1 & !(uid %in% matched_uid)]
  
  matches_t2 <- execute_tier_match(orphaned_treated, "hist_str", controls_pool = d[treat_bin == 0], 
                                   return_type = return_type, tier = 2)
  
  all_matches <- rbind(matches_t1, matches_t2, fill=TRUE)
  
  if (return_type == "summary") return(calculate_summary(d, all_matches))
  return(all_matches)
}

#' Refine matches using Coarsened Exact Matching (CEM)
#'
#' @param dt_match Data.table of raw matches (control_id, treated_id, year).
#' @param dt_full The full panel data containing covariates.
#' @param covariates Character vector. Covariates for refinement.
#' @param cutpoints List. Cutpoints for CEM.
#' @return A data.table with an added 'subclass' column for refinement.
refine_matches <- function(dt_match, dt_full, covariates = c("BirthYear", "exn_all_cum"), 
                           cutpoints = list(BirthYear = 8, exn_all_cum = 8)) {
  
  # 1. Transform matches to long format
  dt_match_long <- rbind(
    dt_match[, .(PeopleID = control_id, year, treat = 0, set = paste0(treated_id, "_", year))],
    dt_match[, head(.SD, 1), treated_id][, .(PeopleID = treated_id, year, treat = 1, set = paste0(treated_id, "_", year))]
  )
  
  # 2. Add covariates
  dt_match_long <- merge(dt_match_long, dt_full, by = c("PeopleID", "year"))
  
  # 3. Apply CEM via matchit
  formula <- as.formula(paste("treat ~ set +", paste(covariates, collapse = " + ")))
  r_match <- matchit(formula, data = dt_match_long, method = "cem", cutpoints = cutpoints)
  
  # 4. Assign subclasses
  dt_match_long[, subclass := r_match$subclass]
  return(dt_match_long)
}

# --- Helper Functions ---

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

execute_tier_match <- function(d, key_col, controls_pool = NULL, return_type = "summary", tier = 1) {
  results <- list()
  years <- sort(unique(d$tid))
  for (yr in years) {
    ctrls <- if(is.null(controls_pool)) d[tid == yr & treat_bin == 0] else controls_pool[tid == yr]
    trtd <- d[tid == yr & treat_bin == 1]
    if (nrow(trtd) == 0 || nrow(ctrls) == 0) next
    setkeyv(trtd, key_col); setkeyv(ctrls, key_col)
    
    if (return_type == "summary") {
      counts <- ctrls[trtd, .(uid = i.uid, N = .N), on = key_col, nomatch = NULL, by = .EACHI]
      results[[as.character(yr)]] <- counts[N > 0]
    } else {
      if (return_type == "matches-info") {
        matched <- ctrls[trtd, list(control_id = uid, treated_id = i.uid, year = yr, match_key = get(key_col), match_tier = tier), 
                         on = key_col, nomatch = NULL, allow.cartesian = TRUE]
      } else {
        matched <- ctrls[trtd, list(control_id = uid, treated_id = i.uid, year = yr), 
                         on = key_col, nomatch = NULL, allow.cartesian = TRUE]
      }
      results[[as.character(yr)]] <- matched
    }
  }
  return(rbindlist(results, fill=TRUE))
}

calculate_summary <- function(d, matches) {
  data.table(
    total_treated = uniqueN(d[treat_bin == 1, uid]),
    treated_matched = uniqueN(matches$uid),
    avg_matches_per_treated = ifelse(uniqueN(matches$uid) > 0, sum(matches$N) / uniqueN(matches$uid), 0),
    total_matches = sum(matches$N)
  )
}
