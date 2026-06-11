library(data.table)
library(MatchIt)
library(cobalt)

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
#'   Must contain match_tier column (i.e., produced with return_type = "matches-info").
#' @param dt_full The full panel data containing covariates.
#' @param covariates Character vector. Covariates for refinement.
#' @param cutpoints List. Cutpoints for CEM.
#' @return A list with:
#'   \item{data}{data.table with subclass column assigned}
#'   \item{matchit}{the raw matchit object for diagnostics}
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
  return(list(data = dt_match_long, matchit = r_match))
}

#' Run the full matching pipeline and return quality diagnostics.
#'
#' @description
#' Convenience wrapper that runs match_treatment_history (in matches-info mode)
#' followed by refine_matches, then calls quality_report_matches. Use this
#' when doing parameter grid search.
#'
#' @param dt Panel data.table.
#' @param params Named list of matching parameters. Expected fields:
#'   L, treatment_col, history_col, unit_id, time_id (passed to match_treatment_history);
#'   covariates, cutpoints (passed to refine_matches).
#' @param dt_full Full panel data for covariate merging in refinement.
#' @param sample_n_treated Optional integer. If set and the number of unique
#'   matched treated units exceeds this, a random sample of that many treated
#'   units (with all their controls) is used for the CEM stage. Stage-1
#'   indicators are always computed on the full match set. Default NULL = no sampling.
#'   Recommended ~1000-2000 for large panels — bal.compute/bal.tab scale poorly
#'   beyond a few million rows.
#' @return A single-row data.table of quality indicators (suitable for rbinding across runs).
run_matching_pipeline <- function(dt, params, dt_full = dt, sample_n_treated = NULL) {

  dt_match <- match_treatment_history(
    dt,
    L             = params$L,
    treatment_col = params$treatment_col,
    history_col   = params$history_col,
    unit_id       = params$unit_id,
    time_id       = params$time_id,
    return_type   = "matches-info"
  )

  if (nrow(dt_match) == 0) {
    warning("run_matching_pipeline: no matches found, returning NA quality row.")
    return(quality_report_matches(dt, dt_match, matchit_obj = NULL,
                                  refined_data = NULL, params = params))
  }

  # subsample treated for CEM stage if requested
  dt_match_for_cem <- dt_match
  if (!is.null(sample_n_treated)) {
    treated_all <- unique(dt_match$treated_id)
    if (length(treated_all) > sample_n_treated) {
      keep <- sample(treated_all, sample_n_treated)
      dt_match_for_cem <- dt_match[treated_id %in% keep]
    }
  }

  refined <- refine_matches(dt_match_for_cem, dt_full,
                            covariates = params$covariates,
                            cutpoints  = params$cutpoints)

  quality_report_matches(dt, dt_match,
                         matchit_obj  = refined$matchit,
                         refined_data = refined$data,
                         params       = params)
}

#' Compute a flat single-row quality report for one parameter combination.
#'
#' @description
#' Extracts stage-1 (exact matching) and stage-2 (CEM) quality indicators into
#' a single data.table row. Designed so that many rows can be rbind-ed into a
#' comparison table when doing parameter grid search.
#'
#' Stage-1 indicators (from dt_match):
#'   n_treated, n_matched, match_rate, pct_tier1, pct_tier2,
#'   avg_controls_per_treated, min_controls_per_treated
#'
#' Stage-2 indicators (from matchit_obj + refined_data via cobalt):
#'   n_after_cem, cem_match_rate, L1_before, L1_after (multivariate L1 via
#'   cobalt::bal.compute, stat = "l1.med"),
#'   and per-covariate SMD and variance ratio before/after (via cobalt::bal.tab
#'   on the covariate frame directly — bypasses the explosion of `set` dummies
#'   that summary.matchit would otherwise report).
#'
#' @param dt_prepared data.table after prepare_data() — used for n_treated denominator.
#'   Alternatively pass the raw panel; the function only needs uniqueN of treated units.
#' @param dt_match data.table from match_treatment_history(..., return_type="matches-info").
#' @param matchit_obj The matchit object from refine_matches()$matchit. NULL if no refinement.
#' @param refined_data The data slot from refine_matches()$data. Required when
#'   matchit_obj is provided (used to compute balance on the real covariates only).
#' @param params Named list of parameters used for this run (appended as columns).
#'   The `covariates` field tells us which covariates to extract balance stats for.
#' @return Single-row data.table.
quality_report_matches <- function(dt_prepared, dt_match, matchit_obj = NULL,
                                   refined_data = NULL, params = list()) {

  # --- Stage 1: exact matching coverage ---
  # dt_prepared may be raw panel or already prepare_data() output; handle both
  if ("treat_bin" %in% names(dt_prepared)) {
    n_treated <- uniqueN(dt_prepared[treat_bin == 1, uid])
  } else {
    t_col <- if (!is.null(params$treatment_col)) params$treatment_col else "exn_pm"
    n_treated <- uniqueN(dt_prepared[get(t_col) > 0, get(if (!is.null(params$unit_id)) params$unit_id else "PeopleID")])
  }

  if (nrow(dt_match) == 0) {
    n_matched  <- 0L
    match_rate <- 0
    pct_tier1  <- NA_real_
    pct_tier2  <- NA_real_
    avg_ctrl   <- NA_real_
    min_ctrl   <- NA_real_
  } else {
    matched_treated <- unique(dt_match$treated_id)
    n_matched  <- length(matched_treated)
    match_rate <- n_matched / n_treated

    # tier breakdown (only available in matches-info mode, which we always use here)
    if ("match_tier" %in% names(dt_match)) {
      tier_counts <- dt_match[, .(n_pairs = .N), by = .(treated_id, match_tier)]
      tier_counts <- tier_counts[, .(tier = match_tier[1]), by = treated_id]  # one row per treated
      pct_tier1 <- mean(tier_counts$tier == 1)
      pct_tier2 <- mean(tier_counts$tier == 2)
    } else {
      pct_tier1 <- NA_real_
      pct_tier2 <- NA_real_
    }

    # controls-per-treated
    ctrl_counts <- dt_match[, .(n_ctrl = uniqueN(control_id)), by = treated_id]
    avg_ctrl <- mean(ctrl_counts$n_ctrl)
    min_ctrl <- min(ctrl_counts$n_ctrl)
  }

  s1 <- data.table(
    n_treated              = n_treated,
    n_matched              = n_matched,
    match_rate             = round(match_rate, 4),
    pct_tier1              = round(pct_tier1, 4),
    pct_tier2              = round(pct_tier2, 4),
    avg_controls_per_treated = round(avg_ctrl, 2),
    min_controls_per_treated = min_ctrl
  )

  # --- Stage 2: CEM balance diagnostics (via cobalt) ---
  if (is.null(matchit_obj) || is.null(refined_data)) {
    s2 <- data.table(
      n_after_cem    = NA_integer_,
      cem_match_rate = NA_real_,
      L1_before      = NA_real_,
      L1_after       = NA_real_
    )
  } else {
    w_after  <- matchit_obj$weights
    treat_v  <- refined_data$treat

    # n after CEM: treated units with positive weight
    n_after_cem    <- sum(treat_v == 1 & w_after > 0)
    # cem_match_rate: of the treated fed INTO cem, what % survived
    n_treated_in_cem <- sum(treat_v == 1)
    cem_match_rate <- if (n_treated_in_cem > 0) round(n_after_cem / n_treated_in_cem, 4) else NA_real_

    # Identify the real covariates (from params; fall back to defaults of refine_matches)
    covs <- if (!is.null(params$covariates)) params$covariates else c("BirthYear", "exn_all_cum")
    covs <- intersect(covs, names(refined_data))

    if (length(covs) == 0) {
      L1_before <- NA_real_
      L1_after  <- NA_real_
      bal_before <- data.table()
      bal_after  <- data.table()
    } else {
      cov_df <- refined_data[, ..covs]

      # Multivariate L1 (l1.med) via cobalt — unweighted vs CEM weights
      L1_before <- tryCatch(
        round(bal.compute(cov_df, treat = treat_v,
                          weights = rep(1, length(treat_v)), stat = "l1.med"), 4),
        error = function(e) NA_real_
      )
      L1_after <- tryCatch(
        round(bal.compute(cov_df, treat = treat_v,
                          weights = w_after, stat = "l1.med"), 4),
        error = function(e) NA_real_
      )


      ## bal.compute(cov_df, treat = treat_v, weights = w_after, stat = "mahalanobis")

      # Per-covariate SMD + variance ratio, before and after
      bt <- tryCatch(
        bal.tab(cov_df, treat = treat_v, weights = w_after, un = TRUE,
                stats = c("mean.diffs", "variance.ratios"),
                binary = "std", continuous = "std"),
        error = function(e) NULL
      )

      if (!is.null(bt) && !is.null(bt$Balance)) {
        bal <- bt$Balance
        # bal is a data.frame with rownames = covariate names, columns:
        #   Diff.Un, V.Ratio.Un, Diff.Adj, V.Ratio.Adj (and Type)

        get_val <- function(rn, col) {
          if (rn %in% rownames(bal) && col %in% colnames(bal)) round(bal[rn, col], 4) else NA_real_
        }
        
        before_list <- lapply(covs, function(v) {
          setNames(data.table(get_val(v, "Diff.Un"), get_val(v, "V.Ratio.Un")),
                   c(paste0("smd_", v, "_before"), paste0("vr_", v, "_before")))
        })
        after_list <- lapply(covs, function(v) {
          setNames(data.table(get_val(v, "Diff.Adj"), get_val(v, "V.Ratio.Adj")),
                   c(paste0("smd_", v, "_after"), paste0("vr_", v, "_after")))
        })
        
        bal_before <- do.call(cbind, before_list)
        bal_after  <- do.call(cbind, after_list)
        
      } else {
        bal_before <- data.table()
        bal_after  <- data.table()
      }
    }

    s2 <- data.table(n_after_cem    = as.integer(n_after_cem),
                     cem_match_rate = cem_match_rate,
                     L1_before      = L1_before,
                     L1_after       = L1_after)
    
    if (nrow(bal_before) > 0) s2 <- cbind(s2, bal_before)
    if (nrow(bal_after)  > 0) s2 <- cbind(s2, bal_after)
  }

  # --- Combine with params ---
  params_dt <- if (length(params) > 0) as.data.table(lapply(params, function(x) {
    if (length(x) > 1) paste(x, collapse = ",") else x
  })) else data.table()

  cbind(params_dt, s1, s2)
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
