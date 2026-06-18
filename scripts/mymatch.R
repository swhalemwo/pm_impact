library(data.table)
library(MatchIt)
library(cobalt)

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Match on treatment history using a two-tier exact matching strategy.
#'
#' @description
#' Performs a hierarchical matching process:
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

  ## Tier 1: exact match on history + Nationality + Gender
  d[, tier1_key := paste(hist_str, Nationality, Gender, sep = "_")]
  matches_t1 <- execute_tier_match(d, "tier1_key", return_type = return_type, tier = 1)

  ## Tier 2: exact match on history only (for orphans from Tier 1)
  matched_uid <- if (return_type == "summary") unique(matches_t1$uid) else unique(matches_t1$treated_id)
  orphaned_treated <- d[treat_bin == 1 & !(uid %in% matched_uid)]

  matches_t2 <- execute_tier_match(orphaned_treated, "hist_str",
                                   controls_pool = d[treat_bin == 0],
                                   return_type = return_type, tier = 2)

  all_matches <- rbind(matches_t1, matches_t2, fill = TRUE)

  if (return_type == "summary") return(calculate_summary(d, all_matches))
  return(all_matches)
}

#' Execute matching for a single year within a tier.
#'
#' @description
#' Given treated and control units in a single year, joins them on the
#' given key column and formats the result according to return_type.
#'
#' @param yr Integer year to process.
#' @param d data.table of all units (or orphaned treated units for Tier 2).
#' @param controls_pool Optional separate pool of controls (used for Tier 2).
#' @param key_col String name of the column to join on.
#' @param return_type One of "summary", "matches", "matches-info".
#' @param tier Integer tier number (1 or 2).
#' @return data.table of matched pairs, or NULL if no matches exist.
ematch_year_single <- function(yr, d, controls_pool, key_col, return_type, tier) {
  ## extract treated and control pools for this year
  ctrls <- if (is.null(controls_pool)) d[tid == yr & treat_bin == 0] else controls_pool[tid == yr]
  trtd <- d[tid == yr & treat_bin == 1]
  if (nrow(trtd) == 0 || nrow(ctrls) == 0) return(NULL)

  if (return_type == "summary") {
    ctrls[trtd, .(uid = i.uid, N = .N), on = key_col, nomatch = NULL, by = .EACHI] %>%
      .[N > 0]
  } else if (return_type == "matches-info") {
    ctrls[trtd, .(control_id = uid, treated_id = i.uid, year = yr,
                   match_key = get(key_col), match_tier = tier),
          on = key_col, nomatch = NULL, allow.cartesian = TRUE]
  } else {
    ctrls[trtd, .(control_id = uid, treated_id = i.uid, year = yr),
          on = key_col, nomatch = NULL, allow.cartesian = TRUE]
  }
}

#' Execute tier matching by applying year-by-year joins.
#'
#' @description
#' Iterates over years and calls ematch_year_single for each.
#' Combines all years into a single data.table.
#'
#' @param d data.table of units (must have tid and treat_bin columns).
#' @param key_col Column name for the exact-match key.
#' @param controls_pool Optional separate control pool.
#' @param return_type One of "summary", "matches", "matches-info".
#' @param tier Integer tier number.
#' @return data.table of matched pairs across all years.
execute_tier_match <- function(d, key_col, controls_pool = NULL,
                               return_type = "summary", tier = 1) {
  years <- sort(unique(d$tid))

  rbindlist(lapply(years, function(yr) {
    ematch_year_single(yr, d, controls_pool, key_col, return_type, tier)
  }), fill = TRUE)
}

#' Calculate summary statistics from matched pairs.
#'
#' @param d Prepared data (must have treat_bin and uid columns).
#' @param matches Summary-mode match output (with uid and N columns).
#' @return data.table with total_treated, treated_matched, avg_matches_per_treated, total_matches.
calculate_summary <- function(d, matches) {
  data.table(
    total_treated = uniqueN(d[treat_bin == 1, uid]),
    treated_matched = uniqueN(matches$uid),
    avg_matches_per_treated = if (uniqueN(matches$uid) > 0) {
      sum(matches$N) / uniqueN(matches$uid)
    } else 0,
    total_matches = sum(matches$N)
  )
}

#' Prepare panel data for matching by constructing lagged history strings.
#'
#' @description
#' Renames unit and time IDs for internal use, creates binary treatment
#' indicator, generates L lags of the cumulative history variable, and
#' concatenates them into a string key for exact matching.
#'
#' @param dt Raw panel data.table.
#' @param L Number of lags.
#' @param t_col Treatment column name.
#' @param h_col History (cumulative) column name.
#' @param uid Unit ID column name.
#' @param tid Time ID column name.
#' @return Modified copy of dt with uid, tid, treat_bin, lag_* and hist_str columns.
prepare_data <- function(dt, L, t_col, h_col, uid, tid) {
  d <- copy(dt)
  setnames(d, c(uid, tid), c("uid", "tid"))

  d[, treat_bin := as.integer(get(t_col) > 0)]

  lag_cols <- paste0("lag_", 1:L)
  d[, (lag_cols) := shift(get(h_col), 1:L), by = uid]
  d <- d[complete.cases(d[, ..lag_cols])]

  d[, hist_str := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_cols]
  d
}

#' Refine matches using Coarsened Exact Matching (CEM).
#'
#' @param dt_match Data.table of raw matches (control_id, treated_id, year).
#'   Must contain match_tier column (from matches-info mode).
#' @param dt_full The full panel data containing covariates.
#' @param covariates Character vector. Covariates for refinement.
#' @param cutpoints List. Cutpoints for CEM.
#' @return A list with:
#'   \item{data}{data.table with subclass column assigned}
#'   \item{matchit}{the raw matchit object for diagnostics}
refine_matches <- function(dt_match, dt_full,
                           covariates = c("BirthYear", "exn_all_cum"),
                           cutpoints = list(BirthYear = 8, exn_all_cum = 8)) {

  ## transform matches to long format (one row per control, one row per treated)
  dt_match_long <- rbind(
    dt_match[, .(PeopleID = control_id, year, treat = 0,
                 set = paste0(treated_id, "_", year))],
    dt_match[, head(.SD, 1), treated_id][, .(PeopleID = treated_id, year, treat = 1,
                                              set = paste0(treated_id, "_", year))]
  )

  ## merge covariates from full panel
  dt_match_long <- merge(dt_match_long, dt_full, by = c("PeopleID", "year"))

  ## apply CEM via matchit
  formula <- as.formula(paste("treat ~ set +", paste(covariates, collapse = " + ")))
  r_match <- matchit(formula, data = dt_match_long, method = "cem",
                     cutpoints = cutpoints)

  ## assign subclasses
  dt_match_long[, subclass := r_match$subclass]
  list(data = dt_match_long, matchit = r_match)
}

#' Run the full matching pipeline and return quality diagnostics.
#'
#' @description
#' Convenience wrapper running match_treatment_history (matches-info mode)
#' then refine_matches, then quality_report_matches. Use for parameter grid search.
#'
#' @param dt Panel data.table.
#' @param params Named list. Expected fields:
#'   L, treatment_col, history_col, unit_id, time_id (passed to match_treatment_history);
#'   covariates, cutpoints (passed to refine_matches).
#' @param dt_full Full panel data for covariate merging in refinement.
#' @param sample_n_treated Optional integer. If set and the number of unique
#'   matched treated units exceeds this, a random sample is used for the CEM stage.
#' @return A single-row data.table of quality indicators.
run_matching_pipeline <- function(dt, params, dt_full = dt,
                                  sample_n_treated = NULL) {

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

  ## subsample treated for CEM stage if requested
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
#' Stage-1 indicators:
#'   n_treated, n_matched, match_rate, pct_tier1, pct_tier2,
#'   avg_controls_per_treated, min_controls_per_treated
#'
#' Stage-2 indicators:
#'   n_after_cem, cem_match_rate, L1_before, L1_after (multivariate L1),
#'   and per-covariate SMD and variance ratio before/after.
#'
#' @param dt_prepared data.table after prepare_data() — used for n_treated denominator.
#' @param dt_match data.table from match_treatment_history(matches-info).
#' @param matchit_obj The matchit object from refine_matches()$matchit. NULL if no refinement.
#' @param refined_data The data slot from refine_matches()$data.
#' @param params Named list of parameters used for this run.
#' @return Single-row data.table.
quality_report_matches <- function(dt_prepared, dt_match, matchit_obj = NULL,
                                   refined_data = NULL, params = list()) {
  ## stage 1: exact matching coverage
  s1 <- qrm_s1(dt_prepared, dt_match, params)

  ## stage 2: CEM balance diagnostics
  s2 <- qrm_s2(matchit_obj, refined_data, params)

  ## combine with params as leading columns
  params_dt <- if (length(params) > 0) {
    as.data.table(lapply(params, function(x) {
      if (length(x) > 1) paste(x, collapse = ",") else x
    }))
  } else data.table()

  cbind(params_dt, s1, s2)
}

#' @rdname quality_report_matches
#' Stage-1 helper: exact matching coverage from dt_match
qrm_s1 <- function(dt_prepared, dt_match, params) {
  ## count treated in the population
  n_treated <- if ("treat_bin" %in% names(dt_prepared)) {
    uniqueN(dt_prepared[treat_bin == 1, uid])
  } else {
    t_col <- params$treatment_col %||% "exn_pm"
    uid_col <- params$unit_id %||% "PeopleID"
    uniqueN(dt_prepared[get(t_col) > 0, get(uid_col)])
  }

  if (nrow(dt_match) == 0) {
    return(data.table(n_treated = n_treated, n_matched = 0L, match_rate = 0,
                      pct_tier1 = NA_real_, pct_tier2 = NA_real_,
                      avg_controls_per_treated = NA_real_,
                      min_controls_per_treated = NA_integer_))
  }

  matched_treated <- unique(dt_match$treated_id)
  n_matched <- length(matched_treated)

  ## tier breakdown: one tier assignment per treated unit
  tier_counts <- if ("match_tier" %in% names(dt_match)) {
    dt_match[, .(tier = match_tier[1]), by = treated_id]
  } else NULL

  pct_tier1 <- if (!is.null(tier_counts)) mean(tier_counts$tier == 1) else NA_real_
  pct_tier2 <- if (!is.null(tier_counts)) mean(tier_counts$tier == 2) else NA_real_

  ## controls-per-treated
  ctrl_counts <- dt_match[, .(n_ctrl = uniqueN(control_id)), by = treated_id]

  data.table(n_treated = n_treated, n_matched = n_matched,
             match_rate = round(n_matched / n_treated, 4),
             pct_tier1 = round(pct_tier1, 4), pct_tier2 = round(pct_tier2, 4),
             avg_controls_per_treated = round(mean(ctrl_counts$n_ctrl), 2),
             min_controls_per_treated = min(ctrl_counts$n_ctrl))
}

#' Build SMD/VR column block for one covariate (helper for qrm_s2).
#'
#' @param v Covariate name (string).
#' @param smd_results data.table from fast_smd_vr with columns cov, smd_un, vr_un, smd_adj, vr_adj.
#' @return Single-row data.table with four columns: smd/vr before/after.
bal_cols_one_cov <- function(v, smd_results) {
  row <- smd_results[cov == v]
  if (nrow(row) == 0) {
    setNames(data.table(NA_real_, NA_real_, NA_real_, NA_real_),
             c(paste0("smd_", v, "_before"), paste0("vr_", v, "_before"),
               paste0("smd_", v, "_after"),  paste0("vr_", v, "_after")))
  } else {
    setNames(data.table(round(row$smd_un, 4), round(row$vr_un, 4),
                         round(row$smd_adj, 4), round(row$vr_adj, 4)),
             c(paste0("smd_", v, "_before"), paste0("vr_", v, "_before"),
               paste0("smd_", v, "_after"),  paste0("vr_", v, "_after")))
  }
}

#' @rdname quality_report_matches
#' Stage-2 helper: CEM balance diagnostics from matchit object
qrm_s2 <- function(matchit_obj, refined_data, params) {
  if (is.null(matchit_obj) || is.null(refined_data)) {
    return(data.table(n_after_cem = NA_integer_, cem_match_rate = NA_real_,
                      L1_before = NA_real_, L1_after = NA_real_))
  }

  w_after <- matchit_obj$weights
  treat_v <- refined_data$treat

  ## n after CEM and match rate
  n_after_cem <- sum(treat_v == 1 & w_after > 0)
  n_treated_in_cem <- sum(treat_v == 1)
  cem_match_rate <- if (n_treated_in_cem > 0) {
    round(n_after_cem / n_treated_in_cem, 4)
  } else NA_real_

  ## identify real covariates (not the set dummies from matchit)
  covs <- intersect(params$covariates %||% c("BirthYear", "exn_all_cum"),
                    names(refined_data))

  if (length(covs) == 0) {
    return(data.table(n_after_cem = as.integer(n_after_cem),
                      cem_match_rate = cem_match_rate,
                      L1_before = NA_real_, L1_after = NA_real_))
  }

  cov_df <- refined_data[, ..covs]

  ## L1 imbalance
  l1_results <- fast_l1_med_both(cov_df, treat_v, w_after,
                                 l1.n = 51L, l1.max.bin = 7L)

  ## per-covariate SMD + variance ratio
  smd_results <- fast_smd_vr(cov_df, treat_v, w_after)

  ## build SMD/VR columns — one per covariate
  bal_parts <- lapply(covs, bal_cols_one_cov, smd_results)

  data.table(n_after_cem = as.integer(n_after_cem),
             cem_match_rate = cem_match_rate,
             L1_before = round(l1_results$before, 4),
             L1_after  = round(l1_results$after, 4)) %>%
    cbind(do.call(cbind, bal_parts))
}

#' Generate one random coarsening and compute its L1 statistic.
#'
#' @description
#' Draws random bin counts for each covariate, bins them using
#' equal-spaced cutpoints, computes L1 under uniform weights.
#' Returns the L1 value and the binned columns for median selection.
#'
#' @param dt_loc data.table with .t (logical treat) and .w_u (uniform weights) columns.
#' @param cov_names Character vector of covariate names.
#' @param ranges List of ranges (min/max) per covariate.
#' @param l1.min.bin Minimum number of bins per coarsening.
#' @param l1.max.bin Maximum number of bins per coarsening.
#' @param nu Integer vector of unique value counts per covariate.
#' @param cut_cols Character vector of column names to store bin assignments.
#' @return Named list with l1_u (numeric) and bin_df (data.table of binned columns).
gen_one_coarsening <- function(dt_loc, cov_names, ranges, l1.min.bin, l1.max.bin,
                               nu, cut_cols) {
  nbins <- vapply(seq_along(cov_names), function(i) {
    rng <- seq(min(l1.min.bin, nu[i]), min(l1.max.bin, nu[i]))
    rng[sample.int(length(rng), 1L)]
  }, integer(1))

  ## assign binned columns in one shot via data.table := with lapply
  dt_loc[, (cut_cols) := lapply(seq_along(cov_names), function(i) {
    if (diff(ranges[[i]]) == 0) return(rep(1L, .N))
    as.integer(findInterval(get(cov_names[i]),
      seq(ranges[[i]][1], ranges[[i]][2], length.out = nbins[i] + 1),
      rightmost.closed = TRUE))
  })]

  ## compute L1 for this coarsening under uniform weights
  agg <- dt_loc[, .(wT = sum(.w_u * .t), wC = sum(.w_u * !.t)), by = cut_cols]
  l1_u <- sum(abs(agg$wT - agg$wC)) / 2

  ## store binned integer columns for median coarsening reuse
  bin_df <- copy(dt_loc[, ..cut_cols])
  dt_loc[, (cut_cols) := NULL]

  list(l1_u = l1_u, bin_df = bin_df)
}


#' Fast multivariate L1 imbalance (before+after), matching cobalt methodology.
#'
#' Implements the exact same algorithm as cobalt::bal.compute(..., stat = "l1.med"):
#' 1. Generate `l1.n` random coarsenings (bin counts per covariate drawn uniformly
#'    from [l1.min.bin, l1.max.bin]).
#' 2. For each coarsening, bin covariates via equal-spaced cutpoints and compute
#'    L1 = sum over cells |weighted_share_treated - weighted_share_control| / 2
#'    /under uniform weights/.
#' 3. Select the coarsening with the median L1 value.
#' 4. Compute both unweighted (before) and CEM-weighted (after) L1 on that single
#'    fixed coarsening.
#'
#' This matches cobalt's methodology exactly (fixed median coarsening).
#' The all-numeric data.table-vectorized inner loop makes this ~5x faster than
#' cobalt, scaling to >10M rows in a few seconds.
#'
#' @param cov_df data.frame or data.table of numeric covariates.
#' @param treat 0/1 or logical vector, length == nrow(cov_df).
#' @param weights_after Weight vector for the "after" (matched) state.
#' @param l1.n Number of random coarsenings to generate (default 51).
#' @param l1.min.bin Min number of bins per covariate per coarsening (default 2).
#' @param l1.max.bin Max number of bins per covariate per coarsening (default 12).
#' @param seed Optional RNG seed for reproducibility.
#' @return Named list with `before` and `after` entries (each numeric in [0,1]).
fast_l1_med_both <- function(cov_df, treat, weights_after,
                             l1.n = 51L, l1.min.bin = 2L, l1.max.bin = 12L,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  dt_loc <- as.data.table(cov_df)
  cov_names <- copy(names(dt_loc))
  n <- length(treat)
  is_t <- as.logical(treat)

  ## uniform weights for finding median coarsening (step 2)
  w_u <- rep(1, n)
  w_u[is_t] <- w_u[is_t] / sum(w_u[is_t])
  w_u[!is_t] <- w_u[!is_t] / sum(w_u[!is_t])

  dt_loc[, `:=`(.t = is_t, .w_u = w_u)]

  ranges <- lapply(cov_names, function(c) range(dt_loc[[c]], na.rm = TRUE))
  nu <- vapply(cov_names, function(c) uniqueN(dt_loc[[c]]), integer(1))
  cut_cols <- paste0(".b_", cov_names)

  ## step 1-2: generate l1.n coarsenings
  l1_results <- lapply(seq_len(l1.n), function(k) {
    gen_one_coarsening(dt_loc, cov_names, ranges, l1.min.bin, l1.max.bin, nu, cut_cols)
  })

  l1s_u <- vapply(l1_results, `[[`, numeric(1), "l1_u")
  bin_ids <- lapply(l1_results, `[[`, "bin_df")

  ## step 3: pick median coarsening
  median_idx <- which(l1s_u == sort(l1s_u, partial = ceiling(l1.n / 2))[ceiling(l1.n / 2)])[1]
  fixed_cuts <- bin_ids[[median_idx]]

  ## step 4: compute before (uniform) and after (weighted) on fixed coarsening
  fixed_cuts[, bin_id := do.call(paste, c(.SD, sep = "_"))]
  w_o <- weights_after
  w_o[is_t] <- w_o[is_t] / sum(w_o[is_t])
  w_o[!is_t] <- w_o[!is_t] / sum(w_o[!is_t])

  dt_fixed <- data.table(bin = fixed_cuts$bin_id, .wu = w_u, .wo = w_o, .t = is_t)
  agg_before <- dt_fixed[, .(wT = sum(.wu * .t), wC = sum(.wu * !.t)), by = bin]
  agg_after  <- dt_fixed[, .(wT = sum(.wo * .t), wC = sum(.wo * !.t)), by = bin]

  list(before = sum(abs(agg_before$wT - agg_before$wC)) / 2,
       after  = sum(abs(agg_after$wT  - agg_after$wC))  / 2)
}


#' Fast multivariate L1 imbalance — single-weight wrapper.
#'
#' Wrapper around fast_l1_med_both for backward compatibility.
#' Returns the entry matching the weight type (before/unweighted or after/weighted).
#'
#' @inheritParams fast_l1_med_both
#' @param weights non-negative weight vector. Use rep(1, n) for unmatched baseline.
#' @return Numeric scalar in [0, 1].
fast_l1_med <- function(cov_df, treat, weights,
                        l1.n = 101L, l1.min.bin = 2L, l1.max.bin = 12L,
                        seed = NULL) {
  both <- fast_l1_med_both(cov_df, treat, weights, l1.n, l1.min.bin, l1.max.bin, seed)
  if (isTRUE(all(weights == 1))) both[["before"]] else both[["after"]]
}


#' Compute SMD and variance ratio for a single covariate.
#'
#' @description
#' Computes unweighted and CEM-weighted standardized mean difference and
#' variance ratio for one covariate. Matches cobalt defaults:
#' s.d.denom = "pooled" (unweighted pooled SD).
#'
#' @param cn Covariate name (string).
#' @param d data.table with .treat (logical) and .w (weights) columns.
#' @return Named list: cn, smd_un, vr_un, smd_adj, vr_adj.
smd_vr_one_cov <- function(cn, d) {
  vT <- d[.treat == TRUE,  get(cn)]
  vC <- d[.treat == FALSE, get(cn)]

  ## unweighted (before)
  mT_u <- mean(vT)
  mC_u <- mean(vC)
  vT_u <- var(vT)
  vC_u <- var(vC)
  sd_p <- sqrt((vT_u + vC_u) / 2)

  ## weighted (after) per cobalt formula:
  ## weighted variance = sum(w_norm * (x - mean)^2) / (1 - sum(w_norm^2))
  wT <- d[.treat == TRUE,  .w]
  wC <- d[.treat == FALSE, .w]
  wTn <- wT / sum(wT)
  wCn <- wC / sum(wC)

  mT_w <- sum(wTn * vT)
  mC_w <- sum(wCn * vC)

  vT_w <- sum(wTn * (vT - mT_w)^2) / (1 - sum(wTn^2))
  vC_w <- sum(wCn * (vC - mC_w)^2) / (1 - sum(wCn^2))

  list(cn,
       smd_un = (mT_u - mC_u) / sd_p,
       vr_un  = vT_u / vC_u,
       smd_adj = (mT_w - mC_w) / sd_p,
       vr_adj  = vT_w / vC_w)
}

#' Fast per-covariate SMD + variance ratio (~17x faster than cobalt::bal.tab).
#'
#' @param cov_df data.frame or data.table of numeric covariates.
#' @param treat 0/1 or logical vector, length == nrow(cov_df).
#' @param weights non-negative weight vector (CEM weights).
#' @return data.table with columns: cov, smd_un, vr_un, smd_adj, vr_adj.
fast_smd_vr <- function(cov_df, treat, weights) {
  d <- as.data.table(cov_df)
  cov_names <- copy(names(d))
  d[, `:=`(.treat = as.logical(treat), .w = weights)]

  ## compute SMD and VR for each covariate via named helper
  out <- lapply(cov_names, smd_vr_one_cov, d)

  rbindlist(out) %>%
    setnames(c("cov", "smd_un", "vr_un", "smd_adj", "vr_adj"))
}
