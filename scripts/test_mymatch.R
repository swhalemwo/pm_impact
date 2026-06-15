library(data.table)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Load small test data (3,000 artists, 58K rows — subset with 1500 PM-treated + 1500 never-treated)
dt <- fread("~/Dropbox/phd/papers/foster/dt_test_smol.csv")

# 1. Test Matching (Summary Mode)
print("--- Testing Summary Mode ---")
summary_res <- match_treatment_history(dt, L = 3, return_type = "summary")
print(summary_res)

# 2. Test Full Match Mode
print("--- Testing Full Match Mode ---")
full_matches <- match_treatment_history(dt, L = 3, return_type = "matches")
print(paste("Full matches rows:", nrow(full_matches)))


if(all(c("control_id", "treated_id", "year") %in% names(full_matches)) && 
   !any(c("match_key", "match_tier") %in% names(full_matches))) {
  print("SUCCESS: Matches mode returns only base columns.")
} else {
  stop("FAILURE: Matches mode returns incorrect columns.")
}

# 3. Consistency Check (Replicate Summary from Full)
print("--- Running Rigorous Consistency Check ---")
full_summary <- data.table(
  total_treated = uniqueN(dt[exn_pm > 0, PeopleID]),
  treated_matched = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = mean(full_matches[, .N, by = treated_id]$N),
  total_matches = nrow(full_matches)
)
stopifnot(all.equal(as.numeric(summary_res$total_matches), as.numeric(full_summary$total_matches)))
stopifnot(all.equal(as.numeric(summary_res$treated_matched), as.numeric(full_summary$treated_matched)))
print("SUCCESS: Consistency check passed.")

# 4. Test Matches-Info Mode
print("--- Testing Matches-Info Mode ---")
matches_info <- match_treatment_history(dt, L = 3, return_type = "matches-info")
stopifnot(all(c("match_key", "match_tier") %in% names(matches_info)))
print("SUCCESS: Matches-Info returns debugging columns.")

# 5. Test Refinement (CEM) on Full Data
print("--- Testing Refinement (CEM) on Full Data ---")
# refine_matches requires matches-info mode (for match_tier column)
refined <- refine_matches(matches_info, dt,
                          covariates = c("BirthYear", "exn_all_cum_ihs"),
                          cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))

if (!is.list(refined) || !all(c("data", "matchit") %in% names(refined))) {
  stop("FAILURE: refine_matches must return a list with 'data' and 'matchit'.")
}
if (!"subclass" %in% names(refined$data)) {
  stop("FAILURE: refined$data missing 'subclass' column.")
}
if (!inherits(refined$matchit, "matchit")) {
  stop("FAILURE: refined$matchit is not a matchit object.")
}
print("SUCCESS: Refinement returned list with data (subclass) and matchit object.")
print(head(refined$data))

# 6. Test number of tiers
dt_tiercnts <- matches_info[, .(nbr_tiers = uniqueN(match_tier)), treated_id]
if (max(dt_tiercnts$nbr_tiers) > 1) {
    stop("some treated are in multiple tiers")
} else {
    print("SUCCESS: all treated are in only one tier")
}

# 7. Test quality_report_matches() directly
# Note: on large match sets (millions of rows) bal.compute/bal.tab are slow.
# We subsample to a manageable number of treated for the CEM stage,
# but keep stage-1 stats computed on the full matches_info.
print("--- Testing quality_report_matches ---")
set.seed(1)
treated_sample <- sample(unique(matches_info$treated_id),
                         min(500, uniqueN(matches_info$treated_id)))
mi_small <- matches_info[treated_id %in% treated_sample]
refined_small <- refine_matches(mi_small, dt,
                                covariates = c("BirthYear", "exn_all_cum_ihs"),
                                cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))

qr <- quality_report_matches(
  dt_prepared  = dt,
  dt_match     = matches_info,        # full set -> stage-1 indicators reflect everything
  matchit_obj  = refined_small$matchit,
  refined_data = refined_small$data,
  params       = list(L = 3, treatment_col = "exn_pm",
                      covariates = c("BirthYear", "exn_all_cum_ihs"))
)

stopifnot(is.data.table(qr))
stopifnot(nrow(qr) == 1)

# stage-1 indicator columns
s1_cols <- c("n_treated", "n_matched", "match_rate",
             "pct_tier1", "pct_tier2",
             "avg_controls_per_treated", "min_controls_per_treated")
missing_s1 <- setdiff(s1_cols, names(qr))
if (length(missing_s1) > 0) stop("FAILURE: quality_report missing stage-1 cols: ", paste(missing_s1, collapse=", "))

# stage-2 indicator columns
s2_cols <- c("n_after_cem", "cem_match_rate", "L1_before", "L1_after")
missing_s2 <- setdiff(s2_cols, names(qr))
if (length(missing_s2) > 0) stop("FAILURE: quality_report missing stage-2 cols: ", paste(missing_s2, collapse=", "))

# per-covariate SMD columns
expected_smd <- c("smd_BirthYear_before", "smd_BirthYear_after",
                  "smd_exn_all_cum_ihs_before", "smd_exn_all_cum_ihs_after")
missing_smd <- setdiff(expected_smd, names(qr))
if (length(missing_smd) > 0) stop("FAILURE: quality_report missing per-covariate SMD cols: ", paste(missing_smd, collapse=", "))

# sanity: match_rate in [0,1], tier pct sum to ~1
stopifnot(qr$match_rate >= 0 && qr$match_rate <= 1)
stopifnot(abs((qr$pct_tier1 + qr$pct_tier2) - 1) < 1e-6)
# n_matched should equal uniqueN(treated_id) of matches_info
stopifnot(qr$n_matched == uniqueN(matches_info$treated_id))

print("SUCCESS: quality_report_matches returns expected schema and sane values.")
print(qr)

# 8. Test fast_smd_vr matches bal.tab exactly
print("--- Testing fast_smd_vr vs cobalt::bal.tab ---")
cdf <- refined$data[, .(BirthYear, exn_pm_cum)]
tv <- refined$data$treat
w <- refined$matchit$weights

bt <- bal.tab(cdf, treat = tv, weights = w, un = TRUE,
              stats = c("mean.diffs", "variance.ratios"),
              binary = "std", continuous = "std")
fs <- fast_smd_vr(cdf, tv, w)

for (i in seq_len(nrow(fs))) {
  cn <- fs$cov[i]
  if (abs(fs$smd_un[i] - bt$Balance[cn, "Diff.Un"])        > 1e-6) stop("FAILURE: smd_un mismatch for ", cn)
  if (abs(fs$smd_adj[i] - bt$Balance[cn, "Diff.Adj"])       > 1e-6) stop("FAILURE: smd_adj mismatch for ", cn)
  if (abs(fs$vr_un[i]  - bt$Balance[cn, "V.Ratio.Un"])      > 1e-6) stop("FAILURE: vr_un mismatch for ", cn)
  if (abs(fs$vr_adj[i] - bt$Balance[cn, "V.Ratio.Adj"])      > 1e-6) stop("FAILURE: vr_adj mismatch for ", cn)
}
print("SUCCESS: fast_smd_vr matches bal.tab exactly.")

print("--- Testing run_matching_pipeline ---")
params <- list(
  L             = 3,
  treatment_col = "exn_pm",
  history_col   = "exn_pm_cum",
  unit_id       = "PeopleID",
  time_id       = "year",
  covariates    = c("BirthYear", "exn_all_cum_ihs"),
  cutpoints     = list(BirthYear = 8, exn_all_cum_ihs = 8)
)
qr2 <- run_matching_pipeline(dt, params, sample_n_treated = 500)

stopifnot(is.data.table(qr2) && nrow(qr2) == 1)
# should contain params columns too
stopifnot("L" %in% names(qr2))
stopifnot(qr2$L == 3)
# core indicators should agree with the manually-assembled run
stopifnot(qr2$n_matched == qr$n_matched)
stopifnot(qr2$n_treated == qr$n_treated)

# verify rbind-ability across runs (the whole point of the flat row design)
combined <- rbind(qr, qr2, fill = TRUE)
stopifnot(nrow(combined) == 2)
print("SUCCESS: run_matching_pipeline produces rbind-able quality rows.")
print(qr2)

print("All tests passed.")
