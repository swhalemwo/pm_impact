library(data.table)
library(tinytest)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Load small test data (3,000 artists, 58K rows)
dt <- fread("~/Dropbox/phd/papers/foster/dt_test_smol.csv")
n_treated_pop <- uniqueN(dt[exn_pm > 0, PeopleID])
expect_true(n_treated_pop == 1500, info = "expected 1500 PM-treated artists in smol dataset")

# === 1. Summary mode ===
summary_res <- match_treatment_history(dt, L = 3, return_type = "summary")
expect_true(is.data.frame(summary_res), info = "summary should be a data.frame")
expect_true(all(c("total_treated", "treated_matched", "total_matches") %in% names(summary_res)),
            info = "summary has expected columns")
expect_true(summary_res$total_treated > 0, info = "non-zero treated count in summary")

# === 2. Full match mode (no debugging columns) ===
full_matches <- match_treatment_history(dt, L = 3, return_type = "matches")
expect_true(nrow(full_matches) > 0, info = "matches should be non-empty")
expect_true(all(c("control_id", "treated_id", "year") %in% names(full_matches)),
            info = "matches has base columns")
expect_false(any(c("match_key", "match_tier") %in% names(full_matches)),
             info = "matches mode should NOT contain match_key or match_tier")

# === 3. Consistency check: summary mode == replicated from full matches ===
full_summary <- data.table(
  total_treated     = n_treated_pop,
  treated_matched   = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = mean(full_matches[, .N, by = treated_id]$N),
  total_matches     = nrow(full_matches)
)
expect_equal(as.numeric(summary_res$total_matches),    as.numeric(full_summary$total_matches),
             info = "total_matches agrees between summary and full mode")
expect_equal(as.numeric(summary_res$treated_matched),  as.numeric(full_summary$treated_matched),
             info = "treated_matched agrees between summary and full mode")

# === 4. Matches-info mode (includes debugging columns) ===
matches_info <- match_treatment_history(dt, L = 3, return_type = "matches-info")
expect_true(all(c("match_key", "match_tier") %in% names(matches_info)),
            info = "matches-info contains match_key and match_tier")

# === 5. Tier exclusivity ===
dt_tiercnts <- matches_info[, .(nbr_tiers = uniqueN(match_tier)), treated_id]
expect_true(max(dt_tiercnts$nbr_tiers) == 1, info = "all treated units are in exactly one tier")

# === 6. Refinement (CEM) ===
refined <- refine_matches(matches_info, dt,
                          covariates = c("BirthYear", "exn_all_cum_ihs"),
                          cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))

expect_true(is.list(refined), info = "refine_matches returns a list")
expect_true(all(c("data", "matchit") %in% names(refined)),
            info = "refine_matches list has data and matchit slots")
expect_true("subclass" %in% names(refined$data),
            info = "refined$data has subclass column")
expect_true(inherits(refined$matchit, "matchit"),
            info = "refined$matchit is a matchit object")
expect_true(nrow(refined$data) > 0, info = "refined data is non-empty")

# === 7. quality_report_matches schema ===
qr <- quality_report_matches(
  dt_prepared  = dt,
  dt_match     = matches_info,
  matchit_obj  = refined$matchit,
  refined_data = refined$data,
  params       = list(L = 3, treatment_col = "exn_pm",
                      covariates = c("BirthYear", "exn_all_cum_ihs"))
)

expect_true(is.data.table(qr), info = "qr is a data.table")
expect_equal(nrow(qr), 1L, info = "qr has exactly one row")

s1_cols <- c("n_treated", "n_matched", "match_rate",
             "pct_tier1", "pct_tier2",
             "avg_controls_per_treated", "min_controls_per_treated")
expect_true(all(s1_cols %in% names(qr)),
            info = sprintf("stage-1 columns present; missing: %s",
                           paste(setdiff(s1_cols, names(qr)), collapse = ", ")))

s2_cols <- c("n_after_cem", "cem_match_rate", "L1_before", "L1_after")
expect_true(all(s2_cols %in% names(qr)),
            info = sprintf("stage-2 columns present; missing: %s",
                           paste(setdiff(s2_cols, names(qr)), collapse = ", ")))

expected_smd <- c("smd_BirthYear_before", "smd_BirthYear_after",
                  "smd_exn_all_cum_ihs_before", "smd_exn_all_cum_ihs_after")
expect_true(all(expected_smd %in% names(qr)),
            info = sprintf("per-covariate SMD columns present; missing: %s",
                           paste(setdiff(expected_smd, names(qr)), collapse = ", ")))

expect_true(qr$match_rate >= 0, info = "match_rate >= 0")
expect_true(qr$match_rate <= 1, info = "match_rate <= 1")
expect_true(abs((qr$pct_tier1 + qr$pct_tier2) - 1) < 1e-6,
            info = "pct_tier1 + pct_tier2 sums to 1")
expect_equal(qr$n_matched, uniqueN(matches_info$treated_id),
             info = "n_matched equals unique treated_id count")

print(qr)

# === 8. fast_smd_vr matches cobalt::bal.tab exactly ===
cdf <- refined$data[, .(BirthYear, exn_all_cum_ihs)]
tv  <- refined$data$treat
w   <- refined$matchit$weights

bt <- bal.tab(cdf, treat = tv, weights = w, un = TRUE,
              stats = c("mean.diffs", "variance.ratios"),
              binary = "std", continuous = "std")
fs <- fast_smd_vr(cdf, tv, w)

for (i in seq_len(nrow(fs))) {
  cn <- fs$cov[i]
  expect_equal(fs$smd_un[i],   bt$Balance[cn, "Diff.Un"],    tolerance = 1e-6,
               info = sprintf("%s smd_un", cn))
  expect_equal(fs$smd_adj[i],  bt$Balance[cn, "Diff.Adj"],    tolerance = 1e-6,
               info = sprintf("%s smd_adj", cn))
  expect_equal(fs$vr_un[i],    bt$Balance[cn, "V.Ratio.Un"],  tolerance = 1e-6,
               info = sprintf("%s vr_un", cn))
  expect_equal(fs$vr_adj[i],   bt$Balance[cn, "V.Ratio.Adj"], tolerance = 1e-6,
               info = sprintf("%s vr_adj", cn))
}

# === 9. run_matching_pipeline end-to-end ===
params <- list(
  L             = 3,
  treatment_col = "exn_pm",
  history_col   = "exn_pm_cum",
  unit_id       = "PeopleID",
  time_id       = "year",
  covariates    = c("BirthYear", "exn_all_cum_ihs"),
  cutpoints     = list(BirthYear = 8, exn_all_cum_ihs = 8)
)
qr2 <- run_matching_pipeline(dt, params)

expect_true(is.data.table(qr2), info = "qr2 is a data.table")
expect_equal(nrow(qr2), 1L, info = "qr2 has exactly one row")
expect_true("L" %in% names(qr2), info = "qr2 includes param columns")
expect_equal(qr2$L, 3L, info = "L param passed through")
expect_equal(qr2$n_matched, qr$n_matched, info = "n_matched consistent")
expect_equal(qr2$n_treated, qr$n_treated, info = "n_treated consistent")

combined <- rbind(qr, qr2, fill = TRUE)
expect_equal(nrow(combined), 2L, info = "rbind-able: combined has 2 rows")

print(qr2)
