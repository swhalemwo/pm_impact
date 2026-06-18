library(data.table)
library(tinytest)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")
library(cobalt)

# --- Timing collector ---
timings <- data.table(label = character(), time_s = numeric())

# Load small test data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test_smol.csv")
## dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv") # not smol data
n_treated_pop <- uniqueN(dt[exn_pm > 0, PeopleID])
expect_true(n_treated_pop == 1500, info = "expected 1500 PM-treated artists in smol dataset")

# === Skip-slow flag: set FOSTER_TEST_SLOW=false to skip slow cobalt cross-checks ===
run_slow <- Sys.getenv("FOSTER_TEST_SLOW", "true") != "false"

# === 1. Summary mode ===
t0 <- Sys.time()
summary_res <- gd_match_history(dt, L = 3, return_type = "summary")
expect_true(is.data.frame(summary_res), info = "summary should be a data.frame")
expect_true(all(c("total_treated", "treated_matched", "total_matches") %in% names(summary_res)),
            info = "summary has expected columns")
expect_true(summary_res$total_treated > 0, info = "non-zero treated count in summary")
timings <- rbind(timings, data.table(label = "1. Summary mode", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 2. Full match mode ===
t0 <- Sys.time()
full_matches <- gd_match_history(dt, L = 3, return_type = "matches")
expect_true(nrow(full_matches) > 0, info = "matches should be non-empty")
expect_true(all(c("control_id", "treated_id", "year") %in% names(full_matches)),
            info = "matches has base columns")
expect_false(any(c("match_key", "match_tier") %in% names(full_matches)),
             info = "matches mode should NOT contain match_key or match_tier")
timings <- rbind(timings, data.table(label = "2. Full match mode", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 3. Consistency check ===
t0 <- Sys.time()
full_summary <- data.table(
  total_treated     = n_treated_pop,
  treated_matched   = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = mean(full_matches[, .N, by = treated_id]$N),
  total_matches     = nrow(full_matches)
)
expect_equal(as.numeric(summary_res$total_matches),   as.numeric(full_summary$total_matches),
             info = "total_matches agrees between summary and full mode")
expect_equal(as.numeric(summary_res$treated_matched), as.numeric(full_summary$treated_matched),
             info = "treated_matched agrees between summary and full mode")
timings <- rbind(timings, data.table(label = "3. Consistency check", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 4. Matches-info mode ===
t0 <- Sys.time()
matches_info <- gd_match_history(dt, L = 3, return_type = "matches-info")
expect_true(all(c("match_key", "match_tier") %in% names(matches_info)),
            info = "matches-info contains match_key and match_tier")
timings <- rbind(timings, data.table(label = "4. Matches-info mode", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 5. Tier exclusivity ===
t0 <- Sys.time()
dt_tiercnts <- matches_info[, .(nbr_tiers = uniqueN(match_tier)), treated_id]
expect_true(max(dt_tiercnts$nbr_tiers) == 1, info = "all treated units are in exactly one tier")
timings <- rbind(timings, data.table(label = "5. Tier exclusivity", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 6. CEM refinement ===
t0 <- Sys.time()
refined <- gd_refine_matches(matches_info, dt,
                             covariates = c("BirthYear", "exn_all_cum_ihs"),
                             cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))
expect_true(is.list(refined), info = "refine_matches returns a list")
expect_true(all(c("data", "matchit") %in% names(refined)),
            info = "refine_matches list has data and matchit slots")
expect_true("subclass" %in% names(refined$data), info = "refined$data has subclass column")
expect_true(inherits(refined$matchit, "matchit"), info = "refined$matchit is a matchit object")
expect_true(nrow(refined$data) > 0, info = "refined data is non-empty")
timings <- rbind(timings, data.table(label = "6. CEM refinement", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 7. gd_quality_report schema ===
t0 <- Sys.time()
qr <- gd_quality_report(
  dt_prepared  = dt,
  dt_match     = matches_info,
  matchit_obj  = refined$matchit,
  refined_data = refined$data,
  params       = list(L = 3, treatment_col = "exn_pm",
                      covariates = c("BirthYear", "exn_all_cum_ihs")))
expect_true(is.data.table(qr), info = "qr is a data.table")
expect_equal(nrow(qr), 1L, info = "qr has exactly one row")
s1_cols <- c("n_treated", "n_matched", "match_rate", "pct_tier1", "pct_tier2",
             "avg_controls_per_treated", "min_controls_per_treated")
expect_true(all(s1_cols %in% names(qr)), info = "stage-1 columns present")
s2_cols <- c("n_after_cem", "cem_match_rate", "L1_before", "L1_after")
expect_true(all(s2_cols %in% names(qr)), info = "stage-2 columns present")
expected_smd <- c("smd_BirthYear_before", "smd_BirthYear_after",
                  "smd_exn_all_cum_ihs_before", "smd_exn_all_cum_ihs_after")
expect_true(all(expected_smd %in% names(qr)), info = "SMD columns present")
expect_true(qr$match_rate >= 0 && qr$match_rate <= 1, info = "match_rate in [0,1]")
expect_true(abs((qr$pct_tier1 + qr$pct_tier2) - 1) < 1e-6, info = "tier pcts sum to 1")
expect_equal(qr$n_matched, uniqueN(matches_info$treated_id), info = "n_matched equals unique treated_id count")
print(qr)
timings <- rbind(timings, data.table(label = "7. gd_quality_report schema", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 8. gd_smd_vr matches bal.tab exactly ===
t0 <- Sys.time()
cdf <- refined$data[, .(BirthYear, exn_all_cum_ihs)]
tv  <- refined$data$treat
w   <- refined$matchit$weights
bt <- bal.tab(cdf, treat = tv, weights = w, un = TRUE,
              stats = c("mean.diffs", "variance.ratios"),
              binary = "std", continuous = "std")
fs <- gd_smd_vr(cdf, tv, w)
for (i in seq_len(nrow(fs))) {
  cn <- fs$cov[i]
  expect_equal(fs$smd_un[i],  bt$Balance[cn, "Diff.Un"],    tolerance = 1e-6, info = sprintf("%s smd_un", cn))
  expect_equal(fs$smd_adj[i], bt$Balance[cn, "Diff.Adj"],    tolerance = 1e-6, info = sprintf("%s smd_adj", cn))
  expect_equal(fs$vr_un[i],   bt$Balance[cn, "V.Ratio.Un"],  tolerance = 1e-6, info = sprintf("%s vr_un", cn))
  expect_equal(fs$vr_adj[i],  bt$Balance[cn, "V.Ratio.Adj"], tolerance = 1e-6, info = sprintf("%s vr_adj", cn))
}
timings <- rbind(timings, data.table(label = "8. gd_smd_vr vs bal.tab", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 9. gd_l1_med_both vs cobalt::bal.compute (l1.med) [SLOW: ~40s] ===
if (run_slow) {
t0 <- Sys.time()
set.seed(42)
l1_results <- list()
l1_params <- list(
  list(cov = c("BirthYear", "exn_all_cum_ihs"), cp = 8),
  list(cov = c("BirthYear", "exn_all_cum_ihs"), cp = 4),
  list(cov = c("BirthYear"),                     cp = 8)
)
for (p in l1_params) {
  cp <- setNames(rep(list(p$cp), length(p$cov)), p$cov)
  r <- gd_refine_matches(matches_info, dt, covariates = p$cov, cutpoints = cp)
  covs <- p$cov
  cdf <- r$data[, ..covs]
  tv  <- r$data$treat
  w   <- r$matchit$weights

  init <- bal.init(cdf, treat = tv, stat = "l1.med")
  co_before <- bal.compute(init, weights = rep(1, length(tv)))
  co_after  <- bal.compute(init, weights = w)

  my <- gd_l1_med_both(cdf, tv, w, l1.n = 101, seed = 42)

  l1_results[[length(l1_results) + 1]] <- data.table(
    combo = paste(p$cov, collapse = "+"),
    co_before, co_after,
    my_before = my$before, my_after = my$after
  )
}
l1_tab <- rbindlist(l1_results)
l1_tab[, `:=`(diff_before = abs(co_before - my_before),
              diff_after  = abs(co_after - my_after),
              sign_agree_before = sign(co_before) == sign(my_before),
              sign_agree_after  = sign(co_after)  == sign(my_after))]

# Basic sanity checks
expect_true(all(l1_tab$sign_agree_before), info = "gd_l1_med_both L1 before sign agrees with cobalt")
expect_true(all(l1_tab$sign_agree_after),  info = "gd_l1_med_both L1 after sign agrees with cobalt")
# Reduction direction must agree
l1_tab[, `:=`(co_reduction = co_before - co_after,
              my_reduction = my_before - my_after)]
expect_true(all(l1_tab$co_reduction > 0),  info = "cobalt says L1 decreases after matching (all combos)")
expect_true(all(l1_tab$my_reduction > 0),  info = "gd_l1_med_both says L1 decreases after matching (all combos)")

# Mean absolute difference for "before" (same methodology, should be very tight)
mad_before <- mean(l1_tab$diff_before)
expect_true(mad_before < 0.005, info = sprintf("L1 before MAE vs cobalt: %.4f < 0.005", mad_before))

# Rank correlation check for both before and after
if (nrow(l1_tab) >= 3) {
  rank_after <- cor(l1_tab$co_after, l1_tab$my_after, method = "spearman")
  expect_true(rank_after > 0.8, info = sprintf("L1 after Spearman rho vs cobalt: %.2f > 0.8", rank_after))

  rank_before <- cor(l1_tab$co_before, l1_tab$my_before, method = "spearman")
  expect_true(rank_before > 0.8, info = sprintf("L1 before Spearman rho vs cobalt: %.2f > 0.8", rank_before))
}

print(l1_tab[, .(combo, co_before, my_before, diff_before, co_after, my_after, diff_after)])
timings <- rbind(timings, data.table(label = "9. L1 vs cobalt bal.compute", time_s = round(as.numeric(Sys.time() - t0), 2)))
}  # end if (run_slow)

# === 10. gd_run_match_pipeline end-to-end ===
t0 <- Sys.time()
params <- list(L = 3, treatment_col = "exn_pm", history_col = "exn_pm_cum",
               unit_id = "PeopleID", time_id = "year",
               covariates = c("BirthYear", "exn_all_cum_ihs"),
               cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))
qr2 <- gd_run_match_pipeline(dt, params)
expect_true(is.data.table(qr2), info = "qr2 is a data.table")
expect_equal(nrow(qr2), 1L, info = "qr2 has exactly one row")
expect_true("L" %in% names(qr2), info = "qr2 includes param columns")
expect_equal(qr2$L, 3L, info = "L param passed through")
expect_equal(qr2$n_matched, qr$n_matched, info = "n_matched consistent")
expect_equal(qr2$n_treated, qr$n_treated, info = "n_treated consistent")
combined <- rbind(qr, qr2, fill = TRUE)
expect_equal(nrow(combined), 2L, info = "rbind-able: combined has 2 rows")
print(qr2)
timings <- rbind(timings, data.table(label = "10. Pipeline e2e", time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 11. gd_l1_med wrapper (before vs after) ===
t0 <- Sys.time()
l1_before <- gd_l1_med(cdf, tv, rep(1, length(tv)), l1.n = 51)
l1_after  <- gd_l1_med(cdf, tv, w, l1.n = 51)
expect_true(is.numeric(l1_before) && l1_before >= 0 && l1_before <= 1,
            info = "gd_l1_med before in [0,1]")
expect_true(is.numeric(l1_after)  && l1_after  >= 0 && l1_after  <= 1,
            info = "gd_l1_med after in [0,1]")
expect_true(l1_before > l1_after, info = "gd_l1_med before > after (balance improves)")
timings <- rbind(timings, data.table(label = "11. gd_l1_med wrapper",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 12. gd_run_match_pipeline — empty matches ===
t0 <- Sys.time()
qr_empty <- gd_run_match_pipeline(dt[FALSE], params)
expect_true(is.data.table(qr_empty), info = "qr_empty is a data.table")
expect_equal(nrow(qr_empty), 1L, info = "qr_empty has one row")
expect_equal(qr_empty$n_matched, 0L, info = "qr_empty n_matched == 0")
expect_true(is.na(qr_empty$L1_before), info = "qr_empty L1 is NA")
timings <- rbind(timings, data.table(label = "12. Pipeline empty matches",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 13. gd_run_match_pipeline — with sample_n_treated ===
t0 <- Sys.time()
qr_smp <- gd_run_match_pipeline(dt, params, sample_n_treated = 10L)
expect_true(is.data.table(qr_smp), info = "qr_smp is a data.table")
expect_equal(nrow(qr_smp), 1L, info = "qr_smp has one row")
expect_equal(qr_smp$n_matched, qr2$n_matched,
             info = "sample_n_treated leaves full match count unchanged")
expect_true(qr_smp$n_after_cem <= 10L,
            info = "sample_n_treated limits CEM stage to ≤10 treated")
timings <- rbind(timings, data.table(label = "13. Pipeline with sampling",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 14. gd_quality_report_s1 — empty dt_match ===
t0 <- Sys.time()
s1_empty <- gd_quality_report_s1(dt, data.table(), params)
expect_true(is.data.table(s1_empty), info = "s1_empty is a data.table")
expect_equal(s1_empty$n_matched, 0L, info = "s1_empty n_matched == 0")
expect_equal(s1_empty$match_rate, 0, info = "s1_empty match_rate == 0")
expect_true(is.na(s1_empty$pct_tier1), info = "s1_empty pct_tier1 is NA")
timings <- rbind(timings, data.table(label = "14. qr_s1 empty match",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 15. gd_quality_report_s2 — NULL matchit + no-cov intersect ===
t0 <- Sys.time()
s2_null <- gd_quality_report_s2(NULL, NULL, params)
expect_true(is.data.table(s2_null), info = "s2_null is a data.table")
expect_true(all(is.na(s2_null$n_after_cem)), info = "s2_null n_after_cem is NA")

s2_nocov <- gd_quality_report_s2(refined$matchit, refined$data,
                                  list(covariates = "nonexistent"))
expect_true(is.data.table(s2_nocov), info = "s2_nocov is a data.table")
expect_true(all(is.na(s2_nocov$L1_before)), info = "s2_nocov L1_before is NA")
timings <- rbind(timings, data.table(label = "15. qr_s2 NULL/no-cov",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 16. gd_bal_cols_cov — no matching cov in smd_results ===
t0 <- Sys.time()
dt_smd <- gd_smd_vr(cdf, tv, w)
bc <- gd_bal_cols_cov("nonexistent", dt_smd)
expect_true(is.data.table(bc), info = "bc is a data.table")
expect_equal(nrow(bc), 1L, info = "bc has one row")
expect_true(all(is.na(unlist(bc))), info = "bc all NA for missing cov")
timings <- rbind(timings, data.table(label = "16. bal_cols no match",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 17. gd_ematch_yr — NULL return (empty year) ===
t0 <- Sys.time()
dt_dummy <- copy(dt)[, `:=`(tid = year, treat_bin = fifelse(exn_pm > 0, 1L, 0L))]
null_res <- gd_ematch_yr(1800, dt_dummy, NULL, "exn_pm", "matches", 1)
expect_true(is.null(null_res), info = "ematch_yr returns NULL for empty year")
timings <- rbind(timings, data.table(label = "17. ematch_yr NULL",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# --- Timing report ---
cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
