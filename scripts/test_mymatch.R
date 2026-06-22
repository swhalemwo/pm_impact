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
expect_true(nrow(matches_info) > 0, info = "matches-info should be non-empty")
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

# === 18. coarse_bin — edge cases ===
t0 <- Sys.time()
expect_equal(coarse_bin(c(1,2,3), 1), c(1L, 1L, 1L), info = "coarse_bin nbins=1 gives all 1")
expect_equal(coarse_bin(rep(5, 10), 8), rep(1L, 10), info = "coarse_bin zero-range gives all 1")
expect_true(all(coarse_bin(1:100, 4) %in% 1:4), info = "coarse_bin 100 values into 4 bins")
expect_true(max(coarse_bin(rnorm(100), 10)) <= 10, info = "coarse_bin max bin ≤ nbins")
expect_true(min(coarse_bin(rnorm(100), 10)) >= 1,  info = "coarse_bin min bin ≥ 1")
timings <- rbind(timings, data.table(label = "18. coarse_bin edges",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 19. gd_refine_matches — subclass and weight invariants ===
t0 <- Sys.time()
## reuse `refined` from section 6
n_s <- refined$data[!is.na(subclass), .N]
n_u <- refined$data[is.na(subclass), .N]
s_ids <- unique(na.omit(refined$data$subclass))

## no duplicate subclass IDs
expect_equal(length(s_ids), length(unique(s_ids)),
             info = "all subclass IDs are unique")
## every matched cell has both treated and control
sub_comp <- refined$data[!is.na(subclass),
                         .(n_treat = sum(treat == 1), n_ctrl = sum(treat == 0)),
                         by = subclass]
expect_true(all(sub_comp$n_treat >= 1),
            info = "every subclass has ≥1 treated")
expect_true(all(sub_comp$n_ctrl >= 1),
            info = "every subclass has ≥1 control")
## weights: 1 for matched, 0 for unmatched
expect_equal(refined$data[!is.na(subclass), range(weight)], c(1, 1),
             info = "all matched weights are 1")
expect_equal(refined$data[is.na(subclass), range(weight)], c(0, 0),
             info = "all unmatched weights are 0")
## matchit weights match
expect_equal(refined$matchit$weights, refined$data$weight,
             info = "matchit$weights = data$weight")
expect_equal(as.integer(refined$matchit$subclass),
             as.integer(refined$data$subclass),
             info = "matchit$subclass matches data$subclass")
## class set
expect_true(inherits(refined$matchit, "matchit"),
            info = "matchit object has class matchit")
timings <- rbind(timings, data.table(label = "19. refine invariants",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 20. gd_refine_matches — single covariate + single bin ===
t0 <- Sys.time()
r_1c <- gd_refine_matches(matches_info, dt,
                           covariates = "BirthYear",
                           cutpoints = list(BirthYear = 8))
expect_true(is.list(r_1c) && all(c("data", "matchit") %in% names(r_1c)),
            info = "single covariate: list structure ok")
expect_true("subclass" %in% names(r_1c$data),
            info = "single covariate: subclass column present")
expect_true(nrow(r_1c$data) == nrow(refined$data),
            info = "single covariate: same row count")

r_1b <- gd_refine_matches(matches_info, dt,
                           covariates = c("BirthYear", "exn_all_cum_ihs"),
                           cutpoints = list(BirthYear = 1, exn_all_cum_ihs = 1))
expect_true(nrow(r_1b$data) == nrow(refined$data),
            info = "single bin: same row count")
## with 1 bin each, coarse covariates collapse to one value;
## subclasses then equal the number of distinct sets (= matched treated)
expect_equal(uniqueN(na.omit(r_1b$data$subclass)),
             uniqueN(r_1b$data[!is.na(subclass) & treat == 1, set]),
             info = "single bin per covariate → one subclass per set")
timings <- rbind(timings, data.table(label = "20. refine single cov/bin",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 21. gd_refine_matches_matchit — structure check ===
t0 <- Sys.time()
r_mi <- gd_refine_matches_matchit(matches_info, dt,
                                   covariates = c("BirthYear", "exn_all_cum_ihs"),
                                   cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8))
expect_true(is.list(r_mi), info = "matchit version returns a list")
expect_true(all(c("data", "matchit") %in% names(r_mi)),
            info = "matchit version has data and matchit slots")
expect_true(inherits(r_mi$matchit, "matchit"),
            info = "matchit version matchit is a matchit object")
expect_equal(nrow(r_mi$data), nrow(refined$data),
             info = "matchit version same row count as DT version")
timings <- rbind(timings, data.table(label = "21. refine_matchit struct",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 22. DT vs matchit CEM — equal-width bins [SLOW] ===
if (run_slow) {
t0 <- Sys.time()
r_dt <- gd_refine_matches(matches_info, dt,
                           covariates = c("BirthYear", "exn_all_cum_ihs"),
                           cutpoints = list(BirthYear = 8L, exn_all_cum_ihs = 8L))
r_mi <- gd_refine_matches_matchit(matches_info, dt,
                                   covariates = c("BirthYear", "exn_all_cum_ihs"),
                                   cutpoints = list(BirthYear = 8L, exn_all_cum_ihs = 8L))
expect_true(gd_match_equivalent(r_dt, r_mi, dt,
                                 c("BirthYear", "exn_all_cum_ihs"),
                                 list(BirthYear = 8L, exn_all_cum_ihs = 8L)),
            info = "equal-width bins: DT and matchit produce equivalent results")
timings <- rbind(timings, data.table(label = "22. DT vs matchit equal-width",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))
}  # end if (run_slow)

# === 23. gd_coarse_breaks — structure and edge cases ===
t0 <- Sys.time()
br8 <- gd_coarse_breaks(dt, c("BirthYear", "exn_all_cum_ihs"),
                         list(BirthYear = 8, exn_all_cum_ihs = 8))
expect_true(is.list(br8), info = "coarse_breaks returns a list")
expect_equal(names(br8), c("BirthYear", "exn_all_cum_ihs"),
             info = "coarse_breaks: names match covariates")
expect_equal(length(br8$BirthYear), 9,
             info = "coarse_breaks: 8 bins → 9 breakpoints")
expect_equal(length(br8$exn_all_cum_ihs), 9,
             info = "coarse_breaks: 8 bins → 9 breakpoints")
expect_true(br8$BirthYear[1] <= min(dt$BirthYear),
            info = "coarse_breaks: first break ≤ min")
expect_true(br8$BirthYear[9] >= max(dt$BirthYear),
            info = "coarse_breaks: last break ≥ max")
expect_true(all(diff(br8$BirthYear) > 0),
            info = "coarse_breaks: strictly increasing")

## edge: 1 bin
br1 <- gd_coarse_breaks(dt, "BirthYear", list(BirthYear = 1))
expect_equal(length(br1$BirthYear), 2,
             info = "coarse_breaks: 1 bin → 2 breakpoints (min, max)")

## edge: fewer unique values than bins
dt_small <- dt[1:10]
br_uniq <- gd_coarse_breaks(dt_small, "BirthYear", list(BirthYear = 20))
expect_true(length(br_uniq$BirthYear) <= 10,
            info = "coarse_breaks: ≤ n unique when fewer than nbins")
timings <- rbind(timings, data.table(label = "23. coarse_breaks structure",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 24. gd_refine_matches with quantile breaks — internal consistency ===
t0 <- Sys.time()
l_breaks <- gd_coarse_breaks(dt, c("BirthYear", "exn_all_cum_ihs"),
                              list(BirthYear = 8, exn_all_cum_ihs = 8))
r_q <- gd_refine_matches(matches_info, dt,
                          covariates = c("BirthYear", "exn_all_cum_ihs"),
                          cutpoints = l_breaks)
expect_true(is.list(r_q), info = "quantile refine: returns a list")
expect_true("subclass" %in% names(r_q$data),
            info = "quantile refine: subclass column present")
expect_true(inherits(r_q$matchit, "matchit"),
            info = "quantile refine: matchit object class")
## check each matched cell has ≥1 treated and ≥1 control
sub_comp_q <- r_q$data[!is.na(subclass),
                        .(n_treat = sum(treat == 1), n_ctrl = sum(treat == 0)),
                        by = subclass]
expect_true(all(sub_comp_q$n_treat >= 1),
            info = "quantile refine: every subclass has ≥1 treated")
expect_true(all(sub_comp_q$n_ctrl >= 1),
            info = "quantile refine: every subclass has ≥1 control")
## weights correct
expect_true(all(r_q$data[!is.na(subclass), weight == 1]),
            info = "quantile refine: all matched have weight 1")
expect_true(all(r_q$data[is.na(subclass), weight == 0]),
            info = "quantile refine: all unmatched have weight 0")
timings <- rbind(timings, data.table(label = "24. refine quantile breaks",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 25. gd_coarse_breaks + refine — fingerprint reproducibility ===
t0 <- Sys.time()
## Running twice with same breaks reproduces identical subclasses
l_breaks2 <- gd_coarse_breaks(dt, c("BirthYear", "exn_all_cum_ihs"),
                               list(BirthYear = 8, exn_all_cum_ihs = 8))
r_q2 <- gd_refine_matches(matches_info, dt,
                           covariates = c("BirthYear", "exn_all_cum_ihs"),
                           cutpoints = l_breaks2)
fp1 <- r_q$data[!is.na(subclass), .(PeopleID, year, treat, s = subclass)][order(PeopleID, year, treat)]
fp2 <- r_q2$data[!is.na(subclass), .(PeopleID, year, treat, s = subclass)][order(PeopleID, year, treat)]
expect_equal(nrow(fp1), nrow(fp2),
             info = "quantile refine: same number of matched rows on re-run")
expect_true(all(fp1$s == fp2$s),
            info = "quantile refine: subclass assignment reproducible")
timings <- rbind(timings, data.table(label = "25. refine quantile reproducible",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))

# === 26. quantile breaks: DT vs matchit CEM [SLOW] ===
if (run_slow) {
t0 <- Sys.time()
l_breaks2 <- gd_coarse_breaks(dt, c("BirthYear", "exn_all_cum_ihs"),
                               list(BirthYear = 8, exn_all_cum_ihs = 8))
r_dt2 <- gd_refine_matches(matches_info, dt,
                            covariates = c("BirthYear", "exn_all_cum_ihs"),
                            cutpoints = l_breaks2)
r_mi2 <- gd_refine_matches_matchit(matches_info, dt,
                                    covariates = c("BirthYear", "exn_all_cum_ihs"),
                                    cutpoints = l_breaks2)
expect_true(gd_match_equivalent(r_dt2, r_mi2, dt,
                                 c("BirthYear", "exn_all_cum_ihs"), l_breaks2),
            info = "quantile breaks: DT and matchit produce equivalent results")
timings <- rbind(timings, data.table(label = "26. quantile DT vs MI",
                                     time_s = round(as.numeric(Sys.time() - t0), 2)))
}  # end if (run_slow)

# --- Timing report ---
cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
