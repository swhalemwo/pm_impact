library(data.table)
library(tinytest)
library(collapse, include.only = c("join", "funique", "set_collapse"))
set_collapse(verbose = 0)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")
library(cobalt)

# --- Timing collector ---
timings <- data.table(label = character(), time_s = numeric())

# Load small test data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test_smol.csv")
dt_med <- fread("~/Dropbox/phd/papers/foster/dt_test_med.csv")
dt_big <- fread("~/Dropbox/phd/papers/foster/dt_test_big.csv")



# === Config definitions ===
configs <- list(
  list(dt = dt, L = 3L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
       uid = "PeopleID", tid = "year",
       covariates = c("BirthYear", "exn_all_cum_ihs"),
       cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       slow = TRUE, label = "smol_L3")
  ## list(dt = dt, L = 5L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
  ##      uid = "PeopleID", tid = "year",
  ##      covariates = c("BirthYear", "exn_all_cum_ihs"),
  ##      cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
  ##      slow = TRUE, label = "smol_L5"),
  ## list(dt = dt_med, L = 5L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
  ##      uid = "PeopleID", tid = "year",
  ##      covariates = c("BirthYear", "exn_all_cum_ihs"),
  ##      cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
  ##      slow = FALSE, label = "med_L5"),
  ## list(dt = dt_big, L = 5L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
  ##      uid = "PeopleID", tid = "year",
  ##      covariates = c("BirthYear", "exn_all_cum_ihs"),
  ##      cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
  ##      slow = FALSE, label = "big_L5")
)

# ==== Test functions ====

gd_test_basic <- function(cfg) {
  ## sections 1-5: summary, matches, consistency, matches-info, tier exclusivity
  t0 <- Sys.time()
  r_s <- gd_match_history(cfg$dt, L = cfg$L, treatment_col = cfg$treatment_col,  ## r_s: match history summary
                           history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid,
                           return_type = "summary")
  expect_true(is.data.frame(r_s), info = sprintf("[%s] summary is data.frame", cfg$label))
  expect_true(all(c("total_treated", "treated_matched", "total_matches") %in% names(r_s)),
              info = sprintf("[%s] summary columns", cfg$label))
  expect_true(r_s$total_treated > 0, info = sprintf("[%s] non-zero treated", cfg$label))

  r_f <- gd_match_history(cfg$dt, L = cfg$L, treatment_col = cfg$treatment_col,  ## r_f: match history full matches
                           history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid,
                           return_type = "matches")
  expect_true(nrow(r_f) > 0, info = sprintf("[%s] matches non-empty", cfg$label))
  expect_true(all(c("control_id", "treated_id", "year") %in% names(r_f)),
              info = sprintf("[%s] matches columns", cfg$label))
  expect_false(any(c("match_key", "match_tier") %in% names(r_f)),
               info = sprintf("[%s] matches no extra cols", cfg$label))

  ## consistency: summary matches full data
  expect_equal(as.numeric(r_s$total_matches), as.numeric(nrow(r_f)),
               info = sprintf("[%s] total_matches consistent", cfg$label))
  expect_equal(as.numeric(r_s$treated_matched), as.numeric(uniqueN(r_f$treated_id)),
               info = sprintf("[%s] treated_matched consistent", cfg$label))

  r_mi <- gd_match_history(cfg$dt, L = cfg$L, treatment_col = cfg$treatment_col,  ## r_mi: match history matches-info
                            history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid,
                            return_type = "matches-info")
  expect_true(nrow(r_mi) > 0, info = sprintf("[%s] matches-info non-empty", cfg$label))
  expect_true(all(c("match_key", "match_tier") %in% names(r_mi)),
              info = sprintf("[%s] matches-info columns", cfg$label))

  ## check tier exclusivity: each treated unit in exactly one tier
  dt_tc <- r_mi[, .(nbr = uniqueN(match_tier)), treated_id]  ## dt_tc: tier counts per treated
  expect_true(max(dt_tc$nbr) == 1, info = sprintf("[%s] one tier per treated", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s basic", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  r_mi
}

gd_test_cem_refine <- function(cfg, mi) {
  ## section 6: CEM refinement structure
  t0 <- Sys.time()
  r <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                          covariates = cfg$covariates, cutpoints = cfg$cutpoints)
  expect_true(is.list(r), info = sprintf("[%s] refine is list", cfg$label))
  expect_true(all(c("data", "matchit") %in% names(r)),
              info = sprintf("[%s] refine slots", cfg$label))
  expect_true("subclass" %in% names(r$data),
              info = sprintf("[%s] subclass column", cfg$label))
  expect_true(inherits(r$matchit, "matchit"),
              info = sprintf("[%s] matchit class", cfg$label))
  expect_true(nrow(r$data) > 0, info = sprintf("[%s] refine non-empty", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s CEM refine", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  r
}

gd_test_quality_report <- function(cfg, mi, refined) {
  ## section 7: quality report schema
  t0 <- Sys.time()
  params <- list(L = cfg$L, treatment_col = cfg$treatment_col, uid = cfg$uid,
                 covariates = cfg$covariates)
  dt_qr <- gd_quality_report(  ## dt_qr: quality report (stage 1 + stage 2 + params)
    dt_prepared  = cfg$dt,
    dt_match     = mi,
    matchit_obj  = refined$matchit,
    refined_data = refined$data,
    params       = params)

  expect_true(is.data.table(dt_qr), info = sprintf("[%s] qr is data.table", cfg$label))
  expect_equal(nrow(dt_qr), 1L, info = sprintf("[%s] qr one row", cfg$label))
  s1_cols <- c("n_treated", "n_matched", "match_rate", "pct_tier1", "pct_tier2",
               "avg_controls_per_treated", "min_controls_per_treated")
  expect_true(all(s1_cols %in% names(dt_qr)), info = sprintf("[%s] stage-1 cols", cfg$label))
  s2_cols <- c("n_after_cem", "cem_match_rate", "L1_before", "L1_after")
  expect_true(all(s2_cols %in% names(dt_qr)), info = sprintf("[%s] stage-2 cols", cfg$label))
  invisible(lapply(cfg$covariates, function(v) {
    expect_true(paste0("smd_", v, "_before") %in% names(dt_qr),
                info = sprintf("[%s] SMD %s before", cfg$label, v))
    expect_true(paste0("smd_", v, "_after") %in% names(dt_qr),
                info = sprintf("[%s] SMD %s after", cfg$label, v))
  }))
  expect_true(dt_qr$match_rate >= 0 && dt_qr$match_rate <= 1,
              info = sprintf("[%s] match_rate in [0,1]", cfg$label))
  expect_true(abs((dt_qr$pct_tier1 + dt_qr$pct_tier2) - 1) < 1e-6,
              info = sprintf("[%s] tier pcts sum to 1", cfg$label))
  expect_equal(dt_qr$n_matched, uniqueN(mi$treated_id),
               info = sprintf("[%s] n_matched correct", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s quality report", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  dt_qr
}

gd_test_smd_vr <- function(cfg, refined) {
  t0 <- Sys.time()
  ## compare gd_smd_vr against cobalt::bal.tab
  dt_cdf <- refined$data[, cfg$covariates, with = FALSE]  ## dt_cdf: covariate data frame
  tv  <- refined$data$treat                                ## tv: treatment vector (0/1)
  w   <- refined$matchit$weights                            ## w: CEM weights

  dt_bt <- bal.tab(dt_cdf, treat = tv, weights = w, un = TRUE,  ## dt_bt: bal.tab output from cobalt
                stats = c("mean.diffs", "variance.ratios"),
                binary = "std", continuous = "std")
  dt_smd <- gd_smd_vr(dt_cdf, tv, w)  ## dt_smd: standardized mean differences + variance ratios
  invisible(lapply(seq_len(nrow(dt_smd)), function(i) {
    cn <- dt_smd$cov[i]
    expect_equal(dt_smd$smd_un[i],  dt_bt$Balance[cn, "Diff.Un"],    tolerance = 1e-6,
                 info = sprintf("[%s] %s smd_un", cfg$label, cn))
    expect_equal(dt_smd$smd_adj[i], dt_bt$Balance[cn, "Diff.Adj"],    tolerance = 1e-6,
                 info = sprintf("[%s] %s smd_adj", cfg$label, cn))
    expect_equal(dt_smd$vr_un[i],   dt_bt$Balance[cn, "V.Ratio.Un"],  tolerance = 1e-6,
                 info = sprintf("[%s] %s vr_un", cfg$label, cn))
    expect_equal(dt_smd$vr_adj[i],  dt_bt$Balance[cn, "V.Ratio.Adj"], tolerance = 1e-6,
                 info = sprintf("[%s] %s vr_adj", cfg$label, cn))
  }))

  ## gd_l1_med wrapper
  l1_before <- gd_l1_med(dt_cdf, tv, rep(1, length(tv)), l1.n = 51)
  l1_after  <- gd_l1_med(dt_cdf, tv, w, l1.n = 51)
  expect_true(is.numeric(l1_before) && l1_before >= 0 && l1_before <= 1,
              info = sprintf("[%s] l1_med before in [0,1]", cfg$label))
  expect_true(is.numeric(l1_after)  && l1_after  >= 0 && l1_after  <= 1,
              info = sprintf("[%s] l1_med after in [0,1]", cfg$label))

  ## gd_bal_cols_cov edge case: no matching cov in smd_results
  dt_bc <- gd_bal_cols_cov("nonexistent", dt_smd)  ## dt_bc: balance columns (smd/vr) for one cov
  expect_true(is.data.table(dt_bc), info = sprintf("[%s] bc data.table", cfg$label))
  expect_equal(nrow(dt_bc), 1L, info = sprintf("[%s] bc one row", cfg$label))
  expect_true(all(is.na(unlist(dt_bc))), info = sprintf("[%s] bc all NA", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s smd_vr", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  list(cdf = dt_cdf, tv = tv, w = w)
}

gd_test_l1_crossval <- function(cfg, mi) {
  ## section 9: L1 cross-validation against cobalt [SLOW]
  if (!run_slow || !cfg$slow) return()
  t0 <- Sys.time()
  set.seed(42)
  l_params <- list(
    list(cov = cfg$covariates,                     cp = 8),
    list(cov = cfg$covariates,                     cp = 6),
    list(cov = cfg$covariates,                     cp = 4),
    list(cov = cfg$covariates,                     cp = 2),
    list(cov = cfg$covariates[1],                   cp = 8),
    list(cov = cfg$covariates[1],                   cp = 4)
  )
  dt_l1 <- rbindlist(lapply(l_params, function(p) {  ## dt_l1: L1 imbalance results across parameter combos
    cp <- setNames(rep(list(p$cp), length(p$cov)), p$cov)
    if (length(cp) == 0) return(NULL)
    r <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                            covariates = p$cov, cutpoints = cp)
    covs <- p$cov
    dt_cov <- r$data[, covs, with = FALSE]
    tv  <- r$data$treat
    w   <- r$matchit$weights

    init <- bal.init(dt_cov, treat = tv, stat = "l1.med")
    co_before <- bal.compute(init, weights = rep(1, length(tv)))
    co_after  <- bal.compute(init, weights = w)
    my <- gd_l1_med_both(dt_cov, tv, w, l1.n = 101, seed = 42)

    data.table(combo = paste(p$cov, collapse = "+"),
               co_before, co_after, my_before = my$before, my_after = my$after)
  }))
  dt_l1[, `:=`(diff_before = abs(co_before - my_before),
                diff_after  = abs(co_after - my_after),
                sign_agree_before = sign(co_before) == sign(my_before),
                sign_agree_after  = sign(co_after) == sign(my_after))]

  expect_true(all(dt_l1$sign_agree_before),
              info = sprintf("[%s] L1 before sign agrees", cfg$label))
  expect_true(all(dt_l1$sign_agree_after),
              info = sprintf("[%s] L1 after sign agrees", cfg$label))
  dt_l1[, `:=`(co_reduction = co_before - co_after,
                my_reduction = my_before - my_after)]
  expect_true(all(dt_l1$co_reduction > 0),
              info = sprintf("[%s] cobalt L1 decreases", cfg$label))
  expect_true(all(dt_l1$my_reduction > 0),
              info = sprintf("[%s] our L1 decreases", cfg$label))
  mad_before <- mean(dt_l1$diff_before)
  expect_true(mad_before < 0.005,
              info = sprintf("[%s] L1 before MAE < 0.005: %.4f", cfg$label, mad_before))
  if (nrow(dt_l1) >= 3) {
    r_after <- cor(dt_l1$co_after, dt_l1$my_after, method = "spearman")
    expect_true(r_after > 0.8,
                info = sprintf("[%s] L1 after rho > 0.8: %.2f", cfg$label, r_after))
  }

  timings <<- rbind(timings, data.table(label = sprintf("%s L1 crossval", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_pipeline <- function(cfg, mi, refined, dt_qr) {
  ## sections 10-13: pipeline e2e, empty matches, sampling
  t0 <- Sys.time()
  params <- list(L = cfg$L, treatment_col = cfg$treatment_col,
                 history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid,
                 covariates = cfg$covariates, cutpoints = cfg$cutpoints)

  dt_qr2 <- gd_run_match_pipeline(cfg$dt, params)  ## dt_qr2: quality report from pipeline run
  expect_true(is.data.table(dt_qr2), info = sprintf("[%s] qr2 data.table", cfg$label))
  expect_equal(nrow(dt_qr2), 1L, info = sprintf("[%s] qr2 one row", cfg$label))
  expect_true("L" %in% names(dt_qr2), info = sprintf("[%s] qr2 has L", cfg$label))
  expect_equal(dt_qr2$L, cfg$L, info = sprintf("[%s] L passed through", cfg$label))
  expect_equal(dt_qr2$n_matched, dt_qr$n_matched,
               info = sprintf("[%s] n_matched consistent", cfg$label))
  expect_equal(dt_qr2$n_treated, dt_qr$n_treated,
               info = sprintf("[%s] n_treated consistent", cfg$label))
  combined <- rbind(dt_qr, dt_qr2, fill = TRUE)
  expect_equal(nrow(combined), 2L, info = sprintf("[%s] rbind-able", cfg$label))

  ## empty matches test: no data → NA results
  dt_qr_empty <- gd_run_match_pipeline(cfg$dt[FALSE], params)  ## dt_qr_empty: quality report, empty input
  expect_true(is.data.table(dt_qr_empty), info = sprintf("[%s] empty data.table", cfg$label))
  expect_equal(nrow(dt_qr_empty), 1L, info = sprintf("[%s] empty one row", cfg$label))
  expect_equal(dt_qr_empty$n_matched, 0L, info = sprintf("[%s] empty n_matched=0", cfg$label))
  expect_true(is.na(dt_qr_empty$L1_before), info = sprintf("[%s] empty L1 NA", cfg$label))

  ## sampling: CEM refinement on 100 treated only
  dt_qr_smp <- gd_run_match_pipeline(cfg$dt, params, sample_n_treated = 100L)  ## dt_qr_smp: quality report with sampling
  expect_true(is.data.table(dt_qr_smp), info = sprintf("[%s] smp data.table", cfg$label))
  expect_equal(nrow(dt_qr_smp), 1L, info = sprintf("[%s] smp one row", cfg$label))
  expect_equal(dt_qr_smp$n_matched, dt_qr2$n_matched,
               info = sprintf("[%s] smp match count unchanged", cfg$label))
  expect_true(dt_qr_smp$n_after_cem <= 100L,
              info = sprintf("[%s] smp ≤100 treated", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s pipeline", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_helpers <- function(cfg, refined) {
  ## sections 14-18: edge cases (quality report, coarse_bin, ematch_yr)
  t0 <- Sys.time()
  params <- list(L = cfg$L, treatment_col = cfg$treatment_col, uid = cfg$uid,
                 covariates = cfg$covariates)

  ## qr_s1: stage 1 — exact matching coverage (empty input)
  dt_s1 <- gd_quality_report_s1(cfg$dt, data.table(), params)
  expect_true(is.data.table(dt_s1), info = "s1_empty data.table")
  expect_equal(dt_s1$n_matched, 0L, info = "s1_empty n_matched=0")
  expect_equal(dt_s1$match_rate, 0, info = "s1_empty match_rate=0")
  expect_true(is.na(dt_s1$pct_tier1), info = "s1_empty pct_tier1 NA")

  ## qr_s2: stage 2 — balance diagnostics (NULL input)
  dt_s2 <- gd_quality_report_s2(NULL, NULL, params)
  expect_true(is.data.table(dt_s2), info = "s2_null data.table")
  expect_true(all(is.na(dt_s2$n_after_cem)), info = "s2_null n_after_cem NA")

  ## s2_nocov: no covariates intersect — tests non-NULL matchit path
  dt_s2_nc <- gd_quality_report_s2(refined$matchit, refined$data,
                                 list(covariates = "nonexistent"))  ## dt_s2_nc: stage 2 no covariates
  expect_true(is.data.table(dt_s2_nc), info = "s2_nocov data.table")
  expect_true(all(is.na(dt_s2_nc$L1_before)), info = "s2_nocov L1 NA")

  ## coarse_bin edges
  expect_equal(coarse_bin(c(1,2,3), 1), c(1L, 1L, 1L), info = "coarse_bin nbins=1")
  expect_equal(coarse_bin(rep(5, 10), 8), rep(1L, 10), info = "coarse_bin zero-range")
  expect_true(all(coarse_bin(1:100, 4) %in% 1:4), info = "coarse_bin 4 bins")
  expect_true(max(coarse_bin(rnorm(100), 10)) <= 10, info = "coarse_bin max ≤ 10")
  expect_true(min(coarse_bin(rnorm(100), 10)) >= 1,  info = "coarse_bin min ≥ 1")

  ## ematch_yr NULL
  dt_tmp <- copy(cfg$dt)[, `:=`(tid = get(cfg$tid),
                                 treat_bin = fifelse(get(cfg$treatment_col) > 0, 1L, 0L))]
  null_res <- gd_ematch_yr(1800, dt_tmp, NULL, "exn_pm", "matches", 1)
  expect_true(is.null(null_res), info = "ematch_yr NULL for empty year")

  timings <<- rbind(timings, data.table(label = sprintf("%s helpers", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_invariants <- function(cfg, mi) {
  ## sections 19-20: subclass/weight invariants, single cov/bin
  t0 <- Sys.time()
  r <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                          covariates = cfg$covariates, cutpoints = cfg$cutpoints)

  ## invariants
  s_ids <- unique(na.omit(r$data$subclass))
  expect_equal(length(s_ids), length(unique(s_ids)),
               info = sprintf("[%s] unique subclasses", cfg$label))
  sc <- r$data[!is.na(subclass), .(nt = sum(treat == 1), nc = sum(treat == 0)), by = subclass]  ## sc: subclass composition (treated/control counts)
  expect_true(all(sc$nt >= 1), info = sprintf("[%s] each subclass has treated", cfg$label))
  expect_true(all(sc$nc >= 1), info = sprintf("[%s] each subclass has control", cfg$label))
  expect_true(all(r$data[!is.na(subclass) & treat == 1, weight == 1]),
              info = sprintf("[%s] treated weight=1", cfg$label))
  expect_true(all(r$data[!is.na(subclass) & treat == 0, weight > 0]),
              info = sprintf("[%s] control weight>0", cfg$label))
  expect_true(all(r$data[is.na(subclass), weight == 0]),
              info = sprintf("[%s] unmatched weight=0", cfg$label))
  expect_equal(r$matchit$weights, r$data$weight,
               info = sprintf("[%s] matchit weights match", cfg$label))
  expect_equal(as.integer(r$matchit$subclass), as.integer(r$data$subclass),
               info = sprintf("[%s] matchit subclass matches", cfg$label))
  expect_true(inherits(r$matchit, "matchit"),
              info = sprintf("[%s] matchit class", cfg$label))

  ## single covariate
  r_1c <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                             covariates = cfg$covariates[1],
                             setNames(list(cfg$cutpoints[[1]]), cfg$covariates[1]))
  expect_true(is.list(r_1c) && all(c("data", "matchit") %in% names(r_1c)),
              info = sprintf("[%s] single cov list", cfg$label))
  expect_true("subclass" %in% names(r_1c$data),
              info = sprintf("[%s] single cov subclass", cfg$label))

  ## single bin per covariate
  cp_1b <- setNames(lapply(cfg$covariates, function(x) 1L), cfg$covariates)
  r_1b <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                             covariates = cfg$covariates, cutpoints = cp_1b)
  expect_equal(uniqueN(na.omit(r_1b$data$subclass)),
               uniqueN(r_1b$data[!is.na(subclass) & treat == 1, set]),
               info = sprintf("[%s] single bin = one per set", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s invariants", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_matchit_structure <- function(cfg, mi) {
  ## section 21: matchit version structure
  t0 <- Sys.time()
  r_mi <- gd_refine_matches_matchit(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                                     covariates = cfg$covariates,
                                     cutpoints = cfg$cutpoints)
  expect_true(is.list(r_mi), info = sprintf("[%s] matchit list", cfg$label))
  expect_true(all(c("data", "matchit") %in% names(r_mi)),
              info = sprintf("[%s] matchit slots", cfg$label))
  expect_true(inherits(r_mi$matchit, "matchit"),
              info = sprintf("[%s] matchit class", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s matchit struct", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_matchit_equiv <- function(cfg, mi) {
  ## sections 22, 26: DT vs matchit equivalence
  if (!run_slow || !cfg$slow) return()
  t0 <- Sys.time()

  ## equal-width bins (integer cutpoints)
  cp_int <- lapply(cfg$cutpoints, function(x) if (is.numeric(x)) as.integer(x) else x)
  r_dt <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                             covariates = cfg$covariates, cutpoints = cp_int)
  r_mi <- gd_refine_matches_matchit(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                                     covariates = cfg$covariates, cutpoints = cp_int)
  expect_true(gd_match_equivalent(r_dt, r_mi, uid = cfg$uid, tid = cfg$tid),
              info = sprintf("[%s] equal-width DT=MI", cfg$label))

  ## quantile breaks
  l_breaks <- gd_coarse_breaks(cfg$dt, cfg$covariates, cfg$cutpoints)
  r_dt2 <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                              covariates = cfg$covariates, cutpoints = l_breaks)
  r_mi2 <- gd_refine_matches_matchit(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                                      covariates = cfg$covariates, cutpoints = l_breaks)
  expect_true(gd_match_equivalent(r_dt2, r_mi2, uid = cfg$uid, tid = cfg$tid),
              info = sprintf("[%s] quantile DT=MI", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s DT=MI", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_coarse_breaks <- function(cfg) {
  ## section 23: coarse_breaks structure
  t0 <- Sys.time()
  br <- gd_coarse_breaks(cfg$dt, cfg$covariates, cfg$cutpoints)
  expect_true(is.list(br), info = sprintf("[%s] coarse_breaks list", cfg$label))
  expect_equal(names(br), cfg$covariates,
               info = sprintf("[%s] coarse_breaks names", cfg$label))
  expect_true(all(diff(br[[1]]) > 0),
              info = sprintf("[%s] coarse_breaks increasing", cfg$label))
  expect_true(br[[1]][1] <= min(cfg$dt[[cfg$covariates[1]]]),
              info = sprintf("[%s] coarse_breaks first ≤ min", cfg$label))

  ## 1 bin
  cp1 <- setNames(list(1), cfg$covariates[1])
  br1 <- gd_coarse_breaks(cfg$dt, cfg$covariates[1], cp1)
  expect_equal(length(br1[[1]]), 2,
               info = sprintf("[%s] 1 bin = 2 breakpoints", cfg$label))

  ## fewer unique values than bins
  br_uniq <- gd_coarse_breaks(cfg$dt[1:10], cfg$covariates[1],
                               setNames(list(20), cfg$covariates[1]))
  expect_true(length(br_uniq[[1]]) <= 10,
              info = sprintf("[%s] ≤n unique breaks", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s coarse_breaks", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_quantile_refine <- function(cfg, mi) {
  ## sections 24-25: quantile breaks internal consistency
  t0 <- Sys.time()
  l_breaks <- gd_coarse_breaks(cfg$dt, cfg$covariates, cfg$cutpoints)
  r_q <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                            covariates = cfg$covariates, cutpoints = l_breaks)
  expect_true(is.list(r_q), info = sprintf("[%s] quantile refine list", cfg$label))
  expect_true("subclass" %in% names(r_q$data),
              info = sprintf("[%s] quantile subclass", cfg$label))
  expect_true(inherits(r_q$matchit, "matchit"),
              info = sprintf("[%s] quantile matchit class", cfg$label))
  sc_q <- r_q$data[!is.na(subclass), .(nt = sum(treat == 1), nc = sum(treat == 0)), by = subclass]  ## sc_q: subclass composition, quantile breaks
  expect_true(all(sc_q$nt >= 1), info = sprintf("[%s] quantile has treated", cfg$label))
  expect_true(all(sc_q$nc >= 1), info = sprintf("[%s] quantile has control", cfg$label))
  expect_true(all(r_q$data[!is.na(subclass) & treat == 1, weight == 1]),
              info = sprintf("[%s] quantile treated weight=1", cfg$label))
  expect_true(all(r_q$data[!is.na(subclass) & treat == 0, weight > 0]),
              info = sprintf("[%s] quantile control weight>0", cfg$label))
  expect_true(all(r_q$data[is.na(subclass), weight == 0]),
              info = sprintf("[%s] quantile unmatched weight=0", cfg$label))

  ## reproducibility: same quantile breaks → identical subclass assignments
  l_breaks2 <- gd_coarse_breaks(cfg$dt, cfg$covariates, cfg$cutpoints)
  r_q2 <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                             covariates = cfg$covariates, cutpoints = l_breaks2)
  fp1 <- r_q$data[!is.na(subclass),
                   .(uid = get(cfg$uid), tid = get(cfg$tid), treat, s = subclass)][order(uid, tid, treat)]
  fp2 <- r_q2$data[!is.na(subclass),
                    .(uid = get(cfg$uid), tid = get(cfg$tid), treat, s = subclass)][order(uid, tid, treat)]
  expect_equal(nrow(fp1), nrow(fp2),
               info = sprintf("[%s] quantile reproducible row count", cfg$label))
  expect_true(all(fp1$s == fp2$s),
              info = sprintf("[%s] quantile reproducible assignments", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s quantile", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

# ==== Orchestrator ====
test_config <- function(cfg) {
  mi   <- gd_test_basic(cfg)
  ref  <- gd_test_cem_refine(cfg, mi)
  dt_qr   <- gd_test_quality_report(cfg, mi, ref)
  gd_test_smd_vr(cfg, ref)
  gd_test_l1_crossval(cfg, mi)
  gd_test_pipeline(cfg, mi, ref, dt_qr)
  gd_test_helpers(cfg, ref)
  gd_test_invariants(cfg, mi)
  gd_test_matchit_structure(cfg, mi)
  gd_test_matchit_equiv(cfg, mi)
  gd_test_coarse_breaks(cfg)
  gd_test_quantile_refine(cfg, mi)
}

# ==== Main execution ====
run_slow <- Sys.getenv("FOSTER_TEST_SLOW", "true") != "false"

for (cfg in configs) {
  test_config(cfg)
}

# --- Timing report ---
cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings, n = 800)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
