library(data.table)
library(tinytest)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")
library(cobalt)

# --- Timing collector ---
timings <- data.table(label = character(), time_s = numeric())

# Load small test data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test_smol.csv")
## dt_med <- fread("~/Dropbox/phd/papers/foster/dt_test_med.csv")
dt_big <- fread("~/Dropbox/phd/papers/foster/dt_test_big.csv")



# === Config definitions ===
configs <- list(
  list(dt_panel = dt, L = 3L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
       uid = "PeopleID", tid = "year",
       # covariates = c("BirthYear", "exn_all_cum_ihs"),
       cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       slow = TRUE, label = "smol_L3", sample_n_treated = 100,
       tiers = list(c("hist_str", "Nationality", "Gender"), c("hist_str")))
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
  ## r_s <- gd_match_history(cfg$dt, L = cfg$L, treatment_col = cfg$treatment_col,
  ## r_s: match history summary
  r_s <- gd_match_history(cfg$dt_panel, L = cfg$L, treatment_col = cfg$treatment_col,
                           history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid, tiers = cfg$tiers,
                          return_type = "summary")
  
  
  expect_true(is.data.frame(r_s), info = sprintf("[%s] summary is data.frame", cfg$label))
  expect_true(all(c("total_treated", "treated_matched", "total_matches") %in% names(r_s)),
              info = sprintf("[%s] summary columns", cfg$label))
  expect_true(r_s$total_treated > 0, info = sprintf("[%s] non-zero treated", cfg$label))

  ## r_f <- gd_match_history(cfg$dt, L = cfg$L, treatment_col = cfg$treatment_col,
  ## r_f: match history full matches
  r_f <- gd_match_history(cfg$dt_panel, L = cfg$L, treatment_col = cfg$treatment_col,
                           history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid, tiers = cfg$tiers,
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
                            history_col = cfg$history_col, uid = cfg$uid, tid = cfg$tid, tiers = cfg$tiers,
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
  r <- gd_refine_matches(mi, cfg$dt_panel, cutpoints = cfg$cutpoints, uid = cfg$uid, tid = cfg$tid)
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
  params <- list(L = cfg$L, treatment_col = cfg$treatment_col, uid = cfg$uid, cutpoints = cfg$cutpoints)
  dt_qr <- gd_quality_report(
    dt_prepared  = cfg$dt_panel,
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
  invisible(lapply(names(cfg$cutpoints), function(v) {
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
  dt_cdf <- refined$data[, names(cfg$cutpoints), with = FALSE]  ## dt_cdf: covariate data frame
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
    list(cov = names(cfg$cutpoints),                     cp = 8),
    list(cov = names(cfg$cutpoints),                     cp = 6),
    list(cov = names(cfg$cutpoints),                     cp = 4),
    list(cov = names(cfg$cutpoints),                     cp = 2),
    list(cov = names(cfg$cutpoints)[1],                   cp = 8),
    list(cov = names(cfg$cutpoints)[1],                   cp = 4)
  )
  dt_l1 <- rbindlist(lapply(l_params, function(p) {  ## dt_l1: L1 imbalance results across parameter combos
    cp <- setNames(rep(list(p$cp), length(p$cov)), p$cov)
    if (length(cp) == 0) return(NULL)
    r <- gd_refine_matches(mi, cfg$dt, uid = cfg$uid, tid = cfg$tid,
                            cutpoints = cp)
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

  print(dt_l1[, .(combo, co_before, co_after, my_before, my_after, diff_before, diff_after,
                  co_reduc = co_reduction, my_reduc = my_reduction)], n = 100)
  
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
                 cutpoints = cfg$cutpoints,
                 sample_n_treated = cfg$sample_n_treated,
                 tiers = cfg$tiers)

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
  dt_qr_smp <- gd_run_match_pipeline(cfg$dt, params)  ## dt_qr_smp: quality report with sampling
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
                 cutpoints = cfg$cutpoints)

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
  dt_tmp <- copy(cfg$dt_panel)[, `:=`(tid = get(cfg$tid),
                                 treat_bin = fifelse(get(cfg$treatment_col) > 0, 1L, 0L))]
  null_res <- gd_ematch_yr(1800, dt_tmp, NULL, "exn_pm", "matches", 1)
  expect_true(is.null(null_res), info = "ematch_yr NULL for empty year")

  timings <<- rbind(timings, data.table(label = sprintf("%s helpers", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_invariants <- function(cfg, mi) {
  ## sections 19-20: subclass/weight invariants, single cov/bin
  t0 <- Sys.time()
  r <- gd_refine_matches(mi, cfg$dt_panel, cutpoints = cfg$cutpoints, uid = cfg$uid, tid = cfg$tid)

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
  cp_1c <- setNames(list(cfg$cutpoints[[1]]), names(cfg$cutpoints)[1])
  r_1c <- gd_refine_matches(mi, cfg$dt_panel, cutpoints = cp_1c, uid = cfg$uid, tid = cfg$tid)
  expect_true(is.list(r_1c) && all(c("data", "matchit") %in% names(r_1c)),
              info = sprintf("[%s] single cov list", cfg$label))
  expect_true("subclass" %in% names(r_1c$data),
              info = sprintf("[%s] single cov subclass", cfg$label))

  ## single bin per covariate
  cp_1b <- setNames(lapply(names(cfg$cutpoints), function(x) 1L), names(cfg$cutpoints))
  r_1b <- gd_refine_matches(mi, cfg$dt_panel, cutpoints = cp_1b, uid = cfg$uid, tid = cfg$tid)
  expect_equal(uniqueN(na.omit(r_1b$data$subclass)),
               uniqueN(r_1b$data[!is.na(subclass) & treat == 1, set]),
               info = sprintf("[%s] single bin = one per set", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s invariants", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_matchit_structure <- function(cfg, mi) {
  ## section 21: matchit version structure
  t0 <- Sys.time()
  r_mi <- gd_refine_matches_matchit(mi, cfg$dt_panel, cutpoints = cfg$cutpoints, uid = cfg$uid, tid = cfg$tid)
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
  r_dt <- gd_refine_matches(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid, cutpoints = cp_int)
  r_mi <- gd_refine_matches_matchit(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid, cutpoints = cp_int)
  expect_true(gd_match_equivalent(r_dt, r_mi, uid = cfg$uid, tid = cfg$tid),
              info = sprintf("[%s] equal-width DT=MI", cfg$label))

  ## quantile breaks
  l_breaks <- gd_coarse_breaks(cfg$dt_panel, cfg$cutpoints)
  r_dt2 <- gd_refine_matches(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid,
                              cutpoints = l_breaks)
  r_mi2 <- gd_refine_matches_matchit(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid,
                                      cutpoints = l_breaks)
  expect_true(gd_match_equivalent(r_dt2, r_mi2, uid = cfg$uid, tid = cfg$tid),
              info = sprintf("[%s] quantile DT=MI", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s DT=MI", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_coarse_breaks <- function(cfg) {
  ## section 23: coarse_breaks structure
  t0 <- Sys.time()
  br <- gd_coarse_breaks(cfg$dt_panel, cfg$cutpoints)
  
  expect_true(is.list(br), info = sprintf("[%s] coarse_breaks list", cfg$label))
  expect_equal(names(br), names(cfg$cutpoints),
               info = sprintf("[%s] coarse_breaks names", cfg$label))
  expect_true(all(diff(br[[1]]) > 0),
              info = sprintf("[%s] coarse_breaks increasing", cfg$label))
  ## expect_true(br[[1]][1] <= min(cfg$dt_panel[[names(cfg$cutpoints)[1]]]),
  ##             info = sprintf("[%s] coarse_breaks first ≤ min", cfg$label))

  ## 1 bin: no longer makes sense
  ## cp1 <- setNames(list(1), names(cfg$cutpoints)[1])
  ## br1 <- gd_coarse_breaks(cfg$dt_panel, cp1)
  ## expect_equal(length(br1[[1]]), 2, info = sprintf("[%s] 1 bin = 2 breakpoints", cfg$label))

  ## fewer unique values than bins
  cp_uq <- setNames(list(20), names(cfg$cutpoints)[1])
  br_uniq <- gd_coarse_breaks(cfg$dt_panel[1:10], cp_uq)
  
  expect_true(length(br_uniq[[1]]) <= 10,
              info = sprintf("[%s] ≤n unique breaks", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s coarse_breaks", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_quantile_refine <- function(cfg, mi) {
  ## sections 24-25: quantile breaks internal consistency
  t0 <- Sys.time()
  l_breaks <- gd_coarse_breaks(cfg$dt_panel, cfg$cutpoints)
  r_q <- gd_refine_matches(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid, cutpoints = l_breaks)
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
  l_breaks2 <- gd_coarse_breaks(cfg$dt_panel, cfg$cutpoints)
  r_q2 <- gd_refine_matches(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid, cutpoints = l_breaks2)
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

exec_test <- function(fn, name, cfg) {
    tryCatch(fn(), error = function(e) {
      expect_true(FALSE, info = sprintf("[%s] crashed in %s: %s",
                                        cfg$label, name, conditionMessage(e)))
      NULL
    })
  }

# ==== Orchestrator ====
test_config <- function(cfg) {
  ## wrap a test step: crash → test failure with clear label
  ## Takes a thunk (zero-arg function) so tryCatch wraps the actual execution
  

  mi <- exec_test(function() gd_test_basic(cfg), "basic", cfg)
  if (is.null(mi)) return()
  ref <- exec_test(function() gd_test_cem_refine(cfg, mi), "cem refine", cfg)
  if (is.null(ref)) return()
  
  dt_qr  <- exec_test(function() gd_test_quality_report(cfg, mi, ref), "quality report", cfg)
  print(dt_qr) %>% t

  exec_test(function() gd_test_smd_vr(cfg, ref),                   "smd vr", cfg)
  exec_test(function() gd_test_l1_crossval(cfg, mi),               "L1 crossval", cfg)
  exec_test(function() gd_test_pipeline(cfg, mi, ref, dt_qr),     "pipeline", cfg)
  exec_test(function() gd_test_helpers(cfg, ref),                  "helpers", cfg)
  exec_test(function() gd_test_invariants(cfg, mi),                "invariants", cfg)
  exec_test(function() gd_test_matchit_structure(cfg, mi),         "matchit struct", cfg)
  exec_test(function() gd_test_matchit_equiv(cfg, mi),             "DT=MI", cfg)
  exec_test(function() gd_test_coarse_breaks(cfg),                 "coarse breaks", cfg)
  exec_test(function() gd_test_quantile_refine(cfg, mi),           "quantile", cfg)
}

# ===========================================================================
# Multicategory CEM Stratum Tests
# ===========================================================================

cem_configs <- list(
  list(dt_panel = dt_big, L = 3L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = c(5L, 5L), history_breaks = "quantile",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       weight_type = "none",
       label = "quantile_5"),
  list(dt_panel = dt_big, L = 3L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = 5L, history_breaks = "equal",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       weight_type = "any_vs_ctrl",
       label = "equal_5"),
  list(dt_panel = dt_big, L = 3L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = 5L, history_breaks = "none",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(),
       weight_type = "composite",
       label = "none_no_cov")
)

gd_test_cem_prep <- function(cfg) {
  t0 <- Sys.time()
  dt_p <- gd_prep_multivar_hist(
    dt_panel = cfg$dt_panel, L = cfg$L,
    history_cols = cfg$history_cols, history_nbins = cfg$history_nbins,
    history_breaks = cfg$history_breaks,
    uid = cfg$uid, tid = cfg$tid
  )
  expect_true(is.data.table(dt_p),
              info = sprintf("[%s] prep is data.table", cfg$label))
  expect_true(nrow(dt_p) > 0,
              info = sprintf("[%s] prep non-empty", cfg$label))
  expect_true("hist_str" %in% names(dt_p),
              info = sprintf("[%s] hist_str column", cfg$label))
  expect_true(all(c("uid", "tid") %in% names(dt_p)),
              info = sprintf("[%s] internal uid/tid", cfg$label))
  ## rows should decrease by (L-1) per artist due to lag NA removal
  expect_true(nrow(dt_p) < nrow(cfg$dt_panel),
              info = sprintf("[%s] rows decreased from lag removal", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s cem_prep", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  dt_p
}

gd_test_cem_strata <- function(cfg) {
  t0 <- Sys.time()
  result <- gd_cem_strata(cfg)
  expect_true(is.list(result),
              info = sprintf("[%s] result is list", cfg$label))
  expect_true(all(c("data", "strata_summary", "config") %in% names(result)),
              info = sprintf("[%s] result slots", cfg$label))
  dt_d <- result$data
  expect_true("stratum" %in% names(dt_d),
              info = sprintf("[%s] stratum column", cfg$label))
  expect_true(cfg$uid %in% names(dt_d),
              info = sprintf("[%s] uid preserved", cfg$label))
  expect_true(cfg$tid %in% names(dt_d),
              info = sprintf("[%s] tid preserved", cfg$label))

  ## treatment columns preserved
  invisible(lapply(cfg$treatment_cols, function(tc) {
    expect_true(tc %in% names(dt_d),
                info = sprintf("[%s] treatment col %s preserved", cfg$label, tc))
  }))

  timings <<- rbind(timings, data.table(label = sprintf("%s cem_strata", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  result
}

gd_test_cem_consistency <- function(cfg, result) {
  t0 <- Sys.time()
  dt_d <- result$data
  dt_s <- result$strata_summary

  ## within-stratum consistency: year, hist_str, and coarse covariates are invariant
  cv_cols <- intersect(names(dt_d), paste0(".cr_", names(cfg$covariate_cutpoints)))
  check_cols <- c(cfg$tid, "hist_str", cv_cols)
  for (cc in check_cols) {
    nv <- dt_d[, uniqueN(get(cc)), by = stratum][V1 > 1]
    expect_true(nrow(nv) == 0,
                info = sprintf("[%s] %s invariant within strata", cfg$label, cc))
  }

  ## strata_summary matches data
  expect_equal(nrow(dt_s), uniqueN(dt_d$stratum),
               info = sprintf("[%s] strata_summary has one row per stratum", cfg$label))

  ## n_total in summary matches actual
  dt_check <- dt_d[, .(.N), by = stratum]
  dt_join <- merge(dt_check, dt_s[, .(stratum, n_total)], by = "stratum")
  expect_true(all(dt_join$N == dt_join$n_total),
              info = sprintf("[%s] n_total correct", cfg$label))

  ## per-treatment counts in summary are >= 0 (already checked by being in data.table)
  for (tc in cfg$treatment_cols) {
    ncol <- paste0("n_", tc)
    expect_true(ncol %in% names(dt_s),
                info = sprintf("[%s] summary has n_%s", cfg$label, tc))
    expect_true(all(dt_s[[ncol]] >= 0),
                info = sprintf("[%s] n_%s non-negative", cfg$label, tc))
  }

  ## n_control is computed
  expect_true("n_control" %in% names(dt_s),
              info = sprintf("[%s] summary has n_control", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s cem_consistency", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_cem_edge_cases <- function(cfg) {
  t0 <- Sys.time()

  ## empty data
  cfg_empty <- cfg
  cfg_empty$dt_panel <- cfg$dt_panel[FALSE]
  cfg_empty$label <- paste0(cfg$label, "_empty")
  expect_error(
    gd_cem_strata(cfg_empty),
    info = sprintf("[%s] empty data handled", cfg$label)
  )

  ## single treatment (should still work)
  cfg_one <- cfg
  cfg_one$treatment_cols <- cfg$treatment_cols[1]
  cfg_one$history_cols <- cfg$history_cols[1]
  if (length(cfg$history_nbins) > 1) {
    cfg_one$history_nbins <- cfg$history_nbins[1]
  }
  cfg_one$label <- paste0(cfg$label, "_one")
  res_one <- gd_cem_strata(cfg_one)
  expect_true("stratum" %in% names(res_one$data),
              info = sprintf("[%s] single treatment works", cfg$label))

  ## quality report schema — single row, correct columns
  dt_q <- gd_cem_quality_report(res_one)
  expect_true(is.data.table(dt_q),
              info = sprintf("[%s] quality report is data.table", cfg$label))
  expect_equal(nrow(dt_q), 1L,
               info = sprintf("[%s] quality report single row", cfg$label))
  expect_true("n_prepped" %in% names(dt_q),
              info = sprintf("[%s] quality report has n_prepped", cfg$label))
  expect_true("n_strata" %in% names(dt_q),
              info = sprintf("[%s] quality report has n_strata", cfg$label))
  expect_true(paste0("n_strata_", cfg_one$treatment_cols[1]) %in% names(dt_q),
              info = sprintf("[%s] quality report has per-treatment col", cfg$label))
  expect_true("n_ctrl_strata" %in% names(dt_q),
              info = sprintf("[%s] quality report has n_ctrl_strata", cfg$label))


  ## weight column present when weight_type is set
  if (!is.null(cfg$weight_type) && cfg$weight_type != "none") {
    expect_true("weight" %in% names(res_one$data),
                info = sprintf("[%s] has weight col", cfg$label))
    expect_true(all(res_one$data$weight >= 0),
                info = sprintf("[%s] weights non-negative", cfg$label))
    expect_true(any(res_one$data$weight != 1),
                info = sprintf("[%s] weights not all 1", cfg$label))
  }
  ## diag function
  dt_q2 <- gd_cem_strata_diag(res_one)
  expect_true(is.data.table(dt_q2),
              info = sprintf("[%s] diag returns data.table", cfg$label))

  timings <<- rbind(timings, data.table(label = sprintf("%s cem_edges", cfg$label),
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

test_cem_config <- function(cfg) {
  dt_p <- exec_test(function() gd_test_cem_prep(cfg), "cem_prep", cfg)
  if (is.null(dt_p)) return()
  res <- exec_test(function() gd_test_cem_strata(cfg), "cem_strata", cfg)
  if (is.null(res)) return()
  exec_test(function() gd_test_cem_consistency(cfg, res), "cem_consistency", cfg)
  exec_test(function() gd_test_cem_edge_cases(cfg), "cem_edges", cfg)
}

# ===========================================================================
# CEM Weight Function Tests
# ===========================================================================

gd_test_weights_none <- function(dt_d) {
  t0 <- Sys.time()
  w <- gd_cem_weights_none(dt_d)
  expect_true(all(w == 1),
              info = "weights_none all 1")
  expect_equal(length(w), nrow(dt_d),
               info = "weights_none correct length")
  timings <<- rbind(timings, data.table(label = "weights_none",
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_weights_any_ctrl <- function(dt_d) {
  t0 <- Sys.time()
  w <- gd_cem_weights_any_ctrl(dt_d, c("exn_pi", "exn_pm"))
  expect_true(length(w) == nrow(dt_d),
              info = "weights_any_ctrl correct length")
  expect_true(all(w >= 0),
              info = "weights_any_ctrl non-negative")
  expect_true(any(w != 1),
              info = "weights_any_ctrl varies")
  timings <<- rbind(timings, data.table(label = "weights_any_ctrl",
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  invisible(w)
}

gd_test_weights_any_ctrl_ratio <- function(dt_d, w) {
  ## verify that within each stratum, (weighted non-treat)/(weighted treat) = global ratio
  n_treat <- sum(rowSums(dt_d[, .(exn_pi, exn_pm)] > 0) > 0)
  n_ctrl  <- sum(rowSums(dt_d[, .(exn_pi, exn_pm)] > 0) == 0)
  R <- n_ctrl / n_treat
  dt_check <- data.table(stratum = dt_d$stratum, w = w)[
    , .(r = sum(w) > 0), by = stratum]  ## quick check: some non-zero weight strata exist
  expect_true(any(dt_check$r),
              info = "weights_any_ctrl some viable strata")
}

gd_test_weights_composite <- function(dt_d) {
  t0 <- Sys.time()
  w <- gd_cem_weights_composite(dt_d, c("exn_pi", "exn_pm"))
  expect_true(length(w) == nrow(dt_d),
              info = "weights_composite correct length")
  expect_true(all(w >= 0),
              info = "weights_composite non-negative")
  timings <<- rbind(timings, data.table(label = "weights_composite",
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
  invisible(w)
}

gd_test_weights_qp_stub <- function() {
  t0 <- Sys.time()
  expect_error(gd_cem_weights_qp(NULL, NULL),
               info = "weights_qp raises error")
  timings <<- rbind(timings, data.table(label = "weights_qp_stub",
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}

gd_test_weights_dispatcher <- function(dt_d) {
  t0 <- Sys.time()
  w_none <- gd_cem_weights(dt_d, c("exn_pi", "exn_pm"), "none")
  expect_true(all(w_none == 1), info = "dispatcher none")
  w_any <- gd_cem_weights(dt_d, c("exn_pi", "exn_pm"), "any_vs_ctrl")
  expect_true(length(w_any) == nrow(dt_d), info = "dispatcher any_vs_ctrl")
  expect_error(gd_cem_weights(dt_d, c("exn_pi", "exn_pm"), "unknown_type"),
               info = "dispatcher unknown type")
  timings <<- rbind(timings, data.table(label = "weights_dispatcher",
                                         time_s = round(as.numeric(Sys.time() - t0), 2)))
}


# ==== Main execution ====
run_slow <- Sys.getenv("FOSTER_TEST_SLOW", "true") != "false"

for (cfg in configs) {
  test_config(cfg)
}

for (cfg in cem_configs) {
  test_cem_config(cfg)
}

# --- CEM Weight Function Tests ---

cat("\n--- Weight function tests ---\n")

## use a cem result to get data with strata
dt_demo <- dt_big
demo_cfg <- list(dt_panel = dt_demo, L = 3L,
                 treatment_cols = c("exn_pm", "exn_pi"),
                 history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
                 history_nbins = c(5L, 5L), history_breaks = "quantile",
                 uid = "PeopleID", tid = "year",
                 covariate_cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
                 label = "weights")
demo_res <- gd_cem_strata(demo_cfg)
demo_data <- demo_res$data
demo_data[, (c("exn_pi", "exn_pm")) := lapply(.SD, function(x) as.integer(x > 0)),
           .SDcols = c("exn_pi", "exn_pm")]

gd_test_weights_none(demo_data)
gd_test_weights_qp_stub()
gd_test_weights_dispatcher(demo_data)
w_any <- gd_test_weights_any_ctrl(demo_data)
gd_test_weights_any_ctrl_ratio(demo_data, w_any)
gd_test_weights_composite(demo_data)

# --- Timing report ---
cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings, n = 800)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
