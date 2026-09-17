library(tinytest)
suppressMessages({
  library(data.table)
  library(purrr)
  library(magrittr)
  library(optweight)
  library(fixest)
  library(broom)
})
source("~/Dropbox/phd/papers/foster/scripts/mySBW.R")
cat("test_mySBW.R loaded\n")

## run via: Rscript -e 'tinytest::run_test_file("scripts/test_mySBW.R")'

# --- Timing collector ---
timings <- data.table(label = character(), time_s = numeric())
t_test <- function(label, expr) {
  t0 <- Sys.time()
  out <- force(expr)
  timings <<- rbind(timings, data.table(label = label,
                                        time_s = round(as.numeric(Sys.time() - t0), 2)))
  out
}

# --- Synthetic panel: complete grid, entry at min(years), history-driven counts ---
gd_syn_dt <- function(n_artists = 2500, years = 2006:2014, seed = 7) {
  set.seed(seed)
  yr_min <- min(years); yr_max <- max(years); n <- yr_max - yr_min + 1
  map_dfr(seq_len(n_artists), function(i) {
    byr <- sample(1940:1990, 1)
    pm <- pi <- oes <- toti <- integer(n)
    pm[1] <- rbinom(1, 1, 0.3); pi[1] <- rpois(1, 0.8); oes[1] <- rpois(1, 1)
    toti[1] <- rbinom(1, 1, 0.1)
    for (j in 2:n) {
      pm[j]   <- rpois(1, 0.05 + 0.55*pm[j-1] + 0.10*pi[j-1])
      pi[j]   <- rpois(1, 0.15 + 0.45*pi[j-1] + 0.10*pm[j-1])
      oes[j]  <- rpois(1, 0.30 + 0.35*oes[j-1] + 0.05*pi[j-1])
      toti[j] <- rbinom(1, 1, pmin(0.6, 0.05 + 0.3*pm[j-1] + 0.2*toti[j-1]))
    }
    ## ar1: artist-specific time-varying numeric covariate (not deterministic
    ## in year, so its lags are informative but not collinear with the level)
    ar1 <- numeric(n); ar1[1] <- rnorm(1)
    for (j in 2:n) ar1[j] <- 0.6 * ar1[j - 1] + rnorm(1)

    data.table(PeopleID = i, year = yr_min:yr_max,
               exn_pm = pm, exn_pi = pi, exn_oes = oes, exn_toti = toti,
               Gender = sample(c("Male", "Female"), 1),
               BirthYear = byr,
               ar1 = round(ar1, 3),
               auc_turnover = cumsum(pm + pi),   ## deterministic "auction" proxy
               exn_all_cum_ihs = asinh(cumsum(pm + pi + oes + toti)))
  })
}

## mimic gl_sbw_proc's preparation chain on a config: panel + factor
## expansion + covariate lags; returns the prepared config (ready for
## gf_sbw_formulas / gd_sbw_years)
gd_prep_cfg <- function(c_sbw) {
  c_sbw <- md_sbw_defaults(c_sbw)
  dt_p <- md_sbw_panel(c_sbw)
  l_exp <- md_sbw_expand_covs(dt_p, c_sbw$covariate_cols)
  c_sbw$covariate_cols <- l_exp$covs
  l_lags <- md_sbw_cov_lags(l_exp$dt, c_sbw)
  c_sbw$covariate_cols <- l_lags$covs
  list(cfg = c_sbw, dt = l_lags$dt)
}

## test cache lives in a per-session temp dir
DIR_CACHE_TEST <- file.path(tempdir(), sprintf("sbw_test_cache_%d", Sys.getpid()))

gd_syn_config <- function(dt, ...) {
  utils::modifyList(list(
    label = "test_syn", dt_panel = dt, uid = "PeopleID", tid = "year",
    treatment_cols = c("exn_pm", "exn_pi", "exn_oes", "exn_toti"),
    covariate_cols = "Gender",
    n_lag = 2L, n_lag_cov = 0L, impact_years = 2L,
    years = 2009:2010, min_n_year = 50L,
    estimand = "ATE", norm = "l2", eps = 1e-3, verbose = FALSE,
    n_cores = 1L, dir_cache = DIR_CACHE_TEST, wcol = "w_sbw"), list(...))
}

# === config validation ===
gd_test_check <- function() {
  dt <- gd_syn_dt(50, 2006:2009, seed = 1)
  c_ok <- gd_syn_config(dt, years = NULL, min_n_year = 5L)

  expect_true(gd_sbw_check(md_sbw_defaults(c_ok)), info = "valid config passes")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(label = "bad/label"))),
               pattern = "label", info = "label regex")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(treatment_cols = "exn_nope"))),
               pattern = "treatment", info = "treatment col exists")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(treatment_cols = "Gender"))),
               pattern = "numeric", info = "treatment numeric")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(estimand = "ATT"))),
               pattern = "ATE", info = "estimand ATE only")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(n_lag = 0))),
               pattern = "n_lag", info = "n_lag >= 1")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(covariate_cols = "exn_pm"))),
               pattern = "overlap", info = "cov/treat disjoint")
  ## mixed numeric + factor covariate set with n_lag_cov > 0 is supported:
  ## numeric covs (and later dummies) are lagged, factors balanced at level t
  expect_true(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(covariate_cols = c("Gender", "ar1"),
                                 n_lag_cov = 2L)))),
    info = "mixed covariate set + n_lag_cov passes check")
  expect_error(gd_sbw_check(utils::modifyList(c_ok, list(n_lag_cov = -1L))),
               pattern = "n_lag_cov", info = "n_lag_cov >= 0")

  ## gap in the grid -> error (note: replace dt_panel by direct assignment;
  ## utils::modifyList recurses into data.table members and breaks)
  c_gap <- c_ok; c_gap$dt_panel <- dt[!(PeopleID == 1 & year == 2007)]
  expect_error(gd_sbw_check(c_gap),
               pattern = "complete per-artist grid", info = "complete grid required")
  ## duplicate rows -> error
  c_dup <- c_ok; c_dup$dt_panel <- rbind(dt, dt[PeopleID == 1 & year == 2006])
  expect_error(gd_sbw_check(c_dup),
               pattern = "duplicate", info = "duplicate (uid, tid) rejected")
}

# === formula generation (after factor-covariate expansion + covariate lags) ===
gd_test_formulas <- function() {
  dt <- gd_syn_dt(50, 2006:2009, seed = 1)
  c_sbw <- gd_syn_config(dt)
  ## mimic gl_sbw_proc: build panel, expand factor covariates, lag covariates
  l_prep <- gd_prep_cfg(c_sbw)
  c_sbw <- l_prep$cfg
  l_fs <- gf_sbw_formulas(c_sbw)

  expect_equal(length(l_fs), 4L, info = "one formula per treatment")
  expect_equal(names(l_fs), c_sbw$treatment_cols, info = "formulas named by treatment")
  rhs <- as.character(l_fs[[1]])[3]
  expect_true(grepl("exn_pm_l1", rhs) && grepl("exn_pm_l2", rhs), info = "own lags in RHS")
  expect_true(grepl("exn_toti_l1", rhs), info = "other-treatment lags in RHS")
  expect_true(grepl("Gender_", rhs), info = "covariates in RHS (expanded dummies)")
  expect_false(grepl("\\bGender\\b", rhs), info = "no bare factor name in RHS")
  ## no same-year treatment on the RHS (only lagged versions)
  expect_false(any(map_lgl(c_sbw$treatment_cols,
                           ~grepl(sprintf("\\b%s\\b(?!_)", .x), rhs, perl = TRUE))),
               info = "no current-year treatments in RHS")

  ## covariate lags enter the RHS for a numeric time-varying covariate
  c_cov <- gd_syn_config(dt, covariate_cols = "ar1", n_lag_cov = 2L)
  l_cov <- gd_prep_cfg(c_cov)
  rhs2 <- as.character(gf_sbw_formulas(l_cov$cfg)[[1]])[3]
  expect_true(grepl("ar1_l1", rhs2) && grepl("ar1_l2", rhs2),
              info = "numeric time-varying covariate lags in RHS")
  expect_true("ar1_l1" %in% l_cov$cfg$covariate_cols &&
              "ar1_l2" %in% l_cov$cfg$covariate_cols,
              info = "numeric covariate lags kept in final covariate set")

  ## time-invariant numeric covariate: lags are redundant and dropped
  c_ti <- gd_syn_config(dt, covariate_cols = "BirthYear", n_lag_cov = 2L)
  l_ti <- gd_prep_cfg(c_ti)
  expect_false(any(grepl("^BirthYear_l", l_ti$cfg$covariate_cols)),
               info = "time-invariant numeric covariate lags dropped")
}

# === factor covariate expansion ===
gd_test_covexp <- function() {
  dt <- gd_syn_dt(50, 2006:2009, seed = 1)
  c_sbw <- gd_syn_config(dt)               ## covariate_cols = "Gender" (character)
  dt_p <- md_sbw_panel(c_sbw)
  l_exp <- md_sbw_expand_covs(dt_p, c_sbw$covariate_cols)

  ## levels sorted c("Female", "Male") -> "Female" is the dropped reference
  expect_true("Gender_Male" %in% l_exp$covs, info = "non-reference dummy created")
  expect_false("Gender_Female" %in% l_exp$covs, info = "reference level dropped")
  expect_false("Gender" %in% l_exp$covs, info = "original factor replaced")
  expect_false("Gender" %in% names(l_exp$dt), info = "factor column dropped from panel")
  expect_true(is.numeric(l_exp$dt$Gender_Male), info = "dummy columns numeric")
  expect_equal(nrow(l_exp$dt), nrow(dt_p), info = "expansion keeps row count")
  expect_equal(sum(l_exp$dt$Gender_Male == 1), sum(dt$Gender == "Male"),
               info = "dummy counts match levels")
  ## numeric covariates untouched
  c_num <- gd_syn_config(dt, covariate_cols = "exn_all_cum_ihs")
  dt_p2 <- md_sbw_panel(c_num)
  l_num <- md_sbw_expand_covs(dt_p2, "exn_all_cum_ihs")
  expect_equal(l_num$covs, "exn_all_cum_ihs", info = "numeric covariate kept as-is")
  expect_true("exn_all_cum_ihs" %in% names(l_num$dt), info = "numeric covariate retained")
}

# === panel builder ===
gd_test_panel <- function() {
  dt_h <- data.table(PeopleID = rep(1:2, each = 5), year = rep(2001:2005, 2),
                     exn_pm = c(1, 2, 3, 4, 5, 0, 0, 0, 0, 0),
                     exn_pi = c(0, 0, 0, 0, 0, 10, 20, 30, 40, 50))
  c_h <- gd_syn_config(dt_h, treatment_cols = c("exn_pm", "exn_pi"),
                       covariate_cols = NULL, years = 2003:2004, label = "hand")
  dt <- md_sbw_panel(c_h)

  expect_true(is.data.table(dt), info = "panel is data.table")
  ## raw future leads are KEPT (usable e.g. for Poisson outcomes and needed
  ## as the outcome-coverage gate in gd_sbw_years)
  expect_true(all(c("exn_pm_f1", "exn_pm_f2") %in% names(dt)), info = "leads retained")

  ## artist 1: pm = 1..5 over 2001..2005
  d1 <- dt[PeopleID == 1]
  expect_equal(d1[year == 2003, exn_pm_l1], 2, info = "lag1 value")
  expect_equal(d1[year == 2003, exn_pm_l2], 1, info = "lag2 value")
  expect_true(is.na(d1[year == 2001, exn_pm_l1]), info = "lag NA at panel start")
  expect_equal(d1[year == 2003, exn_pm_f1], 4, info = "raw lead f1 value")
  expect_equal(d1[year == 2003, exn_pm_d1], 2, info = "d1 = f1 - l1")
  expect_equal(d1[year == 2003, exn_pm_d2], 3, info = "d2 = f2 - l1")
  expect_equal(d1[year == 2003, exn_pm_dm2], -1, info = "dm2 = l2 - l1")
  expect_true(is.na(d1[year == 2005, exn_pm_d1]), info = "d NA at panel end")

  ## artist 2: pi = 10..50
  d2 <- dt[PeopleID == 2]
  expect_equal(d2[year == 2003, exn_pi_dm2], -10, info = "pi dm2")
  expect_equal(d2[year == 2003, exn_pi_d1], 20, info = "pi d1")

  ## n_lag and impact_years are independent: n_lag = 2, impact = 1
  c_ind <- gd_syn_config(dt_h, treatment_cols = c("exn_pm", "exn_pi"),
                         covariate_cols = NULL, years = 2003:2004,
                         label = "indep", n_lag = 2L, impact_years = 1L)
  dt_ind <- md_sbw_panel(c_ind)
  expect_true(all(c("exn_pm_l1", "exn_pm_l2", "exn_pm_dm2") %in% names(dt_ind)),
              info = "lags/dm only up to n_lag (2)")
  expect_true(all(c("exn_pm_f1", "exn_pm_d1") %in% names(dt_ind)),
              info = "leads/d only up to impact_years (1)")
  expect_false(any(grepl("_f2$|_d2$|_dm3$", names(dt_ind))),
               info = "no columns beyond n_lag / impact_years")

  ## clash error when the panel already has pipeline columns
  c_clash <- c_h; c_clash$dt_panel <- copy(dt_h); c_clash$dt_panel[, exn_pm_l1 := 0L]
  expect_error(md_sbw_panel(c_clash),
               pattern = "pipeline-generated", info = "precomputed lag cols rejected")
}

# === covariate lags semantics (numeric time-varying vs time-invariant vs factor) ===
gd_test_covlags <- function() {
  dt <- gd_syn_dt(120, 2006:2012, seed = 4)

  ## numeric time-varying covariate: lags kept
  c_tv <- gd_syn_config(dt, covariate_cols = "ar1", n_lag_cov = 2L)
  l_tv <- gd_prep_cfg(c_tv)
  expect_true(all(c("ar1", "ar1_l1", "ar1_l2") %in% l_tv$cfg$covariate_cols),
              info = "time-varying numeric covariate keeps lags")
  expect_true(all(c("ar1_l1", "ar1_l2") %in% names(l_tv$dt)),
              info = "time-varying numeric lag columns exist in panel")

  ## numeric time-invariant covariate: lags dropped (lag == level value)
  c_ti <- gd_syn_config(dt, covariate_cols = "BirthYear", n_lag_cov = 2L)
  l_ti <- gd_prep_cfg(c_ti)
  expect_equal(l_ti$cfg$covariate_cols, "BirthYear",
               info = "time-invariant numeric covariate: lags dropped")
  expect_false(any(grepl("^BirthYear_l", names(l_ti$dt))),
               info = "time-invariant lag columns removed from panel")

  ## time-invariant factor (Gender): level dummies kept, their lags dropped
  c_g <- gd_syn_config(dt, covariate_cols = "Gender", n_lag_cov = 2L)
  l_g <- gd_prep_cfg(c_g)
  expect_true("Gender_Male" %in% l_g$cfg$covariate_cols, info = "gender dummy kept")
  expect_false(any(grepl("^Gender_Male_l", l_g$cfg$covariate_cols)),
               info = "time-invariant factor dummy lags dropped")

  ## time-varying factor: dummy lags carry information and are kept
  dt_tv <- dt[, .(PeopleID, year)]
  set.seed(11)
  ## status changes for ~30% of artist-years (random walk on 3 levels)
  dt_tv[, status := rbinom(.N, 1, 0.25)]   ## placeholder; rebuilt below per artist
  l_art <- unique(dt_tv$PeopleID)
  dt_status <- map_dfr(l_art, function(id) {
    n <- nrow(dt_tv[PeopleID == id])
    st <- integer(n)
    st[1] <- sample(1:3, 1)
    for (j in 2:n) st[j] <- if (runif(1) < 0.3) sample(1:3, 1) else st[j - 1]
    data.table(PeopleID = id, year = sort(dt_tv[PeopleID == id, year]), status = st)
  })
  dt_tv2 <- merge(dt, dt_status, by = c("PeopleID", "year"))
  dt_tv2[, status := factor(status)]
  c_fv <- gd_syn_config(dt_tv2, covariate_cols = "status", n_lag_cov = 1L,
                        years = 2008:2009)
  l_fv <- gd_prep_cfg(c_fv)
  expect_true(any(grepl("^status_.*_l1$", l_fv$cfg$covariate_cols)),
              info = "time-varying factor dummy lags kept (dummy lag informative)")
  expect_true(all(grepl("^status_", l_fv$cfg$covariate_cols) | l_fv$cfg$covariate_cols == "status"),
              info = "factor replaced by dummy set")

  ## end-to-end fit with a numeric time-varying covariate incl. its lag;
  ## exact balance on a continuous covariate needs sample size, so use a full
  ## pool and a uniform balance tolerance (tols = 1e-3)
  dt_fit_pool <- gd_syn_dt(2500, 2006:2014, seed = 9)
  c_fit <- gd_syn_config(dt_fit_pool, covariate_cols = "ar1", n_lag_cov = 1L,
                         years = 2009:2010, label = "test_covlags", tols = 1e-3)
  r <- gl_sbw_proc(c_fit)
  expect_true(all(r$dt_diag$status == "solved"), info = "covariate-lag fit solves")
  expect_true(all(r$dt_diag$max_cor_w < 0.05), info = "covariate-lag fit balances")
  expect_true(any(grepl("^ar1_l1$", r$config$covariate_cols)),
              info = "covariate lags part of fitted config")
}

# === year eligibility ===
gd_test_years <- function() {
  dt <- gd_syn_dt(300, 2006:2014, seed = 2)  ## smaller pool, all enter 2006
  c_sbw <- gd_syn_config(dt, years = NULL, min_n_year = 50L)
  dt_y <- gd_sbw_years(c_sbw, md_sbw_panel(c_sbw))

  ## lags need years <= 2006 present -> first eligible year is 2008
  ## outcome coverage: t + 2 <= 2014 -> last eligible year is 2012
  expect_equal(dt_y$year, 2008:2012, info = "auto years from coverage + min n")
  expect_true(all(dt_y$n_use == 300), info = "n_use = full pool once lags exist")

  ## explicit years: honored within eligibility, others dropped with warning
  expect_warning(
    dt_y2 <- gd_sbw_years(gd_syn_config(dt, years = c(2009, 2011, 2015)), md_sbw_panel(c_sbw)),
    pattern = "not eligible", info = "warn on ineligible requested years")
  expect_equal(dt_y2$year, c(2009, 2011), info = "explicit eligible years kept")
}


# === fit-population eligibility filter (elig_col / elig_min) ===
gd_test_elig <- function() {
  dt <- copy(gd_syn_dt(600, 2006:2014, seed = 5))
  setorder(dt, PeopleID, year)
  ## monotone "cumulative exhibitions strictly BEFORE t" — the eligibility
  ## use case; the project panel's exn_all_cum has exactly this definition
  dt[, exn_all_cum := cumsum(exn_pm + exn_pi + exn_oes + exn_toti) -
                      (exn_pm + exn_pi + exn_oes + exn_toti), by = PeopleID]
  THR <- 3L

  ## --- validation ---
  expect_error(gd_sbw_check(md_sbw_defaults(gd_syn_config(dt, elig_col = "exn_all_cum"))),
               pattern = "same length", info = "elig_col/elig_min same length")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, elig_col = "not_a_col", elig_min = 1))),
    pattern = "not in dt_panel", info = "elig_col must exist")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, elig_col = "Gender", elig_min = 1))),
    pattern = "numeric", info = "elig_col must be numeric")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, elig_col = "exn_all_cum", elig_min = NA_real_))),
    pattern = "without NAs", info = "elig_min must be complete")

  ## --- mask semantics ---
  c_no <- gd_syn_config(dt, years = 2009:2010, min_n_year = 1L)
  c_el <- gd_syn_config(dt, years = 2009:2010, min_n_year = 1L,
                        elig_col = "exn_all_cum", elig_min = THR)
  cf_no <- gd_prep_cfg(c_no); cf_el <- gd_prep_cfg(c_el)
  b_no <- gd_sbw_fitmask(cf_no$dt, cf_no$cfg)
  b_el <- gd_sbw_fitmask(cf_el$dt, cf_el$cfg)

  ## without elig the mask is exactly the old completeness gate
  l_fitcols <- c(cf_no$cfg$treatment_cols, gl_sbw_balcols(cf_no$cfg),
                 gd_sbw_shiftcols(cf_no$cfg$treatment_cols, cf_no$cfg$impact_years, "f"))
  expect_equal(b_no, complete.cases(cf_no$dt[, .SD, .SDcols = l_fitcols]),
               info = "no elig -> mask == completeness gate")
  ## with elig it is completeness AND (elig_col >= elig_min): elig only removes
  expect_equal(b_el, b_no & cf_el$dt$exn_all_cum >= THR,
               info = "mask == completeness AND elig")
  expect_true(sum(b_el) < sum(b_no), info = "eligibility removes fit rows")

  ## --- KEY property: elig keeps the rows that PRE-FILTERING the panel loses,
  ## because the pre-threshold years stay in the panel as lag donors ---
  dt_pref <- copy(dt)[exn_all_cum >= THR]
  c_pref <- gd_syn_config(dt_pref, years = 2009:2010, min_n_year = 1L,
                          label = "test_elig_pref")
  cf_pref <- gd_prep_cfg(c_pref)
  b_pref <- gd_sbw_fitmask(cf_pref$dt, cf_pref$cfg)
  expect_true(sum(b_pref) < sum(b_el),
              info = "panel pre-filter keeps strictly fewer fit rows than elig")
  ## the pre-filtered fit rows are a subset of the elig fit rows ...
  d_el   <- cf_el$dt[b_el, .(PeopleID, year)]
  d_pref <- cf_pref$dt[b_pref, .(PeopleID, year)]
  expect_equal(nrow(d_pref[!d_el, on = c("PeopleID", "year")]), 0L,
               info = "pre-filtered fit rows subset of elig fit rows")
  ## ... and the extra rows are above the threshold but incomplete after the
  ## pre-filter (NA lags from the truncated panel start)
  d_extra <- d_el[!d_pref, on = c("PeopleID", "year")]
  expect_true(nrow(d_extra) > 0, info = "elig-only rows exist")
  expect_true(all(cf_el$dt[d_extra, on = c("PeopleID", "year"), exn_all_cum] >= THR),
              info = "elig-only rows are above the threshold")
  expect_true(all(is.na(cf_pref$dt[d_extra, on = c("PeopleID", "year"), exn_pm_l2])),
              info = "those rows lose their lags under panel pre-filtering")

  ## --- end-to-end: weights only on eligible rows; elig_min enters the key ---
  r_el <- t_test("elig fit", gl_sbw_proc(c_el))
  expect_true(all(r_el$dt_diag$status == "solved"), info = "eligibility fit solves")
  expect_true(all(r_el$dt_results[!is.na(w_sbw), exn_all_cum] >= THR),
              info = "no weight below the eligibility threshold")
  expect_equal(sum(!is.na(r_el$dt_results$w_sbw)), sum(r_el$dt_years$n_use),
               info = "weighted rows == n_use (elig)")
  ## the DiD runs on the eligible population only
  dt_c <- gd_sbw_coef_did(gl_sbw_did(r_el, l_tps = 1:2, l_dvs = "exn_pm"))
  expect_true(all(is.finite(dt_c$estimate)), info = "DiD on the elig sample is estimable")

  c_el2 <- gd_syn_config(dt, years = 2009:2010, min_n_year = 1L,
                         elig_col = "exn_all_cum", elig_min = THR + 1)
  r_el2 <- gl_sbw_proc(c_el2)
  expect_false(r_el2$cache$dir_run == r_el$cache$dir_run,
               info = "elig_min enters the cache key")
}
# === end-to-end pipeline + caching + registry ===
gd_test_proc <- function() {
  dt <- gd_syn_dt(2500, 2006:2014)
  c_sbw <- gd_syn_config(dt, years = 2009:2010)

  r <- t_test("proc fit", gl_sbw_proc(c_sbw))
  expect_true(all(c("dt_panel", "dt_results", "dt_years", "dt_diag", "config", "cache") %in% names(r)),
              info = "result structure")
  expect_true("w_sbw" %in% names(r$dt_results), info = "weight column added")
  expect_true(all(c("exn_pm_l1", "exn_pm_d1", "exn_pm_dm2") %in% names(r$dt_results)),
              info = "lags and outcome columns present")
  expect_equal(nrow(r$dt_years), 2L, info = "two fit years")
  expect_equal(sum(!is.na(r$dt_results$w_sbw)), sum(r$dt_years$n_use), info = "weighted rows = n_use")
  expect_equal(nrow(r$dt_diag), 2L, info = "diag row per year")
  expect_true(all(r$dt_diag$status == "solved"), info = "all years solved")
  expect_true(all(r$dt_diag$max_cor_w < 0.02), info = "weighted correlations ~ 0")
  expect_true(all(r$dt_diag$ess > 0.1 * r$dt_diag$n_fit), info = "ESS not collapsed")
  expect_true(all(r$dt_diag$max_cor_raw > r$dt_diag$max_cor_w), info = "balance improves over raw")

  ## cache: rerun is a pure reload with identical weights
  r2 <- t_test("proc cached", gl_sbw_proc(c_sbw))
  expect_true(r2$cache$from_cache_all, info = "second run served from cache")
  expect_equal(r$dt_results$w_sbw, r2$dt_results$w_sbw, info = "cached weights identical")
  expect_true(r$cache$dir_run == r2$cache$dir_run, info = "same cache dir")

  ## b_force refits
  r3 <- t_test("proc refit", gl_sbw_proc(c_sbw, b_force = TRUE))
  expect_false(r3$cache$from_cache_all, info = "b_force refits")
  expect_equal(r$dt_results$w_sbw, r3$dt_results$w_sbw, info = "deterministic refit weights")

  ## registry + persisted files
  expect_true(file.exists(file.path(c_sbw$dir_cache, "registry.csv")), info = "registry exists")
  dt_reg <- gd_sbw_registry(c_sbw$dir_cache)
  expect_true(nrow(dt_reg) >= 1, info = "registry has runs")
  expect_true(all(c("label", "hash", "years", "n_years", "n_obs") %in% names(dt_reg)),
              info = "registry columns")
  expect_true(dt_reg$label[1] == "test_syn", info = "registry records label")
  expect_true(file.exists(file.path(r$cache$dir_run, "diag.csv")), info = "diag.csv persisted")
  expect_true(file.exists(file.path(r$cache$dir_run, "config.rds")), info = "config snapshot persisted")
}

# === balance long table ===
gd_test_bal_long <- function() {
  dt <- gd_syn_dt(2500, 2006:2014)
  c_sbw <- gd_syn_config(dt, years = 2009:2010)
  r <- gl_sbw_proc(c_sbw)
  dt_bl <- gd_sbw_bal_long(r$dt_results, r$config)

  ## Gender -> 1 dummy (Gender_Male; Female = dropped reference):
  ## 4 treats x (8 lags + 1 dummy) x 2 years
  expect_equal(nrow(dt_bl), 2 * 4 * 9, info = "long table dimensions")
  expect_true(all(c("year", "treat", "balcol", "cor_raw", "cor_w") %in% names(dt_bl)),
              info = "long table columns")
  expect_true(all(abs(dt_bl$cor_w) < 0.05, na.rm = TRUE), info = "weighted corr small")
  expect_true(any(abs(dt_bl$cor_raw) > 0.1, na.rm = TRUE), info = "raw corr material (data informative)")
  expect_true(any(grepl("^Gender_", dt_bl$balcol)), info = "dummy appears in balance table")
}

# === DiD estimation on the weighted panel ===
gd_test_did <- function() {
  dt <- gd_syn_dt(2500, 2006:2014)
  c_sbw <- gd_syn_config(dt, years = 2009:2010)
  r <- gl_sbw_proc(c_sbw)

  ## formula uses the explicit no-intercept year factor (not | year)
  f_did <- gf_sbw_did(r$config, l_tps = c(-2, 1, 2), l_dvs = "exn_pm")
  expect_true(grepl("0 \\+ factor\\(year\\)", paste(deparse(f_did), collapse = "")),
              info = "formula has 0 + factor(year)")
  expect_false(grepl("\\|", paste(deparse(f_did), collapse = "")),
               info = "formula has no absorbed FE")

  r_did <- t_test("did fit", gl_sbw_did(r))
  expect_true(inherits(r_did, "fixest_multi"), info = "fixest multi-LHS object")
  ## gl_sbw_did attaches the config (with tid) so gd_sbw_coef_did can locate
  ## the year-factor terms without guessing from names
  expect_true(!is.null(attr(r_did, "c_sbw")) && !is.null(attr(r_did, "c_sbw")$tid),
              info = "config attached to the fitted model")

  dt_c <- gd_sbw_coef_did(r_did)
  expect_true(all(c("mdl", "term", "estimate", "std.error", "tp", "dv") %in%
                  names(dt_c)), info = "coef table columns")
  ## per LHS: 4 treatment terms + 1 pooled no_treatment; 4 dvs x 3 tps
  expect_equal(nrow(dt_c), 4 * 3 * 5, info = "4 dvs x 3 tps x (4 treats + no_treatment)")
  expect_true(all(c(-2, 1, 2) %in% unique(dt_c$tp)), info = "pre and post time points")
  ## pre-period (dm2) coefficients ~ 0 once lags are balanced
  expect_true(max(abs(dt_c[tp == -2 & term != "no_treatment", estimate])) < 0.1,
              info = "pre period flat (weighted)")
  ## no year-factor coefficient rows leak into the output
  expect_false(any(grepl("factor\\(", dt_c$term)), info = "no raw factor rows in output")

  ## no_treatment: one row per (dv, tp), finite SE
  dt_nt <- dt_c[term == "no_treatment"]
  expect_equal(nrow(dt_nt), 4 * 3, info = "one no_treatment row per (dv, tp)")
  expect_true(all(is.finite(dt_nt$estimate)), info = "no_treatment estimates finite")
  expect_true(all(is.finite(dt_nt$std.error)), info = "no_treatment SE finite")
  ## explicit config path gives the same table
  dt_c_b <- gd_sbw_coef_did(r_did, c_sbw = r$config)
  expect_equal(as.data.frame(dt_c_b), as.data.frame(dt_c),
               info = "explicit c_sbw equals the attached-config path")

  ## no_treatment equals the equal-weight mean of the year-factor coefficients
  ## (= mean(fixef()) of the absorbed parametrization), and treatment coefs are
  ## unchanged vs the absorbed spec
  l_tr <- r$config$treatment_cols
  f_lhs <- "c(exn_pm_d1, exn_pm_d2)"
  f_tr  <- paste(l_tr, collapse = " + ")
  m_fac <- feols(as.formula(sprintf("%s ~ 0 + factor(year) + %s", f_lhs, f_tr)),
                 r$dt_results, weights = ~w_sbw, cluster = ~PeopleID)
  m_fe  <- feols(as.formula(sprintf("%s ~ %s | year", f_lhs, f_tr)),
                 r$dt_results, weights = ~w_sbw, cluster = ~PeopleID)
  expect_true(all(sapply(seq_along(m_fe), function(i)
    all(abs(coef(m_fe[[i]])[l_tr] - coef(m_fac[[i]])[l_tr]) < 1e-8))),
    info = "treatment coefs identical to absorbed | year spec")
  for (i in 1:2) {
    lvl_mean <- as.numeric(mean(coef(m_fac[[i]])[grep("^factor", names(coef(m_fac[[i]])))]))
    expect_equal(as.numeric(dt_nt[dv == "exn_pm" & tp == i, estimate]), lvl_mean,
                 info = sprintf("no_treatment == mean year level (d%d)", i))
  }

  ## single-year config: no_treatment still identified (factor has one level)
  c_1y <- gd_syn_config(dt, years = 2010, label = "test_did_1y", n_lag = 2L, impact_years = 2L)
  r_1y <- gl_sbw_proc(c_1y)
  dt_c1 <- gd_sbw_coef_did(gl_sbw_did(r_1y, l_tps = 1:2, l_dvs = "exn_pm"))
  expect_true(all(is.finite(dt_c1[term == "no_treatment", std.error])),
              info = "single-year no_treatment SE finite")

  ## robustness: another factor in the model must not be taken for the year
  ## factor (gd_sbw_notreat_one matches the exact "factor(<tid>)" prefix)
  dt_r <- copy(r$dt_results)
  dt_r[, g_fac := factor(PeopleID %% 3)]
  m_twofac <- feols(exn_pm_d1 ~ 0 + factor(year) + g_fac + exn_pm + exn_pi + exn_oes + exn_toti,
                    dt_r, weights = ~w_sbw, cluster = ~PeopleID)
  b <- coef(m_twofac)
  l_year_only <- names(b)[startsWith(names(b), "factor(year)")]
  expect_true(all(startsWith(l_year_only, "factor(year)")), info = "year terms identified")
  expect_false(any(startsWith(l_year_only, "factor(g_fac)")), info = "g_fac not matched as year")
  r_nt <- gd_sbw_notreat_one(m_twofac, tid = "year")
  expect_equal(r_nt$estimate, as.numeric(mean(b[startsWith(names(b), "factor(year)")])),
               info = "notreat uses only the year-level coefficients")
  e_no <- tryCatch(gd_sbw_coef_did(structure(list(m_twofac), names = "lhs: exn_pm_d1",
                                             class = "fixest_multi")),
                   error = function(e) conditionMessage(e))
  expect_true(is.character(e_no) && grepl("no config", e_no),
              info = "informative error when no config is attached")
}

# === small real-data smoke (subset of the foster panel) ===
gd_test_real <- function() {
  dt_f <- fread("~/Dropbox/phd/papers/foster/misc/dt_atyear_foster.csv",
                select = c("PeopleID", "year", "Gender", "exn_pm", "exn_pi", "exn_oes",
                           "auc_turnover"))
  set.seed(3)
  l_ids <- sample(unique(dt_f$PeopleID), 1500)
  dt_fs <- dt_f[PeopleID %in% l_ids & year <= 2015]
  ## raw Gender has ~20 rare/multi-label levels (exact balance infeasible);
  ## recode to the 3-level project convention before balancing
  dt_fs[, Gender := fifelse(Gender %in% c("Male", "Female"), Gender, "Other")]

  c_real <- list(
    label = "real_smoke", dt_panel = dt_fs, uid = "PeopleID", tid = "year",
    treatment_cols = c("exn_pm", "exn_pi", "exn_oes"),
    covariate_cols = "Gender",          ## factor covariate: exercises the dummy expansion
    outcome_cols = "auc_turnover",      ## real additional outcome (auction turnover)
    n_lag = 2L, impact_years = 2L,
    years = 2012:2013, min_n_year = 100L,
    dir_cache = DIR_CACHE_TEST)

  r <- t_test("real fit", gl_sbw_proc(c_real))
  expect_equal(nrow(r$dt_years), 2L, info = "two real fit years")
  expect_true(all(r$dt_diag$status == "solved"), info = "real years solved")
  expect_true(all(r$dt_diag$n_fit >= 100), info = "real n per year")
  expect_true(all(r$dt_diag$max_cor_w < 0.05), info = "real balance ok")
  expect_true(all(r$dt_diag$ess > 0.3 * r$dt_diag$n_fit), info = "real ESS healthy")
  expect_true("w_sbw" %in% names(r$dt_results), info = "real weights attached")

  ## the real outcome got its shift columns and is estimable via gl_sbw_did
  expect_true(all(c("auc_turnover_f1", "auc_turnover_d1", "auc_turnover_dm2") %in%
                  names(r$dt_results)), info = "real outcome shift columns present")
  r_did_o <- t_test("real did outcome", gl_sbw_did(r, l_dvs = "auc_turnover"))
  dt_co <- gd_sbw_coef_did(r_did_o)
  expect_equal(unique(dt_co$dv), "auc_turnover", info = "real outcome-only model")
  expect_true(any(is.finite(dt_co$estimate)), info = "real outcome coefficients estimated")
}

# === cached weight files carry (uid, tid) row identity + PSOCK parallel fit ===
gd_test_cache_format_and_parallel <- function() {
  dt <- gd_syn_dt(2500, 2006:2014)
  c_sbw <- gd_syn_config(dt, years = 2009:2010)

  ## sequential fit; inspect the on-disk weight format
  r <- gl_sbw_proc(c_sbw)
  p_w <- file.path(r$cache$dir_run, "w_2009.rds")
  dt_w <- readRDS(p_w)
  expect_true(is.data.table(dt_w), info = "weight file is a data.table")
  expect_equal(names(dt_w), c("PeopleID", "year", "w"), info = "weight file columns (uid, tid, w)")
  expect_equal(nrow(dt_w), sum(r$dt_years[year == 2009, n_use]), info = "one weight per complete row")
  ## weights match the panel rows on (uid, tid), not on order: shuffle the
  ## panel's 2009 rows and re-attach from the same file
  dt_shuf <- r$dt_results[sample(.N)]
  ## re-slice + read into the shuffled table
  c_sbw_fit <- r$config   ## config after covariate expansion (Gender -> Gender_Male)
  dt_year <- gd_sbw_year_slice(dt_shuf, 2009, c_sbw_fit)
  expect_equal(nrow(dt_year), nrow(dt_w), info = "shuffled slice rows match weight rows")
  ## values attached by join must equal the file weights for those (uid, tid)
  dt_check <- merge(dt_shuf[year == 2009, .(PeopleID, year)],
                    dt_w, by = c("PeopleID", "year"))
  expect_equal(sum(!is.na(dt_check$w)), nrow(dt_w), info = "every file row joins to a panel row")

  ## fork-parallel (mclapply) fit produces the same weights as the sequential fit
  c_par <- gd_syn_config(dt, years = 2009:2010, label = "test_syn_par", n_cores = 2L)
  r_par <- gl_sbw_proc(c_par)
  expect_true(all(r_par$dt_diag$status == "solved"), info = "parallel years solved")
  expect_true(all(r_par$dt_diag$ess > 0.1 * r_par$dt_diag$n_fit), info = "parallel ESS healthy")
  dt_m <- merge(r$dt_results[!is.na(w_sbw), .(PeopleID, year, w_sbw)],
                r_par$dt_results[!is.na(w_sbw), .(PeopleID, year, w_sbw)],
                by = c("PeopleID", "year"), suffixes = c("_seq", "_par"))
  expect_equal(nrow(dt_m), sum(r$dt_years$n_use), info = "sequential and parallel weight sets align")
  expect_true(max(abs(dt_m$w_sbw_seq - dt_m$w_sbw_par)) < 1e-6, info = "sequential == parallel weights")
}

# === additional outcome variables (outcome_cols) ===
gd_test_outcome <- function() {
  dt <- gd_syn_dt(2500, 2006:2014)

  ## config validation
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols = "not_a_col"))),
    pattern = "outcome_cols not in dt_panel", info = "outcome col must exist")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols = "exn_pm"))),
    pattern = "overlap treatment_cols", info = "outcome must not be a treatment")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols = "Gender"))),
    pattern = "outcome_cols must be numeric", info = "outcome must be numeric")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, covariate_cols = "ar1", outcome_cols = "ar1"))),
    pattern = "overlap covariate_cols", info = "outcome must not be a covariate")

  ## hand-checked shift columns for an outcome
  dt_h <- data.table(PeopleID = rep(1:2, each = 6), year = rep(2001:2006, 2),
                     exn_pm = c(1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0),
                     auc = c(10, 12, 15, 19, 24, 30, 100, 90, 81, 73, 66, 60))
  c_h <- gd_syn_config(dt_h, treatment_cols = "exn_pm", outcome_cols = "auc",
                       covariate_cols = NULL, years = 2003:2004, label = "out_hand",
                       n_lag = 2L, impact_years = 2L)
  dt_p <- md_sbw_panel(c_h)
  a1 <- dt_p[PeopleID == 1]
  expect_equal(a1[year == 2003, auc_l1], 12, info = "outcome lag1")
  expect_equal(a1[year == 2003, auc_l2], 10, info = "outcome lag2")
  expect_true(all(c("auc_f1", "auc_f2") %in% names(dt_p)), info = "outcome leads kept")
  expect_equal(a1[year == 2003, auc_f1], 19, info = "outcome lead f1")
  expect_equal(a1[year == 2003, auc_f2], 24, info = "outcome lead f2")
  expect_equal(a1[year == 2003, auc_d1], 7, info = "outcome d1 = f1 - l1")
  expect_equal(a1[year == 2003, auc_d2], 12, info = "outcome d2 = f2 - l1")
  expect_equal(a1[year == 2003, auc_dm2], -2, info = "outcome dm2 = l2 - l1")

  ## weights are invariant to adding an outcome (no cache invalidation)
  ## same label: outcome difference must not change the hash/cache dir
  c0 <- gd_syn_config(dt, years = 2009:2010, label = "test_out")
  c1 <- gd_syn_config(dt, years = 2009:2010, label = "test_out",
                      outcome_cols = "auc_turnover")
  r0 <- gl_sbw_proc(c0)
  r1 <- gl_sbw_proc(c1)
  expect_true(r1$cache$dir_run == r0$cache$dir_run,
              info = "same cache dir with and without outcome (hash unaffected)")
  expect_equal(r0$dt_results$w_sbw, r1$dt_results$w_sbw,
               info = "weights identical with/without outcome")
  expect_equal(r0$dt_years$n_use, r1$dt_years$n_use,
               info = "fit population unchanged by outcome")

  ## outcome in the dt_results + estimable via gl_sbw_did
  expect_true(all(c("auc_turnover_f1", "auc_turnover_d1", "auc_turnover_dm2") %in%
                  names(r1$dt_results)), info = "outcome shift columns in results")
  r_did_all <- gl_sbw_did(r1)      ## default: treatments + outcome dvs
  dt_c <- gd_sbw_coef_did(r_did_all)
  expect_true("auc_turnover" %in% unique(dt_c$dv), info = "outcome dv in default model")
  ## per LHS: 4 treatment terms + 1 no_treatment; 5 dvs x 3 tps
  expect_equal(nrow(dt_c), 5 * 3 * 5, info = "5 dvs x 3 tps x (4 treats + no_treatment)")
  r_did_o <- gl_sbw_did(r1, l_dvs = "auc_turnover")
  dt_co <- gd_sbw_coef_did(r_did_o)
  expect_equal(unique(dt_co$dv), "auc_turnover", info = "subset model to outcome only")
  expect_equal(nrow(dt_co), 1 * 3 * 5, info = "outcome-only model dimensions")
  expect_true(all(is.finite(dt_co[term == "no_treatment", std.error])),
              info = "outcome no_treatment SE finite")
  expect_error(gl_sbw_did(r1, l_dvs = "not_a_dv"), pattern = "unknown dv",
               info = "unknown dv errors")
}

# === cumulative outcomes (outcome_cols_cur: DiD base = t) ===
gd_test_outcome_cur <- function() {
  dt <- copy(gd_syn_dt(2500, 2006:2014))
  ## a genuine CUMULATIVE column built from a treatment's own flow (year t
  ## INCLUDED): exn_pi_cum2(t) = sum_{s<=t} exn_pi(s)
  setorder(dt, PeopleID, year)
  dt[, exn_pi_cum2 := cumsum(exn_pi), by = PeopleID]

  ## --- validation ---
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols_cur = "not_a_col"))),
    pattern = "outcome_cols_cur not in dt_panel", info = "cur outcome must exist")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols_cur = "exn_pm"))),
    pattern = "overlap treatment_cols", info = "cur outcome must not be a treatment")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols_cur = "Gender"))),
    pattern = "outcome_cols_cur must be numeric", info = "cur outcome must be numeric")
  expect_error(gd_sbw_check(md_sbw_defaults(
    gd_syn_config(dt, outcome_cols = "exn_pi_cum2", outcome_cols_cur = "exn_pi_cum2"))),
    pattern = "both outcome_cols", info = "no column in both outcome families")

  ## --- shift semantics: base cur anchors the window at t (excludes year t) ---
  c_l1  <- gd_syn_config(dt, outcome_cols = "exn_pi_cum2", label = "test_curr")
  c_cur <- gd_syn_config(dt, outcome_cols_cur = "exn_pi_cum2", label = "test_curr")
  dt_l1  <- md_sbw_panel(c_l1)
  dt_cur <- md_sbw_panel(c_cur)

  expect_equal(sort(names(dt_l1)), sort(names(dt_cur)),
               info = "identical generated columns for both outcome families")

  ## l1 base: d1 = c(t+1) - c(t-1) = v(t) + v(t+1)
  ## cur base: d1 = c(t+1) - c(t)   = v(t+1)
  ## -> they differ by exactly the year-t flow v(t) = exn_pi
  d_cmp <- na.omit(
    merge(dt_cur[, .(PeopleID, year, cur = exn_pi_cum2_d1, flow_t = exn_pi)],
          dt_l1[, .(PeopleID, year, l1 = exn_pi_cum2_d1)],
          by = c("PeopleID", "year")))
  expect_true(all(abs(d_cmp$l1 - d_cmp$cur - d_cmp$flow_t) < 1e-9),
              info = "l1 base minus cur base = the year-t flow")

  ## cur d2 = sum of the flow over t+1..t+2 (year t excluded), hand-checked
  f_sum <- sum(dt[PeopleID == 1 & year %in% 2010:2011, exn_pi])
  expect_equal(dt_cur[PeopleID == 1 & year == 2009, exn_pi_cum2_d2], f_sum,
               info = "cur d2 = flow summed over t+1..t+2")
  ## dm{k} identical across families (the pre-period already excludes year t)
  expect_equal(dt_l1$exn_pi_cum2_dm2, dt_cur$exn_pi_cum2_dm2,
               info = "dm unchanged by the outcome base")

  ## --- end-to-end: base l1 adds exactly +1 to the OWN coefficient ---
  ## weights depend on treatments/covariates only -> identical across the two
  ## configs (the outcome set is excluded from the cache key)
  r_l1  <- t_test("cur l1 base", gl_sbw_proc(c_l1))
  r_cur <- t_test("cur base", gl_sbw_proc(c_cur))
  expect_true(r_cur$cache$dir_run == r_l1$cache$dir_run,
              info = "outcome family does not change the cache dir")
  expect_equal(r_l1$dt_results$w_sbw, r_cur$dt_results$w_sbw,
               info = "identical weights across outcome families")

  dt_c1 <- gd_sbw_coef_did(gl_sbw_did(r_l1,  l_dvs = "exn_pi_cum2"))
  dt_cc <- gd_sbw_coef_did(gl_sbw_did(r_cur, l_dvs = "exn_pi_cum2"))
  d_m <- merge(dt_c1[, .(dv, tp, term, est_l1 = estimate, se_l1 = std.error)],
               dt_cc[, .(dv, tp, term, est_cur = estimate, se_cur = std.error)],
               by = c("dv", "tp", "term"))
  expect_equal(nrow(d_m), nrow(dt_c1), info = "coef rows align across families")

  ## post rows: the year-t exposure sits inside the l1-base window and is
  ## itself a regressor, so the own coefficient is exactly +1 higher
  expect_true(all(abs(d_m[term == "exn_pi" & tp > 0, est_l1 - est_cur - 1]) < 1e-6),
              info = "cur base removes the mechanical +1 on the own coefficient")
  ## pre (dm) rows are identical: the pre-window already excludes year t
  expect_true(all(abs(d_m[term == "exn_pi" & tp < 0, est_l1 - est_cur]) < 1e-6),
              info = "pre-period own coefficient unchanged by the base")
  ## every other coefficient (incl. no_treatment) is unaffected
  expect_true(all(abs(d_m[term != "exn_pi", est_l1 - est_cur]) < 1e-6),
              info = "all other coefficients unaffected by the base")
  ## SEs identical everywhere (the -1 is exact -> residuals unchanged)
  expect_true(max(abs(d_m$se_l1 - d_m$se_cur)) < 1e-8,
              info = "standard errors unaffected by the base")
}

# === run everything ===
cat("\n--- config checks ---\n");       gd_test_check()
cat("\n--- formulas ---\n");            gd_test_formulas()
cat("\n--- covariate expansion ---\n"); gd_test_covexp()
cat("\n--- panel builder ---\n");       gd_test_panel()
cat("\n--- covariate lags ---\n");     gd_test_covlags()
cat("\n--- year eligibility ---\n");    gd_test_years()
cat("\n--- fit-population eligibility ---\n"); gd_test_elig()
cat("\n--- end-to-end pipeline ---\n"); gd_test_proc()
cat("\n--- balance long ---\n");        gd_test_bal_long()
cat("\n--- DiD ---\n");                 gd_test_did()
cat("\n--- outcomes ---\n");            gd_test_outcome()
cat("\n--- outcomes (cur base) ---\n"); gd_test_outcome_cur()
cat("\n--- cache format + parallel ---\n"); gd_test_cache_format_and_parallel()
cat("\n--- real-data smoke ---\n");     gd_test_real()

cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
