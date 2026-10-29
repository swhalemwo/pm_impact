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
source("~/Dropbox/phd/papers/foster/scripts/aucexp.R")
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
               ## auction turnover: noisy and mostly zero, like the real one.
               ## log1p-scale columns follow the project convention since
               ## 2026-10 (_ihs balance/outcome columns were replaced by _log)
               auc_turnover = rpois(n, 0.15 + 0.5 * pm + 0.3 * pi) * 1e4L,
               auc_turnover_log = log1p(rpois(n, 0.15 + 0.5 * pm + 0.3 * pi) * 1e4),
               aucturn_cum2_log = log1p(cumsum(rpois(n, 0.15 + 0.5 * pm + 0.3 * pi) * 1e4)),
               exn_all_cum_log = log1p(cumsum(pm + pi + oes + toti)))
  })
}

## mimic gl_sbw_proc's preparation chain on a config: panel + factor
## expansion + covariate lags + balance lags; returns the prepared config
## (ready for gf_sbw_formulas / gd_sbw_years)
gd_prep_cfg <- function(c_sbw) {
  l_prep <- md_sbw_prepare(c_sbw)
  list(cfg = l_prep$c_sbw, dt = l_prep$dt_results)
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

## --- fixtures for the aucexp tests --------------------------------
## a minimal SBW config: only the pieces md_auc_sel/gc_auc_euro_cfg read
C_CFG_AUC <- list(uid = "PeopleID",
                  outcome_cols = character(0),
                  outcome_cols_cur = c("x_log", "y_log"))

## a real (tiny) fixest_multi with the lhs names gd_sbw_coef_did expects,
## so the nobs cross-check and the coefficient extraction run on real
## objects rather than on a hand-made list
gd_fake_did <- function(n_obs = 12L) {
  set.seed(4)
  ## exn_pm must cross the year dummies, else fixest drops it as collinear
  dt <- data.table(year = rep(2010:2012, length.out = n_obs),
                   exn_pm = rep(c(0, 1, 2, 0, 2, 1), length.out = n_obs),
                   x_log_d2 = rnorm(n_obs), x_log_dm2 = rnorm(n_obs))
  rx <- feols(c(x_log_d2, x_log_dm2) ~ 0 + factor(year) + exn_pm, dt)
  attr(rx, "c_sbw") <- list(tid = "year")
  rx
}

## hand-built panel + coefficient table for the aucexp unit tests: no SBW,
## no model, so every number is checkable by hand. The dv carries a "_log"
## suffix because gd_auc_lvl_tr() infers the back-transform from it.
gd_syn_auc_dt <- function(n_artists = 10, n_years = 5, seed = 11) {
  set.seed(seed)
  dt <- data.table(PeopleID = rep(seq_len(n_artists), each = n_years),
                   year = rep(seq(2010, by = 1, length.out = n_years), n_artists))
  ## zero-inflated turnover: ~1/3 zeros, rest log-normal, so the weighted
  ## mean level is dominated by a few rows (as in the real data)
  dt[, lvl := fifelse(runif(.N) < 1/3, 0, round(exp(rnorm(.N, 8, 2))))]
  dt[, lvl_log := log1p(lvl)]
  setorder(dt, PeopleID, year)
  ## the panel columns the conversion reads, built by hand
  dt[, x_log := lvl_log]                                  ## level at t
  dt[, x_log_f2 := shift(lvl_log, -2), by = PeopleID]     ## level at t+2
  dt[, x_log_d2 := x_log_f2 - lvl_log]                    ## d2 = v(t+2) - v(t)
  dt[, x_log_l3 := shift(lvl_log, 3), by = PeopleID]      ## level at t-3
  dt[, x_log_dm3 := x_log_l3 - shift(lvl_log, 1), by = PeopleID]
  dt[, x_log_l1 := shift(lvl_log, 1), by = PeopleID]
  dt[, w_sbw := runif(.N, 0.5, 2)]
  dt
}

gd_syn_auc_cfg <- function(dt) {
  list(dt = dt, wcol = "w_sbw", uid = "PeopleID", b_ref = "t+k")
}

## one post-period cell per treatment plus one pre-period (placebo) cell
gd_syn_auc_coef <- function() {
  data.table(dv = "x_log", tp = c(2L, 2L, -3L),
             term = c("exn_pm", "exn_pi", "exn_pm"),
             estimate = c(0.2, -0.1, 0.05),
             std.error = c(0.05, 0.04, 0.03))
}

## built once and reused: a real fixest_multi with lhs names
GD_FAKE_DID_AUC <- NULL
gd_fake_did_get <- function() {
  if (is.null(GD_FAKE_DID_AUC)) GD_FAKE_DID_AUC <<- gd_fake_did()
  GD_FAKE_DID_AUC
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
  c_num <- gd_syn_config(dt, covariate_cols = "exn_all_cum_log")
  dt_p2 <- md_sbw_panel(c_num)
  l_num <- md_sbw_expand_covs(dt_p2, "exn_all_cum_log")
  expect_equal(l_num$covs, "exn_all_cum_log", info = "numeric covariate kept as-is")
  expect_true("exn_all_cum_log" %in% names(l_num$dt), info = "numeric covariate retained")
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

# === extra balancing columns (balance_lags: lags of raw panel columns) ===
gd_test_bal_lags <- function() {
  ## --- spec canonicalization ---
  expect_equal(md_sbw_bal_lags_spec(list(auc = c(2L, 1L, 2L), b = 0L)),
               list(auc = 1:2, b = 0L),
               info = "named list: sorted, deduplicated, integer")
  expect_equal(md_sbw_bal_lags_spec(c("auc_l1", "auc_l3", "b_l2")),
               list(auc = c(1L, 3L), b = 2L),
               info = "character shorthand '{col}_l{k}' parsed")
  expect_equal(md_sbw_bal_lags_spec(NULL), list(), info = "NULL/empty spec")
  expect_equal(md_sbw_bal_lags_spec(list()), list(), info = "empty list spec")
  ## idempotent (md_sbw_defaults runs more than once per pipeline call)
  expect_equal(md_sbw_bal_lags_spec(md_sbw_bal_lags_spec(c("auc_l1"))),
               list(auc = 1L), info = "spec canonicalization is idempotent")
  expect_error(md_sbw_bal_lags_spec(c("auc", "b_l1")), pattern = "must look like",
               info = "unparseable shorthand errors")
  expect_error(md_sbw_bal_lags_spec(list(1:2)), pattern = "NAMED list",
               info = "unnamed list rejected")
  expect_error(md_sbw_bal_lags_spec(list(auc = -1L)), pattern = "non-negative",
               info = "negative distance rejected")
  expect_error(md_sbw_bal_lags_spec(list(auc = 1.5)), pattern = "whole numbers",
               info = "fractional distance rejected")

  ## --- config validation ---
  dt0 <- gd_syn_dt(50, 2006:2009, seed = 1)
  c_ok <- gd_syn_config(dt0, years = NULL, min_n_year = 5L)
  expect_true(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(balance_lags = list(ar1 = 1:2))))),
    info = "valid balance_lags passes check")
  expect_error(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(balance_lags = list(not_a_col = 1L))))),
    pattern = "not in dt_panel", info = "balance base must exist")
  expect_error(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(balance_lags = list(Gender = 1L))))),
    pattern = "must be numeric", info = "balance base must be numeric")
  expect_error(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(balance_lags = list(exn_pm = 1L))))),
    pattern = "treatment columns", info = "treatment base rejected (n_lag covers it)")
  expect_error(gd_sbw_check(md_sbw_defaults(
    utils::modifyList(c_ok, list(balance_cols = "nope")))),
    pattern = "neither dt_panel columns nor generated", info = "unknown balance_cols")

  ## --- generation: hand-checked lag values, per-artist boundary ---
  dt_h <- data.table(PeopleID = rep(1:2, each = 5), year = rep(2001:2005, 2),
                     exn_pm = c(1, 2, 3, 4, 5, 0, 0, 0, 0, 0),
                     auc = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500))
  c_h <- gd_syn_config(dt_h, treatment_cols = "exn_pm", covariate_cols = NULL,
                       years = 2004, label = "bal_hand", n_lag = 2L,
                       impact_years = 1L,
                       balance_lags = list(auc = c(0L, 1L, 2L)))
  l_h <- gd_prep_cfg(c_h)
  dt_hp <- l_h$dt
  expect_equal(l_h$cfg$balance_cols, c("auc", "auc_l1", "auc_l2"),
               info = "balance_cols resolved (level + lags)")
  a1 <- dt_hp[PeopleID == 1]; a2 <- dt_hp[PeopleID == 2]
  expect_equal(a1[year == 2003, auc_l1], 20, info = "balance lag1 value")
  expect_equal(a1[year == 2003, auc_l2], 10, info = "balance lag2 value")
  expect_true(is.na(a1[year == 2001, auc_l1]), info = "balance lag NA at panel start")
  expect_equal(a2[year == 2003, auc_l1], 200, info = "lags do not cross artists")
  expect_equal(a2[year == 2003, auc_l2], 100, info = "lag2 does not cross artists")
  ## the level (k = 0) is the base column itself, not a copy
  expect_true("auc" %in% l_h$cfg$balance_cols && "auc" %in% names(dt_hp),
              info = "k = 0 balances the base column itself")
  ## balance columns are NOT treatments: no formula, no d/dm columns
  expect_equal(names(gf_sbw_formulas(l_h$cfg)), "exn_pm",
               info = "no extra optweightMV formula per balance column")
  expect_false(any(grepl("^auc_d", names(dt_hp))), info = "no d{k} columns for balance vars")
  ## and they appear in the formula RHS exactly once
  rhs_h <- as.character(gf_sbw_formulas(l_h$cfg)[[1]])[3]
  expect_true(grepl("auc_l1", rhs_h) && grepl("auc_l2", rhs_h),
              info = "balance lags in the formula RHS")

  ## --- fit mask: a missing balance value drops that row from the fit ---
  dt_h2 <- copy(dt_h)
  dt_h2[PeopleID == 2 & year == 2003, auc := NA_real_]
  c_m0 <- gd_syn_config(dt_h2, treatment_cols = "exn_pm", covariate_cols = NULL,
                        years = 2001:2004, label = "bal_mask", n_lag = 2L,
                        impact_years = 1L)
  c_m1 <- gd_syn_config(dt_h2, treatment_cols = "exn_pm", covariate_cols = NULL,
                        years = 2001:2004, label = "bal_mask", n_lag = 2L,
                        impact_years = 1L, balance_lags = list(auc = 0L))
  l_m0 <- gd_prep_cfg(c_m0); l_m1 <- gd_prep_cfg(c_m1)
  b_m0 <- gd_sbw_fitmask(l_m0$dt, l_m0$cfg)
  b_m1 <- gd_sbw_fitmask(l_m1$dt, l_m1$cfg)
  l_fitcols <- c(l_m1$cfg$treatment_cols, gl_sbw_balcols(l_m1$cfg),
                 gd_sbw_shiftcols(l_m1$cfg$treatment_cols, l_m1$cfg$impact_years, "f"))
  expect_equal(b_m1, complete.cases(l_m1$dt[, .SD, .SDcols = l_fitcols]),
               info = "mask == completeness over balance cols + leads")
  expect_true(all(b_m1 <= b_m0), info = "balance columns only remove fit rows")
  expect_equal(sum(b_m0) - sum(b_m1), 1L,
               info = "exactly the one row with a missing balance value is dropped")
  idx_bad <- l_m1$dt[, which(PeopleID == 2 & year == 2003)]
  expect_true(length(idx_bad) == 1L && !b_m1[idx_bad],
              info = "the incomplete balance row is not fitted")
  ## a lag that is NA only where the treatment lags are NA already costs nothing
  c_m2 <- gd_syn_config(dt_h, treatment_cols = "exn_pm", covariate_cols = NULL,
                        years = 2003:2004, label = "bal_mask2", n_lag = 2L,
                        impact_years = 1L, balance_lags = list(auc = 2L))
  l_m2 <- gd_prep_cfg(c_m2)
  expect_equal(sum(gd_sbw_fitmask(l_m2$dt, l_m2$cfg)),
               sum(gd_sbw_fitmask(l_m0$dt, l_m0$cfg)),
               info = "balance lag NAs only at panel start cost no fit rows")

  ## --- reuse of prebuilt outcome lags (base is an outcome_col) ---
  c_o <- gd_syn_config(dt_h, treatment_cols = "exn_pm", covariate_cols = NULL,
                       years = 2004, label = "bal_reuse", n_lag = 2L,
                       impact_years = 1L, outcome_cols = "auc",
                       balance_lags = list(auc = 1:2))
  l_o <- gd_prep_cfg(c_o)
  expect_equal(l_o$cfg$balance_cols, c("auc_l1", "auc_l2"),
               info = "outcome lags reused as balance columns")
  expect_equal(l_o$dt$auc_l1, l_h$dt$auc_l1, info = "reused lag values unchanged")
  ## a lag beyond n_lag is generated on demand
  c_o2 <- gd_syn_config(dt_h, treatment_cols = "exn_pm", covariate_cols = NULL,
                        years = 2004, label = "bal_reuse2", n_lag = 1L,
                        impact_years = 1L, outcome_cols = "auc",
                        balance_lags = list(auc = 3L))
  l_o2 <- gd_prep_cfg(c_o2)
  expect_true("auc_l3" %in% names(l_o2$dt), info = "balance lag beyond n_lag built")
  expect_equal(l_o2$dt[PeopleID == 1 & year == 2004, auc_l3], 10,
               info = "generated on-demand lag value")

  ## --- end-to-end: balancing on past auction performance ---
  ## a cumulative, zero-inflated 'auction' variable (like auc_turnover) that is
  ## independent of the treatment histories
  dt_fit <- copy(gd_syn_dt(2500, 2006:2014, seed = 8))
  setorder(dt_fit, PeopleID, year)
  set.seed(31)
  dt_fit[, auc_lvl := asinh(cumsum(rpois(.N, 0.4))), by = PeopleID]

  c_0 <- gd_syn_config(dt_fit, years = 2009:2010, label = "test_bal_fit", tols = 1e-3)
  c_1 <- gd_syn_config(dt_fit, years = 2009:2010, label = "test_bal_fit", tols = 1e-3,
                       balance_lags = list(auc_lvl = 1:2))
  r_0 <- t_test("bal_lags off", gl_sbw_proc(c_0))
  r_1 <- t_test("bal_lags on", gl_sbw_proc(c_1))
  expect_true(all(r_1$dt_diag$status == "solved"), info = "fit with balance_lags solves")
  expect_equal(r_1$config$balance_cols, c("auc_lvl_l1", "auc_lvl_l2"),
               info = "resolved balance_cols in the fitted config")
  expect_true(all(c("auc_lvl_l1", "auc_lvl_l2") %in% names(r_1$dt_results)),
              info = "balance lag columns present in dt_results")
  expect_true(all(c("auc_lvl_l1", "auc_lvl_l2") %in%
                  gd_sbw_bal_long(r_1$dt_results, r_1$config)$balcol),
              info = "balance lags reported in the balance table")
  dt_bl <- gd_sbw_bal_long(r_1$dt_results, r_1$config)
  expect_true(all(abs(dt_bl[grepl("^auc_lvl", balcol), cor_w]) < 0.02),
              info = "weighted corr with past auction performance ~ 0")
  expect_true(any(abs(dt_bl[grepl("^auc_lvl", balcol), cor_raw]) > 0.02),
              info = "raw corr with past auction performance is material")
  ## the extra constraints change the weights, the estimand and the cache key
  expect_false(r_1$cache$dir_run == r_0$cache$dir_run,
               info = "balance_lags enters the cache key")
  expect_true(max(abs(r_1$dt_results$w_sbw - r_0$dt_results$w_sbw), na.rm = TRUE) > 1e-6,
              info = "balance_lags change the weights")
  expect_equal(r_1$dt_years$n_use, r_0$dt_years$n_use,
               info = "no rows lost when the balance lags are observed everywhere")
  ## DiD still estimable on the reweighted panel
  dt_c <- gd_sbw_coef_did(gl_sbw_did(r_1, l_dvs = "exn_pm"))
  expect_true(all(is.finite(dt_c$estimate)), info = "DiD estimable with balance_lags")
  ## the balance columns are RHS-only: they never become DVs
  expect_false("auc_lvl" %in% unique(dt_c$dv), info = "balance vars are not outcomes")

  ## --- balance_cols escape hatch: a precomputed column balanced verbatim ---
  dt_fit[, auc_bin := as.integer(auc_lvl > 0)]
  c_verb <- gd_syn_config(dt_fit, years = 2009:2010, label = "test_bal_verb",
                          tols = 1e-3, balance_cols = "auc_bin")
  r_v <- gl_sbw_proc(c_verb)
  expect_equal(r_v$config$balance_cols, "auc_bin",
               info = "manually supplied balance_cols kept")
  expect_true(any(gd_sbw_bal_long(r_v$dt_results, r_v$config)$balcol == "auc_bin"),
              info = "verbatim balance column balanced")
}

# === column-name resolution: level column, outcome column, back-transform ===
gd_test_auc_names <- function() {
  ## the reference level of a POST-period row lives at t+k, that of a
  ## PRE-period row at t-k: dm{k} = v(t-k) - v(t-1) runs backwards
  expect_equal(gd_auc_lvl_col("x_log", 2L, "t+k"), "x_log_f2",
               info = "post period level is {dv}_f{k}")
  expect_equal(gd_auc_lvl_col("x_log", -3L, "t+k"), "x_log_l3",
               info = "pre period level is {dv}_l{k}")
  expect_equal(gd_auc_lvl_col("x_log", 4L, "t-1"), "x_log_l1",
               info = "t-1 level is {dv}_l1 whatever tp is")
  expect_equal(gd_auc_lvl_col("x_log", -2L, "t"), "x_log",
               info = "t level is the bare dv")
  ## an unknown ref yields NULL; b_ref is validated by match.arg in
  ## gc_auc_euro_cfg, so a bad value cannot reach here through the API
  expect_null(gd_auc_lvl_col("x_log", 2L, 500),
              info = "unknown ref has no level column")
  ## the outcome column is d{k} for post, dm{k} for pre
  expect_equal(gd_auc_ref_cols("x_log", 2L, "t+k"),
               c(lvl_col = "x_log_f2", d_col = "x_log_d2"),
               info = "post ref cols")
  expect_equal(gd_auc_ref_cols("x_log", -3L, "t+k"),
               c(lvl_col = "x_log_l3", d_col = "x_log_dm3"),
               info = "pre ref cols use dm{k}")


  ## back-transform inferred from the dv suffix, and it is a round trip
  vec_raw <- c(0, 1, 5, 1e6)
  expect_equal(gd_auc_lvl_tr("auc_log", log1p(vec_raw)), vec_raw,
               info = "_log -> expm1 round trip")
  expect_equal(gd_auc_lvl_tr("auc_ihs", asinh(vec_raw)), vec_raw,
               info = "_ihs -> sinh round trip")
  expect_equal(gd_auc_lvl_tr("auc_raw", vec_raw), vec_raw,
               info = noquote("no suffix -> identity"))
  ## this is the trap the helper exists for: name-based, not value-based
  expect_false(isTRUE(all.equal(gd_auc_lvl_tr("auc", log1p(vec_raw)), vec_raw)),
               info = "an unsuffixed log column is NOT back-transformed")
}

# === cluster-robust variance of a weighted mean (hand-computed) ===
gd_test_auc_crvar <- function() {
  vec_x <- c(1, 2, 3); vec_w <- c(1, 1, 1)
  ## mean 2, e = (-1, 0, 1); clusters {1,2} and {3} -> rowsums (-1, 1)
  ## sum of squares 2, divided by (sum w)^2 = 9
  expect_equal(gd_auc_crvar(vec_x, vec_w, c(1, 1, 2)), 2/9,
               info = "hand-computed CR0 variance")
  ## one row per cluster is the same arithmetic
  expect_equal(gd_auc_crvar(vec_x, vec_w, 1:3), 2/9,
               info = "no clustering coincides with the iid formula here")
  ## doubling every weight leaves a weighted mean unchanged and must leave
  ## the variance unchanged too (both numerator and denominator scale by 4)
  expect_equal(gd_auc_crvar(vec_x, rep(2, 3), c(1, 1, 2)), 2/9,
               info = "scale invariance in the weights")
  ## fewer clusters can only reduce the sum of squares (aggregation)
  set.seed(5)
  vec_r <- rnorm(200); vec_cl <- sample(1:5, 200, TRUE)
  expect_true(gd_auc_crvar(vec_r, rep(1, 200), vec_cl) <=
              gd_auc_crvar(vec_r, rep(1, 200), 1:200),
              info = "clustering cannot increase the variance")
}

# === config + coefficient-table validation ===
gd_test_auc_cfg <- function() {
  dt <- gd_syn_auc_dt()
  c_ok <- gd_syn_auc_cfg(dt)

  ## r_sbw must be a gl_sbw_proc result
  expect_error(gc_auc_euro_cfg(list(dt_results = dt)),
               pattern = "gl_sbw_proc", info = "missing config rejected")
  expect_error(gc_auc_euro_cfg(list(config = c_ok)),
               pattern = "gl_sbw_proc", info = "missing dt_results rejected")
  ## the weight and the cluster id are required
  expect_error(gc_auc_euro_cfg(list(dt_results = as.list(dt), config = list(uid = "PeopleID"))),
               pattern = "dt_results is not a data.frame",
               info = "non-frame dt_results rejected")
  r_no_w <- list(dt_results = dt[, .(PeopleID, year)],
                 config = list(uid = "PeopleID"))
  expect_error(gc_auc_euro_cfg(r_no_w), pattern = "weight column",
               info = "missing weight column rejected")
  r_no_id <- list(dt_results = dt[, .(w_sbw, year)],
                  config = list(uid = "PeopleID"))
  expect_error(gc_auc_euro_cfg(r_no_id), pattern = "uid column",
               info = "missing uid rejected")
  ## b_ref is matched against the allowed set unless numeric
  r_ok <- list(dt_results = dt, config = C_CFG_AUC)
  expect_error(gc_auc_euro_cfg(r_ok, b_ref = "banana"),
               pattern = "should be one of", info = "bad b_ref rejected")
  expect_error(gc_auc_euro_cfg(r_ok, b_ref = -2),
               pattern = "b_ref must be one of", info = "numeric b_ref rejected")
  ## defaults: wcol falls back to w_sbw, residuals are kept by construction
  expect_equal(gc_auc_euro_cfg(r_ok)$wcol, "w_sbw", info = "default weight column")

  ## coefficient table validation
  expect_error(md_auc_coef_get(NULL, NULL), pattern = "either dt_coef",
               info = "no coef source rejected")
  expect_error(md_auc_coef_get(NULL, data.table(dv = "x")),
               pattern = "missing column", info = "incomplete coef table rejected")
  ## tp arrives from fixest as numeric/double and is coerced to integer,
  ## because the pre/post switch and the column names are integer-based
  dt_c <- md_auc_coef_get(NULL, data.table(dv = "x_log", tp = 2, term = "exn_pm",
                                           estimate = 0.1, std.error = 0.02))
  expect_true(is.integer(dt_c$tp), info = "tp coerced to integer")
  ## a gl_sbw_did result is accepted as the only input: the table is then
  ## rebuilt from the model, so every lhs of the fit appears
  dt_from_did <- md_auc_coef_get(gd_fake_did_get(), NULL)
  expect_true(all(c("x_log", "x_log") %in% dt_from_did$dv),
              info = "coef table rebuilt from r_did")
  expect_true("exn_pm" %in% dt_from_did$term,
              info = "rebuilt coef table carries the treatment")
  expect_true(any(dt_from_did$term == "no_treatment"),
              info = "rebuilt coef table carries the no_treatment row")
  expect_true(all(c(2L, -2L) %in% dt_from_did$tp),
              info = "rebuilt coef table recovers post and pre periods")
}

# === row selection: which coefficients get a euro interpretation ===
gd_test_auc_sel <- function() {
  dt_c <- gd_syn_auc_coef()
  ## add a no_treatment row and a second dv, as gd_sbw_coef_did emits them
  dt_c <- rbind(dt_c,
                data.table(dv = "x_log", tp = 2L, term = "no_treatment",
                           estimate = 1.5, std.error = 0.1),
                data.table(dv = "y_log", tp = 1L, term = "exn_pm",
                           estimate = 0.3, std.error = 0.05))
  c_did <- list(outcome_cols = character(0), outcome_cols_cur = c("x_log", "y_log"))

  ## no_treatment rows are levels, not effects -> always dropped
  dt_sel <- md_auc_sel(dt_c, c_did)
  expect_false("no_treatment" %in% dt_sel$term, info = "no_treatment dropped")
  ## pre and post rows both kept (placebos are diagnostics, see the file header)
  expect_true(any(dt_sel$tp < 0) && any(dt_sel$tp > 0),
              info = "pre and post periods both kept")
  ## default dv set = the config's outcomes that were actually estimated
  expect_true(setequal(unique(dt_sel$dv), c("x_log", "y_log")),
              info = "default dv set from config")
  ## l_tp restricts the periods, negative ones included
  expect_equal(md_auc_sel(dt_c, c_did, l_tp = -3:-2)$tp, -3L,
              info = "l_tp keeps only the requested pre period")
  expect_equal(nrow(md_auc_sel(dt_c, c_did, l_tp = 2L)), 2L,
              info = "l_tp keeps only the requested post period")
  ## l_dv restricts the outcomes
  expect_equal(unique(md_auc_sel(dt_c, c_did, l_dv = "y_log")$dv), "y_log",
              info = "l_dv restricts the dv set")
  ## unknown dv is an error, not a silent empty result
  expect_error(md_auc_sel(dt_c, c_did, l_dv = "nope"), pattern = "not in the coefficient table",
               info = "unknown dv rejected")
  ## filtering everything away is an error
  expect_error(md_auc_sel(dt_c, c_did, l_tp = 99L), pattern = "nothing left",
               info = "empty selection rejected")
  ## a config with no outcomes and no l_dv has nothing to convert
  expect_error(md_auc_sel(dt_c, list(outcome_cols = character(0),
                                     outcome_cols_cur = character(0))),
               pattern = "pass l_dv", info = "no config outcomes rejected")
  ## returned object is a copy (mutating it must not touch the input)
  dt_sel[, term := "zzz"]
  expect_true("exn_pm" %in% dt_c$term, info = "selection is a copy")
}

# === the core conversion, checked algebraically ===
gd_test_auc_calc <- function() {
  dt <- gd_syn_auc_dt()
  c_eu <- gd_syn_auc_cfg(dt)
  dt_sel <- md_auc_sel(gd_syn_auc_coef(), C_CFG_AUC, l_dv = "x_log")
  dt_eu <- gd_auc_euro_calc(dt_sel, c_eu)

  ## with log1p the retransformation is exact:
  ##   expm1(log1p(y) + b) - y = (1 + y) * (exp(b) - 1)
  expect_equal(dt_eu$euro, (exp(dt_eu$estimate) - 1) * (1 + dt_eu$ref_level),
               info = "exact euro AME identity")
  ## the sign of the euro effect follows the sign of beta
  expect_true(all(sign(dt_eu$euro) == sign(dt_eu$estimate)),
              info = "euro effect has the sign of the log coefficient")
  ## the delta-method SE can only exceed the b-only SE (second variance term
  ## is non-negative); on this outcome it is a large share of the total
  expect_true(all(dt_eu$se_euro >= dt_eu$se_euro_b),
              info = "total SE >= b-only SE")
  expect_true(any(dt_eu$se_euro > 1.001 * dt_eu$se_euro_b),
              info = "reference-level variance actually contributes")
  ## reference level = weighted mean of the OBSERVED euro level on the
  ## estimation sample (residuals kept); the mask is lvl & d_out & w
  ## the column is {dv}_f2, i.e. log1p(turnover at t+2) -- NOT the level at t
  msk <- !is.na(dt$x_log_f2) & !is.na(dt$x_log_d2) & !is.na(dt$w_sbw)
  expect_equal(dt_eu[tp == 2, unique(ref_level)],
               c(dt[msk, sum(w_sbw * expm1(x_log_f2)) / sum(w_sbw)]),
               info = "post ref level is the weighted mean euro level at t+k")
  expect_equal(dt_eu[tp == 2, unique(n_used)], sum(msk),
               info = "post n_used is the reconstruction sample size")
  ## the PRE period is anchored at its own level (t-k), which on a
  ## cumulative outcome is SMALLER than the t+k level
  msk_k <- !is.na(dt$x_log_l3) & !is.na(dt$x_log_dm3) & !is.na(dt$w_sbw)
  expect_equal(dt_eu[tp == -3, ref_level],
               dt[msk_k, sum(w_sbw * expm1(x_log_l3)) / sum(w_sbw)],
               info = "pre ref level is the weighted mean euro level at t-k")
  expect_equal(dt_eu[tp == -3, n_used], sum(msk_k),
               info = "pre n_used counts only complete pre rows")
  ## share_of_level is the effect relative to the level it is applied to
  expect_equal(dt_eu$share_of_level,
               (exp(dt_eu$estimate) - 1) * (1 + 1 / dt_eu$ref_level),
               info = "share_of_level = euro / ref_level")
  ## intervals are Wald intervals around the euro estimate
  expect_equal(dt_eu$euro_hi - dt_eu$euro, dt_eu$euro - dt_eu$euro_lo,
               info = "Wald interval is symmetric on the euro scale")
  expect_equal(dt_eu$euro_lo, dt_eu$euro - 1.96 * dt_eu$se_euro,
               info = "interval is 1.96 SE wide")
  ## bookkeeping columns
  expect_true(all(dt_eu$b_ref == "t+k"), info = "b_ref recorded")
  expect_equal(attr(dt_eu, "gnrtdby"), "gd_auc_euro_calc", info = "stamp set")
  expect_equal(unique(dt_eu$dv), "x_log", info = "l_dv honoured")

}

# === nobs cross-check ===
gd_test_auc_check_n <- function() {
  dt_eu <- data.table(dv = "x_log", tp = 2L, term = "exn_pm",
                      estimate = 0.1, euro = 10, se_euro = 1,
                      n_used = as.integer(nobs(gd_fake_did_get()[[1]])))
  ## no model -> the column exists and is NA (kept for structural uniformity)
  expect_true(all(is.na(md_auc_check_n(dt_eu, NULL)$n_model)),
              info = "no r_did -> n_model is NA")
  ## a matching model passes silently
  expect_silent(md_auc_check_n(dt_eu, gd_fake_did_get()))
  expect_equal(md_auc_check_n(dt_eu, gd_fake_did_get())$n_model,
               dt_eu$n_used, info = "n_model read from the model's nobs")
  ## a reconstruction on a different row set must warn
  dt_bad <- copy(dt_eu); dt_bad[, n_used := n_used - 1L]
  expect_warning(md_auc_check_n(dt_bad, gd_fake_did_get()),
                 pattern = "row count mismatch", info = "mismatch warns")
}

# === end to end: synthetic panel -> SBW -> DiD -> coef -> euro ===
gd_test_auc_e2e <- function() {
  dt <- gd_syn_dt(1200, 2008:2014, seed = 21)
  c_sbw <- gd_syn_config(
    dt, years = 2011:2012, n_lag = 3L, impact_years = 2L, label = "test_auc_e2e",
    covariate_cols = c("Gender", "ar1"),
    outcome_cols_cur = c("auc_turnover_log", "aucturn_cum2_log"))

  r_sbw  <- t_test("auc e2e proc", gl_sbw_proc(c_sbw))
  expect_true(all(r_sbw$dt_diag$status == "solved"), info = "auc e2e years solved")

  ## the log outcome columns got the shift family, "cur" base
  expect_true(all(c("auc_turnover_log_f2", "auc_turnover_log_d2",
                    "auc_turnover_log_l3", "auc_turnover_log_dm2") %in%
                  names(r_sbw$dt_results)),
              info = "log outcome shift columns present")
  ## "cur" base: d{k} = v(t+k) - v(t) (exposure-year-free window)
  dt_res <- r_sbw$dt_results
  expect_equal(dt_res[!is.na(auc_turnover_log_d2), auc_turnover_log_d2],
               dt_res[!is.na(auc_turnover_log_d2),
                      auc_turnover_log_f2 - auc_turnover_log],
               info = "outcome_cols_cur uses base t")

  r_did   <- t_test("auc e2e did", gl_sbw_did(r_sbw, l_dvs = c("auc_turnover_log",
                                                              "aucturn_cum2_log")))
  dt_coef <- gd_sbw_coef_did(r_did)

  ## default call: every period, both dvs
  dt_eu <- t_test("auc e2e euro", gd_auc_euro(r_sbw, r_did, dt_coef, verbose = FALSE))
  expect_true(all(c("dv", "tp", "term", "estimate", "euro", "se_euro",
                    "ref_level", "n_used", "n_model") %in% names(dt_eu)),
              info = "euro table has the documented columns")
  expect_true(all(c("auc_turnover_log", "aucturn_cum2_log") %in% dt_eu$dv),
              info = "both log outcomes converted")
  expect_true(any(dt_eu$tp < 0) && any(dt_eu$tp > 0),
              info = "e2e: pre and post periods both converted")
  expect_false("no_treatment" %in% dt_eu$term, info = "e2e: no_treatment excluded")
  ## the identity holds on a real fitted model, not just on the hand panel
  expect_equal(dt_eu$euro, (exp(dt_eu$estimate) - 1) * (1 + dt_eu$ref_level),
               info = "e2e: exact euro identity")
  ## the reconstruction reproduces the model's own sample size
  expect_true(all(dt_eu$n_used == dt_eu$n_model),
              info = "e2e: reconstruction matches model nobs")
  expect_true(all(is.finite(dt_eu$euro)) && all(dt_eu$se_euro > 0),
              info = "e2e: euro estimates and SEs finite and positive")
  ## post-period reference levels are the observed levels, so they must be
  ## non-negative and, on a cumulative outcome, rising in tp
  expect_true(all(dt_eu$ref_level >= 0), info = "e2e: euro levels non-negative")
  ref_cum <- dt_eu[dv == "aucturn_cum2_log" & tp > 0 & term == "exn_pm"][order(tp)]
  expect_true(all(diff(ref_cum$ref_level) > 0),
              info = "e2e: cumulative reference level rises with the horizon")

  ## l_tp / l_dv / b_ref plumbing
  expect_equal(nrow(gd_auc_euro(r_sbw, r_did, dt_coef, l_tp = 1:2, verbose = FALSE)),
               2L * 4L * 2L, info = "e2e: l_tp=1:2 -> (tp x term x dv) rows")
  dt_pre <- gd_auc_euro(r_sbw, r_did, dt_coef, l_tp = -2L, verbose = FALSE)
  expect_true(all(dt_pre$tp == -2L), info = "e2e: l_tp=-2 keeps only the placebo")
  dt_one <- gd_auc_euro(r_sbw, r_did, dt_coef, l_dv = "aucturn_cum2_log",
                        verbose = FALSE)
  expect_equal(unique(dt_one$dv), "aucturn_cum2_log", info = "e2e: l_dv restricts")
  ## b_ref = "t-1" anchors every horizon on the {dv}_l1 column instead, so
  ## the reference level is the weighted mean of THAT column on each
  ## horizon's own mask. It is therefore not constant across horizons, and
  ## for a flow outcome it is not even ordered against the t+k level
  ## (annual turnover fluctuates).
  dt_b1 <- gd_auc_euro(r_sbw, r_did, dt_coef, b_ref = "t-1", verbose = FALSE)
  msk_b1 <- !is.na(dt_res$aucturn_cum2_log_l1) & !is.na(dt_res$aucturn_cum2_log_d2) &
            !is.na(dt_res$w_sbw)
  expect_equal(dt_b1[dv == "aucturn_cum2_log" & tp == 2 & term == "exn_pm", ref_level],
               dt_res[msk_b1, sum(w_sbw * expm1(aucturn_cum2_log_l1)) / sum(w_sbw)],
               info = "e2e: t-1 anchor is the weighted mean of the l1 column")
  ## on a CUMULATIVE outcome that anchor is below the t+k one for the same row
  expect_true(all(dt_b1[dv == "aucturn_cum2_log" & tp > 0, ref_level] <=
                  dt_eu[dv == "aucturn_cum2_log" & tp > 0, ref_level] + 1e-9),
              info = "e2e: t-1 anchor is below t+k on the cumulative outcome")

  ## dt_coef is optional: the coefficient table is rebuilt from r_did
  dt_from_did <- gd_auc_euro(r_sbw, r_did, dt_coef = NULL, verbose = FALSE)
  expect_equal(dt_from_did$euro, dt_eu$euro, info = "e2e: dt_coef is optional")
  expect_equal(attr(dt_eu, "gnrtdby"), "gd_auc_euro", info = "e2e: stamp set")
}


# === run everything ===
cat("\n--- config checks ---\n");       gd_test_check()
cat("\n--- formulas ---\n");            gd_test_formulas()
cat("\n--- covariate expansion ---\n"); gd_test_covexp()
cat("\n--- panel builder ---\n");       gd_test_panel()
cat("\n--- covariate lags ---\n");     gd_test_covlags()
cat("\n--- extra balance lags ---\n");  gd_test_bal_lags()
cat("\n--- year eligibility ---\n");    gd_test_years()
cat("\n--- fit-population eligibility ---\n"); gd_test_elig()
cat("\n--- end-to-end pipeline ---\n"); gd_test_proc()
cat("\n--- balance long ---\n");        gd_test_bal_long()
cat("\n--- DiD ---\n");                 gd_test_did()
cat("\n--- outcomes ---\n");            gd_test_outcome()
cat("\n--- outcomes (cur base) ---\n"); gd_test_outcome_cur()
cat("\n--- cache format + parallel ---\n"); gd_test_cache_format_and_parallel()
cat("\n--- aucexp: names ---\n");        gd_test_auc_names()
cat("\n--- aucexp: crvar ---\n");        gd_test_auc_crvar()
cat("\n--- aucexp: config/coef ---\n");  gd_test_auc_cfg()
cat("\n--- aucexp: selection ---\n");    gd_test_auc_sel()
cat("\n--- aucexp: conversion ---\n");   gd_test_auc_calc()
cat("\n--- aucexp: nobs check ---\n");   gd_test_auc_check_n()
cat("\n--- aucexp: end to end ---\n");   gd_test_auc_e2e()
cat("\n--- real-data smoke ---\n");     gd_test_real()

cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
