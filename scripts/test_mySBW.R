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

  r_did <- t_test("did fit", gl_sbw_did(r))
  expect_true(inherits(r_did, "fixest_multi"), info = "fixest multi-LHS object")
  dt_c <- gd_sbw_coef_did(r_did)
  expect_true(all(c("mdl", "term", "estimate", "std.error", "tp", "dv") %in% names(dt_c)),
              info = "coef table columns")
  expect_equal(nrow(dt_c), 4 * 4 * 3, info = "4 dvs x 3 tps x 4 terms")
  expect_true(all(c(-2, 1, 2) %in% unique(dt_c$tp)), info = "pre and post time points")
  ## pre-period (dm2) coefficients ~ 0 once lags are balanced
  expect_true(max(abs(dt_c[tp == -2, estimate])) < 0.1, info = "pre period flat (weighted)")
}

# === small real-data smoke (subset of the foster panel) ===
gd_test_real <- function() {
  dt_f <- fread("~/Dropbox/phd/papers/foster/misc/dt_atyear_foster.csv",
                select = c("PeopleID", "year", "Gender", "exn_pm", "exn_pi", "exn_oes"))
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
  dt_year <- gd_sbw_year_slice(dt_shuf, 2009, "year", "PeopleID",
                               c(c_sbw_fit$treatment_cols, gl_sbw_balcols(c_sbw_fit)))
  expect_equal(nrow(dt_year), nrow(dt_w), info = "shuffled slice rows match weight rows")
  ## values attached by join must equal the file weights for those (uid, tid)
  dt_check <- merge(dt_shuf[year == 2009, .(PeopleID, year)],
                    dt_w, by = c("PeopleID", "year"))
  expect_equal(sum(!is.na(dt_check$w)), nrow(dt_w), info = "every file row joins to a panel row")

  ## PSOCK parallel fit produces the same weights as the sequential fit
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

# === run everything ===
cat("\n--- config checks ---\n");       gd_test_check()
cat("\n--- formulas ---\n");            gd_test_formulas()
cat("\n--- covariate expansion ---\n"); gd_test_covexp()
cat("\n--- panel builder ---\n");       gd_test_panel()
cat("\n--- covariate lags ---\n");     gd_test_covlags()
cat("\n--- year eligibility ---\n");    gd_test_years()
cat("\n--- end-to-end pipeline ---\n"); gd_test_proc()
cat("\n--- balance long ---\n");        gd_test_bal_long()
cat("\n--- DiD ---\n");                 gd_test_did()
cat("\n--- cache format + parallel ---\n"); gd_test_cache_format_and_parallel()
cat("\n--- real-data smoke ---\n");     gd_test_real()

cat("\n========================================\n")
cat("          Timing Summary\n")
cat("========================================\n")
print(timings)
cat(sprintf("\nTotal: %.2fs\n", sum(timings$time_s)))
