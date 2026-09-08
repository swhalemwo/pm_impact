## =====================================================================
## mySBW.R — SBW pipeline (stable balancing weights, optweight::optweightMV)
##
## Turns an SBW config (single list holding all options) into a weighted
## artist-year panel via per-year optweightMV fits:
##
##   c_sbw <- list(dt_panel       = dt_atyear_foster,  # full panel, complete grid
##                 label          = "fos_k4_l3",        # unique id -> cache folder + registry
##                 uid            = "PeopleID", tid     = "year",
##                 treatment_cols = c("exn_pm", "exn_pi", "exn_oes", "exn_toti"),
##                 covariate_cols = "Gender",           # factor: auto dummy expansion
##                 n_lag          = 3L,                  # lags of all treatments balanced
##                 n_lag_cov      = 0L,                  # lags of the covariate set
##                 impact_years   = 4L,                  # outcome window d1..dK
##                 years          = NULL)                # NULL -> auto (gd_sbw_years)
##
##   r_sbw <- gl_sbw_proc(c_sbw)    # fit weights (cached per year on disk)
##   r_did <- gl_sbw_did(r_sbw)     # pooled DiD on the weighted panel
##
## Stages (each a standalone function):
##   md_sbw_defaults / gd_sbw_check   fill defaults / validate config
##   md_sbw_panel                     treatment lags + leads (kept) + d/dm cols
##   md_sbw_expand_covs               factor covariates -> numeric dummies
##   md_sbw_cov_lags                  covariate lags (n_lag_cov), info-free lags dropped
##   gd_sbw_years                     year eligibility (fit cols + outcome coverage)
##   gf_sbw_formulas                  optweightMV formula list (one per treatment)
##   gf_sbw_hash                      cache key (config spec + years + data content)
##   gl_sbw_weights                   per-year fits (cached as w_<year>.rds),
##                                    fork-parallel (mclapply), explicit weight attach
##   gd_sbw_bal_*                     balance report (weighted vs raw corr, ESS);
##                                    summary + full long table share one code path
##   gl_sbw_proc                      orchestrator -> list(dt_panel, dt_results, ...)
##   gl_sbw_did / gd_sbw_coef_did     DiD outcome model (feols) + coefficient extraction

## Result of gl_sbw_proc:
##   $dt_panel   original input panel, untouched (identical to config$dt_panel)
##   $dt_results working panel: uid/tid + treatment/covariate columns + built
##               lags (l1..) + outcome columns (d{k} / dm{k}) + the weight
##               column (wcol); one row per (uid, tid), weight NA outside the
##               fitted (complete) rows
##   $dt_years   eligible fit years with complete-row counts (n_use)
##   $dt_diag    per-year fit + balance diagnostics (n, ESS, max |corr|, ...)
##   $config     the (defaulted + covariate-expanded) config used
##   $cache      paths/hash of the weight cache used
##
## Cache/registry:
##   Per-year weights are stored as <dir_cache>/<label>/<hash>/w_<year>.rds;
##   each file is a data.table with columns (uid, tid, w) — the row identity
##   travels with the weights, so neither order nor a transient row index is
##   ever relied on (joining back is on (uid, tid)). Whenever fitting
##   happened a run row (label/hash/years/elapsed/timestamp) is appended to
##   <dir_cache>/registry.csv; a rerun with an identical config + panel
##   content (same hash) restores the weights from disk (b_force = TRUE
##   refits). diag.csv + config.rds are persisted next to the weights.
##
## Contracts / gotchas:
## - dt_panel must be a COMPLETE per-artist grid (every year between the
##   artist's first and last panel year) and must NOT contain precomputed
##   lag/lead columns (the pipeline owns the *_l/*_d/*_dm namespace).
## - Covariates: numeric columns and factor/character columns are both
##   supported. Factors are auto-expanded to numeric dummies (all levels
##   minus the first sorted level = reference) — optweightMV is infeasible
##   when passed factors directly. Rare levels (< 0.5% of rows) warn: exact
##   ATE balance with them is usually infeasible (recode). n_lag_cov lags the
##   covariate set (numeric columns AND dummies — a lag of a categorical
##   covariate is a lag of its dummies, so time-varying categorical
##   covariates are lagged meaningfully). Lags that carry no information
##   (time-invariant covariate: lag == level value) are dropped. Mixed
##   numeric + factor covariate sets are fine.
## - Covariates that contain the current-year treatment (e.g. a cumulative
##   count as of year t) cannot be balanced exactly against that treatment —
##   pass the cumulative as of t-1, or rely on the treatment lags.
## - The osqp solver is single-threaded. n_cores > 1 runs the per-year fits
##   with parallel::mclapply (forking) on Unix — dt_working is shared
##   copy-on-write so memory stays near the sequential footprint (a PSOCK
##   cluster serializes every year slice into a fresh R process, which
##   ballooned to 10+ GB with n_cores = 3). OPENBLAS_NUM_THREADS is forced
##   to 1 before forking (forking R while OpenBLAS runs threads deadlocks).
##   On Windows the code falls back to sequential with a warning.
## - eps = osqp solver tolerance (1e-3 ~2x faster than 1e-4, weights
##   essentially identical); tols = *balance* tolerance: NULL = exact (the
##   optweight default), a scalar = uniform tolerance, or a full
##   per-formula tolerance list.
## - estimand is fixed to "ATE" (optweightMV only supports ATE).
## - min_n_year only matters when years = NULL (auto): the minimum number of
##   complete rows a year needs to be kept. Explicitly requested years are
##   kept with as few as one complete row.
## - n_lag (treatment lags / balance horizon) and impact_years (outcome
##   window) are INDEPENDENT: lags go to n_lag, future leads/d to
##   impact_years. Rows are only fitted when their whole outcome window is
##   observed (see gd_sbw_years).
## - Moments/interactions are not hard-wired: precompute e.g. I(x^2) or
##   product columns and list them in covariate_cols.
## =====================================================================


md_sbw_defaults <- function(c_sbw) {
  #' Fill missing config fields with defaults. User-supplied entries always
  #' win. 
  #'
  #' @param c_sbw config list (may be partially specified)
  #' @return the completed config list
  l_dflt <- list(
    covariate_cols = character(0), # names of extra covariates to balance
    n_lag_cov      = 0L,           # lags of the covariate set (numeric + dummies)
    years          = NULL,         # fit years; NULL -> auto (see gd_sbw_years)
    min_n_year     = 10L,          # min complete rows per year for auto years
    estimand       = "ATE",        # optweightMV supports only "ATE"
    norm           = "l2",         # weight norm: l2/l1/linf/entropy
    eps            = 1e-3,         # osqp solver tolerance
    min.w          = 1e-8,         # lower bound on weights
    verbose        = FALSE,        # optweightMV verbosity
    tols           = NULL,         # balance tolerance (NULL = exact)
    n_cores        = 1L,           # years to fit in parallel (PSOCK)
    wcol           = "w_sbw",      # name of the weights column in dt_results
    dir_cache      = "~/Dropbox/phd/papers/foster/misc/sbw_cache")

  for (nm in setdiff(names(l_dflt), names(c_sbw))) c_sbw[[nm]] <- l_dflt[[nm]]
  c_sbw$dir_cache <- path.expand(c_sbw$dir_cache)
  c_sbw
}


gd_sbw_check <- function(c_sbw) {
  #' Validate a config. Errors with an explanatory message on the first
  #' problem found; returns TRUE invisibly when everything is fine.
  #'
  #' @param c_sbw a defaulted config (as returned by md_sbw_defaults)
  #' @return TRUE invisibly

  if (!is.list(c_sbw) || is.null(c_sbw$dt_panel)) {
    stop("config must be a list with at least a dt_panel entry")
  }
  if (!grepl("^[A-Za-z0-9._-]+$", c_sbw$label %||% "")) {
    stop("label must match ^[A-Za-z0-9._-]+$ (it names the cache folder)")
  }
  if (!is.data.table(c_sbw$dt_panel)) stop("dt_panel must be a data.table")

  ## identifiers
  l_need <- c("uid", "tid", "treatment_cols", "n_lag", "impact_years")
  l_miss <- setdiff(l_need, names(c_sbw))
  if (length(l_miss) > 0) stop(sprintf("config missing: %s", paste(l_miss, collapse = ", ")))
  if (!c_sbw$uid %in% names(c_sbw$dt_panel)) stop(sprintf("uid column '%s' not in dt_panel", c_sbw$uid))
  if (!c_sbw$tid %in% names(c_sbw$dt_panel)) stop(sprintf("tid column '%s' not in dt_panel", c_sbw$tid))
  if (c_sbw$uid == c_sbw$tid) stop("uid and tid must be different columns")
  if (anyNA(c_sbw$dt_panel[[c_sbw$uid]]) || anyNA(c_sbw$dt_panel[[c_sbw$tid]])) {
    stop("uid/tid must not contain NAs")
  }

  ## complete grid, no duplicates
  dt_id <- c_sbw$dt_panel[, .(mn = min(get(c_sbw$tid)), mx = max(get(c_sbw$tid)), n = .N),
                          by = c(c_sbw$uid)]
  if (uniqueN(c_sbw$dt_panel, by = c(c_sbw$uid, c_sbw$tid)) != nrow(c_sbw$dt_panel)) {
    stop("duplicate (uid, tid) rows in dt_panel")
  }
  if (!all(dt_id$n == dt_id$mx - dt_id$mn + 1)) {
    stop("dt_panel must be a complete per-artist grid (each artist has rows for ",
         "every year between first and last panel year); fill gaps before passing")
  }

  ## treatments
  l_treats <- c_sbw$treatment_cols
  if (length(l_treats) == 0 || any(!l_treats %in% names(c_sbw$dt_panel))) {
    stop("treatment_cols must be non-empty existing columns of dt_panel")
  }
  if (!all(map_lgl(l_treats, ~is.numeric(c_sbw$dt_panel[[.x]])))) {
    stop("treatment_cols must be numeric (exhibition counts)")
  }

  ## covariates
  l_covs <- c_sbw$covariate_cols
  if (any(!l_covs %in% names(c_sbw$dt_panel))) stop("some covariate_cols not in dt_panel")
  if (length(intersect(l_covs, l_treats)) > 0) stop("covariate_cols must not overlap treatment_cols")
  if (c_sbw$n_lag_cov < 0) stop("n_lag_cov must be >= 0")
  if (c_sbw$wcol %in% names(c_sbw$dt_panel)) {
    stop(sprintf("wcol '%s' already exists in dt_panel", c_sbw$wcol))
  }

  ## numeric knobs
  l_ints <- c(n_lag = c_sbw$n_lag, n_lag_cov = c_sbw$n_lag_cov,
              impact_years = c_sbw$impact_years, min_n_year = c_sbw$min_n_year,
              n_cores = c_sbw$n_cores)
  
  if (!all(map_lgl(l_ints, ~is.numeric(.x) && .x == round(.x)))) {
    stop("n_lag / n_lag_cov / impact_years / min_n_year / n_cores must be whole numbers")
  }
  if (c_sbw$n_lag < 1 || c_sbw$impact_years < 1 || c_sbw$n_lag_cov < 0 ||
      c_sbw$min_n_year < 1 || c_sbw$n_cores < 1) {
    stop("need n_lag >= 1, impact_years >= 1, n_lag_cov >= 0, min_n_year >= 1, n_cores >= 1")
  }
  if (length(c_sbw$years) > 0 && !is.numeric(c_sbw$years)) stop("years must be numeric")
  if (c_sbw$eps <= 0) stop("eps must be positive")

  ## solver options
  if (!identical(c_sbw$estimand, "ATE")) {
    stop("optweightMV only supports estimand = 'ATE'")
  }
  if (!c_sbw$norm %in% c("l2", "l1", "linf", "entropy")) {
    stop("norm must be one of 'l2', 'l1', 'linf', 'entropy'")
  }
  invisible(TRUE)
}


gd_sbw_shiftcols <- function(l_cols, n, dir = "l") {
  #' Column names produced by data.table::shift(.SD, 1:n): col-major, the lag
  #' distance is the fastest-varying index. Used so that generated names match
  #' the columns := creates (treat-major, distance-inner).
  #'
  #' @param l_cols character vector of source column names (e.g. treatments)
  #' @param n      maximum shift distance
  #' @param dir    suffix letter: "l" (lag) or "f" (lead)
  #' @return character vector of length length(l_cols) * n
  paste0(rep(l_cols, each = n), "_", dir, rep(seq_len(n), times = length(l_cols)))
}


gl_sbw_balcols <- function(c_sbw) {
  #' The balance columns: everything each treatment formula balances on —
  #' lagged treatments (t-1..t-n_lag of every treatment type) followed by the
  #' covariate set. IMPORTANT: config$covariate_cols must already be the
  #' FINAL covariate list — md_sbw_expand_covs turns factor columns into
  #' numeric dummies and md_sbw_cov_lags appends their lags — so this
  #' function just adds the names as-is (no lag expansion here; doing it here
  #' would double-lag). 
  #'
  #' @param c_sbw defaulted + covariate-expanded/lagged config
  #' @return character vector of column names
  c(gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
    c_sbw$covariate_cols)
}


gf_sbw_formulas <- function(c_sbw) {
  #' One optweightMV formula per treatment: exn_X ~ <all treatment lags> +
  #' <covariate set>. Same-year other treatments stay out of the RHS so each
  #' fitted coefficient is a partial dose effect conditional on the observed
  #' other treatments of that year.
  #'
  #' @param c_sbw defaulted config; covariate_cols must be the FINAL set
  #'   (numeric covariates + expanded factor dummies + their lags, as
  #'   prepared by md_sbw_expand_covs / md_sbw_cov_lags)
  #' @return named list of formulas, names = treatment_cols

  l_rhs <- paste(c(gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
                   c_sbw$covariate_cols),
                 collapse = " + ")

  map(c_sbw$treatment_cols, ~as.formula(sprintf("%s ~ %s", .x, l_rhs))) %>%
    setNames(c_sbw$treatment_cols)
}


md_sbw_panel <- function(c_sbw) {
  #' Build the working panel from the base dt_panel (a lean copy; the input
  #' is never modified): keeps uid/tid/treatments/covariates and adds
  #'   - treatment lags exn_X_l1..l{n_lag}   (n_lag only; NOT tied to
  #'     impact_years — lags beyond the balance horizon would anyway be
  #'     unbalanced and only serve as placebo pre-periods)
  #'   - treatment leads exn_X_f1..f{K}      (K = impact_years) — KEPT in the
  #'     panel: raw future values are useful in their own right (e.g. Poisson
  #'     outcome models) and gd_sbw_years uses them for outcome-coverage
  #'     checks
  #'   - outcome columns exn_X_d{k}  = lead k - l1  (k = 1..K)
  #'   - outcome columns exn_X_dm{k} = lag k - l1   (k = 2..n_lag; dm1 = 0 is
  #'     the DiD baseline by construction and is not stored)
  #' Covariates are NOT lagged here — factor covariates must first be
  #' expanded to numeric dummies (md_sbw_expand_covs) and their lags are
  #' built afterwards (md_sbw_cov_lags), so that a lag of a factor level is a
  #' lag of a dummy, not a lag of a character column.
  #'
  #' @param c_sbw defaulted config
  #' @return data.table: the working panel

  l_treats <- c_sbw$treatment_cols
  l_covs   <- c_sbw$covariate_cols
  K        <- c_sbw$impact_years
  max_lag  <- c_sbw$n_lag

  ## generated names must not collide with base panel columns (the pipeline
  ## owns the *_l/*_f/*_d/*_dm namespace)
  l_new <- c(gd_sbw_shiftcols(l_treats, max_lag, "l"),  # lags
             gd_sbw_shiftcols(l_treats, K, "f"),        # leads
             sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K)), # forward diffs, backward diffs below
             if (max_lag > 1) sprintf("%s_dm%d", rep(l_treats, each = max_lag - 1), 2:max_lag) else character(0))
  
  l_clash <- intersect(l_new, names(c_sbw$dt_panel))
  if (length(l_clash) > 0) {
    stop(sprintf("dt_panel already contains pipeline-generated columns: %s. ",
                 paste(l_clash, collapse = ", ")),
         "Pass the base panel (no precomputed lags/leads).")
  }

  dt_working <- copy(c_sbw$dt_panel[, c(c_sbw$uid, c_sbw$tid, l_treats, l_covs), with = FALSE])
  setorderv(dt_working, c(c_sbw$uid, c_sbw$tid))

  ## lags of treatments (to max_lag); collapse::L lags within the uid groups
  l_lagcols <- gd_sbw_shiftcols(l_treats, max_lag, "l")
  dt_working[, (l_lagcols) := collapse::L(.SD, n = 1:max_lag, by = get(c_sbw$uid)), .SDcols = l_treats]

  ## leads of treatments (kept in the panel; also used for d{k})
  l_fcols <- gd_sbw_shiftcols(l_treats, K, "f")
  dt_working[, (l_fcols) := collapse::L(.SD, n = -(1:K), by = get(c_sbw$uid)), .SDcols = l_treats]

  ## d{k} = lead{k} - l1
  l_l1cols <- gd_sbw_shiftcols(l_treats, 1, "l")
  l_dcols  <- sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K))
  dt_working[, (l_dcols) := map2(l_fcols, rep(l_l1cols, each = K), ~get(.x) - get(.y))]

  ## dm{k} = lag{k} - l1  (k = 2..max_lag)
  if (max_lag > 1) {
    l_dmcols  <- sprintf("%s_dm%d", rep(l_treats, each = max_lag - 1), 2:max_lag)
    ## source lags: every column except each treatment's _l1 (positions
    ## 1, 1+max_lag, ... in the col-major name order)
    l_dmsrcs  <- l_lagcols[-seq(1, length(l_lagcols), by = max_lag)]
    dt_working[, (l_dmcols) := map2(l_dmsrcs, rep(l_l1cols, each = max_lag - 1), ~get(.x) - get(.y))]
  }

  return(dt_working)
}


md_sbw_dum_for_cov <- function(dt_working, cov_col) {
  #' Dummy matrix for one factor/character covariate: one column per level
  #' minus the first (sorted level order = the dropped reference, consistent
  #' with model.matrix). Column names are make.names-sanitized so they stay
  #' syntactically valid when they appear in a formula (level labels like
  #' "Artist duo, Female" would otherwise break the formula parser).
  #'
  #' @param dt_working data.table containing cov_col
  #' @param cov_col    name of the factor/character column to expand
  #' @return numeric matrix with length(levels) - 1 columns
  x <- dt_working[[cov_col]]
  l_lv <- sort(unique(x[!is.na(x)]))
  if (length(l_lv) < 2) stop(sprintf("covariate '%s' has fewer than 2 levels", cov_col))
  m <- model.matrix(~ factor(x) - 1)   ## one column per level
  colnames(m) <- l_lv
  m <- m[, -1, drop = FALSE]           ## drop the reference level
  colnames(m) <- make.names(colnames(m))
  m
}


gd_sbw_warn_rare <- function(dt_working, cov_col) {
  #' Warn when a covariate level is tiny (< 0.5% of the panel rows). Exact
  #' ATE balance against such a dummy is typically infeasible per year; the
  #' user should recode the covariate into few meaningful levels first (e.g.
  #' the project's Male/Female/Other gender recode).
  #'
  #' @param dt_working data.table containing cov_col
  #' @param cov_col    name of the factor/character covariate

  l_rare <- names(which(table(dt_working[[cov_col]], useNA = "no") < (0.005 * nrow(dt_working))))

  if (length(l_rare) > 0) {
    warning(sprintf("covariate '%s' has rare level(s): %s. ", cov_col,
                    paste(l_rare, collapse = ", ")),
            "Exact ATE balance with such dummies is typically infeasible; ",
            "recode the covariate into a few meaningful levels before balancing.")
  }
}


md_sbw_expand_covs <- function(dt_working, l_covs) {
  #' Expand factor/character covariates into numeric dummies (one per level
  #' minus the first, sorted = reference); numeric covariates stay as-is.
  #' optweightMV is infeasible when passed factor columns directly, so the
  #' pipeline always balances on the expanded dummies.
  #'
  #' @param dt_working working panel (as built by md_sbw_panel)
  #' @param l_covs     covariate names to balance
  #' @return list(dt = dt_working with dummies added and factor columns
  #'              removed, covs = the names actually to balance afterwards)

  l_keep <- l_covs[map_lgl(l_covs, ~is.numeric(dt_working[[.x]]))]
  l_exp  <- setdiff(l_covs, l_keep)
  if (length(l_exp) == 0) return(list(dt = dt_working, covs = l_covs))

  l_dum <- map(l_exp, md_sbw_dum_for_cov, dt_working = dt_working) %>%
    setNames(l_exp)
  l_dum_names <- map2(names(l_dum), l_dum, ~sprintf("%s_%s", .x, colnames(.y))) %>%
    unlist
  if (any(l_dum_names %in% names(dt_working))) {
    stop(sprintf("dummy columns would clash with existing columns: %s",
                 paste(intersect(l_dum_names, names(dt_working)), collapse = ", ")))
  }
  walk(l_exp, gd_sbw_warn_rare, dt_working = dt_working)

  m_dum <- do.call(cbind, unname(l_dum))
  colnames(m_dum) <- l_dum_names

  dt_working <- cbind(dt_working, data.table::as.data.table(m_dum))
  dt_working[, (l_exp) := NULL]        ## factor info now lives in the dummies
  list(dt = dt_working, covs = c(l_keep, l_dum_names))
}


md_sbw_cov_lags <- function(dt_working, c_sbw) {
  #' Build covariate lags. Runs AFTER md_sbw_expand_covs, so every covariate
  #' in config$covariate_cols is numeric (raw numeric columns or factor
  #' dummies) and can be lagged with a plain shift. Adds {cov}_l1..l{n_lag_cov}
  #' to the panel and rewrites config$covariate_cols to the FINAL covariate
  #' set (levels + surviving lags) — gl_sbw_balcols / gf_sbw_formulas then
  #' consume that set directly.
  #'
  #' A lag column is DROPPED when it carries no information: for a
  #' time-invariant covariate (constant within artist, e.g. Gender or
  #' BirthYear) the lag is identical to the level value, so balancing it
  #' would only duplicate an existing constraint. Time-varying covariates
  #' (numeric or categorical dummies) keep their lags.
  #'
  #' @param dt_working working panel with expanded (numeric) covariates
  #' @param c_sbw      defaulted config (covariate_cols already expanded)
  #' @return list(dt = dt_working with covariate lag columns added,
  #'              covs = the final covariate set to balance)
  if (c_sbw$n_lag_cov < 1 || length(c_sbw$covariate_cols) == 0) {
    return(list(dt = dt_working, covs = c_sbw$covariate_cols))
  }

  l_covs <- c_sbw$covariate_cols
  n_lag  <- c_sbw$n_lag_cov
  l_all  <- gd_sbw_shiftcols(l_covs, n_lag, "l")     ## cov-major, lag-inner

  dt_working[, (l_all) := collapse::L(.SD, n = 1:n_lag, by = get(c_sbw$uid)),
             .SDcols = l_covs]

  ## keep only the lag columns that carry within-artist information
  ## (collapse::varying: FALSE for a lag of a time-invariant covariate,
  ## which is identical to the level value and would only duplicate the
  ## constraint)
  l_keep <- collapse::varying(dt_working, by = dt_working[[c_sbw$uid]], cols = l_all)[l_all]

  l_drop <- l_all[!l_keep]
  if (length(l_drop) > 0) dt_working[, (l_drop) := NULL]

  list(dt = dt_working, covs = c(l_covs, l_all[l_keep]))
}


gd_sbw_years <- function(c_sbw, dt_working) {
  #' Eligible fit years. A year is eligible when it has at least one row that
  #' is complete on EVERYTHING the pipeline needs from that row: current
  #' treatments, balance columns (treatment lags + covariate set) AND the
  #' future outcome leads f1..fK (the coverage gate — per row, so artists
  #' whose panel ends within the outcome window are excluded automatically;
  #' no global "last year - K" cutoff is needed). This is the same row
  #' definition gl_sbw_weights later fits on, so every weighted row has its
  #' d1..dK outcomes.
  #' With config years = NULL additionally require n_use >= min_n_year;
  #' explicitly requested years are honored with >= 1 complete row (requested
  #' but ineligible years are dropped with a warning).
  #'
  #' @param c_sbw      defaulted config (covariate set final)
  #' @param dt_working working panel (md_sbw_panel output)
  #' @return data.table with columns (year, n_use) for the eligible years,
  #'   sorted ascending

  tid_col  <- c_sbw$tid
  K         <- c_sbw$impact_years
  l_fitcols <- c(c_sbw$treatment_cols,                       # treatments at t
                 gl_sbw_balcols(c_sbw),                      # lagged treatments + covariates
                 gd_sbw_shiftcols(c_sbw$treatment_cols, K, "f")) # future outcome leads

  dt_years <- dt_working[, .(n_use = sum(complete.cases(.SD))), by = c(tid_col), .SDcols = l_fitcols] %>%
    .[n_use >= 1]

  setnames(dt_years, tid_col, "year")

  v_requested <- c_sbw$years
  if (length(v_requested) > 0) {
    l_miss <- setdiff(v_requested, dt_years$year)
    if (length(l_miss) > 0) {
      warning(sprintf("gd_sbw_years: requested years not eligible (no complete rows or no outcome coverage): %s",
                      paste(l_miss, collapse = ", ")))
    }
    dt_years <- dt_years[year %in% v_requested]
  } else {
    dt_years <- dt_years[n_use >= c_sbw$min_n_year]
  }
  dt_years[order(year)]
  
}


gf_sbw_hash <- function(c_sbw, dt_years, dt_working) {
  #' md5 cache key: the config spec (minus dt_panel/label, which are the data
  #' and the folder name) + the selected years + the content of the used
  #' working-panel columns. Numerics are rounded to 8 digits first to keep
  #' the hash stable across float noise.
  #'
  #' @param c_sbw      defaulted config (covariate set final)
  #' @param dt_years   eligible years table (gd_sbw_years output)
  #' @param dt_working working panel (its lag/lead columns determine the fit)
  #' @return md5 hex string
  c_spec <- c_sbw
  c_spec$dt_panel <- NULL
  c_spec$label <- NULL
  l_spec <- c_spec[order(names(c_spec))]
  spec_str <- map_chr(names(l_spec), ~sprintf("%s=%s", .x,
                                              paste(l_spec[[.x]], collapse = "|"))) %>%
    paste(collapse = ";")

  ## content of all data columns the fit depends on: ids, treatments,
  ## covariate set AND the future outcome leads (they define the fitted row
  ## population — the coverage pattern is part of what the weights are
  ## fitted on)
  l_cols <- unique(c(c_sbw$uid, c_sbw$tid, c_sbw$treatment_cols,
                     c_sbw$covariate_cols,
                     gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$impact_years, "f")))
  dt0 <- dt_working[, ..l_cols]
  data_str <- map_chr(dt0, function(x) {
    if (is.numeric(x)) digest::digest(round(x, 8), algo = "md5")
    else digest::digest(x, algo = "md5", serialize = FALSE)
  }) %>% paste(collapse = "")

  digest::digest(paste(spec_str, paste(dt_years$year, collapse = ","), data_str, sep = "|"),
                 algo = "md5", serialize = FALSE)

}



gd_sbw_fit_year <- function(dt_working, year, p_out, c_sbw, l_fs, l_fit_args) {
  #' Fit optweightMV for ONE year: slice the year's complete fit rows from
  #' the working panel, solve, write the weights file (uid, tid, w) to disk.
  #' Runs either in the master process or in a forked worker.
  #'
  #' @param dt_working  working panel (md_sbw_panel output)
  #' @param year        calendar year to fit
  #' @param p_out       full path of the w_<year>.rds file to write
  #' @param c_sbw       defaulted + covariate-expanded config (uid/tid)
  #' @param l_fs        optweightMV formula list (gf_sbw_formulas output)
  #' @param l_fit_args  optweightMV arguments (identical across years)
  #' @return single-row data.table (n_fit, status, iter, run_time) solver
  #'   diagnostics
  dt_year <- gd_sbw_year_slice(dt_working, year, c_sbw)

  r <- do.call(optweight::optweightMV,
               c(list(formula.list = l_fs, data = dt_year), l_fit_args))
  if (!identical(r$info$status, "solved")) {
    stop(sprintf("optweightMV did not solve for year %s (status: %s). ",
                 year, r$info$status),
         "Try: relax `tols`, drop/recode covariates (factors are auto-expanded ",
         "to dummies; avoid covariates containing the current-year treatment, ",
         "e.g. cumulative counts as of year t), reduce n_lag, or trim extreme values.")
  }

  dt_w <- data.table(uid = dt_year[[c_sbw$uid]], tid = dt_year[[c_sbw$tid]],
                     w = as.numeric(r$weights))
  data.table::setnames(dt_w, c("uid", "tid"), c(c_sbw$uid, c_sbw$tid))
  saveRDS(dt_w, p_out)

  data.table(n_fit = nrow(dt_year), status = r$info$status, iter = r$info$iter,
             run_time = r$info$run_time)
}

gd_sbw_year_slice <- function(dt_working, year, c_sbw) {
  #' One year's rows of the working panel that the fit needs: uid, tid, the
  #' treatments, the balance columns and the future outcome leads (the
  #' coverage gate), complete rows only (optweightMV errors on any NA). The
  #' fit population therefore equals the outcome population: every weighted
  #' row has its d1..dK outcomes.
  #'
  #' @param dt_working working panel (md_sbw_panel output)
  #' @param year       calendar year to slice
  #' @param c_sbw      defaulted config (uid/tid + treatment/impact cols)
  #' @return data.table with columns (uid, tid, treatments, balcols, leads),
  #'   complete cases only
  v_year <- year   ## local scalar: `year` is also a column name in dt_working
  l_cols <- c(c_sbw$uid, c_sbw$tid, c_sbw$treatment_cols,
              gl_sbw_balcols(c_sbw),
              gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$impact_years, "f"))
  l_fitcols <- setdiff(l_cols, c(c_sbw$uid, c_sbw$tid))
  dt_year <- dt_working[get(c_sbw$tid) == v_year, .SD, .SDcols = l_cols]
  dt_year[complete.cases(dt_year[, .SD, .SDcols = l_fitcols])]
}


gd_sbw_read_year <- function(dt_working, year, dir_run, c_sbw) {
  #' Load one year's cached weights from disk and validate them against the
  #' year's complete fit rows (row count + (uid, tid) membership — catches a
  #' stale/mismatched cache instead of silently mis-assigning). PURE: no side
  #' effect on the working panel; gl_sbw_weights attaches all years at once.
  #'
  #' @param dt_working working panel (md_sbw_panel output)
  #' @param year       calendar year
  #' @param dir_run    cache directory holding w_<year>.rds
  #' @param c_sbw      defaulted config (uid/tid)
  #' @return data.table (uid, tid, w) of the validated weights for that year
  dt_year <- gd_sbw_year_slice(dt_working, year, c_sbw)
  dt_w <- readRDS(file.path(dir_run, sprintf("w_%s.rds", year)))
  if (nrow(dt_w) != nrow(dt_year)) {
    stop(sprintf("cached weights for year %s have %d rows but %d complete fit rows",
                 year, nrow(dt_w), nrow(dt_year)))
  }
  if (!identical(sort(dt_w[[c_sbw$uid]]), sort(dt_year[[c_sbw$uid]])) ||
      !identical(sort(dt_w[[c_sbw$tid]]), sort(dt_year[[c_sbw$tid]]))) {
    stop(sprintf("cached (uid, tid) pairs for year %s do not match the current panel", year))
  }
  dt_w[, .SD, .SDcols = c(c_sbw$uid, c_sbw$tid, "w")]
}


gd_sbw_year_weights <- function(dt_working, year, dir_run, c_sbw, l_fs,
                                l_fit_args) {
  #' One worker task: the complete weight pipeline for a single year.
  #' Fits optweightMV when the year is not yet cached (gd_sbw_fit_year),
  #' then reads the weights back (gd_sbw_read_year) and returns them together
  #' with the fit diagnostics. Runs either in the master process or in a
  #' forked worker (dt_working shared copy-on-write).
  #'
  #' @param dt_working  working panel (md_sbw_panel output)
  #' @param year        calendar year
  #' @param dir_run     cache directory holding w_<year>.rds
  #' @param c_sbw       defaulted + covariate-expanded config
  #' @param l_fs        optweightMV formula list
  #' @param l_fit_args  optweightMV arguments (identical across years)
  #' @return list(dt_w = data.table (uid, tid, w) of the year's weights,
  #'              dt_fit = single-row diagnostics (n_fit, status, iter,
  #'              run_time, from_cache))
  p_out <- file.path(dir_run, sprintf("w_%s.rds", year))
  if (!file.exists(p_out)) {
    dt_fit <- gd_sbw_fit_year(dt_working, year, p_out, c_sbw, l_fs, l_fit_args)
    dt_fit[, from_cache := FALSE]
  } else {
    dt_fit <- data.table(n_fit = NA_integer_, status = "cached",
                         iter = NA_real_, run_time = NA_real_,
                         from_cache = TRUE)
  }
  dt_w <- gd_sbw_read_year(dt_working, year, dir_run, c_sbw)
  dt_fit[, n_fit := nrow(dt_w)]
  list(dt_w = dt_w, dt_fit = dt_fit)
}


gl_sbw_fit_missing_years <- function(c_sbw, dt_working, v_miss, dir_run, l_fs,
                                     l_fit_args) {
  #' Fit the per-year weights for all years that are not yet cached. Runs
  #' sequentially, or with parallel::mclapply (fork) when n_cores > 1 on
  #' Unix; OPENBLAS_NUM_THREADS is forced to 1 before forking (forking R
  #' while OpenBLAS runs threads deadlocks). Windows falls back to
  #' sequential with a warning.
  #'
  #' @param c_sbw       defaulted + covariate-expanded config
  #' @param dt_working  working panel (md_sbw_panel output)
  #' @param v_miss      years without cached weights
  #' @param dir_run     cache directory to write w_<year>.rds into
  #' @param l_fs        optweightMV formula list
  #' @param l_fit_args  optweightMV arguments (identical across years)
  #' @return named list (one entry per year in v_miss) of the results of
  #'   gd_sbw_year_weights
  if (length(v_miss) == 0) return(list())

  if (c_sbw$n_cores > 1 && .Platform$OS.type == "unix") {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    l_res <- parallel::mclapply(v_miss, gd_sbw_year_weights,
                                dt_working = dt_working, dir_run = dir_run,
                                c_sbw = c_sbw, l_fs = l_fs,
                                l_fit_args = l_fit_args,
                                mc.cores = c_sbw$n_cores,
                                mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else {
    if (c_sbw$n_cores > 1) {
      warning("n_cores > 1 only parallelizes on Unix (fork); running sequentially")
    }
    l_res <- map(v_miss, gd_sbw_year_weights, dt_working = dt_working,
                 dir_run = dir_run, c_sbw = c_sbw, l_fs = l_fs,
                 l_fit_args = l_fit_args)
  }
  setNames(l_res, v_miss)
}


gl_sbw_weights <- function(c_sbw, dt_working, dt_years, l_fs, b_force = FALSE) {
  #' Per-year optweightMV fits with an on-disk cache (one w_<year>.rds per
  #' year, each storing the weights together with their (uid, tid) rows).
  #' Fits the missing years (gl_sbw_fit_missing_years, sequential or forked
  #' parallel), then loads EVERY year's weights and attaches them to
  #' dt_working in one explicit join on (uid, tid) under config$wcol (by
  #' reference) — there is no per-year side effect, so a failed fit leaves
  #' dt_working untouched. Rows outside the fitted complete rows keep NA.
  #'
  #' @param c_sbw      defaulted + covariate-expanded config
  #' @param dt_working working panel (md_sbw_panel output; must NOT yet
  #'                   contain the wcol column)
  #' @param dt_years   eligible years (gd_sbw_years output)
  #' @param l_fs       formula list (gf_sbw_formulas output)
  #' @param b_force    TRUE: delete cached weights and refit all years
  #' @return list(dt_fit = per-year diagnostics, dir_run, hash,
  #'              from_cache_all, elapsed_s); dt_working gains the wcol
  #'              column as a side effect of the final join
  if (c_sbw$wcol %in% names(dt_working)) {
    stop(sprintf("wcol '%s' already present in dt_working; pass a fresh panel", c_sbw$wcol))
  }

  ## cache location (content-addressed within the label folder)
  hash    <- gf_sbw_hash(c_sbw, dt_years, dt_working)
  dir_run <- file.path(c_sbw$dir_cache, c_sbw$label, substr(hash, 1, 12))
  dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)

  l_paths <- setNames(file.path(dir_run, sprintf("w_%s.rds", dt_years$year)),
                      dt_years$year)
  if (b_force) walk(l_paths, ~if (file.exists(.x)) file.remove(.x))
  v_cached <- map_lgl(l_paths, file.exists)
  v_miss   <- dt_years$year[!v_cached]

  t0 <- Sys.time()

  ## optweightMV arguments that are identical across years
  l_fit_args <- gl_sbw_fit_args(c_sbw, l_fs)

  ## fit the missing years (parallel across years when n_cores > 1)
  l_res <- gl_sbw_fit_missing_years(c_sbw, dt_working, v_miss, dir_run, l_fs,
                                    l_fit_args)

  ## load ALL years' weights (fitted above or from cache) and attach in ONE
  ## join on (uid, tid); v_cached (the original file-existence detection)
  ## must NOT be overwritten — from_cache_all below relies on it
  v_all  <- dt_years$year
  v_from_cache <- setdiff(v_all, v_miss)
  l_res_cached <- map(v_from_cache, gd_sbw_year_weights, dt_working = dt_working,
                      dir_run = dir_run, c_sbw = c_sbw, l_fs = l_fs,
                      l_fit_args = l_fit_args) %>% setNames(v_from_cache)
  l_res <- c(l_res, l_res_cached)[as.character(v_all)]

  dt_all_w <- rbindlist(map(l_res, ~.x$dt_w))
  dt_working[dt_all_w, on = c(c_sbw$uid, c_sbw$tid), (c_sbw$wcol) := i.w]

  dt_fit <- rbindlist(map(l_res, ~.x$dt_fit))
  dt_fit[, year := v_all]
  setorderv(dt_fit, "year")

  list(dt_fit = dt_fit, dir_run = dir_run, hash = hash,
       from_cache_all = all(v_cached),
       elapsed_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2))
}


gl_sbw_fit_args <- function(c_sbw, l_fs) {
  #' optweightMV arguments that are identical across all years.
  #' tols: NULL = exact balance (optweight default), a numeric scalar =
  #' uniform tolerance, otherwise passed through as the full per-formula
  #' tolerance list.
  #'
  #' @param c_sbw defaulted + covariate-expanded config
  #' @param l_fs  optweightMV formula list (number of formulas determines the
  #'              length of a uniform tols repetition)
  #' @return named list of optweightMV arguments (estimand, norm, min.w,
  #'   verbose, eps_rel, eps_abs, tols.list)
  list(estimand = c_sbw$estimand, norm = c_sbw$norm,
       min.w = c_sbw$min.w, verbose = c_sbw$verbose,
       eps_rel = c_sbw$eps, eps_abs = c_sbw$eps,
       tols.list = if (is.null(c_sbw$tols)) {
         list(0)              ## exact balance (optweight default)
       } else if (is.numeric(c_sbw$tols) && length(c_sbw$tols) == 1) {
         rep(list(c_sbw$tols), length(l_fs))  ## uniform tolerance
       } else {
         c_sbw$tols          ## full per-formula tolerance list
       })
}


gd_sbw_cor_w <- function(m, w) {
  #' Weighted correlation matrix of the columns of m.
  #'
  #' Hand-rolled on purpose instead of weights::wtd.cor / Hmisc:
  #'   - weights::wtd.cor ERRORS when a column has zero variance ("0 (non-NA)
  #'     cases") — exactly what happens with a rare dummy or a constant
  #'     covariate within a year;
  #'   - the ~10 lines here are the textbook weighted covariance formula
  #'     (crossprod on mean-centered, weight-scaled columns), vectorized,
  #'     dependency-free, and easy to step through. Zero-variance columns are
  #'     set to NA instead of erroring.
  #'
  #' Row alignment: w must be in the SAME row order as m. The only caller
  #' (gd_sbw_bal_year) builds both from one dt_year without any filtering in
  #' between, so alignment holds by construction (guarded by a length check
  #' there).
  #'
  #' @param m numeric matrix (rows = observations)
  #' @param w numeric weight vector (length nrow(m)); w = rep(1, nrow(m))
  #'   yields the unweighted correlation matrix
  #' @return matrix; NA where a column has (near-)zero variance
  sw  <- sum(w)
  mw  <- colSums(m * w) / sw
  mc  <- m - rep(mw, each = nrow(m))
  covm <- crossprod(mc * w, mc) / sw
  sds <- sqrt(diag(covm))
  r   <- covm / outer(sds, sds)
  r[sds < 1e-12, ] <- NA
  r[, sds < 1e-12] <- NA
  r
}


gd_sbw_mat_col <- function(dt_year, col) {
  #' Matrix form of one column of a (weighted) year: numeric columns as-is;
  #' factor/character columns expanded to one dummy per level (used only by
  #' the balance diagnostics — the fit itself gets pre-expanded dummies).
  #'
  #' @param dt_year one year of the weighted panel
  #' @param col     column name
  #' @return numeric matrix with 1 (numeric) or k (factor) columns
  x <- dt_year[[col]]
  if (is.numeric(x)) {
    m1 <- as.matrix(x)
    colnames(m1) <- col
    m1
  } else {
    m1 <- model.matrix(~ factor(x) - 1)
    colnames(m1) <- make.names(paste0(col, "_", colnames(m1)))
    m1
  }
}


gd_sbw_mat <- function(dt_year, l_cols) {
  #' Column-bind the given columns of dt_year in matrix form (cbind helper
  #' over gd_sbw_mat_col).
  #'
  #' @param dt_year one year of the weighted panel
  #' @param l_cols  column names to include
  #' @return numeric matrix
  l_m <- map(l_cols, gd_sbw_mat_col, dt_year = dt_year)
  ## gd_sbw_mat_col(dt_year, "exn_pm_l1") %>% dim
  do.call(cbind, l_m)
}


gd_sbw_bal_year <- function(dt_year, l_treats, l_balcols, wcol) {
  #' Weighted vs raw correlations (treatment x balance column) for one year
  #' of the weighted panel.
  #'
  #' Transparency of row order: the weights vector and the design matrix are
  #' both extracted from the SAME dt_year (no subsetting in between), so the
  #' k-th row of m corresponds to the k-th element of w by construction; the
  #' length check guards against future edits that break that invariant.
  #'
  #' @param dt_year   one year of the weighted panel (weighted rows only)
  #' @param l_treats  treatment column names
  #' @param l_balcols balance column names (gl_sbw_balcols output)
  #' @param wcol      weight column name
  #' @return data.table with columns (treat, balcol, cor_raw, cor_w)

  w <- dt_year[[wcol]]
  m <- gd_sbw_mat(dt_year, c(l_treats, l_balcols))
  stopifnot(length(w) == nrow(m))

  n_t   <- length(l_treats)
  r_w   <- gd_sbw_cor_w(m, w)
  r_raw <- gd_sbw_cor_w(m, rep(1, nrow(m)))

  data.table(treat   = rep(l_treats, times = ncol(m) - n_t),
             balcol  = rep(colnames(m)[-seq_len(n_t)], each = n_t),
             cor_raw = as.vector(r_raw[seq_len(n_t), -seq_len(n_t), drop = FALSE]),
             cor_w   = as.vector(r_w[seq_len(n_t), -seq_len(n_t), drop = FALSE]))
}


gd_sbw_bal_long_year <- function(dt_results, year, tid_col, wcol, l_treats,
                                 l_balcols) {
  #' Balance rows for one weighted year of dt_results.
  #'
  #' @param dt_results weighted panel (gl_sbw_proc output)
  #' @param year       calendar year to evaluate
  #' @param tid_col    time column name
  #' @param wcol       weight column name
  #' @param l_treats   treatment column names
  #' @param l_balcols  balance column names
  #' @return data.table like gd_sbw_bal_year with an added year column
  v_year <- year   ## local scalar: `year` is also a column name in dt_results
  dt_year <- dt_results[get(tid_col) == v_year & !is.na(get(wcol))]
  dt_b <- gd_sbw_bal_year(dt_year, l_treats, l_balcols, wcol)
  dt_b[, year := v_year]
  dt_b
}


gd_sbw_bal_long <- function(dt_results, c_sbw, v_years = NULL, wcol = NULL) {
  #' Long balance table across years: rows year x treat x balcol with raw and
  #' weighted correlations. Only weighted rows are used; v_years defaults to
  #' all years that have weights.
  #'
  #' @param dt_results weighted panel (gl_sbw_proc output)
  #' @param c_sbw      config (as stored in the result)
  #' @param v_years    optional subset of years; NULL = all weighted years
  #' @param wcol       weight column (defaults to config$wcol)
  #' @return data.table with columns (year, treat, balcol, cor_raw, cor_w)
  c_sbw <- md_sbw_defaults(c_sbw)
  wcol    <- wcol %||% c_sbw$wcol
  tid_col <- c_sbw$tid
  v_years <- v_years %||% dt_results[!is.na(get(wcol)), sort(unique(get(tid_col)))]

  map(v_years, gd_sbw_bal_long_year, dt_results = dt_results, tid_col = tid_col,
      wcol = wcol, l_treats = c_sbw$treatment_cols,
      l_balcols = gl_sbw_balcols(c_sbw)) %>% rbindlist
}


gd_sbw_bal_wstats_year <- function(dt_year, wcol) {
  #' Weight statistics (ESS, mean weight) for one weighted year.
  #'
  #' @param dt_year one year of the weighted panel (weighted rows only)
  #' @param wcol    weight column name
  #' @return single-row data.table (ess, mean_w)
  w <- dt_year[[wcol]]
  data.table(ess = round(sum(w)^2 / sum(w^2), 1),
             mean_w = round(mean(w), 4))
}


gd_sbw_bal_summary <- function(dt_results, c_sbw, v_years = NULL, wcol = NULL,
                               b_return_long = FALSE) {
  #' Per-year balance summary over the weighted panel. Computed from the SAME
  #' long individual-correlation table that gd_sbw_bal_long returns (one code
  #' path — no separate per-year matrix code), aggregated to per-year maxima
  #' and joined with the per-year weight statistics.
  #'
  #' @param dt_results     weighted panel
  #' @param c_sbw          config
  #' @param v_years        optional subset; NULL = all weighted years
  #' @param wcol           weight column (defaults to config$wcol)
  #' @param b_return_long  TRUE: return the full long table (individual
  #'   treat x balcol correlations for every year) instead of the summary
  #' @return data.table: default one row per year
  #'   (year, ess, mean_w, max_cor_raw, max_cor_w); with b_return_long = TRUE
  #'   the long table (year, treat, balcol, cor_raw, cor_w)

  c_sbw <- md_sbw_defaults(c_sbw)
  wcol    <- wcol %||% c_sbw$wcol
  tid_col <- c_sbw$tid
  v_years <- v_years %||% dt_results[!is.na(get(wcol)), sort(unique(get(tid_col)))]

  dt_long <- gd_sbw_bal_long(dt_results, c_sbw, v_years = v_years, wcol = wcol)
  if (b_return_long) return(dt_long)

  ## per-year weight statistics
  ## FIXME: assingment in lambdas
  dt_wstats <- map(v_years, function(year) {
    v_year <- year   ## `year` is also a column name in dt_results
    dt_year <- dt_results[get(tid_col) == v_year & !is.na(get(wcol))]
    dt_ws <- gd_sbw_bal_wstats_year(dt_year, wcol)
    dt_ws[, year := v_year]
    dt_ws
  }) %>% rbindlist

  dt_cors <- dt_long[, .(max_cor_raw = round(max(abs(cor_raw), na.rm = TRUE), 4),
                         max_cor_w = round(max(abs(cor_w), na.rm = TRUE), 4)),
                     by = year]
  merge(dt_wstats, dt_cors, by = "year")[order(year)]
}


gl_sbw_proc <- function(c_sbw, b_force = FALSE) {
  #' Run the SBW pipeline end-to-end for one config (the orchestrator).
  #' Stages: validate config -> build working panel (treatment lags/leads/
  #' d/dm) -> expand factor covariates to numeric dummies -> build covariate
  #' lags (n_lag_cov) -> select eligible years -> build formulas -> fit
  #' per-year weights (cached, parallel via fork) -> attach weights ->
  #' per-year balance diagnostics -> persist diag/config/registry.
  #'
  #' @param c_sbw   config list (see header for the fields and defaults)
  #' @param b_force TRUE to ignore the weight cache and refit all years
  #' @return list(dt_panel, dt_results, dt_years, dt_diag, config, cache):
  #'   dt_panel   = the input panel, untouched (== config$dt_panel)
  #'   dt_results = working panel + weights (the analysis table)
  #'   dt_years   = eligible years with complete-row counts
  #'   dt_diag    = per-year fit/balance diagnostics
  #'   config     = defaulted config (covariate_cols = the FINAL covariate
  #'                set: numeric + expanded dummies + surviving lags)
  #'   cache      = list(dir_cache, dir_run, hash, from_cache_all, elapsed_fit_s)
  
  c_sbw_mod <- md_sbw_defaults(c_sbw)
  gd_sbw_check(c_sbw_mod)

  ## working panel + factor-covariate expansion + covariate lags + eligibility
  dt_results <- md_sbw_panel(c_sbw_mod)
  l_exp <- md_sbw_expand_covs(dt_results, c_sbw_mod$covariate_cols)
  dt_results <- l_exp$dt
  c_sbw_mod$covariate_cols <- l_exp$covs   ## factors now referenced by dummies
  l_lags <- md_sbw_cov_lags(dt_results, c_sbw_mod)
  dt_results <- l_lags$dt

  c_sbw_mod$covariate_cols <- l_lags$covs  ## final covariate set incl. surviving lags
  ## gl_sbw_balcols(c_sbw_mod)

  dt_years <- gd_sbw_years(c_sbw_mod, dt_results)
  if (nrow(dt_years) == 0) {
    stop("no eligible fit years: check impact_years vs panel end, years, and min_n_year")
  }
  l_fs <- gf_sbw_formulas(c_sbw_mod)
  if (c_sbw_mod$verbose) {print(l_fs)}

  ## per-year weights (cached on disk), attached under config$wcol
  r_w <- gl_sbw_weights(c_sbw_mod, dt_results, dt_years, l_fs, b_force = b_force)

  ## per-year fit + balance diagnostics
  dt_bal  <- gd_sbw_bal_summary(dt_results, c_sbw_mod)
  ## gd_sbw_bal_long(dt_results, c_sbw_mod, v_years = 2012, wcol = "w_sbw")
  dt_diag <- merge(r_w$dt_fit, dt_bal, by = "year")
  setorderv(dt_diag, "year")

  if (any(dt_diag$ess < 0.05 * dt_diag$n_fit)) {
    warning("some years have ESS < 5% of n: weights are near-degenerate ",
            "(balance constraints barely feasible). Relax `tols` or simplify ",
            "the balance specification.")
  }

  ## persist diagnostics + config snapshot + one registry row (when fitting)
  fwrite(dt_diag, file.path(r_w$dir_run, "diag.csv"))
  c_spec <- c_sbw_mod
  c_spec$dt_panel <- NULL
  saveRDS(c_spec, file.path(r_w$dir_run, "config.rds"))
  
  if (!r_w$from_cache_all) {
    dt_reg <- data.table(label = c_sbw_mod$label, hash = r_w$hash,
                         dir_run = r_w$dir_run,
                         years = paste(dt_years$year, collapse = "|"),
                         n_years = nrow(dt_years), n_obs = sum(dt_diag$n_fit),
                         elapsed_s = r_w$elapsed_s,
                         time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    fwrite(dt_reg, file.path(c_sbw_mod$dir_cache, "registry.csv"), append = TRUE)
  }

  list(dt_panel = c_sbw_mod$dt_panel, dt_results = dt_results, dt_years = dt_years,
       dt_diag = dt_diag, config = c_sbw_mod,
       cache = list(dir_cache = c_sbw_mod$dir_cache, dir_run = r_w$dir_run,
                    hash = r_w$hash, from_cache_all = r_w$from_cache_all,
                    elapsed_fit_s = r_w$elapsed_s))
}


gd_sbw_registry <- function(dir_cache) {
  #' Read the SBW run registry (most recent run first).
  #'
  #' @param dir_cache cache root (config$dir_cache)
  #' @return data.table of runs; empty data.table if none recorded yet
  p <- file.path(path.expand(dir_cache), "registry.csv")
  if (!file.exists(p)) return(data.table())
  dt_reg <- fread(p)
  dt_reg[order(-seq_len(.N))]
}


gf_sbw_did <- function(c_sbw, l_tps = NULL) {
  #' DiD formula on the weighted panel: c(exn_X_d{k} / exn_X_dm{k} across all
  #' treatments and time points) ~ <same-year treatments> | <tid>.
  #'
  #' Time-point defaults are INDEPENDENT for the pre and post side, matching
  #' the now-independent n_lag and impact_years: pre periods -(2..n_lag)
  #' (l1 is the baseline, dm1 = 0 by construction and never included) and
  #' post periods 1..impact_years. Sorted ascending afterwards.
  #'
  #' @param c_sbw config (as stored in the result)
  #' @param l_tps optional integer vector of relative years (excludes 0)
  #' @return formula
  K       <- c_sbw$impact_years
  max_lag <- c_sbw$n_lag

  if (is.null(l_tps)) l_tps <- sort(c(-(2:max_lag), seq_len(K)))
  if (any(l_tps > K | l_tps < -max_lag | l_tps == 0)) {
    stop(sprintf("time points must be in [%d, %d] \\ {0}: post periods only go to impact_years = %d",
                 -max_lag, K, K))
  }

  l_dvs <- map(c_sbw$treatment_cols, function(treat_col) {
    sprintf("%s_%s%d", treat_col, fifelse(l_tps > 0, "d", "dm"), abs(l_tps))
  }) %>% unlist

  sprintf("c(%s) ~ %s | %s", paste(l_dvs, collapse = ", "),
          paste(c_sbw$treatment_cols, collapse = " + "), c_sbw$tid) %>% as.formula

}


gl_sbw_did <- function(r_sbw, l_tps = NULL, ...) {
  #' Pooled DiD on the SBW-weighted panel:
  #'   feols(c(<d/dm outcomes>) ~ <same-year treatments> | <tid>,
  #'         data = dt_results, weights = ~wcol, cluster = ~uid)
  #' The result is a fixest_multi object; extract coefficients with
  #' gd_sbw_coef_did (year fixed effects make this the meta-analytic average
  #' of the per-year fits).
  #'
  #' @param r_sbw result of gl_sbw_proc
  #' @param l_tps optional time points (gf_sbw_did)
  #' @param ...   further feols arguments (vcov, etc.)
  #' @return fixest_multi model object
  c_sbw <- r_sbw$config

  f_did <- gf_sbw_did(c_sbw, l_tps)
  
  feols(f_did, data = r_sbw$dt_results,
        weights = as.formula(sprintf("~ %s", c_sbw$wcol)),
        cluster = as.formula(sprintf("~ %s", c_sbw$uid)), ...)

}


gd_sbw_coef_did <- function(rx) {
  #' Extract coefficients from a fixest multi-LHS DiD model (gl_sbw_did).
  #'
  #' @param rx fixest_multi object (or any list-like with a tidy method)
  #' @return data.table with columns
  #'   (mdl, term, estimate, std.error, statistic, p.value, dv, tp)
  #'   where dv is the outcome treatment type (exn_pm, ...) and tp is the
  #'   relative year (negative = pre-period dm, positive = post-period d)
  dt_coef <- map(rx, tidy) %>% rbindlist(idcol = "mdl")

  dt_coef[, tp_prep := sub("^lhs: ", "", mdl)] %>%
    .[, `:=`(dv = sub("_(d|dm)[0-9]+$", "", tp_prep),
             tdir = fifelse(grepl("_dm", tp_prep), -1, 1),
             tlen = as.integer(sub("^.*_(d|dm)([0-9]+)$", "\\2", tp_prep)))] %>%
    .[, tp := tdir * tlen] %>%
    .[, `:=`(tdir = NULL, tlen = NULL)]

  return(dt_coef)
}
