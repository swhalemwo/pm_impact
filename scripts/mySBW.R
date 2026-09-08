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
##                 n_lag_cov      = 0L,                  # lags of numeric covariates
##                 impact_years   = 4L,                  # outcome window d1..dK
##                 years          = NULL)                # NULL -> auto (gd_sbw_years)
##
##   r_sbw <- gl_sbw_proc(c_sbw)    # fit weights (cached per year on disk)
##   r_did <- gl_sbw_did(r_sbw)     # pooled DiD on the weighted panel
##
## Stages (each a standalone function):
##   md_sbw_defaults / gd_sbw_check   fill defaults / validate config
##   md_sbw_panel                     lags + leads + d/dm outcome columns
##   gd_sbw_years                     year eligibility (coverage + min n)
##   gf_sbw_formulas                  optweightMV formula list (one per treatment)
##   gf_sbw_hash                      cache key (config spec + years + data content)
##   gl_sbw_weights                   per-year fits, cached as w_<year>.rds, PSOCK parallel
##   gd_sbw_bal_*                     balance report (weighted vs raw corr, ESS)
##   gl_sbw_proc                      orchestrator -> list(dt, dt_years, dt_diag, config, cache)
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
## - Factor/character covariates are auto-expanded to numeric dummies (all
##   levels minus the first sorted level = reference) — optweightMV is
##   infeasible when passed factors directly. Rare levels (< 0.5% or < 100
##   rows) warn: exact ATE balance with them is usually infeasible (recode).
## - Covariates that contain the current-year treatment (e.g. a cumulative
##   count as of year t) cannot be balanced exactly against that treatment —
##   pass the cumulative as of t-1, or rely on the treatment lags.
## - The osqp solver is single-threaded. n_cores > 1 runs years on a PSOCK
##   cluster (the fit function + constants are exported to the workers);
##   set OPENBLAS_NUM_THREADS = 1 (forking + OpenBLAS threads deadlocks).
## - eps = osqp solver tolerance (1e-3 ~2x faster than 1e-4, weights
##   essentially identical); tols = *balance* tolerance: NULL = exact (the
##   optweight default), a scalar = uniform tolerance, or a full
##   per-formula tolerance list.
## - estimand is fixed to "ATE" (optweightMV only supports ATE).
## - min_n_year only matters when years = NULL (auto): the minimum number of
##   complete rows a year needs to be kept. Explicitly requested years are
##   kept with as few as one complete row.
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
    n_lag_cov      = 0L,           # lags of numeric covariates to balance
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
  if (c_sbw$n_lag_cov > 0 && !all(map_lgl(l_covs, ~is.numeric(c_sbw$dt_panel[[.x]])))) {
    stop("n_lag_cov > 0 requires all covariate_cols to be numeric ",
         "(categorical covariates are balanced at level t only)")
  }
  if (c_sbw$wcol %in% names(c_sbw$dt_panel)) {
    stop(sprintf("wcol '%s' already exists in dt_panel", c_sbw$wcol))
  }

  ## numeric knobs
  l_ints <- c(n_lag = c_sbw$n_lag, n_lag_cov = c_sbw$n_lag_cov,
              impact_years = c_sbw$impact_years, min_n_year = c_sbw$min_n_year,
              n_cores = c_sbw$n_cores)
  
  if (!all(map_lgl(l_ints, ~is.integer(.x)))) {
    stop("n_lag / n_lag_cov / impact_years / min_n_year / n_cores must be integers")
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

## FIXME: make sure this handles different covariates (time-varying, categorical) properly
gl_sbw_balcols <- function(c_sbw) {
  #' The balance columns: everything each treatment formula balances on —
  #' lagged treatments (t-1..t-n_lag of every treatment type) followed by the
  #' covariates at year t and (for numeric covariates, n_lag_cov > 0) their
  #' lags. NOTE: returns column *names* (a character vector), not a formula.
  #'
  #' @param c_sbw defaulted config
  #' @return character vector of column names
  c(gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
    c_sbw$covariate_cols,
    if (c_sbw$n_lag_cov > 0) {
      gd_sbw_shiftcols(c_sbw$covariate_cols, c_sbw$n_lag_cov, "l")
    } else {
      character(0)
    })
}


gf_sbw_formulas <- function(c_sbw) {
  #' One optweightMV formula per treatment: exn_X ~ <all treatment lags> +
  #' <covariates>. Same-year other treatments stay out of the RHS so each
  #' fitted coefficient is a partial dose effect conditional on the observed
  #' other treatments of that year.
  #'
  #' @param c_sbw defaulted config (covariate_cols already expanded to
  #'   dummies by md_sbw_expand_covs when factors were involved)
  #' @return named list of formulas, names = treatment_cols

  l_rhs <- paste(c(gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
                   c_sbw$covariate_cols,
                   if (c_sbw$n_lag_cov > 0) {
                     gd_sbw_shiftcols(c_sbw$covariate_cols, c_sbw$n_lag_cov, "l")
                   }),
                 collapse = " + ")

  map(c_sbw$treatment_cols, ~as.formula(sprintf("%s ~ %s", .x, l_rhs))) %>%
    setNames(c_sbw$treatment_cols)
}


md_sbw_panel <- function(c_sbw) {
  #' Build the working panel from the base dt_panel (a lean copy; the input
  #' is never modified): keeps uid/tid/treatments/covariates and adds
  #'   - treatment lags exn_X_l1..l{max_lag}   (max_lag = max(n_lag, impact_years);
  #'     lags beyond n_lag serve as unbalanced placebo pre-periods)
  #'   - outcome columns exn_X_d{k}  = lead k - l1  (k = 1..impact_years)
  #'   - outcome columns exn_X_dm{k} = lag k - l1   (k = 2..L_lag; dm1 = 0 is
  #'     the DiD baseline by construction and is not stored)
  #' Intermediate lead columns are dropped again.
  #'
  #' @param c_sbw defaulted config
  #' @return data.table: the working panel
  l_treats <- c_sbw$treatment_cols
  l_covs   <- c_sbw$covariate_cols
  K        <- c_sbw$impact_years
  ## max_lag    <- max(c_sbw$n_lag, K)
  max_lag    <- c_sbw$n_lag
  

  ## generated names must not collide with base panel columns (the pipeline
  ## owns the *_l/*_f/*_d/*_dm namespace)
  l_new <- c(gd_sbw_shiftcols(l_treats, max_lag, "l"),  # lags
             gd_sbw_shiftcols(l_treats, K, "f"), # leads
             sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K)), # forward diffs
             if (max_lag > 1) sprintf("%s_dm%d", rep(l_treats, each = max_lag - 1), 2:max_lag) else character(0)) # backward diffs
  l_clash <- intersect(l_new, names(c_sbw$dt_panel))
  if (length(l_clash) > 0) {
    stop(sprintf("dt_panel already contains pipeline-generated columns: %s. ",
                 paste(l_clash, collapse = ", ")),
         "Pass the base panel (no precomputed lags/leads).")
  }

  dt_working <- copy(c_sbw$dt_panel[, c(c_sbw$uid, c_sbw$tid, l_treats, l_covs), with = FALSE])
  setorderv(dt_working, c(c_sbw$uid, c_sbw$tid))

  ## lags of treatments (to max_lag)
  l_lagcols <- gd_sbw_shiftcols(l_treats, max_lag, "l")

  ## dt_working[, (l_lagcols) := shift(.SD, 1:max_lag, type = "lag"), by = c(c_sbw$uid), .SDcols = l_treats]
  dt_working[,  (l_lagcols) := collapse::L(.SD, n = 1:max_lag, by = get(c_sbw$uid)), .SDcols = l_treats]
    

  ## covariate lags (numeric covariates only; enforced in gd_sbw_check)
  ## CHECK ME, not fully tested
  ## covariate lagging not sufficiently flexible me thinks: difference between time-invariant and categorical
  if (c_sbw$n_lag_cov > 0) {
    l_covlagcols <- gd_sbw_shiftcols(l_covs, c_sbw$n_lag_cov, "l")
    ## dt_working[, (l_covlagcols) := shift(.SD, 1:c_sbw$n_lag_cov, type = "lag"),
    ##            by = c(c_sbw$uid), .SDcols = l_covs]
    dt_working[, (l_covlagcols) := shift(.SD, n = 1:c_sbw$n_lag_cov, by = get(c_sbw$uid)), .SDcols = l_covs]
  }

  ## leads of treatments (intermediate)
  l_fcols <- gd_sbw_shiftcols(l_treats, K, "f")
  ## dt_working[, (l_fcols) := shift(.SD, 1:K, type = "lead"), by = c(c_sbw$uid), .SDcols = l_treats]
  dt_working[, (l_fcols) := collapse::L(.SD, n=-(1:K), by = get(c_sbw$uid)), .SDcols = l_treats]

  ## d{k} = lead{k} - l1
  l_l1cols <- gd_sbw_shiftcols(l_treats, 1, "l")
  l_dcols  <- sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K))
  dt_working[, (l_dcols) := map2(l_fcols, rep(l_l1cols, each = K), ~get(.x) - get(.y))]

  ## dm{k} = lag{k} - l1  (k = 2..L_lag)
  if (max_lag > 1) {
    l_dmcols  <- sprintf("%s_dm%d", rep(l_treats, each = max_lag - 1), 2:max_lag)
    ## source lags: every column except each treatment's _l1 (positions
    ## 1, 1+max_lag, ... in the col-major name order)
    l_dmsrcs  <- l_lagcols[-seq(1, length(l_lagcols), by = max_lag)]
    dt_working[, (l_dmcols) := map2(l_dmsrcs, rep(l_l1cols, each = max_lag - 1), ~get(.x) - get(.y))]
  }

  ## leads were intermediate -> drop
  ## dt_working[, (l_fcols) := NULL]
  return(dt_working)
}


sbw_dum_for_cov <- function(dt, cl) {
  #' dummy matrix for one factor/character covariate (one column per level
  #' minus the first, sorted level order = reference); column names are
  #' make.names-sanitized level labels so they stay syntactically valid
  x <- dt[[cl]]
  l_lv <- sort(unique(x[!is.na(x)]))
  if (length(l_lv) < 2) stop(sprintf("covariate '%s' has fewer than 2 levels", cl))
  m <- model.matrix(~ factor(x) - 1)   ## one column per level
  colnames(m) <- l_lv
  m <- m[, -1, drop = FALSE]           ## drop the (first sorted) reference level
  colnames(m) <- make.names(colnames(m))
  m
}


sbw_warn_rare <- function(dt, cl) {
  #' warn when a covariate level is tiny (exact balance with such dummies is
  #' typically infeasible per year); the user should recode before balancing
  l_rr <- names(which(table(dt[[cl]], useNA = "no") < max(100, 0.005 * nrow(dt))))
  if (length(l_rr) > 0) {
    warning(sprintf("covariate '%s' has rare level(s): %s. ", cl,
                    paste(l_rr, collapse = ", ")),
            "Exact ATE balance with such dummies is typically infeasible; ",
            "recode the covariate into a few meaningful levels before balancing.")
  }
}


md_sbw_expand_covs <- function(dt, l_covs) {
  #' expand factor/character covariates into numeric dummies (one per level
  #' minus the first, sorted level order = reference, consistent with
  #' model.matrix); numeric covariates stay as-is. optweightMV is infeasible
  #' when passed factor columns directly, hence the expansion.
  #' @return list(dt = dt with dummy columns added / factors removed, covs)

  l_keep <- l_covs[map_lgl(l_covs, ~is.numeric(dt[[.x]]))]
  l_exp <- setdiff(l_covs, l_keep)
  if (length(l_exp) == 0) return(list(dt = dt, covs = l_covs))

  l_dum <- map(l_exp, sbw_dum_for_cov, dt = dt) %>% setNames(l_exp)
  l_nms <- map2(names(l_dum), l_dum, ~sprintf("%s_%s", .x, colnames(.y))) %>% unlist
  if (any(l_nms %in% names(dt))) {
    stop(sprintf("dummy columns would clash with existing columns: %s",
                 paste(intersect(l_nms, names(dt)), collapse = ", ")))
  }
  walk(l_exp, sbw_warn_rare, dt = dt)

  m_dum <- do.call(cbind, unname(l_dum))
  colnames(m_dum) <- l_nms

  dt <- cbind(dt, data.table::as.data.table(m_dum))
  dt[, (l_exp) := NULL]                  ## factor info now lives in the dummies
  list(dt = dt, covs = c(l_keep, l_nms))
}


gd_sbw_years <- function(c_sbw, dt_working) {
  #' Eligible fit years: years with at least one complete row on all fit
  #' columns AND full outcome coverage (year + impact_years <= last panel
  #' year, so every fitted row can be matched with its d1..dK outcomes).
  #' With config years = NULL additionally require n_use >= min_n_year;
  #' explicitly requested years are honored with >= 1 complete row (requested
  #' but ineligible years are dropped with a warning).
  #'
  #' @param c_sbw      defaulted config
  #' @param dt_working working panel (md_sbw_panel output)
  #' @return data.table with columns (year, n_use) for the eligible years,
  #'   sorted ascending

  tid_col  <- c_sbw$tid
  K         <- c_sbw$impact_years
  l_fitcols <- c(c_sbw$treatment_cols, # treatment cols
                 gl_sbw_balcols(c_sbw), # past treatment cols
                 gd_sbw_shiftcols(c_sbw$treatment_cols, K, "f")) # future impact

  dt_years <- dt_working[, .(n_use = sum(complete.cases(.SD))), by = c(tid_col), .SDcols = l_fitcols] %>%
    .[n_use >= 1]

  setnames(dt_years, tid_col, "year")
  ## dt_years <- dt_years[n_use >= 1 & year <= max(dt_working[[tid_col]]) - K]

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
  #' @param c_sbw      defaulted config (covariate_cols expanded)
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

  l_cols <- unique(c(c_sbw$uid, c_sbw$tid, c_sbw$treatment_cols, c_sbw$covariate_cols))
  dt0 <- dt_working[, ..l_cols]
  data_str <- map_chr(dt0, function(x) {
    if (is.numeric(x)) digest::digest(round(x, 8), algo = "md5")
    else digest::digest(x, algo = "md5", serialize = FALSE)
  }) %>% paste(collapse = "")

  digest::digest(paste(spec_str, paste(dt_years$year, collapse = ","), data_str, sep = "|"),
                 algo = "md5", serialize = FALSE)

}


gd_sbw_fit_year <- function(dt_year, p_out, l_fs, l_fit_args, year, uid_col,
                            tid_col) {
  #' Fit optweightMV for one year's complete rows and write the weights to
  #' disk. Runs either in the master process or on a PSOCK worker. The
  #' weights file stores (uid, tid, w) per row, so nothing depends on order.
  #'
  #' @param dt_year    data.table of one year's rows: uid + tid + treatments +
  #'                   balance columns (complete cases only)
  #' @param p_out      full path of the w_<year>.rds file to write
  #' @param l_fs       named list of optweightMV formulas (gf_sbw_formulas)
  #' @param l_fit_args optweightMV arguments (estimand/norm/eps/tols/...)
  #' @param year       the calendar year (for error messages)
  #' @param uid_col    unit id column name (PeopleID)
  #' @param tid_col    time column name (year)
  #' @return list(n_fit, status, iter, run_time) solver diagnostics

  r <- do.call(optweight::optweightMV,
               c(list(formula.list = l_fs, data = dt_year), l_fit_args))
  if (!identical(r$info$status, "solved")) {
    stop(sprintf("optweightMV did not solve for year %s (status: %s). ",
                 year, r$info$status),
         "Try: relax `tols`, drop/recode covariates (factors are auto-expanded ",
         "to dummies; avoid covariates containing the current-year treatment, ",
         "e.g. cumulative counts as of year t), reduce n_lag, or trim extreme values.")
  }

  dt_w <- data.table::data.table(uid = dt_year[[uid_col]],
                                 tid = dt_year[[tid_col]],
                                 w = as.numeric(r$weights))
  data.table::setnames(dt_w, c("uid", "tid"), c(uid_col, tid_col))
  saveRDS(dt_w, p_out)

  list(n_fit = nrow(dt_year), status = r$info$status, iter = r$info$iter,
       run_time = r$info$run_time)
}


gd_sbw_year_slice <- function(dt_working, year, tid_col, uid_col, l_fitcols) {
  #' One year's rows of the working panel, keeping the row identity (uid,
  #' tid) plus the fit columns, complete rows only (optweightMV errors on any
  #' NA). The uid/tid columns travel into the fit so the saved weights carry
  #' their row identity.
  #'
  #' @param dt_working working panel (md_sbw_panel output)
  #' @param year       calendar year to slice
  #' @param tid_col    time column name (e.g. "year")
  #' @param uid_col    unit id column name (e.g. "PeopleID")
  #' @param l_fitcols  columns needed for the fit (treatments + balance cols)
  #' @return data.table with columns (uid, tid, fitcols), complete cases only
  v_year <- year   ## local scalar: `year` is also a column name in dt_working
  dt_year <- dt_working[get(tid_col) == v_year,
                        .SD, .SDcols = c(uid_col, tid_col, l_fitcols)]
  dt_year[complete.cases(dt_year[, .SD, .SDcols = l_fitcols])]
}


gd_sbw_read_year <- function(dt_year, dt_working, dir_run, year, tid_col,
                             uid_col, wcol) {
  #' Load one year's cached weights and attach them to the working panel rows
  #' of that year (matched on (uid, tid), never on row order). Validates that
  #' the cached weight count matches the number of complete fit rows.
  #'
  #' @param dt_year    complete-case rows of the year in dt_working (same
  #'                   definition as gd_sbw_year_slice returns)
  #' @param dt_working working panel to attach weights into (by reference)
  #' @param dir_run    cache directory holding w_<year>.rds
  #' @param year       calendar year
  #' @param tid_col,uid_col column names of time / unit id
  #' @param wcol       name of the weight column to fill
  #' @return number of fitted rows for this year (invisible)
  dt_w <- readRDS(file.path(dir_run, sprintf("w_%s.rds", year)))
  if (nrow(dt_w) != nrow(dt_year)) {
    stop(sprintf("cached weights for year %s have %d rows but %d complete fit rows",
                 year, nrow(dt_w), nrow(dt_year)))
  }
  if (!identical(sort(dt_w[[uid_col]]), sort(dt_year[[uid_col]])) ||
      !identical(sort(dt_w[[tid_col]]), sort(dt_year[[tid_col]]))) {
    stop(sprintf("cached (uid, tid) pairs for year %s do not match the current panel", year))
  }

  dt_w_keep <- dt_w[, .SD, .SDcols = c(uid_col, tid_col, "w")]
  dt_working[dt_w_keep, on = c(uid_col, tid_col), (wcol) := i.w]
  invisible(nrow(dt_w))
}


gl_sbw_weights <- function(c_sbw, dt_working, dt_years, l_fs, b_force = FALSE) {
  #' Per-year optweightMV fits with an on-disk cache (one w_<year>.rds per
  #' year, each storing the weights together with their (uid, tid) rows).
  #' Years are fitted in parallel over a PSOCK cluster when n_cores > 1.
  #' Afterwards every year's weights are attached to dt_working (by
  #' reference) under config$wcol; rows outside the fitted complete rows keep
  #' NA.
  #'
  #' @param c_sbw      defaulted + covariate-expanded config
  #' @param dt_working working panel (md_sbw_panel output)
  #' @param dt_years   eligible years (gd_sbw_years output)
  #' @param l_fs       formula list (gf_sbw_formulas output)
  #' @param b_force    TRUE: delete cached weights and refit all years
  #' @return list(dt_fit = per-year diagnostics, dir_run, hash,
  #'              from_cache_all, elapsed_s)

  uid_col  <- c_sbw$uid
  tid_col  <- c_sbw$tid
  v_years  <- dt_years$year
  l_fitcols <- c(c_sbw$treatment_cols, gl_sbw_balcols(c_sbw))

  ## cache location (content-addressed within the label folder)
  hash    <- gf_sbw_hash(c_sbw, dt_years, dt_working)
  dir_run <- file.path(c_sbw$dir_cache, c_sbw$label, substr(hash, 1, 12))
  dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)

  l_paths  <- setNames(file.path(dir_run, sprintf("w_%s.rds", v_years)), v_years)
  if (b_force) walk(l_paths, ~if (file.exists(.x)) file.remove(.x))
  l_cached <- map_lgl(l_paths, file.exists)

  t0 <- Sys.time()

  ## ---- fit the missing years ----
  v_miss <- v_years[!l_cached]
  if (length(v_miss) > 0) {

    l_year_dt <- map(v_miss, gd_sbw_year_slice, dt_working = dt_working,
                     tid_col = tid_col, uid_col = uid_col,
                     l_fitcols = l_fitcols) %>% setNames(v_miss)

    ## optweightMV arguments that are identical across years
    l_fit_args <- list(estimand = c_sbw$estimand, norm = c_sbw$norm,
                       min.w = c_sbw$min.w, verbose = c_sbw$verbose,
                       eps_rel = c_sbw$eps, eps_abs = c_sbw$eps,
                       tols.list = if (is.null(c_sbw$tols)) {
                         list(0)              ## exact balance (optweight default)
                       } else if (is.numeric(c_sbw$tols) && length(c_sbw$tols) == 1) {
                         rep(list(c_sbw$tols), length(l_fs))  ## uniform tolerance
                       } else {
                         c_sbw$tols          ## full per-formula tolerance list
                       })

    ## one job per year: the data + the output path (year is the list name)
    l_jobs <- map2(l_year_dt, l_paths[as.character(v_miss)],
                   ~list(dt_year = .x, p_out = .y))

    if (c_sbw$n_cores > 1) {
      ## PSOCK cluster across years; export the fit function and the
      ## constants so the workers do not need the package source
      Sys.setenv(OPENBLAS_NUM_THREADS = "1")  ## fork + OpenBLAS deadlocks
      cl <- parallel::makeCluster(c_sbw$n_cores, type = "PSOCK")
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c("gd_sbw_fit_year", "l_fs", "l_fit_args"),
                              envir = environment())
      l_r <- parallel::clusterMap(cl, function(job, year) {
        gd_sbw_fit_year(job$dt_year, job$p_out, l_fs, l_fit_args, year,
                        uid_col, tid_col)
      }, l_jobs, names(l_jobs))
    } else {
      ## sequential: fit years one after another in the master process
      l_r <- map2(l_jobs, v_miss, function(job, year) {
        gd_sbw_fit_year(job$dt_year, job$p_out, l_fs, l_fit_args, year,
                        uid_col, tid_col)
      })
    }

    ## debug: don't delete
    ## gd_sbw_fit_year(l_year_dt$`2012`, p_out = l_jobs$`2012`$p_out, l_fs = l_fs, l_fit_args = l_fit_args, year = 2012, uid_col =uid_col, tid_col = tid_col)
    
    names(l_r) <- v_miss

    dt_fit <- map2(v_miss, l_r, ~data.table(year = .x, n_fit = .y$n_fit,
                                            status = .y$status, iter = .y$iter,
                                            run_time = .y$run_time,
                                            from_cache = FALSE)) %>% rbindlist
  } else {
    dt_fit <- data.table()
  }

  ## ---- attach all years' weights (fresh or cached) to dt_working ----
  l_nfit <- map(v_years, function(year) {
    dt_year <- gd_sbw_year_slice(dt_working, year, tid_col, uid_col, l_fitcols)
    gd_sbw_read_year(dt_year, dt_working, dir_run, year, tid_col, uid_col,
                     c_sbw$wcol)
  }) %>% setNames(v_years)

  ## gd_sbw_year_slice(dt_working, 2012, tid_col, uid_col, l_fitcols)

  ## diagnostics rows for the years that came from cache
  if (any(l_cached)) {
    dt_fit_c <- map(v_years[l_cached], ~data.table(year = .x,
                                                   n_fit = l_nfit[[as.character(.x)]],
                                                   status = "cached", iter = NA_real_,
                                                   run_time = NA_real_, from_cache = TRUE)) %>%
      rbindlist
    dt_fit <- rbind(dt_fit, dt_fit_c)
  }
  setorderv(dt_fit, "year")

  list(dt_fit = dt_fit, dir_run = dir_run, hash = hash,
       from_cache_all = all(l_cached),
       elapsed_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2))
}


gd_sbw_cor_w <- function(m, w) {
  #' Weighted correlation matrix of the columns of m.
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
  #' @param dt_year   one year of the weighted panel (weighted rows only)
  #' @param l_treats  treatment column names
  #' @param l_balcols balance column names (gl_sbw_balcols output)
  #' @param wcol      weight column name
  #' @return data.table with columns (treat, balcol, cor_raw, cor_w)

  w <- dt_year[[wcol]]
  m <- gd_sbw_mat(dt_year, c(l_treats, l_balcols))
  n_t <- length(l_treats)
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


gd_sbw_bal_summary_year <- function(dt_results, year, tid_col, wcol, l_treats,
                                    l_balcols) {
  #' One-row balance summary (ESS, mean weight, max |corr| raw/weighted) for
  #' one weighted year of dt_results.
  #'
  #' @param dt_results weighted panel
  #' @param year       calendar year
  #' @param tid_col    time column name
  #' @param wcol       weight column name
  #' @param l_treats   treatment column names
  #' @param l_balcols  balance column names
  #' @return single-row data.table (year, ess, mean_w, max_cor_raw, max_cor_w)

  v_year <- year   ## local scalar: `year` is also a column name in dt_results
  dt_year <- dt_results[get(tid_col) == v_year & !is.na(get(wcol))]
  w    <- dt_year[[wcol]]
  dt_b <- gd_sbw_bal_year(dt_year, l_treats, l_balcols, wcol)

  data.table(year = v_year,
             ess = round(sum(w)^2 / sum(w^2), 1),
             mean_w = round(mean(w), 4),
             max_cor_raw = round(max(abs(dt_b$cor_raw), na.rm = TRUE), 4),
             max_cor_w = round(max(abs(dt_b$cor_w), na.rm = TRUE), 4))
}


gd_sbw_bal_summary <- function(dt_results, c_sbw, v_years = NULL, wcol = NULL) {
  #' Per-year balance summary over the weighted panel (used by gl_sbw_proc
  #' for the dt_diag table).
  #'
  #' @param dt_results weighted panel
  #' @param c_sbw      config
  #' @param v_years    optional subset; NULL = all weighted years
  #' @param wcol       weight column (defaults to config$wcol)
  #' @return data.table one row per year
  c_sbw <- md_sbw_defaults(c_sbw)
  wcol    <- wcol %||% c_sbw$wcol
  tid_col <- c_sbw$tid
  v_years <- v_years %||% dt_results[!is.na(get(wcol)), sort(unique(get(tid_col)))]

  ## DEBUG
  ## gd_sbw_bal_summary_year(dt_results, 2012, tid_col = tid_col, wcol = wcol, l_treats = c_sbw$treatment_cols, l_balcols = gl_sbw_balcols(c_sbw))

  map(v_years, gd_sbw_bal_summary_year, dt_results = dt_results,
      tid_col = tid_col, wcol = wcol, l_treats = c_sbw$treatment_cols,
      l_balcols = gl_sbw_balcols(c_sbw)) %>% rbindlist

}


gl_sbw_proc <- function(c_sbw, b_force = FALSE) {
  #' Run the SBW pipeline end-to-end for one config (the orchestrator).
  #' Stages: validate config -> build working panel -> expand factor
  #' covariates to dummies -> select eligible years -> build formulas ->
  #' fit per-year weights (cached) -> attach weights -> per-year balance
  #' diagnostics -> persist diag/config/registry.
  #'
  #' @param c_sbw   config list (see header for the fields and defaults)
  #' @param b_force TRUE to ignore the weight cache and refit all years
  #' @return list(dt_panel, dt_results, dt_years, dt_diag, config, cache):
  #'   dt_panel   = the input panel, untouched (== config$dt_panel)
  #'   dt_results = working panel + weights (the analysis table)
  #'   dt_years   = eligible years with complete-row counts
  #'   dt_diag    = per-year fit/balance diagnostics
  #'   config     = defaulted config (covariate_cols = expanded dummy names)
  #'   cache      = list(dir_cache, dir_run, hash, from_cache_all, elapsed_fit_s)
  
  c_sbw_mod <- md_sbw_defaults(c_sbw)
  gd_sbw_check(c_sbw_mod)

  ## working panel + factor-covariate expansion + eligibility
  dt_results <- md_sbw_panel(c_sbw_mod)
  l_exp <- md_sbw_expand_covs(dt_results, c_sbw_mod$covariate_cols)
  dt_results <- l_exp$dt
  c_sbw_mod$covariate_cols <- l_exp$covs   ## factors now referenced by dummies

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
  #' treatments and time points) ~ <same-year treatments> | <tid>. Default
  #' time points: pre -(2..L_lag) (l1 is the baseline, dm1 = 0 by
  #' construction) and post 1..impact_years.
  #'
  #' @param c_sbw config (as stored in the result)
  #' @param l_tps optional integer vector of relative years (excludes 0)
  #' @return formula
  K     <- c_sbw$impact_years
  L_lag <- max(c_sbw$n_lag, K)
  l_tps <- l_tps %||% c(-(2:L_lag), seq_len(K))
  if (any(l_tps > K | l_tps < -L_lag | l_tps == 0)) {
    stop(sprintf("time points must be in [%d, %d] \\ {0}: post periods only go to impact_years = %d",
                 -L_lag, K, K))
  }

  l_dvs <- map(c_sbw$treatment_cols, function(tr) {
    sprintf("%s_%s%d", tr, fifelse(l_tps > 0, "d", "dm"), abs(l_tps))
  }) %>% unlist

  sprintf("c(%s) ~ %s | %s", paste(l_dvs, collapse = ", "),
          paste(c_sbw$treatment_cols, collapse = " + "), c_sbw$tid) %>% as.formula
}


gl_sbw_did <- function(r_sbw, l_tps = NULL, ...) {
  #' pooled DiD on the SBW-weighted panel: feols(multi-LHS d/dm outcomes ~
  #' same-year treatments | tid, weights = ~wcol, cluster = ~uid)
  #' @return fixest_multi object; extract with gd_sbw_coef_did

  c_sbw <- r_sbw$config
  f <- gf_sbw_did(c_sbw, l_tps)
  feols(f, data = r_sbw$dt,
        weights = as.formula(sprintf("~ %s", c_sbw$wcol)),
        cluster = as.formula(sprintf("~ %s", c_sbw$uid)), ...)
}


gd_sbw_coef_did <- function(rx) {
  #' extract coefficients from a fixest multi-LHS DiD model (gl_sbw_did):
  #' columns mdl/term/estimate/std.error/statistic/p.value/dv/tp

  dt_coef <- map(rx, tidy) %>% rbindlist(idcol = "mdl")

  dt_coef[, tp_prep := sub("^lhs: ", "", mdl)] %>%
    .[, `:=`(dv = sub("_(d|dm)[0-9]+$", "", tp_prep),
             tdir = fifelse(grepl("_dm", tp_prep), -1, 1),
             tlen = as.integer(sub("^.*_(d|dm)([0-9]+)$", "\\2", tp_prep)))] %>%
    .[, tp := tdir * tlen] %>%
    .[, `:=`(tdir = NULL, tlen = NULL)]

  dt_coef
}
