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
##
## Cache/registry: weights land in <dir_cache>/<label>/<hash>/w_<year>.rds as
## each year finishes; whenever fitting happened a run row (label/hash/years/
## elapsed) is appended to <dir_cache>/registry.csv. Re-running an identical
## config (same spec + same panel content) restores the weights from disk.
##
## Notes / contracts:
## - dt_panel must be a COMPLETE per-artist grid (rows for every year between
##   the artist's first and last panel year) and must NOT contain precomputed
##   lag/lead columns (the pipeline builds them; clash -> error).
## - Covariates: numeric columns are balanced as-is; factor/character columns
##   are auto-expanded to one numeric dummy per level minus the first (sorted)
##   level — optweightMV is infeasible when passed factors directly. Beware
##   covariates that contain the current-year treatment (e.g. a cumulative
##   count as of year t): balancing corr(treat_t, cum_t) = 0 is structurally
##   infeasible — pass the cumulative as of t-1 instead, or the lag columns.
## - The solver is single-threaded: parallelize across years with a PSOCK
##   cluster (n_cores > 1); forking R + OpenBLAS threads deadlocks.
## - optweightMV only supports estimand = "ATE". eps is the osqp solver
##   tolerance (1e-3 ~2x faster than 1e-4, weights essentially identical);
##   tols is the *balance* tolerance: NULL = exact balance (optweight
##   default), a scalar = uniform tolerance, or a full per-formula list.
## - Every default below can be overridden per config entry; fields not in
##   the config fall back to md_sbw_defaults. min_n_year only matters when
##   years = NULL (auto-selection): it is the minimum number of complete rows
##   a year needs to be kept; explicitly requested years are kept with as few
##   as one complete row.
## - Moments/interactions are not hard-wired: precompute e.g. I(x^2) or
##   product columns and list them in covariate_cols.
## =====================================================================


md_sbw_defaults <- function(c_sbw) {
  #' fill default values for missing config fields; user-supplied entries
  #' always win (plain overlay, no recursion into list members like dt_panel)

  l_dflt <- list(
    covariate_cols = character(0),
    n_lag_cov = 0L,
    years = NULL,
    min_n_year = 10L,
    estimand = "ATE",
    norm = "l2",
    eps = 1e-3,
    min.w = 1e-8,
    verbose = FALSE,
    tols = NULL,
    n_cores = 1L,
    wcol = "w_sbw",
    dir_cache = "~/Dropbox/phd/papers/foster/misc/sbw_cache")

  for (nm in setdiff(names(l_dflt), names(c_sbw))) c_sbw[[nm]] <- l_dflt[[nm]]
  c_sbw$dir_cache <- path.expand(c_sbw$dir_cache)
  c_sbw
}


gd_sbw_check <- function(c_sbw) {
  #' validate an SBW config; returns TRUE invisibly, errors with a message otherwise

  if (!is.list(c_sbw) || is.null(c_sbw$dt_panel)) {
    stop("config must be a list with at least a dt_panel entry")
  }
  if (!grepl("^[A-Za-z0-9._-]+$", c_sbw$label %||% "")) {
    stop("label must match ^[A-Za-z0-9._-]+$ (it names the cache folder)")
  }
  if (!is.data.table(c_sbw$dt_panel)) stop("dt_panel must be a data.table")

  ## ids
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
  if (!all(map_lgl(l_ints, ~is.numeric(.x) && .x == round(.x)))) {
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
  #' shifted column names in data.table::shift(.SD) order: col-major, n inner
  paste0(rep(l_cols, each = n), "_", dir, rep(seq_len(n), times = length(l_cols)))
}


gf_sbw_balcols <- function(c_sbw) {
  #' columns entering the balance constraints (the RHS of every treatment formula)

  c(gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
    c_sbw$covariate_cols,
    if (c_sbw$n_lag_cov > 0) {
      gd_sbw_shiftcols(c_sbw$covariate_cols, c_sbw$n_lag_cov, "l")
    } else {
      character(0)
    })
}


gf_sbw_formulas <- function(c_sbw) {
  #' one optweightMV formula per treatment: exn_X ~ lags of all treatments +
  #' covariates (at t and lagged); same-year other treatments stay out of the
  #' RHS (partial dose effect interpretation)

  c_sbw <- md_sbw_defaults(c_sbw)
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
  #' build working panel (lean copy, original dt_panel untouched): keeps
  #' uid/tid/treats/covs and adds lags (to max(n_lag, impact_years)), plus
  #' DiD outcome columns exn_X_d{k} = f{k} - l1 and exn_X_dm{k} = l{k} - l1
  #' (k >= 2; l1 is the DiD baseline). Leads are intermediate and dropped.

  c_sbw <- md_sbw_defaults(c_sbw)
  l_treats <- c_sbw$treatment_cols
  l_covs <- c_sbw$covariate_cols
  K <- c_sbw$impact_years
  L_lag <- max(c_sbw$n_lag, K)

  ## generated column names must not collide with panel columns (the pipeline
  ## owns the *_l/*_d/*_dm namespace; pass the base panel without these)
  l_new <- c(gd_sbw_shiftcols(l_treats, L_lag, "l"),
             gd_sbw_shiftcols(l_treats, K, "f"),
             sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K)),
             if (L_lag > 1) sprintf("%s_dm%d", rep(l_treats, each = L_lag - 1), 2:L_lag) else character(0))
  l_clash <- intersect(l_new, names(c_sbw$dt_panel))
  if (length(l_clash) > 0) {
    stop(sprintf("dt_panel already contains pipeline-generated columns: %s. ",
                 paste(l_clash, collapse = ", ")),
         "Pass the base panel (no precomputed lags/leads).")
  }

  dt <- copy(c_sbw$dt_panel[, c(c_sbw$uid, c_sbw$tid, l_treats, l_covs), with = FALSE])
  setorderv(dt, c(c_sbw$uid, c_sbw$tid))

  ## lags of treatments (to L_lag; lags beyond n_lag serve as unbalanced placebo pre-periods)
  l_lagcols <- gd_sbw_shiftcols(l_treats, L_lag, "l")
  dt[, (l_lagcols) := shift(.SD, 1:L_lag, type = "lag"),
     by = c(c_sbw$uid), .SDcols = l_treats]

  ## covariate lags (numeric covariates only; enforced in gd_sbw_check)
  if (c_sbw$n_lag_cov > 0) {
    l_covlagcols <- gd_sbw_shiftcols(l_covs, c_sbw$n_lag_cov, "l")
    dt[, (l_covlagcols) := shift(.SD, 1:c_sbw$n_lag_cov, type = "lag"),
       by = c(c_sbw$uid), .SDcols = l_covs]
  }

  ## leads of treatments (intermediate)
  l_fcols <- gd_sbw_shiftcols(l_treats, K, "f")
  dt[, (l_fcols) := shift(.SD, 1:K, type = "lead"),
     by = c(c_sbw$uid), .SDcols = l_treats]

  ## d{k} = f{k} - l1 (k = 1..K)
  l_l1 <- gd_sbw_shiftcols(l_treats, 1, "l")
  l_dcols <- sprintf("%s_d%d", rep(l_treats, each = K), seq_len(K))
  dt[, (l_dcols) := map2(l_fcols, rep(l_l1, each = K), ~get(.x) - get(.y))]

  ## dm{k} = l{k} - l1 (k = 2..L_lag; dm1 = 0 would be the baseline)
  if (L_lag > 1) {
    l_dmcols <- sprintf("%s_dm%d", rep(l_treats, each = L_lag - 1), 2:L_lag)
    ## drop each treatment's _l1 entry (positions 1, 1+L_lag, ... in col-major order)
    l_dmsrcs <- l_lagcols[-seq(1, length(l_lagcols), by = L_lag)]
    dt[, (l_dmcols) := map2(l_dmsrcs, rep(l_l1, each = L_lag - 1), ~get(.x) - get(.y))]
  }

  ## leads are intermediate -> drop
  dt[, (l_fcols) := NULL]
  dt
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


gd_sbw_years <- function(c_sbw, dt) {
  #' eligible fit years: complete rows with all balance columns non-NA and
  #' outcome coverage t + impact_years <= max(panel year). With years = NULL
  #' additionally require n_use >= min_n_year (see header); explicitly
  #' requested years are honored whenever n_use >= 1 (warn about requested
  #' years that drop out).

  c_sbw <- md_sbw_defaults(c_sbw)
  v_tid <- c_sbw$tid
  l_use <- c(c_sbw$treatment_cols, gf_sbw_balcols(c_sbw))
  K <- c_sbw$impact_years

  dt_years <- dt[, .(n_use = sum(complete.cases(.SD))), by = c(v_tid), .SDcols = l_use]
  setnames(dt_years, v_tid, "year")
  dt_years <- dt_years[n_use >= 1 & year <= max(dt[[v_tid]]) - K]

  v_years <- c_sbw$years
  if (length(v_years) > 0) {
    l_miss <- setdiff(v_years, dt_years$year)
    if (length(l_miss) > 0) {
      warning(sprintf("gd_sbw_years: requested years not eligible (no complete rows or no outcome coverage): %s",
                      paste(l_miss, collapse = ", ")))
    }
    dt_years <- dt_years[year %in% v_years]
  } else {
    dt_years <- dt_years[n_use >= c_sbw$min_n_year]
  }
  dt_years[order(year)]
}


gf_sbw_hash <- function(c_sbw, dt_years, dt) {
  #' md5 cache key: config spec (minus dt_panel/label) + selected years +
  #' content of the used working-panel columns (numerics rounded to kill fp
  #' noise). dt is the pipeline-built working panel (lags/leads derived from
  #' it, so its content determines the weights).

  c_spec <- c_sbw
  c_spec$dt_panel <- NULL
  c_spec$label <- NULL
  l_spec <- c_spec[order(names(c_spec))]
  spec_str <- map_chr(names(l_spec), ~sprintf("%s=%s", .x, paste(l_spec[[.x]], collapse = "|"))) %>%
    paste(collapse = ";")

  l_cols <- unique(c(c_sbw$uid, c_sbw$tid, c_sbw$treatment_cols, c_sbw$covariate_cols))
  dt0 <- dt[, ..l_cols]
  data_str <- map_chr(dt0, function(x) {
    if (is.numeric(x)) digest::digest(round(x, 8), algo = "md5")
    else digest::digest(x, algo = "md5", serialize = FALSE)
  }) %>% paste(collapse = "")

  digest::digest(paste(spec_str, paste(dt_years$year, collapse = ","), data_str, sep = "|"),
                 algo = "md5", serialize = FALSE)
}


sbw_fit_year <- function(dty, p_out, l_fs, l_fit_args, y) {
  #' fit one year of optweightMV, write weights to p_out (runs on PSOCK workers
  #' or in the master process); returns solver diagnostics

  r <- do.call(optweight::optweightMV, c(list(formula.list = l_fs, data = dty), l_fit_args))
  if (!identical(r$info$status, "solved")) {
    stop(sprintf("optweightMV did not solve for year %s (status: %s). ",
                 y, r$info$status),
         "Try: relax `tols`, drop/recode covariates (factors are auto-expanded ",
         "to dummies; avoid covariates containing the current-year treatment, ",
         "e.g. cumulative counts as of year t), reduce n_lag, or trim extreme values.")
  }
  w <- as.numeric(r$weights)
  saveRDS(w, p_out)
  list(n_fit = nrow(dty), status = r$info$status, iter = r$info$iter,
       run_time = r$info$run_time)
}


sbw_year_slice <- function(dt, y, v_tid, l_use) {
  #' rows of the working panel for one year with complete balance columns
  #' (keeps .ridx so weights can be re-attached to the full panel later)
  dty <- dt[get(v_tid) == y, .SD, .SDcols = c(".ridx", l_use)]
  dty[complete.cases(dty)]
}


sbw_load_year <- function(dt, y, dir_run, l_use, v_tid) {
  #' load one year's weights from disk and pair them with the .ridx of the
  #' complete rows they belong to; validates length vs complete rows
  dty <- dt[get(v_tid) == y, .SD, .SDcols = c(".ridx", l_use)]
  dty <- dty[complete.cases(dty)]
  w <- readRDS(file.path(dir_run, sprintf("w_%s.rds", y)))
  if (length(w) != nrow(dty)) {
    stop(sprintf("cached weights for year %s have length %d but %d complete rows",
                 y, length(w), nrow(dty)))
  }
  list(dt_w = data.table(.ridx = dty$.ridx, w = w), n_fit = length(w))
}


gd_sbw_attach_w <- function(dt, v_years, dir_run, l_use, v_tid, wcol) {
  #' load per-year weights from disk and attach them to dt (by reference, via
  #' the .ridx column)
  #' @return named vector of per-year n_fit

  l_att <- map(v_years, sbw_load_year, dt = dt, dir_run = dir_run,
               l_use = l_use, v_tid = v_tid) %>% setNames(v_years)

  dt[, (wcol) := NA_real_]
  dt[data.table::rbindlist(map(l_att, ~.x$dt_w)), on = ".ridx", (wcol) := i.w]
  map_dbl(l_att, ~.x$n_fit)
}


gl_sbw_weights <- function(c_sbw, dt, dt_years, l_fs, b_force = FALSE) {
  #' per-year optweightMV fits with on-disk cache (w_<year>.rds per year);
  #' attaches the weights column to dt (by reference, requires a .ridx column);
  #' parallelizes across years with a PSOCK cluster when n_cores > 1.
  #' @return list(dt_fit, dir_run, hash, from_cache_all, elapsed_s)

  c_sbw <- md_sbw_defaults(c_sbw)
  v_tid <- c_sbw$tid
  v_years <- dt_years$year
  l_use <- c(c_sbw$treatment_cols, gf_sbw_balcols(c_sbw))

  ## cache location
  hash <- gf_sbw_hash(c_sbw, dt_years, dt)
  dir_run <- file.path(c_sbw$dir_cache, c_sbw$label, substr(hash, 1, 12))
  dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)

  l_paths <- setNames(file.path(dir_run, sprintf("w_%s.rds", v_years)), v_years)
  if (b_force) walk(l_paths, ~if (file.exists(.x)) file.remove(.x))
  l_cached <- map_lgl(l_paths, file.exists)

  t0 <- Sys.time()

  ## fit the missing years (parallel across years, PSOCK only)
  v_miss <- v_years[!l_cached]
  if (length(v_miss) > 0) {
    l_ycc <- map(v_miss, sbw_year_slice, dt = dt, v_tid = v_tid, l_use = l_use) %>%
      setNames(v_miss)

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

    if (c_sbw$n_cores > 1) {
      Sys.setenv(OPENBLAS_NUM_THREADS = "1") ## fork + OpenBLAS threads deadlocks
      cl <- parallel::makeCluster(c_sbw$n_cores, type = "PSOCK")
      on.exit(parallel::stopCluster(cl))
      l_r <- parallel::clusterMap(cl, sbw_fit_year, l_ycc, l_paths[as.character(v_miss)],
                                  MoreArgs = list(l_fs = l_fs, l_fit_args = l_fit_args,
                                                  y = as.character(v_miss)))
    } else {
      l_r <- map2(l_ycc, l_paths[as.character(v_miss)], sbw_fit_year,
                  l_fs = l_fs, l_fit_args = l_fit_args, y = as.character(v_miss))
    }
    names(l_r) <- v_miss
    dt_fit <- map2(v_miss, l_r, ~data.table(year = .x, n_fit = .y$n_fit,
                                            status = .y$status, iter = .y$iter,
                                            run_time = .y$run_time,
                                            from_cache = FALSE)) %>% rbindlist
  } else {
    dt_fit <- data.table()
  }

  ## load all weights (fresh or cached) and attach them to dt
  l_nfit <- gd_sbw_attach_w(dt, v_years, dir_run, l_use, v_tid, c_sbw$wcol)

  ## diag rows for the cached years
  if (any(l_cached)) {
    dt_fit_c <- map(v_years[l_cached], ~data.table(year = .x, n_fit = l_nfit[as.character(.x)],
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
  #' weighted correlation matrix; NA where a column has (near-)zero variance

  sw <- sum(w)
  mw <- colSums(m * w) / sw
  mc <- m - rep(mw, each = nrow(m))
  covm <- crossprod(mc * w, mc) / sw
  sds <- sqrt(diag(covm))
  r <- covm / outer(sds, sds)
  r[sds < 1e-12, ] <- NA
  r[, sds < 1e-12] <- NA
  r
}


sbw_mat_col <- function(dty, cl) {
  #' matrix form of one column: numeric as-is, factor/character expanded to
  #' one dummy per level (names sanitized like md_sbw_expand_covs)
  x <- dty[[cl]]
  if (is.numeric(x)) {
    m1 <- as.matrix(x)
    colnames(m1) <- cl
    m1
  } else {
    m1 <- model.matrix(~ factor(x) - 1)
    colnames(m1) <- make.names(paste0(cl, "_", colnames(m1)))
    m1
  }
}


gd_sbw_mat <- function(dty, l_cols) {
  #' cbind matrix form of the given columns (used for the balance diagnostics;
  #' the fit itself already receives pre-expanded dummies from md_sbw_expand_covs)

  l_m <- map(l_cols, sbw_mat_col, dty = dty)
  do.call(cbind, l_m)
}


gd_sbw_bal_year <- function(dty, l_treats, l_balcols, wcol) {
  #' weighted vs raw correlations (treatment x balance column) for one year

  w <- dty[[wcol]]
  m <- gd_sbw_mat(dty, c(l_treats, l_balcols))
  n_t <- length(l_treats)
  r_w <- gd_sbw_cor_w(m, w)
  r_raw <- gd_sbw_cor_w(m, rep(1, nrow(m)))

  data.table(treat = rep(l_treats, times = ncol(m) - n_t),
             balcol = rep(colnames(m)[-seq_len(n_t)], each = n_t),
             cor_raw = as.vector(r_raw[seq_len(n_t), -seq_len(n_t), drop = FALSE]),
             cor_w = as.vector(r_w[seq_len(n_t), -seq_len(n_t), drop = FALSE]))
}


sbw_bal_long_year <- function(dt, y, v_tid, wcol, l_treats, l_bal) {
  #' balance rows (treat x balcol correlations) for one weighted year
  dty <- dt[get(v_tid) == y & !is.na(get(wcol))]
  dt_b <- gd_sbw_bal_year(dty, l_treats, l_bal, wcol)
  dt_b[, year := y]
  dt_b
}


gd_sbw_bal_long <- function(dt, c_sbw, v_years = NULL, wcol = NULL) {
  #' long balance table across years: year x treat x balcol with raw/weighted
  #' correlations (weighted rows only; v_years defaults to all weighted years)

  c_sbw <- md_sbw_defaults(c_sbw)
  wcol <- wcol %||% c_sbw$wcol
  v_tid <- c_sbw$tid
  l_treats <- c_sbw$treatment_cols
  l_bal <- gf_sbw_balcols(c_sbw)
  v_years <- v_years %||% dt[!is.na(get(wcol)), sort(unique(get(v_tid)))]

  map(v_years, sbw_bal_long_year, dt = dt, v_tid = v_tid, wcol = wcol,
      l_treats = l_treats, l_bal = l_bal) %>% rbindlist
}


sbw_bal_summary_year <- function(dt, y, v_tid, wcol, l_treats, l_bal) {
  #' one-row balance summary (ESS, mean weight, max |corr| raw/weighted) for
  #' one weighted year
  dty <- dt[get(v_tid) == y & !is.na(get(wcol))]
  w <- dty[[wcol]]
  dt_b <- gd_sbw_bal_year(dty, l_treats, l_bal, wcol)
  data.table(year = y,
             ess = round(sum(w)^2 / sum(w^2), 1),
             mean_w = round(mean(w), 4),
             max_cor_raw = round(max(abs(dt_b$cor_raw), na.rm = TRUE), 4),
             max_cor_w = round(max(abs(dt_b$cor_w), na.rm = TRUE), 4))
}


gd_sbw_bal_summary <- function(dt, c_sbw, v_years = NULL, wcol = NULL) {
  #' per-year balance summary (ESS, mean weight, max |corr| raw/weighted) for
  #' the weighted rows; used by gl_sbw_proc for the diag table

  c_sbw <- md_sbw_defaults(c_sbw)
  wcol <- wcol %||% c_sbw$wcol
  v_tid <- c_sbw$tid
  l_treats <- c_sbw$treatment_cols
  l_bal <- gf_sbw_balcols(c_sbw)
  v_years <- v_years %||% dt[!is.na(get(wcol)), sort(unique(get(v_tid)))]

  map(v_years, sbw_bal_summary_year, dt = dt, v_tid = v_tid, wcol = wcol,
      l_treats = l_treats, l_bal = l_bal) %>% rbindlist
}


gl_sbw_proc <- function(c_sbw, b_force = FALSE) {
  #' run the SBW pipeline end-to-end for a config
  #' @return list(dt, dt_years, dt_diag, config, cache): dt is the working
  #'   panel with the weights column; dt_panel (untouched) is in config$dt_panel

  c_sbw <- md_sbw_defaults(c_sbw)
  gd_sbw_check(c_sbw)

  ## working panel + factor-covariate expansion + per-year eligibility
  dt <- md_sbw_panel(c_sbw)
  l_exp <- md_sbw_expand_covs(dt, c_sbw$covariate_cols)
  dt <- l_exp$dt
  c_sbw$covariate_cols <- l_exp$covs   ## factors now referenced via their dummies
  dt[, .ridx := .I]
  dt_years <- gd_sbw_years(c_sbw, dt)
  if (nrow(dt_years) == 0) {
    stop("no eligible fit years: check impact_years vs panel end, years, and min_n_year")
  }
  l_fs <- gf_sbw_formulas(c_sbw)

  ## per-year weights (cached on disk)
  r_w <- gl_sbw_weights(c_sbw, dt, dt_years, l_fs, b_force = b_force)
  dt[, .ridx := NULL]

  ## per-year fit + balance diag
  dt_bal <- gd_sbw_bal_summary(dt, c_sbw)
  dt_diag <- merge(r_w$dt_fit, dt_bal, by = "year")
  setorderv(dt_diag, "year")
  if (any(dt_diag$ess < 0.05 * dt_diag$n_fit)) {
    warning("some years have ESS < 5% of n: weights are near-degenerate ",
            "(balance constraints barely feasible). Relax `tols` or simplify ",
            "the balance specification.")
  }

  ## persist: diag + config snapshot + registry row (when any fitting happened)
  fwrite(dt_diag, file.path(r_w$dir_run, "diag.csv"))
  c_spec <- c_sbw
  c_spec$dt_panel <- NULL
  saveRDS(c_spec, file.path(r_w$dir_run, "config.rds"))
  if (!r_w$from_cache_all) {
    dt_reg <- data.table(label = c_sbw$label, hash = r_w$hash, dir_run = r_w$dir_run,
                         years = paste(dt_years$year, collapse = "|"),
                         n_years = nrow(dt_years), n_obs = sum(dt_diag$n_fit),
                         elapsed_s = r_w$elapsed_s,
                         time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    fwrite(dt_reg, file.path(c_sbw$dir_cache, "registry.csv"), append = TRUE)
  }

  list(dt = dt, dt_panel = c_sbw$dt_panel, dt_years = dt_years, dt_diag = dt_diag,
       config = c_sbw,
       cache = list(dir_cache = c_sbw$dir_cache, dir_run = r_w$dir_run, hash = r_w$hash,
                    from_cache_all = r_w$from_cache_all, elapsed_fit_s = r_w$elapsed_s))
}


gd_sbw_registry <- function(dir_cache) {
  #' read the SBW run registry (most recent run first); empty data.table if none

  p <- file.path(path.expand(dir_cache), "registry.csv")
  if (!file.exists(p)) return(data.table())
  dt <- fread(p)
  dt[order(-seq_len(.N))]
}


gf_sbw_did <- function(c_sbw, l_tps = NULL) {
  #' DiD formula on the weighted panel: c(<exn_X_d{k} / exn_X_dm{k}> across
  #' treatments) ~ same-year treatments | tid. Default time points: pre
  #' -(2..L_lag) (l1 is the baseline, dm1 = 0 by construction) and post 1..K.

  c_sbw <- md_sbw_defaults(c_sbw)
  K <- c_sbw$impact_years
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
