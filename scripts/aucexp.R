## =====================================================================
## scripts/scrap/aucexp.R — euro-scale translation of SBW DiD coefficients
##
## The auction DiD is estimated on a LOG scale (auc_turnover_log =
## log1p(auc_turnover), aucturn_cum2_log = log1p(aucturn_cum2)), because
## raw turnover cannot be balanced (no arithmetic overlap between the PM
## and non-PM arms: mean turnover at t-1 is EUR 6.4k vs EUR 299k, top 100
## of 211k non-PM rows hold 58% of that arm's turnover -> optweightMV is
## primal-infeasible on raw counts). Log coefficients are hard to read,
## so this file converts them to euro amounts with standard errors.
## Nothing here re-estimates anything; it post-processes the coefficient
## table.
##
## ---------------------------------------------------------------------
## What the conversion is
##
## With log1p as the link, a DiD coefficient b_j on treatment j is a shift
## of the log(1 + y) scale:
##     log1p(y)^(1) = log1p(y)^(0) + b_j
## i.e. a PROPORTIONAL rescaling of the artist's own turnover plus one.
## Back-transforming with expm1 is exact and gives a clean closed form,
##
##     expm1(log1p(y) + b) - y = (1 + y) * (exp(b) - 1)
##
## so the (weighted) average marginal effect in euros is
##     euro_j = (exp(b_j) - 1) * (1 + mean_w(y))
##     mean_w(y) = weighted mean of the LEVEL y at the outcome horizon
## and its delta-method variance is
##     Var(euro_j) = (exp(b_j) * (1 + mean_w(y)))^2 * se(b_j)^2
##                 + (exp(b_j) - 1)^2 * Var(mean_w(y))
## (see "Standard errors" below — the second term is not negligible).
##
## The ranking of the four institution types, the ratios between them and
## the growth rates exp(b) - 1 do NOT depend on any of the choices this
## file makes; only the euro level does.
##
## Pre-period (placebo) coefficients are converted with the same formula,
## anchored at their OWN level: dm{k} = v(t-k) - v(t-1) is a backward
## change, so b_ref = "t+k" resolves to {dv}_l{k} rather than {dv}_f{k}
## (gd_auc_lvl_col). The resulting euro figures are diagnostics -- a
## well-behaved placebo should be ~0 euro as well -- and must NOT be read
## as effects, since the exhibition is dated after the outcome.
##
## ---------------------------------------------------------------------
## Standard errors — do NOT use the b-only delta method
##
## Var(mean_w(y)) is the DOMINANT term for auction turnover: the
## reference level is a weighted mean of a variable where the top 1% of
## rows hold ~97% of total turnover, so its own sampling variance is
## large.  Measured on the full panel (k = 6, exn_pm): b-only se = 592,
## se including the reference variance = 840, artist-cluster bootstrap
## (B = 300) = 918.  The b-only version understates by 30-70% depending
## on the treatment.  Var(mean_w(y)) is CR0 cluster-robust on
## c_sbw$uid; the covariance between b and the reference moment is
## ignored, which leaves the SE ~10% below the cluster bootstrap.
## se_euro_b is still returned so the difference stays visible.
##
## ---------------------------------------------------------------------
## Residuals — the reference level is the OBSERVED one
##
## The reference-level column comes from the pipeline as log1p(y) and is
## back-transformed with expm1() (gd_auc_lvl_tr, keyed off the dv suffix),
## so the euro level used below is the observed turnover at that horizon.
## expm1(log1p(y) + b) then evaluates the shift at the artist's ACTUAL level.
## Equivalently, the regression residual is kept: u(1) = u(0), i.e. the
## treatment moves each artist along the log scale without changing her
## position relative to other artists. This reproduces the data exactly
## when b = 0, and it is the version reported here.
##
## The alternative is to evaluate the shift at the model's FITTED level
## instead (u = 0, the marginaleffects-style "average counterfactual
## prediction": every artist sits exactly on the fitted line). That is a
## different, smaller estimand -- on this outcome the two differ by a
## factor of ~1.5, and 99% of the gap comes from the top 4.6% of rows, so
## it is a tail/Jensen effect, not a zeros effect. To get it, replace
## ref_level in gd_auc_euro_calc by the weighted mean of expm1(fitted log
## level), where the fitted log level is the model prediction (needs
## r_did); report both as bounds rather than silently picking one.
##
## ---------------------------------------------------------------------
## Caveats to state when reporting
##
## 1. Reference level.  b_ref = "t+k" (default) uses the level of y at
##    the outcome horizon, which is where the percentage applies.  That
##    level already contains the estimated treatment effect: measured
##    contamination is 7-9% of the level.  b_ref = "t-1" uses the
##    pre-treatment level; it answers a different question ("what is a b%
##    increase worth on today's turnover?") and is smaller because
##    turnover levels drift upward within the window.
## 2. Tail.  ~97% of total turnover sits in the top 1% of rows, so a euro
##    figure is an average over a very skewed distribution.  Never
##    multiply a per-artist-year figure by the number of artist-years to
##    get a "total market effect".
## 3. The log1p-proportional model cannot represent a nonseller becoming
##    a seller (log1p(0) = 0, so a zero row moves to exp(b) - 1 euro,
##    i.e. almost nothing).  These euro figures therefore reflect the
##    intensive margin; the participation channel has to be added
##    separately from an auc_any LPM (delta P in pp * level among prior
##    sellers).
## 4. The additive-shift-on-a-log-scale is not an effect on E[y]: the
##    choice in "Residuals" is exactly the retransformation problem
##    (Mullahy 1999, Manning & Mullahy 2001).  Duan (1983) smearing is
##    the standard fix and is NOT available here: the smearing factor is
##    mean_w(exp(u)) = 1,989 (one artist contributes 10% of it), because
##    the residual is a Pareto-like tail rather than a tame error
##    distribution.
##
## ---------------------------------------------------------------------
## Usage
##
##   source("scripts/scrap/aucexp.R")
##   dt_euro <- gd_auc_euro(r_sbw_aucexp, r_did_aucexp, dt_coef_aucexp)
##
## Depends on data.table, purrr, magrittr (as mySBW.R does) and on
## `%||%` (base >= 4.4).
## =====================================================================


gc_auc_euro_cfg <- function(r_sbw, b_ref = "t+k") {
  #' Bundle the shared settings of the euro conversion into one config
  #' list (coding_style.md §4.3), so the stage functions take two
  #' arguments instead of six.
  #'
  #' @param r_sbw       result of gl_sbw_proc
  #' @param b_ref       reference level for the proportional rescaling:
  #'                    "t+k" (level at the outcome horizon, default),
  #'                    "t-1" (pre-treatment level), "t", or a number
  #' @return config list (dt, wcol, uid, b_ref)
  if (is.null(r_sbw$dt_results) || is.null(r_sbw$config)) {
    stop("gc_auc_euro_cfg: r_sbw must be the result of gl_sbw_proc ",
         "(needs $dt_results and $config)")
  }
  if (!is.data.frame(r_sbw$dt_results)) {
    stop("gc_auc_euro_cfg: r_sbw$dt_results is not a data.frame")
  }
  if (!is.numeric(b_ref)) b_ref <- match.arg(b_ref, c("t+k", "t-1", "t"))

  ## the weight column and the cluster id drive the whole estimation, so
  ## their absence is a config error, not a data problem -> check now
  c_eu <- list(dt    = r_sbw$dt_results,
               wcol  = r_sbw$config$wcol %||% "w_sbw",
               uid   = r_sbw$config$uid,
               b_ref = b_ref)
  if (!c_eu$wcol %in% names(c_eu$dt)) {
    stop(sprintf("gc_auc_euro_cfg: weight column '%s' not in dt_results", c_eu$wcol))
  }
  if (is.null(c_eu$uid) || !c_eu$uid %in% names(c_eu$dt)) {
    stop("gc_auc_euro_cfg: config has no usable uid column (needed for the ",
         "cluster-robust reference variance)")
  }
  c_eu
}


md_auc_coef_get <- function(r_did, dt_coef = NULL) {
  #' Resolve and validate the coefficient table: use dt_coef when given,
  #' otherwise rebuild it from the DiD result.
  #'
  #' @param r_did   result of gl_sbw_did (may be NULL when dt_coef is given)
  #' @param dt_coef coefficient table from gd_sbw_coef_did, or NULL
  #' @return data.table with at least (dv, tp, term, estimate, std.error)
  if (is.null(dt_coef)) {
    if (is.null(r_did)) {
      stop("md_auc_coef_get: pass either dt_coef (from gd_sbw_coef_did) or r_did")
    }
    dt_coef <- gd_sbw_coef_did(r_did)
  }
  dt_coef <- as.data.table(dt_coef)

  ## the four columns the conversion consumes; fail loudly rather than
  ## producing NAs deeper in the calculation
  l_miss <- setdiff(c("dv", "tp", "term", "estimate", "std.error"), names(dt_coef))
  if (length(l_miss) > 0) {
    stop(sprintf("md_auc_coef_get: dt_coef is missing column(s): %s",
                 paste(l_miss, collapse = ", ")))
  }
  dt_coef[, tp := as.integer(tp)]
  dt_coef
}


md_auc_sel <- function(dt_coef, c_sbw_did, l_dv = NULL, l_tp = NULL) {
  #' Select the rows to convert: TREATMENT coefficients, pre and post
  #' period. no_treatment rows are levels rather than effects and have no
  #' euro interpretation, so they are dropped. Pre-period dm rows are kept
  #' because their euro value is a useful diagnostic (a well-behaved
  #' placebo should be ~0 euro too); they must still NOT be read as
  #' effects, since the "treatment" is in the past.
  #'
  #' @param dt_coef   validated coefficient table
  #' @param c_sbw_did SBW config (for the default dv set)
  #' @param l_dv      outcomes to convert; NULL = the config's auction
  #'                  outcomes that were actually estimated
  #' @param l_tp      relative years; NULL = every period present (pre and
  #'                  post). Pass e.g. -4:-2 for the placebo rows only.
  #' @return data.table (copy), one row per (dv, tp, term)
  l_cfg <- unique(c(c_sbw_did$outcome_cols, c_sbw_did$outcome_cols_cur))

  ## default dv set = the config's outcomes that actually appear in the
  ## coefficient table (gl_sbw_did may have been called on a subset)
  if (is.null(l_dv)) {
    if (length(l_cfg) == 0) {
      stop("md_auc_sel: config has no outcome_cols/outcome_cols_cur; pass l_dv ",
           "(a euro conversion only makes sense for the auction outcomes)")
    }
    l_dv <- intersect(l_cfg, unique(dt_coef$dv))
    if (length(l_dv) == 0) {
      stop(sprintf(paste0("md_auc_sel: none of the config's outcome columns (%s) is ",
                          "in the coefficient table (has: %s); pass l_dv"),
                   paste(l_cfg, collapse = ", "),
                   paste(unique(dt_coef$dv), collapse = ", ")))
    }
  }
  l_bad <- setdiff(l_dv, unique(dt_coef$dv))
  if (length(l_bad) > 0) {
    stop(sprintf("md_auc_sel: dv(s) not in the coefficient table: %s",
                 paste(l_bad, collapse = ", ")))
  }

  ## keep every treatment row, pre and post period alike (tp < 0 are the
  ## dm placebos; they get the same conversion, anchored at their own
  ## level -- see gd_auc_lvl_col); optionally restricted to l_tp
  dt_sel <- dt_coef[dv %in% l_dv & term != "no_treatment"]
  if (!is.null(l_tp)) dt_sel <- dt_sel[tp %in% as.integer(l_tp)]
  if (nrow(dt_sel) == 0) stop("md_auc_sel: nothing left to convert after filtering")
  copy(dt_sel)
}


gd_auc_lvl_col <- function(dv, tp, b_ref) {
  #' Resolve the panel column holding the reference LEVEL of `dv` at
  #' relative year `tp`. The pipeline builds the same shift family for
  #' every dv, so the name depends only on the sign of `tp` and on `b_ref`,
  #' never on the outcome family (d{k} = v(t+k) - v(t-1) for outcome_cols,
  #' v(t+k) - v(t) for outcome_cols_cur, but the denominator is v(t+k)
  #' either way).
  #'
  #' The sign of `tp` matters: post-period d{k} is a FORWARD change from
  #' the baseline, so its level lives at t+k ({dv}_f{k}); pre-period dm{k}
  #' is a BACKWARD change (dm{k} = v(t-k) - v(t-1)), so its level lives at
  #' t-k ({dv}_l{k}). Applying b_ref = "t+k" to a placebo row therefore
  #' means "the level the placebo coefficient is measured against", which
  #' is the only reading that makes the euro figure comparable across the
  #' pre and post periods.
  #'
  #' @param dv    outcome variable name (panel column stem)
  #' @param tp    relative year (non-zero integer)
  #' @param b_ref "t+k", "t-1", "t", or a numeric level
  #' @return column name (character), NA when b_ref is numeric
  if (is.numeric(b_ref)) return(NA_character_)
  switch(b_ref,
         "t+k" = if (tp > 0) sprintf("%s_f%d", dv, as.integer(tp)) else
                              sprintf("%s_l%d", dv, abs(as.integer(tp))),
         "t-1" = sprintf("%s_l1", dv),
         "t"   = dv)
}


gd_auc_crvar <- function(vec_x, vec_w, vec_cl) {
  #' Cluster-robust (CR0) variance of a weighted mean sum(w*x)/sum(w).
  #' Used for the reference level at t+k, whose sampling variance is large
  #' (turnover is extremely heavy tailed) and cannot be ignored in the
  #' euro SE — see the file header.
  #'
  #' @param vec_x  numeric vector
  #' @param vec_w  weights (same length)
  #' @param vec_cl cluster id (same length)
  #' @return scalar variance
  mean_w <- sum(vec_w * vec_x) / sum(vec_w)
  sum(rowsum(vec_w * (vec_x - mean_w), vec_cl, reorder = FALSE)^2) / sum(vec_w)^2
}


gd_auc_ref_cols <- function(dv, tp, b_ref) {
  #' Name the two panel columns the conversion needs for one (dv, tp):
  ##  lvl_col  the reference level (log1p/asinh scale, back-transformed later)
  ##  d_col    the d{k}/dm{k} outcome the DiD was fitted on
  #' @param dv    outcome variable name
  #' @param tp    relative year (non-zero integer)
  #' @param b_ref reference level choice (see gd_auc_lvl_col)
  #' @return character vector (lvl_col, d_col); lvl_col is NA for a
  #'   numeric b_ref
  c(lvl_col = gd_auc_lvl_col(dv, tp, b_ref),
    d_col   = sprintf("%s_%s%d", dv, if (tp > 0) "d" else "dm", abs(as.integer(tp))))
}


gd_auc_lvl_tr <- function(dv, vec) {
  #' Back-transform a reference-level column to euros, inferred from the
  #' dv's suffix: an "_log" dv is stored as log1p(<raw>) and an "_ihs" dv
  #' as asinh(<raw>); anything else is already on the euro scale. Getting
  #' this wrong is silent and catastrophic (the level ends up ~1 instead
  #' of ~10^5), so it is a named function rather than an inline ifelse.
  #'
  #' @param dv  outcome variable name (panel column stem)
  #' @param vec numeric vector (the reference-level column)
  #' @return numeric vector on the euro scale
  if (endsWith(dv, "_log")) return(expm1(vec))
  if (endsWith(dv, "_ihs")) return(sinh(vec))
  vec
}


gd_auc_ref_one <- function(dv, tp, c_eu) {
  #' Reference statistics for ONE (dv, tp): the weighted mean level in
  #' euros, its CR0 cluster-robust variance, and the estimation-sample
  #' size. Factored out of the (dv, tp) grid so each combination is
  #' independently reproducible and nothing large is materialised
  #' (coding_style.md §4.11).
  #'
  #' @param dv   outcome variable name
  #' @param tp   relative year (non-zero integer; negative = placebo dm)
  #' @param c_eu config from gc_auc_euro_cfg
  #' @return one-row data.table (dv, tp, ref_level, var_level, n_used)
  dt_panel <- c_eu$dt
  l_col    <- gd_auc_ref_cols(dv, tp, c_eu$b_ref)

  ## the euro reference level: the OBSERVED level (the residual is kept,
  ## see the file header), taken from the column b_ref points at
  lvl <- if (is.na(l_col[["lvl_col"]])) {
    stop("gd_auc_ref_one: a numeric b_ref needs an explicit level column")
  } else {
    dt_panel[[l_col[["lvl_col"]]]]
  }

  ## estimation sample of the corresponding DiD model: the reference
  ## level, the weight and the d{k}/dm{k} outcome all observed (the level
  ## is a back-transform of the log column, so masking on it also masks
  ## on the log value)
  dt_use <- data.table(lvl   = as.numeric(gd_auc_lvl_tr(dv, lvl)),
                       d_out = as.numeric(dt_panel[[l_col[["d_col"]]]]),
                       w     = as.numeric(dt_panel[[c_eu$wcol]]),
                       cl    = dt_panel[[c_eu$uid]])
  dt_use <- dt_use[!is.na(lvl) & !is.na(d_out) & !is.na(w)]
  if (nrow(dt_use) == 0) {
    stop(sprintf("gd_auc_ref_one: empty estimation sample for dv = %s, tp = %d", dv, tp))
  }

  data.table(dv        = dv,
             tp        = as.integer(tp),
             ref_level = sum(dt_use$w * dt_use$lvl) / sum(dt_use$w),
             var_level = gd_auc_crvar(dt_use$lvl, dt_use$w, dt_use$cl),
             n_used    = nrow(dt_use))
}


gd_auc_euro_calc <- function(dt_coef_sel, c_eu) {
  #' Convert log1p coefficients to euros, vectorised over the coefficient
  #' table: build the small per-(dv, tp) reference table, then compute
  #' every quantity in one data.table expression block.
  #'
  #' @param dt_coef_sel selection from md_auc_sel
  #' @param c_eu        config from gc_auc_euro_cfg
  #' @return long data.table, one row per (dv, tp, term) with columns
  #'   dv, tp, term, estimate, std.error, estimate_exp,
  #'   euro, se_euro, se_euro_b, euro_lo, euro_hi,
  #'   b_ref, ref_level, share_of_level, n_used
  ## one row per (dv, tp) -- a handful of rows, cheap even with a map
  dt_ref <- unique(dt_coef_sel[, .(dv, tp)]) %>%
    .[, map(seq_len(.N), ~gd_auc_ref_one(dv[.x], tp[.x], c_eu)) %>% rbindlist]

  ## cross the coefficients with the reference statistics
  dt_out <- merge(copy(dt_coef_sel), dt_ref, by = c("dv", "tp"), sort = FALSE)

  ## growth factor and the exact euro AME: with log1p,
  ##   expm1(log1p(y) + b) - y = (1 + y) * (exp(b) - 1)
  dt_out[, estimate_exp := exp(estimate)]
  dt_out[, euro := (estimate_exp - 1) * (1 + ref_level)]

  ## delta method on (b, ref_level); the second term dominates here, see
  ## the file header. se_euro_b is the b-only version, kept for contrast.
  dt_out[, se_euro_b := estimate_exp * (1 + ref_level) * std.error]
  dt_out[, se_euro   := sqrt(se_euro_b^2 + (estimate_exp - 1)^2 * var_level)]

  ## intervals, and how large the effect is relative to the level it is
  ## applied to (a sanity ratio, ~ exp(b) - 1 by construction)
  dt_out[, `:=`(euro_lo        = euro - 1.96 * se_euro,
                euro_hi        = euro + 1.96 * se_euro,
                b_ref          = as.character(c_eu$b_ref),
                share_of_level = fifelse(ref_level != 0, euro / ref_level, NA_real_))]

  dt_out[, var_level := NULL]
  md_auc_stamp(setorder(dt_out, dv, tp, term), "gd_auc_euro_calc")
}


gd_auc_euro <- function(r_sbw, r_did = NULL, dt_coef = NULL, l_dv = NULL,
                        l_tp = NULL, b_ref = "t+k",
                        b_check = TRUE, verbose = TRUE) {
  #' Euro-scale average marginal effect of each treatment on an auction
  #' outcome, with a standard error that includes the sampling variance of
  #' the reference level (see the file header).
  #'
  #' @param r_sbw       result of gl_sbw_proc (uses $dt_results and $config)
  #' @param r_did       result of gl_sbw_did; only used for the nobs check
  #'                    and, when dt_coef is NULL, to build the coefficient
  #'                    table
  #' @param dt_coef     coefficient table from gd_sbw_coef_did(r_did)
  #' @param l_dv        outcomes to convert; NULL = the config's auction
  #'                    outcomes that were actually estimated
  #' @param l_tp        relative years; NULL = every period present (pre and
  #'                    post); e.g. -4:-2 for the placebos, 1:6 for the effects
  #' @param b_ref       reference level: "t+k" (default), "t-1", "t", or a number
  #' @param b_check     cross-check the reconstruction against the model nobs
  #' @param verbose     progress messages
  #' @return data.table, one row per (dv, tp, term): see gd_auc_euro_calc
  c_eu      <- gc_auc_euro_cfg(r_sbw, b_ref)
  dt_coef   <- md_auc_coef_get(r_did, dt_coef)
  dt_sel    <- md_auc_sel(dt_coef, r_sbw$config, l_dv, l_tp)
  dt_euro   <- gd_auc_euro_calc(dt_sel, c_eu)

  if (b_check) dt_euro <- md_auc_check_n(dt_euro, r_did)
  if (verbose) {
    cat(sprintf("gd_auc_euro: %d row(s); reference level = %s; residuals kept\n",
                nrow(dt_euro), as.character(c_eu$b_ref)))
  }
  md_auc_stamp(dt_euro, "gd_auc_euro")
}


md_auc_stamp <- function(dt, str_fn) {
  #' Attach the generating-function attribute (coding_style.md §2) so a
  #' data.table in a session can be traced back to the stage that made it.
  #'
  #' @param dt     data.table
  #' @param str_fn generating function name
  #' @return dt, with the attribute set
  attr(dt, "gnrtdby") <- str_fn
  dt
}


md_auc_check_n <- function(dt_euro, r_did) {
  #' Cross-check the reconstruction sample size against the model's own
  #' nobs and warn on a mismatch (a mismatch means the reference level was
  #' computed on different rows than the DiD was fitted on).
  #'
  #' @param dt_euro output of gd_auc_euro_calc
  #' @param r_did   result of gl_sbw_did, or NULL to skip
  #' @return dt_euro (copy) with n_model added
  dt_out <- copy(dt_euro)
  if (is.null(r_did)) return(dt_out[, n_model := NA_integer_])

  ## fixest_multi: [[ takes an INTEGER index only, so match the lhs name
  l_lhs <- names(r_did)
  ## d{k} for the post period, dm{k} for the pre period -- the lhs names
  ## follow the panel column names, so the sign of tp decides the prefix
  dt_out[, n_model := map_int(
           sprintf("lhs: %s_%s%d", dv, ifelse(tp > 0, "d", "dm"), abs(tp)),
           function(str_nm) {
    i_mdl <- match(str_nm, l_lhs)
    if (is.na(i_mdl)) NA_integer_ else as.integer(nobs(r_did[[i_mdl]]))
  })]

  ## one warning is enough: they all share the same cause
  dt_bad <- dt_out[!is.na(n_model) & n_model != n_used]
  if (nrow(dt_bad) > 0) {
    warning(sprintf(paste0("md_auc_check_n: row count mismatch for %s at tp = %d ",
                           "(reconstruction %d, model nobs %d); the reference ",
                           "level is computed on the reconstruction"),
                    dt_bad$dv[1], dt_bad$tp[1], dt_bad$n_used[1], dt_bad$n_model[1]),
            call. = FALSE)
  }
  md_auc_stamp(dt_out, "md_auc_check_n")
}
