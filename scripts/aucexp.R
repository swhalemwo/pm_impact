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
  #'                    "t-1" (pre-treatment level) or "t"
  #' @return config list (dt, wcol, uid, b_ref)
  if (is.null(r_sbw$dt_results) || is.null(r_sbw$config)) {
    stop("gc_auc_euro_cfg: r_sbw must be the result of gl_sbw_proc ",
         "(needs $dt_results and $config)")
  }
  if (!is.data.frame(r_sbw$dt_results)) {
    stop("gc_auc_euro_cfg: r_sbw$dt_results is not a data.frame")
  }
  ## match.arg() would accept this but complain unhelpfully on a non-character
  if (!is.character(b_ref) || length(b_ref) != 1L) {
    stop("gc_auc_euro_cfg: b_ref must be one of 't+k', 't-1', 't'")
  }
  b_ref <- match.arg(b_ref, c("t+k", "t-1", "t"))

  ## the weight column and the cluster id drive the whole estimation, so
  ## their absence is a config error, not a data problem -> check now
  c_eu <- list(dt    = r_sbw$dt_results,
               wcol  = r_sbw$config$wcol %||% "w_sbw",
               uid   = r_sbw$config$uid,
               tid   = r_sbw$config$tid,
               b_ref = b_ref)
  if (!c_eu$wcol %in% names(c_eu$dt)) {
    stop(sprintf("gc_auc_euro_cfg: weight column '%s' not in dt_results", c_eu$wcol))
  }
  if (is.null(c_eu$uid) || !c_eu$uid %in% names(c_eu$dt)) {
    stop("gc_auc_euro_cfg: config has no usable uid column (needed for the ",
         "cluster-robust reference variance)")
  }
  ## c_eu$tid is only needed by gd_auc_euro_sum (the stacked fit's year FE)
  if (is.null(c_eu$tid) || !c_eu$tid %in% names(c_eu$dt)) {
    stop("gc_auc_euro_cfg: config has no usable tid column")
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


md_auc_sel <- function(dt_coef, c_sbw_did, l_dv = NULL, l_tp = NULL,
                       str_dv_pattern = NULL) {
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
  #' @param str_dv_pattern optional regex that trims the DEFAULT dv set (when
  #'                  l_dv is NULL) by name; an explicit l_dv is used as given
  #' @return data.table (copy), one row per (dv, tp, term)
  l_cfg        <- unique(c(c_sbw_did$outcome_cols, c_sbw_did$outcome_cols_cur))
  b_default_dv <- is.null(l_dv)

  ## default dv set = the config's outcomes that actually appear in the
  ## coefficient table (gl_sbw_did may have been called on a subset)
  if (b_default_dv) {
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

  ## Trim the DEFAULT dv set with the pattern: the config also carries the
  ## exhibition counts (exn_*_cum2) and the raw-scale auction columns
  ## (auc_turnover, aucturn_cum2), none of which have a euro meaning. An
  ## explicitly supplied l_dv is taken as given -- the caller has decided.
  if (!is.null(str_dv_pattern) && b_default_dv) {
    l_dv <- l_dv[grepl(str_dv_pattern, l_dv)]
    if (length(l_dv) == 0) {
      stop(sprintf("md_auc_sel: no dv matches the pattern '%s' (available: %s)",
                   str_dv_pattern, paste(unique(dt_coef$dv), collapse = ", ")))
    }
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
  #' @param b_ref "t+k", "t-1" or "t"
  #' @return column name (character)
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
  #' @return character vector (lvl_col, d_col)
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
  lvl <- dt_panel[[l_col[["lvl_col"]]]]

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


gd_auc_ref_grid <- function(dt_coef_sel, c_eu) {
  #' Reference moments for every (dv, tp) in the selection, as one long
  #' data.table -- computed per combination by gd_auc_ref_one, never by an
  #' inline loop over the whole panel.
  #'
  #' @param dt_coef_sel selection from md_auc_sel
  #' @param c_eu        config from gc_auc_euro_cfg
  #' @return data.table (dv, tp, ref_level, var_level, n_used, m_ihs)
  unique(dt_coef_sel[, .(dv, tp)]) %>%
    .[, map(seq_len(.N), ~gd_auc_ref_one(dv[.x], tp[.x], c_eu)) %>% rbindlist] %>%
    setorder(dv, tp) %>%
    md_auc_stamp("gd_auc_ref_grid")
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
  dt_ref <- gd_auc_ref_grid(dt_coef_sel, c_eu)

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
  ## the euro conversion assumes a log1p link, so by default only the log
  ## auction outcomes are converted; pass l_dv to override
  dt_sel    <- md_auc_sel(dt_coef, r_sbw$config, l_dv, l_tp, str_dv_pattern = "^auc.*log$")
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


md_auc_stack_one <- function(str_dv, l_tp, c_eu, l_treat) {
  #' Stack the post-period blocks of ONE dv into a long table and fit them as
  #' a single regression, so that vcov() returns the JOINT covariance of the
  #' horizon-specific coefficients.
  #'
  #' Every regressor is interacted with the horizon, so the stacked fit
  #' decomposes exactly into the separate per-horizon fits (verified to
  #' machine precision); the point of stacking is only to obtain the
  #' off-diagonal covariance blocks, which the separate fits cannot give.
  #' The six horizons share the same artists, weights and overlapping outcome
  #' windows, so their estimates are strongly positively correlated
  #' (measured ~0.6 on the real panel) and a sum of them has a variance well
  #' above the sum of the variances.
  #'
  #' @param str_dv  outcome variable name
  #' @param l_tp    horizons to stack, in accumulation order (positive or
  #'                negative; the sign only drives which columns are read)
  #' @param c_eu    config from gc_auc_euro_cfg
  #' @param l_treat treatment column names (the RHS of the DiD)
  #' @return list(fit = the stacked fixest model, l_tp = horizons used)

  dt_panel <- c_eu$dt

  ## the factor level is the POSITION in l_tp, not tp itself: a pre-period
  ## block would otherwise create factor levels like "-2", which is fragile
  ## in formula terms and in coefficient names
  l_long <- map(seq_along(l_tp), function(i_pos) {
    l_col <- gd_auc_ref_cols(str_dv, l_tp[i_pos], c_eu$b_ref)
    ## keep the panel's OWN column names for uid / tid: the RHS and the
    ## cluster argument are built from them, so renaming would break feols
    d <- data.table(k = i_pos, dv_val = dt_panel[[l_col[["d_col"]]]])
    ## all four original column names must survive: feols resolves the
    ## weights, the cluster and the RHS against this table by name
    d[[c_eu$wcol]] <- dt_panel[[c_eu$wcol]]
    d[[c_eu$tid]]  <- dt_panel[[c_eu$tid]]
    d[[c_eu$uid]]  <- dt_panel[[c_eu$uid]]
    ## the treatment columns travel on the row: the stacked RHS needs them
    for (str_t in l_treat) d[[str_t]] <- dt_panel[[str_t]]
    ## drop rows without an outcome or a weight (the estimation sample)
    d[!is.na(dv_val) & !is.na(get(c_eu$wcol))]
  })
  dt_long <- rbindlist(l_long)

  ## factor(k):<term> for every regressor -- identical to fitting each
  ## horizon separately with 0 + factor(tid) + treatments. With a single
  ## horizon there is nothing to interact with (factor(k) would have one
  ## level and contrasts() would fail), so use the plain formula.
  str_rhs <- if (length(l_tp) > 1L) {
    sprintf("0 + factor(k):factor(%s) + %s", c_eu$tid,
            paste(sprintf("factor(k):%s", l_treat), collapse = " + "))
  } else {
    sprintf("0 + factor(%s) + %s", c_eu$tid, paste(l_treat, collapse = " + "))
  }

  fit <- feols(as.formula(paste("dv_val ~", str_rhs)), dt_long,
               weights = as.formula(paste0("~ ", c_eu$wcol)),
               cluster = as.formula(paste0("~ ", c_eu$uid)))
  list(fit = fit, l_tp = l_tp)
}


md_auc_check_betas <- function(l_stack, dt_blk, l_tp, l_treat,
                               d_tol = 1e-8) {
  #' Assert that the stacked fit reproduces the per-horizon betas it is being
  #' used to attach a covariance to.
  #'
  #' gd_auc_cum_one takes the point estimates from the pipeline's separate
  #' fits (dt_eu) and the covariance from the stacked fit. That is only
  #' coherent while the two agree -- which they do by construction, since the
  #' stacked formula interacts every regressor with the block and therefore
  #' decomposes into the per-block fits. A silent divergence would mix two
  #' models with no visible sign, so it is checked rather than trusted.
  #'
  #' @param l_stack list from md_auc_stack_one
  #' @param dt_blk  per-horizon rows for this dv and block, accumulation order
  #' @param l_tp    horizons in accumulation order (positions in the stack)
  #' @param l_treat treatment names
  #' @param d_tol   tolerance on the absolute difference
  #' @return TRUE, invisibly; errors on a mismatch
  b_multi <- length(l_tp) > 1L
  d_worst <- 0
  for (str_t in l_treat) {
    vec_st <- coef(l_stack$fit)[md_auc_stack_names(seq_along(l_tp), str_t, b_multi)]
    d_blk  <- dt_blk[term == str_t]
    vec_se <- d_blk$estimate[match(l_tp, d_blk$tp)]
    if (length(vec_st) != length(vec_se) || anyNA(vec_se)) next
    d_worst <- max(d_worst, max(abs(vec_st - vec_se)))
  }
  if (d_worst > d_tol) {
    stop(sprintf(paste0("md_auc_check_betas: stacked and separate betas differ ",
                        "by %.3g (tolerance %.1g); the cumulative point estimate ",
                        "would come from a different fit than its SE"),
                 d_worst, d_tol))
  }
  invisible(TRUE)
}


md_auc_stack_names <- function(l_pos, str_term, b_multi) {
  #' Coefficient names of one treatment in the stacked fit.
  #'
  #' With several blocks the regressors are =factor(k)<pos>:<term>=, where
  #' <pos> is the position of the block in the stack (md_auc_stack_one uses
  #' positions, not tp, so pre-period blocks with negative tp work too). With
  #' a single block there is nothing to interact and the plain term name is
  #' used. Getting this wrong is silent (the vcov subset comes back empty),
  #' so it lives in its own function.
  #'
  #' @param l_pos    positions in the stacked fit (1..n)
  #' @param str_term treatment name
  #' @param b_multi  TRUE when the stacked fit has more than one block
  #' @return character vector of coefficient names
  if (b_multi) sprintf("factor(k)%d:%s", as.integer(l_pos), str_term) else str_term
}


gd_auc_cum_one <- function(dt_term, m_vcov, b_multi, str_dv, str_term, b_ref) {
  #' Cumulative euro effects for ONE (dv, term): one row per endpoint, i.e.
  #' the sum over horizons 1..m for every m in the available horizons.
  #'
  #' The sum is S_m = sum_{k<=m} (exp(b_k) - 1)(1 + m_k), so
  #'     dS_m/db_k = exp(b_k)(1 + m_k)   for k <= m, 0 otherwise
  #' and
  #'     Var(S_m) = grad' V_m grad + sum_{k<=m} (exp(b_k) - 1)^2 Var(m_k)
  #' with V_m the matching block of the joint (stacked) covariance. The
  #' second term is the reference-level variance, so with a single horizon
  #' this reproduces gd_auc_euro's se_euro.
  #'
  #' @param dt_term data.table for this (dv, term), ordered by tp, with
  #'                columns tp, estimate, std.error, euro, ref_level,
  #'                var_level
  #' @param m_vcov  joint cluster-robust vcov of the stacked fit
  #' @param b_multi TRUE when the stacked fit spans several horizons
  #' @param str_dv  outcome name (bookkeeping)
  #' @param str_term treatment name (bookkeeping)
  #' @param b_ref   reference level choice (bookkeeping)
  #' @return data.table, one row per cumulative endpoint
  map(seq_len(nrow(dt_term)), function(i_end) {
    dt_sub <- dt_term[seq_len(i_end)]
    l_nm <- md_auc_stack_names(seq_len(i_end), str_term, b_multi)
    l_miss <- setdiff(l_nm, rownames(m_vcov))
    if (length(l_miss) > 0) {
      stop(sprintf("gd_auc_cum_one: coefficient(s) not in the stacked fit: %s",
                   paste(l_miss, collapse = ", ")))
    }
    ## gradient of the sum wrt each horizon's beta, at the fitted values
    vec_grad <- exp(dt_sub$estimate) * (1 + dt_sub$ref_level)
    m_sub    <- m_vcov[l_nm, l_nm, drop = FALSE]
    ## reference-level variance, one term per horizon in the sum
    vec_ref  <- (exp(dt_sub$estimate) - 1)^2 * dt_sub$var_level
    f_var    <- as.numeric(t(vec_grad) %*% m_sub %*% vec_grad) + sum(vec_ref)
    f_sum    <- sum(dt_sub$euro)
    data.table(dv = str_dv, term = str_term,
               period = if (dt_sub$tp[i_end] > 0) "post" else "pre",
               tp_end = dt_sub$tp[i_end],
               l_tp = paste(dt_sub$tp, collapse = ","),
               n_tp = nrow(dt_sub),
               sum_euro = f_sum,
               se_euro = sqrt(f_var),
               euro_lo = f_sum - 1.96 * sqrt(f_var),
               euro_hi = f_sum + 1.96 * sqrt(f_var),
               p_euro = 2 * (1 - pnorm(abs(f_sum / sqrt(f_var)))),
               sum_beta = sum(dt_sub$estimate),
               ## what an independence assumption would have given
               se_naive = sqrt(sum(vec_grad^2 * dt_sub$std.error^2) + sum(vec_ref)),
               b_ref = as.character(b_ref))
  }) %>% rbindlist
}


gd_auc_sum_one_term <- function(str_term, str_dv, dt_eu, m_vcov, b_multi,
                                l_tp_order) {
  #' Cumulative euro effects for one treatment of one outcome.
  #'
  #' @param str_term  treatment name
  #' @param str_dv    outcome name
  #' @param dt_eu     per-horizon euro table (gd_auc_euro_calc + var_level)
  #' @param m_vcov    joint vcov of the stacked fit for this outcome
  #' @param b_multi   TRUE when the stacked fit spans several horizons
  #' @param l_tp_order the block's horizons in accumulation order; the rows
  #'                  are reordered to match, because the stacked fit's
  #'                  coefficient positions follow that order
  #' @return data.table from gd_auc_cum_one, or NULL when the treatment is
  #'   absent from this outcome's model
  dt_term <- dt_eu[dv == str_dv & term == str_term]
  dt_term <- dt_term[match(l_tp_order, tp)]
  if (nrow(dt_term) == 0) return(NULL)
  gd_auc_cum_one(dt_term, m_vcov, b_multi, str_dv, str_term,
                 dt_term$b_ref[1])
}


gd_auc_sum_one_dv <- function(str_dv, dt_eu, c_eu, l_treat,
                              str_period = c("post", "pre")) {
  #' Cumulative euro effects for one outcome and one period block: stack the
  #' block's horizons, take the joint covariance once, then run every
  #' treatment through it.
  #'
  #' @param str_dv    outcome name
  #' @param dt_eu     per-horizon euro table (gd_auc_euro_calc + var_level)
  #' @param c_eu      config from gc_auc_euro_cfg
  #' @param l_treat   treatment names
  #' @param str_period "post" accumulates tp = 1..m; "pre" accumulates
  #'                  tp = -2, -3, ... (i.e. backwards in time, the mirror
  #'                  image of the post series)
  #' @return data.table, one row per (term, cumulative endpoint)

  str_period <- match.arg(str_period)
  dt_dv <- dt_eu[dv == str_dv]
  ## accumulation order: forward for the post block, backward for the pre
  ## block (so that each row adds one more year of history)
  l_tp_dv <- if (str_period == "post") sort(unique(dt_dv[tp > 0, tp])) else
                                      sort(unique(dt_dv[tp < 0, tp]), decreasing = TRUE)
  if (length(l_tp_dv) == 0) return(NULL)
  l_stack <- md_auc_stack_one(str_dv, l_tp_dv, c_eu, l_treat)

  ## dt_coef_did[dv == "auc_turnover2_log" & term == "exn_toti" & tp > 0]
  ## tidy(l_stack$fit) %>% adt %>% .[grepl("exn_toti", term)]
  ## tidy(l_stack$fit) %>% adt %>% print(n=800)
  ## looks good
  
  m_vcov  <- vcov(l_stack$fit)
  b_multi <- length(l_tp_dv) > 1L
  ## The POINT ESTIMATES come from dt_eu (the pipeline's separate fits) while
  ## the SE comes from this stacked fit. That is only legitimate because the
  ## stacked betas are numerically identical to the separate ones -- the
  ## formula interacts every regressor with the block, so it decomposes into
  ## the per-block fits. Check it rather than trust it: if the stacked fit
  ## ever changes (different sample handling, a missing column), the point
  ## estimate and its SE would silently come from two different models.

  md_auc_check_betas(l_stack, dt_dv[tp %in% l_tp_dv], l_tp_dv, l_treat)

  map(l_treat, gd_auc_sum_one_term, str_dv = str_dv,
      dt_eu = dt_dv[tp %in% l_tp_dv], m_vcov = m_vcov, b_multi = b_multi,
      l_tp_order = l_tp_dv) %>%
    rbindlist
  
}


gd_auc_euro_sum <- function(r_sbw, r_did = NULL, dt_coef = NULL, l_dv = NULL,
                            l_tp = NULL, b_ref = "t+k",
                            str_dv_pattern = "^auc.*log$", b_pre = FALSE,
                            b_verbose = TRUE) {
  #' Cumulative euro effects: for every horizon m, the sum of the euro
  #' effects over horizons 1..m, with a delta-method SE that respects the
  #' correlation between horizons.
  #'
  #' Use this (rather than the cumulative-outcome columns) whenever the total
  #' is to be shown alongside the annual effects: euros are additive, so the
  #' total is the sum of the annual amounts by construction, and the reader
  #' never has to reconcile two different log scales.
  #'
  #' The six per-horizon estimates come from the same artists with
  #' overlapping outcome windows, so they correlate at roughly 0.6; the joint
  #' covariance comes from a single stacked regression whose point estimates
  #' are identical to the separate fits (see md_auc_stack_one). Summing their
  #' variances in quadrature would understate the SE by a factor of about two.
  #'
  #' @param r_sbw   result of gl_sbw_proc
  #' @param r_did   result of gl_sbw_did (coef table, optional)
  #' @param dt_coef coefficient table from gd_sbw_coef_did
  #' @param l_dv    outcomes to sum; NULL = the config's auction outcomes
  #' @param l_tp    horizons; NULL = every post period present
  #' @param b_ref   reference level: "t+k" (default), "t-1" or "t"
  #' @param str_dv_pattern regex on the dv names; the default keeps the log
  #'                auction outcomes only -- the euro conversion assumes a
  #'                log1p link, so raw-scale columns (auc_turnover,
  #'                aucturn_cum2) and the exhibition counts (exn_*_cum2)
  #'                are excluded
  #' @param b_pre   also accumulate the PRE-period placebos (tp = -2, -3, ...),
  #'                giving the mirror image of the post series for the left
  #'                half of an event-study plot. The pre and post blocks are
  #'                stacked separately, so each series has its own joint
  #'                covariance; they are not summed across the exhibition.
  #' @param b_verbose progress message
  #' @return data.table, one row per (dv, term, cumulative endpoint) with
  #'   columns dv, term, tp_end, l_tp, n_tp, sum_euro, se_euro, euro_lo,
  #'   euro_hi, p_euro, sum_beta, se_naive, b_ref

  c_eu    <- gc_auc_euro_cfg(r_sbw, b_ref)
  dt_coef <- md_auc_coef_get(r_did, dt_coef)
  dt_sel  <- md_auc_sel(dt_coef, r_sbw$config, l_dv, l_tp, str_dv_pattern)

  ## the post block is always needed; the pre block is optional
  if (nrow(dt_sel[tp > 0]) == 0) {
    stop("gd_auc_euro_sum: no post-period coefficients to sum (tp > 0)")
  }
  if (!b_pre) dt_sel <- dt_sel[tp > 0]

  ## per-horizon euro amounts; the point estimates of the sums must equal
  ## the sums of these, so they are reused rather than recomputed
  dt_eu  <- gd_auc_euro_calc(dt_sel, c_eu)
  ## var_level is the second variance term in Var(S), so it is merged in
  dt_eu  <- merge(dt_eu, gd_auc_ref_grid(dt_sel, c_eu)[, .(dv, tp, var_level)],
                  by = c("dv", "tp"), sort = FALSE)

  ## dt_eu[dv == "auc_turnover2_log" & term == "exn_toti" & tp == 6, .(euro, se_euro)]
  ## dt_euro[dv == "auc_turnover2_log" & term == "exn_toti" & tp == 6, .(euro, se_euro)]

  ## post and (optionally) pre are separate series: each is stacked on its
  ## own, because the two blocks answer different questions and their
  ## cross-block covariance is not needed for either
  dt_out <- map(unique(dt_eu$dv), gd_auc_sum_one_dv, dt_eu = dt_eu,
                c_eu = c_eu, l_treat = r_sbw$config$treatment_cols,
                str_period = "post") %>%
    rbindlist

  ## dt_out[dv == "auc_turnover2_log" & term == "exn_toti" & n_tp == 6][, .(sum_euro, se_euro)]
  
  gd_auc_sum_one_dv("auc_turnover2_log", dt_eu, c_eu, r_sbw$config$treatment_cols, str_period = "post")

  if (b_pre && nrow(dt_eu[tp < 0]) == 0) {
    warning("gd_auc_euro_sum: b_pre = TRUE but the coefficient table has no ",
            "pre-period rows (tp < 0); only the post series is returned",
            call. = FALSE)
  }
  if (b_pre && nrow(dt_eu[tp < 0]) > 0) {
    dt_out <- rbind(dt_out,
                    map(unique(dt_eu$dv), gd_auc_sum_one_dv, dt_eu = dt_eu,
                        c_eu = c_eu, l_treat = r_sbw$config$treatment_cols,
                        str_period = "pre") %>% rbindlist, fill = TRUE)
  }

  if (b_verbose) {
    cat(sprintf("gd_auc_euro_sum: %d row(s) over %d outcome(s); max horizons %d\n",
                nrow(dt_out), uniqueN(dt_out$dv), max(dt_out$n_tp)))
    cat(sprintf("  SE inflation over the independence assumption: %.2f\n",
                mean(dt_out$se_euro / dt_out$se_naive)))
  }
  md_auc_stamp(setorder(dt_out, dv, term, tp_end), "gd_auc_euro_sum")
}
