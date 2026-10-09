# MWE: stable balancing weights (optweightMV) + difference-in-differences for
# the effect of exhibitions by institution type on artists' future exhibitions.
#
# Runtime ~1 minute.
#
# greifer.csv = 3000 randomly sampled artists x the treatment years 2010-2012.
#   PeopleID, year
#   exn_pi, exn_oes, exn_pm, exn_toti        exhibitions in that year by type
#                                            (the 4 treatments: public
#                                            institutions, other exhibition
#                                            spaces, private museums, top-tier)
#   exn_all_cum                              exhibitions strictly before that
#                                            year (sample restriction: >= 10)
#   exn_X_l1 .. exn_X_l4                     same, lagged 1-4 years
#                                            (the balancing variables)
#   exn_X_d1 .. exn_X_d6                     exn_X(t+k) - exn_X(t-1)
#                                            (outcomes, k = 1..6)
#   exn_X_dm2 .. exn_X_dm4                   exn_X(t-k) - exn_X(t-1)
#                                            (pre-period placebo outcomes)
#   Female, career_year                      time-invariant / career covariates
#   exn_all_cum_ihs, exn_gallery_cum_ihs,
#   exn_event_cum_ihs, auc_turnover_cum_ihs  covariates, inverse hyperbolic
#                                            sine transformed, cumulative
#                                            (strictly before t)
#
# All lag/lead/change columns are built from the artist-year panel with
# shift(exn_X, k, type = "lag"/"lead") by PeopleID, i.e. relative to t.

library(data.table)
library(optweight)
library(fixest)
library(parallel)

data_file <- "~/Dropbox/phd/papers/foster/misc/greifer.csv"

dt <- fread(data_file)
dt <- dt[exn_all_cum >= 10]  # established artists only (>= 10 prior exhibitions)

l_treat <- c("exn_pi", "exn_oes", "exn_pm", "exn_toti")
l_lags  <- paste0(rep(l_treat, each = 4), "_l", 1:4)  # 16 lagged treatments
l_cov   <- c("Female", "career_year", "exn_all_cum_ihs", "exn_gallery_cum_ihs",
             "exn_event_cum_ihs", "auc_turnover_cum_ihs")
l_yrs   <- sort(unique(dt$year))

## 1. SBW: each treatment is balanced on the lags of all four treatments plus
##    the covariates. Balance constraints are exact (tols = 0, the default);
##    fitted separately per calendar year, ATE estimand, so each year's targets
##    are that year's covariate means.
##    eps_rel/eps_abs are loosened to 1e-2 only so that this small sample fits
##    in ~20s per year (the paper uses 1e-3, i.e. ~2-3 min per year here).
l_fm  <- lapply(l_treat, function(tv) reformulate(c(l_lags, l_cov), tv))

l_sbw <- mclapply(l_yrs, \(y) optweightMV(l_fm, data = dt[year == y], estimand = "ATE",
                                          eps_rel = 1e-2, eps_abs = 1e-2),
    mc.cores = 3)

dt[, w_sbw := NA_real_]
for (i in seq_along(l_yrs)) dt[year == l_yrs[i], w_sbw := l_sbw[[i]]$weights]
dt[, .(n = .N, ess = round(sum(w_sbw)^2 / sum(w_sbw^2)), max_w = round(max(w_sbw), 1)),
   by = year]

## 2. balance after weighting (treatment vs. lags and covariates), year 2010
cobalt::bal.tab(l_sbw[[1]])

## 3. DiD: change in exhibitions from t-1 to t+6, jointly regressed on the four
##    doses at t, year fixed effects, artist-clustered SEs, SBW-weighted
l_m <- lapply(l_treat, function(tv) {
    feols(reformulate(c("0 + factor(year)", l_treat), paste0(tv, "_d6")),
          data = dt, weights = ~w_sbw, cluster = ~PeopleID)
})
names(l_m) <- paste0(l_treat, "_d6")
etable(l_m, keep = "^exn", digits = 3)

## 4. same on the pre-period change dm2 (= t-2 minus t-1): should be ~0
l_m_pre <- lapply(l_treat, function(tv) {
    feols(reformulate(c("0 + factor(year)", l_treat), paste0(tv, "_dm2")),
          data = dt, weights = ~w_sbw, cluster = ~PeopleID)
})
names(l_m_pre) <- paste0(l_treat, "_dm2")
etable(l_m_pre, keep = "^exn", digits = 3)
