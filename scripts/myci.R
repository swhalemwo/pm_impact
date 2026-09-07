gf_did2 <- function(l_dvs, l_tps, l_iv, l_fe = "stratum") {
  #' can use weights
  #' estimate DVs of multiple forms in one go -> allows still single pass

  l_tps2 <- l_tps[l_tps != 0]

  dt_combs <- expand.grid(dv = l_dvs, tp = l_tps2) %>% adt %>% 
    .[, dv_fmt := ifelse(tp > 0, sprintf("%s_d%s", dv, tp), sprintf("%s_dm%s", dv, abs(tp)))]

  sprintf("c(%s) ~ %s | %s",
          paste0(dt_combs[, dv_fmt], collapse = ","), 
          paste0(l_iv, collapse = " + "),
          l_fe) %>% as.formula

}


gd_reg_pama <- function(r_pama) {

  dt_coef <- r_pama %>% summary %>% adt(keep.rownames = "tp") %>%
    .[, .(tp, estimate, std.error)] %>% 
    .[, statistic := estimate/std.error] %>% .[, p.value := 2* (1-pnorm(abs(statistic)))]

  dt_gof <- r_pama %$% data.table(lag = lag, qoi = qoi, pooled = pooled, se.method = se.method)

  list(dt_coef = dt_coef, dt_gof = dt_gof)
  

}



gd_reg_fixest_multi <- function(r_fixest_multi) {

  dt_coef <- map(r_fixest_multi, tidy) %>% rbindlist(idcol = "mdl")

  dt_gof <- map(r_fixest_multi, glance) %>% rbindlist(idcol = "mdl")

  list(dt_coef = dt_coef, dt_gof = dt_gof)

}

gd_ess <- function(r_weight) {
  dt_ess <- data.table(arm = as.character(r_weight$treat), w = as.numeric(r_weight$weights)) %>%
    .[, .(n_obs = .N, ess = sum(w)^2 / sum(w^2)), by = arm]
  return(dt_ess)
}


check_cors_owmv <- function(r_ow) {
  #' inspect correlations after weighting

  l_treats <- map(r_ow$treat.list, ~attr(.x, "treat.name")) %>% unlist
  ## assume same covs for all treats

  if (isFALSE(Reduce(\(x,y) if(setequal(x,y)) x else return(FALSE), map(r_ow$covs.list, names)))) {
    stop("requires same covariates for each treatment")
  }
  
  dt_covs <- r_ow$covs.list[1] %>% adt

  dt_treats <- r_ow$treat.list %>% adt %>% setnames(l_treats)
  dt_cbn <- cbind(dt_treats, dt_covs)
  
  ## correlation of weighted
  dt_cor_ww <- suppressWarnings(wtd.cor(dt_treats, num_vars(dt_cbn), weight = r_ow$weights)) %>%
    map(~adt(.x, keep.rownames = "vrbl")) %>% ## rbindlist
      ## dt_cor_ww %>%
    rbindlist(idcol = "mes") %>%
    melt(id.vars = c("vrbl", "mes"), measure.vars = l_treats) %>% 
    dcast(vrbl + variable ~ mes) %>% .[, src := "mv"] %>%
    .[, correlation := round(correlation, 3)]

  ## correlation unweighted
  dt_cor_raw <- suppressWarnings(wtd.cor(dt_treats, num_vars(dt_cbn))) %>%
    map(~adt(.x, keep.rownames = "vrbl")) %>% rbindlist(idcol = "mes") %>%
    melt(id.vars = c("vrbl", "mes"), measure.vars = l_treats) %>%
    dcast(vrbl + variable ~ mes) %>% .[, src := "raw"]
  

  rbind(dt_cor_ww, dt_cor_raw) %>%
    .[vrbl != variable] %>%
    ## .[!(vrbl %in% l_treats & variable %in% l_treats)] %>% 
    dcast(vrbl + variable ~ src, value.var = c("correlation", "p.value")) %>%
    .[, corprop_mv_raw := correlation_mv/correlation_raw] # see to which proportion correlation is shrunk


}

check_cors_weightit <- function(rx) {

  dt_treat <- data.table(treat = rx$treat)
  dt_covs <- adt(rx$covs)

  ## l_covs <- names(dt_covs)
  dt_covs_fac <- fact_vars(dt_covs)

  ## expand factor columns
  dt_covs_expnd <- apply(dt_covs_fac, 2, \(x) model.matrix(~x-1)) %>%
    imap(~adt(.x) %>% setnames(old = names(.), new = paste0(.y, names(.)))) %>% cbindlist
  
  dt_covs_cbn <- cbind(num_vars(dt_covs), dt_covs_expnd)

  r_cor_wtd <- wtd.cor(dt_treat, dt_covs_cbn, weight = rx$weight)
  r_cor_uwtd <- wtd.cor(dt_treat, dt_covs_cbn)

  imap(list("wtd" = r_cor_wtd, "uwtd" = r_cor_uwtd),
       ~adt(.x, keep.rownames = "vrbl")[, src := .y]) %>% rbindlist %>%
    dcast(vrbl ~ src, value.var = "correlation") %>% adt %>% 
    .[, corprop_wtd := round(wtd/uwtd,4)]

}

#' extract coefs from DiD model estimated with fixest
#' @param rx DiD model estimated with fixest (multiple LHS)
#' @return data.table with coef estimates, timepoints (tp), terms
gd_coef_did <- function(rx) {

  dt_coef <- map(rx, tidy) %>% rbindlist(idcol = "mdl")

  dt_coef[, tp_prep := gsub("lhs: exn_", "", mdl)] %>%
    .[, c("dv", "tp_prep2") := tstrsplit(tp_prep, "_")] %>% 
    .[, `:=`(tdir = fifelse(grepl("dm", tp_prep2), -1, 1), tlen = as.integer(gsub("[a-z]", "", tp_prep2)))] %>%
    .[, tp := tdir * tlen] %>%
    .[, `:=`(tdir = NULL, tlen = NULL)]

  return(dt_coef)
  
}

#' visualize dt_coef_did
#' @param dt_coef_did data.table of did coefs produced by gd_coef_did
gp_coef_did <- function(dt_coef_did) {
  
  dt_coef_did %>%
    .[grepl("exn_", term)] %>% 
    ggplot(aes(x=tp, y = estimate, color = term,
                    ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error)) +
    geom_point() + geom_errorbar()
}

