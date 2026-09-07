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

