

gp_evst <- function(dt_coef) {
    #' event study plot
    dt_coef[dv %in% l_treats & term != "no_treatment"] %>% gp_coef_did +
        theme(legend.position = "bottom")
}

gp_evst_wintcpt <- function(dt_coef) {
    #' event study with intercept
    
    dt_coef[dv %in% l_treats] %>% gp_coef_did
}


gp_evst_cum <- function(dt_coef) {
    #' event study cumulative

    dt_coef[grepl("cum2", dv) & grepl(paste0(l_treats, collapse = "|"), dv)] %>%
        .[term != "no_treatment"] %>% 
        gp_coef_did
}


gp_evst_cum_wintcpt <- function(dt_coef) {
    #' event study cumulative with intercept
    dt_coef[grepl("cum2", dv) & grepl(paste0(l_treats, collapse = "|"), dv)] %>%
        gp_coef_did
}

gp_compact <- function(dt_coef) {
    dt_coef[term != "no_treatment" & grepl(paste0(l_treats, collapse = "|"), dv)] %>%
        gp_coef_did_compact
}

gp_cpct_intr_cpeqntl <- function(dt_results) {
  dtx <- dt_results[!is.na(w_sbw)] %>%
    .[, cpe_qntl := as.integer(cut(exn_all_cum,
                                   breaks = quantile(exn_all_cum, probs = c(0, 0.5, 0.9, 1)),
                                   include.lowest = T, ordered_result = T)), year] %>%
    .[, cpe_qntl := factor(cpe_qntl)]
  



  f_intr <- sprintf("c(%s) ~ cpe_qntl * (%s) | year",
                    paste0(l_treats, "_cum2_d6", collapse = ", "),
                    paste0(l_treats, collapse = " + ")) %>% as.formula

  r_intr <- feols(f_intr, dtx,  weights = ~w_sbw, cluster = ~PeopleID)



  dt_coefs <- lapply(r_intr, \(x) map(l_treats,
                                      ~gd_slopes_nway(x, .x, l_grp_vars = "cpe_qntl")[, term := .x]) %>%
                                  rbindlist) %>% rbindlist(idcol = "mdl")  %>%
    .[, cpe_qntl := factor(cpe_qntl, labels = c(`1` = "bottom 50%", `2` = "50%-90%", `3` = "top 10%"))]

  

  dt_coefs %>% 
    .[, term_pure := str_extract(term, paste0(c(gsub("exn_", "", l_treats), "no_treatment"), collapse = "|"))] %>%
    .[, dv := str_extract(mdl, paste0(gsub("exn_", "", l_treats), collapse = "|"))] %>% 
    ggplot(aes(x=estimate, y = term_pure, color = term_pure, shape = cpe_qntl, 
               xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.6)) +
    facet_wrap(~dv, scales = "free_x") +
    scale_color_manual(values = l_colors_treat) + 
    guides(color = guide_legend(title = "treatment"), nrow = 2) +
    theme(legend.position=  "bottom", legend.box = "vertical") +
    labs(y = "") +
    scale_shape_manual(values = c(1,2,5))

}


gp_dvauc <- function(dt_coef) {
  dt_coef[term != "no_treatment" & dv == "auc_turnover"] %>% gp_coef_did + facet_wrap(~term) +
    geom_hline(yintercept = 0, linetype = "dashed")
}

gp_dvauc_cum <- function(dt_coef) {
  dt_coef[term != "no_treatment" & dv == "aucturn_cum"] %>% gp_coef_did + facet_wrap(~term) +
    geom_hline(yintercept = 0, linetype = "dashed")
}


gp_weightdist <- function(dt_results) {
  
  ## FIXME: should be in global
  dt_results[, nbr_treats := rowSums(sign(.SD)), .SDcols = l_treats]

  dt_results %>% ggplot(aes(x=w_sbw, color = factor(nbr_treats))) + geom_density() +
    scale_x_continuous(trans = scales::asinh_trans(), breaks = c(0,1,3,5,10,15)) +
    theme(legend.position = "bottom")
  
}

gp_density_cprn <- function(dt_results) {

  l_cols_cpr <- c("exn_all_cum_ihs", "career_year", "aucturn_cum_ihs")

  ## FIXME: should be global
  dt_results[, no_treatment := 1- sign(rowSums(.SD)), .SDcols = l_treats]


  dt_dens_raw <- map(c(l_treats, "no_treatment") %>% setNames(.,.),
                     ~dt_results[!is.na(w_sbw) & get(.x) > 0,
                                            lapply(1:3, \(nbr)
                                                   density(.SD[[nbr]][get(l_cols_cpr[nbr]) > 0]) %$%
                                                   list(x=x,y=y, z = nbr)) %>%
                                            rbindlist, .SDcols = l_cols_cpr]) %>% rbindlist(idcol = "treat") %>%
    .[, src := "raw"]


  dt_dens_sbw <- map(c(l_treats, "no_treatment") %>% setNames(.,.),
                     ~dt_results[!is.na(w_sbw) & get(.x) > 0,
                                            lapply(1:3, \(nbr)
                                                   density(.SD[[nbr]][get(l_cols_cpr[nbr]) > 0],
                                                           weights = (w_sbw/sum(w_sbw))[get(l_cols_cpr[nbr]) > 0]) %$%
                                                   list(x=x,y=y, z = nbr)) %>%
                                            rbindlist, .SDcols = l_cols_cpr]) %>% rbindlist(idcol = "treat") %>%
    .[, src := "sbw"]


  dt_cbn <- rbind(dt_dens_raw, dt_dens_sbw) %>% 
    ## .[!(z ==3 & x < 5)] %>% 
    .[, max_y := max(y), .(z,src, treat)] %>%
    .[, sum_y := sum(y), .(z,src, treat)] %>% 
    .[, y := y/max_y]
  ## .[, y := y/max_y] 

  ## hmm this is still wrong but i can't be bothered to fix it rn

  dt_cbn %>% 
    ggplot(aes(x=x, height=y, y = treat, color = src)) + 
    geom_ridgeline(alpha = 0.3) +
    facet_wrap(~ z, scales = "free_x", nrow = 4, ncol = 3,
               labeller = as_labeller(setNames(l_cols_cpr, 1:3)))
  


}

  
gp_corw <- function(dt_results, sbw_config) {
  dt_bal_summary <- gd_sbw_bal_summary(dt_results, sbw_config, b_return_long = T)

  dt_bal_summary[, .(cor_raw = mean(cor_raw), cor_w = mean(cor_w)), .(treat, balcol)] %>%
    melt(id.vars = c("treat", "balcol")) %>%
    .[variable == "cor_w", variable := "cor_sbw"] %>%
    ## .[variable == "cor_w", treat := "weighted"] %>% 
    ggplot(aes(x=value, y = balcol, shape =variable, color = treat)) +
    geom_point(size = 2, position = position_dodge(width = 0.4)) +
    scale_shape_manual(values = c(1,2))

}

gp_sbw_raw_cprn <- function(dt_results) {

  dt_raw <- map(c(l_treats, "no_treatment") %>% setNames(.,.),
                ## setNames(l_treats, l_treats),
                ~dt_results[!is.na(w_sbw) & get(.x) > 0, colMeans(.SD) %>% as.list,
                            .SDcols = c("exn_all_cum", "career_year", "aucturn_cum_ihs")]) %>%
    rbindlist(idcol = "treat")
  
  dt_sbw <- map(c(l_treats, "no_treatment") %>% setNames(.,.),
                ## setNames(l_treats, l_treats),
                ~dt_results[!is.na(w_sbw) & get(.x) > 0, map(.SD, weighted.mean, w_sbw) %>% as.list,
                            .SDcols = c("exn_all_cum", "career_year", "aucturn_cum_ihs")]) %>%
    rbindlist(idcol = "treat")

  rbind(dt_raw[, src := "raw"], dt_sbw[, src := "sbw"]) %>%
    melt(id.vars = c("treat", "src")) %>%
    ggplot(aes(x=value, y = treat, color = treat, shape = src)) +
    geom_point(size = 3) +
    facet_wrap(~variable, scales = "free_x")


}

# Custom function to strip leading zeros
drop_leading_zero <- function(x) {
  # Format numbers cleanly and replace "0." with "."
  gsub("^0\\.", ".", label_number()(x))
}



gp_ess <- function(dt_results) {

  map(c(l_treats, "no_treatment") %>% setNames(.,.),
      ~dt_results[!is.na(w_sbw) & get(.x) > 0, .(w_sbw = sum(w_sbw), w_unif = .N,
                                                            ess = sum(w_sbw)^2/sum(w_sbw^2))]) %>%
    rbindlist(idcol = "treat") %>%
    .[, `:=`(ratio = ess/w_unif)] %>%
    melt(id.vars = "treat", measure.vars = c("w_unif", "ess", "ratio")) %>%
    .[, panel := ifelse(variable == "ratio", "ratio", "ess")] %>% 
    ggplot(aes(x=value, y =treat, color = treat, shape = variable)) +
    geom_point(size = 3) +
    facet_wrap(~panel, scales = "free_x") +
    scale_x_continuous(labels = drop_leading_zero)

}
