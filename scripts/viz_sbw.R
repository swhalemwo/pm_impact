

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
  dt_coef[term != "no_treatment" & dv == "aucturn_cum2"] %>% gp_coef_did + facet_wrap(~term) +
    geom_hline(yintercept = 0, linetype = "dashed")
}

gp_auc_euro_evst <- function(dt_euro) {
  dt_euro[dv %in% c("auc_turnover_log", "aucturn_cum2_log", "aucturn_cum3_log" )] %>%
    .[, .(term, estimate = euro, std.error = se_euro, dv, tp, euro_hi, euro_lo)] %>% 
    .[, term_pure := str_extract(term, paste0(c(gsub("exn_", "", l_treats), "no_treatment"), collapse = "|"))] %>% 
    ggplot(aes(x=factor(tp), y = estimate, color = term_pure,
               ## ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error)) +
               ymin = euro_lo, ymax = euro_hi)) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_errorbar(position = position_dodge(width = 0.3)) +
    facet_wrap(~dv, scales = "free_y") +
    scale_color_manual(values = l_colors_treat) + 
    guides(color = guide_legend(title = "treatment")) +
    theme(legend.position=  "bottom") +
    labs(x = "timepoint relative to treatment")


}


gp_weightdist <- function(dt_results) {
  
  ## FIXME: should be in global
  dt_results[, nbr_treats := rowSums(sign(.SD)), .SDcols = l_treats]

  dt_results %>% ggplot(aes(x=w_sbw, color = factor(nbr_treats))) + geom_density() +
    scale_x_continuous(trans = scales::asinh_trans(), breaks = c(0,1,3,5,10,15)) +
    theme(legend.position = "bottom")
  
}

gp_density_cprn <- function(dt_results) {

  l_cols_cpr <- c("exn_all_cum_log", "career_year", "aucturn_cum2_log")

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
                            .SDcols = c("exn_all_cum", "career_year", "aucturn_cum2_log")]) %>%
    rbindlist(idcol = "treat")
  
  dt_sbw <- map(c(l_treats, "no_treatment") %>% setNames(.,.),
                ## setNames(l_treats, l_treats),
                ~dt_results[!is.na(w_sbw) & get(.x) > 0, map(.SD, weighted.mean, w_sbw) %>% as.list,
                            .SDcols = c("exn_all_cum", "career_year", "aucturn_cum2_log")]) %>%
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


gp_cpr_range_raw_sbw <- function(dt_results, c_sbw) {

  l_balcols <- gl_sbw_balcols(c_sbw)
  
  l_probs <- c(0.25, 0.5, 0.75)

  dt_sbw <- dt_results[, map(.SD, ~wtd.quantile(.x, weights = w_sbw, probs = l_probs)), .SDcols = l_balcols] %>% t %>%
    adt(keep.rownames = "vrbl") %>% .[, src := "sbw"] %>% 
    setnames(old = paste0("V", 1:3), new = paste0("vlu_", c("low", "med", "high")))

  

  dt_raw <- dt_results[, map(.SD, ~quantile(.x, probs = l_probs)), .SDcols = l_balcols] %>% t %>%
    adt(keep.rownames = "vrbl") %>% .[, src := "raw"] %>% 
    setnames(old = paste0("V", 1:3), new = paste0("vlu_", c("low", "med", "high")))


  rbind(dt_sbw, dt_raw) %>%
    ggplot(aes(x = vlu_med, y = vrbl, color = src, xmin = vlu_low, xmax = vlu_high)) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_errorbar(position = position_dodge(width = 0.3))
    ## facet_wrap(~vrbl, scales = "free")
 
}

gp_cpr_density_rawsbw <- function(dt_results, c_sbw) {

  l_balcols <- gl_sbw_balcols(c_sbw)

  dt_raw <- dt_results[, imap(.SD, ~density(.x) %$% list(x=x,y=y, vrbl = .y)) %>%
                         rbindlist, .SDcols = l_balcols] %>% .[, src := "raw"] 

  ## custom bandwidth
  ## my_bw <- 0.9 * min(weighted_sd, IQR(.x)/1.34) * length(.x)^(-1/5) # Weigh

  dt_sd_w <- dt_results[, imap(.SD, ~list(sd_w = sqrt(Hmisc::wtd.var(.x, weights = w_sbw)),
                                         iqr_adj = IQR(.x)/1.34, lenx = .N, vrbl = .y)) %>% rbindlist, .SDcols = l_balcols]

  dt_sd_w[, min_sdiqr := pmin(sd_w, iqr_adj)] %>%
    .[, bw := 0.9 * min_sdiqr * lenx^(-1/5)]

  
  dt_sbw <- dt_results[, imap(.SD, ~density(.x, weights = w_sbw/sum(w_sbw), bw = dt_sd_w[vrbl == .y, bw]) %$% list(x=x,y=y, vrbl = .y)) %>%
                         rbindlist, .SDcols = dt_sd_w[bw > 0, vrbl]] %>%
    .[, src := "sbw"] 
    

  dt_cbn <- rbind(dt_raw, dt_sbw)

  ## dt_cbn[, .(mean_y = mean(y)), .(vrbl, src)] %>%
  ##   dcast(vrbl ~ src, value.var = "mean_y")

  dt_cbn %>%
    .[y > 0.01] %>% 
    ggplot(aes(x=x, y=y, color = src)) +
    ## geom_line() +
    geom_smooth(se = F, span = 0.5) + 
    facet_wrap(~vrbl, scales = "free")
  



}


gp_loweights <- function(dt_results, c_sbw) {

  ## library(mgcv)

  l_balcols <- gl_sbw_balcols(c_sbw)

  ## dt_results[, .(low_w = mean(w_sbw)), .(exn_pi_l1)] %>%
  ##   ggplot(aes(x=exn_pi_l1, y = low_w)) + # geom_line() +
  ##   geom_smooth(se = F, span = 0.1)

  

  ## l_r_gam <- map(l_balcols[1:21] %>% setNames(.,.),
  ##                ~gam(as.formula(sprintf("w_sbw ~ s(%s, bs = 'cs')", .x)), data = dt_results, method = "REML"))

  ## l_pred <- imap(l_balcols[1:21], ~setnames(data.table(x = unique(dt_results[, get(.x)])), old = "x", new = .x) %>%
  ##                                  .[, pred := predict(l_r_gam[[.x]], newdata = .SD), .SDcols = .x] %>%
  ##                                  .[, .(x = get(.x), src = .x, pred)])
                                   
  ## dt_pred <- rbindlist(l_pred)

  ## dt_pred %>% ggplot(aes(x=x, y=pred)) +
  ##   geom_line() +
  ##   geom_smooth(se = F) + 
  ##   facet_wrap(~src, scales = "free")

  ## idk predicton is just so annoying
  ## 

  dt_cbn <- map(l_balcols, ~dt_results[, .(mean_w = mean(w_sbw), vrbl = .x, .N), .(x = get(.x))]) %>% rbindlist

  dt_cbn %>%
    .[!grepl("gender_rcd", vrbl)] %>% 
    ggplot(aes(x=x, y=mean_w, weight = N)) + # geom_point() +
    facet_wrap(~vrbl, scales = "free") + geom_smooth(se = F, span = 0.3)



  
}

weighted_sd <- function(x,w) sqrt(Hmisc::wtd.var(x,w))

gt_sumstats <- function(r_sbw) {
  #' generate summary table

  l_vrbls_tosum <- c(r_sbw$config$treatment_cols, r_sbw$config$covariate_cols)

  ## functions to aggregate
  l_fs_nowtd <- c(max = max, min = min)
  l_fs_uwtd <- c(mean = mean, sd = sd)
  l_fs_wtd <- c(mean = weighted.mean, sd = weighted_sd)

  dt_res2 <- cbind(r_sbw$dt_res, r_sbw$dt_panel[, .(gender_rcd)])[!is.na(w_sbw)]

  ## summarize those where weighting makes no sense (min, max)
  dt_sum_nowtd <- CJ(f = names(l_fs_nowtd), v = l_vrbls_tosum)
  dt_sum_nowtd[, vlu := apply(dt_sum_nowtd, 1, \(x) l_fs_nowtd[[x[["f"]]]](dt_res2[, get(x[["v"]])]))]

  dt_sum_wtd <- CJ(f = names(l_fs_uwtd), v = l_vrbls_tosum)
  dt_sum_wtd[, vlu_uw := apply(dt_sum_wtd, 1, \(x) l_fs_uwtd[[x[["f"]]]](dt_res2[, get(x[["v"]])]))]
  # idk why this no work
  dt_sum_wtd[, vlu_wtd := apply(dt_sum_wtd, 1, \(x) l_fs_wtd[[x[["f"]]]](dt_res2[, get(x[["v"]])],
    dt_res2[, w_sbw]))]

  dt_sum_wtd2 <- melt(dt_sum_wtd, id.vars = c("f", "v"), value.vars = c("vlu_uw", "vlu_wtd"),
                      variable.name = "b_w")

  ## weighted.mean(dt_res2[, exn_all_cum_log], dt_res2[, w_sbw])
  ## mean(dt_res2[, exn_all_cum_log])
  ## ## huh means are the same: i think this is maintained by SBW
  ## ## but variances are not 
  

  ## rbind(dt_res2[, .(PeopleID, year, w_sbw, exn_all_cum_log, src = "sbw")],
  ##       dt_res2[, .(PeopleID, year, w_sbw = 1, exn_all_cum_log, src = "raw")]) %>% 
  ##   ## .[, weighted.mean(exn_all_cum_log, w_sbw), src] 
  ##   ggplot(aes(x=exn_all_cum_log, color = src, weight = w_sbw)) + geom_density()
                
  dt_sumstats_viz <- join(
    dcast(dt_sum_wtd2, v ~ f + b_w , value.var = "value"),
    dcast(dt_sum_nowtd, v ~ f, value.var = "vlu"),
    on = "v", verbose = 0) %>%
    setnames(old = "v", new = "vrbl")

  dt_vrblinfo <- gc_vvs() %>% chuck("dt_vrblinfo")

  ## order columns as the spanner reads them (min before max, CJ() sorts f alphabetically)
  dt_sumstats_viz2 <- dt_vrblinfo[, .(vrbl, vrbl_lbl)][dt_sumstats_viz, on  = "vrbl"]

  ## select order
  ## first get 
  l_orgforms <- gc_vvs() %>% chuck("l_orgforms")
  l_orgforms_exn <- paste0("exn_", l_orgforms)
  l_orgforms_cum <- keep(setdiff(dt_sumstats_viz2[, vrbl], l_orgforms_exn), ~grepl("exn", .x))
  l_atsubreg <- dt_sumstats_viz2[grepl("atsubreg", vrbl), vrbl]
  l_rest <- setdiff(dt_sumstats_viz2[, vrbl], c(l_orgforms_exn, l_orgforms_cum, l_atsubreg))

  dt_sumstats_viz3 <- map(list(l_orgforms_exn, l_orgforms_cum, l_rest, l_atsubreg),
                          ~dt_sumstats_viz2[vrbl %in% .x]) %>% rbindlist %>%   
                      .[, vrbl := NULL] %>%
                      setcolorder(c("vrbl_lbl", "mean_vlu_uw", "mean_vlu_wtd", "sd_vlu_uw", "sd_vlu_wtd",
                                    "min", "max"))
  

  c_names <- c(vrbl_lbl = "",
               mean_vlu_uw = "Unweighted",
               mean_vlu_wtd = "Weighted",
               sd_vlu_uw = "Unweighted",
               sd_vlu_wtd = "Weighted",
               max = "",
               min = "")
                 
  c_colnames <- gc_colnames(names(dt_sumstats_viz3), c_names, hline_above = F,
                            align = as.list(c("l", rep("r", 6))))

  c_spanner <- gc_spanner(c("Variable", "Mean", "SD", "Min", "Max"), c(1,2,2,1,1))


  c_atr <- list(
    pos = list(-1,0),
    command = c(c_spanner, c_colnames))


  list(
    dt_fmtd = dt_sumstats_viz3,
    align_cfg = c("l", "l", rep("r", 6)),  # xtable: rownames + one entry per column
    hline_after = c(0, dt_sumstats_viz3[, .N]),
    add_to_row = c_atr,
    number_cols = c(F, rep(T, 6)),  # length must equal ncol (subset by position)
    landscape = T)


}


## countrycode(c(145, 151, 154, 155, 21 , 30 , 39 , 419, 53), "un.regionsub.code", "un.regionsub.name")

## dtx <- dt_atyear_foster[, .N, Nationality] %>%
##   .[, `:=`(regionsub.code = countrycode(Nationality, "country.name", "un.regionsub.code"),
##            regionsub.name = countrycode(Nationality, "country.name", "un.regionsub.name"))] %>%
##   .[, .(N = sum(N)), .(regionsub.code, regionsub.name)] %>% .[order(-N)]


ess_w <- function(w) sum(w)^2 / sum(w^2)


fmt_cor <- function(x, digits = 2L) {
  #' format correlation cells for a table
  #' @param x      numeric vector
  #' @param digits number of decimals (scalar)
  #' @return character vector, "--" where x is NA; rounded first so that a
  #'   tiny negative correlation prints "0.000" rather than "-0.000"
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), round(x, digits)))
}


gd_sbw_vrbl_lbl <- function(l_vrbls, dt_vrblinfo = gc_vvs() %>% chuck("dt_vrblinfo")) {
  #' labels for SBW table rows
  #'
  #' Takes the project variable label (gc_vvs) where a variable is known, and
  #' adds the lag distance for the pipeline-generated `_l<k>` columns, so that
  #' `career_year_l2` reads "Career Year (t-2)" and `exn_all_cum2_log_l1` reads
  #' "Prior Exhibitions (all types, log) (t-1)". Unknown columns keep their
  #' column name.
  #'
  #' @param l_vrbls     column names, in table order
  #' @param dt_vrblinfo variable info table (vrbl, vrbl_lbl), from gc_vvs()
  #' @return character vector of labels, same length/order as l_vrbls

  data.table(vrbl = l_vrbls) %>%
    .[, base := sub("_l[0-9]+$", "", vrbl)] %>%
    .[, dist := sub("^.*_l([0-9]+)$", "\\1", vrbl)] %>%
    join(dt_vrblinfo[, .(base = vrbl, base_lbl = vrbl_lbl)], on = "base", how = "left", verbose = 0) %>%
    .[, base_lbl := fifelse(is.na(base_lbl), base, base_lbl)] %>%
    .[, fifelse(vrbl != base, sprintf("%s (t-%s)", base_lbl, dist), base_lbl)]

}


gd_sbw_bal_rows <- function(c_sbw) {
  #' row metadata for the balance table
  #'
  #' The row set is the constraint set itself (gl_sbw_balcols), so the table
  #' cannot drift from what the weights actually balance: treatment lags, then
  #' covariates, then the extra balancing columns (in that order, as returned by
  #' gl_sbw_balcols).
  #'
  #' @param c_sbw config (needs the FINAL covariate_cols / balance_cols, i.e.
  #'   the config stored in a gl_sbw_proc result)
  #' @return data.table (balcol, row_lbl, grp_lbl), in row order

  c_sbw <- md_sbw_defaults(c_sbw)

  l_bal <- gl_sbw_balcols(c_sbw)

  ## three blocks: treatment history, covariates, extra balance columns
  l_grp <- fifelse(l_bal %in% gd_sbw_shiftcols(c_sbw$treatment_cols, c_sbw$n_lag, "l"),
                   "Exposure history (lags)",
                   fifelse(l_bal %in% c_sbw$balance_cols,
                           "Additional balance columns",
                           "Covariates"))

  data.table(balcol = l_bal, row_lbl = gd_sbw_vrbl_lbl(l_bal), grp_lbl = l_grp)
}


gd_sbw_bal_desc <- function(r_sbw, v_years = NULL) {
  #' group descriptives for the SBW report
  #'
  #' Groups are `exn_X > 0` for each treatment, `none` (no exposure of any type)
  #' and `all`. Means are UNWEIGHTED on purpose: because the weights drive the
  #' treatment-covariate correlation to (about) zero, every weighted group mean
  #' collapses onto the overall mean, so a weighted version of this table
  #' carries no treatment-specific information. The weighted side of the story
  #' is the correlation table (gd_sbw_bal_cor).
  #'
  #' All group means are computed in one go as a (group x variable) crossproduct
  #' divided by the group sizes — no per-variable loop.
  #'
  #' @param r_sbw   gl_sbw_proc output
  #' @param v_years optional year subset; NULL = every year that has weights
  #' @return list(dt_mean = long (grp, vrbl, vlu),
  #'              dt_supp = (grp, n, n_artists, ess_ratio))

  c_sbw <- md_sbw_defaults(r_sbw$config)
  dt_res <- r_sbw$dt_results[!is.na(get(c_sbw$wcol))]
  if (!is.null(v_years)) dt_res <- dt_res[get(c_sbw$tid) %in% v_years]

  ## group membership: positive dose per type, none, all
  m_trt <- as.matrix(dt_res[, .SD, .SDcols = c_sbw$treatment_cols])
  m_grp <- cbind(sign(m_trt) > 0, rowSums(sign(m_trt)) == 0, rep(TRUE, nrow(m_trt)))
  colnames(m_grp) <- c(c_sbw$treatment_cols, "none", "all")
  l_grps <- colnames(m_grp)

  ## group means: (group x variable) crossproduct / group size
  m_vlu <- as.matrix(dt_res[, .SD, .SDcols = c(c_sbw$treatment_cols, c_sbw$covariate_cols)])
  m_mean <- crossprod(m_grp * 1L, m_vlu) / colSums(m_grp)

  dt_mean_wtd <- map(l_grps, ~dt_res[which(m_grp[,.x] == T), lapply(.SD, weighted.mean, w_sbw),
                      .SDcols = c(c_sbw$treatment_cols, c_sbw$covariate_cols)]) %>% rbindlist %>%
    .[, grp := l_grps] %>%
    melt(id.vars = "grp", value.vars = c(c_sbw$treatment_cols, c_sbw$covariate_cols),
         variable.name = "vrbl", value.name = "vlu") %>%
    .[, src := "wtd"]

  ## dt_res


  dt_mean_uwtd <- adt(m_mean, keep.rownames = "grp") %>%
    melt(id.vars = "grp", variable.name = "vrbl", value.name = "vlu", variable.factor = FALSE) %>%
    .[, grp := factor(grp, levels = l_grps)] %>%
    .[, src := "uwtd"]

  dt_mean <- rbind(dt_mean_wtd, dt_mean_uwtd)


  ## support: group size, artists (rows repeat within artist), ESS within group
  w <- dt_res[[c_sbw$wcol]]
  v_nart <- map(l_grps, ~dt_res[m_grp[, .x], uniqueN(get(c_sbw$uid))]) %>% unlist
  v_ess <- map(l_grps, ~ess_w(w[m_grp[, .x]])) %>% unlist

  dt_supp <- data.table(grp = l_grps, n = colSums(m_grp), n_artists = v_nart, ess = v_ess,
                        ess_ratio = v_ess / colSums(m_grp))

  list(dt_mean = dt_mean, dt_supp = dt_supp)
}


gd_sbw_bal_cor <- function(r_sbw, v_years = NULL) {
  #' balance correlations for the SBW report table
  #'
  #' The per-year raw and weighted treatment-balance correlations (one pair per
  #' year, treat x balancing column) are averaged with EQUAL WEIGHT over years —
  #' the same way the pooled DiD is an equal-weight fixed-effects average of the
  #' per-year fits. A correlation is not pooled across the stacked panel on
  #' purpose: the covariate distribution drifts over calendar time, so a pooled
  #' correlation would mix within- and between-year variation.
  #'
  #' @param r_sbw   gl_sbw_proc output
  #' @param v_years optional year subset; NULL = every year that has weights
  #' @return list(dt_avg = (treat, balcol, cor_raw, cor_w),  # mean over years
  #'              dt_max = (treat, max_raw, max_w))           # worst year

  c_sbw <- md_sbw_defaults(r_sbw$config)
  dt_long <- gd_sbw_bal_summary(r_sbw$dt_results, c_sbw, v_years = v_years, b_return_long = TRUE)

  ## mean over years; NA where a balancing column had no variance in a year
  dt_avg <- dt_long[, .(cor_raw = mean(cor_raw, na.rm = TRUE),
                        cor_w = mean(cor_w, na.rm = TRUE)), .(treat, balcol)] %>%
    .[is.nan(cor_raw), cor_raw := NA] %>%
    .[is.nan(cor_w), cor_w := NA]

  ## worst case over years and balancing columns, per treatment
  dt_max <- dt_long[, .(max_raw = suppressWarnings(max(abs(cor_raw), na.rm = TRUE)),
                        max_w = suppressWarnings(max(abs(cor_w), na.rm = TRUE))), treat] %>%
    .[!is.finite(max_raw), max_raw := NA] %>%
    .[!is.finite(max_w), max_w := NA]

  list(dt_avg = dt_avg, dt_max = dt_max)
}


gt_sumstats_bytreat <- function(r_sbw, v_years = NULL) {
  #' descriptive table: analysis variables by treatment group (unweighted)
  #'
  #' Rows are the variables (exposures at t, then covariates), columns the groups:
  #' all artist-years, one column per treatment (positive dose), and the
  #' no-exposure group. The footer carries the group size, the number of distinct
  #' artists (rows repeat within artist) and the ESS as a share of the group size.
  #'
  #' @param r_sbw   gl_sbw_proc output
  #' @param v_years optional year subset; NULL = every year that has weights
  #' @return jtls table object (see wtbl)

  c_sbw <- md_sbw_defaults(r_sbw$config)
  l_res <- gd_sbw_bal_desc(r_sbw, v_years = v_years)

  ## row metadata: exposures at t, then covariates, in config order
  dt_rows <- data.table(vrbl = c(c_sbw$treatment_cols, c_sbw$covariate_cols)) %>%
    .[, grp_lbl := c(rep("Exposure at t", length(c_sbw$treatment_cols)),
                 rep("Covariates", length(c_sbw$covariate_cols)))] %>%
    .[, row_lbl := gd_sbw_vrbl_lbl(vrbl)]

  ## wide: rows x groups (all, treatments, none), formatted as means
  l_grps <- c("all", c_sbw$treatment_cols, "none")

  dt_fmt <- l_res$dt_mean %>% copy %>% 
    .[, vlu_fmt := sprintf("%.2f", vlu)] %>%
    dcast(vrbl + grp ~ src , value.var = "vlu_fmt") %>%
    .[, cell_fmt := sprintf("%s (%s)", uwtd, wtd)] %>%
    dcast(vrbl ~ grp, value.var = "cell_fmt") %>%
    .[dt_rows, on = "vrbl"] %>% # reorder rows to the config order
    .[, .SD, .SDcols = c("row_lbl", "grp_lbl", l_grps)]

  ## l_grps2 <- outer(l_grps, c("_uwtd", "_wtd"), paste0) %>% t %>% as.vector

  ## dt_fmt <- dcast(l_res$dt_mean, vrbl ~ grp + src, value.var = "vlu", fill = NA) %>% 
  ##   ## setcolorder(paste0(c("vrbl", l_grps),) %>%
  ##   .[, (l_grps2) := map(.SD, ~sprintf("%.2f", .x)), .SDcols = l_grps2] %>%
  ##   .[dt_rows, on = "vrbl"] %>% # reorder rows to the config order
  ##   .[, .SD, .SDcols = c("row_lbl", "grp_lbl", l_grps2)]

  ## footer: group size, artists, ESS share
  dt_foot <- copy(l_res$dt_supp) %>%
    .[, `:=`(n = format(n, big.mark = " ", trim = TRUE),
             n_artists = format(n_artists, big.mark = " ", trim = TRUE),
             ess = format(round(ess), big.mark = " ", trim = TRUE),
             ess_ratio = sprintf("%.1f", 100 * ess_ratio))] %>%
    melt(id.vars = "grp", variable.name = "stat", value.name = "vlu") %>%
    .[, stat := factor(stat, levels = c("n", "n_artists", "ess", "ess_ratio"))] %>%
    dcast(stat ~ grp, value.var = "vlu") %>%
    .[order(stat)] %>%
    .[, `:=`(row_lbl = c("N (artist-years)", "N (artists)", "ESS", "ESS / N (\\%)")[stat], grp_lbl = "")] %>%
    .[, .SD, .SDcols = c("row_lbl", "grp_lbl", l_grps)]

  dt_viz <- rbind(dt_fmt, dt_foot, fill = TRUE) %>%
    .[, (l_grps) := map(.SD, ~fifelse(is.na(.x), "", .x)), .SDcols = l_grps]

  dt_out <- dt_viz[, .SD, .SDcols = setdiff(names(dt_viz), "grp_lbl")] # grp_lbl is only for the row-block labels

  ## header: reference / treatment group spanner over the column labels
  l_collbls <- c(list(row_lbl = "Variable", all = "All"),
                 setNames(as.list(toupper(gsub("^exn_", "", c_sbw$treatment_cols))), c_sbw$treatment_cols),
                 list(none = "None"))

  c_colnames <- gc_colnames(names(dt_out), l_collbls, hline_above = F,
                            align = as.list(c("l", rep("r", length(l_grps)))))
  c_spanner <- gc_spanner(c("", "", "Exhibition type (positive dose)", ""),
                          c(1, 1, length(c_sbw$treatment_cols), 1))

  ## row block labels (variable rows only, the footer has no block)
  dt_grpstrs <- gc_grpstrs(dt_fmt, grp = "grp_lbl", nbr_cols = ncol(dt_out))

  c_atr <- list(
    pos = c(list(-1, 0, 0), dt_grpstrs$pos),
    command = c(c_spanner, c_colnames, "\\hline \n", dt_grpstrs$grpstr))

  list(
    dt_fmtd = dt_out,
    align_cfg = c("l", "l", rep("r", length(l_grps))), # xtable: rownames + one entry per column
    hline_after = c(nrow(dt_fmt), nrow(dt_viz)), # before the footer, and at the bottom
    add_to_row = c_atr,
    number_cols = c(F, rep(T, length(l_grps))), # length must equal ncol (subset by position)
    landscape = T)

}


gt_bal_bytreat <- function(r_sbw, v_years = NULL, l_digits = NULL) {
  #' balance table: raw vs weighted treatment-balance correlations
  #'
  #' One row per balancing column — treatment lags, covariates and the extra
  #' balancing columns, i.e. exactly the constraint set (gd_sbw_bal_rows) — and
  #' one column per treatment, with the unweighted and the weighted correlation
  #' side by side. Cells are the equal-weighted mean over years (gd_sbw_bal_cor),
  #' the footer holds the worst year x column |correlation| per treatment.
  #'
  #' @param r_sbw     gl_sbw_proc output
  #' @param v_years   optional year subset; NULL = every year that has weights
  #' @param l_digits  optional named vector of decimals per measure
  #'   ("unw"/"wtd"); NULL = 2 unweighted, 3 weighted
  #' @return jtls table object (see wtbl)

  c_sbw <- md_sbw_defaults(r_sbw$config)
  l_res <- gd_sbw_bal_cor(r_sbw, v_years = v_years)
  dt_rows <- gd_sbw_bal_rows(c_sbw)
  l_digits <- l_digits %||% c(unw = 2L, wtd = 3L)

  ## long -> wide: treatment x (unweighted, weighted)
  l_cols <- map(c_sbw$treatment_cols, ~paste0(.x, c("_unw", "_wtd"))) %>% unlist(use.names = FALSE)

  dt_wide <- melt(l_res$dt_avg, id.vars = c("treat", "balcol"),
                  measure.vars = c("cor_raw", "cor_w"), variable.name = "b_w",
                  value.name = "vlu") %>%
    .[, b_w := fcase(b_w == "cor_raw", "unw", default = "wtd")] %>%
    dcast(balcol ~ treat + b_w, value.var = "vlu", fill = NA) %>%
    setcolorder(c("balcol", l_cols))

  ## format before any row reordering (the format step indexes the columns)
  l_dig_cols <- l_digits[sub("^.*_(unw|wtd)$", "\\1", l_cols)] %>% setNames(l_cols)

  dt_fmt <- copy(dt_wide) %>%
    .[, (l_cols) := map(l_cols, ~fmt_cor(dt_wide[[.x]], l_dig_cols[[.x]]))] %>%
    .[dt_rows, on = "balcol"] %>% # order by the constraint set, add labels
    .[, .SD, .SDcols = c("row_lbl", "grp_lbl", l_cols)]

  ## footer: worst case over years and balancing columns
  dt_foot <- melt(l_res$dt_max, id.vars = "treat", variable.name = "b_w", value.name = "vlu") %>%
    .[, b_w := fcase(b_w == "max_raw", "unw", default = "wtd")] %>%
    .[, vlu := fmt_cor(vlu, 3L)] %>%
    dcast(. ~ treat + b_w, value.var = "vlu") %>%
    setcolorder(l_cols) %>%
    .[, `:=`(row_lbl = "Max. abs. cor across years", grp_lbl = "")] %>%
    .[, .SD, .SDcols = c("row_lbl", "grp_lbl", l_cols)]

  dt_viz <- rbind(dt_fmt, dt_foot, fill = TRUE)

  dt_out <- dt_viz[, .SD, .SDcols = setdiff(names(dt_viz), "grp_lbl")] # grp_lbl is only for the row-block labels

  ## header: treatment spanner over each unw/wtd pair
  l_collbls <- c(list(row_lbl = "Variable"),
                 setNames(as.list(rep(c("raw", "SBW"), times = length(c_sbw$treatment_cols))),
                          l_cols))

  c_colnames <- gc_colnames(names(dt_out), l_collbls, hline_above = F,
                            align = as.list(c("l", rep("r", length(l_cols)))))
  c_spanner <- gc_spanner(c("", toupper(gsub("^exn_", "", c_sbw$treatment_cols))),
                          c(1, rep(2, length(c_sbw$treatment_cols))))

  ## row block labels (balancing rows only, the footer has no block)
  dt_grpstrs <- gc_grpstrs(dt_fmt, grp = "grp_lbl", nbr_cols = ncol(dt_out))

  c_atr <- list(
    pos = c(list(-1, 0, 0), dt_grpstrs$pos),
    command = c(c_spanner, c_colnames, "\\hline \n", dt_grpstrs$grpstr))

  list(
    dt_fmtd = dt_out,
    align_cfg = c("l", "l", rep("r", length(l_cols))), # xtable: rownames + one entry per column
    hline_after = c(nrow(dt_fmt), nrow(dt_viz)), # before the footer, and at the bottom
    add_to_row = c_atr,
    number_cols = c(F, rep(T, length(l_cols))), # length must equal ncol (subset by position)
    landscape = T)
}


## map(c(l_treats, "no_treatment"),
##     ~r_sbw$dt_results[get(.x) > 0 & !is.na(w_sbw),
##                       .(mean_uw = mean(career_year), mean_wtd = weighted.mean(career_year, w_sbw),
##                         mean_twtd = weighted.mean(career_year, w_sbw*get(.x)),
##                         treat = .x)]) %>% rbindlist

