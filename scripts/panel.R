gf_foster <- function(dv, iv_type = "cum", transform = NULL, 
                      vars = gc_vvs()$l_orgforms,
                      fe = c("PeopleID", "year"), offset = NULL) {
    
    ivs <- sprintf("exhbN_%s_%s", vars, iv_type)
    if (!is.null(transform)) {
        ivs <- sprintf("%s(%s)", transform, ivs)
    }
    
    fe_str <- paste(fe, collapse = " + ")
    formula_str <- sprintf("%s ~ %s | %s", dv, paste(ivs, collapse = " + "), fe_str)
    
    if (!is.null(offset)) formula_str <- paste0(formula_str, " + offset(", offset, ")")
    return(as.formula(formula_str))
}

## * panel models




gl_mdls <- function(dt_atyear_foster, dt_atyear_ec) {
    gw_fargs(match.call())

    ## outcome auction turnover

    l_mdls <- list(
        r_auc = fepois(gf_foster("auc_turnover", iv_type = "cum", transform = "asinh"),
                       dt_atyear_foster, fixef.rm = "singleton",
                       cluster = "PeopleID"),

        r_auc_s5 = fepois(gf_foster("auc_turnover", iv_type = "s5", transform = "asinh"), dt_atyear_foster, fixef.rm = "singleton",
                       cluster = "PeopleID"),


        r_auc_ihs = fepois(gf_foster("asinh(auc_turnover)", iv_type = "cum", transform = "asinh"), dt_atyear_foster, fixef.rm = "singleton", cluster = "PeopleID"),

        ## outcome exhb in PI: model names describe DV/IV transformation: U (untransformed, I is IHS)
        r_pi_UU = fepois(gf_foster("exhbN_pi", iv_type = "cum"), dt_atyear_foster, fixef.rm = "none",
                       cluster = "PeopleID"),

        r_pi_IU = fepois(gf_foster("asinh(exhbN_pi)", iv_type = "cum", transform = "asinh"), dt_atyear_foster, fixef.rm = "none",
                       cluster = "PeopleID"),

        r_pi_UI = fepois(gf_foster("exhbN_pi", iv_type = "cum", transform = "asinh"), dt_atyear_foster, fixef.rm = "none",
                       cluster = "PeopleID"),

        r_pi_II = fepois(gf_foster("asinh(exhbN_pi)", iv_type = "cum", transform = "asinh"), dt_atyear_foster, fixef.rm = "none",
                       cluster = "PeopleID"),

        r_pi_UU_s5 = fepois(gf_foster("exhbN_pi", iv_type = "s5"), dt_atyear_foster, fixef.rm = "none",
                       cluster = "PeopleID"),

    ## ## outcome PES with interaction
    ## r_pes2 <- fepois(exhbN_pes ~ asinh(exhbN_gallery_cum) + asinh(exhbN_pes_cum)*asinh(exhbN_pm_cum) +
    ##                      asinh(exhbN_event_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

    ## r_pres3 <- fepois(exhbN_pes ~ asinh(exhbN_pes_cum)*asinh(exhbN_pm_cum) +
    ##                       asinh(exhbN_pm_cum)*asinh(exhbN_gallery_cum) + asinh(exhbN_event_cum)*
    ##                       asinh(exhbN_pm_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

    ## evc DV
        r_evc = feols(gf_foster("evc_ay_ttl", iv_type = "cum", transform = "asinh"), dt_atyear_ec, fixef.rm = "none")
        )
    



    ## l_mdls <- list(r_auc = r_auc,
    ##                r_pes = r_pes,
    ##                r_evc = r_evc)

    attr(l_mdls, "gnrtdby") <- as.character(match.call()[[1]])

    return(l_mdls)
}





gp_atyear_evc <- function(dt_atyear_evc) {
    gw_fargs(match.call())

    dt_atyear_evc %>% ggplot(aes(x=evc_ay_ttl + 0.000001)) +
        geom_density() + scale_x_continuous(trans = scales::log10_trans())
}

## gp_pminfl_nonpm <- function(rb1) {

##     tidy(rb1) %>% adt %>%
##         .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>% 
##         ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
##                    color = term)) +
##         geom_point(show.legend = F) +
##         geom_errorbarh(show.legend = F) + 
##         geom_vline(xintercept = 0, linetype = "dashed")
    
## }


## gp_pminfl_nonpm(rb1)


## gp_pminfl_auc <- function(rc1) {

##     tidy(rc1) %>% adt %>%
##         .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>% 
##         ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
##                    color = term)) +
##         geom_point(show.legend = F) +
##         geom_errorbarh(show.legend = F) + 
##         geom_vline(xintercept = 0, linetype = "dashed")
## }

## gp_pminfl_auc(rc1)


gp_pminfl_pes <- function(r_pes) {

    dt_ctgterm_lbls <- gc_vvs() %>% chuck("dt_ctgterm_lbls")

    l_facet_lbls <- c(pes = "Exhibitions (Public Exhb. Space)", auc = "Auction Turnover")
    
    l_cols <- c(pm = "#C77CFF", event = "#F8766D", gallery = "#7CAE00", pes = "#00BFC4")

    imap(list("pes" = r_pes),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        .[, src := factor(src, levels = c("pes", "auc"))] %>% 
        ## .[src == "survival", estimate := estimate *-1] %>%
        ## .[, src := factor(src, levels = c("exhibitions (pes)", "auction turnover"))] %>% 
        ## .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_s5)")[[1]], .I] %>%
        .[order(estimate)] %>%
        .[, term := factor(term, levels = term)] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbar(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        ## facet_grid(~src, scales = "free", labeller = as_labeller(l_facet_lbls)) +
        scale_y_discrete(labels = dt_ctgterm_lbls[, setNames(vlu_lbl, vlu)]) +
        labs(y=NULL, x = "impact") + 
        scale_color_manual(values = l_cols)

    


}

gp_pminfl_auc <- function(r_auc) {

    dt_ctgterm_lbls <- gc_vvs() %>% chuck("dt_ctgterm_lbls")

    l_facet_lbls <- c(pes = "Exhibitions (Public Exhb. Space)", auc = "Auction Turnover")
    
    l_cols <- c(pm = "#C77CFF", event = "#F8766D", gallery = "#7CAE00", pes = "#00BFC4")

    imap(list("auc" = r_auc),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        .[, src := factor(src, levels = c("pes", "auc"))] %>% 
        ## .[src == "survival", estimate := estimate *-1] %>%
        ## .[, src := factor(src, levels = c("exhibitions (pes)", "auction turnover"))] %>% 
        ## .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_s5)")[[1]], .I] %>%
        .[order(estimate)] %>%
        .[, term := factor(term, levels = term)] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbar(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        ## facet_grid(~src, scales = "free", labeller = as_labeller(l_facet_lbls)) +
        scale_y_discrete(labels = dt_ctgterm_lbls[, setNames(vlu_lbl, vlu)]) +
        labs(y=NULL, x = "impact") + 
        scale_color_manual(values = l_cols)

    


}



## gdplt("p_pminfl_pes")



gp_pminfl_cbn <- function(r_pes, r_auc) {

    dt_ctgterm_lbls <- gc_vvs() %>% chuck("dt_ctgterm_lbls")

    l_facet_lbls <- c(pes = "Exhibitions (Public Exhb. Space)", auc = "Auction Turnover")

    imap(list("pes" = r_pes, "auc" = r_auc),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        .[, src := factor(src, levels = c("pes", "auc"))] %>% 
        ## .[src == "survival", estimate := estimate *-1] %>%
        ## .[, src := factor(src, levels = c("exhibitions (pes)", "auction turnover"))] %>% 
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbar(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_grid(~src, scales = "free", labeller = as_labeller(l_facet_lbls)) +
        scale_y_discrete(labels = dt_ctgterm_lbls[, setNames(vlu_lbl, vlu)]) +
        labs(y=NULL, x = "coefficient")

}

## gp_pminfl_cbn(r_pes, r_auc)












gp_coef_trans <- function(l_mdls) {
    #' coefficients of different model configs

    ## specify iv-dv IHS types here, then coef plot, format DV/IV, either I (IHS) or U (untransformed)

    imap(c("UI", "II", "UU", "IU"), ~tidy(l_mdls[[paste0("r_pes", .y)]]) %>% adt %>% .[, src := .x]) %>% rbindlist %>% 
        .[, term := gsub("asinh\\(", "", term)] %>% .[, term := gsub("\\)", "", term)] %>%
        .[, `:=`(dv_trans = substring(src,0, 1), iv_trans = substring(src, 2,2))] %>%
        ## .[dv_trans == "U"] %>% 
        ggplot(aes(x=estimate, y = term, color = dv_trans, shape = iv_trans, xmin = estimate - 1.96*std.error, xmax = estimate + 1.96*std.error)) +
        geom_point(position = position_dodge(width = 0.4), size = 3, stroke = 3) +
        geom_errorbar(position = position_dodge(width = 0.4), linewidth = 1.5) +
        scale_shape_manual(values = c(1,2))  +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme(legend.position = "bottom")
}
    



gp_pred_exhbNcum <- function(dt_atyear_foster) {

    ## get succesful people who started young and have PM
    dt_sample <- dt_atyear_foster[first_year %between% c(2000,2002) & BirthYear > 1970,
                                  .SD[max(exhbN_pm_cum) > 2 & max(exhbN_all_cum) %between% c(80,120)], PeopleID]
        ## .[, exhbN_pes_pred := predict(l_mdls$r_pes1, .SD)]

    dt_sample[, uniqueN(PeopleID)]

    ## l_vpats_yeet <- c("s5", "n5", "cum", "_all", "yeet")

    ## dt_exhbN_long <- melt(dt_sample[PeopleID == 8742], id.vars = c("PeopleID", "year", "career_year"), measure.vars = patterns("exhbN_")) %>%
    ##     .[pall(map(l_vpats_yeet, ~!grepl(.x, variable)))] 
    ##     ## .[pall(map(c(), ~!grepl(.x, variable)))]
    ## dt_exhbN_long[, .N, variable] %>% print(n=80)

    ## i'm reinventing p_descs_panel.. but shittier? 

    ## dt_sample[PeopleID == 8742] %>% ggplot(aes(x = career_year, y = exhbN_all_cum)) +
    ##     geom_step() +
    ##     geom_step(dt_exhbN_long, mapping = aes(x=career_year, y=value, color = variable))


    ## melt(dt_sample[PeopleID == 8742], id.vars = c("PeopleID", "career_year", "year"), measure.vars = c("exhbN_pes", "exhbN_pes_pred")) %>% 
    ##     ggplot(aes(x=career_year, y = value, color = variable)) + geom_line()

    ## dt_sample %>% ggplot(aes(x=exhbN_pm_cum, y = exhbN_pes)) + geom_jitter() + geom_smooth(method = "lm")
    ## dt_sample %>% ggplot(aes(x=exhbN_gallery_cum, y = exhbN_pes)) + geom_jitter() + geom_smooth(method = "lm")



    ## dt_sample %>% ggplot(aes(x=exhbN_gallery_cum, y = exhbN_pes_pred, color = factor(PeopleID))) +
    ##     geom_point(show.legend = F) + geom_line(show.legend = F)

    ## dt_sample2 <- dt_sample[, .SD[sample(1:.N, 1)], PeopleID]

    ## map(0:20, ~dt_sample2 %>% copy %>% .[, exhbN_gallery_cum := exhbN_gallery_cum + .x]) %>% rbindlist %>%
    ##     .[, exhbN_pes_pred := predict(l_mdls$r_pes1, .SD)] %>% 
    ##     ggplot(aes(x=exhbN_gallery_cum, y = exhbN_pes_pred, color = factor(PeopleID))) +
    ##     geom_point(show.legend = F, size = 0.6) +
    ##     geom_line(show.legend = F) + 
    ##     scale_y_continuous(trans = scales::log10_trans()) +
    ##     scale_x_continuous(trans = scales::asinh_trans(), breaks = l_asinh_breaks)




    ## l_vrbls <- sprintf("exhbN_%s_cum", c("pm", "pes", "gallery", "event"))
    l_vrbls <- sprintf("exhbN_%s_s5", c("pm", "pes", "gallery", "event"))

    dt_pred <- map(l_vrbls,
                   ~lapply(0:10, \(t) dt_sample[PeopleID == 8742 & career_year == 0] %>% copy %>% 
                                      ## .[, `:=`(xxx = get(.x) + t, vrbl = .x)] %>%
                                      .[, `:=`(xxx = t, vrbl = .x)] %>%
                                      .[, .SD, .SDcols = setdiff(names(.), .x)] %>% # drop original value
                                      setnames(old = "xxx", new = .x) %>%
                                      .[, vx := get(.x)]
                           ) %>% rbindlist) %>% rbindlist(use.names= T) %>%
        setnafill(fil = 0, cols = l_vrbls) %>% 
        .[, .SD, .SDcols = c("PeopleID", "year", "vrbl", "vx", l_vrbls)] %>% 
        ## .[, pred := predict(l_mdls$r_pes1, .SD)] %>%
        .[, pred := predict(l_mdls$r_pes_s5, .SD)] %>%
        ## .[, vrbl := sub("exhbN_(.*?)_cum", "\\1", vrbl)]
        .[, vrbl := sub("exhbN_(.*?)_s5", "\\1", vrbl)]


    dt_lbls <- merge(dt_pred[, .SD[which.max(vx)], .(vlu = vrbl)],
                     gc_vvs() %>% chuck("dt_ctgterm_lbls"), by = "vlu") %>%
        .[, vx := vx+0.5] %>% 
        .[vlu == "gallery", pred := pred - 0.02]



    dt_pred %>%
        ggplot(aes(x=vx, y = pred, color = vrbl)) +
        geom_line(show.legend = F) +
        geom_point(show.legend = F) +
        ## facet_grid(~vrbl, scales = "free")  + 
        scale_x_continuous(trans = scales::asinh_trans(), breaks = l_asinh_breaks, limits = c(NA, 50))  +
        scale_y_continuous(trans = scales::log_trans()
                           ## breaks = c(, 1, 1.2, 1.4)
                           ) +  # breaks = c(0.5, 1,1.5,2)) +
        theme(legend.position = "bottom")  +
        geom_label(dt_lbls, mapping = aes(label = vlu_lbl, color = vlu), hjust = 0, show.legend = F) +
        labs(x="Number of Exhibitions", y = "Pred. Exhbs (Public Exhb Spaces)")
    
}




gp_grow_perc_arrow <- function(dt_insgrow5) {
    dt_perc <- dt_insgrow5[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)), .(vlu = InstitutionType), .SDcols = "exhbN_perc_n5s5_org"]
    
    dt_perc2 <- rbind(copy(dt_perc)[, x := "t5"],
                      copy(dt_perc)[, `:=`(x = "t0", exhbN_perc_n5s5_org = 0)]) %>%
        merge(gc_vvs() %>% chuck("dt_ctgterm_lbls"), by = "vlu") %>%
        .[, y2 := exhbN_perc_n5s5_org] %>%
        .[vlu == "event", y2 := y2 - 2.5]

    c_lbls_x <- c(t0 = "Before (5 year span)", t5 = "After (5 year span)")

    dt_perc2 %>% 
        ggplot(aes(x=x, y = exhbN_perc_n5s5_org, color = vlu, group = vlu)) +
        geom_point(show.legend = F) +
        geom_line(arrow = arrow(angle = 15, type = "closed"), show.legend = F) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_label(dt_perc2[x == "t5"], mapping = aes(label = vlu_lbl, y = y2),
                   position = position_nudge(x=0.4),
                   hjust = 0, show.legend = F) +
        geom_text(dt_perc2[x == "t5"],
                  mapping = aes(y=y2,label = sprintf("%s%%", format(exhbN_perc_n5s5_org, digits = 2))),
                  show.legend = F,
                  position = position_nudge(x=0.04), hjust = 0) +                   
        scale_x_discrete(expand = expansion(add = c(0.4,1.7)),
                         labels = c_lbls_x) +
        labs(x = NULL, y = "Percentage Change")
}

## * Plotting Coefs

gd_coefs <- function(l_mdls) {
    #' Extract coefficients from list of models
    imap(l_mdls, ~adt(tidy(.x))[, mdl := .y]) %>% rbindlist
}

gp_coefs <- function(dt_coefs, l_facet_lbls = NULL) {
    #' Plot coefficients
    
    dt_ctgterm_lbls <- gc_vvs() %>% chuck("dt_ctgterm_lbls")

    l_orgforms <- sprintf("_%s_", gc_vvs() %>% chuck("l_orgforms")) %>% paste0(collapse = "|")
        
    

    dt_coef_viz <- dt_coefs %>%
        .[grepl("_pi_", mdl)] %>%
        .[, instype := gsub("_", "", str_extract(term, l_orgforms))]
        
        
        

    p <- dt_coef_viz %>% 
        ggplot(aes(x = estimate, y = instype, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = mdl)) +
        geom_point() +
        geom_errorbarh(show.legend = FALSE) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_y_discrete(labels = dt_ctgterm_lbls[, setNames(vlu_lbl, vlu)]) +
        labs(y=NULL, x = "estimate") +
        facet_grid(~mdl, scales = "free_x")
    
    ## if (!is.null(l_facet_lbls)) {
    ##     p <- p + facet_grid(mdl ~ ., scales = "free_y", labeller = as_labeller(l_facet_lbls))
    ## } else {
    ##     p <- p + facet_grid(mdl ~ ., scales = "free_y")
    ## }
    
    return(p)
}
