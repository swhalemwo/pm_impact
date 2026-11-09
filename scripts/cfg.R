
gc_plts <- function() {

    list(
        ## p_pminfl_cbn = list(
        ##     r_pes = quote(l_mdls$r_pi_UI),
        ##     r_auc = quote(l_mdls$r_auc),
        ##     width = 14, height = 7,
        ##     caption = "influence of PM exhibitions"),
        ## p_pminfl_pi = list(
        ##     r_pes = quote(l_mdls$r_pi_UI),
        ##     ## r_pes = quote(l_mdls$r_pes1),
        ##     width = 9, height = 7,
        ##     caption = "influence of PM exhibitions"),
        ## p_pminfl_auc = list(
        ##     r_auc = quote(l_mdls$r_auc_s5),
        ##     ## r_pes = quote(l_mdls$r_pes1),
        ##     width = 13, height = 7,
        ##     caption = "influence of PM exhibitions"),
        ## p_hist_crap = list(
        ##     dt_atyear_foster = quote(dt_atyear_foster),
        ##     width = 14, height = 7,
        ##     caption = "nonsense hist"),
        p_activity_jitter = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height =7,
            caption = "activity jitter"),
        ## p_insgrow5 = list(
        ##     dt_insgrow5 = quote(dt_insgrow5),
        ##     width = 14, height =9,
        ##     caption = "growth relations"),
        ## p_insgrow5_contours = list(
        ##     dt_insgrow5 = quote(dt_insgrow5),
        ##     width = 14, height =9,
        ##     caption = "growth relations"),
        p_cumad_jitter = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "cumulative advantage whatever"),
        ## p_insgrow5_dist = list(
        ##     dt_insgrow5 = quote(dt_insgrow5),
        ##     l_cols_summary = quote(l_cols_summary),
        ##     width = 14, height = 18, caption = "grow5 dists"),
        p_mainseq = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "cumulative advantage"),
        p_dist_exn_pi = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "distribution of yearly PES exhibitions"),
        p_dist_aucturn = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "distribution of yearly aucturn turnover"),
        ## p_coef_trans = list(
        ##     l_mdls = quote(l_mdls),
        ##     width = 14, height = 10,
        ##     caption = "coefs"),
        ## p_pred_exncum = list(
        ##     dt_atyear_foster = quote(dt_atyear_foster),
        ##     width = 14, height = 7,
        ##     caption = "coefs"),
        p_evst = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 12,
            caption = "Event study plot"),
        p_evst_wintcpt = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 12,
            caption = "Event study plot with no-treatment group"),
        p_evst_cum = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 12,
            caption = "Event study plot (cumulative DV)"),
        p_evst_cum_wintcpt = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 12,
            caption = "Event study plot (cumulative DV) with no-treatment group"),
        p_compact = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 12,
            caption = "Event study coefficients (cumulative) at t=6"),
        p_cpct_intr_cpeqntl = list(
            dt_results = quote(r_sbw$dt_results),
            width = 14, height = 12,
            caption = "Event study coefficients (cumulative) at t=6 interacted by prior exhibition count"),
        p_dvauc = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 10,
            caption = "Impact of exhibitions on auction turnover"),
        p_dvauc_cum = list(
            dt_coef = quote(dt_coef_did),
            width = 14, height = 10,
            caption = "Impact of exhibitions on auction turnover (cumulative)"),
        p_weightdist = list(
            dt_results = quote(r_sbw$dt_results),
            width = 14, height = 10,
            caption = "weight distribution by number of distinct treatments"),
        p_density_cprn = list(
            dt_results = quote(r_sbw$dt_results),
            width = 14, height = 10,
            caption = "covariate distribution by treatment before and after SBW"),
        p_corw = list(
            dt_results = quote(r_sbw$dt_results),
            sbw_config = quote(r_sbw$config),
            width = 14, height = 10,
            caption = "treatment-covariate correlations before and after SBW"),
        p_sbw_raw_cprn = list(
            dt_results = quote(r_sbw$dt_results),
            width = 14, height = 7,
            caption = "covariate means by treatment before and after SBW"),
        p_ess = list(
            dt_results = quote(r_sbw$dt_results),
            width = 14, height = 7,
            caption = "effective sample size by treatment"),
        p_auc_euro_evst = list(
            dt_euro = quote(dt_euro),
            width = 14, height = 7,
            caption = "Event study plot of exhibition effects on auction turnover")
        
        
        ## p_grow_arrow = list(#  dont re-render for now
        ##     dt_insgrow5 = quote(dt_insgrow5),
        ##     width = 12, height = 6,
        ##     caption = "descs"),
        ## p_grow_perc_arrow = list(
        ##     dt_insgrow5 = quote(dt_insgrow5),
        ##     width = 12, height = 7,
        ##     caption = "grow"),
        ## p_insevc = list(
        ##     dt_insevc = quote(dt_insevc),
        ##     width = 14, height = 7,
        ##     caption = "EVC distribution")
        )}



gc_tbls <- function() {
    list(
        t_grow5_sumry  = list(
            dt_exhbppl_watyear = quote(dt_exhbppl_watyear),
            dt_insgrow5 = quote(dt_insgrow5),
            l_cols_summary = quote(l_cols_summary),
            caption = "summary stats"),
        t_insgrow5_topperc = list(
            dt_insgrow5 = quote(dt_insgrow5),
            dt_af_instns_rcd = quote(dt_af_instns_rcd),
            l_cols_summary = quote(l_cols_summary),
            caption = "most influential orgs (percentage)"),
        t_deg_info = list(
            g_insflow = quote(g_insflow),
            dt_af_instns_rcd = quote(dt_af_instns_rcd),
            caption = "degree info by form"),
        t_insgrow5_topcnt = list(
            dt_insgrow5 = quote(dt_insgrow5),
            dt_af_instns_rcd = quote(dt_af_instns_rcd),
            l_cols_summary = quote(l_cols_summary),
            caption = "most influential orgs (cnt)"),
        t_sumstats = list(
            r_sbw = quote(r_sbw),
            caption = "sumstats")
    )
}


gc_vvs <- function() {

    l_orgforms <- list(
        l_orgforms = c("pm", "toti",
                       gd_af_inst_clsfcn()[, setdiff(unique(get(c_datargs$clsfcn_to_use)), "yeet")]))



    l_vrbl_lbls <- list(
        InstitutionType         = "Institution Type",
        cnt_prevexhb_cum =   "# prev. Exhibitions",
        exn_oes = "Exhibitions (Other Exhibition Space)",
        exn_pi = "Exhibitions (Public Instiutiton)",
        exn_pm = "Exhibitions (Private Museum)",
        exn_toti = "Exhibitions (Top-tier museum)",
        career_year             = "Career Year",
        gender_rcd              = "Gender",
        gender_rcd_Male              = "Gender (Male)",
        gender_rcd_Other     = "Gender (Other)",
        atsubreg_145 = "Western Asia",
        atsubreg_151 = "Eastern Europe",
        atsubreg_154 = "Northern Europe",
        atsubreg_155 = "Western Europe",
        atsubreg_21  = "Northern America",
        atsubreg_30 = "Eastern Asia",
        atsubreg_39 = "Southern Europe",
        atsubreg_419 = "Latin America and the Caribbean",
        atsubreg_53 = "Australia and New Zealand",
        exn_all_cum_log   = "Prior Exhibitions (all types, log)",  
        exn_event_cum_log  = "Prior Exhibitions (Events,log)",
        exn_gallery_cum_log = "Prior Exhibitions (gallery, log)") %>%
        setNames(names(.), .)  # adjust order to allow easy combinations


    l_vrblgrp_mes <- list(
        num = c("cnt_prevexhb_cum", "exn_oes", "exn_pi", "exn_toti", "exn_pm", "career_year",
                "exn_all_cum_log","exn_event_cum_log","exn_gallery_cum_log"),
        bin = c("gender_rcd_Male", "gender_rcd_Other", "atsubreg_145", "atsubreg_151", "atsubreg_154",
                "atsubreg_155", "atsubreg_21", "atsubreg_30", "atsubreg_39", "atsubreg_419", "atsubreg_53"),
        cat = c("InstitutionType", "gender_rcd"))

    l_vrblinfo <- list(
        vrbl_lbl = l_vrbl_lbls,
        mes = l_vrblgrp_mes)

    dt_vrblinfo_prep <- lapply(seq_along(l_vrblinfo), \(i)
                               imap(l_vrblinfo[[i]], ~data.table(vrbl = .x, xxx = .y)) %>% rbindlist %>%
                               setnames(old = "xxx", new = names(l_vrblinfo)[i])) %>% # head(2) %>% 
        Reduce(\(x,y) join(x,y, on = "vrbl"), .) %>%
        .[, `:=`(vrblgrp = "all", vrblgrp_lbl = "all")]


    l_ctgterm_lbls <- list(
        ## list(vrbl = "InstitutionType", vlu = "nonpm", term_lbl = "Non-PM"),
        list(vrbl = "InstitutionType", vlu = "pm", term_lbl = "Private Museum"),
        list(vrbl = "InstitutionType", vlu = "gallery", term_lbl = "Gallery"),
        list(vrbl = "InstitutionType", vlu = "event", term_lbl = "Event"),
        list(vrbl = "InstitutionType", vlu = "npo", term_lbl = "NPO"),
        list(vrbl = "InstitutionType", vlu = "pes", term_lbl = "Public Exhibition Space"),
        list(vrbl = "InstitutionType", vlu = "pubm", term_lbl = "Public Museum"),
        list(vrbl = "InstitutionType", vlu = "oes", term_lbl = "Other Exhibition Space"),
        list(vrbl = "InstitutionType", vlu = "toti", term_lbl = "Top-Tier Museum"),
        list(vrbl = "InstitutionType", vlu = "pi", term_lbl = "Public Institution"))

    dt_ctgterm_lbls_prep <- rbindlist(l_ctgterm_lbls) %>%
        .[, `:=`(term = paste0(vrbl, vlu), vlu_lbl = term_lbl)] %>%
        dt_vrblinfo_prep[mes == "cat", .(vrbl, vrbl_lbl)][., on = "vrbl"] %>%
        .[, term_lbl := sprintf("%s - %s", vrbl_lbl, term_lbl)]


    ## need artificial merge column (tm = to merge) for cross join??
    dt_instypes <- dt_ctgterm_lbls_prep[vrbl == "InstitutionType", .(tm = 'a', vlu, term, vlu_lbl)]


    dt_ctgterm_interactions <- merge(
        dt_ctgterm_lbls_prep[vrbl %!in% c("InstitutionType", "natreg6", "locreg6")] %>% copy %>% .[, tm := 'a'],
        dt_instypes, allow.cartesian = T, by = 'tm') %>% .[, tm := NULL] %>%
        .[, `:=`(term = sprintf("%s:%s%s", term.y, vrbl, vlu.x),
                 term_lbl = sprintf("%s - %s * Ins. type (%s)", vrbl_lbl, vlu.x, vlu_lbl.y))] %>%
        .[, .(vrbl, vrbl_lbl, vlu = vlu.x, term_lbl, term, vlu_lbl = vlu_lbl.x)]

    dt_ctgterm_lbls <- rbind(dt_ctgterm_lbls_prep, dt_ctgterm_interactions)


    ## generate labels for interactions
    ## interact every variable with every institutiontype, just to have them all
    dt_interactions <- merge(dt_vrblinfo_prep[vrbl %!in% c("InstitutionType", "natreg6", "locreg6")] %>% copy %>%
                             .[, tm := 'a'],
                             dt_instypes, allow.cartesian = T, by = 'tm') %>% .[, tm := NULL] %>% 
        .[, `:=`(vrbl = sprintf("InstitutionType%s:%s", vlu, vrbl),
                 vrbl_lbl = sprintf("%s * Ins. type (%s)", vrbl_lbl, vlu_lbl))] %>%
        .[, .(vrbl, vrbl_lbl, mes, vrblgrp, vrblgrp_lbl)] # level, prom
    
    dt_vrblinfo <- rbind(dt_vrblinfo_prep, dt_interactions)

    dt_gof_cfg <- data.table(
        gof_name = c("nobs", "logLik", "AIC", "BIC", "pseudo.r.squared"),
        digits = c(0, 1, 1, 1, 2),
        gof_lbl = c("N", "logLik", "AIC", "BIC", "pseudo-R^2"))


    c_vvs <- c(## l_orgforms,
        l_orgforms,
        list(dt_vrblinfo = dt_vrblinfo),
        ## list(l_vrbltypes = l_vrbltypes), 
        list(dt_ctgterm_lbls = dt_ctgterm_lbls),
        list(dt_gof_cfg = dt_gof_cfg))

    return(c_vvs)

}


    
