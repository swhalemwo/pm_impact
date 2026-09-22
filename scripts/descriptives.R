
gp_hist_crap <- function(dt_atyear_foster) {
    gw_fargs(match.call())
    ## distribution of exhibitions by form: mostly gallery and pes
    ## fixme: ugly af, idk what's the point


    melt(dt_atyear_foster, id.vars = c("PeopleID", "year"), measure.vars = patterns("exn_.*_s5"),
         variable.name = "vrbl", value.name = "vlu") %>%
        .[grepl("s5", vrbl)] %>% na.omit %>%
        ggplot(aes(x=vlu, fill = vrbl)) +
        ## geom_density(bw = 0.2) +
        geom_histogram(alpha = 0.5) + 
        facet_grid(vrbl~., switch = "y", scales = "free") +
        theme(strip.text.y.left = element_text(angle = 0)) +
        scale_x_continuous(trans = scales::transform_asinh(), breaks = c(0,3,10,30,100,300,600))

}

gp_mainseq <- function(dt_atyear_foster) {
    gw_fargs(match.call())
    #' main-sequence inspired scatterplot: exhb_all_cum ~ career_year, with color of exn_all
    #' basically only shows cumulative advantage tho...
    #' more informative would be tiles
    

    ## descriptive plots
    ## exn_all_cum ~ cy: looks cute, but color is essentially also only cumulative advantage
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
        ggplot(aes(x=career_year, y = exn_all_cum, color = exn_all)) +
        geom_jitter(size = 0.4, height = 0.01) +
        ## scale_y_continuous(trans = scales::log10_trans()) +
        scale_y_continuous(trans = scales::asinh_trans(), breaks = c(0,3, 10, 30,100,300,750)) +
        scale_color_gradient(low = "red", high = "yellow",
                             trans = scales::asinh_trans(), breaks = c(0,3, 10,30,80)) +
        theme(panel.background = element_rect(fill = "black"),
              panel.grid.major = element_line(color = "black"),
              panel.grid.minor = element_line(color = "black")
              )

}


gp_activity_jitter <- function(dt_atyear_foster) {
    gw_fargs(match.call())

    ## exn_all ~ cy: idk what's the point
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>% 
        ggplot(aes(x=career_year, y = exn_all)) + geom_jitter(size = 0.1, height = 0.03) +
        scale_y_continuous(trans = scales::asinh_trans(), breaks = c(0, 3, 10, 30, 60))
    }


gp_cumad_jitter <- function(dt_atyear_foster) {
    gw_fargs(match.call())
    ## exn_all ~ exhnb_all_cum: muh cumulative advantage
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
        ggplot(aes(x=exn_all_cum, y = exn_all)) +
        geom_jitter(size = 0.1, # width = 0.003, if x not log
                    ## width = 0.03, # if x log
                    width = 0.05,
                    height = 0.2) +
        scale_y_continuous(trans = scales::asinh_trans(), breaks = c(0,3,10,30,60)) +
        scale_x_continuous(trans = scales::asinh_trans(), breaks = c(0,3,10,30,100,300,1000))
}


gp_insgrow5 <- function(dt_insgrow5) {
    ## growth count by growth perc: basically shows that PM grow less percent for same count -> orientation towards established



    dt_insgrow5 %>% .[, .SD[sample(1:.N, size = min(.N, 500))], InstitutionType] %>% 
        ggplot(aes(y=exn_grow5_perc, x = exn_grow5_cnt, color = InstitutionType)) +
        geom_point(size = 0.4) + 
        ## geom_smooth(se = F) +
        scale_x_continuous(trans = scales::transform_asinh(), breaks = l_asinh_breaks) +
        scale_y_continuous(trans = scales::transform_asinh(), breaks = l_asinh_breaks) +
        theme(legend.position = "bottom")
    

}


gp_insgrow5_contours <- function(dt_insgrow5) {
    gw_fargs(match.call())
    ## growth scatterplot with contours

    p_insgrow5 <- gp_insgrow5(dt_insgrow5)

    dt_contours <- dt_insgrow5 %>% na.omit %>%
        .[, kde2d_tidy(asinh(exn_grow5_cnt), asinh(exn_grow5_perc),n=150), InstitutionType] %>%
        .[, contourLines_tidy(x,y,z, levels = gd_zlvls(z, c(0.1, 0.2, 0.5, 0.8))[, z]), InstitutionType] %>% 
        .[, level2 := as.integer(factor(level)), InstitutionType] %>% #  .[, .N, level2]
        .[level2 == 4, .(x = sinh(x), y= sinh(y), InstitutionType)]

    p_insgrow5 + 
        geom_path(dt_contours, mapping = aes(x=x, y=y, color = InstitutionType), linewidth = 2) 
}




contourLines_tidy <- function(x,y,z, ...)  {
    ## hacky way of getting contourlines from tidy input

    ## first need to transform z into matrix, this is slow af
    m_z <- data.table(x=x, y=y, z=z) %>% 
        .[, .(xfac = as.integer(factor(x, levels = sort(unique(x)))),
              yfac = as.integer(factor(y, levels = sort(unique(y)))), z)] %>%
        dcast(xfac~yfac, value.var = "z") %>% .[, xfac := NULL] %>% as.matrix

    x_unq <- unique(x) %>% sort
    y_unq <- unique(y) %>% sort
    
    
    contourLines(x_unq,y_unq,m_z, ...) %>% rbindlist
}

gd_zlvls <- function(z, l_lvls) {
    ## get thresholds to pass to contours

    data.table(z=z) %>% .[order(z)] %>%
        .[, zprop := fcumsum(z)/sum(z)] %>% # tail %>% 
        .[J(zprop = l_lvls), on = "zprop", roll= T]
}
        


smooth_it <- function(x, y, n = 1000, method = "natural") {
    #' hacky version to get splines, but not sure if necessary (not used atm)
    t <- seq_along(x)
    new_t <- seq(min(t), max(t), length.out = n)
    new_x <- spline(t, x, xout = new_t, method = method)$y
    new_y <- spline(t, y, xout = new_t, method = method)$y
    data.frame(t = new_t, x = new_x, y = new_y)
}


gt_grow5_sumry <- function(dt_exhbppl_watyear, dt_insgrow5, l_cols_summary) {
    #' summary stats of 5-year growth by form
    gw_fargs(match.call())

    rbind(
        dt_exhbppl_watyear[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)), InstitutionType,
                           .SDcols = intersect(l_cols_summary, names(dt_exhbppl_watyear))] %>%
                               .[, src := "exhblink"] %>% melt(id.vars = c("InstitutionType", "src")), 
        dt_insgrow5[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)), InstitutionType,
                    .SDcols = intersect(names(dt_insgrow5), l_cols_summary)] %>% .[, src := "ins"] %>%
        melt(id.vars = c("InstitutionType", "src"))) %>%
        dcast(InstitutionType + src ~ variable) %>% 
        .[order(InstitutionType), .SD, .SDcols = c("InstitutionType", "src", l_cols_summary)] %>%
        .[, exn_perc_n5s5_org := trimws(format(exn_perc_n5s5_org, digits = 2))] %>% # FIXME: absolute hacky way to get 
        .[exn_perc_n5s5_org == "NA", exn_perc_n5s5_org := ""] %>% 
        ## setnafill(type = "const", fill = 0, cols = l_cols_summary) %>%
        setnames(old = names(.), new = gsub("exn_", "", names(.))) %>% 
        gt_plain
}


gp_insgrow5_dist <- function(dt_insgrow5, l_cols_summary) {
    gw_fargs(match.call())

    ## distribution of percentage count/growth of institutions by institution type
    melt(dt_insgrow5, id.vars = c("InstitutionID", "InstitutionType"), measure.vars = l_cols_summary,         
         variable.name = "mes", value.name = "vlu") %>%
        ggplot(aes(x=vlu, color = InstitutionType)) +
        geom_density()+
        facet_wrap(~mes, ncol =1, scale = "free") +
        scale_x_continuous(trans = scales::asinh_trans(), breaks = l_asinh_breaks)

}
    



gt_insgrow5_topperc <- function(dt_insgrow5, dt_af_instns_rcd, l_cols_summary) {
    gw_fargs(match.call())
    #' most influential orgs by type in terms of percent

    dt_insgrow5[order(-exn_grow5_perc)][, head(.SD), InstitutionType,
                                          .SDcols = c("InstitutionID", l_cols_summary)] %>% 
        dt_af_instns_rcd[, .(InstitutionID, Name, City, Country)][., on = "InstitutionID"] %>%
        .[, `:=`(InstitutionID = NULL)] %>%
        ## .[, .(exn_grow5_cnt, exn_grow5_perc, evc_grow5_cnt, evc_grow5_prop)] %>% 
        setnafill(fill = 999, type = "const", cols = l_cols_summary) %>%
        gt_plain
    
}

gt_insgrow5_topcnt <- function(dt_insgrow5, dt_af_instns_rcd, l_cols_summary) {
    gw_fargs(match.call())

    ## most influential in exhibition count
    dt_insgrow5[order(-exn_grow5_cnt)][, head(.SD), InstitutionType] %>% 
        copy(dt_af_instns_rcd)[, InstitutionType := NULL][., on = "InstitutionID"] %>% 
        .[, `:=`(InstitutionID = NULL)] %>%
        ## .[, .(exn_grow5_cnt, exn_grow5_perc, evc_grow5_cnt, evc_grow5_prop)] %>% 
        setnafill(fill = 999, type = "const", cols = l_cols_summary) %>%
        gt_plain

}


gd_nbrs <- function() {

  l_af_desc_prep <- list(
    list(nbr_name = "N_af_people", nbr_fmt = dt_af_people[, format(.N, big.mark = " ")]),
    list(nbr_name = "N_af_instns", nbr_fmt = dt_af_instns[, format(.N, big.mark = " ")]),
    list(nbr_name = "N_af_exhbs", nbr_fmt = dt_af_exhbs[, format(.N, big.mark = " ")]),
    list(nbr_name = "N_af_instntype", nbr_fmt = dt_af_instns[, uniqueN(InstitutionType)]),
    list(nbr_name = "N_af_exhblinks", nbr_fmt = dt_af_links_exhb_ppl[, format(.N, big.mark = " ")]),
    list(nbr_name = "N_pm_af_instn_matches", nbr_fmt = dt_af_pmdb_matches[AF_IID != "nomatch", .N]),
    list(nbr_name = "perc_pm_af_instn_matches",
         nbr_fmt = format(100*dt_af_pmdb_matches[AF_IID != "nomatch", .N]/dt_pmdb[, .N], digits = 3)),
    list(nbr_name = "perc_af_top10_share",
         nbr_fmt = dt_af_links_exhb_ppl[, .N, PeopleID][order(-N)][, rank := 1:.N] %>% 
           .[, top_10perc := rank < .N * 0.1] %>% .[, .(exhbcnt = sum(N)), top_10perc] %>%
           .[, ttl_cnt := sum(exhbcnt)] %>% .[top_10perc == T, format(100*exhbcnt/ttl_cnt, digits = 2)]
         ))

  dt_af_desc <- rbindlist(l_af_desc_prep)[, grp := "af_desc"]

  ## organization-form exhibition info
  dt_ofexbh <- merge(dt_af_exhbs, dt_af_instns_rcd, by = "InstitutionID") %>%
    .[, .N, InstitutionType] %>% .[, `:=`(N= N*1.0, perc =100*N/sum(N))] %>%
    melt(id.vars = "InstitutionType", variable.name = "mes",) %>%
    .[, .(nbr_name = sprintf("ofexhb_%s_%s", mes, InstitutionType),
          nbr_fmt= format(value, digits = 2), grp = "af_ofexhb_stats"), .I] %>% .[, I := NULL]

  dt_tanp <- gd_tanp_muyr()[order(year, -total)][, head(.SD, 100), year]

  ## instittution type counts
  dt_instype_info <- dt_af_instns_rcd[, .N, InstitutionType] %>%
    .[, .(nbr_name = sprintf("instype_cnt_%s", InstitutionType),
          nbr_fmt = format(N, big.mark = " "), grp = "instype_cnts"), .I] %>% .[, I := NULL]
  
  ## dt_tanp[, .N, .(muci_id, muci_name, city_name)][!dt_tanp_muyr_afid[, .N, muci_id], on = "muci_id"][order(-N)]
  ## tanp coverage
  dt_tanp_info <- tribble(
    ~nbr_name, ~nbr_fmt, 
    "tanp_Nins_ttl", dt_tanp[, uniqueN(muci_id)],
    "tanp_Nins_wafid", dt_tanp_muyr_afid[, uniqueN(muci_id)],
    "tanp_Nins_woafid", dt_tanp[, uniqueN(muci_id)]- dt_tanp_muyr_afid[, uniqueN(muci_id)]) %>% adt %>%
    .[, grp := "tanp"]

  ## dt_exhbppl numbers
  dt_datasize <- merge(
    dt_exhbppl[, .(N_link = .N, N_instn = uniqueN(InstitutionID),
                   N_at = uniqueN(PeopleID), N_exhb = uniqueN(ExhibitionID))] %>%
    melt(measure.vars = names(.), variable.factor = F, value.name = "vlu_exhbppl"),
    dt_exhbppl_wyeet[, .(N_link = .N, N_instn = uniqueN(InstitutionID),
                         N_at = uniqueN(PeopleID), N_exhb = uniqueN(ExhibitionID))] %>%
    melt(measure.vars = names(.), variable.factor = F, value.name = "vlu_wyeet"),
    by = "variable") %>% setnames(old = "variable", new = "vrbl") %>% 
    .[, perc_of_wyeet := 100*(vlu_exhbppl/vlu_wyeet)] %>% 
    melt(id.vars = "vrbl", measure.vars = c("vlu_exhbppl", "perc_of_wyeet")) %>%
    .[, .(nbr_name = sprintf("%s_%s", vrbl, variable),
          nbr_fmt = format(value, big.mark = " ", scientific = F,  trim = T, digits= 3)), variable] %>%
    .[, `:=`(variable  = NULL, grp = "datasize")]

  dt_ynkplt <- gc_ynkplt()
  dt_refplt <- gc_refplt()

  ## dt_ynktbl <- gc_ynktbl()
  dt_reftbl <- gc_reftbl()


  dt_nbrs_cbnd <- rbindlist(list(dt_ynkplt, dt_refplt, dt_reftbl, dt_af_desc, dt_ofexbh, dt_tanp_info,
                                 dt_instype_info, dt_datasize))
  return(dt_nbrs_cbnd)

}


gp_dist_exn_pi <- function(dt_atyear_foster) {

## dt_atyear_foster[, .N, exn_pes] %>%
##     .[, perc := 100*(N/sum(N))] %>% 
##     ggplot(aes(x=exn_pes, y = N)) +
##     ## geom_col(width = 0.3) +
##     geom_step() + 
##     scale_x_continuous(trans = scales::asinh_trans(), breaks = l_asinh_breaks)



    dt_atyear_foster[, .N, cut(exn_pi, unique(l_asinh_breaks))] %>%
        .[, perc := 100*(N/sum(N))] %>% 
        .[, perc_fmt := sprintf("%s%%", format(perc,digits = 1, nsmall = 0)), .I] %>% 
        .[, c("edge_lo", "edge_hi") := jtls::edge(cut), .I] %>%
        .[, edge_lop1 := edge_lo + 1] %>%
        .[, cut_exhb_pi := ifelse(edge_hi == edge_lop1, edge_hi, sprintf("%s-%s", edge_lop1, edge_hi))] %>%
        .[order(cut)] %>% .[, cut_exhb_pi := factor(cut_exhb_pi, cut_exhb_pi)] %>% 
        ggplot(aes(x=cut_exhb_pi, y = N, label = perc_fmt)) +
        geom_col() +
        geom_text(vjust = -0.7) +
        labs(x="Number of exhibitions in Public Institutions per year") +
        scale_y_continuous(labels = scales::label_comma()) +
        coord_cartesian(ylim = c(NA, 5.5e5))

}


gp_dist_aucturn <- function(dt_atyear_foster) {

    dt_atyear_foster %>%
        .[auc_turnover > 0] %>% 
        ggplot(aes(x=auc_turnover)) + geom_density() +
        scale_x_continuous(trans = scales::log10_trans(), breaks = c(0, 10^(c(4,8,seq(1,8,2)))),
                           labels = label_number(scale_cut = scales::cut_short_scale())) +
        labs(x = "Yearly Auction turnover (>0)")
}




gp_grow_arrow <- function(dt_insgrow5) {
    gw_fargs(match.call())
    
    dt_sums <- dt_insgrow5[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)), .(vlu = InstitutionType),
                           .SDcols = c("exn_all_s5", "exn_grow5_cnt", "exn_perc_n5s5_org")] %>%
        .[, `:=`(x = "exn_grow5_cnt", y2 = exn_grow5_cnt)] %>%
        .[vlu == "pes", y2 := y2 + 3]

    dt_sums2 <- merge(dt_sums, gc_vvs() %>% chuck("dt_ctgterm_lbls"), by = "vlu")
    
    c_lbls_x <- c(exn_all_s5 = "Before (5 year span)", exn_grow5_cnt = "After (5 year span)")


    melt(dt_sums, id.vars = "vlu", measure.vars = c("exn_all_s5", "exn_grow5_cnt")) %>%
        ggplot(aes(x=variable, y = value, color = vlu, group = vlu)) +
        geom_point(size = 2, show.legend = F) +
        geom_path(arrow = arrow(angle = 15, type = "closed"), show.legend = F) +
        geom_text_repel(mapping = aes(label = format(value, digits =3)), direction = "y", show.legend= F) + 
        geom_label(dt_sums2, mapping = aes(x=x, y = y2, label = vlu_lbl), 
                   position = position_nudge(x=0.12), hjust = 0, show.legend = F) +
        scale_x_discrete(labels = c_lbls_x,
                         ## expand= c(0.1, 0)) +
                         expand = expansion(add = c(0.4,1.2))) + 
        labs(x = NULL, y = "Number of exhibitions") +
        scale_y_continuous(limits = c(15,NA))
        ## expand_limits(x=c(3,3))

}


