
gp_hist_crap <- function(dt_atyear_foster) {
    gw_fargs(match.call())
    ## distribution of exhibitions by form: mostly gallery and pes
    ## fixme: ugly af, idk what's the point


    melt(dt_atyear_foster, id.vars = c("PeopleID", "year"), measure.vars = patterns("exhbN_.*_s5"),
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
    #' main-sequence inspired scatterplot: exhb_all_cum ~ career_year, with color of exhbN_all
    #' basically only shows cumulative advantage tho...
    #' more informative would be tiles
    

    ## descriptive plots
    ## exhbN_all_cum ~ cy: looks cute, but color is essentially also only cumulative advantage
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
        ggplot(aes(x=career_year, y = exhbN_all_cum, color = exhbN_all)) +
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

    ## exhbN_all ~ cy: idk what's the point
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>% 
        ggplot(aes(x=career_year, y = exhbN_all)) + geom_jitter(size = 0.1, height = 0.03) +
        scale_y_continuous(trans = scales::asinh_trans(), breaks = c(0, 3, 10, 30, 60))
    }


gp_cumad_jitter <- function(dt_atyear_foster) {
    gw_fargs(match.call())
    ## exhbN_all ~ exhnb_all_cum: muh cumulative advantage
    dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
        ggplot(aes(x=exhbN_all_cum, y = exhbN_all)) +
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
        ggplot(aes(y=exhbN_grow5_perc, x = exhbN_grow5_cnt, color = InstitutionType)) +
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
        .[, kde2d_tidy(asinh(exhbN_grow5_cnt), asinh(exhbN_grow5_perc),n=150), InstitutionType] %>%
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
        .[, exhbN_perc_n5s5_org := trimws(format(exhbN_perc_n5s5_org, digits = 2))] %>% # FIXME: absolute hacky way to get 
        .[exhbN_perc_n5s5_org == "NA", exhbN_perc_n5s5_org := ""] %>% 
        ## setnafill(type = "const", fill = 0, cols = l_cols_summary) %>%
        setnames(old = names(.), new = gsub("exhbN_", "", names(.))) %>% 
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

    dt_insgrow5[order(-exhbN_grow5_perc)][, head(.SD), InstitutionType,
                                          .SDcols = c("InstitutionID", l_cols_summary)] %>% 
        dt_af_instns_rcd[, .(InstitutionID, Name, City, Country)][., on = "InstitutionID"] %>%
        .[, `:=`(InstitutionID = NULL)] %>%
        ## .[, .(exhbN_grow5_cnt, exhbN_grow5_perc, evc_grow5_cnt, evc_grow5_prop)] %>% 
        setnafill(fill = 999, type = "const", cols = l_cols_summary) %>%
        gt_plain
    
}

gt_insgrow5_topcnt <- function(dt_insgrow5, dt_af_instns_rcd, l_cols_summary) {
    gw_fargs(match.call())

    ## most influential in exhibition count
    dt_insgrow5[order(-exhbN_grow5_cnt)][, head(.SD), InstitutionType] %>% 
        copy(dt_af_instns_rcd)[, InstitutionType := NULL][., on = "InstitutionID"] %>% 
        .[, `:=`(InstitutionID = NULL)] %>%
        ## .[, .(exhbN_grow5_cnt, exhbN_grow5_perc, evc_grow5_cnt, evc_grow5_prop)] %>% 
        setnafill(fill = 999, type = "const", cols = l_cols_summary) %>%
        gt_plain

}


gd_nbrs <- function() {

    dt_ynkplt <- gc_ynkplt()
    dt_refplt <- gc_refplt()

    ## dt_ynktbl <- gc_ynktbl()
    dt_reftbl <- gc_reftbl()


    dt_nbrs_cbnd <- rbindlist(list(dt_ynkplt, dt_refplt, dt_reftbl))
    return(dt_nbrs_cbnd)

}
