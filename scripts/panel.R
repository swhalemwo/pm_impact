
## * panel models




gl_mdls <- function(dt_atyear_foster, dt_atyear_ec) {
    gw_fargs(match.call())

    ## outcome auction turnover
    r_auc <- fepois(auc_turnover ~ log1p(exhbN_pm_cum) + log1p(exhbN_pes_cum) + log1p(exhbN_gallery_cum) +
                        log1p(exhbN_event_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "singleton",
                    cluster = "PeopleID")

    ## outcome exhb in PES
    r_pes <- fepois(exhbN_pes ~ log1p(exhbN_pm_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_pes_cum) +
                        log1p(exhbN_event_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none",
                    cluster = "PeopleID")

    ## ## outcome PES with interaction
    ## r_pes2 <- fepois(exhbN_pes ~ log1p(exhbN_gallery_cum) + log1p(exhbN_pes_cum)*log1p(exhbN_pm_cum) +
    ##                      log1p(exhbN_event_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

    ## r_pres3 <- fepois(exhbN_pes ~ log1p(exhbN_pes_cum)*log1p(exhbN_pm_cum) +
    ##                       log1p(exhbN_pm_cum)*log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum)*
    ##                       log1p(exhbN_pm_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

    ## evc DV
    r_evc <- feols(evc_ay_ttl ~ log1p(exhbN_pm_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_pes_cum) +
                       log1p(exhbN_event_cum) | PeopleID + year, dt_atyear_evc, fixef.rm = "none")
    


    l_mdls <- list(r_auc = r_auc,
                   r_pes = r_pes,
                   r_evc = r_evc)

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




gp_pminfl_cbn <- function(r_pes, r_auc) {

    imap(list("exhibitions (pes)" = r_pes, "auction turnover" = r_auc),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        ## .[src == "survival", estimate := estimate *-1] %>%
        .[, src := factor(src, levels = c("exhibitions (pes)", "auction turnover"))] %>% 
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbar(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_grid(~src, scales = "free")

}

## gp_pminfl_cbn(r_pes, r_auc)











