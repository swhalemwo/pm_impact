rd1 <- coxph(Surv(career_year, career_year_end, last_exhb) ~ log1p(exhbN_nonpm_cum) + log1p(exhbN_pm_cum) +
                 log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum) + strata(year), dt_atyear_surv[year < 2015])

rc1 <- fepois(auc_turnover ~ log1p(exhbN_pm_cum) + log1p(exhbN_nonpm_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum) | PeopleID + year, dt_atyear, fixef.rm = "singleton")


rb1 <- fepois(exhbN_nonpm ~ log1p(exhbN_pm_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_nonpm_cum) + log1p(exhbN_event_cum) | PeopleID + year, dt_atyear, fixef.rm = "none")


gp_pminfl_nonpm <- function(rb1) {

    tidy(rb1) %>% adt %>%
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed")
    
}

gp_pminfl_nonpm(rb1)


gp_pminfl_auc <- function(rc1) {

    tidy(rc1) %>% adt %>%
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed")
}

gp_pminfl_auc(rc1)




gp_pminfl_surv <- function(rd1) {
    
    tidy(rd1) %>% adt %>%
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        .[, `:=`(estimate = estimate *-1)] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed")

}



gp_pminfl_surv(rd1)


gp_pminfl_cbn <- function(rb1, rc1, rd1) {

    imap(list("exhibitions (nonpm)" = rb1, "auction turnover" = rc1, survival = rd1),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        .[src == "survival", estimate := estimate *-1] %>%
        .[, src := factor(src, levels = c("exhibitions (nonpm)", "auction 1turnover", "survival"))] %>% 
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_grid(~src, scales = "free")

}
