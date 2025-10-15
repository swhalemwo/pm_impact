
## * panel v0,
## idk when they were added.. 
fepois(exhbN_all ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID, dt_atyear)

r_auc_ols <- feols(auc_turnover ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID,
                   dt_atyear, fixef.rm = "singleton")

r_auc_poi <- fepois(auc_turnover ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID,
                    dt_atyear, fixef.rm = "singleton")

screenreg(list(r_auc_ols, r_auc_poi))


## * panel models

## outcome survival: not cute -> probably yeet
rd1 <- coxph(Surv(career_year, career_year_end, last_exhb) ~ log1p(exhbN_pes_cum) + log1p(exhbN_pm_cum) +
                 log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum) + strata(year), dt_atyear_surv[year < 2015])

## outcome auction turnover
rc1 <- fepois(auc_turnover ~ log1p(exhbN_pm_cum) + log1p(exhbN_pes_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum) | PeopleID + year, dt_atyear, fixef.rm = "singleton")

# outcome exhb in PES
rb1 <- fepois(exhbN_pes ~ log1p(exhbN_pm_cum) + log1p(exhbN_gallery_cum) + log1p(exhbN_pes_cum) + log1p(exhbN_event_cum) | PeopleID + year, dt_atyear, fixef.rm = "none")

# outcome PES with interaction
rb2 <- fepois(exhbN_pes ~ log1p(exhbN_gallery_cum) + log1p(exhbN_pes_cum)*log1p(exhbN_pm_cum) + log1p(exhbN_event_cum) | PeopleID + year, dt_atyear, fixef.rm = "none")

rb3 <- fepois(exhbN_pes ~ log1p(exhbN_pes_cum)*log1p(exhbN_pm_cum) + log1p(exhbN_pm_cum)*log1p(exhbN_gallery_cum) + log1p(exhbN_event_cum)*log1p(exhbN_pm_cum) | PeopleID + year, dt_atyear, fixef.rm = "none")

screenreg(list(rb1,rb2,rb3))




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

    imap(list("exhibitions (pes)" = rb1, "auction turnover" = rc1, "survival" = rd1),
         ~adt(tidy(.x))[, src := .y]) %>% rbindlist %>%
        .[src == "survival", estimate := estimate *-1] %>%
        .[, src := factor(src, levels = c("exhibitions (pes)", "auction turnover", "survival"))] %>% 
        .[, term := str_extract_all(term, "(?<=exhbN_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_grid(~src, scales = "free")

}

gp_pminfl_cbn(rb1, rc1, rd1)




## experiment with splines
rb4 <- fepois(exhbN_pes ~ log1p(exhbN_pes_cum)*log1p(exhbN_pm_cum)*bs(exhbN_pes_cum) | PeopleID + year, dt_atyear, fixef.rm = "none")

## do i need this
rs1 <- fepois(exhbN_pes ~ bs(log1p(exhbN_all_cum), df = 10)*log1p(exhbN_all_s5) | PeopleID + year,
              dt_atyear, fixef.rm = "none")

## or this
rs2 <- fepois(exhbN_pes ~ I(bs(log1p(exhbN_all_cum))^2) | PeopleID + year, dt_atyear, fixef.rm = "none")

## or this
rs3 <- fepois(exhbN_pes ~ log1p(exhbN_all_cum)*bs(log1p(exhbN_all_cum)) | PeopleID + year, dt_atyear, fixef.rm = "none")


dt_topred <- expand.grid(exhbN_all_cum = seq(1, 500, length.out = 30),
                         exhbN_all_s5 = seq(1, 15),
                         PeopleID = 11, year = 2015) %>% adt


## density of s5
ggplot(dt_atyear, aes(x = exhbN_all_s5)) + geom_density() +
    scale_x_continuous(transform = scales::log10_trans())


## some hacky way of getting at marginal effect by exhbN_cum 

dt_topred[, `:=`(pred1 = predict(rs3, dt_topred),
                 pred2 = predict(rs3, copy(dt_topred)[, exhbN_all_cum := exhbN_all_cum +1]),
                 pred3 = predict(rs1, dt_topred),
                 pred4 = predict(rs1, copy(dt_topred)[, exhbN_all_s5 := exhbN_all_s5 +1]))] %>%
    .[, `:=`(eff1 = pred2/pred1, eff2 = pred3/pred4)]

## rs3 figure
ggplot(dt_topred, aes(x=exhbN_all_cum, y = eff1)) + geom_point()

## rs1 figure
ggplot(dt_topred, aes(x=exhbN_all_cum, y = eff2, color = exhbN_all_s5)) + geom_point()




