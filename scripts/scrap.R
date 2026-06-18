nlist <- function(...) {
    vx <- substitute(c(...))# [-1L]
    print(vx)
    print(len(vx))
    return(setNames(vx,vx))
}

nlist("hi", "there")
nlist(hi, there)
print(len(c("hi", "there")))

setNames(c("hi", "there"), c("hi", "there"))


dt_atyear_evc[, c(nlist(PeopleID, year), .SD) , .SDcols = patterns("^evc")]


## reversed graph: this would basically be who shows the most established art?
## kinda going back in time? 
g_insflow2 <- dt_edges_wtd[, .(InstitutionID = InstitutionID2, InstitutionID2 = InstitutionID, N)] %>%
    gg_insflow()
     

gd_deg_info(g_insflow) %>% merge(dt_af_instns, by = "InstitutionID") %>%
    .[order(-deg_in)]


## dt_ppl_evcsum <- gd_ppl_evcsum(dt_exhbppl_evc)

## dt_atofyear <- gd_atofyear(dt_atyear_evc)


## ** looking at yuge EVC grow values

dt_insgrow5[evc_grow5_prop > 300000]
dt_insgrow5[!is.infinite(evc_grow5_prop), mean(evc_grow5_prop, na.rm = T), InstitutionType]


## 3 pes                     12958           15.9             164.       -0.0343         7.33e19

dt_exhbppl_watyear[InstitutionID == 12958, .(begin_year, PeopleID, ExhibitionID, evc_s5, evc_n5, evc_grow5_prop)] %>% print(n=80) %>% .[, mean(evc_grow5_prop, na.rm = T)]

merge(
    dt_exhbppl_watyear[PeopleID == 210084, .(begin_year, InstitutionID)],
    dt_insevc[, .(InstitutionID, Name, InstitutionType, evc)], by = "InstitutionID")
## ok so all of 210084's exhibitions for some years were at some shit org %>% average evc super low
## hmm maybe the variation in evc is also so yuge because of sampling....

## hmm if you one artist who makes it somewhere from nowhere, that's a yuuuuuuuuuuge boost
## maybe normalize EVC to some orders of magnitude?


## * experiment with splines
rb4 <- fepois(exn_pes ~ log1p(exn_pes_cum)*log1p(exn_pm_cum)*bs(exn_pes_cum) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

## do i need this
rs1 <- fepois(exn_pes ~ bs(log1p(exn_all_cum), df = 10)*log1p(exn_all_s5) | PeopleID + year,
              dt_atyear_foster, fixef.rm = "none")

## or this
rs2 <- fepois(exn_pes ~ I(bs(log1p(exn_all_cum))^2) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")

## or this
rs3 <- fepois(exn_pes ~ log1p(exn_all_cum)*bs(log1p(exn_all_cum)) | PeopleID + year, dt_atyear_foster, fixef.rm = "none")


dt_topred <- expand.grid(exn_all_cum = seq(1, 500, length.out = 30),
                         exn_all_s5 = seq(1, 15),
                         PeopleID = 11, year = 2015) %>% adt


## density of s5
ggplot(dt_atyear_foster, aes(x = exn_all_s5)) + geom_density() +
    scale_x_continuous(transform = scales::log10_trans())


## some hacky way of getting at marginal effect by exn_cum 

dt_topred[, `:=`(pred1 = predict(rs3, dt_topred),
                 pred2 = predict(rs3, copy(dt_topred)[, exn_all_cum := exn_all_cum +1]),
                 pred3 = predict(rs1, dt_topred),
                 pred4 = predict(rs1, copy(dt_topred)[, exn_all_s5 := exn_all_s5 +1]))] %>%
    .[, `:=`(eff1 = pred2/pred1, eff2 = pred3/pred4)]

## rs3 figure
ggplot(dt_topred, aes(x=exn_all_cum, y = eff1)) + geom_point()

## rs1 figure
ggplot(dt_topred, aes(x=exn_all_cum, y = eff2, color = exn_all_s5)) + geom_point()



## * survival model -> yeet
## outcome survival: not cute -> probably yeet
rd1 <- coxph(Surv(career_year, career_year_end, last_exhb) ~ log1p(exn_pes_cum) + log1p(exn_pm_cum) +
                 log1p(exn_gallery_cum) + log1p(exn_event_cum) + strata(year), dt_atyear_surv[year < 2015])


gp_pminfl_surv <- function(rd1) {
    
    tidy(rd1) %>% adt %>%
        .[, term := str_extract_all(term, "(?<=exn_)[a-zA-Z]+(?=_cum)")[[1]], .I] %>%
        .[, `:=`(estimate = estimate *-1)] %>% 
        ggplot(aes(x = estimate, y = term, xmin = estimate - 1.96 *std.error, xmax = estimate + 1.96 *std.error,
                   color = term)) +
        geom_point(show.legend = F) +
        geom_errorbarh(show.legend = F) + 
        geom_vline(xintercept = 0, linetype = "dashed")

}



gp_pminfl_surv(rd1)




## * look into tweedie

## library(statmod)

re1 <- feols(evc_ay_ttl ~ log1p(exn_pm_cum) + log1p(exn_gallery_cum) + log1p(exn_pes_cum) + log1p(exn_event_cum) | PeopleID + year, dt_atyear_evc, fixef.rm = "none")


re2 <- feglm(evc_ay_ttl ~ log1p(exn_pm_cum) + log1p(exn_gallery_cum) + log1p(exn_pes_cum) + log1p(exn_event_cum) | PeopleID + year, dt_atyear_evc, fixef.rm = "none", family = tweedie())

re3 <- feglm(evc_ay_ttl ~ log1p(exn_pm_cum) + log1p(exn_gallery_cum) + log1p(exn_pes_cum) + log1p(exn_event_cum) | PeopleID + year, dt_atyear_evc, fixef.rm = "none", family = tweedie(var.power = 0.5))

re4 <- feglm(evc_ay_ttl ~ log1p(exn_pm_cum) + log1p(exn_gallery_cum) + log1p(exn_pes_cum) + log1p(exn_event_cum) | PeopleID + year, dt_atyear_evc, fixef.rm = "none", family = tw(link = "log"))



map(seq(0, 5, 0.5), ~data.table(x = rpois(n=100, lambda = .x)) %>% .[, lambda := .x]) %>% rbindlist %>%
    .[, .N, .(x, lambda)] %>%
    ggplot(aes(x=x, y=N, color = factor(lambda))) + geom_line()



## tweedie: used for logged DVs which originally include 0
## https://library.virginia.edu/data/articles/getting-started-tweedie-models-0
## why even log DV? if you don't log DV you can use OLS just fine

library(tweedie)
rtweedie

dt_tweedie_params <- expand.grid(phi = seq(1, 2,0.5), mu = seq(0.5, 2.5), power = seq(1.05, 1.95, 0.2)) %>%
    adt %>% .[, `:=`(n=1000, id = paste0("V", 1:.N))]

apply(dt_tweedie_params[, .(phi, mu, power, n)], 1, \(x) do.call(rtweedie, as.list(x))) %>% adt %>%
    melt(measure.vars = names(.), variable.name = "id") %>%
    merge(dt_tweedie_params) %>%
    ggplot(aes(x=value, color = factor(phi))) +
    geom_density() +
    facet_wrap(power~mu, scales = "free", nrow = 5, ncol = 3)


## * debug differences between cnt_prevexhb and exn_all_cum

dt_exhbppl_watyear_chk <- copy(dt_exhbppl_watyear) %>%
    .[, .(PeopleID, begin_year, cnt_prevexhb, exn_all_cum, diff_cum = cnt_prevexhb - exn_all_cum,
          InstitutionType)]

## dt_exhbppl_watyear_chk[, .N, InstitutionType]

dt_exhbppl_watyear_chk[, .N, diff_cum] %>% .[order(diff_cum)] %>% print(n=80)


dt_exhbppl_watyear_chk[diff_cum < 0]

dt_exhbppl_watyear_chk[PeopleID == 492 & diff_cum < 0]

dt_exhbppl[PeopleID == 492, .(PeopleID, begin_year, cnt_prevexhb)] %>% print(n=85)

dt_atyear_foster[PeopleID == 492, .SD, .SDcols = patterns("PeopleID|^year$|cum$")] %>% print(n=80)

dt_atyear_foster %>% copy %>%
.[PeopleID == 492, .(year, exn_all_cum, exn_all, exn_all_cum2 = fcumsum(exn_all))] %>%
.[, xx := exn_all_cum2 - exn_all] %>% print(n=80)



dt_atyear_raw[PeopleID == 492, .(PeopleID, year, exn_all, xx = fcumsum(exn_all))]
dt_exhbppl[PeopleID == 492 & begin_year < 2006]
dt_exhbppl_wyeet[PeopleID == 492 & begin_year < 2006]

dt_exhbppl_wyeet[PeopleID == 492 & begin_year < 2006, .N, .(PeopleID, begin_year)] %>%
.[, cumx := fcumsum(N)]


## diff between exhbppl and atyear

dt_atyear_rcd_wide[PeopleID == 492] %>% copy %>%
.[, cnt_all := rowSums(.SD, na.rm = T), .SDcols = patterns("exn")]
## here it's still 10 for 2005

## huh here it's 11
dt_cbn_prep1[PeopleID == 492]

dt_atyear_exhbcnt[PeopleID == 492] # here i's 11

dt_cbn[PeopleID == 492 & begin_year == 2005, uniqueN(ExhibitionID)] ## here it's 11

## check out gd_atyear_rcd
dt_atyear_rcd[PeopleID == 492 & begin_year == 2005]

## reproduce manually
dt_m1 <- merge(dt_af_links_exhb_ppl, dt_af_exhbs[, .(ExhibitionID, InstitutionID, begin_year)],
               by = "ExhibitionID")

dt_m2 <- merge(dt_m1, dt_af_instns, by = "InstitutionID", all.x = T)

dt_m2[is.na(InstitutionType)]

## ** more model experimentation

## EE: exhibition effects
## most basic idea: mean by institution type
r_exef1 <- feols(exn_all_n5 ~ InstitutionType, dt_exhbppl_watyear)
dt_exhbppl_watyear[, mean(exn_all_n5, na.rm = T), InstitutionType] # equivalent descriptive

## adding artist FE: now within-artist variation: much much smaller
r_exef2 <- feols(exn_all_n5 ~ InstitutionType | PeopleID, dt_exhbppl_watyear, fixef.rm = "singleton")

## 
r_exef3 <- feols(exn_all_n5 ~ InstitutionType + cnt_prevexhb | PeopleID, dt_exhbppl_watyear, fixef.rm = "singleton")

## testing that with panel data
## not sure if N5 is valid as DV: i feel i'm smoothing/re-using same data multiple times

dt_atyear_foster[PeopleID==44, .(PeopleID, year, exn_all_n5, exn_pm, exn_gallery, exn_pes, exn_all_cum)]

r_exef4 <- fepois(exn_all_n5 ~ exn_pm + exn_gallery + exn_pes + exn_event +
                      log1p(exn_all_cum) | PeopleID, dt_atyear_foster)

r_exef5 <- fepois(exn_all_n5 ~ exn_pm + exn_gallery + exn_pes + exn_event +
                      log1p(exn_all_cum) | PeopleID, dt_atyear_foster, cluster = "PeopleID")







## not sure using dt_exhbppl_watyear is good for estimating effect of exhibition:
## this doesn't include those who didn't experience exhibitions -> those are in dt_atyear



## ** illustrative predictions




             
    

dt_sample[exn_event_cum == exn_gallery_cum & exn_gallery_cum > 10]


dt_atyear_foster[exn_pes_cum < 10, .N]/dt_atyear_foster[, .N] ## 70% of artist-years have less than 10 PES exbhs


merge(dt_sample[, .(PeopleID = unique(PeopleID))],
      vadt(fixef(l_mdls$r_pes1)$PeopleID, name = "PeopleID", vlu = "fixef")[, PeopleID := as.integer(PeopleID)],
      by = "PeopleID")
    ## .[which.min(abs(fixef))]


## ** career-stage variation

dt_atyear_foster2 <- dt_atyear_foster %>% copy

l_qntls <- c(0.25,0.5,0.7, 0.9, 0.95,0.98)
dt_atyear_foster2[, exhbD_all_cum := cut(exn_all_cum, quantile(exn_all_cum, probs = c(0,l_qntls)))]
                                         ## labels = paste0("qn_", c(l_qntls*100))]

dt_atyear_foster2[, .N, exhbD_all_cum]

dt_atyear_foster2[is.na(exhbD_all_cum), .(exn_all_cum)]


## * lola
dt_atyear_foster[, uniqueN(PeopleID)]

fepois(exn_pes ~ asinh(exn_pm_cum), dt_atyear_foster)
fepois(exn_pes ~ asinh(exn_pm_cum) | PeopleID  , dt_atyear_foster)

dt_atyear_foster[, .N, exn_pm_cum]

dt_atyear_foster[, pm_ever := fifelse(exn_pm_cum == 0, 0, 1)]
dt_atyear_foster[, .N, pm_ever]

list(pooled = fepois(exn_pes ~ pm_ever, dt_atyear_foster),
     did = fepois(exn_pes ~ pm_ever | PeopleID, dt_atyear_foster)) %>%
    screenreg


rx_year <- fepois(exn_pes ~ pm_ever * factor(year) | PeopleID, dt_atyear_foster[year %between% c(2000, 2020)])




## * dump lit

dt_lit_pmfoster <- gd_lit_foster("pm_foster")
dt_lit_consecgen <- gd_lit_foster("consec_gen")

fwrite(dt_lit_pmfoster, paste0(c_dirs$misc, "lit_pm_foster.csv"))
fwrite(dt_lit_consecgen, paste0(c_dirs$misc, "lit_consecgen.csv"))


## * slopes

dt_coefs <- gd_coefs(l_mdls)


gp_coefs(dt_coefs)

## * matching

## custom nationality handling for multiple nationality: more compact to avoid long match_strings
dt_natgrid <- dt_af_people[, .(nat = unlist(tstrsplit(Nationality, ","))), PeopleID] %>%
    .[, iso3c := countrycode(nat, "country.name", "iso3c",
                             custom_match = c("Kosovo" = "XKX", "Micronesia" = "MCN"))] %>%
    .[order(PeopleID, iso3c)] %>%
    .[, .(Nationality = paste0(iso3c, collapse = ",")), PeopleID]



dt_test <- dt_atyear_foster[, .SD, .SDcols = c("PeopleID", "year", "Gender", "exn_pi", "exn_pm",
                                               "exn_pm_cum", "exn_pm_cum_ihs",
                                               "exn_pi_cum", "exn_pi_cum_ihs",
                                               "BirthYear", "career_year", "exn_all_cum_ihs",
                                               "exn_all_s5")] %>%
    merge(dt_natgrid, by = "PeopleID") %>%
    .[, any_pm := fifelse(max(exn_pm_cum_ihs) == 0, 0, 1), PeopleID] %>%
    .[, Gender := fifelse(Gender %in% c("Male", "Female"), Gender, "Other")] %>%
    .[any_pm == 1, year_first_pm := year[min(which(exn_pm_cum_ihs != 0))], PeopleID] %>%
    .[, time_to_first_pm := year - year_first_pm] %>% 
    .[year %between% c(2000, END_YEAR)] %>%
    .[order(PeopleID, year)]


fwrite(dt_test, "dt_test.csv")
    
dt_test <- fread("dt_test.csv")

## make some smol data.frame for testing purposes
l_at_wpm <- dt_test[exn_pm > 0, unique(PeopleID)] # artist which have at least one PM exhb
l_at_wopm <- setdiff(dt_test$PeopleID, l_at_wpm) # artist which never have PM

dt_test_smol <- dt_test[PeopleID %in% c(sample(l_at_wpm, 1500), sample(l_at_wopm, 1500))]

fwrite(dt_test_smol, "dt_test_smol.csv")



dt_match <- gd_match_history(dt_test[year >= 2000], L=5, return_type = "matches-info")

dt_match[, .N, match_tier]

dt_match[, uniqueN(treated_id), match_tier]

dt_match[, .N, treated_id]

dt_match[, .(nbr_tiers = uniqueN(match_tier)), treated_id][, max(nbr_tiers)]


t1 <- Sys.time()
l_ref <- gd_refine_matches(dt_match, dt_test, covariates = c("BirthYear", "exn_all_cum_ihs"),
                        cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 16))
t2 <- Sys.time()
t2-t1

dt_ref <- l_ref$data




## dt_ref[!is.na(subclass), map(.SD, mean), treat, .SDcols = c("BirthYear", "exn_all_cum_ihs")]

dt_match[, uniqueN(treated_id)] # 10k
dt_ref[!is.na(subclass), uniqueN(subclass)] # 7.7k -> 2.3k dropped


dt_ref_matched <- dt_ref[!is.na(subclass)] %>% 
    .[treat == 0, weight := 1/.N, subclass] %>%
    .[treat ==1, weight := 1]

## look at balance
dt_ref_matched %>% melt(id.vars = c("set", "PeopleID", "treat", "weight"),
                                  measure.vars = c("BirthYear", "exn_all_cum_ihs", "career_year")) %>%
    ggplot(aes(x=value, color = factor(treat), weight = weight)) + 
    geom_density() +
    facet_wrap(~variable, scales = "free")




params <- list(
  L             = 3,
  treatment_col = "exn_pm",
  history_col   = "exn_pm_cum",
  unit_id       = "PeopleID",
  time_id       = "year",
  covariates    = c("BirthYear", "exn_all_cum_ihs"),
  cutpoints     = list(BirthYear = 8, exn_all_cum_ihs = 8)
)

t1 = Sys.time()
qr2 <- gd_run_match_pipeline(dt_test, params, sample_n_treated = 1500)
t2 = Sys.time()
t2-t1

## ** check gd_l1_med

## need small subset

dt_match_smol <- dt_match[treated_id %in% sample(unique(treated_id), 150)]

dt_match_smol_long <- rbind(
    dt_match_smol[, .(PeopleID = control_id, year, treat = 0, set = paste0(treated_id, "_", year))],
    dt_match_smol[, head(.SD, 1), treated_id][, .(PeopleID = treated_id, year, treat = 1,
                                                  set = paste0(treated_id, "_", year))]) %>%
    merge(dt_test, by = c("PeopleID", "year"))
  
## need to make bal.init: prepare data by binning?

m_bal <- bal.init(dt_match_smol_long[, .SD, .SDcols = params$covariates],
                  treat = dt_match_smol_long$treat, stat = "l1.med", l1.n = 10)
bal.compute(m_bal) # computation is instant? 
m_bal$coarsened.covs %>% table %>% names
## hmm weird: preparation is what takes the time, but only one aggregation is returned
## -> is it choosing the median of aggregations, rather than the median of stats of the aggregations? 

## check variability

dt_balvar <- map(seq(80,100,20), ~lapply(1:10, \(nbr)
                                        bal.init(dt_match_smol_long[, .SD, .SDcols = params$covariates],
                                                 treat = dt_match_smol_long$treat, stat = "l1.med", l1.n = .x) %>%
                                        bal.compute %>%
                                        list(vlu = ., nbr = nbr, l1n = .x)) %>% rbindlist) %>% rbindlist

rbind(dt_balvar, dt_balvar_old_1060) %>% # .[, sd(vlu), l1n]
    ## dt_balvar %>%
    ggplot(aes(x=l1n, y = vlu)) + geom_point()


dt_balvar

## 
gd_l1_med_both(dt_match_smol_long[, .SD, .SDcols = params$covariates], dt_match_smol_long$treat,
            weights.after= rep(1, dt_match_smol_long[, .N]), l1.max.bin = 12, l1.n = 101)

gd_l1_med_both(merge(dt_match, dt_test, by = c("PeopleID", "year"))


## ** matching checks
    
## only look at 2010
dt_atyear_foster[year == 2010 & PeopleID == 17753]

## have both control and treatment in long format
dt_match_long <- rbind(
    dt_match[, .(PeopleID = control_id, year, match_key, match_tier, treat = 0,
                      set = paste0(treated_id, "_", year))],
    dt_match[, head(.SD,1), treated_id][, .(PeopleID= treated_id, year, match_key, match_tier, treat = 1,
                                                 set = paste0(paste0(treated_id, "_", year)))]) %>%
    merge(dt_test, by = c("PeopleID", "year")) #


dt_match_smol <- dt_match_long[year == 2010]

## dt_match_long[set == "13918_2010"]
## dt_match_smol[treated_id == 13918]

## only look at one 
dt_match_tiny <- dt_match_long[set == "13918_2010"] 
    

t1 = Sys.time()
## single obs
r_match <- matchit(treat ~ BirthYear + exn_all_cum,
                   data = dt_match_tiny,
                   method = "cem",
                   cutpoints = list(BirthYear = 5, exn_all_cum = 5))
t2 = Sys.time()

## one year
r_match2 <- matchit(treat ~ set + BirthYear + exn_all_cum,
                    data = dt_match_smol,
                    method = "cem",
                    cutpoints = list(BirthYear = 8, exn_all_cum = 8))

## all years
r_match3 <- matchit(treat ~ set + BirthYear + exn_all_cum,
                    data = dt_match_long,
                    method = "cem",
                    cutpoints = list(BirthYear = 8, exn_all_cum = 8))

## matching entire set means that groups gets bigger -> likely more overlap than if groups were based only on single treat

r_match2 %>% summary

dt_match_long[, `:=`(subclass_all = r_match2$subclass)]


dt_match_long[set == "13918_2010", .(PeopleID, treat, subclass_all, subclass_single = r_match$subclass)] %>% summary
    


## ** just look at effect quickly:


l_leads <- 0:5
l_col_leads <- paste0("exn_pi_n", l_leads)
l_y_diff <- paste0("exn_pi_dil", l_leads)

dt_dv <- dt_atyear_foster[, .(PeopleID, year, exn_pi)] %>%
    .[, exn_pi_l1 := shift(exn_pi, 1), PeopleID] %>%
    .[, (l_col_leads) := map(l_leads, ~shift(exn_pi, -.x)), PeopleID] %>%
    .[, (l_y_diff) := map(l_col_leads, ~get(.x) - exn_pi_l1)]

dt_dv[PeopleID == 11, .(PeopleID, year, exn_pi, exn_pi_dil1, exn_pi_dil2, exn_pi_dil3)] %>% print(n=80)


## hmm but what about weights? -> really should check matchit vignettes
dt_cbn <- merge(dt_ref_matched, dt_dv[, .SD, .SDcols = !"exn_pi"] , by = c("PeopleID", "year"))
    

feols(exn_pi_dil ~ treat | set, dt_cbn, weights = ~weight)
