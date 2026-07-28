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

## add medium sized test
dt_test_med <- dt_test[PeopleID %in% c(sample(l_at_wpm, 3000), sample(l_at_wopm, 3000))]
fwrite(dt_test_med, "dt_test_med.csv")

dt_test_big <- dt_test[PeopleID %in% c(sample(l_at_wpm, 6000), sample(l_at_wopm, 6000))]
fwrite(dt_test_big, "dt_test_big.csv")

dt_med <- fread("~/Dropbox/phd/papers/foster/dt_test_med.csv")
dt_big <- fread("~/Dropbox/phd/papers/foster/dt_test_big.csv")

list(dt = dt_med, L = 5L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
       uid = "PeopleID", tid = "year",
       covariates = c("BirthYear", "exn_all_cum_ihs"),
       cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
     slow = TRUE, label = "med_L5")

list(dt = dt_big, L = 5L, treatment_col = "exn_pm", history_col = "exn_pm_cum",
       uid = "PeopleID", tid = "year",
       covariates = c("BirthYear", "exn_all_cum_ihs"),
       cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       slow = FALSE, label = "big_L5")




dt_match <- gd_match_history(dt_test, L=5, history_col = "exn_pm_cum",
                             treatment_col = "exn_pm", uid = "PeopleID", tid = "year",
                             return_type = "matches-info")

dt_match[, .N, match_tier]
dt_match[, uniqueN(treated_id), match_tier]
dt_match[, .N, treated_id]
dt_match[, .(nbr_tiers = uniqueN(match_tier)), treated_id][, max(nbr_tiers)]



l_ref <- gd_refine_matches(dt_match, dt_test, covariates = c("BirthYear", "exn_all_cum_ihs"),
                           cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 16),
                           uid = "PeopleID", tid = "year")


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
  uid       = "PeopleID",
  tid       = "year",
  covariates    = c("BirthYear", "exn_all_cum_ihs"),
  cutpoints     = list(BirthYear = 8, exn_all_cum_ihs = 8),
  tiers = list(c("hist_str", "Nationality", "Gender"), c("hist_str")),
  sample_n_treated = 100
)


qr2 <- gd_run_match_pipeline(dt_test, params)

## sample n treated grid
## dt_snt_grid <- expand.grid(run = 1:5, snt = c(100, 250, 500, 1000)) %>% adt

## dt_res_snt <- map(dt_snt_grid$snt,
##                   ~gd_run_match_pipeline(dt_test, params, sample_n_treated = .x)[, snt := .x]) %>% rbindlist

## dt_res_snt %>% copy %>%
##     .[, run := 1:.N] %>%
##     .[, L1_diff := L1_before - L1_after] %>%
##     melt(id.vars = c("snt", "run"), measure.vars = c("L1_after", "L1_before", "L1_diff", "cem_match_rate",
##                                                      ## "n_treated", "n_matched", "match_rate", # constant
##                                                      "smd_BirthYear_before", "smd_BirthYear_after"
##                                                      )) %>% 
##     ggplot(aes(x = factor(snt), y= value, color = variable)) + geom_point() +
##     facet_grid(variable ~., scale = "free_y")

## ** systematic testing of sample_n_treated requirements

## general idea:
## no need to run entire pipeline for new cfg
## for stuff like sample_n_treated, i can run those step once, and sample from refined objects


dt_match_full <- gd_match_history(
    dt_test,
    L             = params$L,
    treatment_col = params$treatment_col,
    history_col   = params$history_col,
    uid       = params$uid,
    tid       = params$tid,
    return_type   = "matches-info",
    tiers = list(c("hist_str", "Nationality", "Gender"), c("hist_str", "Nationality"),
                 c("hist_str", "Gender"), "hist_str")
    ## tiers =     params$tiers
    ## tiers = list("hist_str") # big yikes, results in 280m results
)

refined_full <- gd_refine_matches(dt_match_full, dt_test,
                                  covariates = params$covariates,
                                  cutpoints  = params$cutpoints,
                                  uid    = params$uid,
                                  tid = params$tid)



gd_quality_rep_min <- function(dt_panel, refined, params, sample_n_treated) {
    #' sample sample_n_treated and run match stats on dt_ref
    #' matching and refinement is easily doable with entire df, bottleneck is gof
    #' can save the entire matching and refining step, and just re-sample the refined dt? 

    ## get ids of sample of refined
    treated_all <- funique(refined$data$treated_id)
    l_keep <- sample(treated_all, sample_n_treated)
    ids_ref <- which(refined$data$treated_id %in% l_keep)

    ## dt_ref_sample <- dt_match[.(l_keep), on = .(treated_id), nomatch= NULL]

    ## create new refined object (data and matchit (in turn: weights, subclass))
    refined_new <- list(
        data = refined$data[ids_ref],
        matchit = list(weights = refined$matchit$weights[ids_ref],
                       subclass = refined$matchit$subclass[ids_ref]))

    ## generate quality report
    gd_quality_report(dt_test, dt_match,
                    matchit_obj  = refined_new$matchit,
                    refined_data = refined_new$data, 
                    params       = params)

}

## microbenchmark(
##     minrep = gd_quality_rep_min(dt_test, refined_full, params, sample_n_treated = 1000), # 5.2
##     pipeline = gd_run_match_pipeline(dt_test, params, sample_n_treated = 1000), # 7 -> use minrep
##     times = 10)


## set up grid
dt_snt_grid <- expand.grid(run = 1:7, snt = c(30, 100, 300, 1000, 3000)) %>% adt

## actual testing
dt_res_snt <- map(dt_snt_grid$snt,
                  ~gd_quality_rep_min(dt_test, refined_full, params, sample_n_treated = .x)[, snt := .x]) %>%
    rbindlist

## eval
dt_res_snt %>% copy %>%
    .[, run := 1:.N] %>%
    .[, L1_diff := L1_before - L1_after] %>%
    melt(id.vars = c("snt", "run"), measure.vars = c("L1_after", "L1_before", "L1_diff", "cem_match_rate",
                                                     ## "n_treated", "n_matched", "match_rate", # constant
                                                     "smd_BirthYear_before", "smd_BirthYear_after", "tdiff"
                                                     )) %>% 
    ggplot(aes(x = factor(snt), y= value, color = variable)) + geom_point() +
    facet_grid(variable ~., scale = "free_y")


l_cfgs <- gd_build_configs(
  dt = dt_test,
  L = list(3L, 5L),
  treatment_col = list("exn_pm"),
  history_col = list("exn_pm_cum"),
  uid = list("PeopleID"),
  tid = list("year"),
  ## covariates = list(c("BirthYear", "exn_all_cum_ihs"), c("BirthYear")),
  cutpoints = list(list(BirthYear = 8, exn_all_cum_ihs = 8), list(BirthYear = 8)),
  tiers = list(list(c("hist_str", "Nationality", "Gender"), c("hist_str"))),
  sample_n_treated = list(100)
)

l_res <- gd_run_configs(l_cfgs, dt_full = dt_test)

run_match_pipeline(l_cfgs
gd_pipe


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


dt_test <- data.table(stratum = c("a", "a", "b", "b", "b", "c", "c", "c", "d", "d"),
                      treat = factor(c(1,2,1,2,3,1,2,3, 2, 3)),
                      outcome = rnorm(n=10))

feols(outcome ~ treat | stratum, dt_test)
feols(outcome ~ treat | stratum, dt_test[stratum != "d"])


library(lintr)
lint("~/Dropbox/phd/papers/foster/scripts/mymatch.R")




res <- gd_cem_strata(cem_configs[[1]])

gd_cem_quality_report(res) %>% t

res$strata_summary[n_total > 1] %>%
.[, nbr_grps := rowSums(.SD > 0), .SDcols = paste0("n_", c(res$config$treatment_cols, "control"))] %>%
.[nbr_grps > 1, sum(n_total)]

## callgraph (inconvenient af)
options(browser = "surf")
fw <- foodweb(filter = T)
plot(fw)

foodweb(FUN = gd_run_match_pipeline) %>% plot




## some fake data
l_breaks <- gd_coarse_breaks(cfg$dt_panel, cfg$cutpoints)
r_dt2 <- gd_refine_matches(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid,
                              cutpoints = l_breaks)
  r_mi2 <- gd_refine_matches_matchit(mi, cfg$dt_panel, uid = cfg$uid, tid = cfg$tid,
                                      cutpoints = l_breaks)



source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

dt_full <- tribble(
    ~PeopleID,  ~year,  ~BirthYear,  ~exn_all_cum_ihs,  ~treated_id,  ~control_id,  
    1L,         2020L,  1995L,       9,                 1,            1,
    101L,       2020L,  1930L,       1,                 1,            101,
    102L,       2020L,  1950L,       2.5,               1,            102,
    103L,       2020L,  1960L,       4,                 1,            103,
    104L,       2020L,  1975L,       6,                 1,            104,
    105L,       2020L,  1990L,       8,                 1,            105) %>% adt

dt_cem <- dt_full[PeopleID != 1]

## --- Our equal-width (4 bins) ---
cp_int <- list(BirthYear = 4L, exn_all_cum_ihs = 4L)

gd_refine_matches(dt_cem, dt_full, uid = "PeopleID", tid= "year", cutpoints = cp_int)$data[, .(PeopleID, subclass)]

gd_refine_matches_matchit(dt_cem, dt_full, uid = "PeopleID", tid= "year", cutpoints = cp_int)$data[, .(PeopleID, subclass)]



## ** equal spaced breaks

dt <- fread("dt_test_smol.csv")

## MatchIt CEM with same cutpoints
x <- dt$exn_pi_cum_ihs
cp_int <- 5L  ## same as our 5-bin quantile

## run matchit CEM to see what it does internally
form <- as.formula("treat ~ exn_pi_cum_ihs")
dt_tmp <- data.table(treat = sample(0:1, nrow(dt), replace = TRUE), exn_pi_cum_ihs = x)
set.seed(42)
r <- matchit(form, data = dt_tmp, method = "cem", cutpoints = list(exn_pi_cum_ihs = 5L))

## compute our coarse_bin version
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

table(r$subclass, our_bin) # yeah complete overlap just different order

## ** equal-size bins

l_br <- gd_coarse_breaks(dt_tmp, cutpoints = list(exn_pi_cum_ihs = 5))
r_mi_eqsb <- matchit(form, data = dt_tmp, method = "cem", cutpoints = l_br)
r_mi_eqsb$subclass %>% table


## l_br <- gd_coarse_breaks(dt_tmp, cutpoints = list(exn_pi_cum_ihs = 10))
## r_mi_eqsb <- matchit(form, data = dt_tmp, method = "cem", cutpoints = l_br)
## r_mi_eqsb$subclass %>% table
## ## matchit can take custom cutpoints


table(r$subclass, r_mi_eqsb$subclass)

r_coarse_bin <- coarse_bin(x, breaks = l_br$exn_pi_cum_ihs)
table(r_coarse_bin)
table(r_mi_eqsb$subclass)

table(r_coarse_bin, r_mi_eqsb$subclass) 


r_mi_eqsb$subclass %>% table %>% sum # 57971
r_coarse_bin %>% table %>% sum # 57972

x # 57972
dt_tmp # 57972

r_mi_eqsb$subclass[is.na(r_mi_eqsb$subclass)]
r_coarse_bin[1126] # the non-matched in matchit is the singleton in coarse_bin

## try with new coarsening function
l_br_new <- gd_coarse_breaks(dt_tmp, cutpoints = list(exn_pi_cum_ihs = 5))

r_mi_eqsb_new <- matchit(form, data = dt_tmp, method = "cem", cutpoints = l_br_new)


table(r_mi_eqsb$subclass, r_mi_eqsb_new$subclass)# same res

## * check feasability of MI
dt_big_mi <- dt_big %>% copy %>% .[, treat := as.integer(collapse::finteraction(as.integer(exn_pm > 1), as.integer(exn_pi > 1)) != "0.0")]
dt_big_mi[, .N, treat]

r_big <- matchit(treat ~ exn_pm_cum_ihs + Gender + Nationality, dt_big_mi,  method = "cem")
summary(r_big) # still too slow: takes at least 3 min


## * re-assess cem matching


c_cem1 <- list(dt_panel = dt_big, L = 3L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = c(5L, 5L), history_breaks = "quantile",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(BirthYear = 8, exn_all_cum_ihs = 8),
       label = "quantile_5", weight_type = "composite",
       outcome_cols = list(exn_pi = 3, exn_pm = 3))

r_cem <- gd_cem_strata(c_cem1)
dt_res_cem <- gd_cem_quality_report(r_cem) 
gp_cem_quality(r_cem)


## many zero weights... concerning? 
r_cem$data[weight < 5 & !(exn_pm == 0 & exn_pi == 0)] %>% ggplot(aes(x=weight, color = interaction(exn_pm > 0, exn_pi > 0))) + geom_density()


## manually inspect density to assess balance
r_cem$data %>% ggplot(aes(x=BirthYear, color = interaction(as.integer(exn_pm>0),
                                                           as.integer(exn_pi>0)))) + geom_density()

r_cem$config$dt_panel %>% ggplot(aes(x=BirthYear, color = interaction(as.integer(exn_pm>0),
                                                           as.integer(exn_pi>0)))) + geom_density()

r_cem$data %>% .[, exn_pm_bin := as.integer(exn_pm > 0)] %>% 
ggplot(aes(x=exn_pm_cum_ihs, color = factor(exn_pm_bin))) + geom_density()
                                                    


## pre-trend figure
dt_lags <- cbind(r_cem$config$dt_panel[, .(year)],
      collapse::L(r_cem$config$dt_panel, by = ~PeopleID, n=1:3, cols = "exn_pi"))

join(r_cem$data, dt_lags, on = c("PeopleID", "year")) %>% 
.[, map(.SD, ~mean(.x, na.rm = T)), exn_pm_bin, .SDcols = paste0("L", 1:3, ".exn_pi")]


join(r_cem$data, dt_lags, on = c("PeopleID", "year")) %>%
.[year > 2015] %>% 
.[, c(list(N = .N), map(.SD, ~weighted.mean(.x, na.rm = T, weight = weight))), .(exn_pm_bin, stratum),
  .SDcols = paste0("L", 1:3, ".exn_pi")] %>% print(n=10) %>% 
.[, N_treat := uniqueN(exn_pm_bin), stratum] %>%
.[N_treat > 1] %>% na.omit %>% 
.[order(stratum, exn_pm_bin)] %>% print(n=20) %>%
melt(id.vars = c("stratum", "exn_pm_bin"), measure.vars = paste0("L", 1:3, ".exn_pi")) %>%
dcast(stratum + variable ~ exn_pm_bin) %>%
.[, ratio := `0`/`1`] %>%
.[!is.infinite(ratio) & !is.na(ratio)] %>% print(n=10) %>% 
.[, .(value = mean(ratio)), variable]

ggplot(aes(x

dcast(stratum ~ exn_pm_bin, value.var = paste0("L", 1:3, ".exn_pi")) %>%
.[, r_l1 := 


## huh there's something super wrong with career data
## exn_pi_cum_ihs is zero for sooooooo many.. 

join(r_cem$data[stratum == "2016_2_2_2_1_1_1_8_2", .(PeopleID, year)], dt_big, on = c("PeopleID", "year")) %>%
print(n=120) %>% summary

dt_big[exn_pi_cum_ihs == 0]

## hmm maybe it's just that strata? 2016_2_2_2_1_1_1_8_2 mean 2016, cat2 on hist1, cat1 on hist1, cat8 on birthyear (80s), cat 2 on all_ihs 
dt_big[PeopleID == 362685, .(year, exn_all_cum_ihs)]

dt_pretrend_bal <- gd_pre_trend_balance(r_cem, dt_big, outcome_vars = c("exn_pi", "exn_pm", "BirthYear"),
                                        max_lags = 3)

gp_pre_trend(dt_pretrend_bal)


## ** finer coarsening

c_cem_finer <- list(dt_panel = dt_big, L = 3L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = c(10L, 10L), history_breaks = "quantile",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(BirthYear = 10, exn_all_cum_ihs = 10),
       label = "finer", weight_type = "composite")


r_cem_finer <- gd_cem_strata(c_cem_finer)
dt_res_cem_finer <- gd_cem_quality_report(r_cem_finer) 

p_cem1 <- gp_cem_quality(r_cem)

p_cem_finer <- gp_cem_quality(r_cem_finer)
X11()
p_cem_finer


## ** even finer

c_cem_sieve <- list(dt_panel = dt_big, L = 6L,
       treatment_cols = c("exn_pm", "exn_pi"),
       history_cols = c("exn_pm_cum_ihs", "exn_pi_cum_ihs"),
       history_nbins = c(10L, 10L), history_breaks = "quantile",
       uid = "PeopleID", tid = "year",
       covariate_cutpoints = list(BirthYear = 10, exn_all_cum_ihs = 10),
       label = "finer", weight_type = "composite")

r_cem_sieve <- gd_cem_strata(c_cem_sieve)
p_cem_sieve <- gp_cem_quality(r_cem_sieve)




## * weight inspection

dt <- data.table(stratum = c(1,1,1,1, 2,2,2, 2,2,2),
                 treat   = c(1,1,0,0, 1,0,0,0,0,0))


dt_s <- dt[, .(msT = sum(treat), # number of treat in stratum
               msC = sum(treat == 0)), # number of control in stratum
           by = stratum]
mT <- sum(dt_s$msT)  ## 3 treated (nbr total matched)
mC <- sum(dt_s$msC)  ## 4 (7) control (nbr total controls)


## CEM ATT weights (Iacus, King, Porro 2012, eq 6):
##   treated:    w = 1
##   control:    w = (mC / mT) × (msT / msC)
dt[dt_s, on = "stratum", `:=`(msT = i.msT, msC = i.msC)]
dt[treat == 1, w := 1]
dt[treat == 0, w := (mC / mT) * (msT / msC)]

print(dt[, .(stratum, treat, msT, msC, w)])

## Check: within each stratum, sum(treated weight) vs sum(control weight)
dt[, `:=`(Wt = sum(w[treat == 1]), Wc = sum(w[treat == 0])), by = stratum]
unique(dt[, .(stratum, Wt, Wc, ratio = round(Wc / Wt, 2))])


r <- matchit(treat ~ stratum, data = dt, method = "cem",
             cutpoints = list(stratum = NULL))


cbind(dt, weight = r$weights)


dt[, mean(w), treat] # hmm weights by treatment average up to 1 per treatment


dt[, mean(w), .(treat, stratum)]
dt[, sum(w), .(treat, stratum)]

## hmm the weights are only for controls anyways.. and often i won't have any controls... 

## think how this can be expanded to multiple treatment cols

## * benchmarking
## question: how expensive is moving data into and out of C (collapse)
## idea: running loops/maps is slow because of requirement to move data to collapse and back to R
## theory: keeping more stuff in collapse will increase speed
## test: different amount of rows and cols

library(microbenchmark)

dt_big_dup <- map(1:10, ~dt_big) %>% rbindlist

dt_yug <- map(1:50, ~dt_big) %>% rbindlist

dt_bench <- microbenchmark(
    single = dt_big[, collapse::fmean(exn_pm_cum, g = PeopleID)],
    mult =  dt_big[, collapse::fmean(.SD, g = PeopleID),
                   .SDcols = c("exn_pm_cum", "exn_pi_cum", "exn_pm_cum_ihs", "exn_pi_cum_ihs")],
    all = dt_big[, collapse::fmean(.SD, g = PeopleID), .SDcols = collapse::num_vars(dt_big, return = "names")],
    single_dup = dt_big_dup[, collapse::fmean(exn_pm_cum, g = PeopleID)],
    mult_dup =  dt_big_dup[, collapse::fmean(.SD, g = PeopleID),
                   .SDcols = c("exn_pm_cum", "exn_pi_cum", "exn_pm_cum_ihs", "exn_pi_cum_ihs")],
    all_dup = dt_big_dup[, collapse::fmean(.SD, g = PeopleID),
                         .SDcols = collapse::num_vars(dt_big, return = "names")],
    single_yug = dt_yug[, collapse::fmean(exn_pm_cum, g = PeopleID)],
    mult_yug =  dt_yug[, collapse::fmean(.SD, g = PeopleID),
                   .SDcols = c("exn_pm_cum", "exn_pi_cum", "exn_pm_cum_ihs", "exn_pi_cum_ihs")],
    all_yug = dt_yug[, collapse::fmean(.SD, g = PeopleID), .SDcols = collapse::num_vars(dt_big, return = "names")],
    times = 20)

## results
## + Unit: milliseconds
##        expr       min        lq      mean    median        uq        max neval   cld
##      single  2.366411  2.421342  2.525588  2.544428  2.581145   2.828557    50 a    
##        mult  3.876872  3.967331  4.256892  4.088374  4.160908   9.332238    50 a    
##         all  8.530757  8.582639  9.294016  8.730555  8.849160  20.366635    50  b   
##  single_dup 12.654747 12.844385 13.814596 13.491791 13.688410  23.974125    50   c  
##    mult_dup 25.580467 25.906255 27.327870 26.541200 26.836766  37.465952    50    d 
##     all_dup 69.942579 70.991292 74.509336 71.563935 72.848706 164.584789    50     e
 ## single_yug  73.176851  73.388959  78.715957  74.438535  76.338418 115.908775    20     e  
 ##   mult_yug 137.121166 138.109065 141.405213 139.274065 140.992090 162.678472    20      f 
 ##    all_yug 360.097136 360.654044 368.378474 361.921855 372.813254 399.816087    20       g


## both scale better than linear, but columns scale better: 15x data takes only 3.7x for small data and 5.8x for big data
## 10 times the rows take 5.4x for small, 8.2x for large
## hmm the bigger the data, the less clear the benefit? 

dt_bench2 <- dt_bench %>% adt %>% .[, .(time = mean(time)/1e6), expr] %>%
.[, `:=`(width = str_extract(expr, "single|mult|all"),
         length = str_extract(expr, "yug|dup"))] %>%
    .[is.na(length), length := "smol"] %>%
    .[, length := factor(length, levels = c("smol", "dup", "yug"))] %>% 
    .[order(length, -width)] %>%
    .[J(length = c("smol", "dup", "yug"), mult = c(1,10,50)), on = "length"]

join(dt_bench2, dt_bench2[width == "single", .(length, time_single = time)], on = "length") %>%
join(dt_bench2[length == "smol", .(width, time_smol = time)], on = "width") %>%
.[, r_width := time/time_smol] %>%
.[, r_prop := 1-(r_width/mult)]

## not really actually, or rather depends on number of columns: stronger improvements with single/few columns (~40% speedup) than with all columns (15) (~20% speedup)

## * DiD investigation

library(fixest)

r_cem_wdv <- gd_add_outcome_next(r_cem)
r_cem_wdv$data

## basic combined
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pi + exn_pm | stratum, r_cem_wdv$data)

## single 
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pm | stratum, r_cem_wdv$data)
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pi | stratum, r_cem_wdv$data)

## binary
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ as.integer(exn_pm > 0) + as.integer(exn_pi > 0) | stratum, r_cem_wdv$data)

## no stratum
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pi + exn_pm, r_cem_wdv$data)

## interaction
feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pi*Gender + exn_pm*Gender | stratum, r_cem_wdv$data)
## some effects 

feols(c(exn_pi_d1, exn_pi_d2, exn_pi_d3) ~ exn_pi + exn_pm | stratum, r_cem_wdv$data, weights = ~weight)


## check strata variation

collapse::varying(r_cem_wdv$data,  ~stratum)

c_cem1

r_cem_wdv$data[PeopleID == 44] %>% view_xl

r_cem_wdv$config$dt_panel[PeopleID == 44]

## need someone who starts after 2000 but has many exhbs

gp_strata_dist(r_cem_finer)
X11()
gp_strata_dist(r_cem)

