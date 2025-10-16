options(width = 115)

library(pmdata)    # for private museum data
library(jtls)      # custom funcs for tables/figure/number management
library(purrr)     # comfy mapping
library(fixest)    # for fast FE estimation
library(texreg)    # for plotting
library(survival)  # survival models
library(ggsurvfit) # survival models
library(ggplot2)   # plotting
library(countrycode)
library(igraph)    # for flow graph calculations
library(splines)   # for fancy convoluted regression models
library(marginaleffects) # for good predictions
library(Hmisc, include.only = "latexTranslate") # for dealing with nasty names
library(collapse, include.only = c("fcumsum", "funique", "join", ".c"))     # for fast cumsum in data prep
library(scales, include.only = "label_number") # comfy suffix labels (e.g. 2K, 2B )for ranges across magnitudes
library(broom, include.only = "tidy") # model extracting
set.seed(42)


DIR_PROJ <- "~/Dropbox/phd/papers/foster/"
c_dirs <- gc_dirs(dir_proj = DIR_PROJ)
c_dirs <- c(c_dirs, list(tempdata = paste0(DIR_PROJ, "tempdata/")))

l_tbls <- list()
l_plts <- list()
c_tblargs <- list()
c_pltargs <- list()
system(sprintf("rm %s", paste0(c_dirs$misc, "farg_calls.csv")))


PMDATA_LOCS <- gc_pmdata_locs()

END_YEAR <- 2022

l_funcfiles <- c(
    paste0(c_dirs$code, "cfg.R"),
    paste0(c_dirs$code, "data_prep.R")
    ## paste0(c_dirs$code, "regression.R"),
    ## paste0(c_dirs$code, "descriptives.R"), ## FIXME: first into plots/figures
    ## source(paste0(c_dirs$code, "network.R")) # don't need network measures atm
    ## paste0(c_dirs$code, "rem_funcs.R")
    ## paste0(c_dirs$code, "network.R") ## FIXME: add after basic stuff works
)

invisible(lapply(l_funcfiles, \(x) source(x, echo = F, verbose = F, print.eval = F)))


if (interactive()) {stop("it's interactive time")}
## ** main


## FIXME: make independent from infl


dt_pmdb <- gd_pmdb_excl(sel = "any_ever_pm") %>% gd_pmdb

## read in AF data, make IDs a bit more meaningful
dt_af_exhbs <- gd_af_exhbs() %>% setnames(old = "ID", new = "ExhibitionID") %>% 
    .[, begin_year := year(BeginDate)]

dt_af_instns <- gd_af_instns() %>% setnames(old = "ID", new = "InstitutionID") %>%
    .[Country %!in% c("NULL", "Others")]

dt_af_people <- gd_af_people() %>% setnames(old = "ID", new = "PeopleID") %>%
    .[Nationality != "NULL"] # no actual nationality

dt_af_links_exhb_ppl <- gd_af_links_exhb_ppl() 

dt_af_pmdb_matches <- gd_af_pmdb_matches() %>%
    .[dt_pmdb[museum_status %in% c("private museum", "closed"), .(PMDB_ID = ID)], on = "PMDB_ID"]


## set up flags/switches how to convert data
c_datargs <- list(
    clsfcn_to_use = "clsfcn4",
    b_rcd_to_pm = T,
    b_sample_dt_atyear = T)


c_filters_atyear <- list(
    keep_all = expr(T),
    yeet_baby_exhbs = expr(begin_year > BirthYear + 16),
    yeet_babies_completely = expr(first_year > BirthYear + 16))


## auction data
dt_ls_aucres <- gd_ls_aucres() %>% adt
dt_af_ls_matches <- gd_af_ls_matches()
dt_ls_af_aucres <- gd_ls_af_aucres(dt_ls_aucres, dt_af_ls_matches)


dt_af_instns_rcd <- md_instns_rcd(c_datargs$clsfcn_to_use, dt_af_instns, c_datargs$b_rcd_to_pm, dt_af_pmdb_matches)

## --------- artist-year data
dt_atyear_raw <- gd_atyear_raw(dt_af_exhbs, dt_af_people, dt_af_links_exhb_ppl, dt_af_instns_rcd, dt_ls_af_aucres)





## calculate some before/after steps
dt_atyear_foster <- gd_atyear_foster(dt_atyear_raw, c_filters_atyear$yeet_babies_completely)

# dt_exhbppl_sample <- dt_exhbppl[sample(1:.N, size = 1e5)]

## hacky way of integrating EVC, if decided to keep integrate better

source(paste0(c_dirs$code, "network.R"))
dt_edges <- gd_edges(dt_exhbppl)
dt_edges_wtd <- gd_edges_wtd(dt_edges)

g_insflow <- gg_insflow(dt_edges_wtd)
dt_insevc <- gd_insevc(g_insflow, dt_af_instns_rcd)

dt_insevc[order(-evc), .(Name, evc)] # moma, Zentrum für Kunst und Medien Karlsruhe, vienna huh

dt_exhbppl_evc <- gd_exhbppl_evc(dt_exhbppl, dt_insevc)
## dt_ppl_evcsum <- gd_ppl_evcsum(dt_exhbppl_evc)
dt_atyear_evc <- gd_atyear_evc(dt_exhbppl_evc, dt_atyear)
## dt_atofyear <- gd_atofyear(dt_atyear_evc)

## check EVC developments: looks ok? 
dt_atyear_evc[, c(list(PeopleID = PeopleID, year = year), .SD), .SDcols = patterns("^evc")] %>%
    print(n=50)


dt_atyear_foster2 <- merge(dt_atyear_foster,
                           dt_atyear_evc[, c(list(PeopleID = PeopleID, year = year), .SD),
                                         .SDcols =patterns("^evc")],
                           by = c("PeopleID", "year"))




## manually switch in EVC; fuck me that's ugly
dt_exhbppl_watyear <- merge(
    dt_exhbppl,
    ## dt_atyear_foster[, c(list(PeopleID = PeopleID, begin_year = year), .SD), .SDcols = patterns("exhbN")],
    dt_atyear_foster2[, c(list(PeopleID = PeopleID, begin_year = year), .SD), .SDcols = patterns("exhbN|evc")],
    by = c("PeopleID", "begin_year")) %>%
    .[order(PeopleID, begin_year)]

## minor differences between cnt_prevexhb and exhbN_all_cum: exhbN_all_cum is only previous year?
## hmm not quite: diffs can also be negative -> FIXME: make sure there's no ambiguity

dt_exhbppl_watyear[, .(PeopleID, begin_year, cnt_prevexhb, exhbN_all_cum,
                       diff_cum = cnt_prevexhb - exhbN_all_cum)] %>%
    .[, .N, diff_cum] %>% .[order(diff_cum)] %>% print(n=80)

## now can calculate institution-specific effect: total count, relative
# how many new exhibitions are in the next 5 years
dt_exhbppl_watyear[, `:=`(exhbN_grow5_cnt = exhbN_all_cum_n5 - cnt_prevexhb,
                          # by how much does the total exhibition count grow
                          exhbN_grow5_perc =  100*((exhbN_all_cum_n5/cnt_prevexhb)-1),
                          evc_grow5_cnt = evc_n5 - evc_s5,
                          evc_grow5_prop = evc_n5/evc_s5) # hmm do i need log here? 
                          ]

dt_exhbppl_watyear[PeopleID ==44 & begin_year > 2005, .(PeopleID, begin_year, exhbN_all, exhbN_all_cum_n5,
                                                        cnt_prevexhb, exhbN_grow5_cnt, exhbN_grow5_perc,
                                                        evc_grow5_cnt, evc_grow5_perc)] %>%
    print(n=40)


dt_exhbppl_watyear[PeopleID ==44 & begin_year > 2005, .(PeopleID, begin_year, exhbN_all, evc_ay_ttl, evc_s5,
                                                        evc_n5, evc_grow5_cnt, evc_grow5_prop)] %>%
    print(n=80)
    

dt_exhbppl_watyear[, .(mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm = T))]

l_cols_summary <- c("exhbN_grow5_cnt", "exhbN_grow5_perc", "evc_grow5_cnt", "evc_grow5_prop")

dt_exhbppl_watyear[, .(mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm = T),
                       mean_grow5_perc = mean(exhbN_grow5_perc[!is.infinite(exhbN_grow5_perc)], na.rm = T)),
                   InstitutionType]

dt_exhbppl_watyear[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)), InstitutionType, .SDcols = l_cols_summary]

## super weird EVC values: maybe extreme cases?
dt_exhbppl_watyear %>% ggplot(aes(x=evc_grow5_prop, color = InstitutionType)) + geom_density() +
    scale_x_continuous(trans = scales::log10_trans()) +
    coord_cartesian(xlim = c(1e-3,1e3))


dt_exhbppl_watyear[is.infinite(exhbN_grow5_perc), uniqueN(PeopleID)]

dt_atyear_foster

## distribution of institution-level growth factors (ability of individual museums to promote)
dt_insgrow5 <- dt_exhbppl_watyear %>%
    .[, .(mean_grow5_perc = mean(exhbN_grow5_perc[!is.infinite(exhbN_grow5_perc)], na.rm =T),
          mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm =T)), .(InstitutionID, InstitutionType)]

dt_insgrow5 <- dt_exhbppl_watyear[, map(.SD, ~mean(.x[!is.infinite(.x)], na.rm = T)),
                                  .(InstitutionType, InstitutionID), .SDcols = l_cols_summary]

## growth count by growth perc: basically shows that PM grow less percent for same count -> orientation towards established
dt_insgrow5 %>% # [sample(1:.N, size = 5000)] %>%
    ggplot(aes(y=log10(mean_grow5_perc), x = log10(mean_grow5_cnt), color = InstitutionType)) +
    geom_point(size = 0.1, alpha = 0.1) +
    geom_smooth(se = F)


## https://stackoverflow.com/questions/48282989/show-only-high-density-areas-with-ggplot2s-stat-density-2d
dt_insgrow5 %>% ggplot(aes(log10(mean_grow5_perc), x = log10(mean_grow5_cnt), group = InstitutionType)) +
    stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = InstitutionType), bin = 1)
## idk if 2d distribution is necessary, facetted dens plot show it pretty much as well


## distribution of percentage count/growth of institutions by institution type
melt(dt_insgrow5, id.vars = c("InstitutionID", "InstitutionType"), variable.name = "mes", value.name = "vlu") %>%
    ggplot(aes(x=vlu, color = InstitutionType)) + geom_density(alpha = 0.1) +
    facet_wrap(~mes, ncol =1, scale = "free") +
    scale_x_continuous(trans = scales::log10_trans())
    



# most influential orgs by type in terms of percent
dt_insgrow5[order(-mean_grow5_perc)][, head(.SD), InstitutionType] %>% 
    dt_af_instns[, .(InstitutionID, Name, City, Country)][., on = "InstitutionID"] %>% print(n=80)

## most influential in exhibition count
dt_insgrow5[order(-mean_grow5_cnt)][, head(.SD), InstitutionType] %>% 
    copy(dt_af_instns)[, InstitutionType := NULL][., on = "InstitutionID"] %>% print(n=80)



