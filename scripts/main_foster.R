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
library(MASS, include.only = "kde2d")
library(zeallot)
library(stringr)
set.seed(42)

l_asinh_breaks <- c(0, 10^(0:10), 3*10^(0:10)) %>% c(., -1*.) %>% sort # breaks for asinh ggplots

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
    paste0(c_dirs$code, "data_prep.R"),
    paste0(c_dirs$code, "panel.R"),
    paste0(c_dirs$code, "descriptives.R"), 
     ## paste0(c_dirs$code, "rem_funcs.R")
    paste0(c_dirs$code, "network.R") 
)

invisible(lapply(l_funcfiles, \(x) source(x, echo = F, verbose = F, print.eval = F)))


if (interactive()) {stop("it's interactive time")}
## ** main


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


dt_exhbppl_wyeet <- gd_exhbppl(dt_af_people, dt_af_exhbs, dt_af_instns_rcd, dt_ls_af_aucres, dt_af_links_exhb_ppl)

dt_exhbppl <- dt_exhbppl_wyeet[InstitutionType != "yeet"] %>%
    .[eval(c_filters_atyear$yeet_babies_completely)]

# dt_exhbppl_sample <- dt_exhbppl[sample(1:.N, size = 1e5)]

## hacky way of integrating EVC, if decided to keep integrate better



dt_edges <- gd_edges(dt_exhbppl, sample_artists = 0)
dt_edges_wtd <- gd_edges_wtd(dt_edges) %>% .[InstitutionID != InstitutionID2]
g_insflow <- gg_insflow(dt_edges_wtd)


dt_insevc <- gd_insevc(g_insflow, dt_af_instns_rcd)
dt_insevc[order(-evc), .(Name, evc)] # moma, Zentrum für Kunst und Medien Karlsruhe, vienna huh

dt_exhbppl_evc <- gd_exhbppl_evc(dt_exhbppl, dt_insevc)

dt_atyear_evc <- gd_atyear_evc(dt_exhbppl_evc, dt_atyear_foster)


## check EVC developments: looks ok? 
dt_atyear_evc[, c(list(PeopleID = PeopleID, year = year), .SD), .SDcols = patterns("^evc")] %>%
    print(n=50)


dt_exhbppl_watyear <- gd_exhbppl_watyear(dt_exhbppl, dt_atyear_evc)

## "evc_grow5_cnt", "evc_grow5_prop"
l_cols_summary <- c("exhbN_all_cum", "exhbN_all_s5", "exhbN_grow5_cnt", "exhbN_grow5_perc", "exhbN_perc_n5s5", "exhbN_perc_n5s5_org")


dt_insgrow5 <- gd_insgrow5(dt_exhbppl_watyear, l_cols_summary)

l_mdls <- gl_mdls(dt_atyear_foster, dt_atyear_evc)


## dt_atyear_evc[, map(.SD, ~sum(is.na(.x)))] %>% t %>% adt %>% .[, .N, V1]





## super weird EVC values: maybe extreme cases?
## dt_exhbppl_watyear %>% ggplot(aes(x=evc_grow5_prop, color = InstitutionType)) + geom_density() +
##     scale_x_continuous(trans = scales::log10_trans()) +
##     coord_cartesian(xlim = c(1e-3,1e3))


dt_exhbppl_watyear[is.infinite(exhbN_grow5_perc), uniqueN(PeopleID)]


## distribution of institution-level growth factors (ability of individual museums to promote)

    
walk(names(gc_tbls()), ~lapply(c(gtbl, wtbl), \(f) f(.x)))
walk(names(gc_plts()), ~lapply(c(gplt, wplt), \(f) f(.x)))


jtls::gwd_clgrph()
dt_nbrs <- gd_nbrs()
fwrite(dt_nbrs, paste0(c_dirs$misc, "nbrs.csv"), quote = F)



    
gp_pminfl_cbn(l_mdls$r_pes2, l_mdls$r_auc2)
# huh wtf is happening with that PES coef.. really hard to believe


