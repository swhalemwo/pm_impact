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
library(zeallot)   # destructuring assingment, comfy but not very necessary
library(stringr)   # string manip? i guess
library(ggrepel)   # plotting
library(sf) # locating toti
library(rnaturalearth) # locating toti
library(readr, include.only = "read_csv") # read google sheets
library(R.utils, include.only = "withTimeout")
set.seed(42)

refresh_pipeline <- function() {
    invisible(lapply(l_funcfiles, \(x) source(x, echo = F, verbose = F, print.eval = F)))
    message("Pipeline refreshed.")
}

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


if (interactive()) {
    b_min_mdls <- T
    
    ## agent integration for interactive use
    library(mcptools) 
    library(btw)
    options(btw.run_r.enabled = TRUE)

    ## mcptools::mcp_server(type = "http", port = 8080)
    ## Sys.setenv(MCPTOOLS_PORT = "8080")
    ## mcptools::mcp_session()

    x <- mcptools::mcp_session()

} else {
    b_min_mdls <- F
}

b_min_mdls <- T

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
    clsfcn_to_use = "clsfcn5",
    b_rcd_to_pm = T,
    b_sample_dt_atyear = T,
    c_toti = "TIV_10",
    baby = list(
        keep_all = expr(T),
        yeet_baby_exhbs = expr(begin_year > BirthYear + 16),
        yeet_babies_completely = expr(first_year > BirthYear + 16)) %>%
        chuck("yeet_babies_completely"))




## auction data
dt_ls_aucres <- gd_ls_aucres() %>% adt
dt_af_ls_matches <- gd_af_ls_matches()
dt_ls_af_aucres <- gd_ls_af_aucres(dt_ls_aucres, dt_af_ls_matches)

## get the art newspaper museum year data with artfacts id
dt_tanp_muyr_afid <- gd_tanp_muyr_afid()

dt_af_instns_rcd <- md_instns_rcd(c_datargs$clsfcn_to_use, dt_af_instns,
                                  c_datargs$b_rcd_to_pm, dt_af_pmdb_matches) %>%
    gd_instype_toti(dt_tanp_muyr_afid, dt_af_pmdb_matches, c_datargs$c_toti)

## --------- artist-year data
dt_atyear_raw <- gd_atyear_raw(dt_af_exhbs, dt_af_people, dt_af_links_exhb_ppl, dt_af_instns_rcd, dt_ls_af_aucres)


## calculate some before/after steps
dt_atyear_foster <- gd_atyear_foster(dt_atyear_raw, c_datargs$baby)


if (b_min_mdls) { 

    dt_atyear_foster_smol <- dt_atyear_foster[year %between% c(2010, 2012)]
    l_mdls <- gl_mdls(dt_atyear_foster_smol)
} else {
    l_mdls <- gl_mdls(dt_atyear_foster) # , dt_atyear_evc)
}


## dt_atyear_evc[, map(.SD, ~sum(is.na(.x)))] %>% t %>% adt %>% .[, .N, V1]





## super weird EVC values: maybe extreme cases?
## dt_exhbppl_watyear %>% ggplot(aes(x=evc_grow5_prop, color = InstitutionType)) + geom_density() +
##     scale_x_continuous(trans = scales::log10_trans()) +
##     coord_cartesian(xlim = c(1e-3,1e3))


## dt_exhbppl_watyear[is.infinite(exn_grow5_perc), uniqueN(PeopleID)]


## distribution of institution-level growth factors (ability of individual museums to promote)

if (!b_min_mdls) {    
    walk(names(gc_tbls()), ~lapply(c(gtbl, wtbl), \(f) f(.x))) # currently all use grow/insflow objects
}

walk(names(gc_plts()), ~lapply(c(gplt, wplt), \(f) f(.x)))


jtls::gwd_clgrph()
dt_nbrs <- gd_nbrs()
fwrite(dt_nbrs, paste0(c_dirs$misc, "nbrs.csv"), quote = F)



    
## gp_pminfl_cbn(l_mdls$r_pes2, l_mdls$r_auc2)
## huh wtf is happening with that PES coef.. really hard to believe






