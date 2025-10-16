gd_atyear_foster <- function(dt_atyear_raw, c_filters_atyear) {

    gw_fargs(match.call())
    
    ## apply filter
    dt_atyear_raw_fltrd <- dt_atyear_raw[eval(c_filters_atyear)] %>%
        .[, `:=`(age_cut = factor(as.integer(floor(age/5))*5))]
        

    
    l_cols_exhbN <- keep(names(dt_atyear_raw_fltrd),  ~grepl("exhbN_", .x))

    l_cols_exhbN_s5 <- paste0(l_cols_exhbN, "_s5") # previous
    l_cols_exhbN_n5 <- paste0(l_cols_exhbN, "_n5") # next
    
    l_cols_exhbN_cum <- paste0(l_cols_exhbN, "_cum")

    dt_atyear_lagprep1 <- dt_atyear_raw_fltrd %>% copy %>% # .[PeopleID == 44] %>% 
        .[, (l_cols_exhbN_s5) := frollsum(.SD, n=6, fill = NA), PeopleID, .SDcols = l_cols_exhbN] %>%
        .[, (l_cols_exhbN_s5) := map(l_cols_exhbN, ~.SD[, get(paste0(.x, "_s5")) - get(.x)])] %>%
        .[, (l_cols_exhbN_n5) := frollsum(.SD, n = 6, fill = NA, align = "left"), PeopleID,.SDcols=l_cols_exhbN] %>%
        .[, (l_cols_exhbN_n5) := map(l_cols_exhbN, ~.SD[, get(paste0(.x, "_n5")) - get(.x)])] %>%
        .[, (l_cols_exhbN_cum) := fcumsum(.SD), PeopleID, .SDcols = l_cols_exhbN] %>%
        .[, (l_cols_exhbN_cum) := map(l_cols_exhbN, ~.SD[, get(paste0(.x, "_cum")) - get(.x)])]

    dt_atyear_lagprep1[PeopleID ==44, .(PeopleID, age, year, exhbN_all, exhbN_all_s5, exhbN_all_n5, exhbN_all_cum)]
    
    ## get cumsum of future by shifting years
    ## i think shifting by 5 is correct choice, but not 100% sure, maybe also 6?
    dt_cumsum_n5 <- dt_atyear_lagprep1[, c(list(PeopleID = PeopleID, year = year-5), .SD), .SDcols = l_cols_exhbN_cum] %>%
        setnames(old = l_cols_exhbN_cum, new = paste0(l_cols_exhbN_cum, "_n5"))

    dt_atyear_lagprep2 <- merge(dt_atyear_lagprep1, dt_cumsum_n5, by = c("PeopleID", "year"), all.x = T)
    

    dt_atyear_lagprep2[PeopleID ==44, .(PeopleID, age, year, exhbN_all, exhbN_all_s5, exhbN_all_n5, exhbN_all_cum, exhbN_all_cum_n5)] %>% print(n=80)

    attr(dt_atyear_lagprep2, "gnrtdby") <- as.character(match.call()[[1]])

    return(dt_atyear_lagprep2)

}




gd_atyear_rcd <- function(dt_af_exhbs, dt_af_links_exhb_ppl, dt_af_instns_rcd) {
    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    ## expand to artist year based on classification

    #' @param dt_af_exhbs AF exhibition information
    #' @param dt_af_links_exhb_ppl AF links exhibition-people
    #' @param dt_af_instns_rcd data.table of how to aggregately recode AF institutions 
                              
    
    ## first get recoded exhibitions, then link that to people to get atyear with recoded InstitutionType
    dt_atyear_rcd <- merge(dt_af_exhbs, dt_af_instns_rcd[, .(InstitutionID, InstitutionType)],
                           by = "InstitutionID") %>%
        .[, .(ExhibitionID, begin_year, InstitutionType)] %>% 
        merge(dt_af_links_exhb_ppl, ., by = "ExhibitionID") 

    dt_af_inst_clsfcn <- gd_af_inst_clsfcn()
        
    ## cast to wide for merging
    ## if I use different classifications: add indicator to new names here
    dt_atyear_rcd_wide <- dt_atyear_rcd[, .N, .(year = begin_year, PeopleID, InstitutionType)] %>%
        dcast(PeopleID + year ~ InstitutionType, value.var = "N") %>%
        setnames(old = dt_atyear_rcd[, unique(InstitutionType)],
                 new = paste0("exhbN_", dt_atyear_rcd[, unique(InstitutionType)]))

    return(dt_atyear_rcd_wide)
}

gd_atyear_raw <- function(dt_af_exhbs, dt_af_people, dt_af_links_exhb_ppl, dt_af_instns_rcd, dt_ls_af_aucres) {
    #' generate artist-year data.table
    #' @param dt_af_exhbs AF exhibition info
    #' @param dt_af_people AF people info
    #' @param dt_af_links_exhb_ppl AF links exhibition-people
    #' @param dt_af_instns_rcd recoded institutions

    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;
    gw_fargs(match.call())
    

    ## merge exhibition and people
    dt_cbn <- merge(dt_af_links_exhb_ppl, dt_af_people, by = "PeopleID") %>% 
        merge(dt_af_exhbs, by = "ExhibitionID")


    ## get minimum year
    dt_atyear_start <- dt_cbn[, .(first_year = min(begin_year)), .(PeopleID, BirthYear, Nationality, Gender)]
    

    ## expand to grid
    dt_atyear_grid <- dt_atyear_start[!is.na(first_year), .(year = first_year:END_YEAR),
                                      .(PeopleID, BirthYear, Nationality, Gender, first_year)]
    
    ## calculate total exhibition count
    dt_atyear_exhbcnt <- dt_cbn[, .(exhbN_all = .N), .(year = begin_year, PeopleID)]

    ## get the exhibition counts grouped by (reclassified) institutional form
    ## if I want multiple classifications, I can easily get them here
    dt_atyear_rcd_wide <- gd_atyear_rcd(dt_af_exhbs, dt_af_links_exhb_ppl, dt_af_instns_rcd)
        
    ## merge all the exhibition counts together
    dt_cbn_prep1 <- Reduce(\(x,y) merge(x,y, by = c("year", "PeopleID"), all.x = T),
                           list(dt_atyear_grid, dt_atyear_exhbcnt, dt_atyear_rcd_wide))
    
    ## set up exhibition names to refer to them later on
    l_cols_exhbN <- keep(names(dt_cbn_prep1), ~grepl("exhbN_", .x))

    ## fill up no exhibition with zero, calculate last exhibition year
    dt_atyear_raw <- copy(dt_cbn_prep1) %>%
        setnafill("const", cols = l_cols_exhbN, fill = 0) %>%
        .[, `:=`(year_last_exhb = max(year[which(exhbN_all > 0)])), PeopleID] %>% 
        .[, `:=`(last_exhb = (year == year_last_exhb)*1)] %>% # last year
        .[, `:=`(age = year - BirthYear, career_year = year - first_year)] %>% # other career stage data
        .[order(PeopleID, year)]

    ## add ls auction data
    dt_atyear_wls <- merge(dt_atyear_raw, dt_ls_af_aucres,
                           by.x = c("PeopleID", "year"), by.y=  c("AF_PID", "year"), all.x = T) %>%
        setnafill("const", cols = c("count", "price"), fill = 0) %>%
        .[, auc_turnover := count*price]

    
    attr(dt_atyear_wls, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_atyear_wls)
}


gd_ls_af_aucres <- function(dt_ls_af_aucres, dt_af_ls_matches) {
    #' small helper function to get callgraph working
    gw_fargs(match.call())
    dt_ls_af_aucres <- dt_ls_aucres[dt_af_ls_matches, on = .(url =  ls_url)] %>% na.omit
    attr(dt_ls_af_aucres, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_ls_af_aucres)
}

md_instns_rcd <- function(clsfcn_to_use, dt_af_instns, b_rcd_to_pm, dt_af_pmdb_matches) {
    ## recode institutions according to my classification scheme
    ## separate function to handle it comfily in many cases
    #' @param clsfcn_to_use which classification to use
    #' @param dt_af_instns AF institution information
    #' @param rcd_to_pm whether to recode AF institutions to PM
    #' @param dt_af_pmdb_matches which PMDBs are also in AF, to recode those AF to PM
    gw_fargs(match.call())

    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;

    ## set recoding up in a way that later it can be used with multiple classifications 
    ## if that comes up, I can move it to separate function and map it
    
    # clsfcn_to_use <- "clsfcn2" # later map over this
    dt_af_inst_clsfcn <- gd_af_inst_clsfcn() %>%
        setnames(old = clsfcn_to_use, new = "clsfcn_to_use")
    
    ## recode AF instns with my classification system
    dt_af_instns_rcd <- copy(dt_af_instns) %>%
        .[, InstitutionType_AF := InstitutionType] %>% 
        .[dt_af_inst_clsfcn, InstitutionType := i.clsfcn_to_use, on = "InstitutionType"]

    if (b_rcd_to_pm) {
    
        ## apply PM classification: recode our PMs to pm
        dt_af_instns_rcd_pm <- dt_af_instns_rcd %>% copy %>%
            .[dt_af_pmdb_matches[AF_IID != "nomatch", .(AF_IID = as.integer(AF_IID))],
              InstitutionType := "pm" ,on = .(InstitutionID = AF_IID)]
    } else {
        dt_af_instns_rcd_pm <- copy(dt_af_instns_rcd)
    }

    attr(dt_af_instns_rcd_pm, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_af_instns_rcd_pm)
}
