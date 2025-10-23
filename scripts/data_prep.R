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


gd_exhbppl <- function(dt_af_people, dt_af_exhbs, dt_af_instns_rcd, dt_ls_af_aucres, dt_af_links_exhb_ppl) {
    #' 

    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    gw_fargs(match.call())
    

    ## merge people info to exhibition -> organization
    ## dt_af_instns_rcd <- md_instns_rcd(c_datargs$clsfcn_to_use, dt_af_instns, b_rcd_to_pm = c_datargs$b_rcd_to_pm,
    ##                                   dt_af_pmdb_matches) 

    
    ## transform nationality, region and Subregion to lists
    ## tstrsplit and merging on country (not person) to minimize calls to countrycode and list

    dt_natlist <- dt_af_people[, .(Nationality)] %>% funique %>% 
        .[, .(Nationality, nat = trimws(unlist(tstrsplit(Nationality, ","))))] %>% na.omit %>% # somehow tstrsplit adds NAs -> yeet
        .[, .(Nat_list = list(nat)), .(Nationality)]


    ## merge back to people
    dt_af_people_natlist <- merge(dt_af_people, dt_natlist, by = "Nationality")


    dt_exhbppl_prep <- merge(dt_af_exhbs, dt_af_instns_rcd[, .(InstitutionID, InstitutionType, Country)],
                             by = "InstitutionID") %>%
        .[, .(InstitutionID, ExhibitionID, begin_year, BeginDate, EndDate, InstitutionType, Country)] %>% 
        merge(dt_af_links_exhb_ppl, ., by = "ExhibitionID") %>%
        merge(dt_af_people_natlist, by = "PeopleID")

    ## get first year and age, filter down by yeeting babies

    dt_exhbppl_prep2 <- dt_exhbppl_prep[, first_year := min(begin_year), PeopleID] %>%
        .[, age := begin_year - BirthYear] %>%
        .[, gender_rcd := fifelse(Gender %in% c("Male", "Female"), Gender, "Other")] %>% 
        .[, career_year := begin_year - first_year]

    ## %>%
    ##     .[eval(c_filters_atyear$yeet_babies_completely)] 


    c_crycode_reg <- c("Kosovo" = "Europe", "Micronesia" = "Oceania", "Taiwan" = "Asia",
                       "Netherlands Antilles" = "Americas")
    c_crycode_subreg <- c("Kosovo" = "Southern Europe", "Micronesia" = "Melanesia", "Taiwan" = "Eastern Asia",
                          "Netherlands Antilles" = "Latin America and the Caribbean")

    ## whether by exhibition country X artists is national/subregional/regional;
    ## need People ID because link to person required
    dt_insat_nats <- dt_exhbppl_prep2[, .(artist_nat = unlist(Nat_list)), .(Country, PeopleID)] %>%
        .[, `:=`( Subregion = countrycode(Country, "country.name", "un.regionsub.name",
                                          custom_match = c_crycode_subreg),
                 Region = countrycode(Country, "country.name", "un.region.name",
                                      custom_match = c_crycode_reg), 
                 artist_subregion = countrycode(artist_nat, "country.name", "un.regionsub.name",
                                                custom_match = c_crycode_subreg),
                 artist_region = countrycode(artist_nat, 'country.name', 'un.region.name',
                                             custom_match = c_crycode_reg))] %>%
        .[, .(exhb_nat = any(Country == artist_nat), exhb_subreg = any(Subregion == artist_subregion),
              exhb_reg = any(Region == artist_region)), .(PeopleID, Country)]

    dt_exhbppl_prep3 <- merge(dt_exhbppl_prep2, dt_insat_nats, by = c("PeopleID", "Country")) %>%
        .[order(PeopleID, BeginDate)]

    ## add cumulative auction data
    dt_aucres <- dt_ls_af_aucres[order(AF_PID, year)] %>% 
        .[, .(aucturn = count*price,
              begin_year = year+1,
              aucturn_cum = fcumsum(as.numeric(count*price)),
              aucturn_bin = 1), .(PeopleID = AF_PID)] %>%
        .[, `:=`(aucturn_log = log1p(aucturn),
                 aucturn_cum_log = log1p(aucturn_cum))]
    
    dt_aucres_rolled <- dt_aucres[, .(PeopleID, begin_year, aucturn_cum, aucturn_cum_log)] %>%
        .[dt_exhbppl_prep3, on = .(PeopleID, begin_year), roll = Inf, rollend = c(F,T)] %>%
        setnafill(fill = 0, cols = c("aucturn_cum", "aucturn_cum_log"))

    ## dt_aucres[PeopleID == 44]
    ## dt_aucres_rolled[PeopleID == 44 & begin_year > 2002, c(list(PeopleID, begin_year), .SD), .SDcols = patterns("^auc*")] %>% print(n=80)
    ## dt_aucres_rolled[PeopleID == 44, c(list(PeopleID, begin_year), .SD), .SDcols = patterns("^auc*")] %>% print(n=80)

    ## set up and roll auction turnover (lagged)
    
        
    ## find people with gaps in exhb history but auction: is it carried over properly? 
    ## dt_aucres_rolled[begin_year > 2010, .(max_gap = max(begin_year - shift(begin_year), na.rm = T)*1.0), PeopleID] %>%
    ##     .[dt_aucres[, .(PeopleID = unique(PeopleID))], on = "PeopleID"]
        
    
    ## add auction data of years (not cumulative)
    dt_aucres_prepl0 <- dt_ls_af_aucres %>% copy %>%
        .[, auc_turnover_l0 := count*price] %>% .[order(AF_PID, year)] %>%
        .[, .(PeopleID = AF_PID, begin_year = year, auc_turnover_l0)]
    
    ## "lagging" auction values by adding one year (no need for grid)
    dt_aucres_prep_l1 <- dt_aucres_prepl0 %>% copy %>% .[, begin_year := begin_year+1] %>% 
        setnames(old = "auc_turnover_l0", new = "auc_turnover_l1")
    
    ## combine auction success data, calculate growth: 
    ## how much cum auc has changed from end of prev year to end of prevprev year
    dt_exhbppl_prep4 <- merge(dt_aucres_rolled, dt_aucres_prepl0, by = c("PeopleID", "begin_year"), all.x = T) %>%
        merge(dt_aucres_prep_l1, by = c("PeopleID", "begin_year"), all.x = T) %>%
        setnafill(fill = 0 , cols = c("auc_turnover_l0", "auc_turnover_l1")) %>%
        .[, aucturn_cum_grow := 100*((aucturn_cum/(aucturn_cum - auc_turnover_l1))-1)] %>%
        .[is.nan(aucturn_cum_grow), aucturn_cum_grow := 0]
                
    
    ## manually inspecting auction numbers: 29, 44, 31
    ## dt_exhbppl_prep4[PeopleID == 31, c(list(PeopleID=PeopleID, begin_year = begin_year), .SD), .SDcols = patterns("auc*")] %>%
    ##     .[, head(.SD, 1), begin_year] %>% print(n=800)
    ## dt_aucres[PeopleID == 31]
    ## dt_aucres_prepl0[PeopleID == 31] # 106k is in 2011: 2012 should be some growth, 2013 should be no growth (because 2012 has no real growth)


    ## dt_exhbppl_prep4[, .N, aucturn_cum_grow]

    ## dt_exhbppl_prep4[, .(mean_grow = mean(aucturn_cum_grow[aucturn_cum_grow != Inf], na.rm = T),
    ##                      mean_inf = mean(aucturn_cum_grow == Inf, na.rm = T)) , InstitutionType]

    ## dt_exhbppl_prep4[between(aucturn_cum_grow, 1, 2, incbounds = F)] %>% ggplot(aes(x=aucturn_cum_grow, color = InstitutionType)) + geom_density()

    

    ## dt_exhbppl_prep4[PeopleID == 958, .(PeopleID, begin_year, ExhibitionID, auc_turnover_l0, auc_turnover_l1, aucturn_cum)] %>% print(n=80)

    ## create dummies for institution types with dcast
    dt_instype_dummies <- dt_exhbppl_prep4 %>% copy %>% .[, `:=`(x=1, nbr = 1:.N)] %>%
        dcast(nbr + PeopleID ~ paste0("instype_", InstitutionType), value.var = "x", fill = 0) %>%
        .[, nbr := NULL] %>%
        .[, fcumsum(.SD, PeopleID), .SDcols = patterns("^instype")] %>% # use fcumsum for aggregation
        .[, cnt_prevexhb := Reduce(`+`, .SD)-1] # calculate total

    
    ## merge institution type dummies to main data (drop Nationality list, otherwise super slow)
    dt_exhbppl <- cbind(dt_exhbppl_prep4[, Nat_list := NULL], dt_instype_dummies) %>%
        .[begin_year <= END_YEAR]

    ## dt_exhbppl[PeopleID == 879, .(cnt_prevexhb)]
    
    attr(dt_exhbppl, "gnrtdby") <- as.character(match.call()[[1]])

    return(dt_exhbppl)
}
