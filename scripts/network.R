
gd_edges <- function(dt_exhbppl, sample_artists = F) {
    gw_fargs(match.call())

    if (sample_artists > 0) {
        dt_at_sample <- dt_exhbppl[PeopleID %in% sample(unique(PeopleID), size = sample_artists),
                                   .(PeopleID = unique(PeopleID))]
    } else {
        
        dt_at_sample <- dt_exhbppl[, .(PeopleID = unique(PeopleID))]
    }

    dt_exhbppl2 <- dt_exhbppl[dt_at_sample, on = "PeopleID"]
        


    dt_coexhb_prep <- dt_exhbppl2[, .(PeopleID, InstitutionID, BeginDate, BeginDate_chk = BeginDate)]
    dt_coexhb_prep2 <- dt_exhbppl2[, .(PeopleID2 = PeopleID, InstitutionID2 = InstitutionID,
                                      BeginDate2_chk = BeginDate, # some check columns
                                      BeginDate2 = BeginDate)]

    dt_edges <- dt_coexhb_prep[dt_coexhb_prep2,
                               on = .(PeopleID = PeopleID2, BeginDate < BeginDate2), nomatch = NULL] %>%
        .[, tdiff := BeginDate2_chk - BeginDate_chk]

    attr(dt_edges, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_edges)
}


gd_edges_wtd <- function(dt_edges) {
    gw_fargs(match.call())
    
    ## dt_edges[, .N, tdiff]
    ## hmm seems the non-equi join messes up original dates -> use chk columns that aren't used for join for time diff

    ## could implement different tau parameters here, also different years

    dt_edges_wtd <- dt_edges[, .(weight = .N), .(InstitutionID, InstitutionID2)]
    ## 23m unique links, 35k nodes
    attr(dt_edges_wtd, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_edges_wtd)
}


gg_insflow <- function(dt_edges_wtd) {
    gw_fargs(match.call())
    g_insflow <- graph_from_data_frame(dt_edges_wtd, directed = TRUE)

    attr(g_insflow, "gnrtdby") <- as.character(match.call()[[1]])
    return(g_insflow)
    
}

gd_insevc <- function(g_insflow, dt_af_instns_rcd) {
    #' create data.table of institution eigenvector centralities (EVC)
    gw_fargs(match.call())
    
    
    ## EVC; uses weights by default
    eigen_centrality_values <- eigen_centrality(g_insflow, directed = TRUE)
    

    ## Create a data frame of institutions' eigenvector centrality values
    dt_insevc_prep <- data.table(InstitutionID = as.integer(names(eigen_centrality_values$vector)),
                            evc = eigen_centrality_values$vector)

    ## dt_af_instns_rcd <- md_instns_rcd(c_datargs$clsfcn_to_use, dt_af_instns,
    ##                                   c_datargs$b_rcd_to_pm, dt_af_pmdb_matches)

    dt_insevc <- merge(dt_insevc_prep, dt_af_instns_rcd, by = "InstitutionID")
    

    attr(dt_insevc, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_insevc)
}


gd_deg_info <- function(gx) {

    l_degtype <- c("out", "in", "all")

    ## m_deg <- map(l_degtype, ~degree(gx, mode = .x)) %>% Reduce(cbind, .)
    m_deg <- map(l_degtype, ~cbind(degree(gx, mode = .x), strength(gx, mode = .x))) %>% Reduce(cbind, .)
    
    ## colnames(m_deg) <- paste0("deg_", l_degtype)
    colnames(m_deg) <- outer(c("deg_", "str_"), l_degtype, FUN = paste, sep =  "") %>% as.vector
    
    

    dt_deg_info <- adt(m_deg, keep.rownames = "InstitutionID") %>%
        .[, InstitutionID := as.numeric(InstitutionID)]


    return(dt_deg_info)


}


gt_deg_info <- function(g_insflow, dt_af_instns_rcd) {
    gw_fargs(match.call())

    merge(dt_af_instns_rcd[, .(InstitutionID, InstitutionType)],
          gd_deg_info(g_insflow), by = "InstitutionID") %>%
        .[, map(.SD, mean), InstitutionType, .SDcols = patterns("^deg|^str")] %>% gt_plain
}

gd_exhbppl_evc <- function(dt_exhbppl, dt_insevc) {
    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    gw_fargs(match.call())

    ## apply back to people
    dt_exhbppl_evc <- merge(dt_exhbppl, dt_insevc[, .(InstitutionID, evc)],
                            by = "InstitutionID", all.x = T) %>%
        setnafill(fill = 0, cols = "evc") %>%  # make it work in case instns are sampled
        .[, exhb_nbr_artists := .N, ExhibitionID] %>%
        .[, evc_split := evc/exhb_nbr_artists] # 'reward' of institution gets split by person

    
    attr(dt_exhbppl_evc, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_exhbppl_evc)
}


gd_atyear_evc <- function(dt_exhbppl_evc, dt_atyear) {
    if (as.character(match.call()[[1]]) %in% fstd){browser()}
    gw_fargs(match.call())

    ## i want lagged eigv ctrlt: how high are eigenvector centrality when exhibited in PM?
    ## but this needs to be aggregated by year, also maybe some average/sum of past N years? 
    ## link between InstitutionType and Exhibition gets severed, but that sounds kinda realistic:
    ## curators probably have a general picture of how successful an artist is, even without knowing every exhibition
    ## also not just previous exhibition taken into consideration, but overall exhibition profile
    ## dt_exhbppl_evc[, .(mean(evc_split))]

    ## FIXME: evc_ay should be evc_ay_ttl
    dt_atyear_evc_prep <- dt_exhbppl_evc[, .(evc_ay_ttl = sum(evc_split)),
                                         .(PeopleID, year = begin_year)] %>% .[order(PeopleID, year)]

    

    ## apply to artist year
    ## have to do the frollsum for backward and forward separately, apparently can't use `:=` with frollsum
    dt_atyear_evc_prep2 <- merge(dt_atyear, dt_atyear_evc_prep, by = c("PeopleID", "year"), all.x = T) %>%
        setnafill(cols = "evc_ay_ttl", fill = 0) %>%
        .[, evc_s5 := frollsum(.SD, n=6), PeopleID, .SDcols = "evc_ay_ttl"] %>%
        .[, evc_n5 := frollsum(.SD, n=6, align = "left"), PeopleID, .SDcols = "evc_ay_ttl"] %>%
        .[, auc_s5 := frollsum(.SD, n=6), PeopleID, .SDcols = "auc_turnover"] %>%
        .[, auc_n5 := frollsum(.SD, n=6, align = "left"), PeopleID, .SDcols = "auc_turnover"] %>%
        .[, `:=`(evc_s5 = evc_s5 - evc_ay_ttl, evc_n5 = evc_n5 - evc_ay_ttl,
                 auc_s5 = auc_s5 - auc_turnover, auc_n5 = auc_n5 - auc_turnover)]
    
    ## aggregate EVC by institutiontype and year
    dt_atyear_evc_prep3 <- dt_exhbppl_evc[, .(evc = sum(evc_split)),
                                          .(PeopleID, begin_year,
                                            InstitutionType = paste0("evc_ay_", InstitutionType))] %>%
        dcast(PeopleID + begin_year ~ InstitutionType, value.var = "evc")
    
    ## combine s5/n5 with evc-institutiontype-atyear
    dt_atyear_evc <- merge(dt_atyear_evc_prep2, dt_atyear_evc_prep3,
                           by.x = c("PeopleID", "year"), by.y = c("PeopleID", "begin_year"), all.x = T) %>%
        setnafill(fill = 0, cols = keep(names(.), ~grepl("evc_ay_", .x)))


    ## dt_atyear_evc[, .(evc_ay_ttl)] %>% summary
    ## dt_atyear_evc[PeopleID == 11, .(year, evc_ay_tll, evc_s5, evc_n5)] %>% 

    attr(dt_atyear_evc, "gnrtdby") <- as.character(match.call()[[1]])
    return(dt_atyear_evc)
}
