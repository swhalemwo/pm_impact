
dt_exhbppl_wyeet <- gd_exhbppl(dt_af_people, dt_af_exhbs, dt_af_instns_rcd, dt_ls_af_aucres, dt_af_links_exhb_ppl)

dt_exhbppl <- dt_exhbppl_wyeet[InstitutionType != "yeet"] %>%
    .[eval(c_datargs$baby)]

# dt_exhbppl_sample <- dt_exhbppl[sample(1:.N, size = 1e5)]

## hacky way of integrating EVC, if decided to keep integrate better



dt_edges <- gd_edges(dt_exhbppl, sample_artists = 100)
dt_edges_wtd <- gd_edges_wtd(dt_edges) %>% .[InstitutionID != InstitutionID2]
g_insflow <- gg_insflow(dt_edges_wtd)


dt_insevc <- gd_insevc(g_insflow, dt_af_instns_rcd)
dt_insevc[order(-evc), .(Name, evc, InstitutionType)] # moma, Zentrum für Kunst und Medien Karlsruhe, vienna huh

dt_exhbppl_evc <- gd_exhbppl_evc(dt_exhbppl, dt_insevc)

dt_atyear_evc <- gd_atyear_evc(dt_exhbppl_evc, dt_atyear_foster)


## check EVC developments: looks ok? 
dt_atyear_evc[, c(list(PeopleID = PeopleID, year = year), .SD), .SDcols = patterns("^evc")] %>%
    print(n=50)


dt_exhbppl_watyear <- gd_exhbppl_watyear(dt_exhbppl, dt_atyear_evc)

## "evc_grow5_cnt", "evc_grow5_prop"
l_cols_summary <- c("exn_all_cum", "exn_all_s5", "exn_grow5_cnt", "exn_grow5_perc", "exn_perc_n5s5", "exn_perc_n5s5_org")


dt_insgrow5 <- gd_insgrow5(dt_exhbppl_watyear, l_cols_summary)


dt_atyear_evc_smol <- dt_atyear_evc[year %between% c(2010, 2012)]
