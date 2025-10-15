

## ** for now, use dt_atyear from infl
## FIXME: make independent from infl









# add some code

c_filters_atyear <- list(
    keep_all = expr(T),
    yeet_baby_exhbs = expr(begin_year > BirthYear + 16),
    yeet_babies_completely = expr(first_year > BirthYear + 16))


dt_atyear_foster <- gd_atyear_foster(dt_atyear_raw, c_filters_atyear$yeet_babies_completely)

dt_exhbppl_sample <- dt_exhbppl[sample(1:.N, size = 1e5)]




dt_exhbppl_watyear <- merge(
    dt_exhbppl,
    dt_atyear_foster[, c(list(PeopleID = PeopleID, begin_year = year), .SD), .SDcols = patterns("exhbN")],
    by = c("PeopleID", "begin_year")) %>%
    .[order(PeopleID, begin_year)]

## minor differences between cnt_prevexhb and exhbN_all_cum: exhbN_all_cum is only previous year?
## hmm not quite: diffs can also be negative -> FIXME: make sure there's no ambiguity

dt_exhbppl_watyear[, .(PeopleID, begin_year, cnt_prevexhb, exhbN_all_cum, diff_cum = cnt_prevexhb - exhbN_all_cum)] %>%
    .[, .N, diff_cum] %>% .[order(diff_cum)] %>% print(n=80)

## now can calculate institution-specific effect: total count, relative
dt_exhbppl_watyear[, `:=`(exhbN_grow5_cnt = exhbN_all_cum_n5 - cnt_prevexhb, # how many new exhibitions are in the next 5 years
                          exhbN_grow5_perc =  100*((exhbN_all_cum_n5/cnt_prevexhb)-1) # by how much does the total exhibition count grow
                          )]
dt_exhbppl_watyear[PeopleID ==44 & begin_year > 2005, .(PeopleID, begin_year, exhbN_all, exhbN_all_cum_n5, cnt_prevexhb,
                                                        exhbN_grow5_cnt, exhbN_grow5_perc)] %>% print(n=40)
    

dt_exhbppl_watyear[, .(mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm = T))]

dt_exhbppl_watyear[, .(mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm = T),
                       mean_grow5_perc = mean(exhbN_grow5_perc[!is.infinite(exhbN_grow5_perc)], na.rm = T)), InstitutionType]


dt_exhbppl_watyear[is.infinite(exhbN_grow5_perc), uniqueN(PeopleID)]

dt_atyear_foster

## distribution of institution-level growth factors (ability of individual museums to promote)
dt_insgrow5 <- dt_exhbppl_watyear[, .(mean_grow5_perc = mean(exhbN_grow5_perc[!is.infinite(exhbN_grow5_perc)], na.rm =T),
                                      mean_grow5_cnt = mean(exhbN_grow5_cnt, na.rm =T)),
                                  .(InstitutionID, InstitutionType)]

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



