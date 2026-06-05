library(PanelMatch)
library(data.table)
library(fixest)
library(purrr)
library(jtls)
library(broom)



prep_pm_long <- function(pm_object, dt_panel) {
    # pm_object$att[1] %>% names
    
    ## process the sets into one df    
    dt_sets <- keep(pm_object$att, ~len(.x) > 0) %>%
        imap(~attr(.x, "weights") %>% vadt(name = "wbcode2", vlu = "weight") %>%
                 .[, type := "crol"] %>%
                 rbind(data.table(wbcode2 = strsplit(.y, "\\.")[[1]][1], weight = 1, type = "treat")) %>%
                 .[, `:=`(year = strsplit(.y, "\\.")[[1]][2], set = .y)]) %>% rbindlist %>%
        .[, `:=`(wbcode2 = as.integer(wbcode2), year = as.integer(year))]

    l_leads <- attr(pm_object,"lead")
    dv <- attr(pm_object, "outcome")
    l_y_lead <- paste0(dv, "_l" ,l_leads)
    l_y_diff <- paste0(dv, "_dil", l_leads)

    ## construct leading DVs (long-term impact)
    dt_yl <- dt_panel %>% copy %>%
        .[, (l_y_lead) := map(l_leads, ~shift(get(dv), -.x)), wbcode2] %>%
        .[, y_lag := shift(get(dv), 1), wbcode2] %>% 
        .[, (l_y_diff) := map(l_y_lead, ~get(.x) - y_lag)]
        
    ## dt_yl %>% na.omit
    ## dt_yl[, .N, .(wbcode2)]

    dt_sets_wcov <- merge(dt_sets, dt_yl, by= c("wbcode2", "year")) %>%
        .[, setlen := .N, set]

    return(dt_sets_wcov)

}




## most simple OLS: for contemporaneous effect essentially

## wrong approach: uses y
## feols(y  ~ type | set, dt_pm_long, weights = ~weight) # -1.29
## rx <- feols(c(y_l0, y_l1, y_l2, y_l3, y_l4) ~ type | set, dt_pm_long, weights = ~weight) # -1.29 wrong


    


## ## try manual difference
## dt_diff <- dt_pm_long[weight != 0][order(set, -weight)] %>%
##     .[type == "crol", .(y_crol = mean(y-y_lag)), set] %>%
##     merge(dt_pm_long[type == "treat", .(wbcode2, set, y_treat = y - y_lag)], by = "set") %>%
##     .[, y_diff := y_treat - y_crol]


## dt_diff[, mean(y_diff)] # also works



## ## try boostrap
## r_boot <- boot(dt_diff$y_diff, \(data, indices) mean(data[indices]), R = 1000)
## r_boot$t %>% mean # still -1.21

## ## compare estimates: 
## rbind(data.table(vlu = pe_results$bootstrapped.estimates[,1], src = "PM"),
##       data.table(vlu = r_boot$t[,1], src = "mnl")) %>%
##     ggplot(aes(x=vlu, color = src)) + geom_density()
## ## mnl is much wider
    

  
  
## * main

## ** PanelMatch

data(dem)
dem <- as.data.frame(lapply(dem, function(x) if(is.list(x)) unlist(x) else x))

pm_data <- PanelData(panel.data = dem,
                     unit.id = "wbcode2",
                     time.id = "year",
                     treatment = "dem",
                     outcome = "y")

pm_results <- PanelMatch(panel.data = pm_data, lag = 4, 
                         refinement.method = "ps.match", 
                         match.missing = TRUE, 
                         covs.formula = ~tradewb + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                         size.match = 5, qoi = "att",
                         lead = 0:4, 
                         forbid.treatment.reversal = FALSE)

pe_results <- PanelEstimate(sets = pm_results, panel.data = pm_data, se.method = "bootstrap")




## ** custom right approach: uses y diff (ydi)

dt_pm_long <- prep_pm_long(pm_results, adt(dem))


## get treatment indicator
## dt_pm_long[, c("set_cry", "set_year") := tstrsplit(set, "\\.", type.convert = as.integer)]

## dt_pm_long[year == set_year & wbcode2 == set_cry, Dit := 1] %>%
##     setnafill(fill = 0 , cols = "Dit")

## dt_pm_long[order(set), .(wbcode2, year, set, set_year, set_cry, Dit)] %>% print(n=200)



feols(y - y_lag ~ type | set, dt_pm_long, weights = ~weight) # -0.6
rx <- feols(c(y_dil0, y_dil1, y_dil2, y_dil3, y_dil4) ~ type | set, dt_pm_long, weights = ~weight) # DV has to be diff

dt_mnl <- map(rx, ~tidy(.x) %>% adt) %>% rbindlist(idcol = "spec") %>%
    .[, year := as.integer(gsub("lhs: y_dil", "", spec))]


## ** comparison

rbind(
    pe_results %>% summary %>% adt %>% .[, .(estimate, std.error, year = 0:(.N-1), src = "PM")],
    dt_mnl[, .(estimate, std.error, year, src = "MNL")]) %>%
    dcast(year ~ src, value.var = c("estimate", "std.error"))



## ** expansion: decade
dt_pm_long[, decade := 10*floor(year/10)][, decade := factor(decade, levels = sort(unique(decade)))]

r_dcd <- feols(c(y_dil0, y_dil1, y_dil2, y_dil3, y_dil4) ~ type*decade | set, dt_pm_long, weights = ~weight, cluster = ~wbcode2)

## feols(c(y_dil0, y_dil1, y_dil2, y_dil3, y_dil4) ~ type | set + wbcode2, dt_pm_long, weights = ~weight)
## gemini says DiD already accounts for differences from means


dt_dcd <- map(r_dcd, ~tidy(.x) %>% adt) %>% rbindlist(idcol = "spec") %>%
    .[, `:=`(year = gsub("lhs: y_dil", "", spec), decade = gsub("typetreat:decade", "", term))] %>%
    .[decade == "typetreat", decade := levels(dt_pm_long$decade)[1]] %>%
    rbind(copy(dt_mnl)[, `:=`(decade = "all")])


map(r_dcd, ~gd_intr_grp_slopes(.x, "typetreat", "decade")) %>% rbindlist(idcol = "spec") %>%
    .[, `:=`(year = gsub("lhs: y_dil", "", spec))] %>%
    rbind(copy(dt_mnl)[, `:=`(category = "all",p.value = NULL, statistic = NULL)] %>%
          setnames(old = c("std.error", "term"), new = c("se", "vrbl"))) %>%
    ggplot(aes(x=year, y = estimate, group = category, color = category,
               ymin = estimate - 1.96*se, ymax = estimate + 1.96*se, fill = category)) +
    geom_line() +
    geom_ribbon(alpha = 0.2)
    

dt_pm_long[type == "treat", .N, decade]

dt_dcd %>% ggplot(aes(x=year, y = estimate, group = decade, color = decade)) + geom_line()




DisplayTreatment(panel.data = dem.panel, legend.position = "none",
                 xlab = "year", ylab = "Country Code")
