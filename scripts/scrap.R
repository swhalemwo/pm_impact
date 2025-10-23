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
nlist(

dt_atyear_evc[, c(nlist(PeopleID, year), .SD) , .SDcols = patterns("^evc")]


## reversed graph: this would basically be who shows the most established art?
## kinda going back in time? 
g_insflow2 <- dt_edges_wtd[, .(InstitutionID = InstitutionID2, InstitutionID2 = InstitutionID, N)] %>%
    gg_insflow()
     

gd_deg_info(g_insflow) %>% merge(dt_af_instns, by = "InstitutionID") %>%
    .[order(-deg_in)]


## dt_ppl_evcsum <- gd_ppl_evcsum(dt_exhbppl_evc)

## dt_atofyear <- gd_atofyear(dt_atyear_evc)
