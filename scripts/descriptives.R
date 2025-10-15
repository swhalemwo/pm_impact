## distribution of exhibitions by form: mostly gallery and pes
melt(dt_atyear, id.vars = c("PeopleID", "year"), measure.vars = patterns("exhbN_.*_s5"),
     variable.name = "vrbl", value.name = "vlu") %>%
    .[grepl("s5", vrbl)] %>% na.omit %>%
    ggplot(aes(x=log1p(vlu), color = vrbl)) +
    ## geom_density(bw = 0.2) +
    geom_histogram() + 
    facet_grid(vrbl~., switch = "y", scales = "free") +
    theme(strip.text.y.left = element_text(angle = 0))
