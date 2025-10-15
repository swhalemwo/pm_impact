## distribution of exhibitions by form: mostly gallery and pes
melt(dt_atyear, id.vars = c("PeopleID", "year"), measure.vars = patterns("exhbN_.*_s5"),
     variable.name = "vrbl", value.name = "vlu") %>%
    .[grepl("s5", vrbl)] %>% na.omit %>%
    ggplot(aes(x=log1p(vlu), color = vrbl)) +
    ## geom_density(bw = 0.2) +
    geom_histogram() + 
    facet_grid(vrbl~., switch = "y", scales = "free") +
    theme(strip.text.y.left = element_text(angle = 0))


## descriptive plots
## exhbN_all_cum ~ cy: looks cute, but color is essentially also only cumulative advantage
dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
    ggplot(aes(x=career_year, y = exhbN_all_cum, color = exhbN_all)) +
    geom_jitter(size = 0.4, height = 0.01) +
    scale_y_continuous(trans = scales::log10_trans()) +
    scale_color_gradient(low = "red", high = "yellow",
                         trans = scales::log10_trans()) +
    theme(panel.background = element_rect(fill = "black"),
          panel.grid.major = element_line(color = "black"),
          panel.grid.minor = element_line(color = "black")
          )


## exhbN_all ~ cy: idk what's the point
dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
    ggplot(aes(x=career_year, y = exhbN_all + 1)) + geom_jitter(size = 0.1, height = 0.03) +
    scale_y_continuous(trans = scales::log10_trans())

## exhbN_all ~ exhnb_all_cum: muh cumulative advantage
dt_atyear_foster %>% .[sample(1:.N, size = 1e4)] %>%
    ggplot(aes(x=exhbN_all_cum, y = exhbN_all + 1)) +
    geom_jitter(size = 0.1, # width = 0.003, if x not log
                width = 0.03, # if x log
                height = 0.03) +
    scale_y_continuous(trans = scales::log10_trans()) +
    scale_x_continuous(trans = scales::log10_trans()) 
