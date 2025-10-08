

## ** for now, use dt_atyear from infl

fepois(exhbN_all ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID, dt_atyear)

r_auc_ols <- feols(auc_turnover ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID,
                   dt_atyear, fixef.rm = "singleton")

r_auc_poi <- fepois(auc_turnover ~ exhbN_pm_s5 + exhbN_event_s5 + exhbN_gallery_s5 + exhbN_nonpm_s5 | PeopleID,
                    dt_atyear, fixef.rm = "singleton")

screenreg(list(r_auc_ols, r_auc_poi))







# add some code



dt_atyear


