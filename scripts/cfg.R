
gc_plts <- function() {

    list(
        ## FIXME when regression models work properly
        ## p_pminfl_cbn = list(
        ##     rb1 = quote(rb1), rc1 = quote(rc1), rd1 = quote(rd1),
        ##     width = 14, height = 7,
        ##     caption = "influence of PM exhibitions"),
        p_hist_crap = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "nonsense hist"),
        p_activity_jitter = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height =7,
            caption = "activity jitter"),
        p_insgrow5 = list(
            dt_insgrow5 = quote(dt_insgrow5),
            width = 14, height =9,
            caption = "growth relations"),
        p_insgrow5_contours = list(
            dt_insgrow5 = quote(dt_insgrow5),
            width = 14, height =9,
            caption = "growth relations"),
        p_cumad_jitter = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "cumulative advantage whatever"),
        p_insgrow5_dist = list(
            dt_insgrow5 = quote(dt_insgrow5),
            l_cols_summary = l_cols_summary,
            width = 14, height = 7, caption = "grow5 dists"),
        p_mainseq = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "cumulative advantage")
        )}


gc_tbls <- function() {

    list(
        t_grow5_sumry  = list(
            dt_exhbppl_watyear = quote(dt_exhbppl_watyear),
            dt_insgrow5 = quote(dt_insgrow5),
            l_cols_summary = quote(l_cols_summary),
            caption = "summary stats"),
        t_insgrow5_topperc = list(
            dt_insgrow5 = quote(dt_insgrow5),
            l_cols_summary = quote(l_cols_summary),
            caption = "most influential orgs (percentage)"),
        t_deg_info = list(
            g_insflow = quote(g_insflow),
            dt_af_instns_rcd = quote(dt_af_instns_rcd),
            caption = "degree info by form"),
        t_insgrow5_topcnt = list(
            dt_insgrow5 = quote(dt_insgrow5),
            l_cols_summary = quote(l_cols_summary),
            caption = "most influential orgs (cnt)")
        )
}
