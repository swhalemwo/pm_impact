
gc_plts <- function() {

    list(
        p_pminfl_cbn = list(
            rb1 = quote(rb1), rc1 = quote(rc1), rd1 = quote(rd1),
            width = 14, height = 7,
            caption = "influence of PM exhibitions"),
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
        p_mainseq = list(
            dt_atyear_foster = quote(dt_atyear_foster),
            width = 14, height = 7,
            caption = "cumulative advantage")
        )}


gc_tbls <- function() {

    list()
}
