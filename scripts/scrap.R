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
