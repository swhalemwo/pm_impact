
gf_model_formula <- function(dv, ivs, cluster = "PeopleID", fe = c("PeopleID", "year"), offset = NULL) {
    #' Generate formula for fepois/feols models
    #' @param dv Dependent variable (character)
    #' @param ivs Vector of independent variables (character)
    #' @param cluster Clustering variable (character)
    #' @param fe Fixed effects (character vector)
    #' @param offset Optional offset (character)
    
    iv_str <- paste(ivs, collapse = " + ")
    fe_str <- paste(fe, collapse = " + ")
    
    formula_str <- sprintf("%s ~ %s | %s", dv, iv_str, fe_str)
    
    if (!is.null(offset)) {
        formula_str <- paste0(formula_str, " + offset(", offset, ")")
    }
    
    return(as.formula(formula_str))
}
