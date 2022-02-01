#' @title StepBeta internal object
#' @import glue
#' @import stats
#' @export
#' @param the_formula Formula of Beta Regression model
#' @param var_name Names of the variables to keep
#' @return The function updates the formula, it keeps the variables defined by the user

keep_formula_terms <- function(the_formula, var_name){
  if(requireNamespace("glue")){
    var_name <- as.character(glue::glue_collapse(var_name, sep = "+"))

    fmla <- as.formula(paste(paste(the_formula[2],the_formula[1], var_name)))
    return(fmla)
  }
}
