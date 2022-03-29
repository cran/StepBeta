#' @title StepBeta internal object
#' @import stats
#' @export
#' @param model Beta regression model
#' @return It returns the complete formula in a standard form

check_formula_terms <- function(model){
  if(formula(model)[3] == ".()"){
    f <- paste0(attr(model$model,"terms"))
    formula_full_model <- as.formula(paste0(f[2],f[1],f[3]))
  } else {
    formula_full_model <- formula(model)
  }
  return(formula_full_model)
}
