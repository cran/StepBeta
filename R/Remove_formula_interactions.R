#' @title StepBeta internal object
#' @import stats
#' @export
#' @param the_formula Formula of Beta Regression model
#' @return The function returns a reduced form of the formula. It excludes the interactive effects.

remove_formula_interactions <- function(the_formula){
  if(!identical(grep("\\:", as.character(the_formula[3])), integer(0)) | !identical(grep("\\*", as.character(the_formula[3])), integer(0))){
    the_formula <- as.formula(paste(the_formula[2],the_formula[1],gsub("\\|.*","",as.character(the_formula[3]))))
    var_position <- grep(
      ":", attr(terms(the_formula), "term.labels"))
    update(
      the_formula,
      drop.terms(terms(the_formula), var_position,
                 keep.response=TRUE))
  } else{
    return(the_formula)
  }
}
