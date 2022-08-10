#' @title StepBeta internal object
#' @import stats
#' @import combinat
#' @export
#' @param object full model
#' @return The function updates the formula for the dispersion component of the model

dispersion_formula_terms <- function(object){
  Terms <- attr(object$terms$precision,"term.labels")
  main_terms <- Terms[-grep("\\:",Terms)]
    if(requireNamespace("combinat")){
      main_terms_combinations <- unlist(lapply(1:length(Terms),    # Get all combinations
                                                            combinat::combn,
                                                            x = Terms,
                                                            simplify = F),
                                                           recursive = F)
    }
  Output_formula <- list()
  i <- 1
  k <- 1
  while (i < length(main_terms_combinations) + 1 ) {
    if(length(grep("\\:",main_terms_combinations[[i]])) > 0){
      Interactions <- unique(unlist(strsplit(main_terms_combinations[[i]][grep("\\:",main_terms_combinations[[i]])],"\\:")))
    } else{
      Interactions <- NULL
    }
    if(!is.null(Interactions) & all(Interactions %in% main_terms_combinations[[i]][-grep("\\:",main_terms_combinations[[i]])]) == F){
      i <- i + 1
    } else{
      Output_formula[[k]] <- paste("|",gsub(" ","+",paste(main_terms_combinations[[i]] , collapse = " ")))
      i <- i + 1
      k <- k + 1
    }
  }
  return(Output_formula)
}
