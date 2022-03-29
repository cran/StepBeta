###
#' @title Stepwise model selection for Beta Regression
#' @description
#' This function performs a stepwise algorithm to define the best linear predictor according to an user defined criterion (defeault is the Akaike Information Criterion aka AIC).
#' It works for objects of class "betareg". If the object is different from "betareg" class, the function performs the classical "step" function in "stats" package.
#' @author
#' Sergio Garofalo
#'
#' @details
#' StepBeta is different from step (stats) and stepAIC (MASS) functions; for an object of class "betareg" is impossible to use an algorithm which uses the function extractAIC.
#' Starting from a full model it provides a backaward procedure where the scope model is the reduced one.
#'
#' First, StepBeta operates with all the principal effects included in the model; starting from the full model, the algorithm computes all the possible models, it calculates the measure (default is AIC) and it defines as a good predictor the model with lower AIC.
#'
#' Then, based on the previous results, StepBeta operates adding all the possible interactive effects. As in the first passage, the model choosen by the algorithm is the one whose AIC is the lowest.
#'
#' During the procedure, StepBeta considers all the possible models which betareg can fit. There are many cases where betareg function falls into error, in these cases the algorithm does not consider the linear predictor which causes the error and it goes forward.
#'
#' @return
#' The algorithm returns an object of class "betareg"
#'
#' @references
#' Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R. Journal of Statistical Software, 34(2), 1--24. 10.18637/jss.v034.i02
#'
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
#'
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. (has iris3 as iris.)
#'
#' @import stats
#' @importFrom betareg betareg
#' @import glue
#' @export
#' @param object Object of class "betareg". If the class is different the function apply the step function in "stats" package
#' @param k The penalty parameter used for the criterion, e.g. default is k = 2 which identify the classical AIC. BIC can be obtained as k = log(n)
#' @examples
#' ## Starting from a "betareg" model
#'
#' ## Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language.
#' ## Wadsworth & Brooks/Cole. (has iris3 as iris.)
#'
#'
#' ## Prepare the data
#'
#' library(betareg)
#' data <- iris
#' data$Sepal.Length <- data$Sepal.Length/(max(data$Sepal.Length) + 0.01)
#'
#'
#' fullModel <- betareg(Sepal.Length ~ Sepal.Width * Petal.Length *
#'                                     Petal.Width * Species, data = data)
#' reducedModel <- StepBeta(fullModel)
#'
#' summary(reducedModel)
#'
StepBeta <- function(object, k = 2){
  if(requireNamespace("betareg")){
  if(class(object) != "betareg"){

    step(object, direction = "both")
    } else{

    full_model <- object

    formula_full_model <- check_formula_terms(full_model)

    if(!identical(grep("\\|", as.character(formula_full_model[3])), integer(0)) ){
      full_model_VertBar <- paste("|",gsub(".*\\|","",as.character(formula_full_model[3])))
    } else {
      full_model_VertBar <- ""
    }
    full_AIC <- AIC(full_model, k = k)

    formula_NoInt <- remove_formula_interactions(formula_full_model)
    starting_NoInt_AIC <- 0

    Terms <- attr(object$terms$full,"term.labels")
    is.Terms_NoInt <- Terms[-c(grep(":",Terms))]
      if(length(is.Terms_NoInt) > 0){
        Terms_NoInt <- is.Terms_NoInt
      } else {
        Terms_NoInt <- Terms
      }

    Results <- as.data.frame(matrix(ncol=2))
    names(Results) <- c("Model","AIC")
    models <- list()
    starting_variables <- NULL
    diff <- 1

    while(diff > 5e-15){
      i <- 1
      while(i < length(Terms_NoInt) + 1){
        new_formula <- keep_formula_terms(formula_NoInt,c(unique(Terms_NoInt[c(starting_variables,i)])))
        new_formula <- as.formula(paste(new_formula[2], new_formula[1], new_formula[3],full_model_VertBar))
        mod_updated <- try(betareg(new_formula, data = object$model),T)
        if(isTRUE(class(mod_updated) == "try-error")) {
          models[[i]] <- NA
          i <- i + 1
          next
        }

        models[[i]] <- mod_updated
        Results[i,1] <- as.character(new_formula)[3]
        Results[i,2] <- AIC(mod_updated, k = k)
        i <- i + 1
      }
      Var_reduced <- which.min(Results[,2])
      starting_variables <- c(starting_variables ,Var_reduced)
      AIC_reduced <- Results[Var_reduced,2]

      diff <- abs(starting_NoInt_AIC - AIC_reduced)
      starting_NoInt_AIC <- AIC_reduced
    }

    mod_reduced <- models[[which.min(Results[,2])]]

    new_formula <- formula(mod_reduced)
    starting_IntAic <- AIC_reduced

    Terms_Int <- attr(full_model$terms$full,"term.labels")
    Terms_Int <- Terms_Int[-c(1:length(Terms_NoInt))]

    if(identical(Terms_Int, character(0))){
      final_model <- mod_reduced
    } else if(!identical(Terms_Int, character(0)))  {
      Results <- as.data.frame(matrix(ncol=2))
      names(Results) <- c("Model","AIC")
      models <- list()
      starting_variables <- NULL
      diff <- 1

      while(diff > 5e-15){
        i <- 1
        while(i < length(Terms_Int) + 1){
          new_formula <- keep_formula_terms(formula_full_model,c(unique(Terms_Int[c(starting_variables,i)])))
          new_formula <- as.formula(paste(new_formula[2],new_formula[1],as.character(formula(mod_reduced)[3]), "+", new_formula[3],full_model_VertBar))
          mod_updated <- try(betareg(new_formula, data = object$model),T)
          if(isTRUE(class(mod_updated) == "try-error")) {
            i <- i + 1
            models[[i]] <- NA
            next
            }
          models[[i]] <- mod_updated
          Results[i,1] <- as.character(new_formula)[3]
          Results[i,2] <- AIC(mod_updated, k = k)
          i <- i + 1
        }
        Results[i, 1] <- as.character(formula(mod_reduced)[3])
        Results[i, 2] <- AIC_reduced
        models[[i]] <- mod_reduced
        Var_reduced <- which.min(Results[,2])
        starting_variables <- c(starting_variables ,Var_reduced)
        AIC_reduced <- Results[Var_reduced,2]

        diff <- abs(starting_NoInt_AIC - AIC_reduced)
        starting_NoInt_AIC <- AIC_reduced
      }
    }

    final_model <- models[[which.min(Results[,2])]]
    return(final_model)
    }
  }
}

