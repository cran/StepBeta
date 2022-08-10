###
#' @title Stepwise model selection for Beta-Binomial and Negative Binomial Regressions from aod package
#' @description
#' This function performs a stepwise algorithm to define the best linear predictor according to an user defined criterion (defeault is the Akaike Information Criterion aka AIC, but it is also possible to perform the corrected version AICc).
#' It works only for object from betabin function (class "glimML" from "aod" package). If the object is different from "glimMl" class, the function performs the classical step function in "stats" package.
#' @author
#' Sergio Garofalo
#'
#' @details
#' Step_glimML is different from step (stats) and stepAIC (MASS) functions; for an object of class betabin is impossible to use an algorithm which uses the function extractAIC.
#' Starting from a full model it provides a backaward procedure where the scope model is the reduced one.
#'
#' First, Step_glimML operates with all the principal effects included in the model; starting from the full model, the algorithm computes all the possible models, it calculates the measure (default is AIC) and it defines as a good predictor the model with lower AIC.
#'
#' Then, based on the previous results, Step_glimML operates adding all the possible interactive effects. As in the first passage, the model choosen by the algorithm is the one whose AIC is the lowest.
#'
#' During the procedure, Step_glimML considers all the possible models which betabin can fit. There are many cases where betabin function falls into error, in these cases the algorithm does not consider the linear predictor which causes the error and it goes forward.
#'
#' @return
#' The algorithm returns an object of class "glimML"
#'
#' @references
#' Crowder, M.J., 1978. Beta-binomial anova for proportions. Appl. Statist. 27, 34-37.
#'
#' Lawless, J.F., 1987. Negative binomial and mixed Poisson regression. The Canadian Journal of Statistics, 15(3): 209-225.
#'
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
#'
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. (has iris3 as iris.)
#'
#' @import stats
#' @importFrom aod negbin betabin
#' @import glue
#' @import MASS
#' @export
#' @param object Object of class "glimML". If the class is different the function apply step function in "stats" package
#' @param k The penalty parameter used for the criterion, e.g. default is k = 2 which identify the classical AIC. BIC can be obtained as k = log(n)
#' @param overdispersion Provide the stepwise procedure also for the overdispersion component of the model (defined as random) Default is TRUE
#' @param correctAIC Use AICc instead of AIC. Default TRUE is for AICc
#' @examples
#' ## Starting from a "betabinom" model
#'
#' ## Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language.
#' ## Wadsworth & Brooks/Cole. (has iris3 as iris.)
#'
#' ## Prepare the data
#'
#' library(aod)
#' data(iris)
#'
#' ############ Beta Binomial model
#'\dontrun{
#' n <- round(runif(dim(iris)[1],1,50))
#' y <- round(runif(length(n), 1,n))
#' data <- cbind(iris,y,n)
#' fullModel <- betabin(cbind(y, n - y) ~ Sepal.Width * Petal.Length + Petal.Width, ~ Species,
#'                      data = data)
#' reducedModel <- Step_glimML(fullModel)
#' summary(reducedModel)
#'}
#' ############ Negative Binomial model
#' \dontrun{
#' data <- iris
#' data$Sepal.Length <- round(Sepal.length + runif(dim(data)[1],0,1) * 100)
#' fullModel <- negbin(Sepal.Length ~ Sepal.Width * Petal.Length + Petal.Width, ~ Species,
#'                     data = data)
#' reducedModel <-Step_glimML(fullModel)
#' summary(reducedModel)
#' }
#'
Step_glimML <- function(object, k = 2, overdispersion = T, correctAIC = T){
  if(requireNamespace("aod")){
    if(!inherits(object, "glimML")){
      step(object, direction = "both")
    } else{

      full_model <- object

      #Control terms
      ##############
      beta_binomial_control_link <- full_model@link
      ##############
      formula_full_model <- check_formula_terms(full_model)

      #OVERDISPERSION?
      if(as.character(full_model@random)[2] !=1){
        overdispersion_model <- 1
      } else{
        overdispersion_model <- NULL
      }

      #Stepwise both for dispersion parameter
      if(overdispersion == T & !is.null(overdispersion_model)){
        overdispersion_term <- c(paste("~",as.character(full_model@random)[2]),paste("~",1))
      }else{
        overdispersion_term <- paste("~",1)
      }
      if(isS4(AIC(full_model,k = k)) == T){
        if(correctAIC == F){
          full_AIC <- AIC(full_model, k = k)@istats[[2]]
        } else if(correctAIC == T){
          full_AIC <- AIC(full_model, k = k)@istats[[3]]
        }else{
          stop("For models of class glimML it is only possible to abtain AIC or AICc. Value of correctAIC must be TRUE or FALSE")
        }
      } else{
        full_AIC <- NA
      }

      Final_Results <- as.data.frame(matrix(ncol = 3))
      colnames(Final_Results) <- c("Random","Fixed","AIC")

      for(g in 1:length(overdispersion_term)){

        starting_NoInt_AIC <- 0
        Terms <- strsplit(as.character(full_model@formula[3]),'*',fixed = T)[[1]]
        Terms <- unlist(strsplit(Terms,'+',fixed = T))
        Terms <- unlist(strsplit(Terms,':',fixed = T))

        Terms <- gsub(" ","",Terms)
        #is.Terms_NoInt <- Terms[-c(grep(":",Terms))]
        #if(length(is.Terms_NoInt) > 0){
        #  Terms_NoInt <- is.Terms_NoInt
        #} else {
        #Terms_NoInt <- Terms
        #}
        Terms_NoInt <- Combination_Terms(Terms)
        Results <- as.data.frame(matrix(ncol=3))
        names(Results) <- c("Random","Fixed","AIC")
        models <- list()
        starting_variables <- NULL

          i <- 1
          while(i < length(Terms_NoInt) + 1){
            new_formula <- formula(paste(full_model@formula[2], full_model@formula[1], Terms_NoInt[i]))
            if(full_model@method == "BB"){
              mod_updated <-try(betabin(new_formula,formula(overdispersion_term[g]),link = beta_binomial_control_link,
                               data = object@data),silent = T)
              if(isTRUE(class(mod_updated) == "try-error")) {
                models[[i]] <- NA
                i <- i + 1
                next
              }
            } else if(full_model@method == "NB"){
              mod_updated <-try(negbin(new_formula,formula(overdispersion_term[g]),
                                         data = object@data),silent = T)
              if(isTRUE(class(mod_updated) == "try-error")) {
                models[[i]] <- NA
                i <- i + 1
                next
              }
            }
            models[[i]] <- mod_updated
            Results[i,1] <- as.character(mod_updated@random)[2]
            Results[i,2] <- as.character(mod_updated@formula)[3]
            if(isS4(AIC(full_model,k = k)) == T){
              if(correctAIC == F){
                Results[i,3] <- AIC(mod_updated, k = k)@istats[[2]]
              } else if(correctAIC == T){
                Results[i,3] <- AIC(mod_updated, k = k)@istats[[3]]
              }
            } else {
              Results[i,3] <- NA
            }
            i <- i + 1
          }
          Var_reduced <- which.min(Results[,3])
          AIC_reduced <- Results[Var_reduced,3]
          starting_NoInt_AIC <- AIC_reduced

          formula_reduced <- formula(paste(full_model@formula[2],full_model@formula[1],Results[Var_reduced,2]))
          if(full_model@method == "BB"){
            mod_reduced <-try(betabin(new_formula,formula(overdispersion_term[g]),link = beta_binomial_control_link,
                                       data = object@data),silent = T)
            if(isTRUE(class(mod_updated) == "try-error")) {
              models[[i]] <- NA
              i <- i + 1
              next
            }
          } else if(full_model@method == "NB"){
            mod_reduced <-try(negbin(new_formula,formula(overdispersion_term[g]),
                                      data = object@data),silent = T)
            if(isTRUE(class(mod_updated) == "try-error")) {
              models[[i]] <- NA
              i <- i + 1
              next
            }
          }
        new_formula <- mod_reduced@formula
        starting_IntAic <- AIC_reduced
        Terms_new_formula <- unlist(strsplit(as.character(new_formula[3]), "+",fixed = T))
        #####
        #Interactions

        Terms_Int <- Combination_Terms(Terms_new_formula,interaction=T)
        if(length(Terms_Int)> 0){
        if(identical(Terms_Int, character(0))){
          final_model <- mod_reduced
        } else if(!identical(Terms_Int, character(0)))  {
          Results <- as.data.frame(matrix(ncol=3))
          names(Results) <- c("Random","Fixed","AIC")
          models <- list()
          starting_variables <- NULL
          diff <- 1

            i <- 1
            while(i < length(Terms_Int) + 1){

              new_formula <- formula(paste(full_model@formula[2],full_model@formula[1], mod_reduced@formula[3], "+", Terms_Int[i]))
              if(full_model@method == "BB"){
                mod_updated <-try(betabin(new_formula,formula(overdispersion_term[g]),link = beta_binomial_control_link,
                                           data = object@data),silent = T)
                if(isTRUE(class(mod_updated) == "try-error")) {
                  models[[i]] <- NA
                  i <- i + 1
                  next
                }
              } else if(full_model@method == "NB"){
                mod_updated <-try(negbin(new_formula,formula(overdispersion_term[g]),
                                          data = object@data),silent = T)
                if(isTRUE(class(mod_updated) == "try-error")) {
                  models[[i]] <- NA
                  i <- i + 1
                  next
                }
              }
              models[[i]] <- mod_updated
              Results[i,1] <- as.character(mod_updated@random)[2]
              Results[i,2] <- as.character(new_formula)[3]
              if(isS4(AIC(full_model,k = k)) == T){
                if(correctAIC == F){
                  Results[i,3] <- AIC(mod_updated, k = k)@istats[[2]]
                } else if(correctAIC == T){
                  Results[i,3] <- AIC(mod_updated, k = k)@istats[[3]]
                }
              } else {
                Results[i,3] <- NA
              }
              i <- i + 1
            }
            Results[i, 1] <- as.character(mod_reduced@random)[2]
            Results[i, 2] <- as.character(mod_reduced@formula)[3]
            if(isS4(AIC(full_model,k = k)) == T){
              if(correctAIC == F){
                Results[i,3] <- AIC(mod_reduced, k = k)@istats[[2]]
              } else if(correctAIC == T){
                Results[i,3] <- AIC(mod_reduced, k = k)@istats[[3]]
              }
            } else{
              Results[i,3] <- NA
            }
            models[[i]] <- mod_reduced
            Var_reduced <- which.min(Results[,3])

            AIC_reduced <- Results[Var_reduced,3]

            diff <- abs(starting_NoInt_AIC - AIC_reduced)
            starting_NoInt_AIC <- AIC_reduced

        }
        }
        Final_Results <- rbind(Final_Results,Results)
        }

    }
  }
  Final_Results <- Final_Results[-1,]
  final_formula_fixed <- formula(paste(full_model@formula[2],full_model@formula[1],Final_Results$Fixed[which.min(Final_Results$AIC)]))
  final_formula_random <- formula(paste("~",as.numeric(Final_Results$Random[which.min(Final_Results$AIC)])))

  if(full_model@method == "BB"){
    final_model <-aod::betabin(final_formula_fixed,final_formula_random,link = beta_binomial_control_link,
                               data = object@data)
  } else if(full_model@method == "NB"){
    final_model <-aod::negbin(final_formula_fixed,final_formula_random,
                              data = object@data)
  }
 return(final_model)
}


