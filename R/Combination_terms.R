#' @title StepBetaBinomial internal object
#' @import glue
#' @import stats
#' @export
#' @param Terms Variables from the starting model
#' @param interaction Parameter to define which part of linear predictor to operate
#' @return The function create alle possible combination of the linear predictor


Combination_Terms <- function(Terms, interaction = F){
  Combination <- vector()
  if(interaction == F){
    if(length(Terms) > 1){
    for(i in 1:length(Terms)){
      if(i == 1){
        Combination <- c(Combination,t(combn(Terms,i)))
      } else {
        if(requireNamespace("glue")){
          Partial_comb <- apply(t(combn(Terms,i)),1,function(x){glue::glue_collapse(x,sep = "+")})
        }
          Combination <- c(Combination,Partial_comb)
        }
      }
    }else{
      Combination <- Terms
    }
    Combination <- c("1",Combination)
  } else if(interaction == T){
    if(length(Terms) > 1){
   Int_terms <- t(combn(Terms,2))
   if(requireNamespace("glue")){
     Partial_int <- apply(Int_terms,1,function(x){glue::glue_collapse(x,sep = ":")})
   }
   if(length(Partial_int) > 1){
     for(i in 1:length(Partial_int)){
       if(i == 1){
         Combination <- c(Combination,t(combn(Partial_int,i)))
       } else {
         if(requireNamespace("glue")){
           Partial_comb <- apply(t(combn(Partial_int,i)),1,function(x){glue::glue_collapse(x,sep = "+")})
         }
         Combination <- c(Combination,Partial_comb)
       }
     }
   }else{
     Combination <- Partial_int
   }
    }else{
      Combination <- NULL
    }
  }
 return(Combination)
}
