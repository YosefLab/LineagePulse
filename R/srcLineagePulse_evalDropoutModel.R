#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Evaluate logistic drop-out model    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of logistic dropout function given with scale boundaries
#' 
#' Computes value of logistic function for one observation and uses offset
#' to correct model.
#' 
#' @aliases evalDropoutModel_comp
#' 
#' @param vecPiModel: (numeric vector length linear model) Linear model
#'    for drop-out rate in logit space.
#' @param vecPiPredictors: (vector length of predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' 
#' @return scaDropoutRate (scalar): 
#'    Drop-out rate estimate.
#' @export

evalDropoutModel <- function(vecPiModel, vecPiPredictors){
  
  # Set offset parameter
  scaOffset <- 0.001
  
  #vecDropoutRateFit <- 1/(1+exp(-vecLinModelOut))
  #vecDropoutRateFit[vecDropoutRateFit < scaOffset] <- scaOffset
  #vecDropoutRateFit[vecDropoutRateFit > 1-scaOffset] <- 1-scaOffset
  scaDropoutRate <- scaOffset+(1-scaOffset)*1/(1+exp(- (vecPiPredictors %*% vecPiModel )))
  
  return(scaDropoutRate)
}