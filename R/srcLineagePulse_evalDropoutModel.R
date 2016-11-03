#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Evaluate logistic drop-out model    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of impulse function given parameters.
#' 
#' Compute value of impulse function given parameters. Contains
#' log linker for amplitude paramteres. Low impulse model values
#' are caught to avoid mean parameters close to or below zero.
#' 
#' @aliases evalImpulseModel_comp
#' 
#' @param vecPiModel: (numeric vector length linear model) Linear model
#'    for drop-out rate in logit space.
#' @param vecPiPredictors: (vector length of predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' 
#' @return vecY (vec number of vecTimepoints) 
#'    Model expression values of given gene for time points
#' @export

evalDropoutModel <- function(vecPiModel, vecPiPredictors){
  
  # Set offset parameter
  scaOffset <- 0.01
  
  #vecDropoutRateFit <- 1/(1+exp(-vecLinModelOut))
  #vecDropoutRateFit[vecDropoutRateFit < scaOffset] <- scaOffset
  #vecDropoutRateFit[vecDropoutRateFit > 1-scaOffset] <- 1-scaOffset
  scaDropoutRate <- scaOffset+(1-scaOffset)*1/(1+exp(- (vecPiPredictors %*% vecPiModel )))
  
  return(scaDropoutRate)
}