#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++     Evaluate impulse model value    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of impulse function given parameters.
#' 
#' Compute value of impulse function given parameters. Contains
#' log linker for amplitude paramteres. Low impulse model values
#' are caught to avoid mean parameters close to or below zero.
#' 
#' @aliases evalImpulseModel_comp
#' 
#' @seealso Called by \code{evalLogLikImpulse},\code{evalLogLikMean}, 
#'    \code{plotDEGenes}.
#' 
#' @param vecTheta (vector number of parameters) Numerical vector of impulse 
#'    parameters with the order beta, h0, h1, h2, t1, t2.
#' @param vecTimepoints (vector number vecTimepoints) 
#'    Observed vecTimepoints, numeric.
#' 
#' @return vecY (vec number of vecTimepoints) 
#'    Model expression values of given gene for time points
#' @export

evalImpulseModel <- function(vecTheta,vecTimepoints){
  
  # Extract parameters and apply log linker to amplitudes
  beta = vecTheta[1]
  h0 = exp(vecTheta[2])
  h1 = exp(vecTheta[3])
  h2 = exp(vecTheta[4])
  t1 = vecTheta[5]
  t2 = vecTheta[6]
  
  # Compute value of impulse model
  vecY = unlist(lapply(vecTimepoints, function(x) {(1/h1) * 
      (h0 + (h1-h0) * (1/(1+exp(-beta*(x-t1))))) *
      (h2 + (h1-h2) * (1/(1+exp(beta*(x-t2)))))}))
  
  # Catch ver low model values
  vecY[vecY < .Machine$double.eps] <- .Machine$double.eps
  
  return(vecY)
}