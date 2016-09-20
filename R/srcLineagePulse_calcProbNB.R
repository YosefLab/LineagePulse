#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++    Calculate posterior of drop-out  ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Calculate posterior of drop-out
#' 
#' Calculates posterior of observation being a drop-out. This posterior
#' is zero if an observation is non-zero, therefore, the data is read in 
#' as a boolean matrix indicating zero-observations, the actual values
#' are not required. Neighbourhood smoothing can be included.
#' 
#' @seealso Called by \code{fitZINB} and \code{runLineagePulse}.
#' 
#' @param matMu: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param matDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param matboolZero: (bool matrix genes x cells)
#'    Whether observation is zero.
#' @param matboolNotZeroObserved: (bool matrix genes x cells)
#'    Whether observation is real and non-zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return matProbNB: (numeric matrix genes x cells)
#'    Posterior probability of observation not being generated 
#'    by drop-out.
#' @export

calcProbNB <- function( matMu,
  matDispersions,
  matDropout,
  matboolZero,
  matboolNotZeroObserved,
  scaWindowRadius ){
  
  scaNumGenes <- dim(matMu)[1]
  scaNumCells <- dim(matMu)[2]
  
  # Compute probability of zero counts under 
  # negative binomial model.
  matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
  }))
  # Compute posterior of drop-out.
  matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    vecZ <- sapply(seq(1,scaNumCells), function(j){
      if(matboolNotZeroObserved[i,j]){
        scaZ <- 0
      } else if(matboolZero[i,j]) {
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        scaZ <- sum(matDropout[i,j]/(matDropout[i,j] + 
            (1-matDropout[i,j])*matNBZero[i,vecInterval])) *
          1/length(vecInterval)
      } else {
        scaZ <- NA
      }
      return(scaZ)
    })
    return(vecZ)
  }))
  
  return(matZ)
}