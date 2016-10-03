#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++    Calculate posterior of drop-out  ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Calculate posterior of drop-out
#' 
#' Calculates posterior of observation being a drop-out. This posterior
#' is zero if an observation is non-zero, therefore, the data is read in 
#' as a boolean matrix indicating zero-observations, the actual values
#' are not required. Neighbourhood smoothing can be included. Does not use
#' parallelisation and takes vectors: As oppose to \code{calcProbNB_Vector}
#' vectors do not produce format errors and usage within upstream parallelisation
#' works.
#' 
#' @seealso Called by \code{fitZINB}. Same function as \code{calcProbNB_Vector}
#' but for vectors (i.e. one gene) as oppose to matrices: no 
#' parallelisation and no formatting errors.
#' 
#' @param vecMu: (numeric vector samples)
#'    Negative binomial mean parameters of samples.
#' @param vecDispersions: (numeric vector samples)
#'    Negative binomial mean parameters of samples.
#' @param vecDropout: (numeric vector samples)
#'   Drop out rates of samples.
#' @param vecboolZero: (bool matrix genes x cells)
#'    Whether observation is zero.
#' @param vecboolNotZeroObserved: (bool matrix genes x cells)
#'    Whether observation is real and non-zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return vecZ:  (numeric vector samples)
#'    Posterior probability of observation not being generated 
#'    by drop-out.
#' @export

calcPostDrop_Vector <- function( vecMu,
  vecDispersions,
  vecDropout,
  vecboolZero,
  vecboolNotZeroObserved,
  scaWindowRadius ){
  
  scaNumSamples <- length(vecMu)
  
  # Compute probability of zero counts under 
  # negative binomial model.
  vecNBZero <- (vecDispersions/(vecDispersions+vecMu))^vecDispersions
  # Compute posterior of drop-out.
  vecZ <- sapply(seq(1,scaNumSamples), function(j){
    if(vecboolNotZeroObserved[j]){
      scaZ <- 0
    } else if(vecboolZero[j]) {
      if(!is.null(scaWindowRadius)){
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
      } else {
        vecInterval <- j
      }
      scaZ <- sum(vecDropout[j]/(vecDropout[j] + 
          (1-vecDropout[j])*vecNBZero[vecInterval])) *
        1/length(vecInterval)
    } else {
      scaZ <- NA
    }
    return(scaZ)
  })
  
  names(vecZ) <- names(vecMu)
  return(vecZ)
}

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
#' @return matZ: (numeric matrix genes x cells)
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
        if(!is.null(scaWindowRadius)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        } else {
          vecInterval <- j
        }
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
  
  rownames(matZ) <- rownames(matMu)
  colnames(matZ) <- colnames(matMu)
  return(matZ)
}