#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++    Calculate posterior of drop-out  ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Calculate posterior of drop-out
#' 
#' Calculates posterior of observation being a drop-out for a vector. This posterior
#' is zero if an observation is non-zero, therefore, the data is read in 
#' as a boolean matrix indicating zero-observations, the actual values
#' are not required. Neighbourhood smoothing can be included. Does not use
#' parallelisation and takes vectors: As oppose to \code{calcProbNB_Vector}
#' vectors do not produce format errors and usage within upstream parallelisation
#' works.
#' 
#' @seealso Called by \code{fitZINB}. Same function as \code{calcPostDrop_Matrix}
#' but for vectors (i.e. one gene) as oppose to matrices: no 
#' parallelisation and no formatting errors.
#' 
#' @param vecMu: (numeric vector samples)
#'    Negative binomial mean parameters of samples.
#' @param vecDisp: (numeric vector samples)
#'    Negative binomial mean parameters of samples.
#' @param vecDrop: (numeric vector samples)
#'   Drop out rates of samples.
#' @param vecboolZero: (bool vector samples)
#'    Whether observation is zero.
#' @param vecboolNotZero: (bool vector samples)
#'    Whether observation is real and non-zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return vecZ:  (numeric vector samples)
#'    Posterior probability of observation not being generated 
#'    by drop-out.
#'    
#' @author David Sebastian Fischer
#' 
#' @export

calcPostDrop_Vector <- function( vecMu,
  vecDisp,
  vecDrop,
  vecboolZero,
  vecboolNotZero,
  scaWindowRadius ){
  
  scaNumSamples <- length(vecMu)
  
  # Compute probability of zero counts under 
  # negative binomial model.
  vecNBZero <- (vecDisp/(vecDisp+vecMu))^vecDisp
  # Compute posterior of drop-out.
  vecZ <- sapply(seq(1,scaNumSamples), function(j){
    if(vecboolNotZero[j]){
      scaZ <- 0
    } else if(vecboolZero[j]) {
      if(!is.null(scaWindowRadius)){
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumSamples,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
      } else {
        vecInterval <- j
      }
      scaZ <- sum(vecDrop[j]/(vecDrop[j] + 
          (1-vecDrop[j])*vecNBZero[vecInterval])) *
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
#' Calculates posterior of observation being a drop-out for a matrix. This posterior
#' is zero if an observation is non-zero, therefore, the data is read in 
#' as a boolean matrix indicating zero-observations, the actual values
#' are not required. Neighbourhood smoothing can be included.
#' 
#' @seealso Called by \code{fitZINB}. Can be called by user.
#' 
#' @param matMu: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param matDisp: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matDrop: (numeric matrix genes x cells)
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

calcPostDrop_Matrix <- function( matMu,
  matDisp,
  matDrop,
  matboolZero,
  matboolNotZeroObserved,
  scaWindowRadius ){
  
  scaNumGenes <- dim(matMu)[1]
  scaNumCells <- dim(matMu)[2]
  
  # Compute probability of zero counts under 
  # negative binomial model.
  matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    (matDisp[i,]/(matDisp[i,]+matMu[i,]))^matDisp[i,]
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
        scaZ <- sum(matDrop[i,j]/(matDrop[i,j] + 
            (1-matDrop[i,j])*matNBZero[i,vecInterval])) *
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