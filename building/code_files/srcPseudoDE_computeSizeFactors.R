#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++#     Compute Size factors    ++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute size factors for a scRNA-seq dataset
#' 
#' This function computes size factors for each sample
#' in the dataset and expands them to a matrix of the size
#' of the dataset.
#' Size factors scale the negative binomial likelihood
#' model of a gene to the sequencing depth of each sample.
#' Size factors are normalised relative sequencing depth of 
#' each cell. Note that size factors are usually differently computed
#' for bulk data (c.f. ImpulseDE2) but this procedure might
#' be unstable for single-cell data.
#' 
#' @seealso Called by \code{computeNormConst}.
#' 
#' @param matCountsProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCounts}. 
#'    For internal use.
#' 
#' @return vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @export

computeSizeFactors_LineagePulse <- function(matCountsProc){
  
  # Size factors directly represent sequencing depth:
  # Normalised relative sequencing depth.
  vecSeqDepth <- apply(matCountsProc, 2,
    function(cell){ sum(cell, na.rm=TRUE) })
  vecSizeFactors <- vecSeqDepth/sum(vecSeqDepth)*length(vecSeqDepth)
  
  if(any(vecSizeFactors==0)){
    warning("WARNING IN LINEAGEPULSE: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  
  return(vecSizeFactors)
}