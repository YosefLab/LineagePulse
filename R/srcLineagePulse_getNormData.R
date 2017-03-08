#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++    Get depth and batch corrected data  +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Return depth and batch corrected data
#' 
#' The data normalisation is based on the model normalisation used by and inferred
#' by LineagePulse, e.g. for data visualisation.
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
#'    
#' @author David Sebastian Fischer
#' 
#' @export
getNormData <- function(matCounts,
                        lsMuModel){
  
  scaNumCells <- dim(matCounts)[2]
  scaNumConfounders <- length(lsMuModel$lsmatBatchModel)
  matNormData <- do.call(rbind, bplapply(rownames(matCounts), function(id){
    # Extract batch parameters
    vecBatchParam <- array(1, scaNumCells)
    if(!is.null(lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign)){
      for(confounder in seq(1,scaNumConfounders)){
        vecBatchParam <- vecBatchParam*(lsMuModel$lsmatBatchModel[[confounder]][id,][
          lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign[[confounder]]])
      }
    }
    # Normalise counts by depth and batch factors
    vecNormData <- matCounts[id,]/
      lsMuModel$lsMuModelGlobal$vecNormConst/
      vecBatchParam
    return(vecNormData)
  }))
  
  rownames(matNormData) <- rownames(matCounts)
  colnames(matNormData) <- colnames(matCounts)
  return(matNormData)
}