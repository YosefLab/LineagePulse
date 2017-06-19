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
#' @param matCounts: (numeric matrix genes x cells)
#'    Count data.
#' @param lsMuModel: (list) Mean parameter model parameters.
#' @param boolDepth: (bool) [Default TRUE] Whether to normalize for sequencing depth.
#' @param boolBatch: (bool) [Default TRUE] Whether to normalize for batch.
#' 
#' @return matZ: (numeric matrix genes x cells)
#'    Posterior probability of observation not being generated 
#'    by drop-out.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
getNormData <- function(matCounts,
                        lsMuModel,
                        boolDepth = TRUE,
                        boolBatch = TRUE){
    
    if(!boolDepth & !boolBatch) {
        warning("No normalisation chosen")
    }
    
    scaNumCells <- dim(matCounts)[2]
    scaNumConfounders <- length(lsMuModel$lsmatBatchModel)
    matNormData <- do.call(rbind, bplapply(rownames(matCounts), function(id){
        # Extract batch parameters
        vecBatchParam <- array(1, scaNumCells)
        if(!is.null(lsMuModel$lsMuModelGlobal$vecConfounders)){
            for(confounder in seq(1,scaNumConfounders)){
                vecBatchParam <- vecBatchParam*(lsMuModel$lsmatBatchModel[[confounder]][id,][
                    lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign[[confounder]]])
            }
        }
        # Normalise counts by depth and/or batch factors as required
        vecNormData <- matCounts[id,]
        if(boolDepth) {
            vecNormData <-vecNormData/lsMuModel$lsMuModelGlobal$vecNormConst
        }
        if(boolBatch) {
            vecNormData <- vecNormData/vecBatchParam
        }
        return(vecNormData)
    }))
    
    rownames(matNormData) <- rownames(matCounts)
    colnames(matNormData) <- colnames(matCounts)
    return(matNormData)
}