#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++    Calculate posterior of drop-out  +++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Calculate posterior of drop-out
#' 
#' Calculates posterior of observation being a drop-out for a vector.
#' 
#' @seealso For matrices \code{calcPostDrop_Matrix}.
#' 
#' @param vecMu (numeric vector samples)
#' Negative binomial mean parameters of samples.
#' @param vecDisp (numeric vector samples)
#' Negative binomial mean parameters of samples.
#' @param vecDrop (numeric vector samples)
#'   Drop out rates of samples.
#' @param vecboolZero (bool vector samples)
#' Whether observation is zero.
#' @param vecboolNotZero (bool vector samples)
#' Whether observation is real and non-zero.
#' 
#' @return vecZ  (numeric vector samples)
#' Posterior probability of observation not being generated 
#' by drop-out.
#' 
#' @author David Sebastian Fischer
calcPostDrop_Vector <- function( 
    vecMu,
    vecDisp,
    vecDrop,
    vecboolZero,
    vecboolNotZero ){
    
    # Compute probability of zero counts under 
    # negative binomial model.
    vecNB_zero <- (vecDisp/(vecDisp+vecMu))^vecDisp
    # Compute posterior of drop-out
    vecZ_zero <- vecDrop/(vecDrop + (1-vecDrop)*vecNB_zero)
    vecZ <- rep(NA_real_, length(vecMu))
    vecZ[vecboolNotZero] <- 0
    vecZ[vecboolZero] <- vecZ_zero[vecboolZero]

    names(vecZ) <- names(vecMu)
    return(vecZ)
}

#' Calculate posterior of drop-out
#' 
#' Calculates posterior of observation being a drop-out for a matrix.
#' 
#' @seealso Called by \code{plotGene}.
#' 
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' @param vecIDs (vector of strings) [Default NULL]
#' Gene IDs for which posteriors of drop-out are to be computed.
#' 
#' @return matZ (numeric matrix genes x cells)
#' Posterior probability of observation not being generated 
#' by drop-out.
#' 
#' @author David Sebastian Fischer
calcPostDrop_Matrix <- function(
    matCounts,
    lsMuModel,
    lsDispModel, 
    lsDropModel,
    vecIDs=NULL){
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    
    # Compute posterior of drop-out.
    if(!is.null(vecIDs)) {
        vecGenes <- vecIDs
    } else {
        vecGenes <- rownames(matCounts)
    }
    matZ <- do.call(rbind, lapply(vecGenes, function(i){
        # Decompress parameters by gene
        vecMuParam <- decompressMeansByGene(
            vecMuModel=lsMuModel$matMuModel[i,],
            lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, 
                                   function(mat) mat[i,] ),
            lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
            vecInterval=NULL )
        vecDispParam <- decompressDispByGene(
            vecDispModel=lsDispModel$matDispModel[i,],
            lsvecBatchModel=lapply(lsDispModel$lsmatBatchModel, 
                                   function(mat) mat[i,] ),
            lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
            vecInterval=NULL )
        vecPiParam <- decompressDropoutRateByGene(
            matDropModel=lsDropModel$matDropoutLinModel,
            vecMu=vecMuParam,
            vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
            lsDropModelGlobal=lsDropModel$lsDropModelGlobal)
        
        # Evaluate loglikelihood of gene
        vecCounts <- matCounts[i,]
        
        # Compute probability of zero counts under 
        # negative binomial model.
        vecNB_zero <- (vecDispParam/(vecDispParam+vecMuParam))^vecDispParam
        
        # Compute posterior of drop-out
        vecZ_zero <- vecPiParam/(vecPiParam + (1-vecPiParam)*vecNB_zero)
        vecZ <- rep(NA_real_, scaNumCells)
        vecZ[which(!is.na(vecCounts) & vecCounts>0)] <- 0
        vecidxZero <- which(!is.na(vecCounts) & vecCounts==0)
        vecZ[vecidxZero] <- vecZ_zero[vecidxZero]

        return(vecZ)
    }))
    
    rownames(matZ) <- vecGenes
    colnames(matZ) <- colnames(matCounts)
    return(matZ)
}