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
    
    scaNumSamples <- length(vecMu)
    
    # Compute probability of zero counts under 
    # negative binomial model.
    vecNBZero <- (vecDisp/(vecDisp+vecMu))^vecDisp
    # Compute posterior of drop-out.
    vecZ <- sapply(seq(1,scaNumSamples), function(j){
        if(vecboolNotZero[j]){
            scaZ <- 0
        } else if(vecboolZero[j]) {
            scaZ <- sum(vecDrop[j]/(vecDrop[j] + 
                                        (1-vecDrop[j])*vecNBZero[j])) *
                1/length(j)
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
    
    scaNumGenes <- dim(matCounts)[1]
    scaNumCells <- dim(matCounts)[2]
    
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
        vecNBZero <- (vecDispParam/(vecDispParam+vecMuParam))^vecDispParam
        
        vecZ <- sapply(seq(1,scaNumCells), function(j){
            if(!is.na(vecCounts[j]) & vecCounts[j]>0){
                scaZ <- 0
            } else if(!is.na(vecCounts[j]) & vecCounts[j]==0) {
                scaZ <- sum(
                    vecPiParam[j] / 
                        (vecPiParam[j] + (1-vecPiParam[j])*vecNBZero[j])) *
                    1/length(j)
            } else {
                scaZ <- NA
            }
            return(scaZ)
        })
        return(vecZ)
    }))
    
    rownames(matZ) <- vecGenes
    colnames(matZ) <- colnames(matCounts)
    return(matZ)
}