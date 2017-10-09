#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++  Extract model fits and transformed data from fits  +++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Return depth and batch corrected data
#' 
#' The data normalisation is based on the model normalisation used by and inferred
#' by LineagePulse, e.g. for data visualisation.
#' 
#' @seealso Called by \code{fitZINB}. Can be called by user.
#' 
#' @param matCounts (numeric matrix genes x cells)
#' Count data.
#' @param lsMuModel (list) Mean parameter model parameters.
#' @param vecGeneIDs (vector of strings) 
#' Gene IDs for which mean model fits are to be extracted.
#' @param boolDepth (bool) [Default TRUE] Whether to normalize for sequencing depth.
#' @param boolBatch (bool) [Default TRUE] Whether to normalize for batch.
#' 
#' @return (numeric matrix genes x cells)
#' Input data normalized by library size factors (optional) and
#' by inferred batch factors (optional).
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 20,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # Get batch correction on alternative model:
#' # Use H1 model fits.
#' matNormData <- getNormData(
#'      matCounts = lsSimulatedData$counts,
#'      lsMuModel = objLP@lsMuModelH1,
#'      vecGeneIDs = rownames(lsSimulatedData$counts)[1],
#'      boolDepth = TRUE, boolBatch = TRUE)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getNormData <- function(matCounts,
                        lsMuModel,
                        vecGeneIDs,
                        boolDepth = TRUE,
                        boolBatch = TRUE){
    
    ### Check input
    if(!boolDepth & !boolBatch) {
        warning("No normalisation chosen")
    }
    
    if(!all(vecGeneIDs %in% rownames(lsMuModel$matMuModel))) {
        stop(paste0("ERROR getNormData(): Not all vecGeneIDs were ",
                    "rownames of lsMuModel$matMuModel"))
    }
    
    scaNumCells <- dim(matCounts)[2]
    scaNumConfounders <- length(lsMuModel$lsmatBatchModel)
    matNormData <- do.call(rbind, bplapply(vecGeneIDs, function(id){
        # Extract batch parameters
        vecBatchParam <- array(1, scaNumCells)
        if(!is.null(lsMuModel$lsMuModelGlobal$vecConfounders)){
            for(confounder in seq(1,scaNumConfounders)){
                vecBatchParam <- vecBatchParam * 
                    (lsMuModel$lsmatBatchModel[[confounder]][id,][
                    lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign[[confounder]]])
            }
        }
        # Normalise counts by depth and/or batch factors as required
        vecNormData <- matCounts[id,]
        if(boolDepth) {
            vecNormData <- vecNormData/lsMuModel$lsMuModelGlobal$vecNormConst
        }
        if(boolBatch) {
            vecNormData <- vecNormData/vecBatchParam
        }
        return(vecNormData)
    }))
    
    rownames(matNormData) <- vecGeneIDs
    colnames(matNormData) <- colnames(matCounts)
    return(matNormData)
}

#' Get mean model fits
#' 
#' Return mean model fits per gene and cell as matrix for chosen model.
#' 
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param vecGeneIDs (vector of strings) 
#' Gene IDs for which mean model fits are to be extracted.
#' 
#' @return (numeric matrix genes x cells)
#' Mean parameter fits.
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 20,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # Get mean parameter fits on alternative model:
#' # Use H1 model fits.
#' vecMeanFits <- getFitsMean(
#'      lsMuModel = objLP@lsMuModelH1,
#'      vecGeneIDs = rownames(lsSimulatedData$counts)[1])
#' #plot(lsSimulatedData$annot$pseudotime, vecMeanFits)     
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getFitsMean <- function(
    lsMuModel,
    vecGeneIDs=NULL){
    
    ### Check input
    if(!all(vecGeneIDs %in% rownames(lsMuModel$matMuModel))) {
        stop(paste0("ERROR getFitsDropout(): Not all vecGeneIDs ",
                    "were rownames of lsMuModel$matMuModel"))
    }
    
    ### Decompress models into parameter estimates
    matMuParam <- do.call(rbind, lapply(vecGeneIDs, function(i){
        # Decompress parameters by gene
        vecMuParam <- decompressMeansByGene(
            vecMuModel=lsMuModel$matMuModel[i,],
            lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
            lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
            vecInterval=NULL )
        return(vecMuParam)
    }))
    
    rownames(matMuParam) <- vecGeneIDs
    return(matMuParam)
}

#' Get dispersion model fits
#' 
#' Return dispersion model fits per gene and cell as matrix for chosen model.
#' 
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param vecGeneIDs (vector of strings) 
#' Gene IDs for which dispersion model fits are to be extracted.
#' 
#' @return (numeric matrix genes x cells)
#' Dispersion parameter fits.
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 20,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # Get dispersion parameter fits on alternative model:
#' # Use H1 model fits.
#' vecDispersionFits <- getFitsDispersion(
#'      lsDispModel = objLP@lsDispModelH1,
#'      vecGeneIDs = rownames(lsSimulatedData$counts)[1])
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getFitsDispersion <- function(
    lsDispModel,
    vecGeneIDs=NULL){
    
    ### Check input
    if(!all(vecGeneIDs %in% rownames(lsDispModel$matDispModel))) {
        stop(paste0("ERROR getFitsDropout(): Not all vecGeneIDs were ",
                    "rownames of lsDispModel$matDispModel"))
    }
    
    ### Decompress models into parameter estimates
    matDispParam <- do.call(rbind, lapply(vecGeneIDs, function(i){
        # Decompress parameters by gene
        vecDispParam <- decompressDispByGene(
            vecDispModel=lsDispModel$matDispModel[i,],
            lsvecBatchModel=lapply(lsDispModel$lsmatBatchModel, function(mat) mat[i,] ),
            lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
            vecInterval=NULL )
        return(vecDispParam)
    }))
    
    rownames(matDispParam) <- vecGeneIDs
    return(matDispParam)
}

#' Get drop-out model fits
#' 
#' Return drop-out model fits per gene and cell as matrix for chosen models.
#' 
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' @param vecIDs (vector of strings) [Default NULL]
#' Gene IDs for which posteriors of drop-out are to be computed.
#' 
#' @return (numeric matrix genes x cells)
#' Drop-out rate fits.
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 20,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # Get drop-out rate fits on alternative model:
#' # Use H1 model fits.
#' vecDropoutFits <- getFitsDropout(
#'      lsMuModel = objLP@lsMuModelH1,
#'      lsDropModel = objLP@lsDropModel,
#'      vecGeneIDs = rownames(lsSimulatedData$counts)[1])
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getFitsDropout <- function(
    lsMuModel,
    lsDropModel,
    vecGeneIDs=NULL){
    
    ### Check input
    if(!all(vecGeneIDs %in% rownames(lsMuModel$matMuModel))) {
        stop(paste0("ERROR getFitsDropout(): Not all vecGeneIDs were ",
                    "rownames of lsMuModel$matMuModel"))
    }
    
    ### Decompress models into parameter estimates
    matPiParam <- do.call(rbind, lapply(vecGeneIDs, function(i){
        # Decompress parameters by gene
        vecMuParam <- decompressMeansByGene(
            vecMuModel=lsMuModel$matMuModel[i,],
            lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
            lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
            vecInterval=NULL )
        vecPiParam <- decompressDropoutRateByGene(
            matDropModel=lsDropModel$matDropoutLinModel,
            vecMu=vecMuParam,
            vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
            lsDropModelGlobal=lsDropModel$lsDropModelGlobal)
        return(vecPiParam)
    }))
    
    rownames(matPiParam) <- vecGeneIDs
    return(matPiParam)
}

#' Get posteriors of drop-out
#' 
#' Return posteriors of drop-out per gene and cell as matrix for chosen models.
#' 
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' @param vecGeneIDs (vector of strings) [Default NULL]
#' Gene IDs for which posteriors of drop-out are to be computed.
#' 
#' @return (numeric matrix genes x cells)
#' Posterior probability of observation not being generated 
#' by drop-out.
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 20,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # Get posterior of drop-out on alternative model:
#' # Use H1 model fits.
#' vecPosteriorDropoutFits <- getPostDrop(
#'      matCounts = lsSimulatedData$counts,
#'      lsMuModel = objLP@lsMuModelH1,
#'      lsDispModel = objLP@lsDispModelH1,
#'      lsDropModel = objLP@lsDropModel,
#'      vecGeneIDs = rownames(lsSimulatedData$counts)[1])
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getPostDrop <- function(
    matCounts,
    lsMuModel,
    lsDispModel, 
    lsDropModel,
    vecGeneIDs){
    
    ### Check input
    if(!all(vecGeneIDs %in% rownames(matCounts))) {
        stop(paste0("ERROR getFitsDropout(): Not all vecGeneIDs were ",
                    "rownames of matCounts."))
    }
    
    ### Decompress models into parameter estimates
    matPostDrop <- calcPostDrop_Matrix(
        matCounts = matCounts,
        lsMuModel = lsMuModel,
        lsDispModel = lsDispModel, 
        lsDropModel = lsDropModel,
        vecIDs = vecGeneIDs)
    
    rownames(matPostDrop) <- vecGeneIDs
    return(matPostDrop)
}