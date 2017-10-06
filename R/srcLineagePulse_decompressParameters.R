#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++    Decompress parameters: Compute parameter values from model  ++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute mean parameter estimates from mean parameter model for a gene
#' 
#' Takes the model type and computes one mean parameter for each cell for one
#' gene.
#' 
#' @seealso Called by \code{fitZINB}.
#'
#' @param vecMuModel: (numerical vector number of model parameters)
#'    Parameters of mean model for given gene.
#' @param lsMuModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#' @param vecInterval: (integer vector length target cells) [Default NULL]
#'    Positions of cells in ordering, for which parameters are to be 
#'    computed. Default: all cells.
#'    
#' @return vecMu: (numerical vector number of cells)
#'    Mean parameter estimates for given gene given the mean model.
#'    
#' @author David Sebastian Fischer
decompressMeansByGene <- function(
    vecMuModel,
    lsvecBatchModel=NULL,
    lsMuModelGlobal,
    vecInterval=NULL ){
    
    if(lsMuModelGlobal$strMuModel=="constant"){
        if(!is.null(vecInterval)){ 
            vecMu <- rep(vecMuModel, length(vecInterval))
        } else { 
            vecMu <- rep(vecMuModel, lsMuModelGlobal$scaNumCells) 
        }
    } else if(lsMuModelGlobal$strMuModel=="impulse"){
        if(!is.null(vecInterval)){
            vecMu <- evalImpulseModel_comp(
                vecImpulseParam=vecMuModel, 
                vecTimepoints=lsMuModelGlobal$vecPseudotime[vecInterval])
        } else { 
            vecMu <- evalImpulseModel_comp(
                vecImpulseParam=vecMuModel, 
                vecTimepoints=lsMuModelGlobal$vecPseudotime) 
        }
    } else if(lsMuModelGlobal$strMuModel=="splines"){
        if(!is.null(vecInterval)){
            vecMu <- exp(as.vector(lsMuModelGlobal$matSplineBasis[
                vecInterval,,drop=FALSE] %*% vecMuModel))
        } else { 
            vecMu <- exp(as.vector(lsMuModelGlobal$matSplineBasis %*% vecMuModel))
        }
    } else if(lsMuModelGlobal$strMuModel=="groups"){
        if(!is.null(vecInterval)){ 
            vecMu <- vecMuModel[lsMuModelGlobal$vecidxGroups[vecInterval]]
        } else { 
            vecMu <- vecMuModel[lsMuModelGlobal$vecidxGroups]
        }
    } else  if(lsMuModelGlobal$strMuModel=="windows"){
        if(!is.null(vecInterval)){ 
            vecMu <- vecMuModel[vecInterval]
        } else { 
            vecMu <- vecMuModel
        }
    } else if(lsMuModelGlobal$strMuModel=="MM"){
        # Expect mean of ONE MIXTURE component! vecMuModel is scalar
        if(!is.null(vecInterval)){ 
            vecMu <- rep(vecMuModel, length(vecInterval))
        } else { 
            vecMu <- rep(vecMuModel, lsMuModelGlobal$scaNumCells)
        }
    } else {
        stop(paste0("ERROR decompressMeans(): lsMuModelGlobal$strMuModel=", 
                    lsMuModelGlobal$strMuModel, " not recognised."))
    }
    # Scale by batch factors
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        for(conf in seq(1,lsMuModelGlobal$scaNConfounders, by=1)){
            if(!is.null(vecInterval)){
                vecMu <- vecMu*(lsvecBatchModel[[conf]][
                    (lsMuModelGlobal$lsvecidxBatchAssign[[conf]])[vecInterval]])
            } else { 
                vecMu <- vecMu*(lsvecBatchModel[[conf]][
                    lsMuModelGlobal$lsvecidxBatchAssign[[conf]]])
            }
        }
    }
    
    return(vecMu)
}

#' Compute dispersion parameter estimates from mean parameter model for a gene
#' 
#' Takes the model type and computes one dispersion parameter for each cell for one
#' gene.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecDispModel: (numerical vector number of model parameters)
#'    Parameters of dispersion model for given gene.
#' @param lsDispModelGlobal: (list) Global variables for dispersion model,
#'    common to all genes.
#' @param vecInterval: (integer vector length target cells) [Default NULL]
#'    Positions of cells in ordering, for which parameters are to be 
#'    computed. Default: all cells.
#'    
#' @return vecDisp: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#'    
#' @author David Sebastian Fischer
decompressDispByGene <- function(vecDispModel,
                                 lsvecBatchModel=NULL,
                                 lsDispModelGlobal,
                                 vecInterval=NULL ){
    
    if(lsDispModelGlobal$strDispModel=="constant"){
        if(!is.null(vecInterval)) { 
            vecDisp <- rep(vecDispModel, length(vecInterval))
        } else { 
            vecDisp <- rep(vecDispModel, lsDispModelGlobal$scaNumCells)
        }
    } else if(lsDispModelGlobal$strDispModel=="splines"){
        if(!is.null(vecInterval)){
            vecDisp <- exp(as.vector(lsDispModelGlobal$matSplineBasis[
                vecInterval,,drop=FALSE] %*% vecDispModel))
        } else { 
            vecDisp <- exp(as.vector(lsDispModelGlobal$matSplineBasis %*% vecDispModel))
        }
    } else if(lsDispModelGlobal$strDispModel=="groups"){
        if(!is.null(vecInterval)){ 
            vecDisp <- vecDispModel[lsDispModelGlobal$vecidxGroups[vecInterval]]
        } else { 
            vecDisp <- vecDispModel[lsDispModelGlobal$vecidxGroups]
        }
    } else if(lsDispModelGlobal$strDispModel=="MM"){
        # Expect mean of ONE MIXTURE component! vecDispModel is scalar
        if(!is.null(vecInterval)){ 
            vecDisp <- rep(vecDispModel, length(vecInterval))
        } else { 
            vecDisp <- rep(vecDispModel, lsDispModelGlobal$scaNumCells)
        }
    } else {
        stop(paste0("ERROR decompressDispersions(): lsDispModelGlobal$strDispModel=", 
                    lsDispModelGlobal$strDispModel, " not recognised."))
    }
    
    # Scale by batch factors
    if(!is.null(lsDispModelGlobal$lsvecidxBatchAssign)){
        for(conf in seq(1,lsDispModelGlobal$scaNConfounders, by=1)){
            if(!is.null(vecInterval)){
                vecDisp <- vecDisp*(lsvecBatchModel[[conf]][(lsDispModelGlobal$lsvecidxBatchAssign[[conf]])[vecInterval]])
            } else { 
                vecDisp <- vecDisp*(lsvecBatchModel[[conf]][lsDispModelGlobal$lsvecidxBatchAssign[[conf]]])
            }
        }
    }
    
    return(vecDisp)
}

#' Compute dropout rate parameter estimates from dropout rate model for a gene
#' 
#' Compute dropout rate parameter estimates from dropout rate model for a gene.
#' 
#' @seealso Called by \code{fitZINB}.
#'
#' @param vecPiModel: (numerical matrix genes x number of model parameters)
#'    {offset parameter, log(mu) parameter, parameters belonging to
#'    constant predictors}
#'    Parameters of dropout rate model for all cells.
#' @param vecMu: (numerical vector number of genes)
#'    Mean parameter estimates of all cells for given gene.
#' @param vecPiConstPredictors: (numerical vector number of 
#'    constant model predictors) Other model predictors than offset
#'    and the dynamically changing mean parameter. Examples are GC-
#'    content and other gene-specific properties. This would be the 
#'    global parameters as listed in the other decompression
#'    function. Here those are not a list as there is only one object.
#'    
#' @return vecPi: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#'    
#' @author David Sebastian Fischer
#' 
#' @export
decompressDropoutRateByGene <- function(
    matDropModel,
    vecMu=NULL,
    vecPiConstPredictors=NULL,
    lsDropModelGlobal ){
    
    if(lsDropModelGlobal$strDropModel=="logistic"){
        vecPi <- sapply(seq(1,lsDropModelGlobal$scaNumCells), function(j){
            evalDropoutModel_comp(vecPiModel=matDropModel[j,], 
                                  vecPiPredictors=c(1, vecPiConstPredictors))
        })
    } else if(lsDropModelGlobal$strDropModel=="logistic_ofMu"){
        vecPi <- sapply(seq(1,lsDropModelGlobal$scaNumCells), function(j){
            evalDropoutModel_comp(vecPiModel=matDropModel[j,], 
                                  vecPiPredictors=c(1, log(vecMu[j]), vecPiConstPredictors))
        })
    } else {
        stop(paste0("ERROR IN decompressDropoutRateByGene: ",
                    "lsDropModelGlobal$strDropModel=",
                    lsDropModelGlobal$strDropModel,
                    " not recognised."))
    }
    
    return(vecPi)
}

#' Compute dropout rate parameter estimates from dropout rate model for a cell
#' 
#' Compute dropout rate parameter estimates from dropout rate model for a cell.
#' 
#' @seealso Called by \code{fitZINB}.
#'
#' @param vecPiModel: (numerical vector number of model parameters)
#'    {offset parameter, log(mu) paramter, parameters belonging to
#'    constant predictors}
#'    Parameters of dropout rate model for given cell.
#' @param vecMu: (numerical vector number of genes)
#'    Mean parameter estimates of all genes for given cell.
#' @param matPiConstPredictors: (numerical matrix genes x number of 
#'    constant model predictors) Other model predictors than offset
#'    and the dynamically changing mean parameter. Examples are GC-
#'    content and other gene-specific properties. This would be the 
#'    global parameters as listed in the other decompression
#'    function. Here those are not a list as there is only one object.
#'    
#' @return vecPi: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#'    
#' @author David Sebastian Fischer
decompressDropoutRateByCell <- function(
    vecDropModel,
    vecMu=NULL,
    matPiConstPredictors=NULL,
    lsDropModelGlobal){
    
    if(lsDropModelGlobal$strDropModel=="logistic"){
        vecPi <- sapply(seq(1,lsDropModelGlobal$scaNumGenes), function(i){
            evalDropoutModel_comp(vecPiModel=vecDropModel, 
                                  vecPiPredictors=c(1, matPiConstPredictors[i,]))
        })
    } else if(lsDropModelGlobal$strDropModel=="logistic_ofMu"){
        vecPi <- sapply(seq(1,lsDropModelGlobal$scaNumGenes), function(i){
            evalDropoutModel_comp(vecPiModel=vecDropModel, 
                                  vecPiPredictors=c(1, log(vecMu[i]), matPiConstPredictors[i,]))
        })
    } else {
        stop(paste0("ERROR IN decompressDropoutRateByCell: ",
                    "lsDropModelGlobal$strDropModel=",
                    lsDropModelGlobal$strDropModel,
                    " not recognised."))
    }
    
    return(vecPi)
}