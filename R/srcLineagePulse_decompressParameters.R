#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++#++++    Decompress parameters: Compute parameter values from model  ++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# EXAMPLE CODE:

# 1. How to extract parameteres gene-wise?
# Use decompressMeansByGene, decompressDispByGene and 
# decompressDropoutRateByGene. Note that decompressDropoutRateByGene
# performs looping over cell-wise models for you.
#vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
#  lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
#  vecInterval=NULL )
#vecDispParam <- decompressDispersions( vecDispModel=lsDispModel$matDispModel[i,],
#  lsDispModel=lsDispersionModel$lsDispModelGlobal,
#  vecInterval=NULL )
#vecDropoutParam <- decompressDropoutRateByGene( vecDropModel=lsDropModel$matDropModel,
#  vecMu=vecMuParam,
#  vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
      
# 2. How to extract parameters cell-wise?
# Use decompressMeansByGene and decompressDispByGene in loop over genes
# for a single cell and use decompressDropoutRateByCell.
#vecMuParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
#  decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,],
#    lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
#    vecInterval=cell)
#}))
#vecDispParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
#  decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
#    lsDispModelGlobal=lsDispModel$lsMuModelGlobal,
#    vecInterval=cell)
#}))
#vecDropParam <- decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[j,],
#  vecMu=vecMuParam,
#  matPiConstPredictors=lsDropModel$matPiConstPredictors )

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
#'    \itemize{
#'      \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'          "windows"} Name of the mean model
#'      \item scaNumCells: (scalar) [Default NA] Number of cells
#'          for which model is evaluated. Used for constant model.
#'      \item vecPseudotime: (numerical vector number of cells)
#'          [Default NA] Pseudotime coordinates of cells. Used for
#'          impulse model.
#'      \item vecindClusterAssign: (integer vector length number of
#'          cells) [Default NA] Index of cluster assigned to each cell.
#'          Used for clusters model.
#'    }
#' @param vecInterval: (integer vector length target cells) [Default NULL]
#'    Positions of cells in ordering, for which parameters are to be 
#'    computed. Default: all cells.
#'    
#' @return vecMu: (numerical vector number of cells)
#'    Mean parameter estimates for given gene given the mean model.
#' @export

decompressMeansByGene <- function(vecMuModel,
  lsMuModelGlobal,
  vecInterval=NULL ){
  
  # Set interval to entire gene if not given
  # Dont do this to save time at the moment (dont reorder things
  # that are already ordered correctly)
  #if(is.null(vecInterval)){ vecInterval <- seq(1, lsMuModelGlobal$scaNumCells) }
  
  if(lsMuModelGlobal$strMuModel=="constant"){
    if(!is.null(vecInterval)){ 
      vecMu <- rep(vecMuModel, length(vecInterval))
    } else { 
      vecMu <- rep(vecMuModel, lsMuModelGlobal$scaNumCells) 
    }
  } else if(lsMuModelGlobal$strMuModel=="impulse"){
    if(!is.null(vecInterval)){ 
      vecMu <- calcImpulse_comp(vecTheta=vecMuModel, 
        vecTimepoints=lsMuModelGlobal$vecPseudotime[vecInterval])
    } else { 
      vecMu <- calcImpulse_comp(vecTheta=vecMuModel, 
        vecTimepoints=lsMuModelGlobal$vecPseudotime) 
    }
  } else if(lsMuModelGlobal$strMuModel=="clusters"){
    if(!is.null(vecInterval)){ 
      vecMu <- vecMuModel[lsMuModelGlobal$vecindClusterAssign[vecInterval]]
    } else { 
      vecMu <- vecMuModel[lsMuModelGlobal$vecindClusterAssign]
    }
  } else  if(lsMuModelGlobal$strMuModel=="windows"){
    if(!is.null(vecInterval)){ 
      vecMu <- vecMuModel[vecInterval]
    } else { 
      vecMu <- vecMuModel
    }
  } else {
    stop(paste0("ERROR decompressMeans(): lsMuModelGlobal$strMuModel=", 
      lsMuModelGlobal$strMuModel, " not recognised."))
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
#'    \itemize{
#'      \item strDispModel: (str) {"constant"} 
#'    Name of the dispersion model
#'      \item scaNumCells: (scalar) [Default NULL] Number of cells
#'    for which model is evaluated. Used for constant model.
#'      \item vecPseudotime: (numerical vector number of cells)
#'    [Default NULL] Pseudotime coordinates of cells.
#'      \item vecindClusterAssign: (integer vector length number of
#'    cells) [Default NULL] Index of cluster assigned to each cell.
#'    }
#' @param vecInterval: (integer vector length target cells) [Default NULL]
#'    Positions of cells in ordering, for which parameters are to be 
#'    computed. Default: all cells.
#'    
#' @return vecDisp: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#' @export

decompressDispByGene <- function(vecDispModel,
  lsDispModelGlobal,
  vecInterval=NULL ){

  if(lsDispModelGlobal$strDispModel=="constant"){
    if(!is.null(vecInterval)){ scaReps <- length(vecInterval)
    } else { scaReps <- lsDispModelGlobal$scaNumCells }
    vecDisp <- rep(vecDispModel, scaReps)
  } else {
    stop(paste0("ERROR decompressDispersions(): lsDispModelGlobal$strDispModel=", 
      lsDispModelGlobal$strDispModel, " not recognised."))
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
#'    {offset parameter, log(mu) paramter, parameters belonging to
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
#' @export

decompressDropoutRateByGene <- function(matDropModel,
  vecMu,
  vecPiConstPredictors ){
  
  vecPi <- sapply(seq(1,length(vecMu)), function(j){
    1/(1+exp(-( matDropModel[j,] %*% c(1, log(vecMu[j]), vecPiConstPredictors) )))
  })
  
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
#' @export

decompressDropoutRateByCell <- function(vecDropModel,
  vecMu,
  matPiConstPredictors ){
  
  vecPi <- sapply(seq(1,length(vecMu)), function(i){
    1/(1+exp(-( vecDropModel %*% c(1, log(vecMu[i]), matPiConstPredictors[i,]) )))
  })
  
  return(vecPi)
}