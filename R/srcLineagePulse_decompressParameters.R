#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++#++++    Decompress parameters: Compute parameter values from model  ++++++#
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
#' @param strMuModel: (str) {"constant", "impulse", "clusters", 
#'    "windows"} Name of the mean model
#' @param scaNumCells: (scalar) [Default NULL] Number of cells
#'    for which model is evaluated. Used for constant model.
#' @param vecPseudotime: (numerical vector number of cells)
#'    [Default NULL] Pseudotime coordinates of cells. Used for
#'    impulse model.
#' @param vecindClusterAssign: (integer vector length number of
#'    cells) [Default NULL] Index of cluster assigned to each cell.
#'    Used for clusters model.
#'    
#' @return vecMu: (numerical vector number of cells)
#'    Mean parameter estimates for given gene given the mean model.
#' @export

decompressMeans <- function(vecMuModel,
  strMuModel,
  scaNumCells=NULL,
  vecPseudotime=NULL,
  vecindClusterAssign=NULL ){
  
  if(strMuModel=="constant"){
    vecMu <- rep(vecMuModel, scaNumCells)
  } else if(strMuModel=="impulse"){
    vecMu <- calcImpulse_comp(vecMuModel, vecPseudotime)
  } else if(strMuModel=="clusters"){
    vecMu <- vecMuModel[vecindClusterAssign]
  } else  if(strMuModel=="windows"){
    vecMu <- vecMuModel
  } else {
    stop(paste0("ERROR decompressMeans(): strMuModel=", strMuModel, " not recognised."))
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
#' @param strDispModel: (str) {"constant"} 
#'    Name of the dispersion model
#' @param scaNumCells: (scalar) [Default NULL] Number of cells
#'    for which model is evaluated. Used for constant model.
#' @param vecPseudotime: (numerical vector number of cells)
#'    [Default NULL] Pseudotime coordinates of cells.
#' @param vecindClusterAssign: (integer vector length number of
#'    cells) [Default NULL] Index of cluster assigned to each cell.
#'    
#' @return vecDisp: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#' @export

decompressDispersions <- function(vecDispModel,
  strDispModel,
  scaNumCells=NULL,
  vecPseudotime=NULL,
  vecindClusterAssign=NULL ){

  if(strDispModel=="constant"){
    vecDisp <- rep(vecDispModel, scaNumCells)
  } else {
    stop(paste0("ERROR decompressDispersions(): strDispModel=", strDispModel, " not recognised."))
  }

  return(vecDisp)
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
#'    Parameters of dropout rate model for given gene.
#' @param vecMu: (numerical vector number of genes)
#'    Mean parameter estimates of all genes for given cell.
#' @param vecPiConstPredictors: (numerical vector number of 
#'    constant model predictors) Other model predictors than offset
#'    and the dynamically changing mean parameter. Examples are GC-
#'    content and other gene-specific properties.
#'    
#' @return vecPi: (numerical vector number of cells)
#'    Dispersion parameter estimates for given gene 
#'    (one per cell for given gene).
#' @export

decompressDropoutRate <- function(vecPiModel,
  vecMu,
  vecPiConstPredictors ){
  
  vecPi <- sapply(vecMu, function(mu_i){
    1/(1+exp(-( vecPiModel %*% c(1, log(mu_i), vecPiConstPredictors) )))
  })
  
  return(vecPi)
}