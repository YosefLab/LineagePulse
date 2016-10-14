#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++    Cluster expression mean trajectories  ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cluster expression mean trajectories
#' 
#' Cluster expression mean trajectories of given set of genes.
#' 
#' @seealso Called by \code{LineagePulse_main} or by user.
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

clusterGeneTrajectories <- function(vecIDs,
  lsMuModel){
  
}