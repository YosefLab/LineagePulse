#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++    Cluster expression mean trajectories  ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cluster expression mean trajectories
#' 
#' Cluster expression mean trajectories of given set of genes. Clustering
#' has to operate on full mean parameter matrix (genes x cells) 
#' and may therefore be memory intensive. Use with care. Clusters based on
#' correlation of mean trajectories.
#' 
#' @seealso Called by \code{LineagePulse_main} or by user.
#'
#' @param vecIDs: (vector of strings)
#'    Names of genes to cluster.
#' @param lsMuModel: (list length 2)
#'    All objects necessary to compute mean parameters for all
#'    observations.
#'    \itemize{
#'      \item matMuModel: (numerical matrix genes x number of model parameters)
#'    Parameters of mean model for each gene.
#'      \item lsMuModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'        \itemize{
#'          \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'        "windows"} Name of the mean model.
#'          \item scaNumCells: (scalar) [Default NA] Number of cells
#'        for which model is evaluated. Used for constant model.
#'          \item vecPseudotime: (numerical vector number of cells)
#'        [Default NA] Pseudotime coordinates of cells. Used for
#'        impulse model.
#'          \item vecindClusterAssign: (integer vector length number of
#'        cells) [Default NA] Index of cluster assigned to each cell.
#'        Used for clusters model.
#'          \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'        of a gene are simultaneously estiamted as a vector with BFGS
#'        in windows mode.
#'          \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'        for BFGS estimation of impulse model with optim (termination criterium).
#'          \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'        change in objective function for BFGS estimation of impulse 
#'        model with optim (termination criterium).
#'      }
#'    }
#'    
#' @return vecidxClusters: (int vector number of IDs)
#'    Index of cluster for each given ID.
#' @export

clusterGeneTrajectories <- function(vecIDs,
  lsMuModel){
  
  # Check IDs are in model matrix
  if(any(!(vecIDs %in% rownames(lsMuModel$matMuModel)))){
    print("ERROR: Some IDs in vecIDs are not given in model lsMuModel$matMuModel.")
    return(NULL)
  }
  
  # Decompress mean parameters in a memory saving way:
  # Round to integer counts.
  matMuParam <- do.call(rbind, bplapply(vecIDs, function(id){
    round( decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[id,],
      lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
      vecInterval=NULL ) )
  }))
  
  # Cluster genes: Hierarchical cluster based on correlation.
  
  return(vecidxClusters)
  
}