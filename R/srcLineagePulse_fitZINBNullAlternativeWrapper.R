#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# HOWTO debug this: 
# Set bplapply of estimation that throws error to lapply for proper error
# reporting. Or - give nProc=1 and BiocParallel operates in SerialMode with
# proper error reporting -> all estiamtions are sequential though.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit zero-inflated negative binomial model to data
#' 
#' This is the algorithmic core wrapper of LineagePulse that carries out
#' the entire parameter estimation for both H0 and H1.
#'
#' Fit alternative H1 and null H0 zero-inflated negative binomial model 
#' to a data set using cycles of coordinate ascent. The algorithm first
#' fits the either H1 or H0 together with the logistic dropout model by
#' iterating over cell-wise (dropout models) and gene-wise (negative 
#' binomial models) parameters. Subsequently, the remaining model (H0 
#' or H1) is estimated by iterating over zero-inflated negative binomial
#' mean and dispersion parameter estimation condition on the previously
#' estimated logistic drop-out model.
#' 
#' Estimation of H0 and H1 are therefore separate coordinate ascent 
#' procedures yielding different local optima on the overall zero-inflated
#' negative binomial loglikelihood function with different gene-wise constraints.
#' 
#' Convergence is tracked with the the loglikelihood of the entire 
#' data matrix. Every step is a maximum likelihood estimation of the 
#' target parameters conditioned on the remaining parameter estimates. 
#' Therefore, convergence to a local optimum is guaranteed if the algorithm
#' is run until convergence. Parallelisation of each estimation step 
#' is implemented where conditional independences of parameter estimations
#' allow so. 
#' 
#' Convergence can be followed with verbose=TRUE (at each 
#' iteration) or at each step (boolSuperVerbose=TRUE). Variables for the
#' logistic drop-out model are a constant and the estimated mean parameter
#' and other constant gene-specific variables (such as GC-conten) in 
#' matPiConstPredictors. Three modes are available for modelling the mean
#' parameter: As a gene-wise constant (the default null model), by cluster 
#' (this is fast as neighbourhoods don't have to be evaluated), 
#' sliding windows (using neighbourhood smoothing), and as an impulse model.
#' 
#' To save memory, not the entire parameter matrix (genes x cells) but
#' the parmater models are stored in the objects lsMuModel, lsDispModel
#' and lsDropModel. These objects are described in detail in the annotation
#' of the return values of the function. In short, these object contain
#' the gene/cell-wise parameters of the model used to constrain the parameter
#' in question and the predictors necessary to evaluate the parameter model
#' to receive the observation-wise paramter values. Example: Impulse model
#' for the mean parameter: lsMuModel contains the parameter estimates for an
#' impulse model for each gene and pseudotime coordinates. Therefore, the
#' mean parameter for each observation can be computed as the value of the
#' impulse model evaluated at the pseudotime points for each gene.
#' 
#' @seealso Called by \code{runLineagePulse}. Calls parameter estimation
#' wrappers:
#' \code{fitPiZINB}, \code{fitZINBMu}, \code{fitZINBDisp} and
#' \code{fitZINBMuDisp}.
#' Calls \code{evalLogLikMatrix} to follow convergence.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param boolEstimateNoiseBasedOnH0: (bool) [Default: FALSE]
#'    Whether to co-estimate logistic drop-out model with the 
#'    constant null model or with the alternative model. The
#'    co-estimation with the noise model typically extends the
#'    run-time of this model-estimation step strongly. While
#'    the drop-out model is more accurate if estimated based on
#'    a more realistic model expression model (the alternative
#'    model), a trade-off for speed over accuracy can be taken
#'    and the dropout model can be chosen to be estimated based
#'    on the constant null expression model (set to TRUE).
#' @param boolVecWindowsAsBFGS: (bool) [Default FALSE] Whether
#'    mean parameters of a gene are co-estimated in "windows"
#'    mode with BFGS algorithm (optimisation with dimensionality
#'    of number of cells) or estimated one by one, conditioned
#'    one the latest estimates of neighbours. The latter case
#'    (boolVecWindowsAsBFGS=FALSE) is coordinate ascent within the gene
#'    and each mean parameter is optimised once only.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells. Only required if mean model
#'    or dispersion model are fit as a function of pseudotime, 
#'    e.g. impulse model for means.
#' @param scaMaxEstimationCycles: (integer) [Default 20] Maximium number 
#'    of estimation cycles performed in fitZINB(). One cycle
#'    contain one estimation of of each parameter of the 
#'    zero-inflated negative binomial model as coordinate ascent.
#' @param verbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation in high detail with local 
#'    convergence flags and step-by-step loglikelihood computation.
#' 
#' @return (list length 6)
#'    \itemize{
#'      \item lsMuModelH1: (list length 2)
#'    All objects necessary to compute H1 mean parameters for all
#'    observations.
#'      \itemize{
#'        \item matMuModel: (numerical matrix genes x number of model parameters)
#'      Parameters of mean model for each gene.
#'        \item lsMuModelGlobal: (list) Global variables for mean model,
#'      common to all genes.
#'      \itemize{
#'        \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'      "windows"} Name of the mean model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'        \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'      of a gene are simultaneously estiamted as a vector with BFGS
#'      in windows mode.
#'        \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'      for BFGS estimation of impulse model with optim (termination criterium).
#'        \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'      change in objective function for BFGS estimation of impulse 
#'      model with optim (termination criterium).
#'      }
#'    }
#'    \item lsDispModelH1: (list length 2)
#'    All objects necessary to compute H1 dispersion parameters for all
#'    observations.
#'    \itemize{
#'      \item matDispModel: (numerical matrix genes x number of model parameters)
#'    Parameters of dispersion model for each gene.
#'      \item lsDispModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'      \itemize{
#'        \item strDispModel: (str) {"constant"} 
#'      Name of the dispersion model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'      }
#'    }
#'      \item lsMuModelH0: (list length 2)
#'    All objects necessary to compute H0 mean parameters for all
#'    observations.
#'      \itemize{
#'        \item matMuModel: (numerical matrix genes x number of model parameters)
#'      Parameters of mean model for each gene.
#'        \item lsMuModelGlobal: (list) Global variables for mean model,
#'      common to all genes.
#'      \itemize{
#'        \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'      "windows"} Name of the mean model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'        \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'      of a gene are simultaneously estiamted as a vector with BFGS
#'      in windows mode.
#'        \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'      for BFGS estimation of impulse model with optim (termination criterium).
#'        \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'      change in objective function for BFGS estimation of impulse 
#'      model with optim (termination criterium).
#'      }
#'    }
#'    \item lsDispModelH0: (list length 2)
#'    All objects necessary to compute H0 dispersion parameters for all
#'    observations.
#'    \itemize{
#'      \item matDispModel: (numerical matrix genes x number of model parameters)
#'    Parameters of dispersion model for each gene.
#'      \item lsDispModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'      \itemize{
#'        \item strDispModel: (str) {"constant"} 
#'      Name of the dispersion model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'      }
#'    }
#'    \item lsDropModel: (list length 2)
#'    All objects necessary to compute drop-out parameters for all
#'    observations, omitting mean parameters (which are stored in lsMeanModel).
#'      \itemize{
#'        \item matDropoutLinModel: (numeric matrix cells x number of model parameters)
#'      {offset parameter, log(mu) parameter, parameters belonging to
#'      constant predictors}
#'      Parameters of drop-out model for each cell
#'        \item matPiConstPredictors: (numeric matrix genes x number of constant
#'      gene-wise drop-out predictors) Predictors for logistic drop-out 
#'      fit other than offset and mean parameter (i.e. parameters which
#'      are constant for all observations in a gene and externally supplied.)
#'      Is null if no constant predictors are supplied.
#'      }
#'    \item lsFitZINBReporters: (list length 6)
#'    Reporters of behaviour of overall fitting procedure.
#'      \itemize{
#'        \item boolConvergenceH1: (bool) Convergence of
#'      estimation for alternative model H1. Convergence
#'      is evaluated based on the convergence of the 
#'      loglikelihood of the entire data set.
#'        \item boolConvergenceH0: (bool) Convergence of
#'      estimation for null model H0. Convergence
#'      is evaluated based on the convergence of the 
#'      loglikelihood of the entire data set.
#'        \item vecEMLogLikH1: (numeric vector number of 
#'      estimation cycles) Loglikelihood of entire 
#'      data set after each estimation cycle of alternative
#'      model H1.
#'        \item vecEMLogLikH0: (numeric vector number of 
#'      estimation cycles) Loglikelihood of entire 
#'      data set after each estimation cycle of null
#'      model H0.
#'        \item scaKbyGeneH1: (scalar) Degrees of freedom
#'      by gene used in alternative model H1. The logistic
#'      dropout model is ignored as it is shared between
#'      alternative and null model.
#'        \item scaKbyGeneH0: (scalar) Degrees of freedom
#'      by gene used in null model H0. The logistic
#'      dropout model is ignored as it is shared between
#'      alternative and null model.
#'      }
#'    }
#' @export

fitNullAlternative <- function(matCountsProc,
  matPiConstPredictors,
  lsResultsClustering,
  vecNormConst,
  scaWindowRadius=NULL,
  strMuModel="windows",
  strDispModel = "constant",
  boolEstimateNoiseBasedOnH0=TRUE,
  boolVecWindowsAsBFGS=FALSE,
  vecPseudotime=NULL,
  scaMaxEstimationCycles=20,
  boolVerbose=FALSE,
  boolSuperVerbose=FALSE ){
  
  ####################################################
  # Compute degrees of freedom of model for each gene
  # Drop-out model is ignored, would be counted at each gene.
  # 1. Alternative model  H1:
  # Mean model:
  if(strMuModel=="windows"){
    # One mean parameter per cell
    scaKbyGeneH1 <- dim(matCountsProc)[2]
  } else if(strMuModel=="clusters"){
    # One mean parameter per cluster
    scaKbyGeneH1 <- lsResultsClustering$K
  } else if(strMuModel=="MM"){
    # One mean parameter per mixture component
    scaKbyGeneH1 <- dim(matWeights)[2]
  } else if(strMuModel=="impulse"){
    # Six impulse model parameter to model means
    scaKbyGeneH1 <- 6
  } else if(strMuModel=="constant"){
    # One constant mean
    scaKbyGeneH1 <- 1
  } else {
    stop(paste0("ERROR in fitZINB(): strMuModel not recognised: ", strMuModel))
  }
  # Dispersion model
  if(strDispModel=="constant"){
    # One dispersion factor per gene.
    scaKbyGeneH1 <- scaKbyGeneH1 + 1
  } else {
    stop(paste0("ERROR in fitZINB(): strMuModel not recognised: ", strDispModel))
  }
  # 2. Null model H0:
  # One mean per gene and one dispersion factor
  scaKbyGeneH0 <- 1+1
  
  ####################################################
  # Set sequence of model estimation: Which mean model
  # (H0 or H1) is co-estimated with noise model and which
  # is estimated conditioned on the noise model of this estimation.
  if(boolEstimateNoiseBasedOnH0){
    strMuModelA <- "constant"
    strMuModelB <- strMuModel
    strDispModelA <- "constant"
    strDispModelB <- "constant"
    strNameModelA <- paste0("H0: ",strMuModelA)
    strNameModelB <- paste0("H1: ",strMuModelB)
  } else {
    strMuModelA <- strMuModel
    strMuModelB <- "constant"
    strDispModelA <- "constant"
    strDispModelB <- "constant"
    strNameModelA <- paste0("H1: ",strMuModelA)
    strNameModelB <- paste0("H0: ",strMuModelB)
  }
  
  ####################################################
  # Fit model A
  print(paste0("### a) Fit negative binomial model A (",
    strNameModelA,") with noise model."))
  
  tm_cycle <- system.time({
    lsFitsModelA <- fitZINB(matCounts=matCountsProc,
            lsMuModel=lsMuModelA,
            lsDispModel=lsDispModelA,
            lsDropModel=lsDropModel,
            boolFitDrop=TRUE,
            strMuModel=strMuModelA,
            strDispModel=strDispModelB,
            scaWindowRadius=scaWindowRadius,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
  })
  lsMuModelA <- lsFitsModelA$lsMuModel
  lsDispModelA <- lsFitsModelA$lsDispModel
  lsDropModel <- lsFitsModelA$lsDropModel
  boolConvergenceModelA <- lsFisModelA$boolConvergenceModel
  vecEMLogLikModelA <- lsFisModelA$vecEMLogLikModel
  
  print(paste0("Finished fitting zero-inflated negative binomial ",
               "model A with noise model in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  ####################################################
  # Fit model B
  print(paste0("### a) Fit negative binomial model B (",
               strNameModelB,") with noise model."))
  
  tm_cycleB <- system.time({
    lsFitsModelB <- fitZINB(matCounts=matCountsProc,
                            lsDropModel=lsDropModel,
                            boolFitDrop=FALSE,
                            strMuModel=strMuModelB,
                            strDispModel=strDispModelB,
                            scaWindowRadius=scaWindowRadius,
                            boolVerbose=boolVerbose,
                            boolSuperVerbose=boolSuperVerbose)
  })
  lsMuModelB <- lsFitsModelB$lsMuModel
  lsDispModelB <- lsFitsModelB$lsDispModel
  boolConvergenceModelB <- lsFisModelB$boolConvergenceModel
  vecEMLogLikModelB <- lsFisModelB$vecEMLogLikModel
  
  print(paste0("Finished fitting zero-inflated negative binomial ",
               "model B in ", round(tm_cycleB["elapsed"]/60,2)," min."))
  
  # Name rows and columns of parameter matrices
  rownames(lsMuModelA$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelA$matDispModel) <- rownames(matCountsProc)
  rownames(lsMuModelB$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelB$matDispModel) <- rownames(matCountsProc)
  rownames(lsDropModel$matDropoutLinModel) <- colnames(matCountsProc)
  if(!is.null(lsDropModel$matPiConstPredictors)){
    rownames(lsDropModel$matPiConstPredictors) <- rownames(matCountsProc)
  }
  
  if(boolEstimateNoiseBasedOnH0){
    lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelB,
                                boolConvergenceH0=boolConvergenceModelA,
                                vecEMLogLikH1=vecEMLogLikModelB,
                                vecEMLogLikH0=vecEMLogLikModelA,
                                scaKbyGeneH1=scaKbyGeneH1,
                                scaKbyGeneH0=scaKbyGeneH0 )
    objectLineagePulse <- new('LineagePulseObject',
                              dfResults           = NULL,
                              vecDEGenes          = NULL,
                              lsMuModelH1         = lsMuModelB,
                              lsDispModelH1       = lsDispModelB,
                              lsMuModelH0         = lsMuModelA,
                              lsDispModelH0       = lsDispModelA,
                              lsFitZINBReporters  = lsFitZINBReporters,
                              dfAnnotationProc    = NULL,
                              vecNormConst        = vecNormConst,
                              scaNProc            = NULL,
                              scaQThres           = NULL,
                              strReport           = NULL)
  } else {
    lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelA,
                                boolConvergenceH0=boolConvergenceModelB,
                                vecEMLogLikH1=vecEMLogLikModelA,
                                vecEMLogLikH0=vecEMLogLikModelB,
                                scaKbyGeneH1=scaKbyGeneH1,
                                scaKbyGeneH0=scaKbyGeneH0 )
    objectLineagePulse <- new('LineagePulseObject',
                              dfResults           = NULL,
                              vecDEGenes          = NULL,
                              lsMuModelH1         = lsMuModelA,
                              lsDispModelH1       = lsDispModelA,
                              lsMuModelH0         = lsMuModelB,
                              lsDispModelH0       = lsDispModelB,
                              lsFitZINBReporters  = lsFitZINBReporters,
                              dfAnnotationProc    = NULL,
                              vecNormConst        = vecNormConst,
                              scaNProc            = NULL,
                              scaQThres           = NULL,
                              strReport           = NULL)
  }
  
  return(objectLineagePulse)
}