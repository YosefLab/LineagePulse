#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# DEVELOPER COMMENTS ON ALGORITHM BUILDING AND IMPLEMENTATION
# Note: Most problems are numerical and not conceptual problems.

# Previously encountered problems:
# Genes with many zeros and few very high counts: If dispersion is initialised
# too high (=1), the increase in LL cannot always be picked up by dnbinom which
# disables convergence.
# Still didnt see convergence, this might reoccur during fitting. Removed -Inf
# masking from likelihood function for NB part. Seeing convergence of means now.
# Still problems with drop-out rate -> Use more genes to get more reliable 
# traces, still problems. Take out masking of dropout contribution of nonzeros
# in likelihood and force logistic to be decreasing in means - good cnvergence
# until 3rd/4th iter.
# Try replace mean optimise with optim and put masking in likelihood again: Note
# that optim becomes unstable if Inf is returned by objective, while optimise
# seems to be forced out of the "bad" parameter region in which dnbinom yields
# LL=-Inf, optim needs the low value masking to not break down. Optim can be
# initialised outside of the bad parameter region to avoid it. Now only 
# problems in mean estimation (3rd iteration). Remove bounds from BFGS
# optimisation. Cleared problem by setting mean initialisation to log(max+1)
# instead of log(mean+1), was it stuck in bad parameter region again? maybe.
# Initialise mean and drop-out rate to prior value, converging fully!
# Compared BFGS against Nelder-Mead in optim for mean estimation on 200x400
# points: Completed 4 iter with BFGS in 36.03 min on 2 cores 
# (LL=24536057.4015036). Completed 4 iter with Nelder-Mead in 35.22 min on 
# 2 cores (LL=-24536058.6125272). Exact same data set. Both yield 172 DE genes.
# Problem with impulse mode: LL returned from optim does not match LL
# computed from parameters returned by optim. Effect not so heavy in first 
# iteration but very heavy after that. Not due to matLinModel, is it dispersion?
# Dispersion is initialised low, if initialised higher, effects are stronger
# in first iteration already. Is this inaccuracy in reporting parameters? Wouldnt
# be so bad in first iteration with low dispersion param!?
# Arrived at a converging version with ok results -> means are underestimated,
# variance overestimated and drop-out rates are almost all zero if estimating
# on many genes. This could be a bad local optimum of the likelihood, in which
# the non-drop-out componen provides a lot of probability density for the zeros
# and thereby pushing the drop-out component out. This seems to be show in the 
# loglikelihood convergence as a large jump (10%) after the first dispersion
# estimation. Avoid: initialise mu high,
# variance low (dispersion high) and drop-out rate high. This works on a small
# data set. Previously: mu= mean of nonzeros, phi = 0.001, dropout model:
# offset = 0, mu param = -1. Now:  mu= mean of nonzeros, phi = 1, dropout model:
# offset = log(scaPiTarget) - log(1-scaPiTarget) - scaPiLinModelMuParam*
# log(min(vecMuModelInit, na.rm=TRUE)), mu param = -1.

# HOWTO debug this: 
# Set bplapply of estimation that throws error to lapply for proper error
# reporting. Or - give nProc=1 and BiocParallel operates in SerialMode with
# proper error reporting -> all estiamtions are sequential though.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Coordinate mean parameter estimation step
#' 
#' Auxillary function that calls the estimation functions for the
#' different mean models according to their needs. This function
#' only saves space in fitZINB.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecPseudotime: (numerical vector number of cells)
#'    [Default NULL]
#'    Pseudotime coordinates of cells. Only required if mean model
#'    or dispersion model are fit as a function of pseudotime, 
#'    e.g. impulse model for means.
#' @param vecindClusterAssign: (integer vector length number of
#'    cells) [Default NULL] Index of cluster assigned to each cell.
#' @param matMu: (numeric matrix genes x cells) [Default NULL]
#'    Inferred zero inflated negative binomial mean parameters on 
#'    which estimation in this step  is conditioned on. 
#'    Not required for all models.
#' @param  matDispersions: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial dispersion parameters
#'    on which estimation in this step is conditioned on.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates
#'    on which estimation in this step is conditioned on.
#'    These are the observation-wise point estimates, not the
#'    logistic functions.
#' @param matPiConstPredictors: (numeric matrix genes x external 
#'    predictors) [Default NULL]
#' @param matDropoutLinModel: (numeric matrix genes x predictors)
#'    [Default NULL]
#' @param matImpulseParam: (numeric matrix genes x impulse 
#'    parameters 6) Inferred impulse model parameters
#'    if strMuModel is "impulse". NA for all other strMuModel.
#' @param matboolZero: (numeric matrix genes x cells)
#'    [Default NULL] Whether observation is zero.
#' @param matboolNotZeroObserved: (numeric matrix genes x cells)
#'    [Default NULL] Whether observation is non-zero and real.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param boolVecWindowsAsBFGS: (bool) [Default FALSE] Whether
#'    mean parameters of a gene are co-estimated in "windows"
#'    mode with BFGS algorithm (optimisation with dimensionality
#'    of number of cells) or estimated one by one, conditioned
#'    one the latest estimates of neighbours. The latter case
#'    (boolVecWindowsAsBFGS=FALSE) is coordinate ascent within the gene
#'    and each mean parameter is optimised once only.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param MAXIT_BFGS_Impulse: (scalar) [Default 1000]
#'    Maximum number of BFGS iterations used to estiamte
#'    an impulse model for a gene.
#' @param RELTOL_BFGS_Impulse: (scalar) [Default 
#'    sqrt(.Machine$double.eps)] Minimum relativ decrease
#'    in the loglikelihood objective within one iteration
#'    of BFGS impulse model fitting to gene for the estimation
#'    to not complete as converged. The default is the 
#'    default of optim which is overridden by the value set in
#'    fitZINB by default in LineagePulse to save computation 
#'    time.
#' 
#' @return matMuModel: (numeric matrix genes x mu model parameters)
#'    Contains the model parameters according to the used model.
#' @export

fitZINBMu <- function( matCountsProc,
  vecSizeFactors,
  lsMuModel,
  lsDispModel,
  lsDropModel,
  scaWindowRadius ){

  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  if(lsMuModel$lsMuModelGlobal$strMuModel=="windows"){
    # Estimate mean parameter for each cell as ZINB model for cells within pseudotime
    # interval with cell density centred at target cell.
    # Note that this corresponds to maximising smoothed log likelihood but instead
    # of using the implemented cost function evalLogLikSmoothZINB_LinPulse for an entire
    # gene, the optimisation problem is broken up into 1D problems for each mean.
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      # Note: Mean parameter estimates of gene i depend on each other:
      # Either estimate all parameters for gene i together 
      # (quasi-Newton estimation with BFGS: boolVecWindowsAsBFGS=TRUE)
      # or use the latest updates of the remaining parameters during 
      # one-by-one estimation (coordinate ascent, boolVecWindowsAsBFGS=TRUE).
      
      # Decompress parameters
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecInterval=NULL )
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      if(lsMuModel$lsMuModelGlobal$boolVecWindowsAsBFGS){
        vecMu <- fitMuVecWindowsZINB_LinPulse(
          vecCounts=matCountsProc[i,],
          vecMu=vecMuParam,
          vecDisp=vecDispParam,
          vecNormConst=vecSizeFactors,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius )
      } else {
        vecMu <- vecMuParam
        for(j in seq(1,scaNumCells)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          vecMu[j] <- fitMuWindowZINB_LinPulse(vecCounts=matCountsProc[i,vecInterval],
            vecMu=vecMu[vecInterval],
            vecDisp=vecDispParam[vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
            vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
            scaTarget=match(j,vecInterval),
            scaWindowRadius=scaWindowRadius )
        }
      }
      return(vecMu)
    }))
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="clusters"){
    # Estimate mean parameter by cluster. No smoothing is used.
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate mean parameters
      vecMu <- sapply(seq(1,max(lsMuModel$lsMuModelGlobal$vecindClusterAssign)), function(k){
        vecInterval <- lsMuModel$lsMuModelGlobal$vecindClusterAssign==k
        scaMu <- fitMuClusterZINB_LinPulse(
          vecCounts=matCountsProc[i,vecInterval],
          scaMuGuess=lsMuModel$matMuModel[i,k],
          vecDisp=vecDispParam[vecInterval],
          vecNormConst=vecSizeFactors[vecInterval],
          matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
        return(scaMu)
      })
      return(vecMu)
    }))

  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){
    
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecInterval=NULL )
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      vecDropoutParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
        vecMu=vecMuParam,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
      
      #  Compute posterior for parameter
      vecZ <- calcPostDrop_Vector( vecMu=vecMuParam,
        vecDispersions=vecDispParam,
        vecDropout=vecDropoutParam,
        vecboolZero= matCountsProc[i,]==0,
        vecboolNotZeroObserved= !is.na(matCountsProc[i,]) & matCountsProc[i,]>0,
        scaWindowRadius=scaWindowRadius )
      
      # Estimate mean parameters
      lsImpulseFit <- fitMuImpulseZINB_LinPulse(
        vecCounts=matCountsProc[i,],
        vecDisp=vecDispParam,
        vecNormConst=vecSizeFactors,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        vecProbNB=1-vecZ,
        vecPseudotime=lsMuModel$lsMuModelGlobal$vecPseudotime,
        vecImpulseParam=lsMuModel$matMuModel[i,],
        scaWindowRadius=scaWindowRadius,
        MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_Impulse,
        RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_Impulse )
      return(lsImpulseFit$vecBestFitParam)
    }))
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="constant"){
    matMuModel <- do.call(rbind, bplapply( seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate constant mean parameter
      scaMu <- fitMuConstZINB( vecCounts=matCountsProc[i,],
        scaMuGuess=lsMuModel$matMuModel[i,],
        vecDisp=vecDispParam,
        vecNormConst=vecSizeFactors,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        scaWindowRadius=scaWindowRadius )
      return(scaMu)
    }))
  }
  
  return( matMuModel=matMuModel )
}

#' Fit zero-inflated negative binomial model to data
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
#' Convergence is tracked with the theloglikelihood of the entire 
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
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
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
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' @param boolSuperVerbose: (bool) Whether to follow EM-algorithm
#'    progress in high detail with local convergence flags.
#' 
#' @return (list)
#'    \itemize{
#'      \item matMuH1: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial mean parameters
#'        under alternative model H1.
#'      \item matDispersionsH1: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial dispersion parameters
#'        under alternative model H1.
#'      \item matDropoutH1: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial drop out rates.
#'        These are the observation-wise point estimates, not the
#'        logistic functions. Fit under alternative model H1.
#'      \item matMuH0: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial mean parameters
#'        under null model H0.
#'      \item matDispersionsH0: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial dispersion parameters
#'        under null model H0.
#'      \item matDropoutH0: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial drop out rates.
#'        These are the observation-wise point estimates, not the
#'        logistic functions. Fit under null model H0.
#'      \item matDropoutLinModel: (numeric matrix cells x logistic parameters)
#'        Logistic dropout rate model for each cell. Inferred based on
#'        alternative model H1 or null model H0 dependent on 
#'        \code{boolEstimateNoiseBasedOnH0).
#'      \item matImpulseParam: (numeric matrix genes x impulse 
#'        parameters 6) Inferred impulse model parameters
#'        if strMuModel is "impulse". NA for all other strMuModel.
#'      \item lsFitZINBReporters: (list length 6)
#'        \itemize{
#'           \item boolConvergenceH1: (bool) Convergence of
#'            estimation for alternative model H1. Convergence
#'            is evaluated based on the convergence of the 
#'            loglikelihood of the entire data set.
#'          \item boolConvergenceH0: (bool) Convergence of
#'            estimation for null model H0. Convergence
#'            is evaluated based on the convergence of the 
#'            loglikelihood of the entire data set.
#'          \item vecEMLogLikH1: (numeric vector number of 
#'            estimation cycles) Loglikelihood of entire 
#'            data set after each estimation cycle of alternative
#'            model H1.
#'          \item vecEMLogLikH0: (numeric vector number of 
#'            estimation cycles) Loglikelihood of entire 
#'            data set after each estimation cycle of null
#'            model H0.
#'          \item scaKbyGeneH1: (scalar) Degrees of freedom
#'            by gene used in alternative model H1. The logistic
#'            dropout model is ignored as it is shared between
#'            alternative and null model.
#'          \item scaKbyGeneH0: (scalar) Degrees of freedom
#'            by gene used in null model H0. The logistic
#'            dropout model is ignored as it is shared between
#'            alternative and null model.
#'        }
#'    }
#' @export

fitZINB <- function(matCountsProc, 
  lsResultsClustering,
  vecSizeFactors,
  scaWindowRadius=NULL,
  strMuModel="windows",
  strDispModel = "constant",
  boolEstimateNoiseBasedOnH0=TRUE,
  boolVecWindowsAsBFGS=FALSE,
  vecPseudotime=NULL,
  scaMaxEstimationCycles=20,
  nProc=2,
  dirBPLogs,
  verbose=FALSE,
  boolSuperVerbose=FALSE ){
  
  ####################################################
  # Internal Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-6)
  MAXIT_BFGS_Impulse <- 1000 # optim default is 1000
  RELTOL_BFGS_Impulse <- 10^(-4) # optim default is sqrt(.Machine$double.eps)=1e-8
  # Lowering RELTOL_BFGS_IMPULSE gives drastic run time improvements.
  # Set to 10^(-4) to maintain sensible fits without running far into saturation
  # in the objective (loglikelihood).
  
  # Store EM convergence
  vecEMLogLikModelA <- array(NA,scaMaxEstimationCycles)
  vecEMLogLikModelB <- array(NA,scaMaxEstimationCycles)
  scaPredictors <- 2
  matPiConstPredictors <- NULL
  
  # Set smoothing mode
  boolSmoothed <- FALSE
  if(!is.null(scaWindowRadius)){
    if(scaWindowRadius>0){
      boolSmoothed <- TRUE
    } else {
      scaWindowRadius <- NULL
    }
  }
  
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
  # Initialise function
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]  
  
  # Pre-compute constant matrices
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
  matboolZero <- matCountsProc == 0
  
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
  # Initialise model A
  # Initialisation repeated for second estimation to save memory
  print("Initialise parameters model A")
  # Initialise parameters:
  # Mean parameters (mu): Gene-wise mean of non-zero observations.
  # Impulse model: Initialised to constant (mean).
  # Dispersions: Low dispersion factor yielding high variance which makes
  # cost function screening easy in the first iteration.
  # Dropout model: Initialise as offset=0 and log(mu)  parameter which
  # is forced to be negative during fitting, as -1. The parameter corresponding
  # to log(mu) may not be initialised too close to zero, as the cost function 
  # cannot always pick up the signal in such cases, leading to an MLE with this 
  # parameter untouched.
  lsMuModelA <- list( matMuModel=NA,
    lsMuModelGlobal=list( strMuModel=strMuModelA,
      scaNumCells=scaNumCells,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign,
      boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
      MAXIT_BFGS_Impulse=MAXIT_BFGS_Impulse,
      RELTOL_BFGS_Impulse=RELTOL_BFGS_Impulse) )
  vecMuModelInit <- apply(matCountsProc, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
  vecMuModelInit[vecMuModelInit < .Machine$double.eps] <- .Machine$double.eps
  if(strMuModelA=="constant"){
    lsMuModelA$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  } else if(strMuModelA=="impulse"){
    lsMuModelA$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=6)
    lsMuModelA$matMuModel[,2:4] <- log(matrix(vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE))
  } else if(strMuModelA=="clusters"){
    lsMuModelA$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
  } else  if(strMuModelA=="windows"){
    lsMuModelA$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  } else {
    stop(paste0("ERROR fitZINB(): strMuModelA=", strMuModelA, " not recognised."))
  }
  
  lsDispModelA <- list( matDispModel=NA,
    lsDispModelGlobal=list( strDispModel=strDispModelA,
      scaNumCells=scaNumCells,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign) )
  if(strDispModelA=="constant"){
    lsDispModelA$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  } else {
    stop(paste0("ERROR fitZINB(): strDispModelA=", strDispModelA, " not recognised."))
  }
  
  lsDropModel <- list(matDropoutLinModel=NA,
    matPiConstPredictors=matPiConstPredictors)
  # Target initialisation drop-out rate: 0.9, linear model mu
  # parameter = -1 -> solve for offset of linear model:
  scaPiTarget <- 0.99
  scaPiLinModelMuParam <- -1
  scaPiLinModelOffset <- log(scaPiTarget) - log(1-scaPiTarget) - 
    scaPiLinModelMuParam*log(min(vecMuModelInit, na.rm=TRUE))
  lsDropModel$matDropoutLinModel <- cbind(
    rep(scaPiLinModelOffset, scaNumCells), 
    rep(scaPiLinModelMuParam, scaNumCells),
    matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
  
  ####################################################
  # Fit model A
  print(paste0("### a) Fit negative binomial model A (",
    strNameModelA,") with noise model."))
  
  # Evaluate initialisation loglikelihood for model B
  scaLogLikInitA <- evalLogLikMatrix(matCounts=matCountsProc,
    vecSizeFactors=vecSizeFactors,
    lsMuModel=lsMuModelA,
    lsDispModel=lsDispModelA, 
    lsDropModel=lsDropModel,
    scaWindowRadius=scaWindowRadius )
  if(verbose){
    print(paste0("Completed initialisation with ",
      "log likelihood of         ", scaLogLikInitA))
  }
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- scaLogLikInitA
  scaLogLikOld <- NA
  
  tm_cycle <- system.time({
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
      tm_iter <- system.time({
        #####  1. Cell-wise parameter estimation
        # Dropout rate
        # Drop-out estimation is independent between cells and can be parallelised.
        tm_pi <- system.time({
          lsDropModel$matDropoutLinModel <- do.call(rbind, 
            bplapply(seq(1, scaNumCells), function(cell){
              if(!is.null(scaWindowRadius)){
                scaindIntervalStart <- max(1,cell-scaWindowRadius)
                scaindIntervalEnd <- min(scaNumCells,cell+scaWindowRadius)
                vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
              } else {
                vecInterval <- cell
              }
              
              vecDropoutLinModel <- fitPiZINB_LinPulse(
                vecCounts=matCountsProc[,cell],
                vecDropoutLinModel=lsDropModel$matDropoutLinModel[cell,],
                matPiConstPredictors=lsDropModel$matPiConstPredictors,
                lsMuModel=lsMuModelA,
                lsDispModel=lsDispModelA,
                scaNormConst=vecSizeFactors[cell],
                vecInterval=vecInterval,
                scaTarget=match(cell, vecInterval))
              return(vecDropoutLinModel)
            }))
        })
        if(boolSuperVerbose){
          scaLogLikTemp <- evalLogLikMatrix( matCounts=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelA,
            lsDispModel=lsDispModelA, 
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
          print(paste0("# ",scaIter,".1) Drop-out estimation complete: ",
            "loglikelihood of   ", scaLogLikTemp, " in ",
            round(tm_pi["elapsed"]/60,2)," min."))
        }  
        
        ##### 2. Gene-wise parameter estimation: 
        # a) Negative binomial mean parameter
        tm_mu <- system.time({
          lsMuModelA$matMuModel <- fitZINBMu( matCountsProc=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelA,
            lsDispModel=lsDispModelA,
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
        })
        colnames(lsMuModelA$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
        if(boolSuperVerbose){
          scaLogLikTemp <- evalLogLikMatrix( matCounts=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelA,
            lsDispModel=lsDispModelA, 
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
          print(paste0("# ",scaIter, ".2) Mean estimation complete: ",
            "loglikelihood of       ", scaLogLikTemp, " in ",
            round(tm_mu["elapsed"]/60,2)," min."))
        }
        
        # b) Negative binomial dispersion parameter
        # Use MLE of dispersion factor: numeric optimisation of likelihood.
        tm_phi <-system.time({
          if(lsDispModelA$lsDispModelGlobal$strDispModel=="constant"){  
            vecDispFitModelA <- bplapply(seq(1,scaNumGenes), function(i){
              # Decompress parameters
              vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModelA$matMuModel[i,],
                lsMuModelGlobal=lsMuModelA$lsMuModelGlobal,
                vecInterval=NULL )
              vecDropoutParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
                vecMu=vecMuParam,
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
              
              # Estimate constant dispersion factor
              fitDisp <- fitDispZINB_LinPulse(scaDispGuess=lsDispModelA$matDispModel[i,],
                vecCounts=matCountsProc[i,],
                vecSizeFactors=vecSizeFactors,
                vecMuEst=vecMuParam,
                vecDropoutRateEst=vecDropoutParam,
                vecboolNotZeroObserved= !is.na(matCountsProc[i,]) & matCountsProc[i,]>0, 
                vecboolZero= matCountsProc[i,]==0,
                scaWindowRadius=scaWindowRadius,
                boolSmoothed=boolSmoothed )
              return(fitDisp)
            })
          } else {
            #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
            print(paste0("Dispersion parameter model not recognised: ", strDispModel, 
              ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
            stop(paste0("Dispersion parameter model not recognised: ", strDispModel, 
              ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
          }
          vecboolConvergedGLMdispModelA <- sapply(vecDispFitModelA, function(fit) fit["convergence"])
          lsDispModelA$matDispModel <- do.call(rbind, lapply(vecDispFitModelA, function(fit){ fit["par"] }))
          colnames(lsDispModelA$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
        })
        
        # Evaluate Likelihood
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- evalLogLikMatrix( matCounts=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelA,
            lsDispModel=lsDispModelA, 
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
      })
      
      # Iteration complete
      if(boolSuperVerbose){
        if(any(vecboolConvergedGLMdispModelA != 0)){
          print(paste0("Dispersion estimation did not converge in ", 
            sum(vecboolConvergedGLMdispModelA), " cases."))
        }
        print(paste0("# ",scaIter, ".3) Dispersion estimation complete: ",
          "loglikelihood of ", scaLogLikNew, " in ",
          round(tm_phi["elapsed"]/60,2)," min."))
      } else {
        if(verbose){print(paste0("# ",scaIter, ".) complete with ",
          "log likelihood of ", scaLogLikNew, " in ",
          round(tm_iter["elapsed"]/60,2)," min."))}
      }
      vecEMLogLikModelA[scaIter] <- scaLogLikNew
      scaIter <- scaIter+1
    }
  })
  print(paste0("Finished fitting zero-inflated negative binomial ",
    "model A with noise model in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  # Evaluate convergence
  if(all(as.logical(vecboolConvergedGLMdispModelA)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergenceModelA <- TRUE
  } else { boolConvergenceModelA <- FALSE }
  
  ####################################################
  # Initialise model B
  
  print("Initialise parameters model B")
  # Initialise parameters:
  # -> Mean model is re-fit from scratch
  # -> Dispersion estimates from A are used for initialisation
  # -> Drop-out model is kept from model A estimation and
  #   point estimators adjusted to mu initialisation used.
  lsMuModelB <- list( matMuModel=NA,
    lsMuModelGlobal=list( strMuModel=strMuModelB,
      scaNumCells=scaNumCells,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign,
      boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
      MAXIT_BFGS_Impulse=MAXIT_BFGS_Impulse,
      RELTOL_BFGS_Impulse=RELTOL_BFGS_Impulse) )
  # B is alternative model H1
  if(strMuModelB=="constant"){
    lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  } else if(strMuModelB=="impulse"){
    lsMuModelB$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=6)
    lsMuModelB$matMuModel[,2:4] <- log(matrix(vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE))
  } else if(strMuModelB=="clusters"){
    lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
  } else  if(strMuModelB=="windows"){
    lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  } else {
    stop(paste0("ERROR fitZINB(): strMuModelB=", strMuModelB, " not recognised."))
  }
  
  lsDispModelB <- list( matDispModel=NA,
    lsDispModelGlobal=list( strDispModel=strDispModelB,
      scaNumCells=scaNumCells,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign) )
  if(strDispModelB=="constant"){
    if(strDispModelA=="constant"){
      # Use values estimated for model A as initialisation
       lsDispModelB$matDispModel <-  lsDispModelA$matDispModel
    } else {
      # Initialise from scratch if different model used
      lsDispModelB$matDispModel <- matrix(0.001, nrow=scaNumGenes, ncol=1, byrow=FALSE)
    }
  } else {
    stop(paste0("ERROR fitZINB(): strDispModelB=", strDispModelB, " not recognised."))
  }
  
  ####################################################
  # Fit model B
  print(paste0("### b) Fit negative binomial model B (",
    strNameModelB,")."))
  
  # Evaluate initialisation loglikelihood for model B
  scaLogLikInitB <- evalLogLikMatrix(matCounts=matCountsProc,
    vecSizeFactors=vecSizeFactors,
    lsMuModel=lsMuModelB,
    lsDispModel=lsDispModelB, 
    lsDropModel=lsDropModel,
    scaWindowRadius=scaWindowRadius )
  if(verbose){
    print(paste0("Completed initialisation with ",
      "log likelihood of         ", scaLogLikInitB))
  }
  
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- scaLogLikInitB
  scaLogLikOld <- NA
  
  tm_cycle <- system.time({
    while(scaIter == 1 |
        (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
      tm_iter <- system.time({
        ##### Gene-wise parameter estimation: 
        # a) Negative binomial mean parameter
        # Only compute posterior if using closed form estimator for mean:
        # Posterior is not necessary in all other cases, expect for parameter
        # initialisation of impulse model.
        tm_mu <- system.time({
          lsMuModelB$matMuModel <- fitZINBMu( matCountsProc=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelB,
            lsDispModel=lsDispModelB,
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
        })
        colnames(lsMuModelB$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
        if(boolSuperVerbose){
          scaLogLikTemp <- evalLogLikMatrix( matCounts=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelB,
            lsDispModel=lsDispModelB, 
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
          print(paste0("# ",scaIter, ".1) Mean estimation complete: ",
            "loglikelihood of       ", scaLogLikTemp, " in ",
            round(tm_mu["elapsed"]/60,2)," min."))
        }
        
        # b) Negative binomial dispersion parameter
        # Use MLE of dispersion factor: numeric optimisation of likelihood.
        tm_phi <- system.time({
          if(lsDispModelB$lsDispModelGlobal$strDispModel=="constant"){  
            vecDispFitModelB <- bplapply(seq(1,scaNumGenes), function(i){
              # Decompress parameters
              vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModelB$matMuModel[i,],
                lsMuModelGlobal=lsMuModelB$lsMuModelGlobal,
                vecInterval=NULL )
              vecDropoutParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
                vecMu=vecMuParam,
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
              
              # Estimate constant dispersion factor
              fitDisp <- fitDispZINB_LinPulse(scaDispGuess=lsDispModelB$matDispModel[i,],
                vecCounts=matCountsProc[i,],
                vecSizeFactors=vecSizeFactors,
                vecMuEst=vecMuParam,
                vecDropoutRateEst=vecDropoutParam,
                vecboolNotZeroObserved= !is.na(matCountsProc[i,]) & matCountsProc[i,]>0, 
                vecboolZero= matCountsProc[i,]==0,
                scaWindowRadius=scaWindowRadius,
                boolSmoothed=boolSmoothed )
              return(fitDisp)
            })
          } else {
            #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
            print(paste0("Dispersion parameter model not recognised: ", strDispModel, 
              ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
            stop(paste0("Dispersion parameter model not recognised: ", strDispModel, 
              ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
          }
          vecboolConvergedGLMdispModelB <- sapply(vecDispFitModelB, function(fit) fit["convergence"])
          lsDispModelB$matDispModel <- do.call(rbind, lapply(vecDispFitModelB, function(fit){ fit["par"] }))
          colnames(lsDispModelB$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
        })
        
        # Evaluate Likelihood
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- evalLogLikMatrix( matCounts=matCountsProc,
            vecSizeFactors=vecSizeFactors,
            lsMuModel=lsMuModelB,
            lsDispModel=lsDispModelB, 
            lsDropModel=lsDropModel,
            scaWindowRadius=scaWindowRadius )
      })
      
      # Iteration complete
      if(boolSuperVerbose){
        if(any(vecboolConvergedGLMdispModelB != 0)){
          print(paste0("Dispersion estimation did not converge in ", 
            sum(vecboolConvergedGLMdispModelB), " cases."))
        }
        print(paste0("# ",scaIter, ".2) Dispersion estimation complete: ",
          "loglikelihood of ", scaLogLikNew, " in ",
          round(tm_phi["elapsed"]/60,2)," min."))
      } else {
        if(verbose){print(paste0("# ",scaIter, ".) complete with ",
          "log likelihood of ", scaLogLikNew, " in ",
          round(tm_iter["elapsed"]/60,2)," min."))}
      }
      vecEMLogLikModelB[scaIter] <- scaLogLikNew
      scaIter <- scaIter+1
    }
  })
  print(paste0("Finished fitting zero-inflated negative ",
    "binomial model B in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  # Evaluate convergence
  if(all(as.logical(vecboolConvergedGLMdispModelB)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergenceModelB <- TRUE
  } else { boolConvergenceModelB <- FALSE }
  
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
    lsReturn <- list( lsMuModelH1=lsMuModelB,
      lsDispModelH1=lsDispModelB,
      lsMuModelH0=lsMuModelA,
      lsDispModelH0=lsDispModelA,
      lsDropModel=lsDropModel,
      lsFitZINBReporters=lsFitZINBReporters )
  } else {
    lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelA,
      boolConvergenceH0=boolConvergenceModelB,
      vecEMLogLikH1=vecEMLogLikModelA,
      vecEMLogLikH0=vecEMLogLikModelB,
      scaKbyGeneH1=scaKbyGeneH1,
      scaKbyGeneH0=scaKbyGeneH0 )
    lsReturn <- list( lsMuModelH1=lsMuModelA,
      lsDispModelH1=lsDispModelA,
      lsMuModelH0=lsMuModelB,
      lsDispModelH0=lsDispModelB,
      lsDropModel=lsDropModel,
      lsFitZINBReporters=lsFitZINBReporters )
  } 
  
  return(lsReturn)
}