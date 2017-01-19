################################################################################
#######################     LineagePulse package     ###########################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

library(BiocParallel)
#library(BatchJobs)
library(compiler)
library(ggplot2)
library(MASS)
#setwd("/Users/davidsebastianfischer/gitDevelopment/LineagePulse/R")
#setwd("/data/yosef2/users/fischerd/code/LineagePulse/R")
setwd("/home/david/gitDevelopment/code/LineagePulse/R")

source("srcLineagePulse_clusterCellsInPseudotime.R")
source("srcLineagePulse_calcPostDrop.R")
source("srcLineagePulse_calcNormConst.R")
source("srcLineagePulse_classLineagePulseObject.R")
source("srcLineagePulse_decompressParameters.R")
source("srcLineagePulse_estimateMMAssignments.R")
source("srcLineagePulse_evalDropoutModel.R")
source("srcLineagePulse_evalImpulseModel.R")
source("srcLineagePulse_evalLogLikZINB.R")
source("srcLineagePulse_fitZINB_cofitMeanDispersion.R")
source("srcLineagePulse_fitZINB_fitMean.R")
source("srcLineagePulse_fitZINB_fitDispersion.R")
source("srcLineagePulse_fitZINB_fitDropout.R")
source("srcLineagePulse_fitZINB.R")
source("srcLineagePulse_fitZINB_WrapperMixture.R")
source("srcLineagePulse_fitZINB_WrapperLP.R")
source("srcLineagePulse_initialiseImpulseParameters.R")
source("srcLineagePulse_plotComparativeECDF.R")
source("srcLineagePulse_plotGene.R")
source("srcLineagePulse_plotPseudotimeClustering.R")
source("srcLineagePulse_processSCData.R")
source("srcLineagePulse_processSCDataMixture.R")
source("srcLineagePulse_runDEAnalysis.R")
source("srcLineagePulse_simulateDataSet.R")
source("srcLineagePulse_sortGeneTrajectories.R")
source("srcLineagePulse_validateOutput.R")
source("srcLineagePulse_validateOutputSimulation.R")

evalLogLikMuWindowZINB_comp <- cmpfun(evalLogLikMuWindowZINB)
evalLogLikMuVecWindowsZINB_comp <- cmpfun(evalLogLikMuVecWindowsZINB)
evalLogLikMuConstZINB_comp <- cmpfun(evalLogLikMuConstZINB)
evalLogLikMuImpulseZINB_comp <- cmpfun(evalLogLikMuImpulseZINB)
evalLogLikDispConstZINB_comp <- cmpfun(evalLogLikDispConstZINB)

################################################################################
### Main function
################################################################################

#' LineagePulse wrapper: Differential expression analysis in pseudotime
#' 
#' This function performs all steps of longitudinal differential
#' expression analysis in pseudotime for you.
#' 
#' This function is the wrapper function for the LineagePulse algorithm
#' which performs zero-inflated negative binomial model fitting
#' and differential expression analysis as well as
#' output validation. Note that LineagePulse has many input parameters but
#' only few will be relevant for you and you will be able to leave the 
#' remaining ones as their defaults. Read up on specific input parameters
#' in the input parameter annotation of this function. In short,
#' you have to:
#' 
#' 1. Supply data: Count data (matCounts) and pseudotime coordinates
#' of the cells (vecPseudotime). You may decide to also provide cell-wise
#' normalisation factors (such as factors accounting for sequencing depth)
#' (vecNormConstExternal). To ease LineagePulse use when testing
#' multiple lineages within one data set, the input that controls the 
#' set of cells entering LineagePulse is vecPseudotime. Cells not
#' mentioned in vecPseudotime are not included.
#' 2. Decide whether you want to run a test run on a few genes only 
#' (scaSmallRun).
#' 3. Chose the model constraining mean (strMuModel) and dispersion 
#' parameters (strDispModel) for each gene. If you run clusters,
#' you may decide to force the number of clusters (scaKCluster)
#' rather than using internal model selection or use clusters
#' based on true time (time of sampling) (boolClusterInPseudotime).
#' 4. Decide whether you want to use local negative binomial model
#' smooting (scaWindowRadius). 
#' 5. Supply gene-specific drop-out predictors if wanted 
#' (matPiConstPredictors).
#' 6. Set optimisation parameters (boolEstimateNoiseBasedOnH0,
#' boolVecWindowsAsBFGS, boolCoEstDispMean, scaMaxEstimationCycles).
#' 7. Chose the number of processes you want to use (scaNProc), LineagePulse
#' is parallelised on all computation intensive steps. Note that
#' the current parallelisation scheme runs on Unix (MacOS) and Linux but
#' not on windows. Adjust the code section within this wrapper to
#' parallelise on windows.
#' 8. Set the validation output you want to receive to visualise
#' the results (boolValidateZINBfit).
#' 9. Set the level of detail with which you want to follow
#' progress through text printed on the console during a run
#' (verbose, boolSuperVerbose).
#' 
#' Finally, after running LineagePulse, you may continue to work
#' on your data set by:
#' A) Regenerating observation-wise parameters,
#' such as the mean parameter matrix which represents the hidden
#' expression states, with the functions in srcLineagePulse_decompressParameters.R.
#' B) You can also compute the observation-wise probability of 
#' dropout with calcPostDrop.
#' C) Moreover, you can have a closer look at the hidden expression
#' states of the genes with sortGeneTrajectories.
#' D) If you work on simulated data, you can create additional 
#' validation metrics with validateOutputSimulation.
#' 
#' @details Note on optimisation strategies: LineagePulse gives the use
#' control over how the large optimisation problem is broken up.
#' We give the user this choice so that the user can adjust the 
#' accuracy to the available computing infrastructure. Consult
#' with the LineagePulse developers if you are unsure about he settings
#' here, e.g. if you are not sure whether to estimate the noise
#' model on the H1 or H0 expression model.
#' 
#' Note on parameter objects:
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
#' The observation-wise parameter estimates can be recovered with 
#' the functions in: srcLineagePulse_decompressParameters.R
#'
#' The computational complexity of LineagePulse is linear in the
#' number of genes and linear in the number of cells.
#' 
#' @aliases LineagePulse wrapper, main function
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied
#' @param vecNormConstExternal: (numeric vector number of cells) 
#'    Model scaling factors, one per cell. These factors will linearly 
#'    scale the mean model for evaluation of the loglikelihood. 
#'    Must be named according to the column names of matCounts.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param boolClusterInPseudotime: (bool) [Default TRUE]
#'    Whether to cluster cells in pseudotime. If FALSE,
#'    time points supplied in vecPseudotime are treated as clusters.
#'    This requires that time points in vecPseudotime occur 
#'    multiple times and are for example the real time of sampling
#'    of a cell (e.g. how many hours into the experiment).
#' @param scaKCluster: (integer) [Default NULL] Forces number of centroids
#'    in K-means to be K: setting this to an integer (not NULL) skips model
#'    selection in clusterting.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
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
#' @param boolCoEstDispMean: (bool) [Default TRUE]
#'    Whether mean and dispersion parameters are to be co-estimated
#'    (simulatneous optimisation). Only available for certain 
#'    dispersion and mean models:
#'    dispersion models: constant.
#'    mean models: constant, cluster, sliding window vector, impulse.
#'    Note that co-estimation in model estimation B (without drop-
#'    out model estimation) leads to a single step estimation as 
#'    mean and dispersion parameter don't have to be iterated over.
#'    This makes estimation of large data sets with complex H1 mean
#'    model (e.g. impulse) possible, as the drop-out model can be 
#'    estimated based on H0 (boolEstimateNoiseBasedOnH0) so that
#'    the complex model only has to be estimated once (simultaneous
#'    with the dispersion parameters). This may generally lead to better
#'    convergence as the steps in coordinate-ascent are in a larger
#'    space, closer to full gradient ascent. Setting to TRUE
#'    is encouraged. Optimisation routines for individual mean 
#'    and dispersion fitting (if FALSE) exist, but these may be viewed
#'    as non-deprecated parts of an earlier implementation of the
#'    alorithm.
#' @param scaMaxEstimationCycles: (integer) [Default 20] Maximum number 
#'    of estimation cycles performed in fitZINB(). One cycle
#'    contain one estimation of of each parameter of the 
#'    zero-inflated negative binomial model as coordinate ascent.
#' @param scaNProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param boolValidateZINBfit: (bool) [Default TRUE]
#'    Whether to generate evaluation metrics and plots
#'    for parameter values of inferred ZINB model.
#' @param verbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation in high detail with local 
#'    convergence flags and step-by-step loglikelihood computation.
#'    Note: This increases run-time as the loglikelihood is computed
#'    more often. This is usually not a major contributor to runtime
#'    though.
#' 
#' @return dfDEAnalysis: (data frame genes x reported variables) 
#'    Summary of differential expression analysis, sorted by adj.p:
#'    {Gene: gene ID,
#'    p: raw p-value, 
#'    adj.p: BH corrected p-value, 
#'    loglik_full: loglikelihood of alternative model H1,
#'    loglik_red: loglikelihood of null model H0,
#'    deviance: loglikelihood ratio test statistic (the deviance),
#'    mean_H0: inferred gene-wise mean parameter (constant null model),
#'    dispersion_H0: inferred gene-wise dispersion parameter (constant null model)}
#'    
#' @author David Sebastian Fischer
#' 
#' @export
runLineagePulse <- function(matCounts,
  matPiConstPredictors=NULL,
  vecNormConstExternal=NULL,
  vecPseudotime,
  strMuModel="impulse",
  strDispModel="constant",
  boolClusterInPseudotime=TRUE,
  scaKClusters=NULL,
  scaWindowRadius=NULL,
  boolEstimateNoiseBasedOnH0=FALSE,
  boolVecWindowsAsBFGS=FALSE,
  boolCoEstDispMean=TRUE,
  scaMaxEstimationCycles=20,
  scaNProc=1,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE ){
  
  print("LineagePulse v1.0")
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  vecAllGenes <- rownames(matCounts)
  lsProcessedSCData <- processSCData( matCounts=matCounts,
    matPiConstPredictors=matPiConstPredictors,
    vecPseudotime=vecPseudotime,
    vecNormConstExternal=vecNormConstExternal,
    scaSmallRun=scaSmallRun,
    strMuModel=strMuModel,
    strDispModel=strDispModel,
    scaWindowRadius=scaWindowRadius,
    boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
    boolCoEstDispMean=boolCoEstDispMean )
  objectLineagePulse <- lsProcessedSCData$objectLineagePulse
  vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
  matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
  vecPseudotimeProc <- lsProcessedSCData$vecPseudotimeProc
  
  # Clear memory
  rm(matCounts)
  rm(matPiConstPredictors)
  rm(vecPseudotime)
  rm(lsProcessedSCData)
  
  # X. Inialise parallelisation
  print(paste0("Register parallelisation parameters: ", scaNProc, " threads."))
  # Set the parallelisation environment in BiocParallel:
  if(scaNProc > 1){
    # Set worker time out to 60*60*24*7 (7 days)
    # For single machine (FORK) cluster
    register(MulticoreParam(workers=scaNProc)) 
    #timeout=60*60*24*7,
    #log=FALSE, 
    #threshold="INFO", 
    #logdir=dirBPLogs))
    # Use this on windows or if SOCK clusters wanted:
    # For multiple machine (SOCK) cluster
    #register(SnowParam(workers=scaNProc, timeout=60*60*24*7))
  } else {
    # For debugging in serial mode
    register(SerialParam())
  }
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    if(boolClusterInPseudotime){
      # Cluster in pseudotime
      lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotimeProc,
        scaKexternal=scaKClusters)
    } else {
      # Take observation time points as clusters
      print("Chose given grouping (time points, cell types, conditions...).")
      lsResultsClustering <- list()
      lsResultsClustering[[1]] <- match(vecPseudotimeProc, sort(unique(vecPseudotimeProc)))
      lsResultsClustering[[2]] <- sort( unique(vecPseudotimeProc) )
      lsResultsClustering[[3]] <- length(unique(vecPseudotimeProc))
      names(lsResultsClustering) <- c("Assignments","Centroids","K")
    }
    # Plot clustering
    plotPseudotimeClustering(vecPseudotime=vecPseudotimeProc, 
      lsResultsClustering=lsResultsClustering)
  })
  print(paste("Time elapsed during clustering: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  
  # 4. Compute normalisation constants
  print("4. Compute normalisation constants:")
  objectLineagePulse <- calcNormConst(objectLineagePulse=objectLineagePulse,
     vecNormConstExternal=vecNormConstExternalProc)
  
  # 5. Fit ZINB model for both H1 and H0.
  print("5. Fit ZINB model for both H1 and H0.")
  tm_fitmm <- system.time({
    objectLineagePulse <- fitNullAlternative( objectLineagePulse=objectLineagePulse,
      matPiConstPredictors=matPiConstPredictorsProc,
      lsResultsClustering=lsResultsClustering,
      boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
      boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
      boolCoEstDispMean=boolCoEstDispMean,
      strMuModel=strMuModel,
      strDispModel=strDispModel,
      scaMaxEstimationCycles=scaMaxEstimationCycles,
      boolVerbose=boolVerbose,
      boolSuperVerbose=boolSuperVerbose )
  })
  print(paste("Time elapsed during ZINB fitting: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  
  # 6. Differential expression analysis:
  print("6. Differential expression analysis:")
  tm_deanalysis_mf <- system.time({
    objectLineagePulse <- runDEAnalysis( objectLineagePulse=objectLineagePulse )
  })
  print(paste("Time elapsed during differential expression analysis: ",
    round(tm_deanalysis_mf["elapsed"]/60,2)," min",sep=""))
  
  print("LineagePulse complete.")
  return(objectLineagePulse)
}