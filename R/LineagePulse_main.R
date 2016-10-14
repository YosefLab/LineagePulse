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
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/R")
#setwd("/data/yosef2/users/fischerd/code/LineagePulse/R")

source("srcLineagePulse_clusterCellsInPseudotime.R")
source("srcLineagePulse_createAnnotation.R")
source("srcLineagePulse_calcPostDrop.R")
source("srcLineagePulse_calcSizeFactors.R")
source("srcLineagePulse_decompressParameters.R")
source("srcLineagePulse_evalImpulseModel.R")
source("srcLineagePulse_evalLogLikZINB.R")
source("srcLineagePulse_fitZINB.R")
source("srcLineagePulse_fitZINB_fitMean.R")
source("srcLineagePulse_fitZINB_fitDispersion.R")
source("srcLineagePulse_fitZINB_cofitMeanDispersion.R")
source("srcLineagePulse_fitZINB_fitDropout.R")
source("srcLineagePulse_initialiseImpulseParameters.R")
source("srcLineagePulse_plotGene.R")
source("srcLineagePulse_plotPseudotimeClustering.R")
source("srcLineagePulse_plotZINBfits.R")
source("srcLineagePulse_processSCData.R")
source("srcLineagePulse_runDEAnalysis.R")
source("srcLineagePulse_simulateDataSet.R")
source("srcLineagePulse_validateOutput.R")
source("srcLineagePulse_validateOutputSimulation.R")
# Pre-compile function
evalImpulseModel_comp <- cmpfun(evalImpulseModel)
evalLogLikZINB_LinPulse_comp <- cmpfun(evalLogLikZINB_LinPulse)
evalLogLikSmoothZINB_LinPulse_comp <- cmpfun(evalLogLikSmoothZINB_LinPulse)
evalLogLikMuWindowZINB_LinPulse_comp <- cmpfun(evalLogLikMuWindowZINB_LinPulse)
evalLogLikMuVecWindowsZINB_LinPulse_comp <- cmpfun(evalLogLikMuVecWindowsZINB_LinPulse)
evalLogLikMuConstZINB_LinPulse_comp <- cmpfun(evalLogLikMuConstZINB_LinPulse)
evalLogLikMuImpulseZINB_LinPulse_comp <- cmpfun(evalLogLikMuImpulseZINB_LinPulse)
evalLogLikDispConstZINB_LinPulse_comp <- cmpfun(evalLogLikDispConstZINB_LinPulse)
evalLogLikPiZINB_LinPulse_comp <- cmpfun(evalLogLikPiZINB_LinPulse)
evalLogLikDispConstMuConstZINB_LinPulse_comp <- cmpfun(evalLogLikDispConstMuConstZINB_LinPulse)
evalLogLikDispConstMuImpulseZINB_LinPulse_comp <- cmpfun(evalLogLikDispConstMuImpulseZINB_LinPulse)

################################################################################
### Main function
################################################################################

#' LineagePulse wrapper: Differential expression analysis in pseudotime
#' 
#' This function is the wrapper function for the LineagePulse algorithm,
#' which performs data processing, clustering, zero-inflated negative
#' binomial model estimation), differential expression analysis and
#' output validation. Read up on the options you have in LineagePulse
#' in the input parameter annotation of this function. In short,
#' you have to:
#' 1. Supply data: Count data (matCounts) and pseudotime coordinates
#' of the cells (vecPseudotime).
#' 2. Decide whether you want to run a test run on a few genes only 
#' (scaSmallRun).
#' 3. Chose the model constraining mean (strMuModel) and dispersion 
#' parameters (strDispModel) for each gene. If you run clusters,
#' you may decide to force the number of clusters (scaKCluster)
#' rather than using internal model selection or use clusters
#' based on true time (time of sampling) (boolClusterInPseudotime).
#' 4. Decide whether you want to use local negative binomial model
#' smooting (scaWindowRadius). 
#' 5. Supply gene-specific drop-out predictors if wanted.
#' 6. Set optimisation parameters (boolEstimateNoiseBasedOnH0,
#' boolVecWindowsAsBFGS, boolCoEstDispMean, scaMaxEstimationCycles)
#' 7. Chose the number of processes you want to use (scaNProc), LineagePulse
#' is parallelised on all computation intensive steps.
#' 8. Set the validation output you want to receive to visualise
#' the results (boolPlotZINBfits, boolValidateZINBfit).
#' 9. Set the level of detail with which you want to follow
#' progress through text printed on the console during a run
#' (verbose, boolSuperVerbose).
#' 10. Chose the directory into which you want to have 
#' temporary and final output objects to be saved into
#' (dirOut).
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
#' 
#' @details The computational complexity of LineagePulse is linear in the
#' number of genes and linear in the number of cells.
#' 
#' @aliases LineagePulse wrapper
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param boolClusterInPseudotime: (bool) [Default TRUE]
#'    Whether to cluster cells in pseudotime. If FALSE,
#'    time points supplied in vecPseudotime are treated as clusters.
#'    This requires that time points in vecPseudotime occur 
#'    multiple times and are for example the real time of sampling
#'    of a cell (e.g. how many hours into the experiment).
#' @param scaKCluster: (integer) [Default NULL] Forces number of centroids
#'    in K-means to be K: setting this to an integer (not NULL) skips model
#'    selection in clusterting.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
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
#'    dispersion models: constant
#'    mean models: constant, impulse
#'    Note that co-estimation in model estimation B (without drop-
#'    out model estimation) leads to a single step estimation as 
#'    mean and dispersion parameter don't have to be iterated over.
#'    This makes estimation of large data sets with complex H1 mean
#'    model (e.g. impulse) possible, as the drop-out model can be 
#'    estimated based on H0 (boolEstimateNoiseBasedOnH0) so that
#'    the complex model only has to be estimated once (simultaneous
#'    with the dispersion parameters). This may generally lead to better
#'    convergence as the steps in coordinate-ascent are in a larger
#'    space, closer to full gradient ascent.
#' @param scaMaxEstimationCycles: (integer) [Default 20] Maximium number 
#'    of estimation cycles performed in fitZINB(). One cycle
#'    contain one estimation of of each parameter of the 
#'    zero-inflated negative binomial model as coordinate ascent.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param boolPlotZINBfits: (bool) [Default TRUE]
#'    Whether to plot zero-inflated negative binomial fits to selected genes
#'    and clusters.
#' @param boolValidateZINBfit: (bool) [Default TRUE]
#'    Whether to generate evaluation metrics and plots
#'    for parameter values of inferred ZINB model.
#' @param scaNProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param verbose: (bool) [Defaul TRUE]
#'    Whether progress of coordinate ascent within fitZINB is 
#'    reported once per iteration.
#' @param boolSuperVerbose: (bool) [Defaul FALSE]
#'    Whether coordinate ascent within fitZINB is followed
#'    step-by-step rather than reporting once per iteration.
#'    Note: This increases run-time as the loglikelihood is computed
#'    more often. This is usually not a major contributor to runtime
#'    though.
#' @param boolBPlog: (bool) [Default FALSE] Save status of each 
#'    worker into text files in parallelisation. Have to include
#'    library(BatchJobs) above and set the BiocParallel registration
#'    to adjust the reporting. Use only for debugging.
#' @param dirOut: (str directory) [Default NULL]
#'    Directory to which detailed output is saved to.
#'    Defaults to current working directory if NULL.
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
  vecPseudotime,
  scaSmallRun=NULL,
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
  boolPlotZINBfits = FALSE,
  boolValidateZINBfit=TRUE,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE,
  boolBPlog=FALSE,
  dirOut=NULL ){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  lsProcessedSCData <- processSCData( matCounts=matCounts,
    matPiConstPredictors=matPiConstPredictors,
    vecPseudotime=vecPseudotime,
    scaSmallRun=scaSmallRun,
    strMuModel=strMuModel,
    strDispModel=strDispModel,
    scaWindowRadius=scaWindowRadius,
    boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
    boolCoEstDispMean=boolCoEstDispMean,
    dirOut=dirOut )
  matCountsProc <- lsProcessedSCData$matCountsProc
  matCountsProcFull <- lsProcessedSCData$matCountsProcFull
  vecPseudotimeProc <- lsProcessedSCData$vecPseudotimeProc
  dirOut <- lsProcessedSCData$dirOut
  
  save(matCountsProc,file=file.path(dirOut,"LineagePulse_matCountsProc.RData"))
  save(matCountsProcFull,file=file.path(dirOut,"LineagePulse_matCountsProcFull.RData"))
  save(vecPseudotimeProc,file=file.path(dirOut,"LineagePulse_vecPseudotimeProc.RData"))
  # Clear memory
  rm(matCounts)
  rm(vecPseudotime)
  
  # Save input parameters into list (not data files)
  lsInputParam <- list( scaKClusters=scaKClusters,
    scaSmallRun=scaSmallRun,
    boolClusterInPseudotime=boolClusterInPseudotime,
    scaWindowRadius=scaWindowRadius,
    boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
    boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
    strMuModel=strMuModel,
    strDispModel=strDispModel,
    boolPlotZINBfits=boolPlotZINBfits,
    boolValidateZINBfit=boolValidateZINBfit,
    scaMaxEstimationCycles=scaMaxEstimationCycles,
    scaNProc=scaNProc,
    boolVerbose=boolVerbose,
    boolSuperVerbose=boolSuperVerbose,
    dirOut=dirOut )
  save(lsInputParam,file=file.path(dirOut,"LineagePulse_lsInputParam.RData"))
  rm(lsInputParam)
  
  # 2. Inialise parallelisation
  # Create log directory for parallelisation output
  if(boolBPlog){
    dir.create(file.path(dirOut, "BiocParallel_logs"), showWarnings = FALSE)
    dirBPLogs <- file.path(dirOut, "BiocParallel_logs")
  }
  
  print(paste0("Register parallelisation parameters: ", scaNProc, " threads."))
  # Set the parallelisation environment in BiocParallel:
  if(scaNProc > 1){
    # Set worker time out to 60*60*24*7 (7 days)
    # For single machine (FORK) cluster
    register(MulticoreParam(workers=scaNProc)) 
    #timeout=60*60*24*7,
    #log=boolBPlog, 
    #threshold="INFO", 
    #logdir=dirBPLogs))
    # For multiple machine (SOCK) cluster
    #register(SnowParam(workers=scaNProc, timeout=60*60*24*7))
    # For debugging in serial mode
  } else {
    register(SerialParam())
  }
  
  # 3. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    if(boolClusterInPseudotime){
      # Cluster in pseudotime
      lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotimeProc,
        Kexternal=scaKClusters)
    } else {
      # Take observation time points as clusters
      print("Chose grouping by given time points.")
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
  save(lsResultsClustering,file=file.path(dirOut,"LineagePulse_lsResultsClustering.RData"))
  print(paste("Time elapsed during clustering: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  
  # 4. Compute size factors
  print("4. Compute size factors:")
  vecSizeFactors <- calcSizeFactors(matCountsProcFull)
  save(vecSizeFactors,file=file.path(dirOut,"LineagePulse_vecSizeFactors.RData"))
  # Clear memory
  rm(matCountsProcFull)
  
  # 5. Fit ZINB mixture model for both H1 and H0.
  print("5. Fit ZINB mixture model for both H1 and H0.")
  tm_fitmm <- system.time({
    lsZINBFit <- fitZINB( matCountsProc=matCountsProc,
      matPiConstPredictors=matPiConstPredictors,
      lsResultsClustering=lsResultsClustering,
      vecSizeFactors=vecSizeFactors,
      scaWindowRadius=scaWindowRadius,
      boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
      boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
      boolCoEstDispMean=boolCoEstDispMean,
      strMuModel=strMuModel,
      strDispModel=strDispModel,
      vecPseudotime=vecPseudotimeProc,
      scaMaxEstimationCycles=scaMaxEstimationCycles,
      boolVerbose=boolVerbose,
      boolSuperVerbose=boolSuperVerbose )
    lsMuModelH1 <- lsZINBFit$lsMuModelH1
    lsDispModelH1 <- lsZINBFit$lsDispModelH1
    lsMuModelH0 <- lsZINBFit$lsMuModelH0
    lsDispModelH0 <- lsZINBFit$lsDispModelH0
    lsDropModel <- lsZINBFit$lsDropModel
    lsFitZINBReporters <- lsZINBFit$lsFitZINBReporters
  })
  save(lsMuModelH1,file=file.path(dirOut,"LineagePulse_lsMuModelH1.RData"))
  save(lsDispModelH1,file=file.path(dirOut,"LineagePulse_lsDispModelH1.RData"))
  save(lsMuModelH0,file=file.path(dirOut,"LineagePulse_lsMuModelH0.RData"))
  save(lsDispModelH0,file=file.path(dirOut,"LineagePulse_lsDispModelH0.RData"))
  save(lsDropModel,file=file.path(dirOut,"LineagePulse_lsDropModel.RData"))
  save(lsFitZINBReporters,file=file.path(dirOut,"LineagePulse_lsFitZINBReporters.RData"))
  print(paste("Time elapsed during ZINB fitting: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  
  # 6. Compute posterior probability of drop-out under
  # alternative model.
  #matZH1 <- calcProbNB( matMu=matMuH1,
  #  matDispersions=matDispersionsH1,
  #  matDropout=matDropoutH1,
  #  matboolZero= matCountsProc==0,
  #  matboolNotZeroObserved= (!is.na(matCountsProc) & matCountsProc>0),
  #  scaWindowRadius=scaWindowRadius )
  #save(matZH1,file=file.path(dirOut,"LineagePulse_matZH1.RData"))

  # 7. Plot ZINB fits to data.
  #if(boolPlotZINBfits){
  #  tm_plotZINBfits <- system.time({
  #    graphics.off()
  #    print("7. Plot ZINB fits to data.")
  #    vecZINBfitPlotIDs <- rownames(matCountsProc)
  #    vecZINBfitPlotIDs <- vecZINBfitPlotIDs[1:min(100,length(vecZINBfitPlotIDs))]
  #    vecMuH0 <- matMuH0[,1]
  #    vecDispersionsH0 <- matDispersionsH0[,1]
  #    plotZINBfits( vecGeneIDs=vecZINBfitPlotIDs, 
  #      matCounts=matCountsProc,
  #      vecMu=vecMuH0, 
  #      vecDispersions=vecDispersionsH0,
  #      matProbNB=1-matZH1,
  #      scaWindowRadius=scaWindowRadius,
  #      strPDFname="LineagePulse_ZINBfits.pdf" )
  #  })
  #  print(paste("Time elapsed during plotting of ZINP fits: ",
  #    round(tm_plotZINBfits["elapsed"]/60,2)," min",sep=""))
  #}
  
  # 8. Differential expression analysis:
  print("8. Differential expression analysis:")
  tm_deanalysis_mf <- system.time({
    dfDEAnalysis <- runDEAnalysis(
      matCountsProc = matCountsProc,
      vecSizeFactors=vecSizeFactors,
      lsMuModelH1=lsMuModelH1,
      lsDispModelH1=lsDispModelH1,
      lsMuModelH0=lsMuModelH0,
      lsDispModelH0=lsDispModelH0,
      lsDropModel=lsDropModel,
      scaKbyGeneH1=lsFitZINBReporters$scaKbyGeneH1,
      scaKbyGeneH0=lsFitZINBReporters$scaKbyGeneH0,
      scaWindowRadius=scaWindowRadius )
  })
  save(dfDEAnalysis,file=file.path(dirOut,"LineagePulse_dfDEAnalysis.RData"))
  print(paste("Time elapsed during differential expression analysis: ",
    round(tm_deanalysis_mf["elapsed"]/60,2)," min",sep=""))
  
  # 9. Generate validation metrics for inferred ZINB fits
  if(boolValidateZINBfit){
    print("9. Generate ZINB model fit validation metrics:")
    suppressWarnings( validateOutput(dirOutLineagePulse=dirOut,
      dirOutValidation=dirOut) )
  }
  
  print("LineagePulse complete.")
  return(dfDEAnalysis=dfDEAnalysis)
}