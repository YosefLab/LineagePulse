################################################################################
#######################     LineagePulse package     ###########################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

library(BiocParallel)
library(compiler)
library(ggplot2)
library(MASS)
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/R/ImpulseDE2_main.R")
#source("/data/yosef2/users/fischerd/code/ImpulseDE2/R/ImpulseDE2_main.R")

setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/R")
#setwd("/data/yosef2/users/fischerd/code/LineagePulse/R")
source("srcLineagePulse_evalLogLikZINB.R")
source("srcLineagePulse_fitZINB.R")
source("srcLineagePulse_fitZINB_fitMean.R")
source("srcLineagePulse_fitZINB_fitDispersion.R")
source("srcLineagePulse_fitZINB_fitDropout.R")
# Compile function
evalLogLikZINB_LinPulse_comp <- cmpfun(evalLogLikZINB_LinPulse)
evalLogLikSmoothZINB_LinPulse_comp <- cmpfun(evalLogLikSmoothZINB_LinPulse)
evalLogLikMuWindowsZINB_LinPulse_comp <- cmpfun(evalLogLikMuWindowsZINB_LinPulse)
evalLogLikMuVecWindowsZINB_LinPulse_comp <- cmpfun(evalLogLikMuVecWindowsZINB_LinPulse)
evalLogLikMuClustersZINB_LinPulse_comp <- cmpfun(evalLogLikMuClustersZINB_LinPulse)
evalLogLikMuConstZINB_LinPulse_comp <- cmpfun(evalLogLikMuConstZINB_LinPulse)
evalLogLikDispZINB_LinPulse_comp <- cmpfun(evalLogLikDispZINB_LinPulse)
evalLogLikPiZINB_LinPulse_comp <- cmpfun(evalLogLikPiZINB_LinPulse)
source("srcLineagePulse_clusterCellsInPseudotime.R")
source("srcLineagePulse_createAnnotation.R")
source("srcLineagePulse_computeSizeFactors.R")
source("srcLineagePulse_plotZINBfits.R")
source("srcLineagePulse_plotPseudotimeClustering.R")
source("srcLineagePulse_processSCData.R")
source("srcLineagePulse_runDEAnalysis.R")
source("srcLineagePulse_calcProbNB.R")
source("srcLineagePulse_validateOutput.R")
source("srcLineagePulse_validateOutputSimulation.R")
source("srcLineagePulse_plotGene.R")

################################################################################
### Main function
################################################################################

#' LineagePulse wrapper: Differential expression analysis in pseudotime
#' 
#' This function is the wrapper function for the LineagePulse algorithm,
#' which performs data processing, clustering, zero-inflated negative
#' binomial model identification (hyperparameter estimation) and 
#' model-based differential expression analysis with ImpulseDE2 in 
#' the singlecell mode or model-free differential expression analysis.
#' Differential expression is defined as differential expression over
#' time within one condition, PseuoDE does not handle case-control
#' comparisons at the moment.
#' 
#' @details The computational complexity of ImpulseDE2 is linear in the
#' number of genes and linear in the number of cells.
#' \enumerate{
#' \item \textbf{Cluster cells in pseudo-time with K-means:}
#' The number of clusters $K$ is selected based on the gap-statistic.
#' \item \textbf{Hyperparameter estimation:}
#' A zero-inflated negative binomial model is fit to the clusters for each gene with SCONE.
#' Drop-out rates and dispersion factors are retained as hyperparameters.
#' \item \textbf{Differential expression analysis in pseudo time:}
#' \enumerate{
#'    \item \textbf{Model-free:}
#'    A zero-inflated negative binomial model with an overall mean is fit to the data with SCONE (null model).
#'    The fit of the null model is compared to the fit of the alternative model (cluster-wise zero-inflated negative binomial models from hyperparameter estimation)
#'    with a loglikelihood ratio test.
#'    \item \textbf{Model-based:}
#'    ImpulseDE2 (in the batch mode) fits the impulse model to the data based on a zero-inflated negative binomial cost function 
#'    with drop-out rate and dispersion factor set by SCONE.
#'    ImpulseDE2 performs differential expression analysis based on a loglikelihood ratio test.
#'    }
#' }
#' 
#' @aliases LineagePulse
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param scaKCluster: (integer) [Default NULL] Forces number of centroids
#'    in K-means to be K: setting this to an integer (not NULL) skips model
#'    selection in clusterting.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param boolPseudotime: (bool) [Default TRUE]
#'    Whether to treat time as pseudotime or real time. If FALSE,
#'    time points are treated as clusters. This means that the time
#'    of sampling of the single cells is used as their time coordinate:
#'    e.g. 0,24,28,72 hours in the HSMM data set of the Monocle paper.
#' @param boolContPseudotimeFit: (bool) [Default TRUE]
#'    Whether to fit the impulse model to the pseudotime coordinates
#'    of the cells. If false, the pseudotime centroids of the clusters
#'    are chose: Impulse fitting is done based on the same clusters
#'    as hyperparameter estimation.
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
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param boolPlotZINBfits: (bool) [Default TRUE]
#'    Whether to plot zero-inflated negative binomial fits to selected genes
#'    and clusters.
#' @param boolValidateZINBfit: (bool) [Default TRUE]
#'    Whether to generate evaluation metrics and plots
#'    for parameter values of inferred ZINB model.
#' @param scaMaxCycles: (integer) [Default 20] Maximium number 
#'    of estimation cycles performed in fitZINB(). One cycle
#'    contain one estimation of of each parameter of the 
#'    zero-inflated negative binomial model as coordinate ascent.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param verbose: (bool) [Defaul TRUE]
#'    Whether progress of coordinate ascent within fitZINB is 
#'    reported once per iteration.
#' @param boolSuperVerbose: (bool) [Defaul FALSE]
#'    Whether coordinate ascent within fitZINB is followed
#'    step-by-step rather than reporting once per iteration.
#' @param dirOut: (str directory) [Default NULL]
#'    Directory to which detailed output is saved to.
#'    Defaults to current working directory if NULL.
#' @param boolMemorySaving: (bool) [Default FALSE]
#'    Whether LineagePulse is run in memory saving mode:
#'    In memory saving mode, the communication of files between
#'    function is partially exported to writing and reading
#'    files with fixed names from dirOut.
#' 
#' @return (list length 2)
#'    \itemize{
#'      \item   dfDEAnalysis: (data frame genes x reported variables) 
#'    Summary of differential expression analysis, sorted by adj.p:
#'    {Gene: gene ID,
#'    p: raw p-value, 
#'    adj.p: BH corrected p-value, 
#'    loglik_full: loglikelihood of alternative model H1,
#'    loglik_red: loglikelihood of null model H0,
#'    deviance: loglikelihood ratio test statistic (the deviance),
#'    mean_H0: inferred gene-wise mean parameter (constant null model),
#'    dispersion_H0: inferred gene-wise dispersion parameter (constant null model)}
#'      \item   matMu: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      }
#'    
#' @author David Sebastian Fischer
#' 
#' @export

runLineagePulse <- function(matCounts, 
  vecPseudotime,
  scaKClusters=NULL,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit = TRUE,
  scaWindowRadius=NULL,
  boolEstimateNoiseBasedOnH0=FALSE,
  strMuModel="impulse",
  strDispModel = "constant",
  boolPlotZINBfits = FALSE,
  boolValidateZINBfit=TRUE,
  scaMaxCycles=20,
  nProc=1,
  verbose=TRUE,
  boolSuperVerbose=FALSE,
  dirOut=NULL,
  boolMemorySaving=FALSE ){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  lsProcessedSCData <- processSCData( matCounts=matCounts,
    vecPseudotime=vecPseudotime,
    scaSmallRun=scaSmallRun,
    strMuModel=strMuModel,
    scaWindowRadius=scaWindowRadius,
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
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    if(boolPseudotime){
      # Cluster in pseudotime
      lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotimeProc,
        Kexternal=scaKClusters)
      # Plot clustering
      plotPseudotimeClustering(vecPseudotime=vecPseudotimeProc, 
        lsResultsClustering=lsResultsClustering)
    } else {
      # Take observation time points as clusters
      print("Chose grouping by given time points.")
      lsResultsClustering <- list()
      lsResultsClustering[[1]] <- match(vecPseudotimeProc, sort(unique(vecPseudotimeProc)))
      lsResultsClustering[[2]] <- sort( unique(vecPseudotimeProc) )
      lsResultsClustering[[3]] <- length(unique(vecPseudotimeProc))
      names(lsResultsClustering) <- c("Assignments","Centroids","K")
    }
  })
  save(lsResultsClustering,file=file.path(dirOut,"LineagePulse_lsResultsClustering.RData"))
  print(paste("Time elapsed during clustering: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  
  # 3. Create annotation table
  print("3. Create annotation table")
  if(boolContPseudotimeFit){
    # Fit time trace by cell
    dfAnnotation <- createAnnotationByCell(matCounts=matCountsProc,
      vecPseudotime=vecPseudotimeProc)
  } else {
    # Fit time trace by cluster
    dfAnnotation <- createAnnotationByCluster(matCounts=matCountsProc,
      vecPseudotime=vecPseudotimeProc,
      lsResultsClustering=lsResultsClustering)
  }
  save(dfAnnotation,file=file.path(dirOut,"LineagePulse_dfAnnotation.RData"))
  
  # 4. Compute size factors
  print("4. Compute size factors:")
  vecSizeFactors <- computeSizeFactors_LineagePulse(matCountsProcFull)
  save(vecSizeFactors,file=file.path(dirOut,"LineagePulse_vecSizeFactors.RData"))
  # Clear memory
  rm(matCountsProcFull)
  
  # 5. Fit ZINB mixture model for both H1 and H0.
  print("5. Fit ZINB mixture model for both H1 and H0.")
  tm_fitmm <- system.time({
    lsZINBFit <- fitZINB( matCountsProc=matCountsProc, 
      lsResultsClustering=lsResultsClustering,
      vecSizeFactors=vecSizeFactors,
      scaWindowRadius=scaWindowRadius,
      boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
      strMuModel=strMuModel,
      strDispModel=strDispModel,
      vecPseudotime=vecPseudotimeProc,
      nProc=nProc,
      scaMaxCycles=scaMaxCycles,
      verbose=verbose,
      boolSuperVerbose=boolSuperVerbose,
      boolMemorySaving=boolMemorySaving )
    matMuH1  <- lsZINBFit$matMuH1
    matDispersionsH1 <- lsZINBFit$matDispersionsH1
    matDropoutH1 <- lsZINBFit$matDropoutH1
    matMuH0  <- lsZINBFit$matMuH0
    matDispersionsH0 <- lsZINBFit$matDispersionsH0
    matDropoutH0 <- lsZINBFit$matDropoutH0
    matDropoutLinModel <- lsZINBFit$matDropoutLinModel
    matImpulseParam  <- lsZINBFit$matImpulseParam
    lsFitZINBReporters <- lsZINBFit$lsFitZINBReporters
  })
  save(matMuH1,file=file.path(dirOut,"LineagePulse_matMuH1.RData"))
  save(matDispersionsH1,file=file.path(dirOut,"LineagePulse_matDispersionsH1.RData"))
  save(matDropoutH1,file=file.path(dirOut,"LineagePulse_matDropoutH1.RData"))
  save(matMuH0,file=file.path(dirOut,"LineagePulse_matMuH0.RData"))
  save(matDispersionsH0,file=file.path(dirOut,"LineagePulse_matDispersionsH0.RData"))
  save(matDropoutH0,file=file.path(dirOut,"LineagePulse_matDropoutH0.RData"))
  save(matDropoutLinModel,file=file.path(dirOut,"LineagePulse_matDropoutLinModel.RData"))
  save(matImpulseParam,file=file.path(dirOut,"LineagePulse_matImpulseParam.RData"))
  print(paste("Time elapsed during ZINB fitting: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  
  # 6. Compute posterior probability of drop-out under
  # alternative model.
  matZH1 <- calcProbNB( matMu=matMuH1,
    matDispersions=matDispersionsH1,
    matDropout=matDropoutH1,
    matboolZero= matCountsProc==0,
    matboolNotZeroObserved= (!is.na(matCountsProc) & matCountsProc>0),
    scaWindowRadius=scaWindowRadius )
  save(matZH1,file=file.path(dirOut,"LineagePulse_matZH1.RData"))

  # 7. Plot ZINB fits to data.
  if(boolPlotZINBfits){
    tm_plotZINBfits <- system.time({
      graphics.off()
      print("7. Plot ZINB fits to data.")
      vecZINBfitPlotIDs <- rownames(matCountsProc)
      vecZINBfitPlotIDs <- vecZINBfitPlotIDs[1:min(100,length(vecZINBfitPlotIDs))]
      vecMuH0 <- matMuH0[,1]
      vecDispersionsH0 <- matDispersionsH0[,1]
      plotZINBfits( vecGeneIDs=vecZINBfitPlotIDs, 
        matCounts=matCountsProc,
        vecMu=vecMuH0, 
        vecDispersions=vecDispersionsH0,
        matProbNB=1-matZH1,
        scaWindowRadius=scaWindowRadius,
        strPDFname="LineagePulse_ZINBfits.pdf" )
    })
    print(paste("Time elapsed during plotting of ZINP fits: ",
      round(tm_plotZINBfits["elapsed"]/60,2)," min",sep=""))
  }
  
  # 8. Differential expression analysis:
  print("8. Differential expression analysis:")
  tm_deanalysis_mf <- system.time({
    dfDEAnalysis <- runDEAnalysis(
      matCountsProc = matCountsProc,
      vecSizeFactors=vecSizeFactors,
      matMuH1=matMuH1,
      matDispersionsH1=matDispersionsH1,
      matDropoutH1=matDropoutH1,
      matMuH0=matMuH0,
      matDispersionsH0=matDispersionsH0,
      matDropoutH0=matDropoutH0,
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
    suppressWarnings( validateOutput(dirLineagePulseTempFiles=dirOut,
      dirValidationOut=dirOut) )
  }
  
  # TODO: return mu matrix with constant fits and impulse fits
  # if significant.
  
  print("LineagePulse complete.")
  return(list( dfDEAnalysis=dfDEAnalysis,
    matMu=matMuH1 ))
}