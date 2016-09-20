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
#' @param K: (integer) [Default NULL] Forces number of centroids
#'    in K-means to be K: setting this to not NULL skips model
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
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' @param boolOneDispPerGene: (bool) [Default TRUE]
#'    Whether one negative binomial dispersion factor is fitted
#'    per gene or per gene for each cluster.
#' @param boolDEAnalysisImpulseModel: (bool) [Default TRUE]
#'    Whether to perform differential expression analysis with ImpulseDE2.
#' @param boolDEAnalysisModelFree: (bool) [Default FALSE]
#'    Whether to perform model-free differential expression analysis.
#' @param boolPlotZINBfits: (bool) [Default TRUE]
#'    Whether to plot zero-inflated negative binomial fits to selected genes
#'    and clusters.
#' @param scaMaxiterEM: (integer) [Default 20] Maximium number 
#'    of EM-terations performed in fitZINB(), i.e. during
#'    estimation of the zero-inflated negative binomial 
#'    hyperparameters.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' 
#' @return (list length 2)
#'    \itemize{
#'    \item lsImpulseDE2results: (list length 4)
#'    \itemize{
#'      \item vecDEGenes: (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'      \item dfImpulseResults: (data frame) ImpulseDE2 results.
#'      \item lsImpulseFits: (list) List of matrices which
#'        contain parameter fits and model values for given time course for the
#'        case condition (and control and combined if control is present).
#'        Each parameter matrix is called parameter_'condition' and has the form
#'        (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'        converge_H0]) where beta to t2 are parameters of the impulse
#'        model, mu is the single parameter of the mean model, logL are
#'        log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'        is convergence status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model. Each value matrix is called
#'        value_'condition' and has the form (genes x time points) and contains the
#'        counts predicted by the impulse model at the observed time points.
#'      \item dfDESeq2Results: (NULL) DESeq2 results, DESeq2 is not run within
#'        ImpulseDE2 in singlecell mode.
#'    }
#'    \item dfModelFreeDEAnalysis: (data frame) 
#'        Summary of model-free differential expression analysis.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @export

runLineagePulse <- function(matCounts, 
  vecPseudotime,
  K=NULL,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit = TRUE,
  boolOneDispPerGene = TRUE,
  scaWindowRadius=20,
  boolEstimateNoiseBasedOnH0=TRUE,
  strMuModel="impulse",
  boolPlotZINBfits = FALSE,
  scaMaxiterEM=20,
  nProc=1,
  verbose=TRUE,
  boolSuperVerbose=FALSE ){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  lsProcessedSCData <- processSCData(matCounts=matCounts,
    vecPseudotime=vecPseudotime,
    scaSmallRun=scaSmallRun,
    strMuModel=strMuModel,
    scaWindowRadius=scaWindowRadius
    )
  matCountsProc <- lsProcessedSCData$matCountsProc
  matCountsProcFull <- lsProcessedSCData$matCountsProcFull
  vecPseudotimeProc <- lsProcessedSCData$vecPseudotimeProc
  
  save(matCountsProc,file=file.path(getwd(),"LineagePulse_matCountsProc.RData"))
  save(matCountsProcFull,file=file.path(getwd(),"LineagePulse_matCountsProcFull.RData"))
  save(vecPseudotimeProc,file=file.path(getwd(),"LineagePulse_vecPseudotimeProc.RData"))
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    if(boolPseudotime){
      # Cluster in pseudotime
      lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotimeProc,
        Kexternal=K)
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
  save(lsResultsClustering,file=file.path(getwd(),"LineagePulse_lsResultsClustering.RData"))
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
  save(dfAnnotation,file=file.path(getwd(),"LineagePulse_dfAnnotation.RData"))
  
  # 4. Compute size factors
  print("4. Compute size factors:")
  vecSizeFactors <- computeSizeFactors_LineagePulse(matCountsProcFull)
  save(vecSizeFactors,file=file.path(getwd(),"LineagePulse_vecSizeFactors.RData"))
  
  # 5. Fit ZINB mixture model for both H1 and H0.
  print("5. Fit ZINB mixture model for both H1 and H0.")
  tm_fitmm <- system.time({
    lsZINBFit <- fitZINB( matCountsProc=matCountsProc, 
      lsResultsClustering=lsResultsClustering,
      vecSizeFactors=vecSizeFactors,
      vecSpikeInGenes=NULL,
      boolOneDispPerGene=boolOneDispPerGene,
      scaWindowRadius=scaWindowRadius,
      boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
      strMuModel=strMuModel,
      vecPseudotime=vecPseudotimeProc,
      nProc=nProc,
      scaMaxiterEM=scaMaxiterEM,
      verbose=verbose,
      boolSuperVerbose=boolSuperVerbose )
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
  save(matMuH1,file=file.path(getwd(),"LineagePulse_matMuH1.RData"))
  save(matDispersionsH1,file=file.path(getwd(),"LineagePulse_matDispersionsH1.RData"))
  save(matDropoutH1,file=file.path(getwd(),"LineagePulse_matDropoutH1.RData"))
  save(matMuH0,file=file.path(getwd(),"LineagePulse_matMuH0.RData"))
  save(matDispersionsH0,file=file.path(getwd(),"LineagePulse_matDispersionsH0.RData"))
  save(matDropoutH0,file=file.path(getwd(),"LineagePulse_matDropoutH0.RData"))
  save(matDropoutLinModel,file=file.path(getwd(),"LineagePulse_matDropoutLinModel.RData"))
  save(matImpulseParam,file=file.path(getwd(),"LineagePulse_matImpulseParam.RData"))
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
  save(matZH1,file=file.path(getwd(),"LineagePulse_matZH1.RData"))

  # X. Plot ZINB fits to data.
  if(boolPlotZINBfits & FALSE){
    graphics.off()
    print("X. Plot ZINB fits to data.")
    vecZINBfitPlotIDs <- names(sort(apply(matCountsProc, 1, mean),
      decreasing=TRUE)[1:min(20,dim(matCountsProc)[1])])
    plotZINBfits(vecGeneIDs=vecZINBfitPlotIDs, 
      matCounts=matCountsProc,
      matMuCluster=matMuCluster, 
      vecDispersions=vecDispersions,
      matProbNB=matProbNB,
      vecClusterAssignments=lsResultsClustering$Assignments,
      lsResultsClustering=lsResultsClustering,
      dfAnnotation=dfAnnotation, 
      strPDFname="LineagePulse_ZINBfits.pdf" )
  }
  
  # 7. Differential expression analysis:
  print("7. Differential expression analysis:")
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
  save(dfDEAnalysis,file=file.path(getwd(),"LineagePulse_dfDEAnalysis.RData"))
  print(paste("Time elapsed during differential expression analysis: ",
    round(tm_deanalysis_mf["elapsed"]/60,2)," min",sep=""))
  
  print("Completed LineagePulse.")
  return(list( dfDEAnalysis=dfDEAnalysis,
    matMu=matMuH1 ))
}