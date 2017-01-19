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
source("srcLineagePulse_decompressParameters.R")
source("srcLineagePulse_evalImpulseModel.R")
source("srcLineagePulse_evalDropoutModel.R")
source("srcLineagePulse_evalLogLikZINB.R")
source("srcLineagePulse_fitZINB.R")
source("srcLineagePulse_fitZINB_fitMean.R")
source("srcLineagePulse_fitZINB_fitDispersion.R")
source("srcLineagePulse_fitZINB_cofitMeanDispersion.R")
source("srcLineagePulse_fitZINB_fitDropout.R")
source("srcLineagePulse_initialiseImpulseParameters.R")
source("srcLineagePulse_plotComparativeECDF.R")
source("srcLineagePulse_plotGene.R")
source("srcLineagePulse_plotPseudotimeClustering.R")
source("srcLineagePulse_processSCData.R")
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
### Main function: Fit mixture model to data
################################################################################

#' Mixture model wrapper: Differential expression and subpopulation detection
#' 
#' This function performs all steps of mixture model fitting for you.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runMixtureModel <- function(matCounts,
                            vecFixedAssignments=NULL,
                            matPiConstPredictors=NULL,
                            vecNormConstExternal=NULL,
                            strDispModel="constant",
                            scaMaxEstimationCyclesEMlike=20,
                            scaMaxEstimationCyclesDropModel=20,
                            scaNProc=1,
                            boolVerbose=TRUE,
                            boolSuperVerbose=FALSE ){
  
  print("LineagePulse v1.0")
  print("Mixture model v1.0")
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  vecAllGenes <- rownames(matCounts)
  lsProcessedSCData <- processSCDataMixture( matCounts=matCounts,
                                             vecFixedAssignments=vecFixedAssignments,
                                             matPiConstPredictors=matPiConstPredictors,
                                             vecNormConstExternal=vecNormConstExternal,
                                             strDispModel=strDispModel,
                                             scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
                                             scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel )
  matCountsProc <- lsProcessedSCData$matCountsProc
  vecFixedAssignmentsProc <- lsProcessedSCData$vecFixedAssignmentsProc
  vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
  matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
  
  # Clear memory
  rm(matCounts)
  rm(matPiConstPredictors)
  rm(lsProcessedSCData)
  
  # X. Inialise parallelisation
  # Create log directory for parallelisation output
  if(boolBPlog){
    dir.create(file.path(dirOut, "BiocParallel_logs"), showWarnings = FALSE)
    dirBPLogs <- file.path(dirOut, "BiocParallel_logs")
  }
  print(paste0("Parallelisation: ", scaNProc, " threads."))
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
  
  
  # 2. Compute normalisation constants
  print("2. Compute normalisation constants:")
  vecNormConst <- calcNormConst(matCountsProc,
                                vecNormConstExternal=vecNormConstExternalProc)
  
  # 3. Fit ZINB mixture and null model.
  print("3. Fit ZINB mixture and null model.")
  tm_fitmm <- system.time({
    objectLineagePulseMM <- fitMixtureZINBModel(matCounts=matCountsProc,
                                                vecFixedAssignments=vecFixedAssignmentsProc,
                                                vecNormConst=vecNormConst,
                                                scaNMixtures=scaNMixtures,
                                                strDispModel=strDispModel,
                                                scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
                                                scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike)
  })
  objectLineagePulseMM@vecAllGenes <- vecAllGenes
  objectLineagePulseMM@vecFixedAssignments <- vecFixedAssignments
  print(paste("Time elapsed during ZINB mixture and null model fitting: ",round(tm_fitmm["elapsed"]/60,2),
              " min",sep=""))
  
  # 4. Differential expression analysis:
  print("4. Differential expression analysis:")
  tm_deanalysis_mf <- system.time({
    objectLineagePulseMM <- runDEAnalysis(
      matCountsProc = matCountsProc,
      objectLineagePulse=objectLineagePulseMM,
      scaWindowRadius=NULL )
  })
  print(paste("Time elapsed during differential expression analysis: ",
              round(tm_deanalysis_mf["elapsed"]/60,2)," min",sep=""))
  
  print("LineagePulse complete.")
  return(objectLineagePulseMM)
}