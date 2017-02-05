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
setwd("~/gitDevelopment/LineagePulse/R")

source("main_LineagePulse.R")

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
														dfAnnotation,
														vecConfounders,
														boolFixedPopulations,
                            scaNMixtures,
                            matPiConstPredictors=NULL,
                            vecNormConstExternal=NULL,
                            strDispModel="constant",
                            scaMaxEstimationCyclesEMlike=20,
                            scaMaxEstimationCyclesDropModel=20,
                            scaNProc=1,
                            boolVerbose=TRUE,
                            boolSuperVerbose=FALSE ){
  
  STR_VERSION <- "v0.99"
  
  # 1. Data preprocessing
  vecAllGenes <- rownames(matCounts)
  lsProcessedSCData <- processSCDataMixture( matCounts=matCounts,
  																					 dfAnnotation=dfAnnotation,
  																					 vecConfounders=vecConfounders,
  																					 boolFixedPopulations=boolFixedPopulations,
                                             scaNMixtures=scaNMixtures,
                                             matPiConstPredictors=matPiConstPredictors,
                                             vecNormConstExternal=vecNormConstExternal,
                                             strDispModel=strDispModel,
                                             scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
                                             scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
  																					 boolVerbose=boolVerbose,
  																					 boolSuperVerbose=boolSuperVerbose,
  																					 STR_VERSION=STR_VERSION)
  objectLineagePulseMM <- lsProcessedSCData$objectLineagePulse
  vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
  matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
  
  # Clear memory
  rm(matCounts)
  rm(matPiConstPredictors)
  rm(lsProcessedSCData)
  
  # Inialise parallelisation
  # Set the parallelisation environment in BiocParallel:
  if(scaNProc > 1){
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
  strMessage <- paste0("2. Compute normalisation constants:")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  objectLineagePulseMM <- calcNormConst(objectLineagePulseMM,
                                vecNormConstExternal=vecNormConstExternalProc)
  
  # 3. Fit ZINB mixture and null model.
  strMessage <- paste0("3. Fit ZINB mixture and null model.")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  tm_fitmm <- system.time({
    objectLineagePulseMM <- fitMixtureZINBModel(objectLineagePulse=objectLineagePulseMM,
                                                scaNMixtures=scaNMixtures,
                                                strDispModel=strDispModel,
                                                scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
                                                scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
                                                boolVerbose=boolVerbose,
                                                boolSuperVerbose=boolSuperVerbose )
  })
  
  strMessage <- paste0("Time elapsed during ZINB mixture and null model fitting: ",
                       round(tm_fitmm["elapsed"]/60,2), " min")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  # 4. Differential expression analysis:
  strMessage <- paste0("4. Differential expression analysis:")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  tm_deanalysis_mf <- system.time({
    objectLineagePulseMM <- runDEAnalysis(objectLineagePulse=objectLineagePulseMM )
  })
  
  strMessage <- paste0("Time elapsed during differential expression analysis: ",
              round(tm_deanalysis_mf["elapsed"]/60,2)," min")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  strMessage <- paste0("Finished runMixtureModel().")
  objectLineagePulseMM@strReport <- paste0(objectLineagePulseMM@strReport, strMessage)
  if(boolVerbose) print(strMessage)
  
  return(objectLineagePulseMM)
}