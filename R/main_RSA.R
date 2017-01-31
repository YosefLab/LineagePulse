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

source("main_MixtureModel.R")

################################################################################
### Main function: RSA
################################################################################

#' Responsive subset analysis (RSA) wrapper
#' 
#' This function calls the mixture model estimation wrapper in line
#' with the responsive subset analysis (RSA) scheme.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runRSA <- function(matCounts,
									 dfAnnotation,
									 vecConfounders,
                   scaNMixtures,
                   vecPopulationID,
                   matPiConstPredictors=NULL,
                   vecNormConstExternal=NULL,
                   strDispModel="constant",
                   scaMaxEstimationCyclesEMlike=20,
                   scaMaxEstimationCyclesDropModel=20,
                   scaNProc=1,
                   boolVerbose=TRUE,
                   boolSuperVerbose=FALSE ){
  
  # Map mixture IDs to numerical vector
  vecPopulationAssignmentsUnique <-  unique(vecPopulationAssignments[!is.na(vecPopulationAssignments)])
  vecFixedAssignments <- match(vecPopulationAssignments, vecPopulationAssignmentsUnique)
  vecPopulationIDs <- array(NA, scaNMixtures)
  vecPopulationIDs[1:length(vecPopulationAssignmentsUnique)] <- vecPopulationAssignmentsUnique
  if(length(vecPopulationAssignmentsUnique) < scaNMixtures) vecPopulationIDs[(length(vecPopulationAssignmentsUnique)+1):scaNMixtures] <- paste0("Population", (length(vecPopulationAssignmentsUnique)+1):scaNMixtures)
  # Wrap call to runMixtureModel
  objectLineagePulseMM <- runMixtureModel(matCounts=matCounts,
  																				dfAnnotation=dfAnnotation,
  																				vecConfounders=vecConfounders,
                                          scaNMixtures=scaNMixtures,
                                          vecFixedAssignments=vecFixedAssignments,
                                          matPiConstPredictors=matPiConstPredictors,
                                          vecNormConstExternal=vecNormConstExternal,
                                          strDispModel=strDispModel,
                                          scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
                                          scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
                                          scaNProc=scaNProc,
                                          boolVerbose=boolVerbose,
                                          boolSuperVerbose=boolSuperVerbose )
  # Name fixed subpopulations
  colnames(objectLineagePulseMM@matWeights) <- vecPopulationIDs
  colnames(objectLineagePulseMM@lsMuModelH1$matMuModel) <- vecPopulationIDs
  
  return(objectLineagePulseMM)
}