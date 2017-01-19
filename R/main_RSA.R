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

source("main_MixtureModel.R")
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