################################################################################
#######################     LineagePulse package     ###########################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

source("~/gitDevelopment/LineagePulse/R/main_MixtureModel.R")
source("~/gitDevelopment/LineagePulse/R/main_LineagePulse.R")

################################################################################
### Main function: RSA
################################################################################

#' Responsive subset analysis (RSA) wrapper
#' 
#' This function calls the mixture model estimation wrapper in line
#' with the responsive subset analysis (RSA) scheme. An experienced user
#' can directly call runMixtureModel as well.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runRSA <- function(
  matCounts,
  dfAnnotation,
  vecNCentroidsPerPop=NULL,
  vecH0Pop=NULL,
  vecConfounders,
  scaNMixtures,
  matPiConstPredictors=NULL,
  vecNormConstExternal=NULL,
  strDispModel="constant",
  strDropModel="logistic_ofMu",
  strDropFitGroup="PerCell",
  scaMaxEstimationCyclesEMlike=20,
  scaMaxEstimationCyclesDropModel=20,
  scaNProc=1,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE ){
  
  # Wrap call to runMixtureModel
  objectLineagePulseMM <- runMixtureModel(
    matCounts=matCounts,
    dfAnnotation=dfAnnotation,
    vecConfounders=vecConfounders,
    boolFixedPopulations=TRUE,
    vecNCentroidsPerPop=vecNCentroidsPerPop,
    vecH0Pop=vecH0Pop,
    scaNMixtures=scaNMixtures,
    matPiConstPredictors=matPiConstPredictors,
    vecNormConstExternal=vecNormConstExternal,
    strDispModel=strDispModel,
    strDropModel=strDropModel,
    strDropFitGroup=strDropFitGroup,
    scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
    scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
    scaNProc=scaNProc,
    boolVerbose=boolVerbose,
    boolSuperVerbose=boolSuperVerbose )
  
  return(objectLineagePulseMM)
}