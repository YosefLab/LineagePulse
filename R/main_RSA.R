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
#' with the responsive subset analysis (RSA) scheme. An experienced user
#' can directly call runMixtureModel as well.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runRSA <- function(matCounts,
									 dfAnnotation,
									 vecConfounders,
									 scaNMixtures,
									 matPiConstPredictors=NULL,
									 vecNormConstExternal=NULL,
									 strDispModel="constant",
									 scaMaxEstimationCyclesEMlike=20,
									 scaMaxEstimationCyclesDropModel=20,
									 scaNProc=1,
									 boolVerbose=TRUE,
									 boolSuperVerbose=FALSE ){
	
	# Wrap call to runMixtureModel
	objectLineagePulseMM <- runMixtureModel(matCounts=matCounts,
																					dfAnnotation=dfAnnotation,
																					vecConfounders=vecConfounders,
																					boolFixedPopulations=TRUE,
																					scaNMixtures=scaNMixtures,
																					matPiConstPredictors=matPiConstPredictors,
																					vecNormConstExternal=vecNormConstExternal,
																					strDispModel=strDispModel,
																					scaMaxEstimationCyclesEMlike=scaMaxEstimationCyclesEMlike,
																					scaMaxEstimationCyclesDropModel=scaMaxEstimationCyclesDropModel,
																					scaNProc=scaNProc,
																					boolVerbose=boolVerbose,
																					boolSuperVerbose=boolSuperVerbose )
	
	return(objectLineagePulseMM)
}