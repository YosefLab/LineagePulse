calcBIC <- function(){
  
}

compareDropOutModels <- function(
  matCounts,
  dfAnnotation,
  vecConfounders=NULL,
  strMuModel="impulse",
  strDispModel="constant",
  lsstrDropModel="logistic_ofMu",
  lsstrDropFitGroup="PerCell",
  lsmatPiConstPredictors=NULL,
  vecNormConstExternal=NULL,
  scaKClusters=NULL,
  boolEstimateNoiseBasedOnH0=FALSE,
  boolVecWindowsAsBFGS=FALSE,
  boolCoEstDispMean=TRUE,
  scaMaxEstimationCycles=20,
  scaNProc=1,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE){
  
  # Fit list of drop-out models
  # Model selection
  # loglikelihood ratio test
  # BIC
}