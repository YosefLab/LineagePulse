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
  scaMaxEstimationCycles=20,
  scaNProc=1,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE){
  
  # Check input: all drop-out model lists must be of same length
  # and elements must be named in the same way
  if(length(lsstrDropModel)!=length(lsstrDropFitGroup) |
     length(lsstrDropModel)!=length(lsmatPiConstPredictors)){
    stop(paste0("ERROR IN INPUT TO compareDropOutModels: ",
                "Length of all lists defining drop-out models must be the same: ",
                "lsstrDropModel, lsstrDropFitGroup and lsmatPiConstPredictors.",
                "Each position in these lists corresponds to one drop-out model."))
  }
  if(!all(names(lsstrDropModel)==names(lsstrDropFitGroup)) |
     !all(names(lsstrDropModel)==names(lsmatPiConstPredictors))){
    stop(paste0("ERROR IN INPUT TO compareDropOutModels: ",
                "The vector of elments names in all lists",
                "defining drop-out models must be the same: ",
                "lsstrDropModel, lsstrDropFitGroup and lsmatPiConstPredictors.",
                "Each position in these lists corresponds to one drop-out model."))
  }
  strMessage <- paste0("### The following negative binomial model expression model ",
                       "is used for all drop-out models:",
                       "mu: ",strMuModel, ", dispersion: ", strDispModel, ".")
  if(boolVerbose) message(strMessage)
  
  # Fit list of drop-out models
  lsModelsToFit <- names(lsstrDropModel)
  for(m in lsModelsToFit){
    strMessage <- paste0("### Fit drop-out model ", m, 
                         " (", match(m, lsstrDropModel),"/", length(lsstrDropModel), ").")
    if(boolVerbose) message(strMessage)
    
    tm_modelfit <- system.time({
      lsFitsModelA <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                              dfAnnotation=objectLineagePulse@dfAnnotationProc,
                              vecConfounders=vecConfounders,
                              vecNormConst=objectLineagePulse@vecNormConst,
                              lsDropModel=NULL,
                              strMuModel=strMuModel,
                              strDispModel=strDispModel,
                              strDropModel=lsstrDropModel[[m]],
                              strDropFitGroup=lsstrDropFitGroup[[m]],
                              boolVerbose=boolVerbose,
                              boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelA <- lsFitsModelA$lsMuModel
    lsDispModelA <- lsFitsModelA$lsDispModel
    lsDropModel <- lsFitsModelA$lsDropModel
    boolConvergenceModelA <- lsFitsModelA$boolConvergenceModel
    vecEMLogLikModelA <- lsFitsModelA$vecEMLogLikModel
    objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport,
                                           lsFitsModelA$strReport)
    rm(lsFitsModelA)
    
    strMessage <- paste0("Finished fitting zero-inflated negative binomial ",
                         "model with drop-out model", m, 
                         "in ", round(tm_modelfit["elapsed"]/60,2)," min.")
    if(boolVerbose) message(strMessage)
  }
  # Model selection
  # loglikelihood ratio test
  # BIC
}