################################################################################
#######################     LineagePulse package     ###########################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

source("~/gitDevelopment/LineagePulse/R/main_LineagePulse.R")

################################################################################
### Main function: Fit mixture model to data
################################################################################

#' Compute BIC
#' 
#' @param n: (scalar) number of observations
#' @param k: (scalar) number of model parameters
#' @param ll: (scalar) loglikelihood of model
#' 
#' @return (scalar) BIC score
#' 
#' @author David Sebastian Fischer
calcBIC <- function(n,k,ll){
  return(log(n)*k-2*ll)
}

#' Compute AIC
#' 
#' @param k: (scalar) number of model parameters
#' @param ll: (scalar) loglikelihood of model
#' 
#' @return (scalar) AIC score
#'
#' @author David Sebastian Fischer
calcAIC <- function(k,ll){
  return(2*k-2*ll)
}

#' Perform drop-out model selection
#' 
#' Fit multiple drop-out models based on the
#' same mean and dispersion model and compute
#' likelihood-based model selection criteria:
#' P-value from loglikelihood-ratio test (LRT),
#' BIC and AIC.
#'
#' @author David Sebastian Fischer
runDropOutModelSelection <- function(
  matCounts,
  dfAnnotation,
  vecConfoundersMu=NULL,
  vecConfoundersDisp=NULL,
  strMuModel="constant",
  strDispModel="constant",
  lsstrDropModel,
  lsstrDropFitGroup,
  lsmatPiConstPredictors,
  boolCheckLRT=TRUE,
  vecNormConstExternal=NULL,
  scaKClusters=NULL,
  scaMaxEstimationCycles=20,
  scaNProc=1,
  boolVerbose=TRUE,
  boolSuperVerbose=FALSE){
  
  ### 1. Check drop-out model list input: 
  # I) All drop-out model lists must be of same length
  # and elements must be named in the same way
  if(length(lsstrDropModel)!=length(lsstrDropFitGroup) |
     length(lsstrDropModel)!=length(lsmatPiConstPredictors)){
    stop(paste0(
      "ERROR IN INPUT TO compareDropOutModels: ",
      "Length of all lists defining drop-out models must be the same: ",
      "lsstrDropModel, lsstrDropFitGroup and lsmatPiConstPredictors.",
      "Each position in these lists corresponds to one drop-out model."))
  }
  if(!all(names(lsstrDropModel)==names(lsstrDropFitGroup)) |
     !all(names(lsstrDropModel)==names(lsmatPiConstPredictors))){
    stop(paste0(
      "ERROR IN INPUT TO compareDropOutModels: ",
      "The vector of elments names in all lists",
      "defining drop-out models must be the same: ",
      "lsstrDropModel, lsstrDropFitGroup and lsmatPiConstPredictors.",
      "Each position in these lists corresponds to one drop-out model."))
  }
  if(is.null(names(lsstrDropModel))){
    stop(paste0(
      "ERROR IN INPUT TO compareDropOutModels: ",
      "Name lists of drop-out models."))
  }
  # II) Check that drop-out models are nested within each other
  # as they appear in list:
  # Model with minimal degrees of freedom at end of list.
  # This is necessary if LRTs are used.
  if(boolCheckLRT){
    # a) Check lsstrDropFitGroup
    # Find first occurence of ForAllCells and check that all following
    # model are ForAllCells
    if(any(unlist(lsstrDropFitGroup)!="PerCell")){
      idxFirst_ForAllCells <- match("ForAllCells", unlist(lsstrDropFitGroup))
      vecidx_PerCell <- which("PerCell" %in% unlist(lsstrDropFitGroup))
      if(any(vecidx_PerCell>idxFirst_ForAllCells)){
        stop(paste0(
          "ERROR IN INPUT TO compareDropOutModels: ",
          "Drop-out models must be nested within each other ",
          "with the model with the highest number of degrees ",
          "of freedom appearing first. ",
          "That implies that all models with strDropFitGroup=ForAllCells ",
          "must be preceded by all models with strDropFitGroup=PerCell."))
      }
    }
    # b) Check lsstrDropModel
    if(any(unlist(lsstrDropModel)=="logistic_ofMu")){
      idxFirst_logistic <- match("logistic", unlist(lsstrDropFitGroup))
      vecidx_logistic_ofMu <- which("logistic_ofMu" %in% unlist(lsstrDropFitGroup))
      if(any(vecidx_logistic_ofMu>idxFirst_logistic)){
        stop(paste0(
          "ERROR IN INPUT TO compareDropOutModels: ",
          "Drop-out models must be nested within each other ",
          "with the model with the highest number of degrees ",
          "of freedom appearing first. ",
          "That implies that all models with strDropModel=logistic ",
          "must be preceded by all models with strDropModel=logistic_ofMu."))
      }
    }
    # c) Check lsmatPiConstPredictors
    vecboolPiConstPredictorsIsNull <- 
      sapply(lsmatPiConstPredictors, function(m) is.null(m) )
    lsmatPiConstPredictorsNotNull <- 
      lsmatPiConstPredictors[!vecboolPiConstPredictorsIsNull]
    if(any(vecboolPiConstPredictorsIsNull)){
      # Check NULL
      idxFirst_null <- which(vecboolPiConstPredictorsIsNull)[1]
      vecidx_notnull <- which(!vecboolPiConstPredictorsIsNull)
      if(any(vecidx_notnull>idxFirst_null)){
        stop(paste0(
          "ERROR IN INPUT TO compareDropOutModels: ",
          "Drop-out models must be nested within each other ",
          "with the model with the highest number of degrees ",
          "of freedom appearing first. ",
          "That implies that all models with matPiConstPredictors=NULL ",
          "must be preceded by all models with matPiConstPredictors!=NULL"))
      }
      # Check columns
      if(sum(!vecboolPiConstPredictorsIsNull)>1){
        for(i1 in seq(1,(length(lsmatPiConstPredictorsNotNull)-1))){
          for(i2 in seq(i1+1,length(lsmatPiConstPredictorsNotNull))){
            for(c2 in seq(1, dim(lsmatPiConstPredictorsNotNull[[i2]])[2])){
              boolColMatch <- any(sapply(
                seq(1, dim(lsmatPiConstPredictorsNotNull[[i2]])[2]), 
                function(c1) all(c1==c2) ))
              if(!boolColMatch){
                if(any(vecidx_notnull>idxFirst_null)){
                  stop(paste0(
                    "ERROR IN INPUT TO compareDropOutModels: ",
                    "Drop-out models must be nested within each other ",
                    "with the model with the highest number of degrees ",
                    "of freedom appearing first. lsmatPiConstPredictorsNotNull: ",
                    "Column ", c2, " of model ",
                    vecidx_notnull[i2], " does not occur ",
                    "in model ", vecidx_notnull[i1], "."))
                }
              }
            }
          }
        }
      }
    }
  }
  
  ### 2. Data preprocessing
  vecAllGenes <- rownames(matCounts)
  lsProcessedSCData <- processSCData(
    matCounts=matCounts,
    dfAnnotation=dfAnnotation,
    vecConfoundersMu=vecConfoundersMu,
    vecConfoundersDisp=vecConfoundersDisp,
    matPiConstPredictors=lsmatPiConstPredictors[[1]],
    vecNormConstExternal=vecNormConstExternal,
    strMuModel=strMuModel,
    strDispModelFull=strDispModelFull,
    strDispModelRed=NULL,
    scaMaxEstimationCycles=scaMaxEstimationCycles,
    boolVerbose=boolVerbose,
    boolSuperVerbose=boolSuperVerbose)
  objectLineagePulse <- lsProcessedSCData$objectLineagePulse
  vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
  # Have to handle matPiConstPredictorsProc individually
  lsmatPiConstPredictorsProc <- lapply(lsmatPiConstPredictors, function(mat){
    return(mat[rownames(objectLineagePulse@matCountsProc),,drop=FALSE])
  })
  
  strMessage <- paste0("The following negative binomial model expression model \n",
                       "is used for all drop-out models: ",
                       "mu=",strMuModel, ", dispersion=", strDispModel, ".")
  objectLineagePulse@strReport <- 
    paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) message(strMessage)
  
  # Clear memory
  rm(matCounts)
  rm(dfAnnotation)
  rm(lsProcessedSCData)
  gc()
  
  # Inialise parallelisation
  # Set the parallelisation environment in BiocParallel:
  if(scaNProc > 1){
    register(MulticoreParam(workers=scaNProc)) 
    #register(SnowParam(workers=scaNProc, timeout=60*60*24*7))
  } else {
    # For debugging in serial mode
    register(SerialParam())
  }
  
  ### 3. Compute normalisation constants
  strMessage <- paste0("--- Compute normalisation constants:")
  objectLineagePulse@strReport <- 
    paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  objectLineagePulse <- calcNormConst(
    objectLineagePulse=objectLineagePulse,
    vecNormConstExternal=vecNormConstExternalProc)
  
  ### 4. Fit list of drop-out models
  scaNGenes <- dim(objectLineagePulse@matCountsProc)[1]
  scaNCells <- dim(objectLineagePulse@matCountsProc)[2]
  scaNModels <- length(lsstrDropModel)
  
  vecModelsToFit <- names(lsstrDropModel)
  lsFits <- list()
  for(m in vecModelsToFit){
    strMessage <- paste0(
      "### Fit drop-out model ", m, 
      " (", match(m, vecModelsToFit),"/", length(vecModelsToFit), ").")
    objectLineagePulse@strReport <- 
      paste0(objectLineagePulse@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_modelfit <- system.time({
      lsFitsModelM <- fitZINB(
        matCounts=objectLineagePulse@matCountsProc,
        dfAnnotation=objectLineagePulse@dfAnnotationProc,
        vecConfoundersDisp=vecConfoundersDisp,
        vecConfoundersMu=vecConfoundersMu,
        vecNormConst=objectLineagePulse@vecNormConst,
        lsDropModel=NULL,
        strMuModel=strMuModel,
        strDispModel=strDispModel,
        strDropModel=lsstrDropModel[[m]],
        strDropFitGroup=lsstrDropFitGroup[[m]],
        matPiConstPredictors=lsmatPiConstPredictorsProc[[m]],
        boolVerbose=boolVerbose,
        boolSuperVerbose=boolSuperVerbose)
      
      vecLogLikFit <- evalLogLikMatrix(
        matCounts=objectLineagePulse@matCountsProc,
        lsMuModel=lsFitsModelM$lsMuModel,
        lsDispModel=lsFitsModelM$lsDispModel, 
        lsDropModel=lsFitsModelM$lsDropModel,
        matWeights=NULL )
    })
    strMessage <- paste0(
      "Fitted ZINB model with drop-out model ", m, 
      " with ll of  ", sum(vecLogLikFit), 
      " in ", round(tm_modelfit["elapsed"]/60,2)," min.")
    objectLineagePulse@strReport <- 
      paste0(lsFitsModelM$strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    lsFits[[m]]$lsMuModel <- lsFitsModelM$lsMuModel
    lsFits[[m]]$lsDispModel <- lsFitsModelM$lsDispModel
    lsFits[[m]]$lsDropModel <- lsFitsModelM$lsDropModel
    lsFits[[m]]$strReport <- lsFitsModelM$strReport
    lsFits[[m]]$vecLogLikFit <- vecLogLikFit
    if(lsstrDropFitGroup[[m]]=="PerCell"){
      lsFits[[m]]$scaDF_AllDropModels <- scaNCells*
        dim(lsFitsModelM$lsDropModel$matDropoutLinModel)[2]
    } else if(lsstrDropFitGroup[[m]]=="ForAllCells"){
      lsFits[[m]]$scaDF_AllDropModels <- 
        dim(lsFitsModelM$lsDropModel$matDropoutLinModel)[2]
    } else {
      stop(paste0("ERROR IN runDropOutModelSelection: ",
                  "lsstrDropFitGroup[[m]]=", lsstrDropFitGroup[[m]], 
                  " not recognised."))
    }
    
    rm(lsFitsModelM)
  }
  
  ### 5. Model selection
  strMessage <- paste0("### Model selection:")
  objectLineagePulse@strReport <- 
    paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  # Compute total degrees of freedom
  scaDF_MuModel <- lsFits[[1]]$lsMuModel$lsMuModelGlobal$scaDegFreedom
  scaDF_DispModel <- lsFits[[1]]$lsDispModel$lsDispModelGlobal$scaDegFreedom
  
  # a) loglikelihood ratio test
  vecDF_AllDropModels  <- sapply(lsFits, function(lsFit) lsFit$scaDF_AllDropModels)
  vecDF <- scaNGenes*(scaDF_MuModel+scaDF_DispModel)+vecDF_AllDropModels
  vecLL <- sapply(lsFits, function(lsFit) sum(lsFit$vecLogLikFit) )
  matLL <- do.call(cbind, lapply(lsFits, function(lsFit) sum(lsFit$vecLogLikFit) ))
  if(boolCheckLRT){
    strMessage <- paste0("# 1. Loglikelihood ratio tests")
    objectLineagePulse@strReport <- 
      paste0(objectLineagePulse@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    matPvalLRT <- matrix(NA, nrow=scaNModels, ncol=scaNModels)
    for(i1 in seq(1,scaNModels)){
      for(i2 in seq(i1,scaNModels)){
        scaDeviance <- 2*(vecLL[i1] - vecLL[i2])
        scaDeltaDF <- vecDF[i1] - vecDF[i2]
        matPvalLRT[i1,i2] <- pchisq(scaDeviance,
                                    df=scaDeltaDF,
                                    lower.tail=FALSE)
      }
    }
    rownames(matPvalLRT) <- names(lsFits)
    colnames(matPvalLRT) <- names(lsFits)
  } else { matPvalLRT <- NULL }
  
  # b) BIC
  strMessage <- paste0("# 2. BIC")
  objectLineagePulse@strReport <- 
    paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  vecBIC <- sapply(seq(1,scaNModels), function(i){
    calcBIC(n=scaNGenes*scaNCells,
            k=vecDF[i],
            ll=vecLL[i])
  })
  
  # c) AIC
  strMessage <- paste0("# 3. AIC")
  objectLineagePulse@strReport <- 
    paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  vecAIC <- sapply(seq(1,scaNModels), function(i){
    calcAIC(k=vecDF[i],
            ll=vecLL[i])
  })
  
  return(list(
    vecDF=vecDF,
    vecLL=vecLL,
    matLL=matLL,
    matPvalLRT=matPvalLRT,
    vecBIC=vecBIC,
    vecAIC=vecAIC,
    objectLineagePulse=objectLineagePulse,
    lsFits=lsFits
  ))
}