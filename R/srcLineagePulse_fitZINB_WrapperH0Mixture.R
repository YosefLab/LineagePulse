#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++  M-step to fit mixture model conditioned on assignments  ++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit ZINB H0 mixture model
#' 
#' This function performs the H0 model fitting necessary for differential
#' expression analysis in the context of RSA.
#' If the H0 is a single mixture, this function returns the input obejct
#' because the constant model used to estimate the drop-out model is already
#' the MLE H0!
#' Takes previously fit drop-out model (same as H1).
#'  
#' @export
fitH0MixtureZINBModel <- function(objectLineagePulse,
                                  vecidxMixturesH0=NULL,
                                  boolVerbose=TRUE,
                                  boolSuperVerbose=FALSE){
  
  if(length(vecidxMixturesH0)>1){
    # Save constant fit to separate slot and free null model slot
    objectLineagePulse@lsMuModelConst <- objectLineagePulse@lsMuModelH0
    objectLineagePulse@lsDispModelConst <- objectLineagePulse@lsDispModelH0
    objectLineagePulse@lsMuModelH0 <- NULL
    objectLineagePulse@lsDispModelH0 <- NULL
    
    strMessage <- paste0("### a) Fit H0 mixture models.")
    objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    # Set all parameter subspace combinations to fit
    scaNMixtures <- dim(objectLineagePulse@matWeights)[2]
    scaKH0 <- length(vecidxMixturesH0)
    vecidxMixturesH1 <- setdiff(seq(1,scaNMixtures), vecidxMixturesH0)
    scaKH1 <- length(vecidxMixturesH1)
    matMapsH1 <- matrix(1, nrow=scaKH0^scaKH1, ncol=scaKH1)
    for(i in seq(1,scaKH0^scaKH1)){
      matMapsH1[i,1] <- (i-1) %% (scaKH0) +1
      if(scaKH1 > 1){
        for(k1 in seq(2, scaKH1)){
          matMapsH1[i,2:(k1-1)] <- (matMapsH1[i,2:(k1-1)]-1) %% (scaKH0) +1
          matMapsH1[i,k1] <- (i-1) %/% (scaKH0^(k1-1)) +1
        }
      }
    }
    matMixMaps <- matrix(seq(1, scaNMixtures), nrow=scaKH0^scaKH1, ncol=scaNMixtures, byrow=TRUE)
    matMixMaps[,vecidxMixturesH1] <- matMapsH1
    print(matMixMaps)
    
    scaNGenes <- dim(objectLineagePulse@matCountsProc)[1]
    scaNCells <- dim(objectLineagePulse@matCountsProc)[2]
    scaNH0SubSpaces <- dim(matMixMaps)[1]
    
    # Initialise objects which save fits:
    # Vector which keeps current LL of gene-wise MLE out of the
    # subspaces covered to the given point in the loop.
    vecLLRef <- array(-Inf, scaNGenes)
    # Objects which carry the current MLE fits
    lsMuModelRed <- objectLineagePulse@lsMuModelH1
    lsMuModelRed$matMuModel <- matrix(NA, nrow=scaNGenes, ncol=dim(objectLineagePulse@lsMuModelH1$matMuModel)[2])
    lsMuModelRed$lsmatBatchModel <- lapply(objectLineagePulse@lsMuModelH1$lsmatBatchModel, function(mat){
      matrix(NA, nrow=scaNGenes, ncol=dim(mat)[2])
    })
    lsDispModelRed <- objectLineagePulse@lsDispModelH1
    lsDispModelRed$matDispModel <- matrix(NA, nrow=scaNGenes, ncol=dim(objectLineagePulse@lsDispModelH1$matDispModel)[2])
    
    # Initialise parameters according to reference mixtures mapped
    matMuModelInit <- objectLineagePulse@lsMuModelH1$matMuModel[,vecidxMixturesH0]
    lsmatBatchModelInit <- objectLineagePulse@lsMuModelH1$lsmatBatchModel
    if(objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel=="constant"){
      matDispModelInit <- objectLineagePulse@lsDispModelH1$matDispModel
    } else {
      stop(paste0("Dispersion  model ", 
                  objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel
                  ," not yet coded in fitH0MixtureZINBModel."))
    }
    
    tm_H0fit <- system.time({
      ### A) Perform one optimisation per sub space
      # Update best gene-wise estimator and throw away others
      # to save memory.
      for(l in seq(1,scaNH0SubSpaces)){
        matWeightsRed <- do.call(cbind, lapply(vecidxMixturesH0, function(m){
          apply(as.matrix(objectLineagePulse@matWeights[,matMixMaps[l,]==m]),1,sum)
        }))
        # M-like step: Estimate mixture model parameters
        tm_mstep <- system.time({
          lsZINBFitsRed <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                                   dfAnnotation=objectLineagePulse@dfAnnotationProc,
                                   vecConfounders=objectLineagePulse@vecConfounders,
                                   vecNormConst=objectLineagePulse@vecNormConst,
                                   matWeights=matWeightsRed,
                                   matPiConstPredictors=NULL,
                                   lsDropModel=objectLineagePulse@lsDropModel,
                                   matMuModelInit=matMuModelInit,
                                   lsmatBatchModelInit=lsmatBatchModelInit,
                                   matDispModelInit=matDispModelInit,
                                   strMuModel="MM",
                                   strDispModel=objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel,
                                   scaMaxEstimationCycles=1,
                                   boolVerbose=FALSE,
                                   boolSuperVerbose=FALSE)
          objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport,
                                                 lsZINBFitsRed$strReport)
          vecLLNew <- evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                       lsMuModel=lsZINBFitsRed$lsMuModel,
                                       lsDispModel=lsZINBFitsRed$lsDispModel, 
                                       lsDropModel=objectLineagePulse@lsDropModel,
                                       matWeights=matWeightsRed )
        })
        strMessage <- paste0("# M-step in subspace ",l," complete: ",
                             "loglikelihood of  ", sum(vecLLNew), " in ",
                             round(tm_mstep["elapsed"]/60,2)," min.")
        objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
        if(boolSuperVerbose) print(strMessage)
        
        if(lsZINBFitsRed$boolConvergenceModel !=0){
          strMessage <- paste0("Model estimation did not converge.")
          objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
          if(boolSuperVerbose) print(strMessage)
        }
        # Update gene-wise model fits to MLE
        vecidxNewMLE <- vecLLNew > vecLLRef
        vecLLRef[vecidxNewMLE] <- vecLLNew[vecidxNewMLE]
        if(any(vecidxNewMLE)){
          lsMuModelRed$matMuModel[vecidxNewMLE,] <- lsZINBFitsRed$lsMuModel$matMuModel[vecidxNewMLE,matMixMaps[l,]]
          for(i in seq(1, length(lsZINBFitsRed$lsMuModel$lsmatBatchModel))){
            lsMuModelRed$lsmatBatchModel[[i]][vecidxNewMLE,] <- lsZINBFitsRed$lsMuModel$lsmatBatchModel[[i]][vecidxNewMLE,]
          }
          if(objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel=="constant"){
            lsDispModelRed$matDispModel[vecidxNewMLE,] <- lsZINBFitsRed$lsDispModel$matDispModel[vecidxNewMLE,]
          } else {
            stop(paste0("Dispersion  model ", 
                        objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel
                        ," not yet coded in fitH0MixtureZINBModel."))
          }
        }
      }
      # Update degrees of freedom of null model
      lsMuModelRed$lsMuModelGlobal$scaDegFreedom <- length(vecidxMixturesH0)
      if(objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel=="constant"){
        lsDispModelRed$lsDispModelGlobal$scaDegFreedom <- 1
      } else {
        stop(paste0("Dispersion  model ", 
                    objectLineagePulse@lsDispModelH1$lsDispModelGlobal$strDispModel
                    ," not yet coded in fitH0MixtureZINBModel."))
      }
      # Update fit objects in LineagePulseObject
      objectLineagePulse@lsMuModelH0 <- lsMuModelRed
      objectLineagePulse@lsDispModelH0 <- lsDispModelRed
    })
    strMessage <- paste0("### Finished fitting H0 mixture ZINB model ",
                         "in ", round(tm_H0fit["elapsed"]/60,2)," min.")
    objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
  }
  
  return(objectLineagePulse)
}
