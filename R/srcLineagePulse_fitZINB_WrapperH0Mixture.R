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
fitH0MixtureZINBModel <- function(
    objLP,
    vecidxMixturesH0=NULL,
    boolVerbose=TRUE,
    boolSuperVerbose=FALSE){
    
    if(length(vecidxMixturesH0)>1){
        # Save constant fit to separate slot and free null model slot
        objLP@lsMuModelConst <- objLP@lsMuModelH0
        objLP@lsDispModelConst <- objLP@lsDispModelH0
        objLP@lsMuModelH0 <- NULL
        objLP@lsDispModelH0 <- NULL
        
        strMessage <- paste0("### a) Fit H0 mixture models.")
        objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
        if(boolVerbose) print(strMessage)
        
        # Set all parameter subspace combinations to fit
        scaNMixtures <- dim(objLP@matWeights)[2]
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
        
        scaNGenes <- dim(objLP@matCountsProc)[1]
        scaNCells <- dim(objLP@matCountsProc)[2]
        scaNH0SubSpaces <- dim(matMixMaps)[1]
        
        # Initialise objects which save fits:
        # Vector which keeps current LL of gene-wise MLE out of the
        # subspaces covered to the given point in the loop.
        vecLLRef <- array(-Inf, scaNGenes)
        # Objects which carry the current MLE fits
        lsMuModelRed <- objLP@lsMuModelH1
        lsMuModelRed$matMuModel <- matrix(NA, nrow=scaNGenes, ncol=dim(objLP@lsMuModelH1$matMuModel)[2])
        lsMuModelRed$lsmatBatchModel <- lapply(objLP@lsMuModelH1$lsmatBatchModel, function(mat){
            matrix(NA, nrow=scaNGenes, ncol=dim(mat)[2])
        })
        lsDispModelRed <- objLP@lsDispModelH1
        lsDispModelRed$matDispModel <- matrix(NA, nrow=scaNGenes, ncol=dim(objLP@lsDispModelH1$matDispModel)[2])
        
        # Initialise parameters according to reference mixtures mapped
        matMuModelInit <- objLP@lsMuModelH1$matMuModel[,vecidxMixturesH0]
        lsmatBatchModelInitMu <- objLP@lsMuModelH1$lsmatBatchModel
        matDispModelInit <- objLP@lsDispModelH1$matDispModel
        lsmatBatchModelInitDisp <- objLP@lsDispModelH1$lsmatBatchModel
        
        tm_H0fit <- system.time({
            ### A) Perform one optimisation per sub space
            # Update best gene-wise estimator and throw away others
            # to save memory.
            for(l in seq(1,scaNH0SubSpaces)){
                matWeightsRed <- do.call(cbind, lapply(vecidxMixturesH0, function(m){
                    apply(as.matrix(objLP@matWeights[,matMixMaps[l,]==m]),1,sum)
                }))
                # M-like step: Estimate mixture model parameters
                tm_mstep <- system.time({
                    lsZINBFitsRed <- fitZINB(
                        matCounts=objLP@matCountsProc,
                        dfAnnotation=objLP@dfAnnotationProc,
                        vecConfounders=objLP@vecConfounders,
                        vecConfoundersDisp=objLP@vecConfoundersDisp,
                        vecNormConst=objLP@vecNormConst,
                        matWeights=matWeightsRed,
                        matPiConstPredictors=NULL,
                        lsDropModel=objLP@lsDropModel,
                        matMuModelInit=matMuModelInit,
                        lsmatBatchModelInitMu=lsmatBatchModelInitMu,
                        matDispModelInit=matDispModelInit,
                        lsmatBatchModelInitDisp=lsmatBatchModelInitDisp,
                        strMuModel="MM",
                        strDispModel=objLP@lsDispModelH1$lsDispModelGlobal$strDispModel,
                        scaMaxEstimationCycles=1,
                        boolVerbose=FALSE,
                        boolSuperVerbose=FALSE)
                    objLP@strReport <- paste0(objLP@strReport,
                                              lsZINBFitsRed$strReport)
                    vecLLNew <- evalLogLikMatrix(matCounts=objLP@matCountsProc,
                                                 lsMuModel=lsZINBFitsRed$lsMuModel,
                                                 lsDispModel=lsZINBFitsRed$lsDispModel, 
                                                 lsDropModel=objLP@lsDropModel,
                                                 matWeights=matWeightsRed )
                })
                strMessage <- paste0("# M-step in subspace ",l," complete: ",
                                     "loglikelihood of  ", sum(vecLLNew), " in ",
                                     round(tm_mstep["elapsed"]/60,2)," min.")
                objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                if(boolSuperVerbose) print(strMessage)
                
                if(lsZINBFitsRed$boolConvergenceModel !=0){
                    strMessage <- paste0("Model estimation did not converge.")
                    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
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
                    if(objLP@lsDispModelH1$lsDispModelGlobal$strDispModel=="constant"){
                        lsDispModelRed$matDispModel[vecidxNewMLE,] <- lsZINBFitsRed$lsDispModel$matDispModel[vecidxNewMLE,]
                    } else {
                        stop(paste0("Dispersion  model ", 
                                    objLP@lsDispModelH1$lsDispModelGlobal$strDispModel
                                    ," not yet coded in fitH0MixtureZINBModel."))
                    }
                }
            }
            # Update degrees of freedom of null model
            lsMuModelRed$lsMuModelGlobal$scaDegFreedom <- length(vecidxMixturesH0)
            if(objLP@lsDispModelH1$lsDispModelGlobal$strDispModel=="constant"){
                lsDispModelRed$lsDispModelGlobal$scaDegFreedom <- 1
            } else {
                stop(paste0("Dispersion  model ", 
                            objLP@lsDispModelH1$lsDispModelGlobal$strDispModel
                            ," not yet coded in fitH0MixtureZINBModel."))
            }
            # Update fit objects in LineagePulseObject
            objLP@lsMuModelH0 <- lsMuModelRed
            objLP@lsDispModelH0 <- lsDispModelRed
        })
        strMessage <- paste0("### Finished fitting H0 mixture ZINB model ",
                             "in ", round(tm_H0fit["elapsed"]/60,2)," min.")
        objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
        if(boolVerbose) print(strMessage)
    }
    
    return(objLP)
}
