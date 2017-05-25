###  EM-like iteration for mixture model

#' Fit ZINB mixture model
#' 
#' Structure of code:
#' A: Fit a constant model to all cells. The resulting drop-out model
#'    is used for the EM-like iteration. This is also the null model
#'    for the LRT as the full model is fir based on this drop-out model.
#' B: EM-like iteration to fit H1 mixture model.
#'  
#' @export
fitMixtureZINBModel <- function(
    objLP,
    scaNMixtures,
    strDispModel="constant",
    strDropModel="logistic_ofMu",
    strDropFitGroup="PerCell",
    scaMaxEstimationCyclesDropModel=20,
    scaMaxEstimationCyclesEMlike=20,
    boolVerbose=TRUE,
    boolSuperVerbose=FALSE){
    
    ####################################################
    # Internal Numerical Estimation Parameters:
    # Minimim fractional liklihood increment necessary to
    # continue EM-iterations of assignment estimation.
    scaPrecEMAssignments <- 1-10^(-4)
    # Numerical optmisation of impulse model hyperparameters
    MAXIT_BFGS_MM <- 1000 # optim default is 1000
    RELTOL_BFGS_MM <- 10^(-4) # optim default is sqrt(.Machine$double.eps)=1e-8
    # Lowering RELTOL_BFGS_IMPULSE gives drastic run time improvements.
    # Set to 10^(-4) to maintain sensible fits without running far into saturation
    # in the objective (loglikelihood).
    ####################################################
    
    scaNGenes <- dim(objLP@matCountsProc)[1]
    scaNCells <- dim(objLP@matCountsProc)[2]
    
    ### (A) Pre-estimate drop-out model to speed up mixture model fitting
    strMessage <- paste0("### a) Fit constant ZINB model and set drop-out model.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    tm_cycle <- system.time({
        lsZINBFitsRed <- fitZINB(
            matCounts=objLP@matCountsProc,
            dfAnnotation=objLP@dfAnnotationProc,
            vecConfoundersMu=objLP@vecConfoundersMu,
            vecConfoundersDisp=objLP@vecConfoundersDisp,
            vecNormConst=objLP@vecNormConst,
            matWeights=NULL,
            matPiConstPredictors=NULL,
            lsDropModel=NULL,
            strMuModel="constant",
            strDispModel="constant",
            strDropModel=strDropModel,
            strDropFitGroup=strDropFitGroup,
            scaMaxEstimationCycles=scaMaxEstimationCyclesDropModel,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
        lsMuModelRed <- lsZINBFitsRed$lsMuModel
        lsDispModelRed <- lsZINBFitsRed$lsDispModel
        lsDropModel <- lsZINBFitsRed$lsDropModel
        boolConvergenceModelRed <- lsZINBFitsRed$boolConvergenceModel
        vecEMLogLikModelRed <- lsZINBFitsRed$vecEMLogLikModel
        objLP@strReport <- paste0(objLP@strReport,
                                  lsZINBFitsRed$strReport)
        rm(lsZINBFitsRed)
    })
    
    strMessage <- paste0("### Finished fitting H0 constant ZINB model ",
                         "in ", round(tm_cycle["elapsed"]/60,2)," min.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    ### (B) EM-like estimation cycle: fit mixture model
    strMessage <- paste0("### b) EM-like iteration: Fit H1 mixture model")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    # (I) Initialise estimation
    # Initialise expression models as null model estimate:
    # Correct fixed weights (RSA)
    # Set weights to uniform distribution
    matWeights <- matrix(1/scaNMixtures, 
                         nrow=scaNCells, ncol=scaNMixtures)
    if(objLP@boolFixedPopulations){
        vecFixedAssignments <- objLP@dfAnnotationProc$populations
        matWeights[!is.na(vecFixedAssignments),] <- 0
        scaNCentroidsUsed <- 0
        lsvecFixedCentrByPop <- list()
        for(p in seq(1, length(objLP@vecNCentroidsPerPop))){
            vecidxCentroidsOfPop <- seq(scaNCentroidsUsed+1, scaNCentroidsUsed+objLP@vecNCentroidsPerPop[p])
            matWeights[which(objLP@dfAnnotationProc$populations %in% 
                                 names(objLP@vecNCentroidsPerPop[p])), 
                       vecidxCentroidsOfPop] <- 1/objLP@vecNCentroidsPerPop[p]
            lsvecFixedCentrByPop[[match(objLP@vecNCentroidsPerPop[p], 
                                        objLP@vecNCentroidsPerPop)]] <- vecidxCentroidsOfPop
            scaNCentroidsUsed <- scaNCentroidsUsed+p
        }
        names(lsvecFixedCentrByPop) <- names(objLP@vecNCentroidsPerPop)
    } else {
        lsvecFixedCentrByPop <- NULL 
    }
    # Use cells to initialise centroids: Chose cells so that 
    # the overall distance is maximised.
    lsMuModelFull <- lsMuModelRed
    lsMuModelFull$lsMuModelGlobal$strMuModel <- "MM"
    vecidxCellsForCentroids <- initialiseCentroidsFromCells(
        matCounts=objLP@matCountsProc,
        lsvecFixedCentrByPop=lsvecFixedCentrByPop,
        vecAssignPop=objLP@dfAnnotationProc$populations,
        scaN=scaNMixtures)
    # Add pseudo count to not generate error in optim in log space
    lsMuModelFull$matMuModel <- as.matrix(objLP@matCountsProc[
        ,as.double(vecidxCellsForCentroids),sparse=FALSE]) # not a sparseMatrix
    lsMuModelFull$matMuModel[lsMuModelFull$matMuModel < 10^(-5)] <- 10^(-5)
    ###
    #lsMuModelFull$matMuModel <- matrix(
    #    as.vector(lsMuModelFull$matMuModel), 
    #    nrow = dim(lsMuModelFull$matMuModel)[1],
    #    ncol = scaNMixtures,
    #    byrow = FALSE)
    if(!is.null(lsMuModelFull$lsMuModelGlobal$scaNConfounders)){
        for(conf in seq(1, lsMuModelFull$lsMuModelGlobal$scaNConfounders, by=1)) {
            lsMuModelFull$lsmatBatchModel[[conf]] <- 
                matrix(1, nrow=scaNGenes,ncol=lsMuModelFull$lsMuModelGlobal$vecNBatches[conf])
        }
    }
    ###
    
    lsDispModelFull <- lsDispModelRed
    lsDispModelFull$lsDispModelGlobal$scaNMix <- scaNMixtures
    if(strDispModel == "MM") {
        lsDispModelFull$lsDispModelGlobal$strDispModel <- "MM"
        lsDispModelFull$matDispModel <- matrix(
            as.vector(lsDispModelFull$matDispModel), 
            nrow = dim(lsDispModelFull$matDispModel)[1],
            ncol = scaNMixtures,
            byrow = FALSE)
    }
    if(!is.null(lsDispModelFull$lsDispModelGlobal$scaNConfounders)) {
        for(conf in seq(1, lsDispModelFull$lsDispModelGlobal$scaNConfounders, by=1)) {
            lsDispModelFull$lsmatBatchModel[[conf]] <- 
                matrix(1, nrow=scaNGenes,ncol=lsDispModelFull$lsDispModelGlobal$vecNBatches[conf])
        }
    }
    
    # (II) Estimation iteration on full model
    # Set iteration reporters
    scaIter <- 1
    scaLogLikNew <- sum(evalLogLikMatrix(
        matCounts=objLP@matCountsProc,
        lsMuModel=lsMuModelFull,
        lsDispModel=lsDispModelFull, 
        lsDropModel=lsDropModel,
        matWeights=matWeights ))
    scaLogLikOld <- NA
    vecEMLogLikModelFull <- array(NA, scaMaxEstimationCyclesEMlike)
    
    strMessage <- paste0("#  .   Initialised MM: ",
                         "ll   ", scaLogLikNew)
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolSuperVerbose) print(strMessage)
    
    boolMixtureReset <- FALSE
    tm_RSAcycle <- system.time({
        while(scaIter == 1 | 
              boolMixtureReset| 
              (scaLogLikNew > scaLogLikOld*scaPrecEMAssignments & scaIter <= scaMaxEstimationCyclesEMlike)){
            boolMixtureReset <- FALSE
            # E-like step: Estimation of mixture assignments
            tm_estep <- system.time({
                lsWeightFits <- estimateMMAssignmentsMatrix(
                    matCounts=objLP@matCountsProc,
                    dfAnnotation=objLP@dfAnnotationProc,
                    boolFixedPopulations=objLP@boolFixedPopulations,
                    lsvecFixedCentrByPop=lsvecFixedCentrByPop,
                    vecNormConst=objLP@vecNormConst,
                    lsMuModel=lsMuModelFull,
                    lsDispModel=lsDispModelFull,
                    lsDropModel=lsDropModel,
                    matWeights=matWeights,
                    MAXIT_BFGS_MM=MAXIT_BFGS_MM,
                    RELTOL_BFGS_MM=RELTOL_BFGS_MM )
                matWeights <- lsWeightFits$matWeights
            })
            
            scaLogLikTemp <- sum(evalLogLikMatrix(
                matCounts=objLP@matCountsProc,
                lsMuModel=lsMuModelFull,
                lsDispModel=lsDispModelFull, 
                lsDropModel=lsDropModel,
                matWeights=matWeights ))
            
            strMessage <- paste0("# ", scaIter,".   E-step complete: ",
                                 "ll  ", scaLogLikTemp, " in ",
                                 round(tm_estep["elapsed"]/60,2)," min.")
            objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
            if(boolSuperVerbose) print(strMessage)
            
            strMessage <- paste0("# ", scaIter,".   Cumulative mixture weights  ", 
                                 paste(round(apply(matWeights,2,sum),2), collapse=" "))
            objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
            if(boolSuperVerbose) print(strMessage)
            
            if(any(lsWeightFits$vecConvergence !=0 )){
                strMessage <- paste0("Weight estimation did not convergen in ",
                                     sum(lsWeightFits$vecConvergence !=0), " cases.")
                objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                if(boolSuperVerbose) print(strMessage)
            }
            
            # Catch mixture drop-out
            #vecidxAssignedMixture <- apply(matWeights, 1, which.max)
            #vecboolMixtureDropped <- sapply(seq(1,scaNMixtures), function(mixture){
            #  sum(vecidxAssignedMixture==mixture) <= 1
            #})
            vecboolMixtureDropped <- array(FALSE, scaNMixtures)
            if(scaIter>1 & any(vecboolMixtureDropped)){
                boolMixtureReset <- TRUE
                # Reinitialise one mixture component
                # Get badly described cell: Chose via likelihood under of non-dropout mixtures.
                # Don't chose minimum to not catch outlier cells.
                #vecidxBadlyDescrCells <- sort(apply(matWeights, 1, max), index.return=TRUE, decreasing=FALSE)$ix
                matLLCellgivenMixture <- do.call(cbind, lapply(which(!vecboolMixtureDropped), function(m){
                    sapply(seq(1,scaNCells), function(j){
                        vecMuParam <- do.call(c, lapply(seq(1,scaNGenes), function(i){
                            decompressMeansByGene(
                                vecMuModel=lsMuModelFull$matMuModel[i,m],
                                lsvecBatchModel=lapply(lsMuModelFull$lsmatBatchModel, function(mat) mat[i,] ),
                                lsMuModelGlobal=lsMuModelFull$lsMuModelGlobal,
                                vecInterval=j)
                        }))
                        vecDropParam <- decompressDropoutRateByCell(
                            vecDropModel=lsDropModel$matDropoutLinModel[j,],
                            vecMu=vecMuParam,
                            matPiConstPredictors=lsDropModel$matPiConstPredictors,
                            lsDropModelGlobal=lsDropModel$lsDropModelGlobal )
                        if(lsDispModelFull$lsDispModelGlobal$strDispModel=="constant"){
                            vecDispParam <- as.vector(lsDispModelFull$matDispModel)
                        } else {
                            stop(paste0("ERROR estimateMMAssignmentsMatrix(): strDispModel=", lsDispModel$lsDispModelGlobal$strDispModel, " not recognised."))
                        }
                        vecCounts <- objLP@matCountsProc[,j]
                        scaLogLik <- evalLikZINB_comp(
                            vecCounts=vecCounts,
                            vecMu=vecMuParam*lsMuModelFull$lsMuModelGlobal$vecNormConst[j],
                            vecDisp=vecDispParam, 
                            vecPi=vecDropParam,
                            vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
                            vecboolZero= !is.na(vecCounts) & vecCounts==0 )
                        return(scaLogLik)
                    })
                }))
                vecidxBadlyDescrCells <- sort(apply(matLLCellgivenMixture,1,max), decreasing=FALSE, index.return=TRUE)$ix
                idxCellForCentroidsReset <- vecidxBadlyDescrCells[round(length(vecidxBadlyDescrCells)/5)]
                vecCentroid <- objLP@matCountsProc[,idxCellForCentroidsReset]
                # Add pseudo count to not generate error in optim in log space
                # Do not generate disadvantage through pseudocount for modeling
                # of zeros for new centroid: Pseudocount in minimum modelled value.
                vecCentroid[vecCentroid==0] <- min(lsMuModelFull$matMuModel)
                lsMuModelFull$matMuModel[,which(vecboolMixtureDropped)[1]] <- vecCentroid
                # Compute new likelihood
                scaLogLikOld <- scaLogLikNew
                scaLogLikNew <- sum(evalLogLikMatrix(
                    matCounts=objLP@matCountsProc,
                    lsMuModel=lsMuModelFull,
                    lsDispModel=lsDispModelFull, 
                    lsDropModel=lsDropModel,
                    matWeights=matWeights ))
                
                strMessage <- paste0("# ", scaIter,".   E-Reset complete: ",
                                     "ll ", scaLogLikNew, 
                                     " (Mixture ", which(vecboolMixtureDropped)[1],").")
                objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                if(boolSuperVerbose) print(strMessage)
            } else {
                # M-like step: Estimate mixture model parameters
                tm_mstep <- system.time({
                    lsZINBFitsFull <- fitZINB(
                        matCounts=objLP@matCountsProc,
                        dfAnnotation=objLP@dfAnnotationProc,
                        vecConfoundersMu=objLP@vecConfoundersMu,
                        vecConfoundersDisp=objLP@vecConfoundersDisp,
                        vecNormConst=objLP@vecNormConst,
                        matWeights=matWeights,
                        matPiConstPredictors=NULL,
                        lsDropModel=lsDropModel,
                        matMuModelInit=lsMuModelFull$matMuModel,
                        lsmatBatchModelInitMu=lsMuModelFull$lsmatBatchModel,
                        matDispModelInit=lsDispModelFull$matDispModel,
                        lsmatBatchModelInitDisp=lsDispModelFull$lsmatBatchModel,
                        strMuModel="MM",
                        strDispModel=strDispModel,
                        scaMaxEstimationCycles=1,
                        boolVerbose=FALSE,
                        boolSuperVerbose=FALSE)
                    lsMuModelFull <- lsZINBFitsFull$lsMuModel
                    lsDispModelFull <- lsZINBFitsFull$lsDispModel
                    boolConvergenceModelFull <- lsZINBFitsFull$boolConvergenceModel
                    objLP@strReport <- paste0(objLP@strReport,
                                              lsZINBFitsFull$strReport)
                    rm(lsZINBFitsFull)
                })
                vecLLNew <- evalLogLikMatrix(
                    matCounts=objLP@matCountsProc,
                    lsMuModel=lsMuModelFull,
                    lsDispModel=lsDispModelFull, 
                    lsDropModel=lsDropModel,
                    matWeights=matWeights )
                # Check that no gene has a worse model fit than the null model,
                # re-initialise centroid coordinates of this gene to constant if this 
                # is the case. Under arbitrary mixture assignments, setting all
                # mixtures for this gene to the inferred H0 will yield the H0
                # loglikelihood and will therefore result in an increase in 
                # the overall loglikelihood. Therefore:
                # 1. Convergence is not broken.
                # 2. The model get s a chance to leave a bad local minimum.
                # Note: Also update dispersion to H0.
                # Note: Due to numerical LL thresholding, the null model evaluated
                # as a constant model for all cells and evaluated in a mixture model
                # in which every centroid is the constant model, may differ.
                # The null model on the mixture model LL is the value we have to
                # compare against to guarantee convergence.
                vecLLH0MM <- evalLogLikMatrix(
                    matCounts=objLP@matCountsProc,
                    lsMuModel=lsMuModelRed,
                    lsDispModel=lsDispModelRed, 
                    lsDropModel=lsDropModel,
                    matWeights=matWeights,
                    boolConstModelOnMMLL=TRUE)
                vecboolBadGeneModel <- vecLLNew < vecLLH0MM
                if(any(vecboolBadGeneModel)){
                    lsMuModelFull$matMuModel[vecboolBadGeneModel,] <- matrix(
                        lsMuModelRed$matMuModel[vecboolBadGeneModel,],
                        nrow=sum(vecboolBadGeneModel), 
                        ncol=dim(lsMuModelFull$matMuModel)[2],
                        byrow=FALSE)
                    if(!is.null(lsMuModelFull$lsMuModelGlobal$scaNConfounders)) {
                        for(conf in seq(1, lsMuModelFull$lsMuModelGlobal$scaNConfounders, by=1)) {
                            lsMuModelFull$lsmatBatchModel[[conf]][vecboolBadGeneModel,] <- 
                                lsMuModelRed$lsmatBatchModel[[conf]][vecboolBadGeneModel,]
                        }
                    }
                    lsDispModelFull$matDispModel[vecboolBadGeneModel,] <- matrix(
                        lsDispModelRed$matDispModel[vecboolBadGeneModel,],
                        nrow=sum(vecboolBadGeneModel), 
                        ncol=dim(lsDispModelFull$matDispModel)[2],
                        byrow=FALSE)
                    if(!is.null(lsDispModelFull$lsDispModelGlobal$scaNConfounders)) {
                        for(conf in seq(1, lsDispModelFull$lsDispModelGlobal$scaNConfounders, by=1)) {
                            lsDispModelFull$lsmatBatchModel[[conf]][vecboolBadGeneModel,] <- 
                                lsDispModelRed$lsmatBatchModel[[conf]][vecboolBadGeneModel,]
                        }
                    }
                    vecLLNew <- evalLogLikMatrix(
                        matCounts=objLP@matCountsProc,
                        lsMuModel=lsMuModelFull,
                        lsDispModel=lsDispModelFull, 
                        lsDropModel=lsDropModel,
                        matWeights=matWeights )
                    strMessage <- paste0("# ",scaIter,".   M-step resulted in ",
                                         sum(vecboolBadGeneModel), 
                                         " gene models which are worse than a constant model.")
                    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                    if(boolSuperVerbose) print(strMessage)
                }
                scaLogLikOld <- scaLogLikNew
                #vecLogLikIter <- lsZINBFitsFull$vecEMLogLikModel
                #scaLogLikNew <- vecLogLikIter[sum(!is.na(vecLogLikIter))]
                scaLogLikNew <- sum(vecLLNew)
                
                strMessage <- paste0("# ",scaIter,".   M-step complete: ",
                                     "ll  ", scaLogLikNew, " in ",
                                     round(tm_mstep["elapsed"]/60,2)," min.")
                objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                if(boolSuperVerbose) print(strMessage)
                
                if(boolConvergenceModelFull !=0){
                    strMessage <- paste0("Model estimation did not converge.")
                    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
                    if(boolSuperVerbose) print(strMessage)
                }
            }
            
            strMessage <- paste0("# ",scaIter,". iteration complete: ",
                                 "ll ", scaLogLikNew, " in ",
                                 round(tm_estep["elapsed"]/60,2)," min.")
            objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
            if(boolVerbose & !boolSuperVerbose) print(strMessage)
            
            vecEMLogLikModelFull[scaIter] <- scaLogLikNew
            scaIter <- scaIter+1
        }
    })
    
    strMessage <- paste0("### Finished fitting H1 mixture ZINB model ",
                         "in ", round(tm_cycle["elapsed"]/60,2)," min.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    # Name rows and columns of parameter matrices
    rownames(matWeights) <- colnames(objLP@matCountsProc)
    
    lsFitZINBReporters <- list(
        boolConvergenceH1=boolConvergenceModelFull,
        boolConvergenceConst=boolConvergenceModelRed,
        vecEMLogLikH1=vecEMLogLikModelFull,
        vecEMLogLikConst=vecEMLogLikModelRed )
    
    objLP@lsMuModelH1        <- lsMuModelFull
    objLP@lsDispModelH1      <- lsDispModelFull
    objLP@lsMuModelH0        <- lsMuModelRed
    objLP@lsDispModelH0      <- lsDispModelRed
    objLP@lsDropModel        <- lsDropModel
    objLP@matWeights         <- matWeights
    objLP@lsFitZINBReporters <- lsFitZINBReporters
    
    # Fit non-constant null mixture model if necessary
    # Only possible if there s a reference set of mixtures (RSA scenario)
    if(objLP@boolFixedPopulations){
        objLP <- fitH0MixtureZINBModel(
            objLP=objLP,
            vecidxMixturesH0=unlist(lsvecFixedCentrByPop[objLP@vecH0Pop]),
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    }
    
    return(objLP)
}