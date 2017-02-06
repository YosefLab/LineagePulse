#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  EM-like iteration for mixture model  ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit ZINB mixture model
#' 
#' Structure of code:
#' A: Fit a constant model to all cells. The resulting drop-out model
#'    is used for the EM-like iteration. This is also the null model
#'    for the LRT as the full model is fir based on this drop-out model.
#' B: EM-like iteration to get mixture assignments.
#' C: Fit full model based on mixture assignments as MLE to 
#'    do LRT later.
#'    
#' @export
fitMixtureZINBModel <- function(objectLineagePulse,
                                scaNMixtures,
                                strDispModel="const",
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
  
  scaNGenes <- dim(objectLineagePulse@matCountsProc)[1]
  scaNCells <- dim(objectLineagePulse@matCountsProc)[2]
  
  ### (A) Pre-estimate drop-out model to speed up mixture model fitting
  strMessage <- paste0("### a) Fit H0 constant ZINB model and set drop-out model.")
  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  tm_cycle <- system.time({
    lsZINBFitsRed <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                             dfAnnotation=objectLineagePulse@dfAnnotationProc,
                             vecConfounders=objectLineagePulse@vecConfounders,
                             vecNormConst=objectLineagePulse@vecNormConst,
                             matWeights=NULL,
                             matPiConstPredictors=NULL,
                             scaWindowRadius=NULL,
                             boolVecWindowsAsBFGS=FALSE,
                             lsDropModel=NULL,
                             strMuModel="constant",
                             strDispModel="constant",
                             scaMaxEstimationCycles=scaMaxEstimationCyclesDropModel,
                             boolVerbose=boolVerbose,
                             boolSuperVerbose=boolSuperVerbose)
    lsMuModelRed <- lsZINBFitsRed$lsMuModel
    lsDispModelRed <- lsZINBFitsRed$lsDispModel
    lsDropModel <- lsZINBFitsRed$lsDropModel
    boolConvergenceModelRed <- lsZINBFitsRed$boolConvergenceModel
    vecEMLogLikModelRed <- lsZINBFitsRed$vecEMLogLikModel
    objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport,
                                           lsZINBFitsRed$strReport)
  })
  
  strMessage <- paste0("### Finished fitting H0 constant ZINB model ",
                               "in ", round(tm_cycle["elapsed"]/60,2)," min.")
  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  ### (B) EM-like estimation cycle: fit mixture model
  strMessage <- paste0("### b) EM-like iteration: Fit mixture model")
  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  # (I) Initialise estimation
  # Initialise expression models as null model estimate:
  # Use cells to initialise centroids: Chose cells so that 
  # the overall distance is maximised.
  lsMuModelFull <- lsMuModelRed
  lsMuModelFull$lsMuModelGlobal$strMuModel <- "MM"
  vecidxCellsForCentroids <- initialiseCentroidsFromCells(matCounts=objectLineagePulse@matCountsProc,
                                                          vecFixedCells=objectLineagePulse@dfAnnotationProc$populations,
                                                          scaN=scaNMixtures)
  # Add pseudo count to not generate error in optim in log space
  lsMuModelFull$matMuModel <- objectLineagePulse@matCountsProc[,vecidxCellsForCentroids]
  lsMuModelFull$matMuModel[lsMuModelFull$matMuModel==0] <- 10^(-5)
  lsDispModelFull <- lsDispModelRed
  # Set weights to uniform distribution
  matWeights <- matrix(1/scaNMixtures, 
                       nrow=scaNCells, ncol=scaNMixtures)
  # Correct fixed weights (RSA)
  if(objectLineagePulse@boolFixedPopulations){
    vecFixedAssignments <- objectLineagePulse@dfAnnotationProc$populations
    matWeights[!is.na(vecFixedAssignments),] <- 0
    matWeights[!is.na(vecFixedAssignments),vecFixedAssignments[!is.na(vecFixedAssignments)]] <- 1
  }
  
  # (II) Estimation iteration on full model
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- sum(evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                       lsMuModel=lsMuModelFull,
                                       lsDispModel=lsDispModelFull, 
                                       lsDropModel=lsDropModel,
                                       matWeights=matWeights,
                                       scaWindowRadius=NULL ))
  scaLogLikOld <- NA
  vecEMLogLikModelFull <- array(NA, scaMaxEstimationCyclesEMlike)
  
  strMessage <- paste0("#  .   Initialised MM: ",
                       "loglikelihood of   ", scaLogLikNew)
  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolSuperVerbose) print(strMessage)
  
  boolResetDropout <- FALSE
  tm_RSAcycle <- system.time({
    while(scaIter == 1 | 
          boolResetDropout| 
          (scaLogLikNew > scaLogLikOld*scaPrecEMAssignments & scaIter <= scaMaxEstimationCyclesEMlike)){
      boolResetDropout <- FALSE
      # E-like step: Estimation of mixture assignments
      tm_estep <- system.time({
        lsWeightFits <- estimateMMAssignmentsMatrix(matCounts=objectLineagePulse@matCountsProc,
                                                    dfAnnotation=objectLineagePulse@dfAnnotationProc,
                                                    boolFixedPopulations=objectLineagePulse@boolFixedPopulations,
                                                    vecNormConst=objectLineagePulse@vecNormConst,
                                                    lsMuModel=lsMuModelFull,
                                                    lsDispModel=lsDispModelFull,
                                                    lsDropModel=lsDropModel,
                                                    matWeights=matWeights,
                                                    MAXIT_BFGS_MM=MAXIT_BFGS_MM,
                                                    RELTOL_BFGS_MM=RELTOL_BFGS_MM )
        matWeights <- lsWeightFits$matWeights
      })
      
      scaLogLikTemp <- sum(evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                            lsMuModel=lsMuModelFull,
                                            lsDispModel=lsDispModelFull, 
                                            lsDropModel=lsDropModel,
                                            matWeights=matWeights,
                                            scaWindowRadius=NULL ))
      
      strMessage <- paste0("# ", scaIter,".   E-step complete: ",
                   "loglikelihood of  ", scaLogLikTemp, " in ",
                   round(tm_estep["elapsed"]/60,2)," min.")
      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
      if(boolSuperVerbose) print(strMessage)
      
      if(any(lsWeightFits$vecConvergence !=0 )){
        strMessage <- paste0("Weight estimation did not convergen in ",
                             sum(lsWeightFits$vecConvergence !=0), " cases.")
        objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
        if(boolSuperVerbose) print(strMessage)
      }
      
      # Catch mixture drop-out
      vecidxAssignedMixture <- apply(matWeights, 1, which.max)
      vecboolMixtureDropped <- sapply(seq(1,scaNMixtures), function(mixture){
        sum(vecidxAssignedMixture==mixture) <= 1
      })
      if(scaIter >1 & any(vecboolMixtureDropped)){
        # Reinitialise one mixture component
        # Get badly described cell:
        # Don't chose minimum to not catch outlier cells.
        vecidxMaxWeightsSort <- sort(apply(matWeights, 1, max), index.return=TRUE, decreasing=FALSE)$ix
        idxCellForCentroidsReset <- vecidxMaxWeightsSort[round(length(vecidxMaxWeightsSort)/5)]
        vecCentroid <- objectLineagePulse@matCountsProc[,idxCellForCentroidsReset]
        # Add pseudo count to not generate error in optim in log space
        # Do not generate disadvantage through pseudocount for modeling
        # of zeros for new centroid: Pseudocount in minimum modelled value.
        vecCentroid[vecCentroid==0] <- min(lsMuModelFull$matMuModel)
        lsMuModelFull$matMuModel[,which(vecboolMixtureDropped)[1]] <- vecCentroid
        # Compute new likelihood
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- sum(evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                             lsMuModel=lsMuModelFull,
                                             lsDispModel=lsDispModelFull, 
                                             lsDropModel=lsDropModel,
                                             matWeights=matWeights,
                                             scaWindowRadius=NULL ))
        
        strMessage <- paste0("# ", scaIter,".   E-Reset complete: ",
                             "loglikelihood of ", scaLogLikNew, 
                             " (Mixture ", which(vecboolMixtureDropped)[1],").")
        objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
        if(boolSuperVerbose) print(strMessage)
      }
      
      # M-like step: Estimate mixture model parameters
      tm_mstep <- system.time({
        lsZINBFitsFull <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                                  dfAnnotation=objectLineagePulse@dfAnnotationProc,
                                  vecConfounders=objectLineagePulse@vecConfounders,
                                  vecNormConst=objectLineagePulse@vecNormConst,
                                  matWeights=matWeights,
                                  matPiConstPredictors=NULL,
                                  scaWindowRadius=NULL,
                                  boolVecWindowsAsBFGS=FALSE,
                                  lsDropModel=lsDropModel,
                                  matMuModelInit=lsMuModelFull$matMuModel,
                                  lsmatBatchModelInit=lsMuModelFull$lsmatBatchModel,
                                  matDispModelInit=lsDispModelFull$matDispModel,
                                  strMuModel="MM",
                                  strDispModel="constant",
                                  scaMaxEstimationCycles=1,
                                  boolVerbose=TRUE,
                                  boolSuperVerbose=TRUE)
        lsMuModelFull <- lsZINBFitsFull$lsMuModel
        lsDispModelFull <- lsZINBFitsFull$lsDispModel
        boolConvergenceModelFull <- lsZINBFitsFull$boolConvergenceModel
        objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport,
                                               lsZINBFitsFull$strReport)
      })
      scaLogLikOld <- scaLogLikNew
      vecLogLikIter <- lsZINBFitsFull$vecEMLogLikModel
      scaLogLikNew <- vecLogLikIter[sum(!is.na(vecLogLikIter))]
      
      strMessage <- paste0("# ",scaIter,".   M-step complete: ",
                           "loglikelihood of  ", scaLogLikNew, " in ",
                           round(tm_mstep["elapsed"]/60,2)," min.")
      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
      if(boolSuperVerbose) print(strMessage)
      
      if(lsZINBFitsFull$boolConvergenceModel !=0){
        strMessage <- paste0("Model estimation did not converge.")
        objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
        if(boolSuperVerbose) print(strMessage)
      }
      
      strMessage <- paste0("# ",scaIter,". iteration complete: ",
                           "loglikelihood of ", scaLogLikNew, " in ",
                           round(tm_estep["elapsed"]/60,2)," min.")
      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
      if(boolVerbose & !boolSuperVerbose) print(strMessage)
      
      vecEMLogLikModelFull[scaIter] <- scaLogLikNew
      scaIter <- scaIter+1
    }
  })
  
  strMessage <- paste0("### Finished fitting H1 mixture ZINB model ",
                       "in ", round(tm_cycle["elapsed"]/60,2)," min.")
  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  # Name rows and columns of parameter matrices
  rownames(matWeights) <- colnames(objectLineagePulse@matCountsProc)
  
  lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelFull,
                              boolConvergenceH0=boolConvergenceModelRed,
                              vecEMLogLikH1=vecEMLogLikModelFull,
                              vecEMLogLikH0=vecEMLogLikModelRed )
  
  objectLineagePulse@lsMuModelH1        <- lsMuModelFull
  objectLineagePulse@lsDispModelH1      <- lsDispModelFull
  objectLineagePulse@lsMuModelH0        <- lsMuModelRed
  objectLineagePulse@lsDispModelH0      <- lsDispModelRed
  objectLineagePulse@lsDropModel        <- lsDropModel
  objectLineagePulse@matWeights         <- matWeights
  objectLineagePulse@lsFitZINBReporters <- lsFitZINBReporters
  
  return(objectLineagePulse)
}