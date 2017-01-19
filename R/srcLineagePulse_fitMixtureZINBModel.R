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

fitMixtureZINBModel <- function(matCounts,
                                vecFixedAssignments=NULL,
                                vecNormConst,
                                scaNMixtures,
                                strDispModel="const",
                                scaMaxEstimationCyclesDropModel=20,
                                scaMaxEstimationCyclesEMlike=20){
  
  scaNGenes <- dim(matCounts)[1]
  scaNCells <- dim(matCounts)[2]
  
  ### (A) Pre-estimate drop-out model to speed up mixture model fitting
  if(boolVerbose) print(paste0("### a) Fit H0 constant ZINB model and set drop-out model."))
  
  tm_cycle <- system.time({
    lsZINBFitsRed <- fitZINB(matCounts=matCounts,
                             vecNormConst=vecNormConst,
                             vecPseudotime=NULL,
                             lsResultsClustering=NULL,
                             matWeights=matWeights,
                             matPiConstPredictors=NULL,
                             scaWindowRadius=NULL,
                             boolVecWindowsAsBFGS=NULL,
                             lsDropModel=NULL,
                             strMuModel="MM",
                             strDispModel="constant",
                             scaMaxEstimationCycles=scaMaxEstimationCyclesDropModel,
                             boolVerbose=TRUE,
                             boolSuperVerbose=TRUE)
    lsZINBFitsRed <- lsZINBFitsRed$lsMuModel
    lsZINBFitsRed <- lsZINBFitsRed$lsDispModel
    lsDropModel <- lsZINBFitsRed$lsDropModel
    boolConvergenceModelRed <- lsZINBFitsRed$boolConvergenceModel
    vecEMLogLikModelRed <- lsZINBFitsRed$vecEMLogLikModel
  })
  
  if(boolVerbose) print(paste0("### Finished fitting H0 constant ZINB model ",
                               "in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  ### (B) EM-like estimation cycle: fit mixture model
  if(boolVerbose) print(paste0("### b) EM-like iteration: Fit mixture assignments."))
  
  # (I) Initialise estimation: Weights
  # Set weights to uniform distribution
  matWeights <- matrix(1/scaNMixtures, 
                       nrow=scaNCells, ncols=scaNMixtures)
  # Correct fixed weights (RSA)
  if(!is.null(vecFixedAssignments)){
    matWeights[!is.na(vecFixedAssignments),] <- 0
    matWeights[!is.na(vecFixedAssignments),vecFixedAssignments[!is.na(vecFixedAssignments)]] <- 1
  }
  
  # (II) Estimation iteration on full model
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- scaLogLikInitA
  scaLogLikOld <- NA
  
  tm_RSAcycle <- system.time({
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCyclesEMlike)){
      # E-like step: Estimation of mixture assignments
      lsWeightFits <- estimateMMAssignmentsMatrix(matCounts=matCounts,
                                                  vecFixedAssignments=vecFixedAssignments,
                                                  lsMuModel=lsMuModelFull,
                                                  lsDispModel=lsDispModelFull,
                                                  lsDropModel=lsDropModelFull,
                                                  matWeights=matWeights )
      matWeights <- lsWeightFits$matWeights
      
      # M-like step: Estimate mixture model parameters
      lsZINBFitsFull <- fitZINB(matCounts=matCounts,
                                vecNormConst=vecNormConst,
                                vecPseudotime=NULL,
                                lsResultsClustering=NULL,
                                matWeights=matWeights,
                                matPiConstPredictors=NULL,
                                scaWindowRadius=NULL,
                                boolVecWindowsAsBFGS=NULL,
                                lsDropModel=lsDropModel,
                                strMuModel="MM",
                                strDispModel="constant",
                                scaMaxEstimationCycles=1,
                                boolVerbose=TRUE,
                                boolSuperVerbose=TRUE)
      lsMuModelFull <- lsZINBFitsFull$lsMuModel
      lsDispModelFull <- lsZINBFitsFull$lsDispModel
      lsDropModelFull <- lsZINBFitsFull$lsDropModel
    }
  })
  
  if(boolVerbose) print(paste0("### Finished fitting assignments ",
                               "in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  ### (C) Set degrees of freedom
  # Compute degrees of freedom of model for each gene
  # Drop-out model is ignored, would be counted at each gene.
  # 1. Alternative model  H1:
  # One mean parameter per mixture component
  scaKbyGeneH1 <- dim(matWeights)[2]
  # Dispersion model
  if(strDispModel=="constant"){
    # One dispersion factor per gene.
    scaKbyGeneH1 <- scaKbyGeneH1 + 1
  } else {
    stop(paste0("ERROR in fitMixtureZINBModel(): strMuModel not recognised: ", strDispModel))
  }
  # 2. Null model H0:
  # One mean per gene and one dispersion factor
  scaKbyGeneH0 <- 1+1
  
  # Name rows and columns of parameter matrices
  rownames(lsMuModelFull$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelFull$matDispModel) <- rownames(matCountsProc)
  rownames(lsMuModelRed$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelRed$matDispModel) <- rownames(matCountsProc)
  rownames(lsDropModel$matDropoutLinModel) <- colnames(matCountsProc)
  rownames(matWeights) <- colnames(matCountsProc)
  
  lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelFull,
                              boolConvergenceH0=boolConvergenceModelRed,
                              vecEMLogLikH1=vecEMLogLikModelFull,
                              vecEMLogLikH0=vecEMLogLikModelRed,
                              scaKbyGeneH1=scaKbyGeneH1,
                              scaKbyGeneH0=scaKbyGeneH0 )
  
  return( new('LineagePulseObject',
              dfResults           = NULL,
              matCounts           = matCountsProc,
              vecFixedAssignments = NULL,
              vecAllGenes         = NULL,
              lsMuModelH1         = lsMuModelFull,
              lsDispModelH1       = lsDispModelFull,
              lsMuModelH0         = lsMuModelRed,
              lsDispModelH0       = lsDispModelRed,
              lsDropModel         = lsDropModel,
              matWeights          = matWeights,
              lsFitZINBReporters  = lsFitZINBReporters,
              dfAnnotationProc    = NULL,
              vecNormConst        = vecNormConst,
              strReport           = NULL) )
}