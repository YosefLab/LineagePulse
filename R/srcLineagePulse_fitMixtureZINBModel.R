#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++  EM-like iteration for RSA  +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

fitMixtureZINBModel <- function(matCounts,
                                vecNormConst,
                                scaNMixtures,
                                scaMaxEstimationCycles=10){
  
  scaNGenes <- dim(matCounts)[1]
  scaNCells <- dim(matCounts)[2]
  # (I) Initialise estimation: Weights
  # Set weights to uniform distribution
  matWeights <- matrix(1/scaNMixtures, 
                       nrow=scaNGenes, ncols=scaNCells)
  
  # (II) Estimation iteration on full model
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- scaLogLikInitA
  scaLogLikOld <- NA
  
  tm_RSAcycle <- system.time({
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
      # M-like step: Estimate mixture model parameters
      lsZINBFitsWeights <- fitZINB(matCounts=matCounts,
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
                                   scaMaxEstimationCycles=20,
                                   boolVerbose=TRUE,
                                   boolSuperVerbose=TRUE)
      lsMuModelWeights <- lsZINBFitsWeights$lsMuModel
      lsDispModelWeights <- lsZINBFitsWeights$lsDispModel
      lsDropModelWeights <- lsZINBFitsWeights$lsDropModel
      
      # E-like step: Estimation of mixture assignments
      lsWeightFits <- estimateMMAssignmentsMatrix(matCounts=matCounts,
                                                  lsMuModel=lsMuModelWeights,
                                                  lsDispModel=lsDispModelWeights,
                                                  lsDropModel=lsDropModelWeights,
                                                  matWeights=matWeights )
      matWeights <- lsWeightFits$matWeights
    }
  })
  
  # Do full and alternative model fits
  ####################################################
  # Compute degrees of freedom of model for each gene
  # Drop-out model is ignored, would be counted at each gene.
  # 1. Alternative model  H1:
  # Mean model:
  if(strMuModel=="windows"){
    # One mean parameter per cell
    scaKbyGeneH1 <- dim(matCountsProc)[2]
  } else if(strMuModel=="clusters"){
    # One mean parameter per cluster
    scaKbyGeneH1 <- lsResultsClustering$K
  } else if(strMuModel=="MM"){
    # One mean parameter per mixture component
    scaKbyGeneH1 <- dim(matWeights)[2]
  } else if(strMuModel=="impulse"){
    # Six impulse model parameter to model means
    scaKbyGeneH1 <- 6
  } else if(strMuModel=="constant"){
    # One constant mean
    scaKbyGeneH1 <- 1
  } else {
    stop(paste0("ERROR in fitMixtureZINBModel(): strMuModel not recognised: ", strMuModel))
  }
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
  
  # Fit full model
  if(boolVerbose) print(paste0("### a) Fit H1 mixture ZINB model."))
  
  tm_cycle <- system.time({
    lsZINBFitsFull <- fitZINB(matCounts=matCounts,
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
                              scaMaxEstimationCycles=20,
                              boolVerbose=TRUE,
                              boolSuperVerbose=TRUE)
    lsMuModelFull <- lsZINBFitsFull$lsMuModel
    lsDispModelFull <- lsZINBFitsFull$lsDispModel
    lsDropModelFull <- lsZINBFitsFull$lsDropModel
    boolConvergenceModelFull <- lsZINBFitsFull$boolConvergenceModel
    vecEMLogLikModelFull <- lsZINBFitsFull$vecEMLogLikModel
  })
  
  if(boolVerbose) print(paste0("### Finished fitting H1 mixture ZINB model ",
                               "model in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  ####################################################
  # Fit model B
  if(boolVerbose) print(paste0("### b) Fit H0 constant ZINB model."))
  
  tm_cycleB <- system.time({
    lsFitsModelRed <- fitZINB(matCounts=matCounts,
                              vecNormConst=vecNormConst,
                              vecPseudotime=NULL,
                              lsResultsClustering=NULL,
                              matWeights=matWeights,
                              matPiConstPredictors=NULL,
                              scaWindowRadius=NULL,
                              boolVecWindowsAsBFGS=NULL,
                              lsDropModel=lsDropModelFull,
                              strMuModel="constant",
                              strDispModel="constant",
                              scaMaxEstimationCycles=20,
                              boolVerbose=TRUE,
                              boolSuperVerbose=TRUE)
    lsMuModelRed <- lsZINBFitsRed$lsMuModel
    lsDispModelRed <- lsZINBFitsRed$lsDispModel
    lsDropModelRed <- lsZINBFitsRed$lsDropModel
    boolConvergenceModelRed <- lsZINBFitsRed$boolConvergenceModel
    vecEMLogLikModelRed <- lsZINBFitsRed$vecEMLogLikModel
  })
  
  if(boolVerbose) print(paste0("### Finished fitting H0 constant ZINB model ",
                               "model in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  # Name rows and columns of parameter matrices
  rownames(lsMuModelFull$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelFull$matDispModel) <- rownames(matCountsProc)
  rownames(lsMuModelRed$matMuModel) <- rownames(matCountsProc)
  rownames(lsDispModelRed$matDispModel) <- rownames(matCountsProc)
  rownames(lsDropModel$matDropoutLinModel) <- colnames(matCountsProc)
  if(!is.null(lsDropModelFull$matPiConstPredictors)){
    rownames(lsDropModelFull$matPiConstPredictors) <- rownames(matCountsProc)
  }
  
  lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelFull,
                              boolConvergenceH0=boolConvergenceModelRed,
                              vecEMLogLikH1=vecEMLogLikModelFull,
                              vecEMLogLikH0=vecEMLogLikModelRed,
                              scaKbyGeneH1=scaKbyGeneH1,
                              scaKbyGeneH0=scaKbyGeneH0 )
  lsReturn <- list( lsMuModelH1=lsMuModelFull,
                    lsDispModelH1=lsDispModelFull,
                    lsMuModelH0=lsMuModelRed,
                    lsDispModelH0=lsDispModelRed,
                    lsDropModel=lsDropModelFull,
                    lsFitZINBReporters=lsFitZINBReporters )
  
  return(lsReturn)
}