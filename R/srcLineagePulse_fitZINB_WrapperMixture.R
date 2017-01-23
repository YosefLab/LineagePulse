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
  if(boolVerbose) print(paste0("### a) Fit H0 constant ZINB model and set drop-out model."))
  
  tm_cycle <- system.time({
    lsZINBFitsRed <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                             vecNormConst=objectLineagePulse@vecNormConst,
                             vecPseudotime=NULL,
                             lsResultsClustering=NULL,
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
  })
  
  if(boolVerbose) print(paste0("### Finished fitting H0 constant ZINB model ",
                               "in ", round(tm_cycle["elapsed"]/60,2)," min."))
  
  ### (B) EM-like estimation cycle: fit mixture model
  if(boolVerbose) print(paste0("### b) EM-like iteration: Fit mixture model"))
  
  # (I) Initialise estimation
  # Initialise expression models as null model estimate
  lsMuModelFull <- lsMuModelRed
  lsDispModelFull <- lsDispModelRed
  # Set weights to uniform distribution
  matWeights <- matrix(1/scaNMixtures, 
                       nrow=scaNCells, ncol=scaNMixtures)
  # Correct fixed weights (RSA)
  if(!is.null(objectLineagePulse@vecFixedAssignments)){
    matWeights[!is.na(objectLineagePulse@vecFixedAssignments),] <- 0
    matWeights[!is.na(objectLineagePulse@vecFixedAssignments),objectLineagePulse@vecFixedAssignments[!is.na(objectLineagePulse@vecFixedAssignments)]] <- 1
  }
  # Set expression models to NULL
  lsMuModelFull <- NULL
  lsDispModelFull <- NULL
  
  # (II) Estimation iteration on full model
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- -Inf
  scaLogLikOld <- NA
  vecEMLogLikModelFull <- array(NA, scaMaxEstimationCyclesEMlike)
  
  tm_RSAcycle <- system.time({
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEMAssignments & scaIter <= scaMaxEstimationCyclesEMlike)){
      # M-like step: Estimate mixture model parameters
      tm_mstep <- system.time({
        lsZINBFitsFull <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                                  vecNormConst=objectLineagePulse@vecNormConst,
                                  vecPseudotime=NULL,
                                  lsResultsClustering=NULL,
                                  matWeights=matWeights,
                                  matPiConstPredictors=NULL,
                                  scaWindowRadius=NULL,
                                  boolVecWindowsAsBFGS=FALSE,
                                  lsDropModel=lsDropModel,
                                  matMuModelInit=lsMuModelFull$matMuModel,
                                  matDispModelInit=lsDispModelFull$matDispModel,
                                  strMuModel="MM",
                                  strDispModel="constant",
                                  scaMaxEstimationCycles=1,
                                  boolVerbose=FALSE,
                                  boolSuperVerbose=FALSE)
        lsMuModelFull <- lsZINBFitsFull$lsMuModel
        lsDispModelFull <- lsZINBFitsFull$lsDispModel
      })
      if(boolSuperVerbose){
        vecLogLikIter <- lsZINBFitsFull$vecEMLogLikModel
        scaLogLikTemp <- vecLogLikIter[sum(!is.na(vecLogLikIter))]
        print(paste0("# ",scaIter,".   M-step complete: ",
                     "loglikelihood of ", scaLogLikTemp, " in ",
                     round(tm_mstep["elapsed"]/60,2)," min."))
        if(lsZINBFitsFull$boolConvergenceModel !=0) print(paste0("Model estimation did not converge."))
      }
      
      # E-like step: Estimation of mixture assignments
      tm_estep <- system.time({
        lsWeightFits <- estimateMMAssignmentsMatrix(matCounts=objectLineagePulse@matCountsProc,
                                                    vecFixedAssignments=objectLineagePulse@vecFixedAssignments,
                                                    lsMuModel=lsMuModelFull,
                                                    lsDispModel=lsDispModelFull,
                                                    lsDropModel=lsDropModel,
                                                    vecNormConst=objectLineagePulse@vecNormConst,
                                                    matWeights=matWeights,
                                                    MAXIT_BFGS_MM=MAXIT_BFGS_MM,
                                                    RELTOL_BFGS_MM=RELTOL_BFGS_MM )
        matWeights <- lsWeightFits$matWeights
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- sum(lsWeightFits$vecLL, na.rm=TRUE)
      })
      if(boolSuperVerbose){
        print(paste0("# ",scaIter,".   E-step complete: ",
                     "loglikelihood of ", scaLogLikNew, " in ",
                     round(tm_estep["elapsed"]/60,2)," min."))
        if(any(lsWeightFits$vecConvergence !=0 )) print(paste0("Weight estimation did not convergen in ",
                                                               sum(lsWeightFits$vecConvergence !=0), " cases."))
      } else if(boolVerbose){
        print(paste0("# ",scaIter,". iteration complete: ",
                     "loglikelihood of ", scaLogLikNew, " in ",
                     round(tm_estep["elapsed"]/60,2)," min."))
      }
      
      vecEMLogLikModelFull[scaIter] <- scaLogLikNew
      scaIter <- scaIter+1
    }
  })
  # Final model estimation to get MLE
  tm_finalmstep <- system.time({
    lsZINBFitsFull <- fitZINB(matCounts=objectLineagePulse@matCountsProc,
                              vecNormConst=objectLineagePulse@vecNormConst,
                              vecPseudotime=NULL,
                              lsResultsClustering=NULL,
                              matWeights=matWeights,
                              matPiConstPredictors=NULL,
                              scaWindowRadius=NULL,
                              boolVecWindowsAsBFGS=FALSE,
                              lsDropModel=lsDropModel,
                              matMuModelInit=lsMuModelFull$matMuModel,
                              matDispModelInit=lsDispModelFull$matDispModel,
                              strMuModel="MM",
                              strDispModel="constant",
                              scaMaxEstimationCycles=1,
                              boolVerbose=FALSE,
                              boolSuperVerbose=FALSE)
    lsMuModelFull <- lsZINBFitsFull$lsMuModel
    lsDispModelFull <- lsZINBFitsFull$lsDispModel
    boolConvergenceModelFull <- lsZINBFitsFull$boolConvergenceModel
  })
  if(boolSuperVerbose){
    vecLogLikIter <- lsZINBFitsFull$vecEMLogLikModel
    scaLogLikTemp <- vecLogLikIter[length(vecLogLikIter)]
    print(paste0("# Final M-step complete: ",
                 "loglikelihood of ", scaLogLikTemp, " in ",
                 round(tm_finalmstep["elapsed"]/60,2)," min."))
    if(lsZINBFitsFull$boolConvergenceModel !=0) print(paste0("Model estimation did not converge."))
  }
  
  if(boolVerbose) print(paste0("### Finished fitting H1 mixture ZINB model ",
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
  rownames(lsMuModelFull$matMuModel) <- rownames(objectLineagePulse@matCountsProc)
  rownames(lsDispModelFull$matDispModel) <- rownames(objectLineagePulse@matCountsProc)
  rownames(lsMuModelRed$matMuModel) <- rownames(objectLineagePulse@matCountsProc)
  rownames(lsDispModelRed$matDispModel) <- rownames(objectLineagePulse@matCountsProc)
  rownames(lsDropModel$matDropoutLinModel) <- colnames(objectLineagePulse@matCountsProc)
  rownames(matWeights) <- colnames(objectLineagePulse@matCountsProc)
  
  lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceModelFull,
                              boolConvergenceH0=boolConvergenceModelRed,
                              vecEMLogLikH1=vecEMLogLikModelFull,
                              vecEMLogLikH0=vecEMLogLikModelRed,
                              scaKbyGeneH1=scaKbyGeneH1,
                              scaKbyGeneH0=scaKbyGeneH0 )
  
  objectLineagePulse@lsMuModelH1        <- lsMuModelFull
  objectLineagePulse@lsDispModelH1      <- lsDispModelFull
  objectLineagePulse@lsMuModelH0        <- lsMuModelRed
  objectLineagePulse@lsDispModelH0      <- lsDispModelRed
  objectLineagePulse@lsDropModel        <- lsDropModel
  objectLineagePulse@matWeights         <- matWeights
  objectLineagePulse@lsFitZINBReporters <- lsFitZINBReporters
  objectLineagePulse@strReport          <- NULL
  
  return(objectLineagePulse)
}