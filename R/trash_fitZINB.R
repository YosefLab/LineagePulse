} else {
  # Estimate mean and dispersion parameters sequentially.
  # a) Negative binomial mean parameter
  tm_mu <- system.time({
    lsFitMu <- fitZINBMu( matCounts=matCounts,
                          vecNormConst=vecNormConst,
                          lsMuModel=lsMuModel,
                          lsDispModel=lsDispModel,
                          lsDropModel=lsDropModel,
                          scaWindowRadius=scaWindowRadius )
  })
  lsMuModel$matMuModel <- lsFitMu$matMuModel
  colnames(lsMuModel$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
  vecboolMuEstConverged <- lsFitMu$vecConvergence
  
  if(boolSuperVerbose){
    if(any(vecboolMuEstConverged != 0)){
      print(paste0("Mean estimation did not converge in ", 
                   sum(vecboolDispEstConverged), " cases."))
    }
    scaLogLikTemp <- evalLogLikMatrix( matCounts=matCounts,
                                       vecNormConst=vecNormConst,
                                       lsMuModel=lsMuModel,
                                       lsDispModel=lsDispModel, 
                                       lsDropModel=lsDropModel,
                                       scaWindowRadius=scaWindowRadius )
    print(paste0("# ",scaIter, ".2) Mean estimation complete: ",
                 "loglikelihood of         ", scaLogLikTemp, " in ",
                 round(tm_mu["elapsed"]/60,2)," min."))
  }
  
  # b) Negative binomial dispersion parameter
  # Use MLE of dispersion factor: numeric optimisation of likelihood.
  tm_phi <-system.time({
    lsFitDispModel <- fitZINBDisp( matCounts=matCounts,
                                   vecNormConst=vecNormConst,
                                   lsMuModel=lsMuModel,
                                   lsDispModel=lsDispModel,
                                   lsDropModel=lsDropModel,
                                   scaWindowRadius=scaWindowRadius )
  })
  lsDispModel$matDispModel <- lsFitDispModel$matDispModel
  colnames(lsDispModel$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
  vecboolDispEstConverged <- lsDispModel$vecConvergence
}

==============
  
####################################################
# Initialise model B

print("Initialise parameters model B")
# Initialise parameters:
# -> Mean model is re-fit from scratch
# -> Dispersion estimates from A are used for initialisation
# -> Drop-out model is kept from model A estimation and
#   point estimators adjusted to mu initialisation used.
lsMuModelB <- list( matMuModel=NA,
                    lsMuModelGlobal=list( strMuModel=strMuModelB,
                                          scaNumCells=scaNumCells,
                                          vecPseudotime=vecPseudotime,
                                          vecindClusterAssign=lsResultsClustering$Assignments,
                                          boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
                                          MAXIT_BFGS_Impulse=MAXIT_BFGS_Impulse,
                                          RELTOL_BFGS_Impulse=RELTOL_BFGS_Impulse) )
# B is alternative model H1
if(strMuModelB=="constant"){
  lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
} else if(strMuModelB=="impulse"){
  lsMuModelB$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=6)
  lsMuModelB$matMuModel[,2:4] <- log(matrix(vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE))
} else if(strMuModelB=="clusters"){
  lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
} else  if(strMuModelB=="windows"){
  lsMuModelB$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
} else {
  stop(paste0("ERROR fitZINB(): strMuModelB=", strMuModelB, " not recognised."))
}

lsDispModelB <- list( matDispModel=NA,
                      lsDispModelGlobal=list( strDispModel=strDispModelB,
                                              scaNumCells=scaNumCells,
                                              vecPseudotime=vecPseudotime,
                                              vecindClusterAssign=lsResultsClustering$Assignments) )
if(strDispModelB=="constant"){
  if(strDispModelA=="constant"){
    # Use values estimated for model A as initialisation
    lsDispModelB$matDispModel <-  lsDispModelA$matDispModel
  } else {
    # Initialise from scratch if different model used
    lsDispModelB$matDispModel <- matrix(0.001, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  }
} else {
  stop(paste0("ERROR fitZINB(): strDispModelB=", strDispModelB, " not recognised."))
}

####################################################
# Fit model B
print(paste0("### b) Fit negative binomial model B (",
             strNameModelB,")."))

tm_cycle <- system.time({
  lsFitsModelA <- fitZINB(matCounts=matCountsProc,
                          lsMuModel=lsMuModelB,
                          lsDispModel=lsDispModelB,
                          lsDropModel=lsDropModel,
                          boolFitDrop=FALSE,
                          scaWindowRadius=scaWindowRadius,
                          boolVerbose=boolVerbose,
                          boolSuperVerbose=boolSuperVerbose)
})
lsMuModelB <- lsFitsModelA$lsMuModel
lsDispModelB <- lsFitsModelA$lsDispModel
boolConvergenceModelB <- lsFisModelA$boolConvergenceModel

print(paste0("Finished fitting zero-inflated negative binomial ",
             "model A with noise model in ", round(tm_cycle["elapsed"]/60,2)," min."))
# Evaluate initialisation loglikelihood for model B
scaLogLikInitB <- evalLogLikMatrix(matCounts=matCountsProc,
                                   vecNormConst=vecNormConst,
                                   lsMuModel=lsMuModelB,
                                   lsDispModel=lsDispModelB, 
                                   lsDropModel=lsDropModel,
                                   scaWindowRadius=scaWindowRadius )
if(boolVerbose){
  print(paste0("Completed initialisation with ",
               "log likelihood of         ", scaLogLikInitB))
}

# Set iteration reporters
scaIter <- 1
scaLogLikNew <- scaLogLikInitB
scaLogLikOld <- NA

tm_cycle <- system.time({
  if(boolCoEstDispMean){
    tm_dispmu <- system.time({
      # Estimate mean and dispersion parameters simultaneously.
      # a/b) Negative binomial mean AND dispersion parameter.
      lsFitMuDisp <- fitZINBMuDisp(matCountsProc=matCountsProc,
                                   vecNormConst=vecNormConst,
                                   lsMuModel=lsMuModelB,
                                   lsDispModel=lsDispModelB,
                                   lsDropModel=lsDropModel,
                                   scaWindowRadius=scaWindowRadius)
      lsDispModelB$matDispModel <- lsFitMuDisp$matDispModel
      colnames(lsDispModelB$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
      lsMuModelB$matMuModel <- lsFitMuDisp$matMuModel
      colnames(lsMuModelB$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
      
      vecboolMuEstConvergedB <- lsFitMuDisp$vecConvergence
      vecboolDispEstConvergedB <- lsFitMuDisp$vecConvergence
      
      # Evaluate Likelihood
      scaLogLik <- evalLogLikMatrix( matCounts=matCountsProc,
                                     vecNormConst=vecNormConst,
                                     lsMuModel=lsMuModelB,
                                     lsDispModel=lsDispModelB, 
                                     lsDropModel=lsDropModel,
                                     scaWindowRadius=scaWindowRadius )
    })
    
    # Iteration complete
    if(boolSuperVerbose){
      if(any(vecboolDispEstConvergedB != 0)){
        print(paste0("Mean-Dispersion estimation did not converge in ", 
                     sum(vecboolMuEstConvergedB), " cases [codes: ",
                     paste(unique(vecboolMuEstConvergedB[vecboolMuEstConvergedB!=0])), "]."))
      }
    } 
    if(boolVerbose){print(paste0("Mean+Disp estimation complete: ",
                                 "log likelihood of        ", scaLogLik, " in ",
                                 round(tm_dispmu["elapsed"]/60,2)," min."))
    }
  } else {
    # Estimate mean and dispersion parameters sequentially.
    while(scaIter == 1 |
          (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
      tm_iter <- system.time({
        ##### Gene-wise parameter estimation: 
        # a) Negative binomial mean parameter
        # Only compute posterior if using closed form estimator for mean:
        # Posterior is not necessary in all other cases, expect for parameter
        # initialisation of impulse model.
        tm_mu <- system.time({
          lsFitMuModelB <- fitZINBMu( matCountsProc=matCountsProc,
                                      vecNormConst=vecNormConst,
                                      lsMuModel=lsMuModelB,
                                      lsDispModel=lsDispModelB,
                                      lsDropModel=lsDropModel,
                                      scaWindowRadius=scaWindowRadius )
        })
        lsMuModelB$matMuModel <- lsFitMuModelB$matMuModel
        colnames(lsMuModelB$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
        vecboolMuEstConvergedB <- lsFitMuModelB$vecConvergence
        
        if(boolSuperVerbose){
          if(any(vecboolMuEstConvergedB != 0)){
            print(paste0("Mean estimation did not converge in ", 
                         sum(vecboolMuEstConvergedB), " cases."))
          }
          scaLogLikTemp <- evalLogLikMatrix( matCounts=matCountsProc,
                                             vecNormConst=vecNormConst,
                                             lsMuModel=lsMuModelB,
                                             lsDispModel=lsDispModelB, 
                                             lsDropModel=lsDropModel,
                                             scaWindowRadius=scaWindowRadius )
          print(paste0("# ",scaIter, ".1) Mean estimation complete: ",
                       "loglikelihood of       ", scaLogLikTemp, " in ",
                       round(tm_mu["elapsed"]/60,2)," min."))
        }
        
        # b) Negative binomial dispersion parameter
        # Use MLE of dispersion factor: numeric optimisation of likelihood.
        tm_phi <- system.time({
          lsFitDispModelB <- fitZINBDisp( matCountsProc=matCountsProc,
                                          vecNormConst=vecNormConst,
                                          lsMuModel=lsMuModelB,
                                          lsDispModel=lsDispModelB,
                                          lsDropModel=lsDropModel,
                                          scaWindowRadius=scaWindowRadius )
        })
        lsDispModelB$matDispModel <- lsFitDispModelB$matDispModel
        colnames(lsDispModelB$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
        vecboolDispEstConvergedB <- lsDispModelB$vecConvergence
        
        # Evaluate Likelihood
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- evalLogLikMatrix( matCounts=matCountsProc,
                                          vecNormConst=vecNormConst,
                                          lsMuModel=lsMuModelB,
                                          lsDispModel=lsDispModelB, 
                                          lsDropModel=lsDropModel,
                                          scaWindowRadius=scaWindowRadius )
      })
      
      # Iteration complete
      if(boolSuperVerbose){
        if(any(vecboolDispEstConvergedB != 0)){
          print(paste0("Dispersion estimation did not converge in ", 
                       sum(vecboolDispEstConvergedB), " cases [codes: ",
                       paste(unique(vecboolDispEstConvergedB[vecboolDispEstConvergedB!=0])), "]."))
        }
        print(paste0("# ",scaIter, ".2) Dispersion estimation complete: ",
                     "loglikelihood of ", scaLogLikNew, " in ",
                     round(tm_phi["elapsed"]/60,2)," min."))
      } else {
        if(boolVerbose){print(paste0("# ",scaIter, ".) complete with ",
                                     "log likelihood of ", scaLogLikNew, " in ",
                                     round(tm_iter["elapsed"]/60,2)," min."))}
      }
      vecEMLogLikModelB[scaIter] <- scaLogLikNew
      scaIter <- scaIter+1
    }
  }
})
print(paste0("Finished fitting zero-inflated negative ",
             "binomial model B in ", round(tm_cycle["elapsed"]/60,2)," min."))

# Evaluate convergence
if(all(as.logical(vecboolDispEstConvergedB)) &
   all(as.logical(vecboolMuEstConvergedB)) &
   scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
  boolConvergenceModelB <- TRUE
} else { boolConvergenceModelB <- FALSE }
