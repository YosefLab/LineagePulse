fitZINB <- function(matCounts,
                    vecNormConst,
                    vecPseudotime=NULL,
                    lsResultsClustering=NULL,
                    matWeights=NULL,
                    matPiConstPredictors=NULL,
                    scaWindowRadius=NULL,
                    lsDropModel,
                    strMuModel,
                    strDispModel,
                    scaMaxEstimationCycles=20,
                    boolVerbose=TRUE,
                    boolSuperVerbose=TRUE){
  
  ####################################################
  # Internal Numerical Estimation Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-4)
  # Numerical optmisation of impulse model hyperparameters
  MAXIT_BFGS_Impulse <- 1000 # optim default is 1000
  RELTOL_BFGS_Impulse <- 10^(-4) # optim default is sqrt(.Machine$double.eps)=1e-8
  # Lowering RELTOL_BFGS_IMPULSE gives drastic run time improvements.
  # Set to 10^(-4) to maintain sensible fits without running far into saturation
  # in the objective (loglikelihood).
  # Numerical optmisation of dropout model hyperparameters
  MAXIT_BFGS_Pi <- 1000 
  RELTOL_BFGS_Pi <- 10^(-4)
  ####################################################
  
  vecEMLogLikModel <- array(NA, scaMaxEstimationCycles)
  scaNumGenes <- dim(matCounts)[1]
  scaNumCells <- dim(matCounts)[2]  
  boolFitDrop <- is.null(lsDropModel)
  
  # Initialise
  if(boolVerbose) print("Initialise parameter models.")
  
  # a) Mu model
  # Mean parameters (mu): Gene-wise mean of non-zero observations.
  # Impulse model: Initialised to constant (mean).
  lsMuModel <- list( matMuModel=NA,
                     lsMuModelGlobal=list( strMuModel=strMuModel,
                                           scaNumCells=scaNumCells,
                                           vecPseudotime=vecPseudotime,
                                           vecindClusterAssign=lsResultsClustering$Assignments,
                                           boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
                                           MAXIT_BFGS_Impulse=MAXIT_BFGS_Impulse,
                                           RELTOL_BFGS_Impulse=RELTOL_BFGS_Impulse) )
  vecMuModelInit <- apply(matCounts, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
  vecMuModelInit[vecMuModelInit < .Machine$double.eps] <- .Machine$double.eps
  if(strMuModel=="constant"){
    lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  } else if(strMuModel=="impulse"){
    lsMuModel$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=6)
    lsMuModel$matMuModel[,2:4] <- log(matrix(vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE))
  } else if(strMuModel=="clusters"){
    lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
  } else if(strMuModel=="MM"){
    lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=dim(matWeights)[2], byrow=FALSE)
  } else  if(strMuModel=="windows"){
    lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  } else {
    stop(paste0("ERROR fitZINB(): strMuModel=", strMuModel, " not recognised."))
  }
  
  # b) Dispersion model
  # Dispersions: Low dispersion factor yielding high variance which makes
  # cost function screening easy in the first iteration.
  lsDispModel <- list( matDispModel=NA,
                       lsDispModelGlobal=list( strDispModel=strDispModel,
                                               scaNumCells=scaNumCells,
                                               vecPseudotime=vecPseudotime,
                                               vecindClusterAssign=lsResultsClustering$Assignments) )
  if(strDispModel=="constant"){
    lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=1, byrow=FALSE)
  } else {
    stop(paste0("ERROR fitZINB(): strDispModel=", strDispModel, " not recognised."))
  }
  
  # c) Drop-out model: Only if this is ot given.
  # Dropout model: Initialise as offset=0 and log(mu)  parameter which
  # is forced to be negative during fitting, as -1. The parameter corresponding
  # to log(mu) may not be initialised too close to zero, as the cost function 
  # cannot always pick up the signal in such cases, leading to an MLE with this 
  # parameter untouched.
  if(is.null(lsDropModel)){
    lsDropModel <- list(matDropoutLinModel=NA,
                        matPiConstPredictors=matPiConstPredictors,
                        lsPiOptimHyperparam=list(
                          MAXIT_BFGS_Pi=MAXIT_BFGS_Pi,
                          RELTOL_BFGS_Pi=RELTOL_BFGS_Pi))
    # Target initialisation drop-out rate: 0.99, linear model mu
    # parameter = -1 -> solve for offset of linear model:
    scaPiTarget <- 0.99
    scaPiLinModelMuParam <- -1
    scaPiLinModelOffset <- log(scaPiTarget) - log(1-scaPiTarget) - 
      scaPiLinModelMuParam*log(min(vecMuModelInit, na.rm=TRUE))
    scaPredictors <- 2
    if(!is.null(matPiConstPredictors)){
      scaPredictors <- scaPredictors + dim(matPiConstPredictors)[2]
    }
    lsDropModel$matDropoutLinModel <- cbind(
      rep(scaPiLinModelOffset, scaNumCells), 
      rep(scaPiLinModelMuParam, scaNumCells),
      matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
  }
  
  # Evaluate initialisation loglikelihood
  scaLogLikInitA <- evalLogLikMatrix(matCounts=matCounts,
                                     vecNormConst=vecNormConst,
                                     lsMuModel=lsMuModel,
                                     lsDispModel=lsDispModel, 
                                     lsDropModel=lsDropModel,
                                     scaWindowRadius=scaWindowRadius )
  if(boolVerbose){
    print(paste0("Initialisation has  ",
                 "log likelihood of           ", scaLogLikInitA))
  }
  
  # Set iteration reporters
  scaIter <- 1
  scaLogLikNew <- scaLogLikInitA
  scaLogLikOld <- NA
  while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
    tm_iter <- system.time({
      if(boolFitDrop){
        #####  1. Cell-wise parameter estimation
        # Dropout rate
        # Drop-out estimation is independent between cells and can be parallelised.
        tm_pi <- system.time({
          lsFitsPi <- bplapply(seq(1, scaNumCells), function(cell){
            if(!is.null(scaWindowRadius)){
              scaindIntervalStart <- max(1,cell-scaWindowRadius)
              scaindIntervalEnd <- min(scaNumCells,cell+scaWindowRadius)
              vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
            } else {
              vecInterval <- cell
            }
            
            lsFitPi <- fitPiZINB(
              vecCounts=matCounts[,cell],
              vecDropoutLinModel=lsDropModel$matDropoutLinModel[cell,],
              matPiConstPredictors=lsDropModel$matPiConstPredictors,
              lsPiOptimHyperparam=lsDropModel$lsPiOptimHyperparam,
              lsMuModel=lsMuModel,
              lsDispModel=lsDispModel,
              scaNormConst=vecNormConst[cell],
              vecInterval=vecInterval,
              scaTarget=match(cell, vecInterval))
            return(lsFitPi)
          })
          lsDropModel$matDropoutLinModel <-  do.call(rbind, lapply(lsFitsPi, function(cell) cell$vecLinModel))
          vecboolPiEstConverged <- sapply(lsFitsPi, function(cell) cell$scaConvergence)  
        })
        colnames(lsDropModel$matDropoutLinModel) <- NULL # Want this so that column names dont grow to par.par.par...
        if(boolSuperVerbose){
          scaLogLikTemp <- evalLogLikMatrix( matCounts=matCounts,
                                             vecNormConst=vecNormConst,
                                             lsMuModel=lsMuModel,
                                             lsDispModel=lsDispModel, 
                                             lsDropModel=lsDropModel,
                                             scaWindowRadius=scaWindowRadius )
          if(any(vecboolPiEstConverged != 0)){
            print(paste0("Dropout estimation did not converge in ", 
                         sum(vecboolPiEstConverged), " cases [codes: ",
                         paste(unique(vecboolPiEstConverged[vecboolPiEstConverged!=0])), "]."))
          }
          if(any(vecboolPiEstConverged==1001)){
            print(paste0("Fatal dropout estimation error in ", 
                         sum(vecboolPiEstConverged==1001), " cases."))
          }
          print(paste0("# ",scaIter,".   Drop-out estimation complete: ",
                       "loglikelihood of     ", scaLogLikTemp, " in ",
                       round(tm_pi["elapsed"]/60,2)," min."))
        }  
      }
      
      ##### 2. Gene-wise parameter estimation:
        # Estimate mean and dispersion parameters simultaneously.
        # a/b) Negative binomial mean AND dispersion parameter.
        tm_mudisp <- system.time({
          lsFitMuDisp <- fitZINBMuDisp(matCounts=matCounts,
                                       vecNormConst=vecNormConst,
                                       lsMuModel=lsMuModel,
                                       lsDispModel=lsDispModel,
                                       lsDropModel=lsDropModel,
                                       scaWindowRadius=scaWindowRadius,
                                       matWeights=matWeights)
        })
        lsDispModel$matDispModel <- lsFitMuDisp$matDispModel
        colnames(lsDispModel$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
        lsMuModel$matMuModel <- lsFitMuDisp$matMuModel
        colnames(lsMuModel$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
        
        vecboolMuEstConverged <- lsFitMuDisp$vecConvergence
        vecboolDispEstConverged <- lsFitMuDisp$vecConvergence
        
      # Evaluate Likelihood
      scaLogLikOld <- scaLogLikNew
      scaLogLikNew <- evalLogLikMatrix( matCounts=matCounts,
                                        vecNormConst=vecNormConst,
                                        lsMuModel=lsMuModel,
                                        lsDispModel=lsDispModel, 
                                        lsDropModel=lsDropModel,
                                        scaWindowRadius=scaWindowRadius )
    })
    
    # Iteration complete
    if(boolSuperVerbose){
      if(any(vecboolDispEstConverged != 0)){
        print(paste0("(Mean-) Dispersion estimation did not converge in ", 
                     sum(vecboolDispEstConverged), " cases [codes: ",
                     paste(unique(vecboolDispEstConverged[vecboolDispEstConverged!=0])), "]."))
      }
      print(paste0("# ",scaIter, ".   Mean+Disp co-estimation complete: ",
                     "loglikelihood of ", scaLogLikNew, " in ",
                     round(tm_mudisp["elapsed"]/60,2)," min."))
    } else {
      if(boolVerbose){print(paste0("# ",scaIter, ".) complete with ",
                                   "log likelihood of   ", scaLogLikNew, " in ",
                                   round(tm_iter["elapsed"]/60,2)," min."))}
    }
    vecEMLogLikModel[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  
  # Evaluate convergence
  if(all(as.logical(vecboolDispEstConverged)) &
     all(as.logical(vecboolMuEstConverged)) &
     scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergenceModel <- TRUE
  } else { boolConvergenceModel <- FALSE }
  
  return(list(
    lsMuModel=lsMuModel,
    lsDispModel=lsDispModel,
    lsDropModel=lsDropModel,
    boolConvergenceModel=boolConvergenceModel,
    vecEMLogLikModel=vecEMLogLikModel ))
}