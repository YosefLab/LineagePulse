#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++    Estimate mixture assignments  ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model
#' weight estimation for mixture models on one cell
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of the mixture assignment probability vector for one cell. 
#' 
#' The probabilities are fit in logit space for numerical 
#' reasons.
#' 
#' @param vecTheta: (numeric vector length numer of
#'    mixture components) 
#'    Weight (/mixture assignment) vector of probabilities
#'    in logit space.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @parm matMu (numeric matrix genes x mixture components)
#'    Mean parameter estimates in current model.
#' @parm matDisp (numeric matrix genes x mixture components)
#'    Dispersion parameter estimates in current model.
#' @parm matMu (numeric matrix genes x mixture components)
#'    Drop-out rate estimates in current model.
#' @param scaNormConst: (scalar) 
#'    Model scaling factor of current cell.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial log-likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikWeightsMMZINB <- function(vecTheta,
                                    vecCounts,
                                    matMu,
                                    matDispParam,
                                    matDropParam,
                                    scaNormConst,
                                    vecboolNotZero,
                                    vecboolZero){ 
  
  # (I) Linker functions
  # Logistic linker for probability
  vecWeights <- 1/(1+exp(-vecTheta))
  
  # Loop over models:
  scaNCells <- length(vecCounts)
  
  scaLogLik <- evalLogLikCellMM_comp(vecCounts=vecCounts,
                                     matMuModel=matMu,
                                     matDispParam=matDispParam,
                                     matDropParam=matDropParam,
                                     scaNormConst=scaNormConst,
                                     vecWeights=vecWeights,
                                     vecboolNotZero= !is.na(vecCounts) & vecCounts>=0, 
                                     vecboolZero= !is.na(vecCounts) & vecCounts==0 )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

evalLogLikWeightsMMZINB_comp <- cmpfun(evalLogLikWeightsMMZINB)

fitMMAssignmentsCell <- function(vecCounts,
                                 matMu,
                                 matDispParam,
                                 matDropParam,
                                 scaNormConst,
                                 vecWeights){
  
  # Project weights into logit space (linker) for fitting
  vecWeightInit <- -log(1/vecWeights-1)
  
  fitWeights <- tryCatch({
    optim(    
      par=vecWeightInit,
      fn=evalLogLikWeightsMMZINB_comp,
      vecCounts=vecCounts,
      matMu=matMu,
      matDispParam=matDispParam,
      matDropParam=matDropParam,
      scaNormConst=scaNormConst,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) &vecCounts==0,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial ",
                 "mixture model weights: fitMMAssignmentsCell(). ",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # Extract weights from linker space values
  vecWeights <- 1/(1+exp(-fitWeights$par))
  
  return(list( vecWeights=vecWeights,
               scaLL=fitWeights$value,
               scaConvergence=fitWeights$convergence ))
}

estimateMMAssignmentsMatrix <- function(matCounts,
                                        vecFixedAssignments=NULL,
                                        lsMuModel,
                                        lsDispModel,
                                        lsDropModel,
                                        vecNormConst,
                                        matWeights ){
  
  scaNumCells <- dim(matCounts)[2]
  scaNumMixtures <- dim(matWeights)[2]
  
  # Select non a-priori fixed cells
  if(!is.null(vecFixedAssignments)) vecidxToFit <- which(!is.na(vecFixedAssignments))
  else vecidxToFit <- seq(1, scaNumCells) # Select all cells
  
  # Parallelise weight estimation over cells
  lsFitsWeights <- bplapply( vecidxToFit, function(j){
    # Parameter decompression: Cell-specific
    matDropParam <- do.call(cbind, lapply(seq(1,scaNumMixtures), function(m){
      decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[j,],
                                  vecMu=lsMuModel$matMuModel[,m],
                                  matPiConstPredictors=lsDropModel$matPiConstPredictors )
    }))
    if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){
      vecDispParam <- as.vector(lsDispModel$matDispModel)
      matDispParam <- do.call(cbind, lapply(seq(1,scaNumMixtures), function(m) vecDispParam ))
    } else {
      stop(paste0("ERROR estimateMMAssignmentsMatrix(): strDispModel=", lsDispModel$lsDispModelGlobal$strDispModel, " not recognised."))
    }
    fitWeights <- fitMMAssignmentsCell(vecCounts=matCounts[,j],
                                       matMu=lsMuModel$matMuModel,
                                       matDispParam=matDispParam,
                                       matDropParam=matDropParam,
                                       scaNormConst=vecNormConst[j],
                                       vecWeights=matWeights[j,])
    return(fitWeights)
  })
  matWeightsToFit <- do.call(rbind, lapply(lsFitsWeights, function(fit) fit$vecWeights ))
  matWeights[vecidxToFit,] <- matWeightsToFit
  vecLLFit <- sapply(lsFitsWeights, function(fit) fit$scaLL )
  vecConvcergenceFit <- sapply(lsFitsWeights, function(fit) fit$scaConvergence )
  
  return(list( matWeights=matWeights,
               vecLL=vecLLFit,
               vecConvcergence=vecConvcergenceFit ))
}