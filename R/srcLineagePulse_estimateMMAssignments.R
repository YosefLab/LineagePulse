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
																		matMuParam,
                                    matDispParam,
                                    matDropParam,
                                    scaNormConst,
                                    vecboolNotZero,
                                    vecboolZero){
  
  # (I) Prevent parameter shrinkage/explosion
  vecTheta[vecTheta < -log(10)*10] <- -log(10)*10
  vecTheta[vecTheta > log(10)*10] <- log(10)*10
  
  # (II) Linker functions
  # Logistic linker for probability
  vecWeights <- 1/(1+exp(-vecTheta))
  # Enforce normalisation
  vecWeights <- vecWeights/sum(vecWeights)
  
  # Loop over models:
  scaNCells <- length(vecCounts)
  
  scaLogLik <- evalLogLikCellMM_comp(vecCounts=vecCounts,
  																	 matMuParam=matMuParam,
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
																 matMuParam,
                                 matDispParam,
                                 matDropParam,
                                 scaNormConst,
                                 vecWeights,
                                 MAXIT=1000,
                                 RELTOL=10^(-8) ){
  
  # Project weights into logit space (linker) for fitting
  vecWeightInit <- -log(1/vecWeights-1)
  
  fitWeights <- tryCatch({
    optim(    
      par=vecWeightInit,
      fn=evalLogLikWeightsMMZINB_comp,
      vecCounts=vecCounts,
      matMuParam=matMuParam,
      matDispParam=matDispParam,
      matDropParam=matDropParam,
      scaNormConst=scaNormConst,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      method="BFGS",
      control=list(maxit=MAXIT,
                   reltol=RELTOL,
                   fnscale=-1) )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial ",
                 "mixture model weights: fitMMAssignmentsCell(). ",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    print(paste0("vecWeightInit ", paste(vecWeightInit,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts=vecCounts,
                               vecWeightInit=vecWeightInit,
    													 matMuParam=matMuParam, 
                               matDispParam=matDispParam,
                               scaNormConst=scaNormConst)
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # Impose sensitivity bounds
  vecWeightsInLogitSpace <- fitWeights$par
  vecWeightsInLogitSpace[vecWeightsInLogitSpace < -log(10)*10] <- -log(10)*10
  vecWeightsInLogitSpace[vecWeightsInLogitSpace > log(10)*10] <- log(10)*10
  # Extract weights from linker space values
  vecWeights <- 1/(1+exp(-vecWeightsInLogitSpace))
  # Enforce normalisation
  vecWeights <- vecWeights/sum(vecWeights)
  
  return(list( vecWeights=vecWeights,
               scaLL=fitWeights$value,
               scaConvergence=fitWeights$convergence ))
}

estimateMMAssignmentsMatrix <- function(matCounts,
																				dfAnnotation,
																				boolFixedPopulations,
                                        lsMuModel,
                                        lsDispModel,
                                        lsDropModel,
                                        vecNormConst,
                                        matWeights,
                                        MAXIT_BFGS_MM=1000,
                                        RELTOL_BFGS_MM=10^(-8 ) ){
  
  
	scaNumGenes <- dim(matCounts)[1]
	scaNumCells <- dim(matCounts)[2]
  scaNumMixtures <- dim(matWeights)[2]
  
  # Select non a-priori fixed cells
  if(boolFixedPopulations) vecidxToFit <- which(is.na(dfAnnotation$populations))
  else vecidxToFit <- seq(1, scaNumCells) # Select all cells
  
  # Parallelise weight estimation over cells
  lsFitsWeights <- bplapply( vecidxToFit, function(j){
  	
  	matMuParam <- do.call(cbind, lapply(seq(1,dim(matWeights)[2]), function(m){
  		do.call(c, lapply(seq(1,scaNumGenes), function(i){
  			decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,m],
  														lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
  														lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
  														vecInterval=j)
  		}))
  	}))
  	
    # Parameter decompression: Cell-specific
    matDropParam <- do.call(cbind, lapply(seq(1,scaNumMixtures), function(m){
      decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[j,],
                                  vecMu=matMuParam[,m],
                                  matPiConstPredictors=lsDropModel$matPiConstPredictors )
    }))
    if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){
      vecDispParam <- as.vector(lsDispModel$matDispModel)
      matDispParam <- do.call(cbind, lapply(seq(1,scaNumMixtures), function(m) vecDispParam ))
    } else {
      stop(paste0("ERROR estimateMMAssignmentsMatrix(): strDispModel=", lsDispModel$lsDispModelGlobal$strDispModel, " not recognised."))
    }
    fitWeights <- fitMMAssignmentsCell(vecCounts=matCounts[,j],
    																	 matMuParam=matMuParam,
                                       matDispParam=matDispParam,
                                       matDropParam=matDropParam,
                                       scaNormConst=vecNormConst[j],
                                       vecWeights=matWeights[j,],
                                       MAXIT=MAXIT_BFGS_MM,
                                       RELTOL=RELTOL_BFGS_MM )
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