#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++     Fit dispersion parameters of ZINB model    +++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# File divided into:
# (I) OBJECTIVES - loglikelihood returning functions called within optim in II.
# (II) FITTING COORDINATORS - functions coordinating the specific model fitting,
# including calling optim and performing error handling.
# (III) Top level auxillary function - called by fitZINB and calls II.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (I) OBJECTIVES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#' Cost function for zero-inflated negative binomial 
#' dispersion parameter fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial overdispersion on single gene given
#' the drop-out rate and negative binomial mean parameter. The
#' dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' and to growth above a threshold to avoid shrinkage of the 
#' dispersion factor to zero/ expansion to infinity.
#' 
#' @seealso Called by \code{fitDispZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuEst: (vector number of cells) Negative binomial
#'    mean parameter estimate of clusters to which cells
#'    belong.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell. 
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikDispConstZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  # Log linker function to fit positive dispersions
  scaDispEst <- exp(scaTheta)
  
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  scaDispEst[scaDispEst < .Machine$double.eps] <- .Machine$double.eps
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  scaDispEst[scaDispEst > 1/.Machine$double.eps] <- 1/.Machine$double.eps
  
  vecDispersions <- rep(scaDispEst, length(vecCounts))
  
  if(!is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikSmoothZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecMuEst,
      vecSizeFactors=vecSizeFactors,
      vecDispEst=vecDispersions, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
  } else {
    scaLogLik <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecMuEst*vecSizeFactors,
      vecDispEst=vecDispersions, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero)
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) FITTING COORDINATORS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit zero-inflated negative binomial dispersion parameter as constant
#' 
#' Fits single zero-inflated negative binomial dispersion parameter
#' as maximum likelihood estimator to a gene as a constant.
#' Guards against too large and too small dispersion parameters.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuEst: (vector number of cells) Negative binomial
#'    mean parameter estimate of clusters to which cells
#'    belong.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell. 
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitDispZINB <- function( scaDispGuess,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  fitDisp <- tryCatch({ 
    unlist(optim(
      par=log(scaDispGuess),
      fn=evalLogLikDispConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecMuEst=vecMuEst,
      vecSizeFactors=vecSizeFactors,
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000, fnscale=-1)
    )[c("par","convergence")] )
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial dispersion parameter: fitDispZINB().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecMuEst ", paste(vecMuEst,collapse=" ")))
    print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecMuEst, log(scaDispGuess), vecSizeFactors, vecDropoutRateEst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecMuEst", "logscaDispEst", "vecSizeFactors", "vecDropoutRateEst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  scaDisp <- exp(fitDisp[1])
  
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < .Machine$double.eps){scaDisp <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){scaDisp <- 1/.Machine$double.eps}
  
  scaConvergence <- fitDisp[2]
  
  return(list(scaDisp=scaDisp,
    scaConvergence=scaConvergence))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Top level auxillary function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Chose mode of dispersion parameter estimation
#' 
#' Auxillary function that calls the estimation functions for the
#' different dispersion models according to their needs.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' 
#' @return 
#' @export

fitZINBDisp <- function( matCountsProc,
  vecSizeFactors,
  lsMuModel,
  lsDispModel,
  lsDropModel,
  scaWindowRadius ){
  
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  
  if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){  
    lsFitDisp <- bplapply(seq(1,scaNumGenes), function(i){
      # Decompress parameters
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecInterval=NULL )
      vecDropoutParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
        vecMu=vecMuParam,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
      
      # Estimate constant dispersion factor
      fitDisp <- fitDispZINB(scaDispGuess=lsDispModel$matDispModel[i,],
        vecCounts=matCountsProc[i,],
        vecSizeFactors=vecSizeFactors,
        vecMuEst=vecMuParam,
        vecDropoutRateEst=vecDropoutParam,
        vecboolNotZeroObserved= !is.na(matCountsProc[i,]) & matCountsProc[i,]>0, 
        vecboolZero= matCountsProc[i,]==0,
        scaWindowRadius=scaWindowRadius )
      return(fitDisp)
    })
    matDispModel <- do.call(rbind, lapply(lsFitDisp, function(i){ i$scaDisp }))
    vecConvergence <- sapply(lsFitDisp, function(i) i$scaConvergence)
    
  } else {
    #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
    print(paste0("Dispersion parameter model not recognised: ", 
      lsDispModel$lsDispModelGlobal$strDispModel, 
      ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
    stop(paste0("Dispersion parameter model not recognised: ", 
      lsDispModel$lsDispModelGlobal$strDispModel, 
      ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
  }
  
  return( list(matDispModel=matDispModel,
    vecConvergence=vecConvergence) )
}