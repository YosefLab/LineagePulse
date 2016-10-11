#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++     Fit dispersion parameters of ZINB model    +++++++++++++++++++#
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
#' @seealso Called by \code{fitDispZINB_LinPulse}.
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
#' @param boolSmoothed: (bool) [Default FALSE]
#'    Whether to used scaWindowRadius smoothed cost function
#'    or plain cost function.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikDispZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL,
  boolSmoothed=FALSE ){ 
  
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
  
  if(boolSmoothed){
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

#' Fit zero-inflated negative binomial dispersion parameter
#' 
#' Fits single zero-inflated negative binomial dispersion parameter
#' as maximum likelihood estimator to a gene.
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

fitDispZINB_LinPulse <- function( scaDispGuess,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL,
  boolSmoothed=FALSE ){ 
  
  fitDisp <- tryCatch({ unlist(optim(
    par=log(scaDispGuess),
    fn=evalLogLikDispZINB_LinPulse_comp,
    vecCounts=vecCounts,
    vecMuEst=vecMuEst,
    vecSizeFactors=vecSizeFactors,
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero,
    scaWindowRadius=scaWindowRadius,
    boolSmoothed=boolSmoothed,
    method="BFGS",
    control=list(maxit=1000, fnscale=-1)
  )[c("par","convergence")] )
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial dispersion parameter: fitDispZINB_LinPulse().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecMuEst ", paste(vecMuEst,collapse=" ")))
    print(paste0("vecDropout ", paste(vecDropout,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecMuEst, log(scaDispGuess), vecSizeFactors, vecDropout)
    names(lsErrorCausingGene) <- c("vecCounts", "vecMuEst", "logscaDispEst", "vecSizeFactors", "vecDropout")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })

  fitDisp["par"] <- exp(fitDisp["par"])
  
  # Catch boundary of likelihood domain on dispersion space
  if(fitDisp["par"] < .Machine$double.eps){fitDisp["par"] <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(fitDisp["par"] > 1/.Machine$double.eps){fitDisp["par"] <- 1/.Machine$double.eps}
  
  return(fitDisp)
}