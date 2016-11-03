#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Fit mean parameters of ZINB model    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# File divided into:
# (I) OBJECTIVES - loglikelihood returning functions called within optim in II.
# (II) FITTING COORDINATORS - functions coordinating the specific model fitting,
# including calling optim and performing error handling.
# (III) Top level auxillary function - called by fitZINB and calls II.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (I) OBJECTIVES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter.
#' The mean is modelled by cell and constrained to fit a sliding 
#' window of cells. This is the objective for the constant mean model,
#' either used for all cells of a gene or all cells withing a cluster
#' of a gene.
#' The mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of mean parameter estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDispEst: (scalar vector number of samples) 
#'    Negative binomial dispersion  parameter for given 
#'    gene and observations.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuConstZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  # (I) Linker functions
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  # (II) Compute drop-out rates
  vecDropoutRateEst <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
    vecMu=rep(scaMu, dim(matDropoutLinModel)[1]),
    vecPiConstPredictors=vecPiConstPredictors )
  #vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
  #  sum(matDropoutLinModel[cell,] * c(1,log(scaMu),vecPiConstPredictors))
  #})
  #vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (III) Evaluate loglikelihood (this is the cost function) 
  if(is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=scaMu*vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    scaLogLik <- evalLogLikSmoothZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=rep(scaMu, length(vecCounts)),
      vecSizeFactors=vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' under sliding windows mean model for single mean (coordinate ascent)
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter.
#' The mean is modelled by cell and constrained to fit a sliding 
#' window of cells. This function is used if the means for 
#' observations of a gene are fit sequentially (coordinate ascent).
#' 
#' The mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' mean parameter shrinking beyond a numerical threshold to zero.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of mean parameter estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecMu: (vector number of cells in neighbourhood) Mean
#'    parameter estimates in neighourhood of target cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter. 
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed) and has non-zero count.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaTarget: (integer) Position of target cell
#'    (whose mean is to be estimated) within given interval.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuWindowZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecMu,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaTarget,
  scaWindowRadius){ 
  
  # (I) Linker functions
  scaMuNew <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMuNew < .Machine$double.eps){ scaMuNew <- .Machine$double.eps }
  
  scaN <- length(vecCounts)
  # (II) Compute drop-out rates
  vecMuNeighbourhoodCurrent <- vecMu # Means which are not updated
  vecMuNeighbourhoodCurrent[scaTarget] <- scaMuNew # Mean which is updated
  vecDropoutRateEst <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
    vecMu=vecMuNeighbourhoodCurrent,
    vecPiConstPredictors=vecPiConstPredictors )
  #vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
  #  sum(matDropoutLinModel[cell,] * c(1,log(vecMuNeighbourhoodCurrent[cell]),vecPiConstPredictors))
  #})
  #vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (III) Evaluate loglikelihood (this is the cost function) 
  # Compute loglikelihood terms with contribution from new mean parameter
  # Cis terms: From neighbourhood around target parameter
  scaLogLikCis <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
    vecMu=scaMuNew*vecNormConst,
    vecDispEst=rep(vecDisp[scaTarget],scaN), 
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero )
  # Trans terms: From other neighbourhoods via drop-out rate as function
  # of target mean parameter.
  if(scaN>1){
    if(scaTarget==1){ vecindTransTerms <- seq(2,scaN)
    } else if(scaTarget==scaN){ vecindTransTerms <- seq(1,scaN-1)
    } else { vecindTransTerms <- c(seq(1,scaTarget-1),seq(scaTarget+1,scaN)) }
    scaLogLikTrans <- evalLogLikZINB_LinPulse_comp( vecCounts=rep(vecCounts[scaTarget],length(vecindTransTerms)),
      vecMu=vecMu[vecindTransTerms]*vecNormConst[scaTarget],
      vecDispEst=vecDisp[vecindTransTerms], 
      vecDropoutRateEst=rep(vecDropoutRateEst[scaTarget],length(vecindTransTerms)),
      vecboolNotZeroObserved=rep(vecboolNotZeroObserved[scaTarget],length(vecindTransTerms)), 
      vecboolZero=rep(vecboolZero[scaTarget],length(vecindTransTerms)) )
  } else {
    scaLogLikTrans<- 0
  }
  scaLogLik <- scaLogLikCis + scaLogLikTrans
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' under sliding windows mean model for multiple means (BFGS)
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter.
#' The mean is modelled by cell and constrained to fit a sliding 
#' window of cells. This function is used if the means for 
#' observations of a gene are fit simulatneously (BFGS).
#' 
#' The mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' mean parameter shrinking beyond a numerical threshold to zero.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecTheta: (numeric vector number of means to be estimated) 
#'    Log of mean parameter estimates.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter. 
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed) and has non-zero count.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuVecWindowsZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecDisp,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius){ 
  
  # (I) Linker functions
  vecMu <- exp(vecTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  # (II) Compute drop-out rates
  vecDropoutRateEst <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
    vecMu=vecMu,
    vecPiConstPredictors=vecPiConstPredictors )
  #vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
  #  sum(matDropoutLinModel[cell,] * c(1,log(vecMu[cell]),vecPiConstPredictors))
  #})
  #vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (III) Evaluate loglikelihood (this is the cost function) 
  scaLogLik <- evalLogLikSmoothZINB_LinPulse_comp(
    vecCounts=vecCounts,
    vecMu=vecMu,
    vecSizeFactors=vecNormConst,
    vecDispEst=vecDisp, 
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero,
    scaWindowRadius=scaWindowRadius)
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model mean co-estimation under constant dispersion and constant mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate/model. The mean parameter is modelled with 
#' an impulse model for the given gene.
#' The impulse model amplitude parameters are fit in log space and are
#' therefore fit as positive scalars. The cost function is insensitive to the
#' amplitude parameters shrinking beyond a numerical threshold to zero which 
#' may cause numerical errors.
#' 
#' @aliases evalLogLikMuImpulseZINB_LinPulse_comp 
#' 
#' @seealso Called by
#' \code{fitDispConstMuImpulseOneInitZINB}.
#' 
#' @param vecTheta: (numeric vector impulse parameters (6)) 
#'    Impulse model parameter estimates.
#' @param vecCounts: (integer vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
#'    Unique (pseudo)time coordinates of cells.
#' @param vecDisp: (numeric vector number of cells)
#'    Dispersion parameter estimates for observations of a gene.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to cell in list of sorted
#'    time points. vecTimepoints[vecindTimepointAssign]==vecPseudotime
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikMuImpulseZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecTimepoints,
  vecDisp,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecindTimepointAssign, 
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL){  
  
  # (I) Linker functions
  # Log linker for amplitudes and catching of low model values is in evalImpulseModel_comp
  vecImpulseValue <- evalImpulseModel_comp(vecTheta,vecTimepoints)[vecindTimepointAssign]
  
  # (II) Compute drop-out rates
  vecDropoutRateEst <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
    vecMu=vecImpulseValue,
    vecPiConstPredictors=vecPiConstPredictors )
  #vecLinModelOut <- sapply(seq(1,length(vecImpulseValue)), function(cell){
  #  sum(matDropoutLinModel[cell,] * c(1,log(vecImpulseValue[cell]),vecPiConstPredictors))
  #})
  #vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (III) Evaluate loglikelihood of estimate
  if(is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=vecImpulseValue*vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    scaLogLik <- evalLogLikSmoothZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=vecImpulseValue,
      vecSizeFactors=vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
  }
  
  return(scaLogLik)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) FITTING COORDINATORS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Fits single negative binomial mean parameters numerically as 
#' maximum likelihood estimator to a gene: Constant mean model.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaMuGuess: (scalar) Initialisation for mean parameter
#'    to be estimated.
#' @param scaDispEst: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitMuConstZINB <- function(vecCounts,
  scaMuGuess,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaWindowRadius){ 

  # Numerical maximum likelihood estimator
  fitMu <- tryCatch({
    unlist(optim(
      par=log(scaMuGuess),
      evalLogLikMuConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuConstZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("scaMeanGuess ",scaMeanGuess))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaMu <- exp(fitMu[1])
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  scaConvergence <- fitMu[2]
  
  return(list(scaMu=scaMu,
    scaConvergence=scaConvergence))
}

#' Fit negative binomial mean to single cell in 
#' an interval (cluster or sliding windowmode, coordinate ascent)
#' 
#' Fits negative binomial mean parameter numerically as maximum likelihood estimator
#' to a cell in an interval. The routine distinguishes fitting to a cluster
#' and fitting to a sliding window. Moreover, the routine allows for the
#' use of the closed form maximum likelihood estimator if all normalisation
#' constants are one and the drop-out rate is not fit dynamically (i.e. is
#' handled as a constant).
#' Note that this fitting routine is meant for fitting a single mean
#' to a set of observations and can be used in wrappers which compute
#' every mean separately (assuming independence).
#' 
#' @seealso Called by \code{fitZINB}. Alternative to simultaneous
#' maximum likelihood estimation of mean parameters of all cells in
#' interval used in \code{fitMuVecZINB_LinPulse}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMu: (vector number of cells in neighbourhood)
#'    [Defaul NULL] 
#'    Mean parameter estimates in neighourhood of target cell.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (numeric vector number of cells) Posterior
#'    of observation not being drop-out for closed-form MLE.
#' @param scaTarget: (integer) Position of target cell
#'    (whose mean is to be estimated) within given interval.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuWindowZINB <- function(vecCounts,
  vecMu,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaTarget,
  scaWindowRadius=NULL ){ 
  
  fitMu <- tryCatch({
    unlist(optim(
      par=log(vecMu[scaTarget]),
      evalLogLikMuWindowZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecMu=vecMu,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaTarget=scaTarget,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuWindowZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("log(vecMu[scaTarget]) ",log(vecMu[scaTarget])))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecMu", paste(vecMu,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    print(paste0("scaTarget ", paste(scaTarget,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaMu <- exp(fitMu[1])
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  scaConvergence <- fitMu[2]
  
  return(list(scaMu=scaMu,
    scaConvergence=scaConvergence))
}

#' Fit negative binomial means to all cells in an interval (BFGS)
#' 
#' Fits negative binomial mean parameters numerically as maximum likelihood estimator
#' to all cells in an interval simultaneously using BFGS.
#' 
#' @seealso Called by \code{fitZINB}. Alternative to cell-wise sequential
#' maximum likelihood estimation used in \code{fitMuZINB}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuGuess: (vector number of cells in neighbourhood)
#'    [Defaul NULL] 
#'    Mean parameter estimates in neighourhood of target cell.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (numeric vector number of cells) Posterior
#'    of observation not being drop-out for closed-form MLE.
#' @param scaTarget: (integer) Position of target cell
#'    (whose mean is to be estimated) within given interval.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuVecWindowsZINB<- function(vecCounts,
  vecMuGuess,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaWindowRadius=NULL ){
  
  fitMu <- tryCatch({
    unlist(optim(
      par=log(vecMuGuess),
      evalLogLikMuVecWindowsZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuWindowZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("log(vecMuGuess) ",log(vecMuGuess)))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecMuGuess, vecDisp, 
      matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecMuGuess","vecDisp",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  vecMu <- exp(fitMu[1:length(vecMuGuess)])
  # Catch boundary of likelihood domain on mu space
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  scaConvergence <- fitMu[length(vecMu)+1]
  
  return(list(vecMu=vecMu,
    scaConvergence=scaConvergence))
}

#' Fit an impulse model to observations of a gene
#' 
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the impulse model and returns
#' the fitted (maximum likelihood) model.
#' 
#' @seealso Called by \code{fitImpulse_gene_LP}. This function
#' performs optimisation of one impulse model initialisation,
#' \code{fitImpulse_gene_LP} coordinates the overall fitting
#' of an impulse model to a gene.
#' 
#' @param scaDispGuess: (scalar) 
#'    Dispersion parameter estimate for given gene.
#' @param vecParamGuessPeak (numeric vector number of parameters [6]) 
#'    Impulse model parameter estimates for given gene.
#' @param vecCounts: (integer vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecTimepoints: (numerical vector number of unique time points)
#'    Unique (pseudo)time coordinates of cells.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to cell in list of sorted
#'    time points. vecTimepoints[vecindTimepointAssign]==vecPseudotime
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' @param MAXIT: (integer) [Default 1000] Maximum number of 
#'    estimation iterations optim.
#' @param RELTOL: (scalar) [Default sqrt(.Machine$double.eps)]
#'    Relative tolerance for optim.
#'   
#' @return list: (length 3)
#'    \itemize{
#'      \item vecImpulseParam: (numeric vector [beta, h0, h1, h2, t1, t2])
#'        Impulse model parameter estimates.
#'      \item scaLL: (scalar) Loglikelihood of fit.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#' @export

fitMuImpulseOneInitZINB <- function(vecImpulseParamGuess,
  vecCounts,
  vecTimepoints,
  vecDisp,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecindTimepointAssign,
  scaWindowRadius=NULL,
  MAXIT=1000,
  RELTOL=sqrt(.Machine$double.eps),
  trace=0,
  REPORT=10 ){
  
  fitDispImpulse <- tryCatch({
    unlist( optim(
      par=vecImpulseParamGuess, 
      fn=evalLogLikMuImpulseZINB_LinPulse_comp, 
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS", 
      control=list(maxit=MAXIT,
        reltol=RELTOL,
        fnscale=-1, 
        trace=trace, 
        REPORT=REPORT)
    )[c("par","value","convergence")] )
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting impulse model: fitDispConstMuImpulseZINB().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecImpulseParamGuess ", paste(vecImpulseParamGuess,collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecTimepoints ", paste(vecTimepoints,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("matDropoutLinModel ", paste(matDropoutLinModel,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
    print(paste0("scaWindowRadius", paste(scaWindowRadius,collapse=" ")))
    print(paste0("MAXIT ", MAXIT))
    lsErrorCausingGene <- list(, vecImpulseParamGuess, vecCounts, vecTimepoints, 
      vecDisp, vecNormConst, vecindTimepointAssign, scaWindowRadius, MAXIT)
    names(lsErrorCausingGene) <- c("vecImpulseParamGuess", "vecCounts", "vecTimepoints",
      "vecDisp", "vecNormConst", "vecindTimepointAssign", "scaWindowRadius", "MAXIT")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries  
  vecImpulseParam <- fitDispImpulse[1:6]
  scaLL <- fitDispImpulse[7]
  scaConvergence <- fitDispImpulse[8]
  
  return( list(vecImpulseParam=vecImpulseParam,
    scaLL=scaLL,
    scaConvergence=scaConvergence) )
}

#' Fit means as impulse model and dispersions as constant
#' 
#' Computes impulse parameter initialisation for valley
#' and peak model and uses both and the prior parameter fit
#' in three separate optimisation runs to obtain the best 
#' impulse model fit to the data, simultaneous with fitting a 
#' constant dispersion factor.
#' 
#' @seealso Called by \code{fitZINB}. Calls optimisation wrapper
#' \code{fitDispConstMuImpulseOneInitZINB} for each initialisation.
#' Code similar to \code{ImpulseDE2::fitImpulse_gene}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecImpulseParamGuess: (numerical vector impulse parameters)
#'    Previous impulse model parameter values.
#' @param lsMuModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'    \itemize{
#'      \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'    "windows"} Name of the mean model.
#'      \item scaNumCells: (scalar) [Default NA] Number of cells
#'    for which model is evaluated. Used for constant model.
#'      \item vecPseudotime: (numerical vector number of cells)
#'    [Default NA] Pseudotime coordinates of cells. Used for
#'    impulse model.
#'      \item vecindClusterAssign: (integer vector length number of
#'    cells) [Default NA] Index of cluster assigned to each cell.
#'    Used for clusters model.
#'      \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'    of a gene are simultaneously estiamted as a vector with BFGS
#'    in windows mode.
#'      \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'    for BFGS estimation of impulse model with optim (termination criterium).
#'      \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'    change in objective function for BFGS estimation of impulse 
#'    model with optim (termination criterium).
#'    }
#' @param vecDisp: (numerical vector) Negative binomial
#'    dispersion parameter estimates for each cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item vecImpulseParam: (numeric vector [beta, h0, h1, h2, t1, t2])
#'        Impulse model parameter estimates.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#' @export

fitMuImpulseZINB <- function(vecCounts,
  vecImpulseParamGuess,
  lsMuModelGlobal,
  vecDisp,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  scaWindowRadius=NULL,
  MAXIT=1000,
  RELTOL=sqrt(.Machine$double.eps) ){
  
  # Set reporter parameters for optim BFGS optimisation
  trace <- 0 # Report typ: 0 none, 2 yes
  REPORT <- 1 # Frequency of reporting in iterations
  # Try new peak and valley initialisations?
  # Increases time complexity of mean estimation by factor 3
  # but seems to make a difference on simulated data.
  boolUseNewInits <- TRUE
  
  # (I) Get unique time points to reduce number of evaluations
  # Compute time point specifc parameters
  vecTimepoints <- unique( lsMuModelGlobal$vecPseudotime )
  # Get vector of numeric time point assignment indices:
  vecindTimepointAssign <- match(lsMuModelGlobal$vecPseudotime, vecTimepoints)
  
  # (II) Initialise impulse model
  # Decompress parameters for initialisation
  vecMuParam <- decompressMeansByGene( vecMuModel=vecImpulseParamGuess,
    lsMuModelGlobal=lsMuModelGlobal,
    vecInterval=NULL )
  vecDropoutParam <- decompressDropoutRateByGene( matDropModel=matDropoutLinModel,
    vecMu=vecMuParam,
    vecPiConstPredictors=vecPiConstPredictors )
  # The previous parameter estiamte is kept as a reference and
  # used as an initialisation
  # Compute initialisations for peak and valley
  lsParamGuesses <- initialiseImpulseParametes(vecCounts=vecCounts,
    lsMuModelGlobal=lsMuModelGlobal,
    vecMu=vecMuParam,
    vecDisp=vecDisp,
    vecDrop=vecDropoutParam,
    vecNormConst=vecNormConst,
    scaWindowRadius=scaWindowRadius)
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # (III) Compute new parameters
  # 1. Initialisation: Prior best fit
  lsFitPrior <- fitMuImpulseOneInitZINB(
    vecImpulseParamGuess=vecImpulseParamGuess,
    vecCounts=vecCounts,
    vecTimepoints=vecTimepoints,
    vecDisp=vecDisp,
    matDropoutLinModel=matDropoutLinModel,
    vecPiConstPredictors=vecPiConstPredictors,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    scaWindowRadius=scaWindowRadius,
    MAXIT=lsMuModelGlobal$MAXIT, 
    RELTOL=lsMuModelGlobal$RELTOL, 
    trace=trace, REPORT=REPORT )  
  if(boolUseNewInits){
    # 2. Initialisation: Peak
    lsFitPeak <- fitMuImpulseOneInitZINB(
      vecImpulseParamGuess=vecParamGuessPeak,
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=lsMuModelGlobal$MAXIT, 
      RELTOL=lsMuModelGlobal$RELTOL,  
      trace=trace, REPORT=REPORT )
    # 3. Initialisation: Valley
    lsFitValley <- fitMuImpulseOneInitZINB(
      vecImpulseParamGuess=vecParamGuessValley,
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=lsMuModelGlobal$MAXIT, 
      RELTOL=lsMuModelGlobal$RELTOL,  
      trace=trace, REPORT=REPORT )
    
    # (IV) Find best fit
    lsFits <- list(lsFitPeak, lsFitValley, lsFitPrior)
    vecLL <- sapply(lsFits , function(fit) fit$scaLL)
    indMaxLL <- match(max(vecLL), vecLL)
    lsFitBest <- lsFits[[indMaxLL]]
  } else {
    lsFitBest <- lsFitPrior
  }
  
  vecImpulseParam <- lsFitBest$vecImpulseParam
  scaConvergence <- lsFitBest$scaConvergence
  
  # THIS CODE IS ONLY FOR DEVELOPERS TO DEBUG IMPULSE FITTING
  # Follow the choice of model and loglikelihoods
  # Last time this code flagged a convergence problem, a parameter
  # was assigned wrongly in the function calls above (vecPseudotime).s
  if(FALSE){
    # Compute predicted means
    vecImpulseValue <- evalImpulseModel_comp(vecImpulseParam,vecTimepoints)[vecindTimepointAssign]
    vecImpulseValueOld <- evalImpulseModel_comp(vecImpulseParamGuess,vecTimepoints)[vecindTimepointAssign]
    
    # Compute best set from last iteration
    vecLinModelOutOld <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValueOld[cell])) * matDropoutLinModel[cell,])
    })
    vecDropoutOld <- 1/(1+exp(-vecLinModelOutOld))
    scaLLOld <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecImpulseValueOld*vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutOld,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0, 
      vecboolZero= !is.na(vecCounts) & vecCounts==0 )
    
    # report all new parame
    vecLinModelOut <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValue[cell])) * matDropoutLinModel[cell,])
    })
    vecDropout <- 1/(1+exp(-vecLinModelOut))
    scaLLRef <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecImpulseValue*vecNormConst,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropout,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0, 
      vecboolZero= !is.na(vecCounts) & vecCounts==0 )
    print(paste0("Old:", scaLLOld, " ,New recomputed: ", 
      scaLLRef, " ,New from optim: ", lsFitBest$scaLL))
  }
  
  return(list(vecImpulseParam=vecImpulseParam,
    scaConvergence=scaConvergence))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Top level auxillary function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Chose mode of mean parameter estimation
#' 
#' Auxillary function that calls the estimation functions for the
#' different mean models according to their needs.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param lsMuModel: (list length 2)
#'    All objects necessary to compute mean parameters for all
#'    observations.
#'    \itemize{
#'      \item matMuModel: (numerical matrix genes x number of model parameters)
#'    Parameters of mean model for each gene.
#'      \item lsMuModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'        \itemize{
#'          \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'        "windows"} Name of the mean model.
#'          \item scaNumCells: (scalar) [Default NA] Number of cells
#'        for which model is evaluated. Used for constant model.
#'          \item vecPseudotime: (numerical vector number of cells)
#'        [Default NA] Pseudotime coordinates of cells. Used for
#'        impulse model.
#'          \item vecindClusterAssign: (integer vector length number of
#'        cells) [Default NA] Index of cluster assigned to each cell.
#'        Used for clusters model.
#'          \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'        of a gene are simultaneously estiamted as a vector with BFGS
#'        in windows mode.
#'          \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'        for BFGS estimation of impulse model with optim (termination criterium).
#'          \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'        change in objective function for BFGS estimation of impulse 
#'        model with optim (termination criterium).
#'      }
#'    }
#' @param lsDispModel: (list length 2)
#'    All objects necessary to compute dispersion parameters for all
#'    observations.
#'    \itemize{
#'      \item matDispModel: (numerical matrix genes x number of model parameters)
#'    Parameters of dispersion model for each gene.
#'      \item lsDispModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'        \itemize{
#'          \item strDispModel: (str) {"constant"} 
#'        Name of the dispersion model.
#'          \item scaNumCells: (scalar) [Default NA] Number of cells
#'        for which model is evaluated. Used for constant model.
#'          \item vecPseudotime: (numerical vector number of cells)
#'        [Default NA] Pseudotime coordinates of cells. Used for
#'        impulse model.
#'          \item vecindClusterAssign: (integer vector length number of
#'        cells) [Default NA] Index of cluster assigned to each cell.
#'        Used for clusters model.
#'      }
#'    }
#' @param lsDropModel: (list length 2)
#'    All objects necessary to compute drop-out parameters for all
#'    observations, omitting mean parameters (which are stored in lsMeanModel).
#'    \itemize{
#'      \item matDropoutLinModel: (numeric matrix cells x number of model parameters)
#'    {offset parameter, log(mu) parameter, parameters belonging to
#'    constant predictors}
#'    Parameters of drop-out model for each cell
#'      \item matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied.
#'    }
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' 
#' @return list (length 2)
#'    \itemize{
#'      \item matMuModel:
#'      \item vecConvergence:
#'    }
#' @export

fitZINBMu <- function( matCountsProc,
  vecSizeFactors,
  lsMuModel,
  lsDispModel,
  lsDropModel,
  scaWindowRadius ){

  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  if(lsMuModel$lsMuModelGlobal$strMuModel=="windows"){
    # Estimate mean parameter for each cell as ZINB model for cells within pseudotime
    # interval with cell density centred at target cell.
    # Note that this corresponds to maximising smoothed log likelihood but instead
    # of using the implemented cost function evalLogLikSmoothZINB_LinPulse for an entire
    # gene, the optimisation problem is broken up into 1D problems for each mean.
    lsFitMu <- bplapply(seq(1,scaNumGenes), function(i){
      # Note: Mean parameter estimates of gene i depend on each other:
      # Either estimate all parameters for gene i together 
      # (quasi-Newton estimation with BFGS: boolVecWindowsAsBFGS=TRUE)
      # or use the latest updates of the remaining parameters during 
      # one-by-one estimation (coordinate ascent, boolVecWindowsAsBFGS=TRUE).
      
      # Decompress parameters
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecInterval=NULL )
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      if(lsMuModel$lsMuModelGlobal$boolVecWindowsAsBFGS){
        fitMu <- fitMuVecWindowsZINB(
          vecCounts=matCountsProc[i,],
          vecMu=vecMuParam,
          vecDisp=vecDispParam,
          vecNormConst=vecSizeFactors,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius )
      } else {
        vecMu <- vecMuParam
        scaConvergence <- 0
        for(j in seq(1,scaNumCells)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          fitMu <- fitMuWindowZINB(vecCounts=matCountsProc[i,vecInterval],
            vecMu=vecMu[vecInterval],
            vecDisp=vecDispParam[vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
            vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
            scaTarget=match(j,vecInterval),
            scaWindowRadius=scaWindowRadius )
          # Update parameter vector for next step of coordinate ascent
          # within this mu vector estimation.
          vecMu[j] <- fitMu$scaMu
          if(fitMu$scaConvergence){ scaConvergence <- fitMu$scaConvergence }
        }
        fitMu <- list(vecMu=vecMu,
          scaConvergence=scaConvergence)
      }
      return(fitMu)
    })
    matMuModel <- do.call(rbind, lapply(lsFitMu, function(i) i$vecMu ))
    vecConvergence <- sapply(lsFitMu,  function(i) i$scaConvergence)
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="clusters"){
    # Estimate mean parameter by cluster. No smoothing is used.
    lsFitMu <- bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate mean parameters
      fitMu <- lapply(unique(lsMuModel$lsMuModelGlobal$vecindClusterAssign), function(k){
        vecInterval <- lsMuModel$lsMuModelGlobal$vecindClusterAssign==k
        fitMuCluster <- fitMuConstZINB(
          vecCounts=matCountsProc[i,vecInterval],
          scaMuGuess=lsMuModel$matMuModel[i,k],
          vecDisp=vecDispParam[vecInterval],
          matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          vecNormConst=vecSizeFactors[vecInterval],
          scaWindowRadius=NULL )
        return(fitMuCluster)
      })
      vecMu <- sapply(fitMu, function(k) k$scaMu )
      scaConvergence <- max(sapply(fitMu, function(k) k$scaConvergence ))

      return(list(vecMu=vecMu,
          scaConvergence=scaConvergence))
    })
    matMuModel <- do.call(rbind, lapply(lsFitMu, function(i) i$vecMu ))
    vecConvergence <- sapply(lsFitMu,  function(i) i$scaConvergence)

  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){
    
    lsFitMu <- bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate mean parameters
      fitMu <- fitMuImpulseZINB(
        vecCounts=matCountsProc[i,],
        vecImpulseParamGuess=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecDisp=vecDispParam,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        vecNormConst=vecSizeFactors,
        scaWindowRadius=scaWindowRadius )
      return(fitMu)
    })
    matMuModel <- do.call(rbind, lapply(lsFitMu, function(i) i$vecImpulseParam ))
    vecConvergence <- sapply(lsFitMu,  function(i) i$scaConvergence)
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="constant"){
    lsFitMu <- bplapply( seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate constant mean parameter
      fitMu <- fitMuConstZINB( vecCounts=matCountsProc[i,],
        scaMuGuess=lsMuModel$matMuModel[i,],
        vecDisp=vecDispParam,
        vecNormConst=vecSizeFactors,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        scaWindowRadius=scaWindowRadius )
      return(fitMu)
    })
    matMuModel <- do.call(rbind, lapply(lsFitMu, function(i) i$scaMu ))
    vecConvergence <- sapply(lsFitMu,  function(i) i$scaConvergence)
  }
  
  return( list(matMuModel=matMuModel,
    vecConvergence=vecConvergence) )
}