#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++     Co-Fit mean and dispersion parameters of ZINB model    ++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Note that this file handles co-estimation of mean and dispersion model, as oppose
# to sequential estimation as handled by the respective files.
# File divided into:
# (I) OBJECTIVES - loglikelihood returning functions called within optim in II.
# (II) FITTING COORDINATORS - functions coordinating the specific model fitting,
# including calling optim and performing error handling.
# (III) Top level auxillary function - called by fitZINB and calls II.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (I) OBJECTIVES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and constant mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersionmean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled as a constant for the given gene.
#' 
#' The mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' mean factor shrinking beyond a numerical threshold to zero which 
#' may cause numerical errors.
#' The dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' and to growth above a threshold to avoid shrinkage of the 
#' dispersion factor to zero/ expansion to infinity.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecTheta: (numeric vector length 2) Log of dispersion parameter 
#'    estimate and log of mean parameter estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
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

evalLogLikDispConstMuConstZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  scaMu <- exp(vecTheta[2])
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < .Machine$double.eps){ scaDisp <- .Machine$double.eps }
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){ scaDisp <- 1/.Machine$double.eps }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  # (III) Compute drop-out rates
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(scaMu),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (IV) Evaluate loglikelihood of estimate
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
#' @param vecTheta: (numeric vector 1 + number of cells) 
#'    Log of dispersion and log of mean parameter estimates.
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

evalLogLikDispConstMuVecWindowsZINB_LinPulse <- function(vecTheta,
  vecCounts,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2:length(vecTheta)])
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < .Machine$double.eps){ scaDisp <- .Machine$double.eps }
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){ scaDisp <- 1/.Machine$double.eps }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  # (III) Compute drop-out rates
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(vecMu[cell]),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (IV) Evaluate loglikelihood (this is the cost function) 
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

#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and cluster mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersionmean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled as by cluster for the given gene.
#' 
#' The mean parameters are fit in log space and is therefore fit
#' as positive scalars. The cost function is insensitive to the
#' mean factors shrinking beyond a numerical threshold to zero which 
#' may cause numerical errors.
#' The dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' and to growth above a threshold to avoid shrinkage of the 
#' dispersion factor to zero/ expansion to infinity.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates which are optimised.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
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

evalLogLikDispConstMuClustersZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecboolNotZeroObserved,
  vecboolZero,
  vecindClusterAssign ){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2:length(vecTheta)])
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < .Machine$double.eps){ scaDisp <- .Machine$double.eps }
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){ scaDisp <- 1/.Machine$double.eps }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  vecMuParam <- vecMu[vecindClusterAssign]
  
  # (III) Compute drop-out rates
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(vecMuParam[cell]),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (IV) Evaluate loglikelihood of estimate
  scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
    vecMu=vecMuParam*vecNormConst,
    vecDispEst=vecDisp, 
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and constant mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersionmean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled as a constant for the given gene.
#' The impulse model amplitude parameters are fit in log space and are
#' therefore fit as positive scalars. The cost function is insensitive to the
#' amplitude parameters shrinking beyond a numerical threshold to zero which 
#' may cause numerical errors.
#' The dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' and to growth above a threshold to avoid shrinkage of the 
#' dispersion factor to zero/ expansion to infinity.
#' 
#' @aliases evalLogLikDispConstMuImpulseZINB_LinPulse_comp
#' 
#' @seealso Called by \code{fitImpulse}::\code{fitImpulse_matrix}::
#' \code{fitImpulse_gene}::\code{optimiseImpulseModelFit}.
#' Calls \code{evalImpulseModel} and \code{evalLogLikZINB_comp}.
#' 
#' @param vecTheta: (numeric vector dispersion (1) and impulse parameters (6)) 
#'    Dispersion model and impulse model parameter estimates.
#' @param vecCounts: (integer vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
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
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikDispConstMuImpulseZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecTimepoints,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecindTimepointAssign, 
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL){  
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker for amplitudes and catching of low model values is in evalImpulseModel_comp
  vecImpulseValue <- evalImpulseModel_comp(vecTheta[2:7],vecTimepoints)[vecindTimepointAssign]
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < .Machine$double.eps){ scaDisp <- .Machine$double.eps }
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){ scaDisp <- 1/.Machine$double.eps }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # (III) Compute drop-out rates
  vecLinModelOut <- sapply(seq(1,length(vecImpulseValue)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(vecImpulseValue[cell]),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
  # (IV) Evaluate loglikelihood of estimate
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

#' Cost function zero-inflated negative binomial model for mean fitting under
#' constant dispersion and constant mean model
#' 
#' Fits single negative binomial mean and dispersion parameter numerically as 
#' maximum likelihood estimators to a gene: Constant dispersion and
#' constant mean model.
#' 
#' @seealso Called by \code{fitZINBMuDisp}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param scaMuGuess: (scalar) Initialisation for mean parameter
#'    to be estimated.
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
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitDispConstMuConstZINB <- function(vecCounts,
  scaDispGuess,
  scaMuGuess,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaWindowRadius){ 
  
  # (I) Numerical maximum likelihood estimation
  fitDispMu <- tryCatch({
    unlist(optim(
      par=c(log(scaDispGuess), log(scaMuGuess)),
      evalLogLikDispConstMuConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuConstZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("scaMeanGuess ",scaMeanGuess))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < .Machine$double.eps){scaDisp <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){scaDisp <- 1/.Machine$double.eps}
  
  scaMu <- exp(fitDispMu[2])
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  scaConvergence <- fitDispMu[3]
  
  return(list(scaDisp=scaDisp,
    scaMu=scaMu,
    scaConvergence=scaConvergence))
}

#' Fit negative binomial means and dispersion to all observations of a gene
#' 
#' Fits negative binomial mean parameters and dispersion parameter
#' to all observations of a gene. The mean parameter is modelled
#' by local smoothing (windows) and the dispersion parameter as a constant.
#' 
#' @seealso Called by \code{fitZINBMuDisp}. Alternative to cell-wise sequential
#' maximum likelihood estimation used in \code{fitDispConstMuWindowZINB}.
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

fitDispConstMuVecWindowsZINB<- function(vecCounts,
  scaDispGuess,
  vecMuGuess,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaWindowRadius=NULL ){
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=c(log(scaDispGuess), log(vecMuGuess)),
      evalLogLikDispConstMuVecWindowsZINB_LinPulse_comp,
      vecCounts=vecCounts,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuVecWindowsZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("log(vecMuGuess) ",log(vecMuGuess)))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecMuGuess, scaDispGuess, 
      matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecMuGuess","scaDispGuess",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < .Machine$double.eps){scaDisp <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){scaDisp <- 1/.Machine$double.eps}
  
  vecMu <- exp(fitDispMu[2:(length(vecMuGuess)+1)])
  # Catch boundary of likelihood domain on mu space
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  scaConvergence <- fitDispMu[length(fitDispMu)]
  
  return(list(scaDisp=scaDisp,
    vecMu=vecMu,
    scaConvergence=scaConvergence))
}

#' Fit negative binomial mean parameter and dispersion parameter
#' to a cluster of cells
#' 
#' Fit negative binomial mean parameter and dispersion parameter
#' to a cluster of cells. This is the wrapper performing numerical
#' optimisation with optim. Note that this cannot be performed with
#' fitDispConstMuConstZINB for each cluster as the dispersion
#' parameter is assumed to be shared between clusters.
#' 
#' @seealso Called by \code{fitZINBMu}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param scaMuGuess: (scalar) 
#'    Mean parameter estimate for cluster.
#' @param scaDispGuess: (scalar) 
#'    Dispersion parameter estimate for given gene.
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
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Dispersion parameter estimate. 
#'      \item vecMu: (numeric vector number of clusters)
#'        Mean parameter estimates for each cluster.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#' @export

fitDispConstMuClusterZINB <- function(vecCounts,
  vecMuGuess,
  scaDispGuess,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecindClusterAssign){ 
  
  # scaWindowRadius is set to NULL because smoothing
  # within clusters does't make sense - the clusters already impose
  # a constraint on the means.
  fitDispMu <- tryCatch({
    unlist(optim(    
      par=c(log(scaDispGuess), log(vecMuGuess)),
      evalLogLikDispConstMuClustersZINB_LinPulse_comp,
      vecCounts=vecCounts,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindClusterAssign=vecindClusterAssign,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )[c("par","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuClusterZINB().",
      " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(paste0("log(vecMuGuess) ",log(log(vecMuGuess))))
    print(paste0("scaDispGuess ", paste(scaDispGuess,collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print("matDropoutLinModel")
    print(matDropoutLinModel)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, scaDispGuess, matDropoutLinModel, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "scaDispGuess",
      "matDropoutLinModel", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < .Machine$double.eps){scaDisp <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){scaDisp <- 1/.Machine$double.eps}
  
  vecMu <- exp(fitDispMu[2:(length(vecMuGuess)+1)])
  # Catch boundary of likelihood domain on mu space
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  scaConvergence <- fitDispMu[length(fitDispMu)]
  
  return(list(scaDisp=scaDisp,
    vecMu=vecMu,
    scaConvergence=scaConvergence))
}

#' Fit an impulse modeland a constant dispersion parameter
#' to a gene
#' 
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the impulse model and a constant
#' dispersion parameter and returns the fitted (maximum likelihood) model.
#' This is the wrapper that calls optim.
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
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
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
#' @return list: (length 4)
#'    \itemize{
#'      \item scaDisp: (scalar) Dispersion parameter estimate. 
#'      \item vecImpulseParam: (numeric vector [beta, h0, h1, h2, t1, t2])
#'        Impulse model parameter estimates.
#'      \item scaLL: (scalar) Loglikelihood of fit.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#' @export

fitDispConstMuImpulseOneInitZINB <- function(scaDispGuess,
  vecImpulseParamGuess,
  vecCounts,
  vecTimepoints, 
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
      par=c(log(scaDispGuess), vecImpulseParamGuess), 
      fn=evalLogLikDispConstMuImpulseZINB_LinPulse_comp, 
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= vecCounts==0,
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
    print(paste0("vecParamGuess ", paste(c(scaDispGuess, vecImpulseParamGuess),collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecTimepoints ", paste(vecTimepoints,collapse=" ")))
    print(paste0("matDropoutLinModel ", paste(matDropoutLinModel,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
    print(paste0("scaWindowRadius", paste(scaWindowRadius,collapse=" ")))
    print(paste0("MAXIT ", MAXIT))
    lsErrorCausingGene <- list(c(scaDispGuess, vecImpulseParamGuess), vecCounts, vecTimepoints, 
      vecNormConst, vecindTimepointAssign, scaWindowRadius, MAXIT)
    names(lsErrorCausingGene) <- c("vecParamGuess", "vecCounts", "vecTimepoints",
      "vecNormConst", "vecindTimepointAssign", "scaWindowRadius", "MAXIT")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispImpulse[1])
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < .Machine$double.eps){scaDisp <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){scaDisp <- 1/.Machine$double.eps}
  
  vecImpulseParam <- fitDispImpulse[2:7]
  scaLL <- fitDispImpulse[8]
  scaConvergence <- fitDispImpulse[9]
  
  return( list(scaDisp=scaDisp,
    vecImpulseParam=vecImpulseParam,
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
#' @param scaDispGuess: (scalar) Negative binomial
#'    dispersion parameter estimate.
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
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (vector number of cells) Posterior probability
#'    of not being drop-out for all observations. Used for impulse
#'    model initialisation (not for likelihood evaluation!).
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
#'      \item scaDisp: (scalar) Dispersion parameter estimate. 
#'      \item vecImpulseParam: (numeric vector [beta, h0, h1, h2, t1, t2])
#'        Impulse model parameter estimates.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#' @export

fitDispConstMuImpulseZINB <- function(vecCounts, 
  scaDispGuess,
  vecImpulseParamGuess,
  lsMuModelGlobal,
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  scaWindowRadius ){
  
  # Set reporter parameters for optim BFGS optimisation
  trace <- 0 # Report typ: 0 none, 2 yes
  REPORT <- 1 # Frequency of reporting in iterations
  # Try new peak and valley initialisations?
  # Increases time complexity of mean estimation by factor 3
  # but seems to make a difference on simulated data.
  boolUseNewInits <- TRUE
  
  # (I) Process data
  # Compute time point specifc parameters
  vecTimepoints <- sort(unique( lsMuModelGlobal$vecPseudotime ))
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
    vecDisp=rep(scaDispGuess, length(vecCounts)),
    vecDrop=vecDropoutParam,
    vecNormConst=vecNormConst,
    scaWindowRadius=scaWindowRadius)
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # (II) Compute new parameters
  # 1. Initialisation: Prior best fit
  lsFitPrior <- fitDispConstMuImpulseOneInitZINB(
    scaDispGuess=scaDispGuess,
    vecImpulseParamGuess=vecImpulseParamGuess,
    vecCounts=vecCounts,
    vecTimepoints=vecTimepoints, 
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
    lsFitPeak <- fitDispConstMuImpulseOneInitZINB(
      scaDispGuess=scaDispGuess,
      vecImpulseParamGuess=vecParamGuessPeak,
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints, 
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=lsMuModelGlobal$MAXIT, 
      RELTOL=lsMuModelGlobal$RELTOL, 
      trace=trace, REPORT=REPORT )
    # 3. Initialisation: Valley
    lsFitValley <- fitDispConstMuImpulseOneInitZINB(
      scaDispGuess=scaDispGuess,
      vecImpulseParamGuess=vecParamGuessValley,
      vecCounts=vecCounts,
      vecTimepoints=vecTimepoints, 
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
  
  scaDisp <- lsFitBest$scaDisp
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
      vecDispEst=rep(scaDispGuess, length(vecCounts)), 
      vecDropoutRateEst=vecDropoutOld,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0, 
      vecboolZero=vecCounts==0 )
    
    # report all new parame
    vecLinModelOut <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValue[cell])) * matDropoutLinModel[cell,])
    })
    vecDropout <- 1/(1+exp(-vecLinModelOut))
    scaLLRef <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecImpulseValue*vecNormConst,
      vecDispEst=rep(scaDisp, length(vecCounts)), 
      vecDropoutRateEst=vecDropout,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0, 
      vecboolZero=vecCounts==0 )
    print(paste0("Old:", scaLLOld, " ,New recomputed: ", 
      scaLLRef, " ,New from optim: ", lsFitBest$scaLL))
  }
  
  return(list(scaDisp=scaDisp,
    vecImpulseParam=vecImpulseParam,
    scaConvergence=scaConvergence))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Top level auxillary function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Coordinate mean and dispersion parameter co-estimation step
#' 
#' Auxillary function that calls the estimation functions for the
#' different mean models according to their needs. Note that one
#' function has to be coded for each combination of mean and dispersion
#' models.
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
#' @return list (length 3)
#'    \itemize{
#'      \item matMuModel: (numeric matrix genes x mu model parameters)
#'        Contains the model parameters according to the used model.
#'      \item matMuModel: (numeric matrix genes x mu model parameters)
#'        Contains the model parameters according to the used model.
#'      \item vecConvergence: (numeric vector number of genes) 
#'        Convergence status of optim for each gene.
#'    }
#' @export

fitZINBMuDisp <- function( matCountsProc,
  vecSizeFactors,
  lsMuModel,
  lsDispModel,
  lsDropModel,
  scaWindowRadius ){
  
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  
  # Nest mean models within dispersion models
  if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){
    if(lsMuModel$lsMuModelGlobal$strMuModel=="windows" & 
        lsMuModel$lsMuModelGlobal$boolVecWindowsAsBFGS){
      # The dispersion parameter co-estimation couples all LL
      # terms of a gene so that there is no computational 
      # advantage in not estimating all mean parameters
      # simultaneously (boolVecWindowsAsBFGS=TRUE).
      lsFitDispMu <- bplapply(seq(1,scaNumGenes), function(i){
        # Decompress parameters
        vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          vecInterval=NULL )
        scaDispParam <- lsDispModel$matDispModel[i,]
        
        fitDispMu <- fitDispConstMuVecWindowsZINB(
          vecCounts=matCountsProc[i,],
          scaDispGuess=scaDispParam,
          vecMu=vecMuParam,
          vecNormConst=vecSizeFactors,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius )
        return(fitDispMu)
      })
      matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$scaDisp))
      matMuModel <- do.call(rbind, lapply(lsFitDispMu, function(i) i$vecMu ))
      vecConvergence <- sapply(lsFitDispMu,  function(i) i$scaConvergence)
      
    } else if(lsMuModel$lsMuModelGlobal$strMuModel=="clusters"){
      # Estimate mean parameter by cluster. No smoothing is used.
      lsFitsDispMu <- bplapply(seq(1,scaNumGenes), function(i){        
        # Estimate mean parameters
        fitDispMu <-fitDispConstMuClusterZINB(
          vecCounts=matCountsProc[i,],
          vecMuGuess=lsMuModel$matMuModel[i,],
          scaDispGuess=lsDispModel$matDispModel[i,],
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          vecNormConst=vecSizeFactors,
          vecindClusterAssign=lsMuModel$lsMuModelGlobal$vecindClusterAssign )
        
        return(fitDispMu)
      })
      matMuModel <- do.call(rbind, lapply(lsFitsDispMu, function(i) i$vecMu ))
      matDispModel <- do.call(rbind, lapply(lsFitsDispMu, function(i) i$scaDisp ))
      vecConvergence <- sapply(lsFitsDispMu,  function(i) i$scaConvergence)
      
    } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){      
      lsFitDispMu <- bplapply(seq(1,scaNumGenes), function(i){ 
        # Decompress parameters
        vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
          lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
          vecInterval=NULL)
        
        # Estimate mean parameters
        fitDispMu <- fitDispConstMuImpulseZINB(
          vecCounts=matCountsProc[i,],
          scaDispGuess=lsDispModel$matDispModel[i,],
          vecImpulseParamGuess=lsMuModel$matMuModel[i,],
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          vecNormConst=vecSizeFactors,
          scaWindowRadius=scaWindowRadius )
        return(fitDispMu)
      })
      matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$scaDisp))
      matMuModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$vecImpulseParam))
      vecConvergence <- sapply(lsFitDispMu,  function(i) i$scaConvergence)
      
    } else if(lsMuModel$lsMuModelGlobal$strMuModel=="constant"){
      lsFitDispMu <- bplapply( seq(1,scaNumGenes), function(i){
        
        # Decompress parameters
        # Don't need to decompress any
        
        # Estimate constant mean parameter
        fitDispMu <- fitDispConstMuConstZINB( vecCounts=matCountsProc[i,],
          scaDispGuess=lsDispModel$matDispModel[i,],
          scaMuGuess=lsMuModel$matMuModel[i,],
          vecNormConst=vecSizeFactors,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius )
        return(fitDispMu)
      })
      matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$scaDisp))
      matMuModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$scaMu))
      vecConvergence <- sapply(lsFitDispMu,  function(i) i$scaConvergence)
    }else {
      #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
      print(paste0("Mean parameter model not recognised for co-estimation with dispersion: ", 
        lsMuModel$lsMuModelGlobal$strMuModel, 
        ". Only constant and impulse model implemented. ",
        "Use sequential estimation."))
      stop(paste0("Mean parameter model not recognised for co-estimation with dispersion: ", 
        lsMuModel$lsMuModelGlobal$strMuModel, 
        ". Only constant and impulse model implemented. ",
        "Use sequential estimation."))
    }  
  } else {
    #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
    print(paste0("Dispersion parameter model not recognised: ", 
      lsDispModel$lsDispModelGlobal$strDispModel, 
      ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
    stop(paste0("Dispersion parameter model not recognised: ", 
      lsDispModel$lsDispModelGlobal$strDispModel, 
      ". Only constant model implemented. Contact david.seb.fischer@gmail.com for alternatives."))
  }
  
  return( list(matMuModel=matMuModel,
    matDispModel=matDispModel,
    vecConvergence=vecConvergence) )
}