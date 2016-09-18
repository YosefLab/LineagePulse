#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Fit mean parameters of ZINB model    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecMu: (vector number of cells in neighbourhood) Mean
#'    parameter estimates in neighourhood of target cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
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
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuWindowsZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecMu,
  vecDropoutRateEst=NULL,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaTarget,
  scaWindowRadius,
  boolDynamicPi=FALSE){ 
  
  # Log linker function to fit positive means
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  scaN <- length(vecCounts)
  # Correct drop-out rate estimate for mean estimate
  # Restimate drop-out rate based on new mean parameter
  # which is a predictor of the drop-out rate.
  # Reestimate one drop-out rate per cell of neighbourhood
  if(boolDynamicPi){
    matPredictorsPi <- matrix(vecPredictorsPi, nrow=scaN, 
      ncol= length(vecPredictorsPi), byrow=TRUE)
    matPredictorsPi[,2] <- log(vecMu)
    matPredictorsPi[scaTarget,2] <- scaTheta
    vecLinModelOut <- sapply(seq(1,scaN), function(cell){
      matLinModelPi[cell,] %*% matPredictorsPi[cell,]
    })
    vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  }
  
  # Compute loglikelihood terms with contribution from new mean parameter
  # Cis terms: From neighbourhood around target parameter
  scaLogLikCis <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
    vecMu=scaMu*vecNormConst,
    vecDispEst=rep(vecDisp[scaTarget],scaN), 
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero )
  # Trans terms: From other neighbourhoods via drop-out rate as function
  # of target mean parameter.
  if(scaN>1 & boolDynamicPi){
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
#' @param vecTheta: (numeric vector number of means to be estimated) 
#'    Log of mean parameter estimates.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
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
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuVecWindowsZINB_LinPulse <- function(vecTheta,
  vecCounts,
  vecDisp,
  vecDropoutRateEst=NULL,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius,
  boolDynamicPi=FALSE){ 
  
  # Log linker function to fit positive means
  vecMu <- exp(vecTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  scaN <- length(vecCounts)
  # Correct drop-out rate estimate for mean estimate
  # Restimate drop-out rate based on new mean parameter
  # which is a predictor of the drop-out rate.
  # Reestimate one drop-out rate per cell of neighbourhood
  if(boolDynamicPi){
    matPredictorsPi <- matrix(vecPredictorsPi, nrow=scaN, 
      ncol= length(vecPredictorsPi), byrow=TRUE)
    matPredictorsPi[,2] <- log(vecMu)
    vecLinModelOut <- sapply(seq(1,scaN), function(cell){
      matLinModelPi[cell,] %*% matPredictorsPi[cell,]
    })
    vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  }
  
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

#' Cost function zero-inflated negative binomial model for mean fitting
#' under cluster mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter.
#' The mean is modelled by cell and constrained to fit a sliding 
#' window of cells. This function is used if the means for 
#' observations of a gene are fit by cluster.
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
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
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
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuClustersZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecDropoutRateEst=NULL,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero,
  boolDynamicPi=FALSE){ 
  
  # Log linker function to fit positive means
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  # Correct drop-out rate estimate for mean estimate
  # This is not used in the current framework
  # Restimate drop-out rate based on new mean parameter
  # which is a predictor of the drop-out rate.
  # Reestimate one drop-out rate per cell of neighbourhood
  if(boolDynamicPi){
    vecPredictorsPi[2] <- scaTheta
    vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
      sum(matLinModelPi[cell,] * vecPredictorsPi)
    })
    vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  }

  scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
    vecMu=scaMu*vecNormConst,
    vecDispEst=vecDisp, 
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter.
#' The mean is modelled by cell and constrained to fit a sliding 
#' window of cells. This function is used if a constant mean is 
#' assumed.
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
#' @param vecY: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDispEst: (scalar vector number of samples) 
#'    Negative binomial dispersion  parameter for given 
#'    gene and observations.
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
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuConstZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  # Log linker function to fit positive means
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  # Adjust drop-out rates if pi is dynamic in mu
  if(!is.null(matLinModelPi)){
    vecPredictors <- c(1,log(scaMu))
    vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
      sum(matLinModelPi[cell,] * vecPredictors)
    })
    vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  }
  
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

#' Fit impulse model as mean parameter estimation
#' 
#' Computes impulse parameter initialisation for valley
#' and peak model and uses both and the prior parameter fit
#' in three separate optimisation runs to obtain the best 
#' impulse model fit to the data. To avoid convergence and
#' numerical errors, the best fit is compared against the prior
#' fit and the better one kept. Low impulse model values
#' are masked.
#' This function contains parts of \code{fitImpulse_gene} from
#' ImpulseDE2.
#' 
#' @seealso Called by \code{fitZINB}. Code similar to \code{
#' ImpulseDE2::fitImpulse_gene}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param scaDispersionEstimate: (numerical scalar) Negative binomial
#'    dispersion parameter estimate.
#' @param vecDropoutRate: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x 2) Logistic linear
#'    model parameters of the dropout rate as a function of the mean.
#' @param vecProbNB: (vector number of cells) Posterior probability
#'    of not being drop-out for all observations. Used for impulse
#'    model initialisation (not for likelihood evaluation!).
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not NA and has nonzero count.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
#' @param vecImpulseParam: (numerical vector impulse parameters)
#'    Previous impulse model parameter values.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return \enumerate{ 
#'  \item vecBestFitParam: (numerical vector impulse parameters)
#'  Parameter values of best impulse model fit found.
#'  \item vecImpulseValue: (scalar) Loglikelihood of best 
#'  impulse model fit
#'  }
#' @export

fitImpulse_gene_LP <- function(vecCounts, 
  scaDispersionEstimate,
  vecDropoutRate=NULL,
  matLinModelPi=NULL,
  vecProbNB,
  vecNormConst,
  vecboolZero,
  vecboolNotZeroObserved,
  vecPseudotime,
  vecImpulseParam,
  scaWindowRadius=NULL){
  
  # (I) Process data
  # Compute time point specifc parameters
  vecTimepoints <- sort(unique( vecPseudotime ))
  # Get vector of numeric time point assignment indices:
  vecindTimepointAssign <- match(vecPseudotime, vecTimepoints)
  
  # (II) Fit Impulse model
  # The previous parameter estiamte is kept as a reference and
  # used as an initialisation
  vecParamGuessPrior <- vecImpulseParam
  # Compute initialisations
  lsParamGuesses <- estimateImpulseParam(
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    vecDropoutRate=vecDropoutRate, 
    vecProbNB=vecProbNB,
    vecTimepointAssign=vecPseudotime, 
    vecNormConst=vecNormConst,
    strMode="singlecell",
    strSCMode="continuous" )
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # (I) Previous parameters
  scaLLOld <- evalLogLikImpulseSC_comp(
    vecTheta=vecParamGuessPrior,
    vecX=vecTimepoints,
    vecY=vecCounts, 
    scaDispEst=scaDispersionEstimate,
    vecDropoutRateEst=NA,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero,
    scaWindowRadius=scaWindowRadius )
    
  # (II) Compute new parameters
  # 1. Initialisation: Prior best fit
  vecFitPrior <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessPrior,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=NA,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolZero=vecboolZero, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode="singlecell", 
    MAXIT=1000)  
  # 2. Initialisation: Peak
  vecFitPeak <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessPeak,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=NA,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolZero=vecboolZero, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode="singlecell", 
    MAXIT=1000)
  # 3. Initialisation: Valley
  vecFitValley <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessValley,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=NA,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolZero=vecboolZero, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode="singlecell", 
    MAXIT=1000)
  
  # (IV) Process fits
  dfFitsByInitialisation <- cbind(vecFitPeak, vecFitValley, vecFitPrior)
  
  # Select best fit and report fit type
  # Report mean fit objective value as null hypothesis, too.
  # match() selects first hit if maximum occurs multiple times
  indBestFit <- match(max(dfFitsByInitialisation["value",]),dfFitsByInitialisation["value",])
  
  if(scaLLOld < dfFitsByInitialisation["value",indBestFit]){
    # Select best of the new models if achieve higher loglikelihood than
    # previous model
    vecBestFitParam <- c(dfFitsByInitialisation[1:6,indBestFit])
  } else {
    # Keep old parameters to guarantee convergence
    vecBestFitParam <- vecParamGuessPrior
  }
  
  # Compute predicted means
  vecImpulseValue <- calcImpulse_comp(vecBestFitParam,vecTimepoints)[vecindTimepointAssign]
  vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps

  if(FALSE){
    # Follow choice of model
    # Valley is chosen very often (?) and reported likelihoods are wrong for valley
    # Was reporting observed values and not zeros for valley initialisation
    # to optimisation... routine -> seems to have solved the problem
    print(dfFitsByInitialisation["value",])
    print(paste0(scaLLOld, " ", dfFitsByInitialisation["value",indBestFit]))
    vecLinModelOut <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValue[cell])) * matLinModelPi[cell,])
    })
    vecDropout <- 1/(1+exp(-vecLinModelOut))
    scaLLRef <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecImpulseValue*vecNormConst,
      vecDispEst=rep(scaDispersionEstimate, length(vecCounts)), 
      vecDropoutRateEst=vecDropout,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
    print(scaLLRef)
  }
  
  return(list(vecBestFitParam=vecBestFitParam,
    vecImpulseValue=vecImpulseValue))
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
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (numeric vector number of cells) Posterior
#'    of observation not being drop-out for closed-form MLE.
#' @param scaTarget: (integer) Position of target cell
#'    (whose mean is to be estimated) within given interval.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuZINB_LinPulse <- function(vecCounts,
  vecMu=NULL,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecProbNB=NULL,
  scaTarget=NULL,
  scaWindowRadius=NULL,
  strMuModel="windows",
  boolDynamicPi=FALSE ){ 
  
  if(all(vecNormConst==1) & !boolDynamicPi){
    # Closed form maximum likelihood estimator of observed, not joint likelihood
    scaMu <- sum(vecCounts*vecProbNB, na.rm=TRUE)/sum(vecProbNB, na.rm=TRUE)
  } else {
    # Numerical maximum likelihood estimator
    scaMu <- tryCatch({
      if(strMuModel=="windows"){
        exp(unlist(optim(
          par=log(vecMu[scaTarget]),
          evalLogLikMuWindowsZINB_LinPulse_comp,
          vecCounts=vecCounts,
          vecMu=vecMu,
          vecDisp=vecDisp,
          vecDropoutRateEst=vecDropoutRateEst,
          matLinModelPi=matLinModelPi,
          vecPredictorsPi=vecPredictorsPi,
          vecNormConst=vecNormConst,
          vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
          vecboolZero=vecCounts==0,
          scaTarget=scaTarget,
          scaWindowRadius=scaWindowRadius,
          boolDynamicPi=boolDynamicPi,
          method="BFGS",
          control=list(maxit=1000,fnscale=-1) )["par"]))
      } else if(strMuModel=="clusters"){
        exp(unlist(optim(
          par=log(vecMu[1]),
          evalLogLikMuClustersZINB_LinPulse_comp,
          vecCounts=vecCounts,
          vecDisp=vecDisp,
          vecDropoutRateEst=vecDropoutRateEst,
          matLinModelPi=matLinModelPi,
          vecPredictorsPi=vecPredictorsPi,
          vecProbNB=vecProbNB,
          vecNormConst=vecNormConst,
          vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
          vecboolZero=vecCounts==0,
          boolDynamicPi=boolDynamicPi,
          method="BFGS",
          control=list(maxit=1000,fnscale=-1) )["par"]))
      } else {
        stop(paste0("Unrecognised strMuModel in fitMuZINB_LinPulse: ", strMuModel))
      }
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuZINB_LinPulse().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
      print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
      print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
      lsErrorCausingGene <- list(vecCounts, vecDisp, vecDropoutRateEst, vecNormConst)
      names(lsErrorCausingGene) <- c("vecCounts", "vecDisp", "vecDropoutRateEst","vecNormConst")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      print(strErrorMsg)
      stop(strErrorMsg)
    })
  }
  
  # Catch boundary of likelihood domain on mu space
  if(is.na(scaMu) | scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
}

#' Fit negative binomial means to all cells in an interval (BFGS)
#' 
#' Fits negative binomial mean parameters numerically as maximum likelihood estimator
#' to all cells in an interval simultaneously using BFGS.
#' 
#' 
#' @seealso Called by \code{fitZINB}. Alternative to cell-wise sequential
#' maximum likelihood estimation used in \code{fitMuZINB_LinPulse}.
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
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (numeric vector number of cells) Posterior
#'    of observation not being drop-out for closed-form MLE.
#' @param scaTarget: (integer) Position of target cell
#'    (whose mean is to be estimated) within given interval.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuVecZINB_LinPulse<- function(vecCounts,
  vecMu,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  scaWindowRadius=NULL,
  boolDynamicPi=FALSE ){
  
  vecMu <- tryCatch({
    exp(unlist(optim(
      par=log(vecMu),
      evalLogLikMuVecWindowsZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      vecDropoutRateEst=vecDropoutRateEst,
      matLinModelPi=matLinModelPi,
      vecPredictorsPi=vecPredictorsPi,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      boolDynamicPi=boolDynamicPi,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuVecZINB_LinPulse().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecDropoutRateEst, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp", "vecDropoutRateEst","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Catch boundary of likelihood domain on mu space
  vecMu[is.na(vecMu) | vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  return(vecMu)
}

#' Fit negative binomial means as impulse model 
#' to all cells in an interval
#' 
#' Fits negative binomial mean parameters numerically as maximum likelihood estimator
#' constrained to follow the impulse model to all cells in an interval simultaneously.
#' 
#' 
#' @seealso Called by \code{fitZINB}. Alternative to cell-wise sequential
#' maximum likelihood estimation used in \code{fitMuZINB_LinPulse}.
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
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matLinModelPi: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPredictorsPi: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecProbNB: (numeric vector number of cells) Posterior
#'    of observation not being drop-out for closed-form MLE.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
#' @param vecImpulseParam: (numerical vector impulse parameters)
#'    Previous impulse model parameter values.
#' @param boolDynamicPi: (bool) Whether drop-out rate is
#'    treated as function of mean or as a constant.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuImpulseZINB_LinPulse <- function(vecCounts,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecProbNB=NULL,
  scaWindowRadius=NULL,
  vecPseudotime=NULL,
  vecImpulseParam=NULL,
  boolDynamicPi=FALSE ){
  
  lsImpulseFit <- tryCatch({
    fitImpulse_gene_LP(vecCounts=vecCounts, 
      scaDispersionEstimate=vecDisp[1],
      vecDropoutRate=vecDropoutRateEst,
      matLinModelPi=matLinModelPi,
      vecProbNB=vecProbNB,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      vecPseudotime=vecPseudotime,
      vecImpulseParam=vecImpulseParam,
      scaWindowRadius=scaWindowRadius)
    
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuImpulseZINB_LinPulse().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecDropoutRateEst, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp", "vecDropoutRateEst","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  return(lsImpulseFit)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Fits single negative binomial mean parameters numerically as 
#' maximum likelihood estimator to a gene: Constant mean model.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaDispEst: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitMuConstZINB <- function(vecCounts,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  scaWindowRadius){ 

  # Numerical maximum likelihood estimator
  scaMeanGuess <- log(mean(vecCounts, na.rm=TRUE)+1)
  scaMu <- tryCatch({
    exp(unlist(optim(
      par=scaMeanGuess,
      evalLogLikMuConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      vecDropoutRateEst=vecDropoutRateEst,
      matLinModelPi=matLinModelPi,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuZINB().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("scaMeanGuess ",scaMeanGuess))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
    print("matLinModelPi")
    print(matLinModelPi)
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecDropoutRateEst, 
      matLinModelPi, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp", "vecDropoutRateEst",
      "matLinModelPi", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
}