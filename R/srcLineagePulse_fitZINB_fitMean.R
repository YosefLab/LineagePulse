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
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(scaMu),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
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
  vecLogMuNeighbourhoodCurrent <- log(vecMu) # Means which are not updated
  vecLogMuNeighbourhoodCurrent[scaTarget] <- log(scaMuNew) # Mean which is updated
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,vecLogMuNeighbourhoodCurrent[cell],vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
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
  vecLinModelOut <- sapply(seq(1,length(vecCounts)), function(cell){
    sum(matDropoutLinModel[cell,] * c(1,log(vecMu[cell]),vecPiConstPredictors))
  })
  vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  
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
#' @param matDropoutLinModel: (matrix number of cells x 2) Logistic linear
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
  matDropoutLinModel=NULL,
  vecProbNB,
  vecNormConst,
  vecboolZero,
  vecboolNotZeroObserved,
  vecPseudotime,
  vecImpulseParam,
  scaWindowRadius=NULL,
  MAXIT=1000,
  RELTOL=sqrt(.Machine$double.eps) ){
  
  # Set reporter parameters for optim BFGS optimisation
  trace <- 0 # Report typp: 0 none, 2 yes
  REPORT <- 1 # Frequency of reporting in iterations
  # Try new peak and valley initialisations?
  # Increases time complexity of mean estimation by factor 3
  # but seems to make a difference on simulated data.
  boolUseNewInits <- TRUE
  
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
  
  # (II) Compute new parameters
  # 1. Initialisation: Prior best fit
  vecFitPrior <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessPrior,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=NA,
    matLinModelPi=matDropoutLinModel,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolZero=vecboolZero, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode="singlecell", 
    MAXIT=MAXIT, RELTOL=RELTOL, 
    trace=trace, REPORT=REPORT )  
  if(boolUseNewInits){
    # 2. Initialisation: Peak
    vecFitPeak <- optimiseImpulseModelFit(
      vecParamGuess=vecParamGuessPeak,
      vecTimepoints=vecTimepoints, 
      vecCounts=vecCounts,
      scaDispersionEstimate=scaDispersionEstimate,
      vecDropoutRate=NA,
      matLinModelPi=matDropoutLinModel,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      vecboolZero=vecboolZero, 
      vecboolNotZeroObserved=vecboolNotZeroObserved,
      scaWindowRadius=scaWindowRadius,
      strMode="singlecell", 
      MAXIT=MAXIT, RELTOL=RELTOL,
      trace=trace, REPORT=REPORT)
    # 3. Initialisation: Valley
    vecFitValley <- optimiseImpulseModelFit(
      vecParamGuess=vecParamGuessValley,
      vecTimepoints=vecTimepoints, 
      vecCounts=vecCounts,
      scaDispersionEstimate=scaDispersionEstimate,
      vecDropoutRate=NA,
      matLinModelPi=matDropoutLinModel,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      vecboolZero=vecboolZero, 
      vecboolNotZeroObserved=vecboolNotZeroObserved,
      scaWindowRadius=scaWindowRadius,
      strMode="singlecell", 
      MAXIT=MAXIT, RELTOL=RELTOL, 
      trace=trace, REPORT=REPORT)
    
    # (IV) Process fits
    dfFitsByInitialisation <- cbind(vecFitPeak, vecFitValley, vecFitPrior)
    
    # Select best fit and report fit type
    # Report mean fit objective value as null hypothesis, too.
    # match() selects first hit if maximum occurs multiple times
    indBestFit <- match(max(dfFitsByInitialisation["value",]),dfFitsByInitialisation["value",])
    
    vecBestFitParam <- c(dfFitsByInitialisation[1:6,indBestFit])
  } else {
    vecBestFitParam <- vecFitPrior[1:6]
  }
  
  # THIS CODE IS ONLY FOR DEVELOPERS TO DEBUG IMPULSE FITTING
  # Follow the choise of model and loglikelihoods
  if(TRUE){
    # Compute predicted means
    vecImpulseValue <- calcImpulse_comp(vecBestFitParam,vecTimepoints)[vecindTimepointAssign]
    vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
    vecImpulseValueOld <- calcImpulse_comp(vecParamGuessPrior,vecTimepoints)[vecindTimepointAssign]
    vecImpulseValueOld[vecImpulseValueOld < .Machine$double.eps] <- .Machine$double.eps
    
    vecLinModelOutOld <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValueOld[cell])) * matDropoutLinModel[cell,])
    })
    vecDropoutOld <- 1/(1+exp(-vecLinModelOutOld))
    scaLLOld <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts,
      vecMu=vecImpulseValueOld*vecNormConst,
      vecDispEst=rep(scaDispersionEstimate, length(vecCounts)), 
      vecDropoutRateEst=vecDropoutOld,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
    print(paste0("Old:", scaLLOld))
    
    # report all new parame
    print(dfFitsByInitialisation["value",]) # switch for below if only prior initialisation is used
    #print(vecFitPrior["value"])
    vecLinModelOut <- sapply(seq(1, length(vecCounts)), function(cell){
      sum(c(1,log(vecImpulseValue[cell])) * matDropoutLinModel[cell,])
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
  
  return(list(vecBestFitParam=vecBestFitParam))
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
  scaMu <- tryCatch({
    exp(unlist(optim(
      par=log(scaMuGuess),
      evalLogLikMuConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
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
  
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
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

fitMuWindowZINB_LinPulse <- function(vecCounts,
  vecMu,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaTarget,
  scaWindowRadius=NULL ){ 
  
  scaMu <- tryCatch({
    exp(unlist(optim(
      par=log(vecMu[scaTarget]),
      evalLogLikMuWindowZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecMu=vecMu,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      scaTarget=scaTarget,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuWindowZINB_LinPulse().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
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

fitMuVecWindowsZINB_LinPulse<- function(vecCounts,
  vecMu,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors,
  scaWindowRadius=NULL ){
  
  vecMu <- tryCatch({
    exp(unlist(optim(
      par=log(vecMu),
      evalLogLikMuVecWindowsZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuVecZINB_LinPulse().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Catch boundary of likelihood domain on mu space
  vecMu[is.na(vecMu) | vecMu < .Machine$double.eps] <- .Machine$double.eps
  
  return(vecMu)
}


#' Fit negative binomial mean parameter to a cluster of cells
#' 
#' Fits negative binomial mean parameter numerically as maximum likelihood estimator
#' to a cell in an interval.
#' 
#' @seealso Called by \code{fitZINBMu}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param scaMuGuess: (scalar) 
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
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuClusterZINB_LinPulse <- function(vecCounts,
  scaMuGuess,
  vecDisp,
  vecNormConst,
  matDropoutLinModel,
  vecPiConstPredictors ){ 
  
  # scaWindowRadius is set to NULL because smoothing
  # within clusters does't make sense - the clusters already impose
  # a constraint on the means.
  scaMu <- tryCatch({
    exp(unlist(optim(    
      par=log(scaMuGuess),
      evalLogLikMuConstZINB_LinPulse_comp,
      vecCounts=vecCounts,
      vecDisp=vecDisp,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      scaWindowRadius=NULL,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1) )["par"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuClusterZINB_LinPulse().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Catch boundary of likelihood domain on mu space
  if(is.na(scaMu) | scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
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
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
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
#' 
#' @return scaMu: (scalar) MLE of mean parameter.
#' @export

fitMuImpulseZINB_LinPulse <- function(vecCounts,
  vecDisp,
  vecNormConst,
  matDropoutLinModel=NULL,
  vecPiConstPredictors=NULL,
  vecProbNB=NULL,
  scaWindowRadius=NULL,
  vecPseudotime=NULL,
  vecImpulseParam=NULL,
  MAXIT=1000,
  RELTOL=sqrt(.Machine$double.eps) ){
  
  lsImpulseFit <- tryCatch({
    fitImpulse_gene_LP(vecCounts=vecCounts, 
      scaDispersionEstimate=vecDisp[1],
      matDropoutLinModel=matDropoutLinModel,
      vecProbNB=vecProbNB,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      vecPseudotime=vecPseudotime,
      vecImpulseParam=vecImpulseParam,
      scaWindowRadius=scaWindowRadius,
      MAXIT=MAXIT, RELTOL=RELTOL )
    
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuImpulseZINB_LinPulse().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecDisp, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "vecDisp","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  return(lsImpulseFit)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Top level auxillary function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Coordinate mean parameter estimation step
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
#' @param vecPseudotime: (numerical vector number of cells)
#'    [Default NULL]
#'    Pseudotime coordinates of cells. Only required if mean model
#'    or dispersion model are fit as a function of pseudotime, 
#'    e.g. impulse model for means.
#' @param vecindClusterAssign: (integer vector length number of
#'    cells) [Default NULL] Index of cluster assigned to each cell.
#' @param matMu: (numeric matrix genes x cells) [Default NULL]
#'    Inferred zero inflated negative binomial mean parameters on 
#'    which estimation in this step  is conditioned on. 
#'    Not required for all models.
#' @param  matDispersions: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial dispersion parameters
#'    on which estimation in this step is conditioned on.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates
#'    on which estimation in this step is conditioned on.
#'    These are the observation-wise point estimates, not the
#'    logistic functions.
#' @param matPiConstPredictors: (numeric matrix genes x external 
#'    predictors) [Default NULL]
#' @param matDropoutLinModel: (numeric matrix genes x predictors)
#'    [Default NULL]
#' @param matImpulseParam: (numeric matrix genes x impulse 
#'    parameters 6) Inferred impulse model parameters
#'    if strMuModel is "impulse". NA for all other strMuModel.
#' @param matboolZero: (numeric matrix genes x cells)
#'    [Default NULL] Whether observation is zero.
#' @param matboolNotZeroObserved: (numeric matrix genes x cells)
#'    [Default NULL] Whether observation is non-zero and real.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param boolVecWindowsAsBFGS: (bool) [Default FALSE] Whether
#'    mean parameters of a gene are co-estimated in "windows"
#'    mode with BFGS algorithm (optimisation with dimensionality
#'    of number of cells) or estimated one by one, conditioned
#'    one the latest estimates of neighbours. The latter case
#'    (boolVecWindowsAsBFGS=FALSE) is coordinate ascent within the gene
#'    and each mean parameter is optimised once only.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param MAXIT_BFGS_Impulse: (scalar) [Default 1000]
#'    Maximum number of BFGS iterations used to estiamte
#'    an impulse model for a gene.
#' @param RELTOL_BFGS_Impulse: (scalar) [Default 
#'    sqrt(.Machine$double.eps)] Minimum relativ decrease
#'    in the loglikelihood objective within one iteration
#'    of BFGS impulse model fitting to gene for the estimation
#'    to not complete as converged. The default is the 
#'    default of optim which is overridden by the value set in
#'    fitZINB by default in LineagePulse to save computation 
#'    time.
#' 
#' @return matMuModel: (numeric matrix genes x mu model parameters)
#'    Contains the model parameters according to the used model.
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
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
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
        vecMu <- fitMuVecWindowsZINB_LinPulse(
          vecCounts=matCountsProc[i,],
          vecMu=vecMuParam,
          vecDisp=vecDispParam,
          vecNormConst=vecSizeFactors,
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius )
      } else {
        vecMu <- vecMuParam
        for(j in seq(1,scaNumCells)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          vecMu[j] <- fitMuWindowZINB_LinPulse(vecCounts=matCountsProc[i,vecInterval],
            vecMu=vecMu[vecInterval],
            vecDisp=vecDispParam[vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
            vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
            scaTarget=match(j,vecInterval),
            scaWindowRadius=scaWindowRadius )
        }
      }
      return(vecMu)
    }))
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="clusters"){
    # Estimate mean parameter by cluster. No smoothing is used.
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate mean parameters
      vecMu <- sapply(seq(1,max(lsMuModel$lsMuModelGlobal$vecindClusterAssign)), function(k){
        vecInterval <- lsMuModel$lsMuModelGlobal$vecindClusterAssign==k
        scaMu <- fitMuClusterZINB_LinPulse(
          vecCounts=matCountsProc[i,vecInterval],
          scaMuGuess=lsMuModel$matMuModel[i,k],
          vecDisp=vecDispParam[vecInterval],
          vecNormConst=vecSizeFactors[vecInterval],
          matDropoutLinModel=lsDropModel$matDropoutLinModel[vecInterval,],
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
        return(scaMu)
      })
      return(vecMu)
    }))

  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){
    
    matMuModel <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
        vecInterval=NULL )
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      vecDropoutParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
        vecMu=vecMuParam,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
      
      #  Compute posterior for parameter
      vecZ <- calcPostDrop_Vector( vecMu=vecMuParam,
        vecDispersions=vecDispParam,
        vecDropout=vecDropoutParam,
        vecboolZero= matCountsProc[i,]==0,
        vecboolNotZeroObserved= !is.na(matCountsProc[i,]) & matCountsProc[i,]>0,
        scaWindowRadius=scaWindowRadius )
      
      # Estimate mean parameters
      lsImpulseFit <- fitMuImpulseZINB_LinPulse(
        vecCounts=matCountsProc[i,],
        vecDisp=vecDispParam,
        vecNormConst=vecSizeFactors,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        vecProbNB=1-vecZ,
        vecPseudotime=lsMuModel$lsMuModelGlobal$vecPseudotime,
        vecImpulseParam=lsMuModel$matMuModel[i,],
        scaWindowRadius=scaWindowRadius,
        MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_Impulse,
        RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_Impulse )
      return(lsImpulseFit$vecBestFitParam)
    }))
    
  } else if(lsMuModel$lsMuModelGlobal$strMuModel=="constant"){
    matMuModel <- do.call(rbind, bplapply( seq(1,scaNumGenes), function(i){
      
      # Decompress parameters
      vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
        lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
        vecInterval=NULL)
      
      # Estimate constant mean parameter
      scaMu <- fitMuConstZINB( vecCounts=matCountsProc[i,],
        scaMuGuess=lsMuModel$matMuModel[i,],
        vecDisp=vecDispParam,
        vecNormConst=vecSizeFactors,
        matDropoutLinModel=lsDropModel$matDropoutLinModel,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
        scaWindowRadius=scaWindowRadius )
      return(scaMu)
    }))
  }
  
  return( matMuModel=matMuModel )
}