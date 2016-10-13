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
#' Calls \code{calcImpulse} and \code{evalLogLikZINB_comp}.
#' 
#' @param vecTheta: (numeric vector dispersion (1) and impulse parameters (6)) 
#'    Dispersion model and impulse model parameter estimates.
#' @param vecCounts: (integer vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
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
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecPseudotime).
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
  vecPseudotime,
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
  # Log linker for amplitudes is in calcImpulse_comp.
  vecImpulseValue <- calcImpulse_comp(vecTheta[2:7],vecPseudotime)[vecindTimepointAssign]
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < .Machine$double.eps){ scaDisp <- .Machine$double.eps }
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/.Machine$double.eps){ scaDisp <- 1/.Machine$double.eps }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # Parameter shrinkage for the impulse model is caught at the model value
  # evaluation, not at the parameteres themselves! Because model valuate
  # can deviate from amplitudes if transition times are switched.
  vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
  
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
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
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
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecPseudotime).
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' @param MAXIT: (scalar) [Default 100] Number of iterations, which are performed 
#'    to fit the impulse model to the clusters.
#'   
#' @return list: (length 2)
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
  vecPseudotime, 
  matDropoutLinModel,
  vecPiConstPredictors,
  vecNormConst,
  vecindTimepointAssign,
  scaWindowRadius=NULL,
  MAXIT=100,
  RELTOL=sqrt(.Machine$double.eps),
  trace=0,
  REPORT=10 ){
  
  fitDispImpulse <- tryCatch({
    unlist( optim(
      par=c(log(scaDispGuess), vecImpulseParamGuess), 
      fn=evalLogLikDispConstMuImpulseZINB_LinPulse_comp, 
      vecCounts=vecCounts,
      vecPseudotime=vecPseudotime,
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
    print(paste0("vecPseudotime ", paste(vecPseudotime,collapse=" ")))
    print(paste0("matDropoutLinModel ", paste(matDropoutLinModel,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
    print(paste0("scaWindowRadius", paste(scaWindowRadius,collapse=" ")))
    print(paste0("MAXIT ", MAXIT))
    lsErrorCausingGene <- list(c(scaDispGuess, vecImpulseParamGuess), vecCounts, vecPseudotime, 
      vecNormConst, vecindTimepointAssign, scaWindowRadius, MAXIT)
    names(lsErrorCausingGene) <- c("vecParamGuess", "vecCounts", "vecPseudotime",
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
#' @param scaDispGuess: (numerical scalar) Negative binomial
#'    dispersion parameter estimate.
#' @param vecImpulseParamGuess: (numerical vector impulse parameters)
#'    Previous impulse model parameter values.
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
  matDropoutLinModel,
  vecPiConstPredictors,
  vecProbNB,
  vecNormConst,
  vecPseudotime,
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
  
  # (I) Process data
  # Compute time point specifc parameters
  vecTimepoints <- sort(unique( vecPseudotime ))
  # Get vector of numeric time point assignment indices:
  vecindTimepointAssign <- match(vecPseudotime, vecTimepoints)
  
  # (II) Fit Impulse model
  # The previous parameter estiamte is kept as a reference and
  # used as an initialisation
  # Compute initialisations for peak and valley
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
  lsFitPrior <- fitDispConstMuImpulseOneInitZINB(
    scaDispGuess=scaDispGuess,
    vecImpulseParamGuess=vecImpulseParamGuess,
    vecCounts=vecCounts,
    vecPseudotime=vecPseudotime, 
    matDropoutLinModel=matDropoutLinModel,
    vecPiConstPredictors=vecPiConstPredictors,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    scaWindowRadius=scaWindowRadius,
    MAXIT=MAXIT, RELTOL=RELTOL, 
    trace=trace, REPORT=REPORT )  
  if(boolUseNewInits){
    # 2. Initialisation: Peak
    lsFitPeak <- fitDispConstMuImpulseOneInitZINB(
      scaDispGuess=scaDispGuess,
      vecImpulseParamGuess=vecParamGuessPeak,
      vecCounts=vecCounts,
      vecPseudotime=vecPseudotime, 
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=MAXIT, RELTOL=RELTOL, 
      trace=trace, REPORT=REPORT )
    # 3. Initialisation: Valley
    lsFitValley <- fitDispConstMuImpulseOneInitZINB(
      scaDispGuess=scaDispGuess,
      vecImpulseParamGuess=vecParamGuessValley,
      vecCounts=vecCounts,
      vecPseudotime=vecPseudotime, 
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecNormConst=vecNormConst,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=MAXIT, RELTOL=RELTOL, 
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
  # Follow the choise of model and loglikelihoods
  if(FALSE){
    # Compute predicted means
    vecImpulseValue <- calcImpulse_comp(vecImpulseParam,vecTimepoints)[vecindTimepointAssign]
    vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
    vecImpulseValueOld <- calcImpulse_comp(vecImpulseParamGuess,vecTimepoints)[vecindTimepointAssign]
    vecImpulseValueOld[vecImpulseValueOld < .Machine$double.eps] <- .Machine$double.eps
    
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
    if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){      
      lsFitDispMu <- bplapply(seq(1,scaNumGenes), function(i){ 
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
        fitDispMu <- fitDispConstMuImpulseZINB(
          vecCounts=matCountsProc[i,],
          scaDispGuess=lsDispModel$matDispModel[i,],
          vecImpulseParamGuess=lsMuModel$matMuModel[i,],
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          vecProbNB=1-vecZ,
          vecNormConst=vecSizeFactors,
          vecPseudotime=lsMuModel$lsMuModelGlobal$vecPseudotime,
          scaWindowRadius=scaWindowRadius,
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_Impulse,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_Impulse )
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