#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Note: Most problems are numerical and not conceptual problems.

# Previously encountered problems:
# Genes with many zeros and few very high counts: If dispersion is initialised
# too high (=1), the increase in LL cannot always be picked up by dnbinom which
# disables convergence.
# Still didnt see convergence, this might reoccur during fitting. Removed -Inf
# masking from likelihood function for NB part. Seeing convergence of means now.
# Still problems with drop-out rate -> Use more genes to get more reliable 
# traces, still problems. Take out masking of dropout contribution of nonzeros
# in likelihood and force logistic to be decreasing in means - good cnvergence
# until 3rd/4th iter.
# Try replace mean optimise with optim and put masking in likelihood again: Note
# that optim becomes unstable if Inf is returned by objective, while optimise
# seems to be forced out of the "bad" parameter region in which dnbinom yields
# LL=-Inf, optim needs the low value masking to not break down. Optim can be
# initialised outside of the bad parameter region to avoid it. Now only 
# problems in mean estimation (3rd iteration). Remove bounds from BFGS
# optimisation. Cleared problem by setting mean initialisation to log(max+1)
# instead of log(mean+1), was it stuck in bad parameter region again? maybe.
# Initialise mean and drop-out rate to prior value, converging fully!
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of logistic drop-out paramater model on single gene given
#' the negative binomial mean and dispersion parameters.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecTheta: (numeric vector length linear model) Linear model
#'    for drop-out rate in logit space.
#' @param matPredictorsPi: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecCounts: (numeric vector number of genes) Observed expression values 
#'    of gene in target cell.
#' @param vecMu: (vector number of cells) Negative binomial
#'    mean parameter estimate.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param scaNormConst: (scalar) 
#'    Model scaling factors for cell which takes
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikPiZINB_LinPulse <- function(vecTheta,
  matPredictorsPi,
  vecCounts,
  matMu,
  matDisp,
  vecNormConst,
  vecboolNotZeroObserved,
  vecboolZero ){ 
  
  # Estimate drop-out rate parameters according to proposed
  # linear model:
  # Force link to mean to be negative
  vecTheta[2] <- -exp(vecTheta[2])
  vecLinModelOut <- matPredictorsPi %*% vecTheta
  vecDropoutRateFit <- 1/(1+exp(-vecLinModelOut))
  
  # Loglikelihood is evaluated on each window which was has
  # target cell in its neighbourhood. Note that the negative binomial
  # parameters therefore change for each iteration as these represent
  # SEPARATE windows, and not a smoothing of the target cell to its 
  # neighbours.
  scaLogLik <- sum(sapply(seq(1,length(vecNormConst)), function(cell){
    evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=matMu[,cell]*vecNormConst[cell],
      vecDispEst=matDisp[,cell], 
      vecDropoutRateEst=vecDropoutRateFit,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  }))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Optimisation function for drop-out model fitting
#' 
#' This function fits a logistic drop-out model to a cell based
#' on given gene-specific predictors (which enter the linear model).
#' Parameter estimation of the linear model is performed by maximum
#' likelihood based on the overall likelihood.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param matPredictorsPi: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecCounts: (numeric vector number of genes) Observed expression values 
#'    of gene in target cell.
#' @param matMu: (matrix genes x neighbourhood) Negative binomial
#'    mean parameter estimate.
#' @param matDisp: (matrix genes x neighbourhood) Negative binomial
#'    dispersion parameter estimate.
#' @param scaNormConst: (numeric vector number of cells in neighbourhood) 
#'    Model scaling factors for cell which takes
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return vecLinModel: (numeric vector length linear model) 
#'    Linear model for drop-out rate in logit space for given cell.
#' @export

fitPiZINB_LinPulse <- function( vecLinModelPi,
  matPredictorsPi,
  vecCounts,
  matMu,
  matDisp,
  vecNormConst ){ 
  
  # Catch matrix format error
  if(length(vecNormConst)==1){
    matMu <- matrix(matMu,nrow=length(matMu),ncol=1, byrow=FALSE)
    matDisp <- matrix(matDisp,nrow=length(matDisp),ncol=1, byrow=FALSE)
  }
  # Numerical maximum likelihood estimaton of linear model
  # Initialise optimisation of linear model:
  #vecParamGuess <- rep(1, dim(matPredictorsPi)[2])
  vecParamGuess <- vecLinModelPi
  vecParamGuess[2] <- log(-vecParamGuess[2])
  lsLinModelFit <- tryCatch({
    optim(
      par=vecParamGuess,
      evalLogLikPiZINB_LinPulse_comp,
      matPredictorsPi=matPredictorsPi,
      vecCounts=vecCounts,
      matMu=matMu,
      matDisp=matDisp,
      vecNormConst=vecNormConst,
      vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
      vecboolZero=vecCounts==0,
      method="BFGS",
      control=list(maxit=1000,fnscale=-1)
    )[c("par","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting logistic drop-out model: fitPiZINB_LinPulse().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("matDisp ", paste(matDisp,collapse=" ")))
    print(paste0("matMu ", paste(matMu,collapse=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, matDisp, matMu, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "matDisp", "matMu","vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  vecLinModel <- unlist(lsLinModelFit["par"])
  vecLinModel[2] <- -exp(vecLinModel[2])
  scaConvergence <- unlist(lsLinModelFit["convergence"])
  
  return(vecLinModel)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
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
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matDropoutLinMod: (matrix number of cells x 2) Logistic linear
#'    model parameters of the dropout rate as a function of the mean.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuZINB_LinPulse <- function(scaTheta,
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
  # This is not used in the current framework
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
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' Note that this fitting routine is meant for fitting a single mean
#' to a set of observations and can be used in wrappers which compute
#' every mean separately (assuming independence).
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitMuZINB_LinPulse <- function(vecCounts,
  vecMu=NULL,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecPredictorsPi=NULL,
  vecProbNB,
  scaTarget=NULL,
  scaWindowRadius,
  boolDynamicPi=FALSE ){ 
  
  if(all(vecNormConst==1) & !boolDynamicPi){
    # Closed form maximum likelihood estimator of observed, not joint likelihood
    scaMu <- sum(vecCounts*vecProbNB, na.rm=TRUE)/sum(vecProbNB, na.rm=TRUE)
  } else {
    # Numerical maximum likelihood estimator
    scaMu <- tryCatch({
      exp(unlist(optim(
        par=log(vecMu[scaTarget]),
        evalLogLikMuZINB_LinPulse_comp,
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
        #lower = -36,
        #upper = 16,
        method="BFGS",
        control=list(maxit=1000,fnscale=-1) )["par"]))
        #maximum = TRUE)["maximum"]))
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
      stop(strErrorMsg)
    })
  }
  
  # Catch boundary of likelihood domain on mu space
  if(is.na(scaMu) | scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
}

#' Cost function zero-inflated negative binomial model for dispersion fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial overdispersion on single gene given
#' the drop-out rate and negative binomial mean parameter. The
#' dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
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
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
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
  scaWindowRadius=NULL){ 
  
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
  
  if(is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=vecMuEst*vecSizeFactors,
      vecDispEst=vecDispersions, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    scaLogLik <- evalLogLikSmoothZINB_LinPulse( vecCounts=vecCounts,
      vecMu=vecMuEst,
      vecSizeFactors=vecSizeFactors,
      vecDispEst=vecDispersions,
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius )
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
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
  scaWindowRadius=NULL){ 
  
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
  
  # Catch boundary of likelihood domain on mu space
  if(fitDisp["par"] < .Machine$double.eps){fitDisp["par"] <- .Machine$double.eps}
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(fitDisp["par"] > 1/.Machine$double.eps){fitDisp["par"] <- 1/.Machine$double.eps}
  
  return(fitDisp)
}

#' Fit zero-inflated negative binomial model to data
#' 
#' Fit zero-inflated negative binomial model to data: One mean per cluster 
#' and either one dispersion parameter across all observations in all cluster
#' for a gene or one dispersion parameter per cluster per gene. Dropout rate
#' and dispersion factor inferred here are used as hyperparamters in the
#' impulse fitting stage of PseudoDE.
#' To parallelise this code, replace lapply by bplapply in dispersion
#' factor and drop-out rate estimation and uncomment the BiocParallel
#' register() command at the beginning of this function.
#' Counts can be imputed under the ZINB model as:
#' matCountsProcImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecSpikeInGenes: (string vector) Names of genes
#'    which correspond to external RNA spike-ins. Currently
#'    not used.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' @param boolOneDispPerGene: (bool) [Default TRUE]
#'    Whether one negative binomial dispersion factor is fitted
#'    per gene or per gene for each cluster.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' @param boolSuperVerbose: (bool) Whether to follow EM-algorithm
#'    progress in high detail with local convergence flags. 
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' 
#' @return (list length 6)
#'    \itemize{
#'      \item vecDispersions: (numeric matrix genes x clusters)
#'        Inferred negative binomial dispersions.
#'      \item matDropout: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial drop out rates.
#'      \item matProbNB: (numeric matrix genes x cells)
#'        Inferred probabilities of every observation to be generated
#'        from the negative binomial component of the zero-inflated
#'        negative binomial mixture model.
#'      \item matCountsProcImputed: (numeric matrix genes x cells)
#'        Data predicted with inferred zero-inflated negative binomial
#'        model.
#'      \item matMuCluster: (numeric matrix genes x clusters)
#'        Inferred negative binomial cluster means.
#'      \item boolConvergence: (bool) Convergence of EM algorithm.
#'    }
#' @export

fitZINB <- function(matCountsProc, 
  lsResultsClustering,
  vecSizeFactors,
  vecSpikeInGenes=NULL,
  boolOneDispPerGene=TRUE,
  scaWindowRadius=NULL,
  scaMaxiterEM=100,
  verbose=FALSE,
  boolSuperVerbose=FALSE,
  nProc=1){
  
  # Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-6)
  
  # Store EM convergence
  vecEMLogLik <- array(NA,scaMaxiterEM)
  scaPredictors <- 2
  boolDynamicPi <- TRUE
  matConstPredictorsPi <- NULL
  
  # Set number of processes to be used for parallelisation
  # This function is currently not parallelised to reduce memory usage.
  # Read function description for further information.
  register(MulticoreParam(nProc))
  
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]  
  
  # (I) Initialisation
  # The initial model is estimated under the assumption that all zero-counts
  # are drop-outs.
  
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
  matboolZero <- matCountsProc == 0
  
  vecMu <- apply(matCountsProc, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
  matMu <- matrix(vecMu, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  matDispersions <- matrix(0.001, nrow=scaNumGenes, ncol=scaNumCells)
  matDropout <- matrix(0.5, nrow=scaNumGenes, ncol=scaNumCells)
  #matDropout[matboolNotZeroObserved] <- 0
  matLinModelPi <- cbind(rep(0, scaNumCells), rep(-10^(-10), scaNumCells), 
    matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
  if(boolSuperVerbose){
    scaLogLikInit <- sum(unlist(
      bplapply( seq(1,scaNumGenes), function(i){
        evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
          vecMu=matMu[i,],
          vecSizeFactors=vecSizeFactors,
          vecDispEst=matDispersions[i,], 
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius)
      })
    ))
    print(paste0("Completed initialisation with log likelihood of ", scaLogLikInit))
  }
  
  # (II) EM itertion
  scaIter <- 1
  scaLogLikNew <- 0
  scaLogLikOld <- 0
  
  while(scaIter == 1 | scaIter == 2 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){
    ##### 1. Gene-wise parameter estimation: 
    # a) Negative binomial mean parameter
    # Only compute posterior if using closed form estimator for mean:
    # Posterior is not necessary in all other cases!
    if(all(vecSizeFactors==1) & !boolDynamicPi){    
      matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
        (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
      }))
      matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
        vecZ <- sapply(seq(1,scaNumCells), function(j){
          if(matboolNotZeroObserved[i,j]){
            scaZ <- 0
          } else {
            scaindIntervalStart <- max(1,j-scaWindowRadius)
            scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
            vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
            scaZ <- sum(matDropout[i,j]/(matDropout[i,j] + 
                (1-matDropout[i,j])*matNBZero[i,vecInterval])) *
              1/length(vecInterval)
          }
          return(scaZ)
        })
        return(vecZ)
      }))
    }
    if(!is.null(scaWindowRadius)){
      # Estimate mean parameter for each cell as ZINB model for cells within pseudotime
      # interval with cell density centred at target cell.
      # Note that this corresponds to maximising smoothed log likelihood but instead
      # of using the implemented cost function evalLogLikSmoothZINB_LinPulse for an entire
      # gene, the optimisation problem is broken up into 1D problems for each mean.
      matMu <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
        # Note: estimation if confounded, always use latest update
        # of means when estimating new parameters
        vecMu <- matMu[i,]
        for(j in seq(1,scaNumCells)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          vecMu[j] <- fitMuZINB_LinPulse(vecCounts=matCountsProc[i,vecInterval],
            vecMu=vecMu[vecInterval],
            vecDisp=matDispersions[i,vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            vecDropoutRateEst=matDropout[i,vecInterval],
            vecProbNB=1-matZ[i,vecInterval],
            vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,]),
            matLinModelPi=matLinModelPi[vecInterval,],
            scaTarget=match(j,vecInterval),
            scaWindowRadius=scaWindowRadius,
            boolDynamicPi=boolDynamicPi )
          # vecDropout is re-estimated in mean estimation based on the new mean.
        }
        return(vecMu)
      }))
    } else {
      # Estimate mean parameter by cluster. No smoothing is used.
      matMuCluster <- do.call(cbind, lapply(seq(1,max(vecindClusterAssign)), function(k){
        vecMu <- unlist(bplapply(seq(1,scaNumGenes), function(i){
          scaMu <- fitMuZINB_LinPulse(vecCounts=matCountsProc[i,vecindClusterAssign==k],
            vecDisp=matDispersions[i,vecindClusterAssign==k],
            vecNormConst=vecSizeFactors[vecindClusterAssign==k],
            vecDropoutRateEst=matDropout[i,vecindClusterAssign==k],
            vecProbNB=1-matZ[i,vecindClusterAssign==k])
          return(scaMu)
        }))
        return(vecMu)
      }))
      matMu <- matMuCluster[,vecindClusterAssign]
    }
    # These udates are done during mean estimation too
    # Reestimate drop-out rates based on new means
    matDropout <- do.call(cbind, lapply(seq(1, scaNumCells), function(cell){
      vecLinModelOut <- cbind(1,log(matMu[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
      vecDropout <- 1/(1+exp(-vecLinModelOut))
      return(vecDropout)
    }))
    if(boolSuperVerbose){
      scaLogLikTemp <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){          
          evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecSizeFactors=vecSizeFactors,
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,],
            scaWindowRadius=scaWindowRadius)
        })
      ))
      print(paste0("1a) Mean estimation complete: loglikelihood of ", scaLogLikTemp))
    }
    
    # b) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood.
    if(boolOneDispPerGene){  
      vecDispFit <- bplapply(seq(1,scaNumGenes), function(i){
        fitDispZINB_LinPulse(scaDispGuess=matDispersions[i,1],
          vecCounts=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius )
      })
    } else {
      #  not coded
      stop("not coded")
    }
    vecboolConvergedGLMdisp <- sapply(vecDispFit, function(fit) fit["convergence"])
    vecDispersions <- sapply(vecDispFit, function(fit){ fit["par"] })
    matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    if(boolSuperVerbose){
      print(paste0("Dispersion estimation did not converge in ", 
        sum(vecboolConvergedGLMdisp), " cases."))
      scaLogLikTemp <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecSizeFactors=vecSizeFactors,
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,],
            scaWindowRadius=scaWindowRadius)
        })
      ))
      print(paste0("1b) Dispersion estimation complete: loglikelihood of ", scaLogLikTemp))
    }  
    
    #####  2. Cell-wise parameter estimation
    # Dropout rate
    # Fit dropout rate with GLM
    # Drop-out estimation is independent between cells and can be parallelised.
    matLinModelPi <- do.call(rbind, 
      bplapply(seq(1, scaNumCells), function(cell){
        scaindIntervalStart <- max(1,cell-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,cell+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        
        vecLinModelPi <- fitPiZINB_LinPulse(
          vecLinModelPi=matLinModelPi[cell,],
          matPredictorsPi=cbind(1,log(matMu[,cell]),matConstPredictorsPi),
          vecCounts=matCountsProc[,cell],
          matMu=matMu[,vecInterval],
          matDisp=matDispersions[,vecInterval],
          vecNormConst=vecSizeFactors[vecInterval])
        vecLinModelOut <- cbind(1,log(matMu[,cell]),matConstPredictorsPi) %*% vecLinModelPi
        vecDropout <- 1/(1+exp(-vecLinModelOut))
        return(vecLinModelPi)
      }))
    matDropout <- do.call(cbind, 
      lapply(seq(1, scaNumCells), function(cell){
        vecLinModelOut <- cbind(1,log(matMu[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
        vecDropout <- 1/(1+exp(-vecLinModelOut))
        return(vecDropout)
      }))
    if(boolSuperVerbose){
      print(paste0("Drop-out rate estimation did not converge in ", 
        sum(vecboolConvergedGLMdisp), " cases."))
      scaLogLikTemp <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecSizeFactors=vecSizeFactors,
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,],
            scaWindowRadius=scaWindowRadius)
        })
      ))
      print(paste0("2) Drop-out estimation complete: loglikelihood of ", scaLogLikTemp))
    }  
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    if(is.null(scaWindowRadius)){
      scaLogLikNew <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,]*vecSizeFactors,
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,])
        })
      ))
    } else {
      scaLogLikNew <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecSizeFactors=vecSizeFactors,
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,],
            scaWindowRadius=scaWindowRadius)
        })
      ))
    }
    
    # EM-iteration complete
    if(verbose){print(paste0("Completed iteration ", scaIter, " with log likelihood of ", scaLogLikNew))}
    vecEMLogLik[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  # Evaluate convergence
  #if(all(as.logical(vecboolConvergedGLMdisp)) & all(as.logical(vecboolConvergedGLMpi)) &
  if(all(as.logical(vecboolConvergedGLMdisp)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergence <- TRUE
  } else {
    boolConvergence <- FALSE
  }
  
  # Compute mixture probabilities and imputed counts  
  matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
  }))
  matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    vecZ <- sapply(seq(1,scaNumCells), function(j){
      if(matboolNotZeroObserved[i,j]){
        scaZ <- 0
      } else {
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        scaZ <- sum(matDropout[i,j]/(matDropout[i,j] + 
            (1-matDropout[i,j])*matNBZero[i,vecInterval])) *
          1/length(vecInterval)
      }
      return(scaZ)
    })
    return(vecZ)
  }))
  matProbNB <- 1 - matZ
  
  # Name rows and columns of output
  rownames(matDropout) <- rownames(matCountsProc)
  colnames(matDropout) <- colnames(matCountsProc)
  rownames(matMu) <- rownames(matCountsProc)
  colnames(matMu) <- colnames(matCountsProc)
  names(vecDispersions) <- rownames(matCountsProc)
  rownames(matDispersions) <- rownames(matCountsProc)
  colnames(matDispersions) <- colnames(matCountsProc)
  rownames(matProbNB) <- rownames(matCountsProc)
  colnames(matProbNB) <- colnames(matCountsProc)
  if(is.null(scaWindowRadius)){
    rownames(matMuCluster) <- rownames(matCountsProc)
  } else {
    matMuCluster <- NA
  }
  
  # Check dispersions
  if(any(is.na(matDispersions) | !is.finite(matDispersions))){
    matDispersions[is.na(matDispersions) | !is.finite(matDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions. Set to 1.")
  }
  
  return(list( vecDispersions=vecDispersions,
    matDropout=matDropout,
    matDropoutLinModel=matLinModelPi,
    matProbNB=matProbNB,
    matMuCluster=matMuCluster,
    matMu=matMu,
    boolConvergence=boolConvergence,
    vecEMLogLik=vecEMLogLik ))
}