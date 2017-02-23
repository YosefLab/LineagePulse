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
#' of negative binomial dispersion and mean paramater on single gene given
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
#' @aliases evalLogLikDispConstMuConstZINB_comp
#'
#' @seealso Called by fitting wrapper:
#' \code{fitDispConstMuConstZINB}.
#' Calls \code{evalLogLikGene}.
#' 
#' @param vecTheta: (numeric vector length 2) 
#'    Log of dispersion parameter 
#'    and log of mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuConstZINB <- function(vecTheta,
                                           vecCounts,
                                           lsMuModelGlobal,
                                           matDropoutLinModel,
                                           vecPiConstPredictors,
                                           vecboolNotZero,
                                           vecboolZero,
                                           scaWindowRadius=NULL ){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2])
  scaNParamUsed <- 1+length(vecMu)
  
  # (II) Prevent parameter shrinkage/explosion
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10) 
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  vecMu[vecMu < 10^(-10)] <- 10^(-10)
  vecMu[vecMu > 10^(10)] <- 10^(10)
  
  vecMuParam <- rep(vecMu, length(vecCounts))
  
  # Extract batch factors
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchParam <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxBatchAssign]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
      vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
      
      vecMuParam <- vecMuParam*vecBatchParam
    }
    #vecMuParam[vecMuParam < 10^(-10)] <- 10^(-10)
  }
  
  # (III) Compute drop-out rates
  vecPi <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
                                       vecMu=vecMuParam,
                                       vecPiConstPredictors=vecPiConstPredictors )
  
  # (IV) Evaluate loglikelihood of estimate
  scaLogLik <- evalLogLikGene(vecCounts=vecCounts,
                              vecMu=vecMuParam,
                              vecNormConst=lsMuModelGlobal$vecNormConst,
                              vecDisp=vecDisp, 
                              vecPi=vecPi,
                              vecboolNotZero=vecboolNotZero, 
                              vecboolZero=vecboolZero,
                              scaWindowRadius=scaWindowRadius )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compiled function: evalLogLikDispConstMuConstZINB
#' 
#' Pre-compile heavily used function.
#' Refer to \link{evalLogLikDispConstMuConstZINB}.
#' 
#' @seealso \link{evalLogLikDispConstMuConstZINB}
#' 
#' @param vecTheta: (numeric vector length 2) 
#'    Log of dispersion parameter 
#'    and log of mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuConstZINB_comp <- cmpfun(evalLogLikDispConstMuConstZINB)

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
#' @seealso Called by fitting wrapper:
#' \code{fitDispConstMuVecWindowsZINB}.
#' Calls \code{evalLogLikGene}.
#' Compiled function: \link{evalLogLikDispConstMuVecWindowsZINB_comp}.
#' 
#' @param vecTheta: (numeric vector 1 + number of cells) 
#'    Log of dispersion and log of mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuVecWindowsZINB <- function(vecTheta,
                                                vecCounts,
                                                lsMuModelGlobal,
                                                matDropoutLinModel,
                                                vecPiConstPredictors,
                                                vecboolNotZero,
                                                vecboolZero,
                                                scaWindowRadius){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2:(1+length(vecCounts))])
  scaNParamUsed <- 1+length(vecMu)
  
  # (II) Prevent parameter shrinkage/explosion
  if(scaDisp < 10^(-10)){ scaDisp <- 10^(-10) }
  if(scaDisp > 1/10^(-10)){ scaDisp <- 1/10^(-10) }
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  vecMu[vecMu < 10^(-10)] <- 10^(-10)
  vecMu[vecMu > 10^(10)] <- 10^(10)
  
  vecMuParam <- vecMu
  
  # Extract batch factors
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchParam <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxBatchAssign]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
      vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
      
      vecMuParam <- vecMuParam*vecBatchParam
    }
    #vecMuParam[vecMuParam < 10^(-10)] <- 10^(-10)
  }
  
  # (III) Compute drop-out rates
  vecPi <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
                                       vecMu=vecMuParam,
                                       vecPiConstPredictors=vecPiConstPredictors )
  
  # (IV) Evaluate loglikelihood (this is the cost function) 
  scaLogLik <- evalLogLikSmoothZINB_comp(
    vecCounts=vecCounts,
    vecMu=vecMuParam,
    vecNormConst=lsMuModelGlobal$vecNormConst,
    vecDisp=vecDisp, 
    vecPi=vecPi,
    vecboolNotZero=vecboolNotZero, 
    vecboolZero=vecboolZero,
    scaWindowRadius=scaWindowRadius)
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compiled function: evalLogLikDispConstMuVecWindowsZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to: \link{evalLogLikDispConstMuVecWindowsZINB}.
#' 
#' @param vecTheta: (numeric vector 1 + number of cells) 
#'    Log of dispersion and log of mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuVecWindowsZINB_comp <- cmpfun(evalLogLikDispConstMuVecWindowsZINB)


#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and cluster mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersion and mean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled by cluster for the given gene.
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
#' @seealso Called by fitting wrapper:
#' \code{fitDispConstMuClusterZINB}.
#' Calls \code{evalLogLikGene}.
#' Compiled function: \link{evalLogLikDispConstMuClustersZINB_comp}.
#' 
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuClustersZINB <- function(vecTheta,
                                              vecCounts,
                                              lsMuModelGlobal,
                                              matDropoutLinModel,
                                              vecPiConstPredictors,
                                              vecboolNotZero,
                                              vecboolZero ){ 
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2:(max(lsMuModelGlobal$vecidxClusterAssign)+1)])
  scaNParamUsed <- 1+length(vecMu)
  
  # (II) Prevent parameter shrinkage/explosion
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  vecMu[vecMu < 10^(-10)] <- 10^(-10)
  vecMu[vecMu > 10^(10)] <- 10^(10)
  
  vecMuParam <- vecMu[lsMuModelGlobal$vecidxClusterAssign]
  
  # Extract batch factors
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchParam <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxBatchAssign]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
      vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
      
      vecMuParam <- vecMuParam*vecBatchParam
    }
    #vecMuParam[vecMuParam < 10^(-10)] <- 10^(-10)
  }
  
  # (III) Compute drop-out rates
  vecPi <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
                                       vecMu=vecMuParam,
                                       vecPiConstPredictors=vecPiConstPredictors )
  
  # (IV) Evaluate loglikelihood of estimate
  scaLogLik <- evalLogLikGene( vecCounts=vecCounts,
                               vecMu=vecMuParam,
                               vecNormConst=lsMuModelGlobal$vecNormConst,
                               vecDisp=vecDisp, 
                               vecPi=vecPi,
                               vecboolNotZero=vecboolNotZero, 
                               vecboolZero=vecboolZero )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compiled function: evalLogLikDispConstMuClustersZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikDispConstMuClustersZINB}.
#' 
#' @seealso \link{evalLogLikDispConstMuClustersZINB}.
#'
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuClustersZINB_comp <- cmpfun(evalLogLikDispConstMuClustersZINB)

#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and mixture model mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersion and mean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled as a mixture model (one constant per
#' mixture) for the given gene.
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
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param matWeights: (probability matrix mixture components x cells)
#'    Mixture assignments of each cell to each mixture component.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial log-likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuMMZINB <- function(vecTheta,
                                        vecCounts,
                                        lsMuModelGlobal,
                                        matDropoutLinModel,
                                        vecPiConstPredictors,
                                        vecboolNotZero,
                                        vecboolZero,
                                        matWeights){ 
  
  scaNCells <- length(vecCounts)
  scaNMixtures <- dim(matWeights)[1]
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker function to fit positive mean
  vecMu <- exp(vecTheta[2:(scaNMixtures+1)])
  scaNParamUsed <- 1+length(vecMu)
  
  # (II) Prevent parameter shrinkage/explosion
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10) 
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10) 
  
  vecMu[vecMu < 10^(-10)] <- 10^(-10)
  vecMu[vecMu > 10^(10)] <- 10^(10)
  
  matMuParam <- matrix(vecMu, nrow=scaNCells, ncol=scaNMixtures, byrow=TRUE)
  
  # Extract batch factors
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchParam <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxBatchAssign]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
      vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
      
      matMuParam <- matMuParam*matrix(vecBatchParam, nrow=scaNCells, ncol=scaNMixtures, byrow=FALSE)
    }
  }
  
  # Decompress parameters
  matDropParam <- do.call(cbind, lapply(seq(1,scaNMixtures), function(m){
    decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
                                vecMu=matMuParam[,m],
                                vecPiConstPredictors=vecPiConstPredictors )
  }))
  scaLogLik <- evalLogLikGeneMM(vecCounts=vecCounts,
                                matMuParam=matMuParam,
                                vecNormConst=lsMuModelGlobal$vecNormConst,
                                matDispParam=matrix(scaDisp, nrow=scaNCells, ncol=scaNMixtures),
                                matDropParam=matDropParam,
                                matWeights=matWeights,
                                vecboolNotZero= !is.na(vecCounts) & vecCounts>=0, 
                                vecboolZero= !is.na(vecCounts) & vecCounts==0 )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compiled function: evalLogLikDispConstMuMMZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikDispConstMuMMZINB}.
#' 
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the log mean parameter. 
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param matWeights: (probability matrix mixture components x cells)
#'    Mixture assignments of each cell to each mixture component.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial log-likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuMMZINB_comp <- cmpfun(evalLogLikDispConstMuMMZINB)

#' Cost function zero-inflated negative binomial model for dispersion
#' and mean co-estimation under constant dispersion and constant mean model
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow simultaneous numerical optimisation
#' of negative binomial dispersionmean paramater on single gene given
#' the drop-out rate/model. The dispersion parameter is modelled as a constant
#' and the mean parameter is modelled by the impulse model for the given gene.
#' 
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
#' @seealso Called by fitting wrapper:
#' \code{fitDispConstMuImpulseOneInitZINB}.
#' Calls \code{evalImpulseModel} and \code{evalLogLikGene}.
#' Compiled function: \link{evalLogLikDispConstMuImpulseZINB_comp}.
#' 
#' @param vecTheta: (numeric vector dispersion (1) and impulse parameters (6)) 
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
#'    Unique (pseudo)time coordinates of cells.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to cell in list of sorted
#'    time points. vecTimepoints[vecindTimepointAssign]==vecPseudotime
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuImpulseZINB <- function(vecTheta,
                                             vecCounts,
                                             lsMuModelGlobal,
                                             vecTimepoints,
                                             vecindTimepointAssign,
                                             matDropoutLinModel,
                                             vecPiConstPredictors,
                                             vecboolNotZero, 
                                             vecboolZero,
                                             scaWindowRadius=NULL){  
  
  # (I) Linker functions
  # Log linker function to fit positive dispersion factor
  scaDisp <- exp(vecTheta[1])
  # Log linker for amplitudes
  vecImpulseParam <- vecTheta[2:8]
  vecImpulseParam[3:5] <- exp(vecImpulseParam[3:5])
  scaNParamUsed <- 1+length(vecImpulseParam)
  
  # (II) Prevent parameter shrinkage/explosion
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  vecDisp <- rep(scaDisp, length(vecCounts))
  
  # Prevent amplitude shrinkage/explosion
  vecImpulseParam[3:5][vecImpulseParam[3:5] < 10^(-10)] <- 10^-10
  vecImpulseParam[3:5][vecImpulseParam[3:5] > 10^(10)] <- 10^10
  
  vecMuParam <- evalImpulseModel_comp(vecImpulseParam=vecImpulseParam,
                                      vecTimepoints=vecTimepoints)[vecindTimepointAssign]
  
  # Extract batch factors
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchParam <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxBatchAssign]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
      vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
      
      vecMuParam <- vecMuParam*vecBatchParam
    }
    #vecMuParam[vecMuParam < 10^(-10)] <- 10^(-10)
  }
  
  # (III) Compute drop-out rates
  vecPi <- decompressDropoutRateByGene(matDropModel=matDropoutLinModel,
                                       vecMu=vecMuParam,
                                       vecPiConstPredictors=vecPiConstPredictors )
  
  # (IV) Evaluate loglikelihood of estimate
  scaLogLik <- evalLogLikGene(vecCounts=vecCounts,
                              vecMu=vecMuParam,
                              vecNormConst=lsMuModelGlobal$vecNormConst,
                              vecDisp=vecDisp, 
                              vecPi=vecPi,
                              vecboolNotZero=vecboolNotZero, 
                              vecboolZero=vecboolZero,
                              scaWindowRadius=scaWindowRadius )
  
  return(scaLogLik)
}

#' Compiled function: evalLogLikDispConstMuImpulseZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikDispConstMuImpulseZINB}.
#' 
#' @seealso \link{evalLogLikDispConstMuImpulseZINB}
#' 
#' @param vecTheta: (numeric vector dispersion (1) and impulse parameters (6)) 
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
#'    Unique (pseudo)time coordinates of cells.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to cell in list of sorted
#'    time points. vecTimepoints[vecindTimepointAssign]==vecPseudotime
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikDispConstMuImpulseZINB_comp <- cmpfun(evalLogLikDispConstMuImpulseZINB)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) FITTING COORDINATORS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Numerical fitting wrapper for constant dispersion and constant mean model
#' 
#' Fits single negative binomial mean and dispersion parameter numerically as 
#' maximum likelihood estimators to a gene: Constant dispersion and
#' constant mean model.
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper \code{fitZINBMuDisp}.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikDispConstMuConstZINB}.
#' 
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param scaDispGuess: (scalar) Initialisation of dispersion parameters.
#' @param scaMuGuess: (scalar) Initialisation of mean parameters.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item scaMu: (scalar) Negative binomial mean 
#'        parameter estimate.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @export

fitDispConstMuConstZINB <- function(vecCounts,
                                    scaMuGuess,
                                    lsvecBatchParamGuess,
                                    lsMuModelGlobal,
                                    scaDispGuess,
                                    matDropoutLinModel,
                                    vecPiConstPredictors,
                                    scaWindowRadius,
                                    MAXIT=1000,
                                    RELTOL=10^(-8) ){ 
  
  # (I) Numerical maximum likelihood estimation
  if(!is.null(lsvecBatchParamGuess)) vecTheta <- c(log(scaDispGuess), log(scaMuGuess), log(unlist(lsvecBatchParamGuess)))
  else vecTheta <- c(log(scaDispGuess), log(scaMuGuess))
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=vecTheta,
      evalLogLikDispConstMuConstZINB_comp,
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=MAXIT, 
                   reltol=RELTOL,
                   fnscale=-1) )[c("par","value","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuConstZINB().",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    scaLLInit <- evalLogLikDispConstMuConstZINB_comp(
      vecTheta=c(log(scaDispGuess), log(scaMuGuess)),
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius)
    print(paste0("c(log(scaDispGuess), log(scaMuGuess)) ", 
                 paste(c(log(scaDispGuess), log(scaMuGuess)),collapse=" ")))
    print(paste0("scaLLInit ", scaLLInit))
    lsErrorCausingGene <- list(
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      vecParamGuess=c(log(scaDispGuess), log(scaMuGuess)),
      matDropoutLinModel=matDropoutLinModel, 
      vecPiConstPredictors=vecPiConstPredictors,
      scaWindowRadius=scaWindowRadius,
      scaLLInit=scaLLInit)
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  # Catch boundary of likelihood domain on dispersion space
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(scaDisp > 10^(10)) scaDisp <- 10^(10)
  
  vecMuModel <- exp(fitDispMu[2])
  # Catch boundary of likelihood domain on mu space
  vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
  vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
  
  scaNParamUsed <- 2
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    lsvecBatchFactors <- lapply(lsMuModelGlobal$lsvecidxBatchAssign, function(vecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactors <- c(1, exp(fitDispMu[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
      vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
      return(vecBatchFactors)
    })
  } else { lsvecBatchFactors <- NA }
  
  scaLL <- fitDispMu["value"]
  scaConvergence <- fitDispMu["convergence"]
  
  #Test
  vecMuParam <- decompressMeansByGene( vecMuModel=vecMuModel,
                                       lsvecBatchModel=lsvecBatchFactors,
                                       lsMuModelGlobal=lsMuModelGlobal,
                                       vecInterval=NULL )
  vecDispParam <- rep(scaDisp, length(vecCounts))
  vecPiParam <- decompressDropoutRateByGene( matDropModel=matDropoutLinModel,
                                             vecMu=vecMuParam,
                                             vecPiConstPredictors=vecPiConstPredictors )
  
  scaLogLik <- evalLogLikGene(vecCounts=vecCounts,
                              vecMu=vecMuParam,
                              vecNormConst=lsMuModelGlobal$vecNormConst,
                              vecDisp=vecDispParam, 
                              vecPi=vecPiParam,
                              vecboolNotZero=vecCounts>0, 
                              vecboolZero=vecCounts==0,
                              scaWindowRadius=scaWindowRadius )
  
  return(list(scaDisp=scaDisp,
              vecMuModel=vecMuModel,
              lsvecBatchFactors=lsvecBatchFactors,
              scaConvergence=scaConvergence,
              scaLL=scaLL))
}

#' Numerical fitting wrapper for constant dispersion and sliding window 
#' mean model
#' 
#' Fits single negative binomial mean and dispersion parameter numerically as 
#' maximum likelihood estimators to a gene: local smoothing (sliding windows) 
#' mean model and constant dispersion model.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper \code{fitZINBMuDisp}.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikDispConstMuVecWindowsZINB}.
#' 
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param vecMuGuess: (vector number of cells)
#'    Initialisation for mean parameters to be estimated.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item vecMu: (numeric vector number of cells)
#'        Negative binomial mean parameter estimates.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
fitDispConstMuVecWindowsZINB<- function(vecCounts,
                                        vecMuGuess,
                                        lsvecBatchParamGuess,
                                        lsMuModelGlobal,
                                        scaDispGuess,
                                        matDropoutLinModel,
                                        vecPiConstPredictors,
                                        scaWindowRadius=NULL,
                                        MAXIT=1000,
                                        RELTOL=10^(-8) ){
  
  if(!is.null(lsvecBatchParamGuess)) vecTheta <- c(log(scaDispGuess), log(vecMuGuess), log(unlist(lsvecBatchParamGuess)))
  else vecTheta <- c(log(scaDispGuess), log(vecMuGuess))
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=vecTheta,
      evalLogLikDispConstMuVecWindowsZINB_comp,
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) &vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS",
      control=list(maxit=MAXIT, 
                   reltol=RELTOL,
                   fnscale=-1) )[c("par","value","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuVecWindowsZINB().",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    scaLLInit <- evalLogLikDispConstMuVecWindowsZINB_comp(
      vecTheta=c(log(scaDispGuess), log(vecMuGuess)),
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0)
    print(paste0("c(log(scaDispGuess), log(vecMuGuess)) ", paste(c(scaDispGuess, vecMuGuess),collapse=" ")))
    print(paste0("scaLLInit ", scaLLInit))
    lsErrorCausingGene <- list(
      vecParamGuess=c(scaDispGuess, vecMuGuess), 
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel, 
      vecPiConstPredictors=vecPiConstPredictors,
      scaLLInit=scaLLInit )
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  
  vecMuModel <- exp(fitDispMu[2:(length(vecMuGuess)+1)])
  # Catch boundary of likelihood domain on mu space
  vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
  vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
  
  scaNParamUsed <- length(vecMuGuess)+1
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    lsvecBatchFactors <- lapply(lsMuModelGlobal$lsvecidxBatchAssign, function(vecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactors <- c(1, exp(fitDispMu[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
      vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
      return(vecBatchFactors)
    })
  } else { lsvecBatchFactors <- NA }
  
  scaLL <- fitDispMu["value"]
  scaConvergence <- fitDispMu["convergence"]
  
  return(list(scaDisp=scaDisp,
              vecMuModel=vecMuModel,
              lsvecBatchFactors=lsvecBatchFactors,
              scaConvergence=scaConvergence,
              scaLL=scaLL))
}

#' Numerical fitting wrapper for constant dispersion
#' cluster-wise mean model
#' 
#' Fits single negative binomial mean and dispersion parameter numerically as 
#' maximum likelihood estimators to a gene: Constant dispersion and cluster-wise
#' mean model.
#' Note that this cannot be performed with
#' fitDispConstMuConstZINB for each cluster as the dispersion
#' parameter is assumed to be shared between clusters.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper \code{fitZINBMuDisp}.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikDispConstMuClustersZINB}.
#' 
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param vecMuGuess: (vector number of clusters)
#'    Initialisation for mean parameters to be estimated.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecClusterAssign: (integer vector length number of
#'    cells) Index of cluster assigned to each cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item vecMu: (numeric vector number of clusters)
#'        Negative binomial mean parameter estimates for each cluster.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
fitDispConstMuClusterZINB <- function(vecCounts,
                                      vecMuGuess,
                                      lsvecBatchParamGuess,
                                      lsMuModelGlobal,
                                      scaDispGuess,
                                      matDropoutLinModel,
                                      vecPiConstPredictors,
                                      MAXIT=1000,
                                      RELTOL=10^(-8) ){ 
  
  # scaWindowRadius is set to NULL because smoothing
  # within clusters does't make sense - the clusters already impose
  # a constraint on the means.
  if(!is.null(lsvecBatchParamGuess)) vecTheta <- c(log(scaDispGuess), log(vecMuGuess), log(unlist(lsvecBatchParamGuess)))
  else vecTheta <- c(log(scaDispGuess), log(vecMuGuess))
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=vecTheta,
      evalLogLikDispConstMuClustersZINB_comp,
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) &vecCounts==0,
      method="BFGS",
      control=list(maxit=MAXIT, 
                   reltol=RELTOL,
                   fnscale=-1) )[c("par","value","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuClusterZINB().",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    scaLLInit <- evalLogLikDispConstMuClustersZINB_comp(
      vecTheta=c(log(scaDispGuess), log(vecMuGuess)),
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecClusterAssign=vecClusterAssign,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0)
    print(paste0("c(log(scaDispGuess), log(vecMuGuess)) ", paste(c(log(scaDispGuess), log(vecMuGuess)),collapse=" ")))
    print(paste0("scaLLInit ", scaLLInit))
    lsErrorCausingGene <- list(
      vecParamGuess=c(scaDispGuess, vecMuGuess), 
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel, 
      vecPiConstPredictors=vecPiConstPredictors,
      vecClusterAssign=vecClusterAssign, 
      scaLLInit=scaLLInit )
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  
  vecMuModel <- exp(fitDispMu[2:(length(vecMuGuess)+1)])
  # Catch boundary of likelihood domain on mu space
  vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
  vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
  
  scaNParamUsed <- length(vecMuGuess)+1
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    lsvecBatchFactors <- lapply(lsMuModelGlobal$lsvecidxBatchAssign, function(vecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactors <- c(1, exp(fitDispMu[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
      vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
      return(vecBatchFactors)
    })
  } else { lsvecBatchFactors <- NA }
  
  scaLL <- fitDispMu["value"]
  scaConvergence <- fitDispMu["convergence"]
  
  return(list(scaDisp=scaDisp,
              vecMuModel=vecMuModel,
              lsvecBatchFactors=lsvecBatchFactors,
              scaConvergence=scaConvergence,
              scaLL=scaLL))
}

#' Numerical fitting wrapper for constant dispersion
#' mixture-model mean model
#' 
#' Fits single negative binomial mean and dispersion parameter numerically as 
#' maximum likelihood estimators to a gene: Constant dispersion and mixture-model
#' mean model.
#' Note that this cannot be performed with
#' fitDispConstMuConstZINB for each cluster as the dispersion
#' parameter is assumed to be shared between clusters.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper \code{fitZINBMuDisp}.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikDispConstMuClustersZINB}.
#' 
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param vecMuGuess: (vector number of clusters)
#'    Initialisation for mean parameters to be estimated.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matWeights: (probability matrix mixture components x cells)
#'    Mixture assignments of each cell to each mixture component.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item vecMu: (numeric vector number of clusters)
#'        Negative binomial mean parameter estimates for each cluster.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
fitDispConstMuMMZINB <- function(vecCounts,
                                 vecMuGuess,
                                 lsvecBatchParamGuess,
                                 lsMuModelGlobal,
                                 scaDispGuess,
                                 matDropoutLinModel,
                                 vecPiConstPredictors,
                                 matWeights,
                                 MAXIT=1000,
                                 RELTOL=10^(-8) ){ 
  
  # scaWindowRadius is set to NULL because smoothing
  # within clusters does't make sense - the clusters already impose
  # a constraint on the means.
  if(!is.null(lsvecBatchParamGuess)) vecTheta <- c(log(scaDispGuess), log(vecMuGuess), log(unlist(lsvecBatchParamGuess)))
  else vecTheta <- c(log(scaDispGuess), log(vecMuGuess))
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=vecTheta,
      evalLogLikDispConstMuMMZINB_comp,
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      matWeights=matWeights,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) &vecCounts==0,
      method="BFGS",
      control=list(maxit=MAXIT, 
                   reltol=RELTOL,
                   fnscale=-1) )[c("par","value","convergence")])
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitDispConstMuMMZINB().",
                 " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    scaLLInit <- evalLogLikDispConstMuMMZINB_comp(
      vecTheta=c(log(scaDispGuess), log(vecMuGuess)),
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      matWeights=matWeights,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0)
    print(paste0("c(log(scaDispGuess), log(vecMuGuess)) ", paste(c(log(scaDispGuess), log(vecMuGuess)),collapse=" ")))
    print(paste0("scaLLInit ", scaLLInit))
    lsErrorCausingGene <- list(
      vecParamGuess=c(scaDispGuess, vecMuGuess), 
      vecCounts=vecCounts, 
      lsMuModelGlobal=lsMuModelGlobal,
      matDropoutLinModel=matDropoutLinModel, 
      vecPiConstPredictors=vecPiConstPredictors,
      matWeights=matWeights, 
      scaLLInit=scaLLInit )
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  if(scaDisp > 10^(10)) scaDisp <- 10^(10)
  
  vecMuModel <- exp(fitDispMu[2:(length(vecMuGuess)+1)])
  # Catch boundary of likelihood domain on mu space
  vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
  vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
  
  scaNParamUsed <- length(vecMuGuess)+1
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    lsvecBatchFactors <- lapply(lsMuModelGlobal$lsvecidxBatchAssign, function(vecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactors <- c(1, exp(fitDispMu[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
      vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
      return(vecBatchFactors)
    })
  } else { lsvecBatchFactors <- NA }
  
  scaLL <- fitDispMu["value"]
  scaConvergence <- fitDispMu["convergence"]
  
  return(list(scaDisp=scaDisp,
              vecMuModel=vecMuModel,
              lsvecBatchFactors=lsvecBatchFactors,
              scaConvergence=scaConvergence,
              scaLL=scaLL))
}

#' Numerical fitting wrapper for constant dispersion
#' impulse mean model for single initialisation
#' 
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the impulse model and a constant
#' dispersion parameter and returns the fitted (maximum likelihood) model.
#' This is the wrapper that calls optim.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called by \code{fitDispConstMuImpulseZINB}. This function
#' performs optimisation of one impulse model initialisation,
#' \code{fitDispConstMuImpulseZINB} coordinates the overall fitting
#' of an impulse model to a gene.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikDispConstMuImpulseZINB}.
#' 
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param vecParamGuessPeak (numeric vector number of parameters [6]) 
#'    Initialisation for impulse model for mean parameters.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecTimepoints: (numerical vector number of unique time coordinates)
#'    Unique (pseudo)time coordinates of cells.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to cell in list of sorted
#'    time points. vecTimepoints[vecindTimepointAssign]==vecPseudotime
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' @param MAXIT: (integer) [Default 1000] Maximum number of 
#'    estimation iterations optim.
#' @param RELTOL: (scalar) [Default sqrt(10^(-10))]
#'    Relative tolerance for optim.
#'  @param trace: (integer) optim control parameter for reporting.
#'  @param REPORT: (integer) optim control parameter for reporting.
#'   
#' @return list: (length 4)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item vecImpulseParam: (numeric vector length 6)
#'        {beta, log(h0), log(h1), log(h2), t1, t2}
#'        Impulse model parameter estimates.
#'      \item scaLL: (scalar) Loglikelihood of fit.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
fitDispConstMuImpulseOneInitZINB <- function(vecCounts,
                                             vecImpulseParamGuess,
                                             lsvecBatchParamGuess,
                                             lsMuModelGlobal,
                                             scaDispGuess,
                                             vecTimepoints,
                                             vecindTimepointAssign, 
                                             matDropoutLinModel,
                                             vecPiConstPredictors,
                                             scaWindowRadius=NULL,
                                             MAXIT=1000,
                                             RELTOL=10^(-8) ){
  
  if(!is.null(lsvecBatchParamGuess)) vecTheta <- c(log(scaDispGuess), vecImpulseParamGuess, log(unlist(lsvecBatchParamGuess)))
  else vecTheta <- c(log(scaDispGuess), vecImpulseParamGuess)
  
  fitDispMu <- tryCatch({
    unlist(optim(
      par=vecTheta,			
      fn=evalLogLikDispConstMuImpulseZINB_comp, 
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      vecTimepoints=vecTimepoints,
      vecindTimepointAssign=vecindTimepointAssign,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius,
      method="BFGS", 
      control=list(maxit=MAXIT, 
                   reltol=RELTOL,
                   fnscale=-1)
    )[c("par","value","convergence")] )
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting impulse model: fitDispConstMuImpulseZINB().",
                 " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(strErrorMsg)
    scaLLInit <- evalLogLikDispConstMuImpulseZINB_comp(
      vecTheta=c(log(scaDispGuess), vecImpulseParamGuess),
      vecCounts=vecCounts,
      lsMuModelGlobal=lsMuModelGlobal,
      vecTimepoints=vecTimepoints,
      vecindTimepointAssign=vecindTimepointAssign,
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      scaWindowRadius=scaWindowRadius)
    print(paste0("c(log(scaDispGuess), vecImpulseParamGuess) ", paste(c(log(scaDispGuess), vecImpulseParamGuess),collapse=" ")))
    print(paste0("scaLLInit ", scaLLInit))
    lsErrorCausingGene <- list(
      vecParamGuess=c(scaDispGuess, vecImpulseParamGuess), 
      vecCounts=vecCounts, 
      vecTimepoints=vecTimepoints, 
      vecindTimepointAssign=vecindTimepointAssign, 
      matDropoutLinModel=matDropoutLinModel, 
      vecPiConstPredictors=vecPiConstPredictors,
      scaWindowRadius=scaWindowRadius, 
      scaLLInit=scaLLInit )
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  scaDisp <- exp(fitDispMu[1])
  if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
  if(scaDisp > 1/10^(-10)) scaDisp <- 1/10^(-10)
  
  vecMuModel <- fitDispMu[2:8]
  vecMuModel[3:5] <- exp(vecMuModel[3:5])
  vecMuModel[3:5][vecMuModel[3:5] < 10^(-10)] <- 10^(-10)
  vecMuModel[3:5][vecMuModel[3:5] > 10^(10)] <- 10^(10)
  
  scaNParamUsed <- 1+length(vecMuModel)
  if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
    lsvecBatchFactors <- lapply(lsMuModelGlobal$lsvecidxBatchAssign, function(vecidxBatchAssign){
      scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactors <- c(1, exp(fitDispMu[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
      vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
      return(vecBatchFactors)
    })
  } else { lsvecBatchFactors <- NA }
  
  scaLL <- fitDispMu["value"]
  scaConvergence <- fitDispMu["convergence"]
  
  return( list(scaDisp=scaDisp,
               vecMuModel=vecMuModel,
               lsvecBatchFactors=lsvecBatchFactors,
               scaLL=scaLL,
               scaConvergence=scaConvergence) )
}

#' Numerical fitting wrapper for constant dispersion
#' impulse mean model for multiple initialisation
#' 
#' Computes impulse parameter initialisation for valley
#' and peak model and uses both and the prior parameter fit
#' in three separate optimisation runs to obtain the best 
#' impulse model fit to the data, simultaneous with fitting a 
#' constant dispersion factor.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper \code{fitZINBMuDisp}.
#' Calls optimisation wrapper \code{fitDispConstMuImpulseOneInitZINB} 
#' for each initialisation.
#' 
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param scaDispGuess: (scalar) Initialisation for dispersion parameter
#'    to be estimated.
#' @param vecImpulseParamGuess (numeric vector number of parameters [6]) 
#'    Initialisation for impulse model for mean parameters.
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
#'      \item vecClusterAssign: (integer vector length number of
#'    cells) [Default NA] Index of cluster assigned to each cell.
#'    Used for clusters model.
#'      \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'    of a gene are simultaneously estiamted as a vector with BFGS
#'    in windows mode.
#'      \item MAXIT_BFGS_MuDisp: (int) Maximum number of iterations
#'    for BFGS estimation of impulse model with optim (termination criterium).
#'      \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'    change in objective function for BFGS estimation of impulse 
#'    model with optim (termination criterium).
#'    }
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param matDropoutLinModel: (matrix number of cells x number of predictors)
#'    Logistic linear model parameters of the dropout rate 
#'    as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors: (numeric vector constant gene-wise coefficients)
#'    Constant gene-wise coeffiecients, i.e. predictors which are not
#'    the offset and not the mean parameter.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#' 
#' @return list (length 3)
#'    \itemize{
#'      \item scaDisp: (scalar) Negative binomial dispersion 
#'        parameter estimate. 
#'      \item vecImpulseParam: (numeric vector length 6)
#'        {beta, log(h0), log(h1), log(h2), t1, t2}
#'        Impulse model parameter estimates.
#'      \item scaConvergence: (scalar) Convergence status of optim.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @export
fitDispConstMuImpulseZINB <- function(vecCounts, 
                                      vecImpulseParamGuess,
                                      lsvecBatchParamGuess,
                                      lsMuModelGlobal,
                                      scaDispGuess,
                                      matDropoutLinModel,
                                      vecPiConstPredictors,
                                      scaWindowRadius,
                                      MAXIT=1000,
                                      RELTOL=10^(-8) ){
  
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
                                       lsvecBatchModel=lsvecBatchParamGuess,
                                       lsMuModelGlobal=lsMuModelGlobal,
                                       vecInterval=NULL )
  vecPiParam <- decompressDropoutRateByGene( matDropModel=matDropoutLinModel,
                                             vecMu=vecMuParam,
                                             vecPiConstPredictors=vecPiConstPredictors )
  # The previous parameter estiamte is kept as a reference and
  # used as an initialisation
  # Compute initialisations for peak and valley
  lsParamGuesses <- initialiseImpulseParametes(vecCounts=vecCounts,
                                               lsMuModelGlobal=lsMuModelGlobal,
                                               vecMu=vecMuParam,
                                               vecDisp=rep(scaDispGuess, length(vecCounts)),
                                               vecDrop=vecPiParam,
                                               scaWindowRadius=scaWindowRadius)
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # (II) Compute new parameters
  # 1. Initialisation: Prior best fit
  vecParamGuessPrior <- vecImpulseParamGuess
  vecParamGuessPrior[3:5] <- log(vecParamGuessPrior[3:5])
  lsFitPrior <- fitDispConstMuImpulseOneInitZINB(
    vecCounts=vecCounts,
    vecImpulseParamGuess=vecParamGuessPrior,
    lsvecBatchParamGuess=lsvecBatchParamGuess,
    lsMuModelGlobal=lsMuModelGlobal,
    scaDispGuess=scaDispGuess,
    vecTimepoints=vecTimepoints, 
    matDropoutLinModel=matDropoutLinModel,
    vecPiConstPredictors=vecPiConstPredictors,
    vecindTimepointAssign=vecindTimepointAssign,
    scaWindowRadius=scaWindowRadius,
    MAXIT=MAXIT, 
    RELTOL=RELTOL)
  if(boolUseNewInits){
    # 2. Initialisation: Peak
    lsFitPeak <- fitDispConstMuImpulseOneInitZINB(
      vecCounts=vecCounts,
      vecImpulseParamGuess=vecParamGuessPeak,
      lsvecBatchParamGuess=lsvecBatchParamGuess,
      lsMuModelGlobal=lsMuModelGlobal,
      scaDispGuess=scaDispGuess,
      vecTimepoints=vecTimepoints, 
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=MAXIT, 
      RELTOL=RELTOL)
    # 3. Initialisation: Valley
    lsFitValley <- fitDispConstMuImpulseOneInitZINB(
      vecCounts=vecCounts,
      vecImpulseParamGuess=vecParamGuessValley,
      lsvecBatchParamGuess=lsvecBatchParamGuess,
      lsMuModelGlobal=lsMuModelGlobal,
      scaDispGuess=scaDispGuess,
      vecTimepoints=vecTimepoints, 
      matDropoutLinModel=matDropoutLinModel,
      vecPiConstPredictors=vecPiConstPredictors,
      vecindTimepointAssign=vecindTimepointAssign,
      scaWindowRadius=scaWindowRadius,
      MAXIT=MAXIT, 
      RELTOL=RELTOL)
    
    # (IV) Find best fit
    lsFits <- list(lsFitPeak, lsFitValley, lsFitPrior)
    vecLL <- sapply(lsFits , function(fit) fit$scaLL)
    if(all(!is.na(vecLL))){
      # If any impulse fitting (of the three initialisations)
      # was successful.
      # Chose best value
      indMaxLL <- match(max(vecLL, na.rm=TRUE), vecLL)
      lsFitBest <- lsFits[[indMaxLL]]
    } else if(is.na(vecLL[3])){
      # If optimisation of previous fit was not successfull:
      # Make sure new value is better than previous
      scaLLGuess <- evalLogLikGene(vecCounts=vecCounts,
                                   vecMu=vecMuParam*vecNormConst,
                                   vecDisp=rep(scaDispGuess, length(vecCounts)), 
                                   vecPi=vecPiParam,
                                   vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
                                   vecboolZero= !is.na(vecCounts) &vecCounts==0,
                                   scaWindowRadius=scaWindowRadius)
      indMaxLL <- match(max(vecLL, na.rm=TRUE), vecLL)
      if(vecLL[indMaxLL] < scaLLGuess){
        lsFitBest <- list(scaDisp=scaDispGuess,
                          vecImpulseParam=vecImpulseParamGuess,
                          scaConvergence=1002)
      } else {
        lsFitBest <- lsFits[[indMaxLL]]
      }
    } else {
      # If none of the three initilisations was successful:
      # Use prior paramter values
      lsFitBest <- list(scaDisp=scaDispGuess,
                        vecImpulseParam=vecImpulseParamGuess,
                        scaConvergence=1001)
    }
  } else {
    # Make sure new value is better than previous
    scaLLGuess <- evalLogLikGene(vecCounts=vecCounts,
                                 vecMu=vecMuParam*lsMuModelGlobal$vecNormConst,
                                 vecDisp=rep(scaDispGuess, length(vecCounts)), 
                                 vecPi=vecPiParam,
                                 vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
                                 vecboolZero= !is.na(vecCounts) &vecCounts==0,
                                 scaWindowRadius=scaWindowRadius)
    if(lsFitPrior$scaLL < scaLLGuess){
      lsFitBest <- list(scaDisp=scaDispGuess,
                        vecImpulseParam=vecImpulseParamGuess,
                        scaConvergence=1002)
    } else {
      lsFitBest <- lsFitPrior
    }
  }
  
  scaDisp <- lsFitBest$scaDisp
  vecMuModel <- lsFitBest$vecMuModel
  lsvecBatchFactors <- lsFitBest$lsvecBatchFactors
  scaConvergence <- lsFitBest$scaConvergence
  scaLL <- lsFitBest$scaLL
  
  return(list(scaDisp=scaDisp,
              vecMuModel=vecMuModel,
              lsvecBatchFactors=lsvecBatchFactors,
              scaConvergence=scaConvergence,
              scaLL=scaLL))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Top level auxillary function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Coordinate mean and dispersion parameter co-estimation step
#' 
#' Auxillary function that calls the estimation functions for the
#' different mean and dispersion models according to their needs. Note that one
#' function has to be coded for each combination of mean and dispersion
#' models.
#' 
#' @seealso Called by \code{fitZINB}. Calls fitting wrappers:
#' \code{fitDispConstMuConstZINB},
#' \code{fitDispConstMuClustersZINB},
#' \code{fitDispConstMuVecWindowsZINB} and
#' \code{fitDispConstMuImpulseZINB}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Observed read counts, not observed are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
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
#'          \item vecClusterAssign: (integer vector length number of
#'        cells) [Default NA] Index of cluster assigned to each cell.
#'        Used for clusters model.
#'          \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'        of a gene are simultaneously estiamted as a vector with BFGS
#'        in windows mode.
#'          \item MAXIT_BFGS_MuDisp: (int) Maximum number of iterations
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
#'          \item vecClusterAssign: (integer vector length number of
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
#'        Contains the mean model parameters according to the used model.
#'      \item matDispModel: (numeric matrix genes x disp model parameters)
#'        Contains the dispersion model parameters according to the used model.
#'      \item vecConvergence: (numeric vector number of genes) 
#'        Convergence status of optim for each gene.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @export
fitZINBMuDisp <- function( matCountsProc,
                           vecNormConst,
                           lsMuModel,
                           lsDispModel,
                           lsDropModel,
                           scaWindowRadius,
                           matWeights){
  
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  
  lsFitDispMu <- bplapply(seq(1,scaNumGenes), function(i){
    if(!is.null(lsMuModel$lsMuModelGlobal$vecConfounders)) lsvecBatchParamGuess <- lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,2:dim(mat)[2]] )
    else lsvecBatchParamGuess <- NULL
    
    # Nest mean models within dispersion models
    if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){
      if(lsMuModel$lsMuModelGlobal$strMuModel=="constant"){
        # Decompress parameters
        # Don't need to decompress any
        
        fitDispMu <- fitDispConstMuConstZINB(
          vecCounts=matCountsProc[i,],
          scaMuGuess=lsMuModel$matMuModel[i,],
          lsvecBatchParamGuess=lsvecBatchParamGuess,
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          scaDispGuess=lsDispModel$matDispModel[i,],
          matDropoutLinModel=lsDropModel$matDropoutLinModel,
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius,
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
        return(fitDispMu)
      } else if(lsMuModel$lsMuModelGlobal$strMuModel=="windows" & 
                lsMuModel$lsMuModelGlobal$boolVecWindowsAsBFGS){
        # The dispersion parameter co-estimation couples all LL
        # terms of a gene so that there is no computational 
        # advantage in not estimating all mean parameters
        # simultaneously (boolVecWindowsAsBFGS=TRUE).
        # Decompress parameters
        vecMuParam <- decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,],
                                            lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                                            vecInterval=NULL )
        scaDispParam <- lsDispModel$matDispModel[i,]
        
        fitDispMu <- fitDispConstMuVecWindowsZINB(
          vecCounts=matCountsProc[i,],
          vecMu=vecMuParam, # Mu model
          lsvecBatchParamGuess=lsvecBatchParamGuess,
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          scaDispGuess=scaDispParam, # Disp model
          matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius,
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp  )
        return(fitDispMu)
      } else if(lsMuModel$lsMuModelGlobal$strMuModel=="clusters"){
        # Estimate mean parameter by cluster. No smoothing is used.
        # Estimate mean parameters
        fitDispMu <-fitDispConstMuClusterZINB(
          vecCounts=matCountsProc[i,],
          vecMuGuess=lsMuModel$matMuModel[i,], # Mu model
          lsvecBatchParamGuess=lsvecBatchParamGuess,
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          scaDispGuess=lsDispModel$matDispModel[i,], # Disp model
          matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
        
        return(fitDispMu)
      } else if(lsMuModel$lsMuModelGlobal$strMuModel=="MM"){
        # Estimate mean parameter by cluster. No smoothing is used.
        
        # Estimate mean parameters
        fitDispMu <-fitDispConstMuMMZINB(
          vecCounts=matCountsProc[i,],
          vecMuGuess=lsMuModel$matMuModel[i,], # Mu model
          lsvecBatchParamGuess=lsvecBatchParamGuess,
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          scaDispGuess=lsDispModel$matDispModel[i,], # Disp model
          matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          matWeights=matWeights,
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
        
        return(fitDispMu)
      } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){      
        # Decompress parameters
        vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
                                             lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                                             vecInterval=NULL)
        
        # Estimate mean parameters
        fitDispMu <- fitDispConstMuImpulseZINB(
          vecCounts=matCountsProc[i,],
          vecImpulseParamGuess=lsMuModel$matMuModel[i,], # Mu model
          lsvecBatchParamGuess=lsvecBatchParamGuess,
          lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
          scaDispGuess=lsDispModel$matDispModel[i,], # Disp model
          matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
          vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
          scaWindowRadius=scaWindowRadius,
          MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
          RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
        return(fitDispMu)
      } else {
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
  })
  
  matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$scaDisp))
  matMuModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$vecMuModel))
  lsmatBatchModel <- lapply(seq(1,length(lsMuModel$lsmatBatchModel)), function(confounder){
    do.call(rbind, lapply(lsFitDispMu,  function(i) i$lsvecBatchFactors[[confounder]] ))
  })
  vecConvergence <- sapply(lsFitDispMu,  function(i) i$scaConvergence)
  vecLL <- sapply(lsFitDispMu,  function(i) i$scaLL)
  
  return( list(matMuModel=matMuModel,
               lsmatBatchModel=lsmatBatchModel,
               matDispModel=matDispModel,
               vecConvergence=vecConvergence,
               vecLL=vecLL) )
}