#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++     Fit dropout parameters of ZINB model    +++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of logistic drop-out paramater model on single gene given
#' the negative binomial mean and dispersion parameters.
#'
#' @seealso Called by fitting wrapper:
#' \code{fitPiZINB}.
#' Calls \code{evalLogLikZINB}.
#' Compiled version: \link{evalLogLikPiZINB_comp}.
#' 
#' @param vecTheta: (numeric vector length linear model) 
#'    Parameter estimates for logit linear model for drop-out rate.
#' @param vecCounts (count vector number of genes)
#'    Observed read counts, not observed are NA.
#' @param vecMu: (vector number of cells) Negative binomial
#'    mean parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param matPiPredictors: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikPiZINB_SingleCell <- function(
  vecTheta,
  vecCounts,
  vecMu,
  scaNormConst,
  vecDisp,
  matPiAllPredictors,
  strDropModel,
  vecboolNotZero,
  vecboolZero ){ 
  
  # (I) Linker functions
  # Force mean parameter to be negative
  if(strDropModel=="logistic_ofMu"){
    vecTheta[2] <- -exp(vecTheta[2])
  }
  
  # (II) Prevent parameter shrinkage/explosion
  vecTheta[vecTheta < -10^(10)] <- -10^(10)
  vecTheta[vecTheta > 10^(10)] <- 10^(10)
  if(strDropModel=="logistic_ofMu"){
    if(vecTheta[2] > -10^(-10)) vecTheta[2] <- -10^(-10)
  }
  
  vecPiEst <- sapply(seq(1,dim(matPiAllPredictors)[1]), function(i){
    evalDropoutModel_comp(
      vecPiModel=vecTheta, 
      vecPiPredictors=matPiAllPredictors[i,] )
  })
  
  # (III) Evaluate loglikelihood of estimate
  scaLogLik <- evalLogLikZINB_comp(
    vecCounts=vecCounts,
    vecMu=vecMu*scaNormConst,
    vecDisp=vecDisp, 
    vecPi=vecPiEst,
    vecboolNotZero=vecboolNotZero, 
    vecboolZero=vecboolZero )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

evalLogLikPiZINB_SingleCell_comp <- cmpfun(evalLogLikPiZINB_SingleCell)

# Adapted for memory usage
evalLogLikPiZINB_ManyCells <- function(
  vecTheta,
  matCounts,
  lsMuModel,
  lsDispModel,
  lsDropModel){ 
  
  scaNumGenes <- dim(matCounts)[1]
  # (I) Linker functions
  # Force mean parameter to be negative
  if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    vecTheta[2] <- -exp(vecTheta[2])
  }
  
  # (II) Prevent parameter shrinkage/explosion
  vecTheta[vecTheta < -10^(10)] <- -10^(10)
  vecTheta[vecTheta > 10^(10)] <- 10^(10)
  if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    if(vecTheta[2] > -10^(-10)) vecTheta[2] <- -10^(-10)
  }
  lsDropModel$matDropoutLinModel <- matrix(
    vecTheta,
    nrow=lsDropModel$lsDropModelGlobal$scaNumCells,
    ncol=length(vecTheta),
    byrow = TRUE)
  
  # (III) Evaluate loglikelihood of estimate
  # Have to decompress parameters in loop to keep memory usage down.
  # This slows as this could be pre-computed outside of optim!
  scaLogLik <- sum(evalLogLikMatrix(
    matCounts=matCounts,
    lsMuModel=lsMuModel,
    lsDispModel=lsDispModel, 
    lsDropModel=lsDropModel,
    matWeights=NULL ))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}


#' Compiled function: evalLogLikPiZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikPiZINB}.
#' 
#' @seealso \link{evalLogLikPiZINB}
#' 
#' @param vecTheta: (numeric vector length linear model) 
#'    Parameter estimates for logit linear model for drop-out rate.
#' @param vecCounts (count vector number of genes)
#'    Observed read counts, not observed are NA.
#' @param vecMu: (vector number of cells) Negative binomial
#'    mean parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param matPiPredictors: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikPiZINB_ManyCells_comp <- cmpfun(evalLogLikPiZINB_ManyCells)

#' Optimisation function for drop-out model fitting
#' 
#' This function fits a logistic drop-out model to a cell based
#' on given gene-specific predictors (which enter the linear model).
#' Parameter estimation of the linear model is performed by maximum
#' likelihood based on the overall likelihood.
#'
#' @seealso Called by drop-out estimation wrapper code in \code{fitZINB}.
#' Calls loglikelihood wrapper inside of optim:
#' \code{evalLogLikPiZINB}.
#' 
#' @param vecDropoutLinModel: (numeric vector length linear model)
#'    Previous parameterisation of linear model for drop-out 
#'    rate in logit space for given cell.
#' @param matPiConstPredictors: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecCounts (count vector number of genes)
#'    Observed read counts, not observed are NA.
#' @param lsMuModel:
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecidxInterval: (integer vector neighbourhood)
#'    Positions of cells within smooting interval (neighbourhood)
#'    of target cell.
#' @param scaTarget: (integer) Index of target cell in interval.
#' 
#' @return vecLinModel: (numeric vector length linear model) 
#'    Linear model for drop-out rate in logit space for given cell.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
fitPi_SingleCell <- function(
  vecParamGuess,
  vecCounts,
  vecMuParam,
  scaNormConst,
  vecDispParam,
  matPiAllPredictors,
  lsDropModelGlobal ){ 
  
  # (I) Numerical maximum likelihood estimation of linear model
  lsLinModelFit <- tryCatch({
    optim(
      par=vecParamGuess,
      evalLogLikPiZINB_SingleCell_comp,
      matPiAllPredictors=matPiAllPredictors,
      vecCounts=vecCounts,
      vecMu=vecMuParam,
      vecDisp=vecDispParam,
      scaNormConst=scaNormConst,
      strDropModel=lsDropModelGlobal$strDropModel,
      vecboolNotZero=!is.na(vecCounts) & vecCounts>0,
      vecboolZero= !is.na(vecCounts) & vecCounts==0,
      method="BFGS",
      control=list(maxit=lsDropModelGlobal$MAXIT_BFGS_Pi,
                   reltol=lsDropModelGlobal$RELTOL_BFGS_Pi,
                   fnscale=-1)
    )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting logistic drop-out model: fitPi_SingleCell()."))
    print(strErrorMsg)
    stop()
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  vecLinModel <- unlist(lsLinModelFit[["par"]])
  if(lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    vecLinModel[2] <- -exp(vecLinModel[2])
  }
  # # Catch boundary of likelihood domain on parameter space
  vecLinModel[vecLinModel < -10^(10)] <- -10^(10)
  vecLinModel[vecLinModel > 10^(10)] <- 10^(10)
  if(lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    if(vecLinModel[2] > -10^(-10)) vecLinModel[2] <- -10^(-10)
  }
  
  scaConvergence <- unlist(lsLinModelFit[["convergence"]])
  scaLL <- unlist(lsLinModelFit[["value"]])
  
  return( list(vecLinModel=vecLinModel,
               scaConvergence=scaConvergence,
               scaLL=scaLL) )
}

fitPi_ManyCells <- function(
  vecParamGuess,
  matCounts,
  lsMuModel,
  lsDispModel,
  lsDropModel){ 
  
  # (I) Numerical maximum likelihood estimation of linear model
  lsLinModelFit <- tryCatch({
    optim(
      par=vecParamGuess,
      evalLogLikPiZINB_ManyCells_comp,
      matCounts=matCounts,
      lsMuModel=lsMuModel,
      lsDispModel=lsDispModel,
      lsDropModel=lsDropModel,
      method="BFGS",
      control=list(maxit=lsDropModel$lsDropModelGlobal$MAXIT_BFGS_Pi,
                   reltol=lsDropModel$lsDropModelGlobal$RELTOL_BFGS_Pi,
                   fnscale=-1)
    )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting logistic drop-out model: fitPi_ManyCells()."))
    print(strErrorMsg)
    print(paste0("vecParamGuess ", paste(vecParamGuess, collapse=" ")))
    stop()
  })
  
  # (II) Extract results and correct for sensitivity boundaries
  vecLinModel <- unlist(lsLinModelFit[["par"]])
  if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    vecLinModel[2] <- -exp(vecLinModel[2])
  }
  # # Catch boundary of likelihood domain on parameter space
  vecLinModel[vecLinModel < -10^(10)] <- -10^(10)
  vecLinModel[vecLinModel > 10^(10)] <- 10^(10)
  if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
    if(vecLinModel[2] > -10^(-10)) vecLinModel[2] <- -10^(-10)
  }
  
  scaConvergence <- unlist(lsLinModelFit[["convergence"]])
  scaLL <- unlist(lsLinModelFit[["value"]])
  
  return( list(vecLinModel=vecLinModel,
               scaConvergence=scaConvergence,
               scaLL=scaLL) )
}

fitZINBPi <- function(matCounts,
                      lsMuModel,
                      lsDispModel,
                      lsDropModel){
  
  scaNumGenes <- dim(matCounts)[1]
  scaNumCells <- dim(matCounts)[2]
  if(lsDropModel$lsDropModelGlobal$strDropFitGroup=="PerCell"){
    lsFitsPi <- bplapply(seq(1,scaNumCells), function(j){
      # Decompress parameters
      vecMuParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
        decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,],
                              lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
                              lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                              vecInterval=j)
      }))
      vecDispParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
        decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
                             lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                             vecInterval=j)
      }))
      if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
        matPiAllPredictors<-cbind(rep(1,scaNumGenes),
                                  log(vecMuParam),
                                  lsDropModel$matPiConstPredictors)
      } else if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic"){
        matPiAllPredictors<-cbind(rep(1,scaNumGenes),
                                  lsDropModel$matPiConstPredictors)
      } else {
        stop(paste0("ERROR fitZINBPi: ",
                    "lsDropModel$lsDropModelGlobal$strDropModel=",
                    lsDropModel$lsDropModelGlobal$strDropModel,
                    " not recognised."))
      }
      # Initialise optimisation of linear model:
      vecParamGuess <- lsDropModel$matDropoutLinModel[j,]
      if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
        vecParamGuess[2] <- log(-vecParamGuess[2])
      }
      lsFitPi <- fitPi_SingleCell(
        vecParamGuess=vecParamGuess,
        vecCounts=matCounts[,as.double(j)],
        vecMuParam=vecMuParam,
        scaNormConst=lsMuModel$lsMuModelGlobal$vecNormConst[j],
        vecDispParam=vecDispParam,
        matPiAllPredictors=matPiAllPredictors,
        lsDropModelGlobal=lsDropModel$lsDropModelGlobal )
      return(lsFitPi)
    })
    matDropoutLinModel <- do.call(rbind, lapply(lsFitsPi, function(cell) cell$vecLinModel))
    vecboolConverged <- sapply(lsFitsPi, function(cell) cell$scaConvergence)
    vecLL <- sapply(lsFitsPi, function(cell) cell$scaLL) 
  } else if(lsDropModel$lsDropModelGlobal$strDropFitGroup=="ForAllCells"){
    # Initialise optimisation of linear model:
    vecParamGuess <- lsDropModel$matDropoutLinModel[1,] # All rows are the same!
    if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
      vecParamGuess[2] <- log(-vecParamGuess[2])
    }
    
    lsFitPi <- fitPi_ManyCells(
      vecParamGuess=vecParamGuess,
      matCounts=matCounts,
      lsMuModel=lsMuModel,
      lsDispModel=lsDispModel,
      lsDropModel=lsDropModel)
    matDropoutLinModel <- matrix(lsFitPi$vecLinModel,
                                 nrow=scaNumCells,
                                 ncol=length(lsFitPi$vecLinModel),
                                 byrow = TRUE)
    vecboolConverged <- lsFitPi$scaConvergence
    vecLL <- lsFitPi$scaLL
  } else {
    stop(paste0("ERROR fitZINBPi: ",
                "lsDropModel$lsDropModelGlobal$strDropFitGroup=",
                lsDropModel$lsDropModelGlobal$strDropFitGroup,
                " not recognised."))
  }
  
  return(list(
    matDropoutLinModel = matDropoutLinModel,
    vecboolConverged   = vecboolConverged,
    vecLL              = vecLL 
  ))
}