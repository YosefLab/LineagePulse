#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++     Fit dropout parameters of ZINB model    +++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# File divided into:
# ++++ Model cost functions -evalLogLikZINB (not in this file)
# +++ Fitting cost functions: return loglik of proposed parameters
# +++ - evalLogLikPiZINB_SingleCell fast if model is fit for single cell
# +++ - evalLogLikPiZINB_ManyCells memory efficient on many cells
# ++ Optim wrappers for single fits
# ++ - fitPi_SingleCell wrapper for fit on single cell
# ++ - fitPi_ManyCells wrapper for fit on multiple cells
# + Overall model fitting wrapper: Top level auxillary function called by fitZINB.
# + - fitPi wrapper for all drop-out models to be fit
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (I) Fitting cost functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of logistic drop-out paramater model on single gene given
#' the negative binomial mean and dispersion parameters.
#'
#' @seealso Called by fitting wrapper: \code{fitPi_SingleCell}.
#' Calls \code{evalLogLikZINB}.
#' Compiled version: \link{evalLogLikPiZINB_SingleCell_comp}.
#' 
#' @param vecTheta (numeric vector length linear model) 
#' Parameter estimates for logit linear model for drop-out rate.
#' @param vecCounts (count vector number of genes)
#' Observed read counts, not observed are NA.
#' @param vecMu (vector number of cells) Negative binomial
#' mean parameter estimate.
#' @param scaNormConst (scalar) 
#' Model scaling factors, one per cell.
#' @param vecDisp (vector number of cells) Negative binomial
#' dispersion parameter estimate.
#' @param matPiAllPredictors (matrix genes x predictors) Predictors of
#' the drop-out rate in the linear model. Minimum are a constant
#' offset and log of the negative binomial mean parameter. 
#' Other gene-specific predictors can be added.
#' @param strDropModel (str) {"logistic_ofMu", "logistic"}
#' [Default "logistic_ofMu"] Definition of drop-out model.
#' "logistic_ofMu" - include the fitted mean in the linear model
#' of the drop-out rate and use offset and matPiConstPredictors.
#' "logistic" - only use offset and matPiConstPredictors.
#' @param vecidxNotZero (bool vector number of cells)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of cells)
#' Whether observation is zero.
#' 
#' @return scaLogLik (scalar) Value of cost function:
#' zero-inflated negative binomial likelihood.
#' 
#' @author David Sebastian Fischer
evalLogLikPiZINB_SingleCell <- function(
    vecTheta,
    vecCounts,
    vecMu,
    scaNormConst,
    vecDisp,
    matPiAllPredictors,
    strDropModel,
    vecidxNotZero,
    vecidxZero){ 
    
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
    
    vecPiEst <- sapply(seq_len(nrow(matPiAllPredictors)), function(i){
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
        vecidxNotZero=vecidxNotZero, 
        vecidxZero=vecidxZero )
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' 
#' Refer to evalLogLikPiZINB_SingleCell().
#'
#' @seealso Compiled version of \code{evalLogLikPiZINB_SingleCell}.
#' 
#' @param vecTheta (numeric vector length linear model) 
#' Parameter estimates for logit linear model for drop-out rate.
#' @param vecCounts (count vector number of genes)
#' Observed read counts, not observed are NA.
#' @param vecMu (vector number of cells) Negative binomial
#' mean parameter estimate.
#' @param scaNormConst (scalar) 
#' Model scaling factors, one per cell.
#' @param vecDisp (vector number of cells) Negative binomial
#' dispersion parameter estimate.
#' @param matPiAllPredictors (matrix genes x predictors) Predictors of
#' the drop-out rate in the linear model. Minimum are a constant
#' offset and log of the negative binomial mean parameter. 
#' Other gene-specific predictors can be added.
#' @param strDropModel (str) {"logistic_ofMu", "logistic"}
#' [Default "logistic_ofMu"] Definition of drop-out model.
#' "logistic_ofMu" - include the fitted mean in the linear model
#' of the drop-out rate and use offset and matPiConstPredictors.
#' "logistic" - only use offset and matPiConstPredictors.
#' @param vecidxNotZero (bool vector number of cells)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of cells)
#' Whether observation is zero.
#' 
#' @return scaLogLik (scalar) Value of cost function:
#' zero-inflated negative binomial likelihood.
#' 
#' @author David Sebastian Fischer
evalLogLikPiZINB_SingleCell_comp <- compiler::cmpfun(evalLogLikPiZINB_SingleCell)

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' for many cells
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of logistic drop-out paramater model on multiple cells given
#' the negative binomial mean and dispersion parameters.
#' This function is optimised for memory usage vs evalLogLikPiZINB_SingleCell
#' at the cost of computation time: The parameter models are not
#' kept as a gene x cell matrix but as the raw model objects which
#' are evaluated within the cost function.
#'
#' @seealso Called by fitting wrapper: \code{fitPi_ManyCells}.
#' Calls \code{evalLogLikMatrix}.
#' Compiled version: \link{evalLogLikPiZINB_ManyCells_comp}.
#' 
#' @param vecTheta (numeric vector length linear model) 
#' Parameter estimates for logit linear model for drop-out rate.
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' 
#' @return scaLogLik (scalar) Value of cost function:
#' zero-inflated negative binomial likelihood.
#' 
#' @author David Sebastian Fischer
evalLogLikPiZINB_ManyCells <- function(
    vecTheta,
    matCounts,
    lsMuModel,
    lsDispModel,
    lsDropModel){ 
    
    scaNumGenes <- nrow(matCounts)
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

#' Cost function zero-inflated negative binomial model for drop-out fitting
#' 
#' Refer to \link{evalLogLikPiZINB_ManyCells}.
#' 
#' @seealso Compiled version of \link{evalLogLikPiZINB_ManyCells}
#' 
#' @param vecTheta (numeric vector length linear model) 
#' Parameter estimates for logit linear model for drop-out rate.
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' 
#' @return scaLogLik (scalar) Value of cost function:
#' zero-inflated negative binomial likelihood.
#' 
#' @author David Sebastian Fischer
evalLogLikPiZINB_ManyCells_comp <- compiler::cmpfun(evalLogLikPiZINB_ManyCells)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) Optin wrappers
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Optim wrapper for drop-out model fitting on single cell
#' 
#' This function fits a logistic drop-out model to a cell based
#' on given gene-specific predictors (which enter the linear model).
#' Parameter estimation of the linear model is performed by maximum
#' likelihood based on the overall likelihood.
#'
#' @seealso Called by drop-out estimation wrapper code in \code{fitPi}.
#' Calls fitting cost function:
#' \code{evalLogLikPiZINB_SingleCell_comp}.
#' 
#' @param vecParamGuess (numeric vector length linear model) 
#' Parameter estimates for logit linear model for drop-out rate.
#' @param vecCounts (count vector number of genes)
#' Observed read counts, not observed are NA.
#' @param vecMuParam (vector number of cells) Negative binomial
#' mean parameter estimate.
#' @param scaNormConst (scalar) 
#' Model scaling factors, one per cell.
#' @param vecDispParam (vector number of cells) Negative binomial
#' dispersion parameter estimate.
#' @param matPiAllPredictors (matrix genes x predictors) Predictors of
#' the drop-out rate in the linear model. Minimum are a constant
#' offset and log of the negative binomial mean parameter. 
#' Other gene-specific predictors can be added.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' 
#' @return vecLinModel (numeric vector length linear model) 
#' Linear model for drop-out rate in logit space for given cell.
#' 
#' @author David Sebastian Fischer
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
            vecidxNotZero= which(!is.na(vecCounts) & vecCounts>0),
            vecidxZero= which(!is.na(vecCounts) & vecCounts==0),
            method="BFGS",
            control=list(maxit=lsDropModelGlobal$MAXIT_BFGS_Pi,
                         reltol=lsDropModelGlobal$RELTOL_BFGS_Pi,
                         fnscale=-1)
        )[c("par","value","convergence")]
    }, error=function(strErrorMsg){
        message(paste0("ERROR: Fitting logistic drop-out model: fitPi_SingleCell()."))
        message(strErrorMsg)
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

#' Optim wrapper for drop-out model fitting on many cells
#' 
#' This function fits a logistic drop-out model to a set of cells based
#' on given gene-specific predictors (which enter the linear model).
#' Parameter estimation of the linear model is performed by maximum
#' likelihood based on the overall likelihood.
#'
#' @seealso Called by drop-out estimation wrapper code in \code{fitPi}.
#' Calls fitting cost function:
#' \code{evalLogLikPiZINB_ManyCells_comp}.
#' 
#' @param vecParamGuess (numeric vector length linear model) 
#' Initial parameter estimates for logit linear model for drop-out rate.
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' 
#' @return vecLinModel (numeric vector length linear model) 
#' Linear model for drop-out rate in logit space for given cell.
#' 
#' @author David Sebastian Fischer
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
        message(paste0("ERROR: Fitting logistic drop-out model: fitPi_ManyCells()."))
        message(strErrorMsg)
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Overall model fitting wrapper
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Global wrapper for fitting of all drop-out models
#' 
#' This function fits all drop-out models required on this data set.
#'
#' @seealso Called by ZINB fitting wrapper \code{fitModel}.
#' Calls optim wrapper \code{fitPi_ManyCells} or \code{fitPi_SingleCell}.
#' 
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' 
#' @return vecLinModel (numeric vector length linear model) 
#' Linear model for drop-out rate in logit space for given cell.
#' 
#' @author David Sebastian Fischer
fitPi <- function(
    matCounts,
    lsMuModel,
    lsDispModel,
    lsDropModel){
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    if(lsDropModel$lsDropModelGlobal$strDropFitGroup=="PerCell"){
        lsFitsPi <- bplapply(seq_len(scaNumCells), function(j){
            # Decompress parameters
            if(lsMuModel$lsMuModelGlobal$strMuModel=="MM") {
                stop("ERROR: fitZINBPi not adapted to strMuModel=MM yet. Talk to David.")
            }
            vecMuParam <- do.call(c, lapply(seq_len(scaNumGenes), function(i){
                decompressMeansByGene(
                    vecMuModel=lsMuModel$matMuModel[i,],
                    lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
                    lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                    vecInterval=j)
            }))
            if(lsDispModel$lsDispModelGlobal$strDispModel=="MM") {
                stop("ERROR: fitZINBPi not adapted to strDispModel=MM yet. Talk to David.")
            }
            vecDispParam <- do.call(c, lapply(seq_len(scaNumGenes), function(i){
                decompressDispByGene(
                    vecDispModel=lsDispModel$matDispModel[i,],
                    lsvecBatchModel=lapply(lsDispModel$lsmatBatchModel, function(mat) mat[i,] ),
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
                stop("ERROR fitZINBPi: ",
                     "lsDropModel$lsDropModelGlobal$strDropModel=",
                     lsDropModel$lsDropModelGlobal$strDropModel,
                     " not recognised.")
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
        vecConverged <- sapply(lsFitsPi, function(cell) cell$scaConvergence)
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
        matDropoutLinModel <- matrix(
            lsFitPi$vecLinModel,
            nrow=scaNumCells,
            ncol=length(lsFitPi$vecLinModel),
            byrow = TRUE)
        vecConverged <- lsFitPi$scaConvergence
        vecLL <- lsFitPi$scaLL
    } else {
        stop("ERROR fitZINBPi: ",
             "lsDropModel$lsDropModelGlobal$strDropFitGroup=",
             lsDropModel$lsDropModelGlobal$strDropFitGroup,
             " not recognised.")
    }
    
    # Name model matrix dimensions
     rownames(matDropoutLinModel) <- colnames(matCounts)
    
    return(list(
        matDropoutLinModel = matDropoutLinModel,
        vecConverged       = vecConverged,
        vecLL              = vecLL 
    ))
}