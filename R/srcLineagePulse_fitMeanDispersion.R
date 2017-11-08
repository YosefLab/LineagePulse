#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++     Co-Fit mean and dispersion parameters of ZINB model    +++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# File divided into:
# ++++ Model cost functions -evalLogLikZINB (not in this file)
# +++ Fitting cost functions: return loglik of proposed parameters
# +++ - evalLogLikContinuousZINB 
# ++ Optim wrappers for single fits
# ++ - fitContinuousZINB wrapper for single initialisation on one gene
# ++ - fitImpulseZINB wrapper for multiple initialisation of impulse model
# + Overall model fitting wrapper: 
# + Top level auxillary function called by fitZINB.
# + - fitMeanDisp wrapper for all mean and dispersion models to be fit
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (I) Fitting cost function
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function (zero-inflated) negative binomial model for mean and
#' dispersion model fitting
#' 
#' Log likelihood of (zero inflated) negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of mean and dispersion paramater model on single gene given
#' the drop-out model.
#' 
#' @param vecTheta (numeric vector dispersion (1) and impulse parameters (6))
#' Dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'Observed read counts, not observed are NA.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param lsDispModelGlobal (list)
#' Object containing meta-data of gene-wise dispersion parameter models.
#' @param vecTimepoints (numerical vector number of unique time coordinates)
#' Unique (pseudo)time coordinates of cells.
#' @param vecindTimepointAssign (numeric vector number samples) 
#' Index of time point assigned to cell in list of sorted
#' time points. 
#' @param matDropoutLinModel (matrix number of cells x number of predictors)
#' Logistic linear model parameters of the dropout rate 
#' as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors 
#' (numeric vector constant gene-wise coefficients)
#' Constant gene-wise coeffiecients, i.e. predictors which are not
#' the offset and not the mean parameter.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' @param vecPiParam (numeric vector number of observations)
#' Pre-evaluated drop-out model if model is not a function on the mean
#' parameter to be fit.
#' @param vecidxNotZero (index vector vector number of cells)
#' Observation which are larger than zero.
#' @param vecidxZero (index vector number of cells)
#' Observation which are zero.
#'
#' @return scaLogLik (scalar) 
#' Value of cost function (likelihood) for given gene.
#'
#' @author David Sebastian Fischer
evalLogLikMuDispGeneFit <- function(
    vecTheta,
    vecCounts,
    lsMuModelGlobal,
    lsDispModelGlobal,
    vecTimepoints,
    vecindTimepointAssign,
    matDropoutLinModel,
    vecPiConstPredictors,
    lsDropModelGlobal,
    vecPiParam=NULL,
    vecidxNotZero, 
    vecidxZero){  
    
    scaNParamUsed <- 0
    # Disp
    if(lsDispModelGlobal$strDispModel == "constant") {
        vecDispModel <- exp(vecTheta[scaNParamUsed + 1])
        scaNParamUsed <- length(vecDispModel)
        
        vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10) 
        vecDispModel[vecDispModel > 10^(10)] <- 10^(10) 
        
        vecDispParam <- rep(vecDispModel, lsDispModelGlobal$scaNumCells)
    } else if(lsDispModelGlobal$strDispModel == "groups") {
        vecDispModel <- exp(vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsDispModelGlobal$scaNGroups,
                by = 1)])
        scaNParamUsed <- length(vecDispModel)
        
        vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10) 
        vecDispModel[vecDispModel > 10^(10)] <- 10^(10)
        
        vecDispParam <- vecDispModel[lsDispModelGlobal$vecidxGroups]
    } else if(lsDispModelGlobal$strDispModel == "splines") {
        vecDispModel <- vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsDispModelGlobal$scaNSplines,
                by = 1)]
        scaNParamUsed <- length(vecDispModel)
        
        vecDispModel[vecDispModel < -10*log(10)] <- -10*log(10) 
        vecDispModel[vecDispModel > 10*log(10)] <- 10*log(10)
        
        vecDispParam <- exp(as.vector(lsDispModelGlobal$matSplineBasis %*% 
                                          vecDispModel))
    } else {
        stop("evalLogLikMuDispGeneFit()")
    }
    
    # Extract batch factors
    if(!is.null(lsDispModelGlobal$vecConfounders)){
        for(b in seq_len(lsDispModelGlobal$scaNConfounders)){
            vecBatchParamDisp <- c( 1, exp(vecTheta[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )[ lsDispModelGlobal$lsvecidxBatchAssign[[b]] ]
            scaNParamUsed <- scaNParamUsed + 
                lsDispModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchParamDisp[vecBatchParamDisp < 10^(-10)] <- 10^(-10)
            vecBatchParamDisp[vecBatchParamDisp > 10^(10)] <- 10^(10)
            
            vecDispParam <- vecDispParam*vecBatchParamDisp
        }
    }
    
    # Mu
    if(lsMuModelGlobal$strMuModel == "constant") {
        vecMuModel <- exp(vecTheta[scaNParamUsed + 1])
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
        
        vecMuParam <- rep(vecMuModel, lsMuModelGlobal$scaNumCells)
    } else if(lsMuModelGlobal$strMuModel == "groups") {
        vecMuModel <- exp(vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsMuModelGlobal$scaNGroups,
                by = 1)])
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
        
        vecMuParam <- vecMuModel[lsMuModelGlobal$vecidxGroups]
    } else if(lsMuModelGlobal$strMuModel == "splines") {
        vecMuModel <- vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsMuModelGlobal$scaNSplines,
                by = 1)]
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < -10^(10)] <- -10^(10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
        
        vecMuParam <- exp(as.vector(lsMuModelGlobal$matSplineBasis %*% 
                                        vecMuModel))
    } else if(lsMuModelGlobal$strMuModel == "impulse") {
        vecImpulseParam <- vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + 7,
                by = 1)]
        vecImpulseParam[3:5] <- exp(vecImpulseParam[3:5]) 
        # Log linker for amplitudes
        scaNParamUsed <- scaNParamUsed + length(vecImpulseParam)
        
        vecImpulseParam[1:2][vecImpulseParam[1:2] < 10^(-10)] <- 10^(-10)
        vecImpulseParam[3:5][vecImpulseParam[3:5] < 10^(-10)] <- 10^(-10)
        vecImpulseParam[3:5][vecImpulseParam[3:5] > 10^(10)] <- 10^(10)
        
        vecMuParam <- evalImpulseModel_comp(
            vecImpulseParam=vecImpulseParam,
            vecTimepoints=lsMuModelGlobal$vecContinuousCovar)
    } else {
        stop("evalLogLikImpulseZINB()")
    }
    
    # Extract batch factors
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
            scaNBatchFactors <- max(vecidxBatchAssign)-1 
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchParam <- c(1, exp(vecTheta[
                seq(from = scaNParamUsed+1,
                    to = scaNParamUsed+scaNBatchFactors,
                    by = 1)]))[vecidxBatchAssign]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchParam[vecBatchParam < 10^(-10)] <- 10^(-10)
            vecBatchParam[vecBatchParam > 10^(10)] <- 10^(10)
            
            vecMuParam <- vecMuParam*vecBatchParam
        }
    }
    
    # (III) Compute drop-out rates
    if(is.null(vecPiParam) & lsDropModelGlobal$strDropModel != "none"){
        vecPiParam <- decompressDropoutRateByGene(
            matDropModel=matDropoutLinModel,
            vecMu=vecMuParam,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal)
    }
    
    # (IV) Evaluate loglikelihood of estimate
    scaLogLik <- evalLogLikGene(
        vecCounts=vecCounts,
        vecMu=vecMuParam,
        vecNormConst=lsMuModelGlobal$vecNormConst,
        vecDisp=vecDispParam, 
        vecPi=vecPiParam,
        vecidxNotZero=vecidxNotZero, 
        vecidxZero=vecidxZero )
    
    return(scaLogLik)
}

#' Compiled function: evalLogLikMuDispGeneFit
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikMuDispGeneFit}.
#' 
#' @seealso \link{evalLogLikMuDispGeneFit}
#' 
#' @param vecTheta (numeric vector dispersion (1) and impulse parameters (6)) 
#' Dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'Observed read counts, not observed are NA.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param lsDispModelGlobal (list)
#' Object containing meta-data of gene-wise dispersion parameter models.
#' @param vecTimepoints (numerical vector number of unique time coordinates)
#'Unique (pseudo)time coordinates of cells.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'Index of time point assigned to cell in list of sorted
#'time points. 
#' @param matDropoutLinModel (matrix number of cells x number of predictors)
#' Logistic linear model parameters of the dropout rate 
#' as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors 
#' (numeric vector constant gene-wise coefficients)
#' Constant gene-wise coeffiecients, i.e. predictors which are not
#' the offset and not the mean parameter.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' @param vecPiParam (numeric vector number of observations)
#' Pre-evaluated drop-out model if model is not a function on the mean
#' parameter to be fit.
#' @param vecidxNotZero (index vector vector number of cells)
#' Observation which are larger than zero.
#' @param vecidxZero (index vector number of cells)
#' Observation which are zero.
#'
#' @return scaLogLik (scalar) Value of cost function (likelihood) for given gene.
#'
#' @author David Sebastian Fischer
evalLogLikMuDispGeneFit_comp <- compiler::cmpfun(evalLogLikMuDispGeneFit)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) Optim wrappers
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Optim wrapper for gene-wise models other than mixture model.
#' 
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the supplied mean and dispersion model.
#' This is the wrapper that calls optim.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @seealso Called once for each gene by \code{fitZINBMuDisp}
#' or within wrapper \code{fitZINBImpulse} 
#' once for each initalisation of each gene.
#' Calls fitting likelihood functions: 
#' \code{evalLogLikContinuousZINB_comp}.
#' 
#' @param vecCounts (count vector number of cells)
#'Observed read counts, not observed are NA.
#' @param vecMuModelGuess 
#' (numeric vector number of mean model parameters) 
#' Initialisation for impulse model.
#' @param lsvecBatchParamGuessMu (list) 
#' Object containing initialisation for mean parameter batch correction model.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param vecDispGuess 
#' (numeric vector number of dispersion model parameters)
#' Initialisation for dispersion model.
#' @param lsvecBatchParamGuessDisp (list) 
#' Object containing initialisation for 
#' dispersion parameter batch correction model.
#' @param lsDispModelGlobal (list)
#' Object containing meta-data of gene-wise dispersion parameter models.
#' @param matDropoutLinModel 
#' (matrix number of cells x number of predictors)
#' Logistic linear model parameters of the dropout rate 
#' as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors 
#' (numeric vector constant gene-wise coefficients)
#' Constant gene-wise coeffiecients, i.e. predictors which are not
#' the offset and not the mean parameter.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' @param vecPiParam (numeric vector number of observations)
#' Pre-evaluated drop-out model if model is not a function on the mean
#' parameter to be fit.
#' 
#' @return list
#'\itemize{
#' \item vecMuModel (numeric vector number of mu model parameters)
#' Contains the mean model parameters according to the used model.
#' \item lsvecBatchModelMu (list) 
#' Fit of batch correction models for mean parameter to given gene.
#' \item vecDispModel 
#' (numeric vector number of dispersion model parameters)
#' Contains the dispersion model parameters according to the used model.
#' \item lsvecBatchModelDisp (list) 
#' Fit of batch correction models for dispersion parameter to given gene.
#' \item scaConvergence (numeric vector number of genes) 
#' Convergence status of optim for given gene.
#' \item scaLL (numeric vector number of genes) 
#' Likelihood of model fit for given gene.
#'}
#'
#' @author David Sebastian Fischer
fitMuDispGene <- function(
    vecCounts,
    vecMuModelGuess,
    lsvecBatchParamGuessMu,
    lsMuModelGlobal,
    vecDispGuess,
    lsvecBatchParamGuessDisp,
    lsDispModelGlobal,
    matDropoutLinModel,
    vecPiConstPredictors,
    lsDropModelGlobal,
    vecPiParam){
    
    if(lsDispModelGlobal$strDispModel %in% c("splines")) {
        vecTheta <- vecDispGuess # don't use log space
    } else {
        vecTheta <- log(vecDispGuess)
    }
    if(!is.null(lsvecBatchParamGuessDisp)){
        vecTheta <- c(vecTheta, log(unlist(lsvecBatchParamGuessDisp)))
    }
    if(lsMuModelGlobal$strMuModel %in% c("impulse", "splines")) {
        vecTheta <- c(vecTheta, vecMuModelGuess) # don't use log space
    } else {
        vecTheta <- c(vecTheta, log(vecMuModelGuess))
    }
    if(!is.null(lsvecBatchParamGuessMu)){
        vecTheta <- c(vecTheta, log(unlist(lsvecBatchParamGuessMu)))
    }
    
    fitDispMu <- tryCatch({
        unlist(optim(
            par=vecTheta,			
            fn=evalLogLikMuDispGeneFit_comp, 
            vecCounts=vecCounts,
            lsMuModelGlobal=lsMuModelGlobal,
            lsDispModelGlobal=lsDispModelGlobal,
            matDropoutLinModel=matDropoutLinModel,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal,
            vecPiParam=vecPiParam,
            vecidxNotZero= which(!is.na(vecCounts) & vecCounts>0),
            vecidxZero= which(!is.na(vecCounts) & vecCounts==0),
            method="BFGS", 
            control=list(maxit=lsMuModelGlobal$MAXIT, 
                         reltol=lsMuModelGlobal$RELTOL,
                         fnscale=-1)
        )[c("par","value","convergence")] )
    }, error=function(strErrorMsg){
        message(paste0("ERROR: Fitting impulse model: fitContinuousZINB()."))
        message(strErrorMsg)
        message(paste0("vecTheta ", paste0(vecTheta, collapse = " ")))
        stop(strErrorMsg)
    })
    
    # (II) Extract results and correct for sensitivity boundaries
    scaNParamUsed <- 0
    # Disp
    if(lsDispModelGlobal$strDispModel == "constant") {
        vecDispModel <- exp(fitDispMu[scaNParamUsed + 1])
        scaNParamUsed <- 1
        
        vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10)
        vecDispModel[vecDispModel > 10^(10)] <- 10^(10)
    } else if(lsDispModelGlobal$strDispModel == "groups") {
        vecDispModel <- exp(fitDispMu[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsDispModelGlobal$scaNGroups,
                by = 1)])
        scaNParamUsed <- scaNParamUsed + length(vecDispModel)
        
        vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10)
        vecDispModel[vecDispModel > 10^(10)] <- 10^(10)
    } else if(lsDispModelGlobal$strDispModel == "splines") {
        vecDispModel <- fitDispMu[seq(
            from = scaNParamUsed + 1,
            to = scaNParamUsed + lsDispModelGlobal$scaNSplines,
            by = 1)]
        scaNParamUsed <- scaNParamUsed + length(vecDispModel)
        
        vecDispModel[vecDispModel < -10*log(10)] <- -10*log(10)
        vecDispModel[vecDispModel > 10*log(10)] <- 10*log(10)
    } else {
        stop("fitImpulseOneInitZINB()")
    }
    
    # Extract batch factors
    if(!is.null(lsDispModelGlobal$vecConfounders)){
        lsvecBatchFactorsDisp <- list()
        for(b in seq_len(lsDispModelGlobal$scaNConfounders)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1,
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + 
                lsDispModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            
            lsvecBatchFactorsDisp[[b]] <- vecBatchFactors
        }
    } else { 
        lsvecBatchFactorsDisp <- NULL 
    }
    
    # Mu
    if(lsMuModelGlobal$strMuModel == "constant") {
        vecMuModel <- exp(fitDispMu[scaNParamUsed + 1])
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
    } else if(lsMuModelGlobal$strMuModel == "groups") {
        vecMuModel <- exp(fitDispMu[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + lsMuModelGlobal$scaNGroups,
                by = 1)])
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
    } else if(lsMuModelGlobal$strMuModel == "impulse") {
        vecMuModel <- fitDispMu[seq(from = scaNParamUsed + 1,
                                    to = scaNParamUsed + 7,
                                    by = 1)]
        vecMuModel[3:5] <- exp(vecMuModel[3:5])
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[1:2][vecMuModel[1:2] < 10^(-10)] <- 10^(-10)
        vecMuModel[3:5][vecMuModel[3:5] < 10^(-10)] <- 10^(-10)
        vecMuModel[3:5][vecMuModel[3:5] > 10^(10)] <- 10^(10)
    } else if(lsMuModelGlobal$strMuModel == "splines") {
        vecMuModel <- fitDispMu[seq(
            from = scaNParamUsed + 1,
            to = scaNParamUsed + lsMuModelGlobal$scaNSplines,
            by = 1)]
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < -10^(10)] <- -10^(10)
        vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
    } else {
        stop("fitMuDispGene()")
    }
    
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        lsvecBatchFactorsMu <- list()
        for(b in seq_len(lsMuModelGlobal$scaNConfounders)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1,
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + 
                lsMuModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            
            lsvecBatchFactorsMu[[b]] <- vecBatchFactors
        }
    } else { 
        lsvecBatchFactorsMu <- NULL 
    }
    
    scaLL <- fitDispMu["value"]
    scaConvergence <- fitDispMu["convergence"]
    
    return( list(vecDispModel=vecDispModel,
                 lsvecBatchFactorsDisp=lsvecBatchFactorsDisp,
                 vecMuModel=vecMuModel,
                 lsvecBatchFactorsMu=lsvecBatchFactorsMu,
                 scaLL=scaLL,
                 scaConvergence=scaConvergence) )
}

#' Multiple initilalisation wrapper for impulse mean model
#' 
#' Multiple initialisation are tried for the impulse model. Therefore,
#' this wrapper sits ontop of fitContinuousZINB() in the fitting hierarchy
#' and wraps multiple initialisations at the level of one gene.
#' 
#' @seealso Called by mean-dispersion co-estimation wrapper 
#' \code{fitZINBMuDisp}.
#' Calls optimisation wrapper \code{fitContinuousZINB} 
#' for each initialisation.
#' 
#' @param vecCounts (count vector number of cells)
#'Observed read counts, not observed are NA.
#' @param vecImpulseParamGuess 
#' (numeric vector number of impulse model parameters) 
#' Initialisation for impulse model.
#' @param lsvecBatchParamGuessMu (list) 
#' Object containing initialisation for mean parameter batch correction model.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param vecDispGuess 
#' (numeric vector number of dispersion model parameters)
#' Initialisation for dispersion model.
#' @param lsvecBatchParamGuessDisp (list) 
#' Object containing initialisation for 
#' dispersion parameter batch correction model.
#' @param lsDispModelGlobal (list)
#' Object containing meta-data of gene-wise dispersion parameter models.
#' @param matDropoutLinModel
#' (matrix number of cells x number of predictors)
#' Logistic linear model parameters of the dropout rate 
#' as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors 
#' (numeric vector constant gene-wise coefficients)
#' Constant gene-wise coeffiecients, i.e. predictors which are not
#' the offset and not the mean parameter.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' @param vecPiParam (numeric vector number of observations)
#' Pre-evaluated drop-out model if model is not a function on the mean
#' parameter to be fit.
#' 
#' @return list
#'\itemize{
#' \item vecMuModel (numeric vector number of mu model parameters)
#' Contains the mean model parameters according to the used model.
#' \item lsvecBatchModelMu (list) 
#' Fit of batch correction models for mean parameter to given gene.
#' \item vecDispModel (numeric vector number of dispersion model parameters)
#' Contains the dispersion model parameters according to the used model.
#' \item lsvecBatchModelDisp (list) 
#' Fit of batch correction models for dispersion parameter to given gene.
#' \item scaConvergence (numeric vector number of genes) 
#' Convergence status of optim for given gene.
#' \item scaLL (numeric vector number of genes) 
#' Likelihood of model fit for given gene.
#'}
#'
#' @author David Sebastian Fischer
fitMuDispGeneImpulse <- function(
    vecCounts, 
    vecImpulseParamGuess,
    lsvecBatchParamGuessMu,
    lsMuModelGlobal,
    vecDispGuess,
    lsvecBatchParamGuessDisp,
    lsDispModelGlobal,
    matDropoutLinModel,
    vecPiConstPredictors,
    lsDropModelGlobal,
    vecPiParam ){
    
    # Try new peak and valley initialisations?
    # Increases time complexity of mean estimation by factor 3
    # but seems to make a difference on simulated data.
    boolUseNewInits <- TRUE
    
    # (I) Initialise impulse model
    # Decompress parameters for initialisation
    vecMuParam <- decompressMeansByGene(
        vecMuModel=vecImpulseParamGuess,
        lsvecBatchModel=lsvecBatchParamGuessMu,
        lsMuModelGlobal=lsMuModelGlobal,
        vecInterval=NULL )
    vecDispParam <- decompressDispByGene(
        vecDispModel=vecDispGuess,
        lsvecBatchModel=lsvecBatchParamGuessDisp,
        lsDispModelGlobal=lsDispModelGlobal,
        vecInterval=NULL )
    if(is.null(vecPiParam) & lsDropModelGlobal$strDropModel != "none"){
        vecPiParamInit <- decompressDropoutRateByGene(
            matDropModel=matDropoutLinModel,
            vecMu=vecMuParam,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal)
    } else { 
        vecPiParamInit <- vecPiParam 
    }
    
    # The previous parameter estiamte is kept as a reference and
    # used as an initialisation
    # Compute initialisations for peak and valley
    lsParamGuesses <- initialiseImpulseParameters(
        vecCounts=vecCounts,
        lsMuModelGlobal=lsMuModelGlobal)
    vecParamGuessPeak <- lsParamGuesses$peak
    vecParamGuessValley <- lsParamGuesses$valley
    
    # (II) Compute new parameters
    # 1. Initialisation: Prior best fit
    vecParamGuessPrior <- vecImpulseParamGuess
    vecParamGuessPrior[3:5] <- log(vecParamGuessPrior[3:5])
    lsFitPrior <- fitMuDispGene(
        vecCounts=vecCounts,
        vecMuModelGuess=vecParamGuessPrior,
        lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
        lsMuModelGlobal=lsMuModelGlobal,
        vecDispGuess=vecDispGuess,
        lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
        lsDispModelGlobal=lsDispModelGlobal,
        matDropoutLinModel=matDropoutLinModel,
        vecPiConstPredictors=vecPiConstPredictors,
        lsDropModelGlobal=lsDropModelGlobal,
        vecPiParam=vecPiParam)
    if(boolUseNewInits){
        # 2. Initialisation: Peak
        lsFitPeak <- fitMuDispGene(
            vecCounts=vecCounts,
            vecMuModelGuess=vecParamGuessPeak,
            lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
            lsMuModelGlobal=lsMuModelGlobal,
            vecDispGuess=vecDispGuess,
            lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
            lsDispModelGlobal=lsDispModelGlobal,
            matDropoutLinModel=matDropoutLinModel,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal,
            vecPiParam=vecPiParam)
        # 3. Initialisation: Valley
        lsFitValley <- fitMuDispGene(
            vecCounts=vecCounts,
            vecMuModelGuess=vecParamGuessValley,
            lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
            lsMuModelGlobal=lsMuModelGlobal,
            vecDispGuess=vecDispGuess,
            lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
            lsDispModelGlobal=lsDispModelGlobal,
            matDropoutLinModel=matDropoutLinModel,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal,
            vecPiParam=vecPiParam)
        
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
            scaLLGuess <- evalLogLikGene(
                vecCounts=vecCounts,
                vecMu=vecMuParam*lsMuModelGlobal$vecNormConst,
                vecDisp=vecDispParam, 
                vecPi=vecPiParam,
                vecidxNotZero= !is.na(vecCounts) & vecCounts>0, 
                vecidxZero= !is.na(vecCounts) &vecCounts==0)
            indMaxLL <- match(max(vecLL, na.rm=TRUE), vecLL)
            if(vecLL[indMaxLL] < scaLLGuess){
                lsFitBest <- list(vecDispGuess=vecDispGuess,
                                  vecImpulseParam=vecImpulseParamGuess,
                                  scaConvergence=1002)
            } else {
                lsFitBest <- lsFits[[indMaxLL]]
            }
        } else {
            # If none of the three initilisations was successful:
            # Use prior paramter values
            lsFitBest <- list(vecDispGuess=vecDispGuess,
                              vecImpulseParam=vecImpulseParamGuess,
                              scaConvergence=1001)
        }
    } else {
        # Make sure new value is better than previous
        scaLLGuess <- evalLogLikGene(
            vecCounts=vecCounts,
            vecMu=vecMuParam*lsMuModelGlobal$vecNormConst,
            vecDisp=vecDispParam, 
            vecPi=vecPiParam,
            vecidxNotZero= !is.na(vecCounts) & vecCounts>0, 
            vecidxZero= !is.na(vecCounts) &vecCounts==0)
        if(lsFitPrior$scaLL < scaLLGuess){
            lsFitBest <- list(vecDispGuess=vecDispGuess,
                              vecImpulseParam=vecImpulseParamGuess,
                              scaConvergence=1002)
        } else {
            lsFitBest <- lsFitPrior
        }
    }
    
    return(list(vecDispModel=lsFitBest$vecDispModel,
                lsvecBatchFactorsDisp=lsFitBest$lsvecBatchFactorsDisp,
                vecMuModel=lsFitBest$vecMuModel,
                lsvecBatchFactorsMu=lsFitBest$lsvecBatchFactorsMu,
                scaConvergence=lsFitBest$scaConvergence,
                scaLL=lsFitBest$scaLL))
}

#' Optim wrapper for gene-wise models other than mixture model.
#' 
#' NOT YET SUPPORTED.
#' LINEAGEPULSE CODE WILL BE EXTENDED BY MODULAR FUNCTIONALITIES
#' AND THIS IS ONE INSTANCE OF A PLACEHOLDER USED FOR DEVELOPING.
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the supplied mean and dispersion model.
#' This is the wrapper that calls optim.
#' 
#' This function performs error handling of the numerical fitting procedure.
#' This function corrects for the likelihood sensitivity bounds used in the 
#' cost function.
#' 
#' @param vecCounts (count vector number of cells)
#'Observed read counts, not observed are NA.
#' @param vecMuGuess (numeric vector number of mean model parameters) 
#' Initialisation for impulse model.
#' @param lsvecBatchParamGuessMu (list) 
#' Object containing initialisation for mean parameter batch correction model.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param vecDispGuess (numeric vector number of dispersion model parameters) 
#' Initialisation for dispersion model.
#' @param lsvecBatchParamGuessDisp (list) 
#' Object containing initialisation for 
#' dispersion parameter batch correction model.
#' @param lsDispModelGlobal (list)
#' Object containing meta-data of gene-wise dispersion parameter models.
#' @param matDropoutLinModel (matrix number of cells x number of predictors)
#' Logistic linear model parameters of the dropout rate 
#' as a function of the mean and constant gene-wise coefficients.
#' @param vecPiConstPredictors 
#' (numeric vector constant gene-wise coefficients)
#' Constant gene-wise coeffiecients, i.e. predictors which are not
#' the offset and not the mean parameter.
#' @param lsDropModelGlobal (list)
#' Object containing meta-data of cell-wise drop-out parameter models.
#' @param vecPiParam (numeric vector number of observations)
#' Pre-evaluated drop-out model if model is not a function on the mean
#' parameter to be fit.
#' @param matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' @param MAXIT (scalar) 
#' Maximum number of BFGS iterations handed to optim().
#' @param RELTOL (scalar) 
#' Cost function convergence criterium handed to optim().
#' 
#' @return list
#'\itemize{
#' \item vecMuModel (numeric vector number of mu model parameters)
#' Contains the mean model parameters according to the used model.
#' \item lsvecBatchModelMu (list) 
#' Fit of batch correction models for mean parameter to given gene.
#' \item vecDispModel (numeric vector number of dispersion model parameters)
#' Contains the dispersion model parameters according to the used model.
#' \item lsvecBatchModelDisp (list) 
#' Fit of batch correction models for dispersion parameter to given gene.
#' \item scaConvergence (numeric vector number of genes) 
#' Convergence status of optim for given gene.
#' \item scaLL (numeric vector number of genes) 
#' Likelihood of model fit for given gene.
#'}
#'
#' @author David Sebastian Fischer
fitMuDispGeneMM <- function(
    vecCounts,
    vecMuGuess,
    lsvecBatchParamGuessMu,
    lsMuModelGlobal,
    vecDispGuess,
    lsvecBatchParamGuessDisp,
    lsDispModelGlobal,
    matDropoutLinModel,
    vecPiConstPredictors,
    lsDropModelGlobal,
    vecPiParam,
    matWeights,
    MAXIT,
    RELTOL) {
    
    stop("LineagePulse ERROR: fitMMZINB() not yet supported.",
         " Contact developer.")
    return(NULL)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (III) Overall model fitting wrapper
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Coordinate mean and dispersion parameter co-estimation step
#' 
#' Auxillary function that calls the estimation functions for the
#' different mean and dispersion models according to their needs. Note that one
#' function has to be coded for each combination of mean and dispersion
#' models.
#' 
#' @seealso Called by \code{fitModel}. Calls fitting wrappers:
#' \code{fitImpulseZINB},
#' \code{fitContinuousZINB}.
#' 
#' @param matCountsProc (matrix genes x cells)
#'Observed read counts, not observed are NA.
#' @param vecNormConst (numeric vector number of cells) 
#' Model scaling factors, one per cell.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' @param matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' 
#' @return list
#'\itemize{
#' \item matMuModel (numeric matrix genes x mu model parameters)
#' Contains the mean model parameters according to the used model.
#' \item lsmatBatchModelMu (list) 
#' Fit of batch correction models for mean parameter.
#' \item matDispModel (numeric matrix genes x disp model parameters)
#' Contains the dispersion model parameters according to the used model.
#' \item lsmatBatchModelDisp (list) 
#' Fit of batch correction models for dispersion parameter.
#' \item vecConvergence (numeric vector number of genes) 
#' Convergence status of optim for each gene.
#' \item vecLL (numeric vector number of genes) Likelihood of model fit.
#'}
#'
#' @author David Sebastian Fischer
fitMuDisp <- function(
    matCountsProc,
    vecNormConst,
    lsMuModel,
    lsDispModel,
    lsDropModel,
    matWeights){
    
    scaNumGenes <- nrow(matCountsProc)
    scaNumCells <- ncol(matCountsProc)
    
    # Parallelize fitting across gene.
    lsFitDispMu <- bplapply(seq_len(scaNumGenes), function(i){
        # Extract batch factors in format required for fitting.
        if(!is.null(lsMuModel$lsMuModelGlobal$vecConfounders)) {
            lsvecBatchParamGuessMu <- 
                lapply(lsMuModel$lsmatBatchModel, function(mat) {
                    mat[i,2:dim(mat)[2],drop=FALSE] 
                })
        } else {
            lsvecBatchParamGuessMu <- NULL
        }
        if(!is.null(lsDispModel$lsDispModelGlobal$vecConfounders)) {
            lsvecBatchParamGuessDisp <- 
                lapply(lsDispModel$lsmatBatchModel, function(mat) {
                    mat[i,2:dim(mat)[2],drop=FALSE] 
                })
        } else {
            lsvecBatchParamGuessDisp <- NULL
        }
        
        # Decompress drop-out rates if not function of mean
        if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
            vecPiParam <- NULL
        } else if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic") {
            vecPiParam <- decompressDropoutRateByGene(
                matDropModel=lsDropModel$matDropoutLinModel,
                vecMu=NULL,
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal)
        } else if (lsDropModel$lsDropModelGlobal$strDropModel == "none") {
            vecPiParam <- NULL
        }
        
        # Call fitting wrapper.
        if(lsMuModel$lsMuModelGlobal$strMuModel=="MM"){
            fitDispMu <- fitMuDispGeneMM(
                vecCounts=matCountsProc[as.double(i),], 
                vecMuGuess=lsMuModel$matMuModel[i,], 
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,], 
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, 
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam,
                matWeights=matWeights,
                MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
                RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
            
            return(fitDispMu)
        } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){      
            # Estimate mean parameters
            fitDispMu <- fitMuDispGeneImpulse(
                vecCounts=matCountsProc[as.double(i),], 
                vecImpulseParamGuess=lsMuModel$matMuModel[i,], 
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,],
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, 
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam )
            return(fitDispMu)
        } else if(lsMuModel$lsMuModelGlobal$strMuModel %in% 
                  c("splines", "groups", "constant")){      
            # Estimate mean parameters
            fitDispMu <- fitMuDispGene(
                vecCounts=matCountsProc[as.double(i),], 
                vecMuModelGuess=lsMuModel$matMuModel[i,],
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,], 
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, 
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam )
            return(fitDispMu)
        } else {
            #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
            message(paste0(
                "Mean parameter model not recognised for ",
                "co-estimation with dispersion: ",
                lsMuModel$lsMuModelGlobal$strMuModel, 
                ". Only constant and impulse model implemented. ",
                "Use sequential estimation."))
            stop("Mean parameter model not recognised for ",
                 " co-estimation with dispersion: ", 
                 lsMuModel$lsMuModelGlobal$strMuModel, 
                 ". Only constant and impulse model implemented. ",
                 "Use sequential estimation.")
        }  
    })
    
    # Extract model matrices from fit list
    matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) { 
        i$vecDispModel 
    }))
    matMuModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) {
        i$vecMuModel 
    }))
    if(!is.null(lsMuModel$lsMuModelGlobal$scaNConfounders)) {
        lsmatBatchModelMu <- 
            lapply(seq_len(lsMuModel$lsMuModelGlobal$scaNConfounders), 
                   function(confounder){
                       do.call(rbind, lapply(lsFitDispMu,  function(i) {
                           i$lsvecBatchFactorsMu[[confounder]] 
                       }))
                   })
    } else {
        lsmatBatchModelMu <- NULL
    }
    if(!is.null(lsDispModel$lsDispModelGlobal$scaNConfounders)) {
        lsmatBatchModelDisp <- 
            lapply(seq_len(lsDispModel$lsDispModelGlobal$scaNConfounders),
                   function(confounder){
                       do.call(rbind, lapply(lsFitDispMu,  function(i) {
                           i$lsvecBatchFactorsDisp[[confounder]] 
                       }))
                   })
    } else {
        lsmatBatchModelDisp <- NULL
    }
    
    ## name model matrices dimensions
    # rows
    rownames(matMuModel) <- rownames(matCountsProc)
    if(!is.null(lsmatBatchModelMu)) {
        for(i in seq_along(lsmatBatchModelMu)) {
            rownames(lsmatBatchModelMu[[i]]) <- rownames(matCountsProc)
        }
    }
    rownames(matDispModel) <- rownames(matCountsProc)
    if(!is.null(lsmatBatchModelDisp)) {
        for(i in seq_along(lsmatBatchModelDisp)) {
            rownames(lsmatBatchModelDisp[[i]]) <- rownames(matCountsProc)
        }
    }
    # columns
    if(lsMuModel$lsMuModelGlobal$strMuModel=="groups") {
        colnames(lsMuModel$matMuModel) <- vecGroups
    } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse") {
        colnames(lsMuModel$matMuModel) <- 
            c("beta1", "beta2", "h0", "h1", "h2", "t1", "t2")
    } else {
        colnames(lsMuModel$matMuModel) <- NULL
    }
    if(!is.null(lsmatBatchModelMu)) {
        for(i in seq_along(lsmatBatchModelMu)) {
            colnames(lsmatBatchModelMu) <- NULL
        } 
    }
    colnames(matDispModel) <- NULL 
    if(!is.null(lsmatBatchModelDisp)) {
        for(i in seq_along(lsmatBatchModelDisp)) {
            colnames(lsmatBatchModelDisp[[i]]) <- NULL 
        }
    }
    
            
    vecConvergence <- sapply(lsFitDispMu,  function(i) i$scaConvergence)
    vecLL <- sapply(lsFitDispMu,  function(i) i$scaLL)
    
    return( list(matMuModel=matMuModel,
                 lsmatBatchModelMu=lsmatBatchModelMu,
                 matDispModel=matDispModel,
                 lsmatBatchModelDisp=lsmatBatchModelDisp,
                 vecConvergence=vecConvergence,
                 vecLL=vecLL) )
}