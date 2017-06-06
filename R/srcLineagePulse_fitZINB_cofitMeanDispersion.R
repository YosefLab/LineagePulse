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
#' and mean co-estimation under mixture model
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
#' @param lsMuModelGlobal: (list) 
#'    List with global (not gene-specific) mean model parameters.
#' @param lsDispModelGlobal: (list) 
#'    List with global (not gene-specific) dispersion model parameters.
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
evalLogLikFitMMZINB <- function(
    vecTheta,
    vecCounts,
    lsMuModelGlobal,
    lsDispModelGlobal,
    matDropoutLinModel,
    vecPiConstPredictors,
    lsDropModelGlobal,
    vecPiParam=NULL,
    vecboolNotZero,
    vecboolZero,
    matWeights){ 
    
    #scaNCells <- length(vecCounts)
    #scaNMixtures <- dim(matWeights)[2]
    
    # Disp
    if(lsDispModelGlobal$strDispModel == "constant") {
        vecDispModel <- exp(vecTheta[1])
    } else if(lsDispModelGlobal$strDispModel == "MM") {
        vecDispModel <- exp(vecTheta[seq(1,lsDispModelGlobal$scaNMix,by=1)])
    }
    scaNParamUsed <- length(vecDispModel)
    
    vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10) 
    vecDispModel[vecDispModel > 10^(10)] <- 10^(10) 
    
    matDispParam <- matrix(vecDispModel, nrow=lsDispModelGlobal$scaNumCells, 
                           ncol=lsDispModelGlobal$scaNMix, byrow=TRUE)
    
    # Extract batch factors
    if(!is.null(lsDispModelGlobal$vecConfounders)){
        for(b in seq(1, lsDispModelGlobal$scaNConfounders, by=1)){
            vecBatchParamDisp <- c( 1, exp(vecTheta[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )[ lsDispModelGlobal$lsvecidxBatchAssign[[b]] ]
            scaNParamUsed <- scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchParamDisp[vecBatchParamDisp < 10^(-10)] <- 10^(-10)
            vecBatchParamDisp[vecBatchParamDisp > 10^(10)] <- 10^(10)
            
            matDispParam <- matDispParam*matrix(
                vecBatchParamDisp, nrow=lsDispModelGlobal$scaNumCells, 
                ncol=lsDispModelGlobal$scaNMix, byrow=FALSE)
        }
    }
    
    # Mu
    vecMuModel <- exp(vecTheta[seq(from = scaNParamUsed+1, 
                                   to = scaNParamUsed + lsMuModelGlobal$scaNMix, 
                                   by = 1)])
    scaNParamUsed <- scaNParamUsed+length(vecMuModel)
    
    vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
    vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
    
    matMuParam <- matrix(vecMuModel, nrow=lsMuModelGlobal$scaNumCells, 
                         ncol=lsMuModelGlobal$scaNMix, byrow=TRUE)
    
    # Extract batch factors
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        for(b in seq(1, lsMuModelGlobal$scaNConfounders, by=1)){
            #scaNBatchFactors <- lsMuModelGlobal$vecNBatches[b] - 1 #oldmax(vecidxBatchAssign)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchParamMu <- c( 1, exp(vecTheta[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )[ lsMuModelGlobal$lsvecidxBatchAssign[[b]] ]
            scaNParamUsed <- scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchParamMu[vecBatchParamMu < 10^(-10)] <- 10^(-10)
            vecBatchParamMu[vecBatchParamMu > 10^(10)] <- 10^(10)
            
            matMuParam <- matMuParam*matrix(
                vecBatchParamMu, nrow=lsMuModelGlobal$scaNumCells, 
                ncol=lsMuModelGlobal$scaNMix, byrow=FALSE)
        }
    }
    
    # Decompress parameters
    if(is.null(vecPiParam)){
        matPiParam <- do.call(cbind, lapply(seq(1,lsMuModelGlobal$scaNMix), function(m){
            decompressDropoutRateByGene(
                matDropModel=matDropoutLinModel,
                vecMu=matMuParam[,m],
                vecPiConstPredictors=vecPiConstPredictors,
                lsDropModelGlobal=lsDropModelGlobal)
        }))
    } else {
        matPiParam <- matrix(vecPiParam,
                             nrow=lsMuModelGlobal$scaNumCells,
                             ncol=lsMuModelGlobal$scaNMix,
                             byrow=FALSE)
    }
    
    scaLogLik <- evalLogLikGeneMM(
        vecCounts=vecCounts,
        matMuParam=matMuParam,
        vecNormConst=lsMuModelGlobal$vecNormConst,
        matDispParam=matDispParam,
        matDropParam=matPiParam,
        matWeights=matWeights,
        vecboolNotZero=vecboolNotZero, 
        vecboolZero=vecboolZero)
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}

#' Compiled function: evalLogLikDispMMMuMMZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikDispMMMuMMZINB}.
#' 
#' @param vecTheta: (numeric vector length 1+number of clusters) 
#'    {Log dispersion parameter, log mean parameters}
#'    Log dispersion and mean parameter estimates.
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param lsMuModelGlobal: (list) 
#'    List with global (not gene-specific) mean model parameters.
#' @param lsDispModelGlobal: (list) 
#'    List with global (not gene-specific) dispersion model parameters.
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
evalLogLikFitMMZINB_comp <- cmpfun(evalLogLikFitMMZINB)


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
#' \code{fitImpulseOneInitZINB}.
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
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikContinuousZINB <- function(
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
    vecboolNotZero, 
    vecboolZero){  
    
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
        stop("evalLogLikImpulseZINB()")
    }
    
    # Extract batch factors
    if(!is.null(lsDispModelGlobal$vecConfounders)){
        for(b in seq(1, lsDispModelGlobal$scaNConfounders, by=1)){
            vecBatchParamDisp <- c( 1, exp(vecTheta[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )[ lsDispModelGlobal$lsvecidxBatchAssign[[b]] ]
            scaNParamUsed <- scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1
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
        
        vecMuModel[vecMuModel < -10*log(10)] <- -10*log(10)
        vecMuModel[vecMuModel > 10*log(10)] <- 10*log(10)
        
        vecMuParam <- exp(as.vector(lsMuModelGlobal$matSplineBasis %*% 
                                        vecMuModel))
    } else if(lsMuModelGlobal$strMuModel == "impulse") {
        vecImpulseParam <- vecTheta[
            seq(from = scaNParamUsed + 1,
                to = scaNParamUsed + 7,
                by = 1)]
        vecImpulseParam[3:5] <- exp(vecImpulseParam[3:5]) # Log linker for amplitudes
        scaNParamUsed <- scaNParamUsed + length(vecImpulseParam)
        
        vecImpulseParam[3:5][vecImpulseParam[3:5] < 10^(-10)] <- 10^(-10)
        vecImpulseParam[3:5][vecImpulseParam[3:5] > 10^(10)] <- 10^(10)
        
        vecMuParam <- evalImpulseModel_comp(
            vecImpulseParam=vecImpulseParam,
            vecTimepoints=lsMuModelGlobal$vecPseudotime)
    } else {
        stop("evalLogLikImpulseZINB()")
    }
    
    # Extract batch factors
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        for(vecidxBatchAssign in lsMuModelGlobal$lsvecidxBatchAssign){
            scaNBatchFactors <- max(vecidxBatchAssign)-1 # Batches are counted from 1
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
    if(is.null(vecPiParam)){
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
        vecboolNotZero=vecboolNotZero, 
        vecboolZero=vecboolZero )
    
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
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikContinuousZINB_comp <- cmpfun(evalLogLikContinuousZINB)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# (II) FITTING COORDINATORS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
fitMMZINB <- function(
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
    MAXIT=1000,
    RELTOL=10^(-8) ){ 
    
    vecTheta <- log(vecDispGuess)
    if(!is.null(lsvecBatchParamGuessDisp)){
        vecTheta <- c(vecTheta, log(unlist(lsvecBatchParamGuessDisp)))
    }
    vecTheta <- c(vecTheta, log(vecMuGuess))
    if(!is.null(lsvecBatchParamGuessMu)){
        vecTheta <- c(vecTheta, log(unlist(lsvecBatchParamGuessMu)))
    }
    
    fitDispMu <- tryCatch({
        unlist(optim(
            par=vecTheta,
            evalLogLikFitMMZINB_comp,
            vecCounts=vecCounts,
            lsMuModelGlobal=lsMuModelGlobal,
            lsDispModelGlobal=lsDispModelGlobal,
            matDropoutLinModel=matDropoutLinModel,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal,
            vecPiParam=vecPiParam,
            matWeights=matWeights,
            vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
            vecboolZero= !is.na(vecCounts) & vecCounts==0,
            method="BFGS",
            control=list(maxit=MAXIT, 
                         reltol=RELTOL,
                         fnscale=-1) )[c("par","value","convergence")])
    }, error=function(strErrorMsg){
        print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMMZINB().",
                     " Wrote report into LinagePulse_lsErrorCausingGene.RData"))
        print(strErrorMsg)
        stop(strErrorMsg)
    })
    
    # (II) Extract results and correct for sensitivity boundaries
    scaNParamUsed <- 0
    # Disp
    if(lsDispModelGlobal$strDispModel == "constant") {
        vecDispModel <- exp(fitDispMu[scaNParamUsed + 1])
        scaNParamUsed <- 1
    } else if(lsDispModelGlobal$strDispModel == "MM") {
        vecDispModel <- exp(fitDispMu[seq(
            from = scaNParamUsed + 1,
            to = scaNParamUsed + lsDispModelGlobal$scaNMix,
            by=1)])
        scaNParamUsed <- lsDispModelGlobal$scaNMix
    }
    
    vecDispModel[vecDispModel < 10^(-10)] <- 10^(-10) 
    vecDispModel[vecDispModel > 10^(10)] <- 10^(10)  
    
    # Extract batch factors
    if(!is.null(lsDispModelGlobal$vecConfounders)){
        lsvecBatchFactorsDisp <- list()
        for(b in seq(1, lsDispModelGlobal$scaNConfounders, by=1)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            
            lsvecBatchFactorsDisp[[b]] <- vecBatchFactors
        }
    } else { 
        lsvecBatchFactorsDisp <- NULL 
    }
    
    # Mu
    vecMuModel <- exp(fitDispMu[seq(from = scaNParamUsed + 1,
                                    to = scaNParamUsed + lsMuModelGlobal$scaNMix,
                                    by = 1)])
    scaNParamUsed <- scaNParamUsed + lsMuModelGlobal$scaNMix
    
    vecMuModel[vecMuModel < 10^(-10)] <- 10^(-10)
    vecMuModel[vecMuModel > 10^(10)] <- 10^(10)
    
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        lsvecBatchFactorsMu <- list()
        for(b in seq(1, lsMuModelGlobal$scaNConfounders, by=1)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1
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
    
    return(list(vecDispModel=vecDispModel,
                vecMuModel=vecMuModel,
                lsvecBatchFactorsDisp=lsvecBatchFactorsDisp,
                lsvecBatchFactorsMu=lsvecBatchFactorsMu,
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
fitContinuousZINB <- function(
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
            fn=evalLogLikContinuousZINB_comp, 
            vecCounts=vecCounts,
            lsMuModelGlobal=lsMuModelGlobal,
            lsDispModelGlobal=lsDispModelGlobal,
            matDropoutLinModel=matDropoutLinModel,
            vecPiConstPredictors=vecPiConstPredictors,
            lsDropModelGlobal=lsDropModelGlobal,
            vecPiParam=vecPiParam,
            vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
            vecboolZero= !is.na(vecCounts) & vecCounts==0,
            method="BFGS", 
            control=list(maxit=lsMuModelGlobal$MAXIT, 
                         reltol=lsMuModelGlobal$RELTOL,
                         fnscale=-1)
        )[c("par","value","convergence")] )
    }, error=function(strErrorMsg){
        print(paste0("ERROR: Fitting impulse model: fitContinuousZINB()."))
        print(strErrorMsg)
        print(paste0("vecTheta ", paste0(vecTheta, collapse = " ")))
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
        for(b in seq(1, lsDispModelGlobal$scaNConfounders, by=1)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + lsDispModelGlobal$vecNBatches[b] - 1
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
        
        vecMuModel[3:5][vecMuModel[3:5] < 10^(-10)] <- 10^(-10)
        vecMuModel[3:5][vecMuModel[3:5] > 10^(10)] <- 10^(10)
    } else if(lsMuModelGlobal$strMuModel == "splines") {
        vecMuModel <- fitDispMu[seq(
            from = scaNParamUsed + 1,
            to = scaNParamUsed + lsMuModelGlobal$scaNSplines,
            by = 1)]
        scaNParamUsed <- scaNParamUsed + length(vecMuModel)
        
        vecMuModel[vecMuModel < -10*log(10)] <- -10*log(10)
        vecMuModel[vecMuModel > 10*log(10)] <- 10*log(10)
    } else {
        stop("fitImpulseOneInitZINB()")
    }
    
    if(!is.null(lsMuModelGlobal$lsvecidxBatchAssign)){
        lsvecBatchFactorsMu <- list()
        for(b in seq(1, lsMuModelGlobal$scaNConfounders, by=1)){
            vecBatchFactors <- c( 1, exp(fitDispMu[
                seq(from = scaNParamUsed+1, 
                    to = scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1, 
                    by = 1)]) )
            scaNParamUsed <- scaNParamUsed + lsMuModelGlobal$vecNBatches[b] - 1
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
#' Calls optimisation wrapper \code{fitImpulseOneInitZINB} 
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
fitImpulseZINB <- function(
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
        lsvecBatchModel=lsvecBatchParamGuess,
        lsMuModelGlobal=lsMuModelGlobal,
        vecInterval=NULL )
    vecDispParam <- decompressDispByGene(
        vecDispModel=vecDispGuess,
        lsvecBatchModel=lsvecBatchParamGuessDisp,
        lsDispModelGlobal=lsDispModelGlobal,
        vecInterval=NULL )
    if(is.null(vecPiParam)){
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
        lsMuModelGlobal=lsMuModelGlobal,
        vecMu=vecMuParam,
        vecDisp=vecDispParam,
        vecDrop=vecPiParamInit)
    vecParamGuessPeak <- lsParamGuesses$peak
    vecParamGuessValley <- lsParamGuesses$valley
    
    # (II) Compute new parameters
    # 1. Initialisation: Prior best fit
    vecParamGuessPrior <- vecImpulseParamGuess
    vecParamGuessPrior[3:5] <- log(vecParamGuessPrior[3:5])
    lsFitPrior <- fitContinuousZINB(
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
        lsFitPeak <- fitContinuousZINB(
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
        lsFitValley <- fitContinuousZINB(
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
                vecMu=vecMuParam*vecNormConst,
                vecDisp=vecDispParam, 
                vecPi=vecPiParam,
                vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
                vecboolZero= !is.na(vecCounts) &vecCounts==0)
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
        scaLLGuess <- evalLogLikGene(
            vecCounts=vecCounts,
            vecMu=vecMuParam*lsMuModelGlobal$vecNormConst,
            vecDisp=rep(scaDispGuess, length(vecCounts)), 
            vecPi=vecPiParam,
            vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
            vecboolZero= !is.na(vecCounts) &vecCounts==0)
        if(lsFitPrior$scaLL < scaLLGuess){
            lsFitBest <- list(scaDisp=scaDispGuess,
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
fitZINBMuDisp <- function(
    matCountsProc,
    vecNormConst,
    lsMuModel,
    lsDispModel,
    lsDropModel,
    matWeights){
    
    scaNumGenes <- dim(matCountsProc)[1]
    scaNumCells <- dim(matCountsProc)[2]
    
    lsFitDispMu <- bplapply(seq(1,scaNumGenes), function(i){
        if(!is.null(lsMuModel$lsMuModelGlobal$vecConfounders)) {
            lsvecBatchParamGuessMu <- lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,2:dim(mat)[2],drop=FALSE] )
        } else {
            lsvecBatchParamGuessMu <- NULL
        }
        if(!is.null(lsDispModel$lsDispModelGlobal$vecConfounders)) {
            lsvecBatchParamGuessDisp <- lapply(lsDispModel$lsmatBatchModel, function(mat) mat[i,2:dim(mat)[2],drop=FALSE] )
        } else {
            lsvecBatchParamGuessDisp <- NULL
        }
        
        # Decompress drop-out rates if not function of mean
        if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
            vecPiParam <- NULL
        } else {
            vecPiParam <- decompressDropoutRateByGene(
                matDropModel=lsDropModel$matDropoutLinModel,
                vecMu=NULL,
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal)
        } 
        
        if(lsMuModel$lsMuModelGlobal$strMuModel=="MM"){
            fitDispMu <- fitMMZINB(
                vecCounts=matCountsProc[as.double(i),], # sparseMatrix index must be double
                vecMuGuess=lsMuModel$matMuModel[i,], # Mu model
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,], # Disp model
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam,
                matWeights=matWeights,
                MAXIT=lsMuModel$lsMuModelGlobal$MAXIT_BFGS_MuDisp,
                RELTOL=lsMuModel$lsMuModelGlobal$RELTOL_BFGS_MuDisp )
            
            return(fitDispMu)
        } else if(lsMuModel$lsMuModelGlobal$strMuModel=="impulse"){      
            # Estimate mean parameters
            fitDispMu <- fitImpulseZINB(
                vecCounts=matCountsProc[as.double(i),], # sparseMatrix index must be double
                vecImpulseParamGuess=lsMuModel$matMuModel[i,], # Mu model
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,], # Disp model
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam )
            return(fitDispMu)
        } else if(lsMuModel$lsMuModelGlobal$strMuModel %in% c("splines", "groups", "constant")){      
            # Estimate mean parameters
            fitDispMu <- fitContinuousZINB(
                vecCounts=matCountsProc[as.double(i),], # sparseMatrix index must be double
                vecMuModelGuess=lsMuModel$matMuModel[i,], # Mu model
                lsvecBatchParamGuessMu=lsvecBatchParamGuessMu,
                lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                vecDispGuess=lsDispModel$matDispModel[i,], # Disp model
                lsvecBatchParamGuessDisp=lsvecBatchParamGuessDisp,
                lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                matDropoutLinModel=lsDropModel$matDropoutLinModel, # Pi model
                vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                lsDropModelGlobal=lsDropModel$lsDropModelGlobal,
                vecPiParam=vecPiParam )
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
    })
    
    matDispModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$vecDispModel))
    matMuModel <- do.call(rbind, lapply(lsFitDispMu,  function(i) i$vecMuModel))
    if(!is.null(lsMuModel$lsMuModelGlobal$scaNConfounders)) {
        lsmatBatchModelMu <- lapply(seq(1, lsMuModel$lsMuModelGlobal$scaNConfounders, by=1), function(confounder){
            do.call(rbind, lapply(lsFitDispMu,  function(i) i$lsvecBatchFactorsMu[[confounder]] ))
        })
    } else {
        lsmatBatchModelMu <- NULL
    }
    if(!is.null(lsDispModel$lsDispModelGlobal$scaNConfounders)) {
        lsmatBatchModelDisp <- lapply(seq(1, lsDispModel$lsDispModelGlobal$scaNConfounders, by=1), function(confounder){
            do.call(rbind, lapply(lsFitDispMu,  function(i) i$lsvecBatchFactorsDisp[[confounder]] ))
        })
    } else {
        lsmatBatchModelDisp <- NULL
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