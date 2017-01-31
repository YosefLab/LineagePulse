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
evalLogLikPiZINB <- function(vecTheta,
														 vecCounts,
														 matMu,
														 scaNormConst,
														 matDisp,
														 matPiPredictors,
														 vecboolNotZero,
														 vecboolZero ){ 
	
	# (I) Linker functions
	# Force mean parameter to be negative
	vecTheta[2] <- -exp(vecTheta[2])
	
	# (II) Prevent parameter shrinkage/explosion
	vecTheta[vecTheta < -10^(10)] <- -10^(10)
	vecTheta[vecTheta > 10^(10)] <- 10^(10)
	if(vecTheta[2] > -.Machine$double.eps){ vecTheta[2] <- -.Machine$double.eps }
	
	vecPiEst <- sapply(seq(1,dim(matPiPredictors)[1]), function(i){
		evalDropoutModel_comp(vecPiModel=vecTheta, 
													vecPiPredictors=matPiPredictors[i,])
	})
	
	# (III) Evaluate loglikelihood of estimate
	# Loglikelihood is evaluated on each window which was has
	# target cell in its neighbourhood. Note that the negative binomial
	# parameters therefore change for each iteration as these represent
	# SEPARATE windows, and not a smoothing of the target cell to its 
	# neighbours. ("Trans-terms")
	scaLogLik <- sum(sapply(seq(1,dim(matMu)[2]), function(cell){
		evalLogLikZINB_comp( vecCounts=vecCounts,
												 vecMu=matMu[,cell]*scaNormConst,
												 vecDisp=matDisp[,cell], 
												 vecPi=vecPiEst,
												 vecboolNotZero=vecboolNotZero, 
												 vecboolZero=vecboolZero )
	}))
	
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
evalLogLikPiZINB_comp <- cmpfun(evalLogLikPiZINB)

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

fitPiZINB <- function(vecCounts,
											lsMuModel,
											lsDispModel,
											lsDropModel,
											vecidxInterval,
											idxTarget){ 
	
	scaNumGenes <- length(vecCounts)
	# Decompress parameters
	matMuParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
		decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,],
													lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
													lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
													vecInterval=vecidxInterval)
	}))
	matDispParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
		decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
												 lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
												 vecInterval=vecidxInterval)
	}))
	matPiPredictors <- cbind(1, log(matMuParam[,match(idxTarget,vecidxInterval)]), 
													 matPiConstPredictors)
	
	# (I) Numerical maximum likelihood estimation of linear model
	# Initialise optimisation of linear model:
	#vecParamGuess <- rep(1, dim(matPiPredictors)[2])
	vecParamGuess <- vecDropoutLinModel
	vecParamGuess[2] <- log(-vecParamGuess[2])
	boolError <- FALSE
	lsLinModelFit <- tryCatch({
		optim(
			par=vecParamGuess,
			evalLogLikPiZINB_comp,
			matPiPredictors=lsDropModel$lsDropModelGlobal$matPiPredictors,
			vecCounts=vecCounts,
			matMu=matMuParam,
			matDisp=matDispParam,
			scaNormConst=lsMuModel$lsMuModelGlobal$vecNormConst[idxTarget],
			vecboolNotZero=!is.na(vecCounts) & vecCounts>0,
			vecboolZero= !is.na(vecCounts) & vecCounts==0,
			method="BFGS",
			control=list(maxit=lsPiOptimHyperparam$MAXIT_BFGS_Pi,
									 reltol=lsPiOptimHyperparam$RELTOL_BFGS_Pi,
									 fnscale=-1)
		)[c("par","convergence")]
	}, error=function(strErrorMsg){
		print(paste0("ERROR: Fitting logistic drop-out model: fitPiZINB().",
								 " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
		print(strErrorMsg)
		scaLLInit <- evalLogLikPiZINB_comp(
			vecTheta=vecParamGuess,
			matPiPredictors=lsDropModel$lsDropModelGlobal$matPiPredictors,
			vecCounts=vecCounts,
			matMu=matMuParam,
			matDisp=matDispParam,
			scaNormConst=lsMuModel$lsMuModelGlobal$vecNormConst[idxTarget],
			vecboolNotZero= !is.na(vecCounts) & vecCounts>0,
			vecboolZero= !is.na(vecCounts) & vecCounts==0)
		print(paste0("scaLLInit", scaLLInit))
		print(paste0("vecParamGuess ", paste(vecParamGuess, collapse=" ")))
		lsErrorCausingGene <- list( vecParamGuess=vecParamGuess,
																vecCounts=vecCounts, 
																matPiPredictors=lsDropModel$lsDropModelGlobal$matPiPredictors,
																matMuParam=matMuParam, 
																matDispParam=matDispParam,
																scaNormConst=lsMuModel$lsMuModelGlobal$vecNormConst[idxTarget],
																scaLLInit=scaLLInit )
		save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
		boolError <- TRUE
		# Return intialisation
		return(vecParamGuess)
	})
	
	# (II) Extract results and correct for sensitivity boundaries
	vecLinModel <- unlist(lsLinModelFit["par"])
	vecLinModel[2] <- -exp(vecLinModel[2])
	# # Catch boundary of likelihood domain on parameter space
	vecLinModel[vecLinModel < -10^(10)] <- -10^(10)
	vecLinModel[vecLinModel > 10^(10)] <- 10^(10)
	if(vecLinModel[2] > -.Machine$double.eps){ vecLinModel[2] <- -.Machine$double.eps }
	
	if(boolError){
		scaConvergence <- 1001
	} else {
		scaConvergence <- unlist(lsLinModelFit["convergence"])    
	}
	
	return( list(vecLinModel=vecLinModel,
							 scaConvergence=scaConvergence) )
}