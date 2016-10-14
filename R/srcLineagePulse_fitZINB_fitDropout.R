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
#' @seealso Called by \code{fitPiZINB_LinPulse}.
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
  scaNormConst,
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
  # neighbours. ("Trans-terms")
  scaLogLik <- sum(sapply(seq(1,dim(matMu)[2]), function(cell){
    evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=matMu[,cell]*scaNormConst,
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
#' @param vecDropoutLinModel: (numeric vector length linear model)
#'    Previous parameterisation of linear model for drop-out 
#'    rate in logit space for given cell.
#' @param matPiConstPredictors: (matrix genes x predictors) Predictors of
#'    the drop-out rate in the linear model. Minimum are a constant
#'    offset and log of the negative binomial mean parameter. 
#'    Other gene-specific predictors can be added.
#' @param vecCounts: (numeric vector number of genes) Observed expression values 
#'    of gene in target cell.
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
#' @param scaNormConst: (numeric vector number of cells in neighbourhood) 
#'    Model scaling factors for cell which takes
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecInterval: (integer vector neighbourhood)
#'    Positions of cells within smooting interval (neighbourhood)
#'    of target cell.
#' @param scaTarget: (integer) Index of target cell in interval.
#' 
#' @return vecLinModel: (numeric vector length linear model) 
#'    Linear model for drop-out rate in logit space for given cell.
#' @export

fitPiZINB_LinPulse <- function( vecDropoutLinModel,
  matPiConstPredictors,
  vecCounts,
  lsMuModel,
  lsDispModel,
  scaNormConst,
  vecInterval,
  scaTarget){ 
  
  scaNumGenes <- length(vecCounts)
  # Decompress parameters
  matMuParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
    decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,],
      lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
      vecInterval=vecInterval)
  }))
  matDispParam <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
    decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
      lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
      vecInterval=vecInterval)
  }))
  matPiPredictors <- cbind(1, log(matMuParam[,scaTarget]), 
    matPiConstPredictors)
  
  # Numerical maximum likelihood estimaton of linear model
  # Initialise optimisation of linear model:
  #vecParamGuess <- rep(1, dim(matPredictorsPi)[2])
  vecParamGuess <- vecDropoutLinModel
  vecParamGuess[2] <- log(-vecParamGuess[2])
  lsLinModelFit <- tryCatch({
    optim(
      par=vecParamGuess,
      evalLogLikPiZINB_LinPulse_comp,
      matPredictorsPi=matPiPredictors,
      vecCounts=vecCounts,
      matMu=matMuParam,
      matDisp=matDispParam,
      scaNormConst=scaNormConst,
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
    print(paste0("scaNormConst ", paste(scaNormConst,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, matDisp, matMu, scaNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "matDisp", "matMu","scaNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  vecLinModel <- unlist(lsLinModelFit["par"])
  vecLinModel[2] <- -exp(vecLinModel[2])
  scaConvergence <- unlist(lsLinModelFit["convergence"]) # Not used atm
  
  return(vecLinModel)
}