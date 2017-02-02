#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++    Likelihood ZINB    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute likelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). This
#' is the core function used for every likelihood evaluation 
#' in LineagePulse, such as maximum likelihood-based estimation. 
#' It operates on a vector of counts, such as observations of a gene. Note that
#' for the sake of numerical stability, lower bounds on loglikelihood terms
#' are implemented.
#' 
#' @seealso Called directly by likelihood wrappers 
#' \code{evalLogLikSmoothZINB} and \code{evalLogLikGene}
#' Called directly by \code{evalLogLikPiZINB},
#' \code{evalLogLikDispConstMuClustersZINB}
#' and \code{evalLogLikMuWindowZINB}
#' Moreover called by \code{plotGene}.
#' Compiled verison: \link{evalLikZINB_comp}
#' Loglikelihood: \link{evalLogLikZINB}
#'
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (scalar vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#'    
#' @return vecLik: (vector number of samples) 
#'    Likelihood under zero-inflated
#' 	  negative binomial model of each sample.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLikZINB <- function(vecCounts,
                        vecMu,
                        vecDisp, 
                        vecPi, 
                        vecboolNotZero, 
                        vecboolZero){  
  
  # Note on handling very low probabilities: 
  # The following scalar is used as a sensitivity bound for
  # log likelihoods:
  # scaPrecLim <- -323*log(10)
  # Accordingle, this is the bound for likelihoods
  scaPrecLim <- 10^(-323)
  
  # Likelihood of zero counts:
  vecLikZeros <- (1-vecPi[vecboolZero])*
    (vecDisp[vecboolZero]/(vecDisp[vecboolZero] + 
                             vecMu[vecboolZero]))^vecDisp[vecboolZero] +
    vecPi[vecboolZero]
  # Likelihood of non-zero counts:
  vecLikNonzeros <- (1-vecPi[vecboolNotZero])*dnbinom(
    vecCounts[vecboolNotZero], 
    mu=vecMu[vecboolNotZero], 
    size=vecDisp[vecboolNotZero], 
    log=FALSE)
  # Compute likelihood of all data:
  vecLik <- array(NA, length(vecCounts))
  vecLik[vecboolZero] <- vecLikZeros 
  vecLik[vecboolNotZero] <- vecLikNonzeros
  # Catch low likehood observations
  vecLik[vecLik < scaPrecLim] <- scaPrecLim #?
  
  # Return vector of likelihoods for weighting.
  return(vecLik)
}

#' Compiled function: evalLikZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLikZINB}.
#' 
#' @seealso \link{evalLikZINB}
#' 
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (scalar vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#'    
#' @return scaLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLikZINB_comp <- cmpfun(evalLikZINB)

#' Compute loglikelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). This
#' is the core function used for every likelihood evaluation 
#' in LineagePulse, such as maximum likelihood-based estimation. 
#' It operates on a vector of counts, such as observations of a gene. Note that
#' for the sake of numerical stability, lower bounds on loglikelihood terms
#' are implemented.
#' 
#' @seealso Called directly by likelihood wrappers 
#' \code{evalLogLikSmoothZINB} and \code{evalLogLikGene}
#' Called directly by \code{evalLogLikPiZINB},
#' \code{evalLogLikDispConstMuClustersZINB}
#' and \code{evalLogLikMuWindowZINB}
#' Moreover called by \code{plotGene}.
#' Compiled verison: \link{evalLogLikZINB_comp}
#'
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (scalar vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikZINB <- function(vecCounts,
                           vecMu,
                           vecDisp, 
                           vecPi, 
                           vecboolNotZero, 
                           vecboolZero){  
  
  # Note on handling very low probabilities: vecLikZeros
  # typically does not have zero elements as it has the 
  # the summand drop-out rate. Also the log cannot be
  # broken up over the sum to dnbinom. In contrast to that,
  # the log is taken in dnbinom for vecLikNonzeros to avoid 
  # zero probabilities. Zero probabilities are handled
  # through substitution of the minimal probability under
  # machine precision.
  scaLogPrecLim <- -323*log(10)
  
  # Likelihood of zero counts:
  vecLogLikZeros <- log((1-vecPi[vecboolZero])*
                          (vecDisp[vecboolZero]/(vecDisp[vecboolZero] + 
                                                   vecMu[vecboolZero]))^vecDisp[vecboolZero] +
                          vecPi[vecboolZero])
  vecLogLikZeros[vecLogLikZeros < scaLogPrecLim] <- scaLogPrecLim
  scaLogLikZeros <- sum(vecLogLikZeros)
  # Likelihood of non-zero counts:
  vecLogLikNonzerosDropout <- log(1-vecPi[vecboolNotZero])
  vecLogLikNonzerosNB <- dnbinom(
    vecCounts[vecboolNotZero], 
    mu=vecMu[vecboolNotZero], 
    size=vecDisp[vecboolNotZero], 
    log=TRUE)
  vecLogLikNonzerosDropout[vecLogLikNonzerosDropout < scaLogPrecLim] <- scaLogPrecLim
  vecLogLikNonzerosNB[vecLogLikNonzerosNB < scaLogPrecLim] <- scaLogPrecLim
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( vecLogLikNonzerosDropout+vecLogLikNonzerosNB )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compiled function: evalLogLikZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikZINB}
#' 
#' @seealso \link{evalLogLikZINB}
#'
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (scalar vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikZINB_comp <- cmpfun(evalLogLikZINB)


#' Compute smoothed loglikelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' This likelihood function is a wrapper for computing the likelihood with
#' neighbourhood smoothing for a vector of counts.
#' 
#' @seealso Called directly by \code{evalLogLikGene} and by
#' loglikelihood wrappers that operate on sliding window mean
#' estimates: \code{evalLogLikMuVecWindowsZINB} and
#' \code{evalLogLikDispConstMuVecWindowsZINB}. Calls
#' \code{evalLogLikZINB} on neighbourhoods.
#' Compiled version: \link{evalLogLikSmoothZINB}
#'
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (numeric vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (numeric vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikSmoothZINB <- function(vecCounts,
                                 vecMu,
                                 vecNormConst,
                                 vecDisp, 
                                 vecPi, 
                                 vecboolNotZero, 
                                 vecboolZero,
                                 scaWindowRadius=NULL ){
  
  scaNumCells <- length(vecCounts)
  scaLogLik <- sum(sapply(seq(1,scaNumCells), 
                          function(j){
                            scaindIntervalStart <- max(1,j-scaWindowRadius)
                            scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
                            vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
                            scaLogLikCell <- evalLogLikZINB_comp(vecCounts=vecCounts[vecInterval],
                                                                 vecMu=vecMu[j]*vecNormConst[vecInterval],
                                                                 vecDisp=rep(vecDisp[j], length(vecInterval)), 
                                                                 vecPi=vecPi[vecInterval], 
                                                                 vecboolNotZero=vecboolNotZero[vecInterval], 
                                                                 vecboolZero=vecboolZero[vecInterval])
                            return(scaLogLikCell)
                          }))
  return(scaLogLik)
}

#' Compiled function: evalLogLikSmoothZINB
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikSmoothZINB}
#' 
#' @seealso \link{evalLogLikSmoothZINB}
#'
#' @param vecCounts (count vector number of samples)
#'    Observed read counts, not observed are NA.
#' @param vecMu (numeric vector number of samples) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (numeric vector number of samples) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of samples) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of samples)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikSmoothZINB_comp <- cmpfun(evalLogLikSmoothZINB)


#' Wrapper for log likelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' This likelihood function is a wrapper which choses whether to 
#' compute smoothed or non-smoothed likelihood. Calling the 
#' non-smoothed routing avoids the use of a for loop in the case
#' of no smoothing.
#' 
#' 
#' @seealso Called directly by \code{evalLogLikMatrix} and by
#' loglikelihood wrappers for various models:
#' \code{evalLogLikMuConstZINB},
#' \code{evalLogLikDispConstZINB},
#' \code{evalLogLikDispConstMuConstZINB},
#' \code{evalLogLikMuImpulseZINB} and
#' \code{evalLogLikDispConstMuImpulseZINB}.
#' Called by fitting wrappers to ensure convergence:
#' \code{fitMuImpulseZINB} and
#' \code{fitDispConstMuImpulseZINB}.
#' Additionally by \code[plotGene}. Calls either 
#' \code{evalLogLikZINB} or \code{evalLogLikSmoothZINB}.
#'
#' @param vecCounts (count vector number of cells)
#'    Observed read counts, not observed are NA.
#' @param vecMu (numeric vector number of cells) 
#'    Negative binomial mean parameter.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#' @param vecDisp: (numeric vector number of cells) 
#'    Negative binomial dispersion parameters.
#' @param vecPi: (probability vector number of cells) 
#'    Drop-out rate estimates.
#' @param vecboolNotZero: (bool vector number of cells)
#'    Whether observation is larger than zero.
#' @param vecboolZero: (bool vector number of cells)
#'    Whether observation is zero.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikGene <- function(vecCounts,
                           vecMu,
                           vecNormConst,
                           vecDisp, 
                           vecPi,
                           vecboolNotZero, 
                           vecboolZero,
                           scaWindowRadius=NULL ){
  
  if(!is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikSmoothZINB_comp(vecCounts=vecCounts,
                                           vecMu=vecMu,
                                           vecNormConst=vecNormConst,
                                           vecDisp=vecDisp, 
                                           vecPi=vecPi,
                                           vecboolNotZero=vecboolNotZero, 
                                           vecboolZero=vecboolZero,
                                           scaWindowRadius=scaWindowRadius)
  } else {
    scaLogLik <- evalLogLikZINB_comp(vecCounts=vecCounts,
                                     vecMu=vecMu*vecNormConst,
                                     vecDisp=vecDisp, 
                                     vecPi=vecPi,
                                     vecboolNotZero=vecboolNotZero, 
                                     vecboolZero=vecboolZero )
  }
  return(scaLogLik)
}

evalLogLikGeneMM <- function(vecCounts,
                             matMuParam, # every column is rep of one scalar (one model, evaluated at each cell!)
														 vecNormConst,
														 matDispParam,
														 matDropParam,
														 matWeights,
                             vecboolNotZero, 
                             vecboolZero ){
  
  # Loop over models:
  matLikSum <- do.call(cbind, lapply(seq(1, dim(matWeights)[2]), function(m){
  	evalLikZINB_comp(vecCounts=vecCounts,
                     vecMu=matMuParam[,m]*vecNormConst,
                     vecDisp=matDispParam[,m], 
                     vecPi=matDropParam[,m],
                     vecboolNotZero=vecboolNotZero, 
                     vecboolZero=vecboolZero ) *matWeights[,m]
  }))
  scaLogLik <- sum(log(apply(matLikSum,1,sum)), na.rm=TRUE)
  
  return(scaLogLik)
}

evalLogLikGeneMM_comp <- cmpfun(evalLogLikGeneMM)

evalLogLikCellMM <- function(vecCounts,
                             matMuParam,
                             matDispParam,
                             matDropParam,
                             scaNormConst,
                             vecWeights,
                             vecboolNotZero, 
                             vecboolZero ){
  
  scaNGenes <- length(vecCounts)
  # Initialise sum of likelihoods, this is the mixture model sum under the log
  vecLikSum <- array(0, scaNGenes)
  # Loop over models:
  matLikSum <- do.call(cbind, lapply(seq(1, length(vecWeights)), function(m){
    # Evaluate loglikelihood of observations of cell under current model
    evalLikZINB_comp(vecCounts=vecCounts,
                     vecMu=matMuParam[,m]*scaNormConst,
                     vecDisp=matDispParam[,m], 
                     vecPi=matDropParam[,m],
                     vecboolNotZero=vecboolNotZero, 
                     vecboolZero=vecboolZero )*vecWeights[m]
  }))
  scaLogLik <- sum(log(apply(matLikSum,1,sum)), na.rm=TRUE)
  
  return(scaLogLik)
}

evalLogLikCellMM_comp <- cmpfun(evalLogLikCellMM)

#' Wrapper for log likelihood of zero-inflated negative binomial model
#' for a matrix of counts (parallelised).
#' 
#' This likelihood function is a wrapper computes loglikelihood
#' of entire data set by parallelising loglikelihood computation
#' over genes.
#' 
#' @seealso Called directly by \code{fitZINB} to track
#' convergence of estimation iteration on entire data set.
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
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval radius.
#'
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
evalLogLikMatrix <- function(matCounts,
                             lsMuModel,
                             lsDispModel, 
                             lsDropModel,
                             matWeights=NULL,
                             scaWindowRadius=NULL ){
  
  scaNGenes <- dim(matCounts)[1]
  scaNCells <- dim(matCounts)[2]
  
  scaLogLik <- sum(unlist(
    bplapply( seq(1,scaNGenes), function(i){
      if(lsMuModel$lsMuModelGlobal$strMuModel=="MM"){
        
        matMuParam <- do.call(cbind, lapply(seq(1,dim(matWeights)[2]), function(m){
      		decompressMeansByGene(vecMuModel=lsMuModel$matMuModel[i,m],
      													lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
      													lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
      													vecInterval=NULL)
      	}))
      	
        if(lsDispModel$lsDispModelGlobal$strDispModel=="constant"){
          vecDispParam <- decompressDispByGene(vecDispModel=lsDispModel$matDispModel[i,],
                                               lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                                               vecInterval=NULL)
          matDispParam <- do.call(cbind, lapply(seq(1,dim(matWeights)[2]), function(m) vecDispParam ))
        } else {
          stop(paste0("ERROR evalLogLikMatrix(): strDispModel=", strDispModel, " not recognised."))
        }
      	
      	matDropParam <- do.call(cbind, lapply(seq(1,dim(matWeights)[2]), function(m){
      		decompressDropoutRateByGene(matDropModel=lsDropModel$matDropoutLinModel,
      																							vecMu=matMuParam[,m],
      																							vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
      	}))
      	
        vecCounts <- matCounts[i,]
        scaLL <- evalLogLikGeneMM(vecCounts=vecCounts,
                                  matMuParam=matMuParam,
        													vecNormConst=lsMuModel$lsMuModelGlobal$vecNormConst,
        													matDispParam=matDispParam,
        													matDropParam=matDropParam,
        													matWeights=matWeights,
                                  vecboolNotZero= !is.na(vecCounts) & vecCounts>=0, 
                                  vecboolZero= !is.na(vecCounts) & vecCounts==0 )
        
      } else {
        # Decompress parameters by gene
        vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
        																		 lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, function(mat) mat[i,] ),
                                             lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                                             vecInterval=NULL )
        vecDispParam <- decompressDispByGene( vecDispModel=lsDispModel$matDispModel[i,],
                                              lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                                              vecInterval=NULL )
        vecPiParam <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
                                                   vecMu=vecMuParam,
                                                   vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
        
        # Evaluate loglikelihood of gene
        vecCounts <- matCounts[i,]
        scaLL <- evalLogLikGene(vecCounts=vecCounts,
                                vecMu=vecMuParam,
                                vecNormConst=lsMuModel$lsMuModelGlobal$vecNormConst,
                                vecDisp=vecDispParam, 
                                vecPi=vecPiParam,
                                vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
                                vecboolZero= !is.na(vecCounts) & vecCounts==0,
                                scaWindowRadius=scaWindowRadius)
      }
      return(scaLL)
    })
  ))
  return(scaLogLik)
}