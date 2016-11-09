#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++    Likelihood ZINB    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute log likelihood of zero-inflated negative binomial model for one gene
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). This
#' is the core function used for every maximum likelihood-based estimation
#' in LineagePulse.
#' 
#' @aliases evalLogLikZINB_comp
#' 
#' @seealso Called directly by likelihood wrappers 
#'    \code{evalLogLikSmoothZINB} and \code{evalLogLikGene}.
#'    Called directly by \code{evalLogLikDispConstMuClustersZINB}
#'    and 
#'    Moreover by \code{plotGene}.
#'
#' @param vecCounts (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param vecDisp: (scalar vector number of samples) Negative binomial dispersion 
#'    parameter for given gene and observations.
#' @param vecPi: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#'    
#' @return scaLogLik: (scalar) Likelihood under zero-inflated
#' 	  negative binomial model.
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

#' Compute smoothed log likelihood of zero-inflated negative binomial model for one gene
#' 
#' This liklihood function is a wrapper for computing the likelihood with
#' neighbourhood smoothing.
#' 
#' @aliases evalLogLikSmoothZINB_comp
#' 
#' @seealso Called directly by \code{evalLogLikGene} and by
#' loglikelihood wrappers that operate on sliding window mean
#' estimates: \code{evalLogLikMuVecWindowsZINB} and
#' \code{evalLogLikDispConstMuVecWindowsZINB}.
#'
#' @param vecCounts (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDisp: (scalar vector number of samples) Negative binomial dispersion 
#'    parameter for given gene and observations.
#' @param vecPi: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZero: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
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

evalLogLikMatrix <- function(matCounts,
  vecNormConst,
  lsMuModel,
  lsDispModel, 
  lsDropModel,
  scaWindowRadius=NULL ){
  
  scaLogLik <- sum(unlist(
    bplapply( seq(1,dim(matCounts)[1]), function(i){
      # Decompress parameters by gene
      vecMuParam <- decompressMeansByGene( vecMuModel=lsMuModel$matMuModel[i,],
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
        vecNormConst=vecNormConst,
        vecDisp=vecDispParam, 
        vecPi=vecPiParam,
        vecboolNotZero= !is.na(vecCounts) & vecCounts>0, 
        vecboolZero= !is.na(vecCounts) & vecCounts==0,
        scaWindowRadius=scaWindowRadius)
      return(scaLL)
    })
  ))
  return(scaLogLik)
}