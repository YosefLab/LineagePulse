#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++    Likelihood ZINB    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute log likelihood of zero-inflated negative binomial model for one gene
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq).
#' 
#' @aliases evalLogLikZINB_comp
#' 
#' @seealso Called by \code{fitZINB} and
#' \code{runModelFreeDEAnalysis}.
#'
#' @param vecCounts (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param vecDispEst: (scalar vector number of samples) Negative binomial dispersion 
#'    parameter for given gene and observations.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikZINB_LinPulse <- function(vecCounts,
  vecMu,
  vecDispEst, 
  vecDropoutRateEst, 
  vecboolNotZeroObserved, 
  vecboolZero){  

  # Note on handling very low probabilities: vecLikZeros
  # typically does not have zero elements as it has the 
  # the summand drop-out rate. Also the log cannot be
  # broken up over the sum to dnbinom. In contrast to that,
  # the log is taken in dnbinom for vecLikNonzeros to avoid 
  # zero probabilities. Zero probabilities are handled
  # through substitution of the minimal probability under
  # machine precision.
  
  # Likelihood of zero counts:
  vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
      (vecDispEst[vecboolZero]/(vecDispEst[vecboolZero] + 
          vecMu[vecboolZero]))^vecDispEst[vecboolZero] +
      vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros > .Machine$double.eps]) +
      sum(vecLikZeros <= .Machine$double.eps)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLogLikNonzeros <- log(1-vecDropoutRateEst[vecboolNotZeroObserved]) +
    dnbinom(
      vecCounts[vecboolNotZeroObserved], 
      mu=vecMu[vecboolNotZeroObserved], 
      size=vecDispEst[vecboolNotZeroObserved], 
      log=TRUE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( vecLogLikNonzeros[vecLogLikNonzeros > log(.Machine$double.eps)] ) +
      sum(vecLogLikNonzeros <= log(.Machine$double.eps))*log(.Machine$double.eps)
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compute smoothed log likelihood of zero-inflated negative binomial model for one gene
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). It includes
#' a smoothin penalty on the means.
#' 
#' @aliases evalLogLikSmoothZINB_comp
#' 
#' @seealso Called by \code{fitZINB} and
#' \code{runModelFreeDEAnalysis}.
#'
#' @param vecCounts (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param vecDispEst: (scalar vector number of samples) Negative binomial dispersion 
#'    parameter for given gene and observations.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikSmoothZINB_LinPulse <- function(vecCounts,
  vecMu,
  vecDispEst, 
  vecDropoutRateEst, 
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL ){
  
  scaNumCells <- length(vecCounts)
  scaLogLik <- sum(lapply(seq(1,scaNumCells), 
    function(j){
      scaindIntervalStart <- max(1,j-scaWindowRadius)
      scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
      vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
      scaLogLikCell <- evalLogLikZINB_LinPulse_comp(vecCounts=vecCounts[vecInterval],
        vecMu=rep(vecMu[j], length(vecInterval)),
        vecDispEst=rep(vecDispEst[j], length(vecInterval)), 
        vecDropoutRateEst=vecDropoutRateEst[vecInterval], 
        vecboolNotZeroObserved[vecInterval], 
        vecboolZero[vecInterval])
      return(scaLogLikCell)
    }))
  return(scaLogLik)
}