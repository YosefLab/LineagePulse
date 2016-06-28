#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++    Likelihood ZINB    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute log likelihood of negative binomial model for a vector of 
#' zero observations
#' 
#' This is a custom version of dnbinom only valid for x=0 which is hopefully
#' faster. The likihood is computed for zero observations, with the number
#' of observations assumed to be the number of reported means and dispersions.
#' Guards against numerical errors.
#' 
#' @aliases evalLogLikNBZero_comp
#' 
#' @seealso Called by \code{fitZINB} and
#' \code{runModelFreeDEAnalysis}.
#'
#' @param vecMu (vector number of observations) Negative binomial
#'    mean parameter for each observation in gene.
#' @param vecDispEst: (vector number of observations) Negative binomial  
#'    dispersion parameter for each observation in gene.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikNBZero <- function(vecMu,
  vecDispEst){
  
  vecLogLik <- vecDispEst*(log(vecDispEst) - log(vecDispEst + vecMu))
  scaLogLik <- sum( log(vecLogLik[vecLogLik!=0]) +
      sum(vecLogLik==0)*log(.Machine$double.eps) )
  
  return(scaLogLik)
}

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
#' @param vecY (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param scaDispEst: (scalar) Negative binomial dispersion 
#'    parameter for given gene.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikZINB_PseudoDE <- function(vecY,
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
  #vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
  #  dnbinom(
  #    vecY[vecboolZero], 
  #    mu=vecMu[vecboolZero], 
  #    size=vecDispEst[vecboolZero], 
  #    log=FALSE) +
  #  vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  #scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
  #    sum(vecLikZeros==0)*log(.Machine$double.eps) )
  scaLogLikZeros <- evalLogLikNBZero_comp(vecMu=vecMu[vecboolZero],
    vecDispEst=vecDispEst[vecboolZero])
  # Likelihood of non-zero counts:
  vecLikNonzeros <- log(1-vecDropoutRateEst[vecboolNotZeroObserved]) +
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMu[vecboolNotZeroObserved], 
      size=vecDispEst[vecboolNotZeroObserved], 
      log=TRUE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( vecLikNonzeros[is.finite(vecLikNonzeros)]) +
      sum(!is.finite(vecLikNonzeros))*log(.Machine$double.eps)
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}