#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++    Model-free DE analysis  ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Differential expression analysis
#' 
#' Performs differential expression analysis. This functiion is
#' broken up into:
#' (I) Fit null model.
#' (II) Compute likelihood of full and reduced model.
#' (II) Differential expression analysis.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param vecDispersions: (numeric matrix genes x clusters)
#'    Inferred negative binomial dispersions.
#' @param matMuCluster: (numeric matrix genes x clusters)
#'    Inferred negative binomial cluster means.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param boolConvergenceZINBH1: (bool) Convergence of EM algorithm for
#'    full model.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes) 
#'    Summary of model-free differential expression analysis.
#' @export

runDEAnalysis <- function(matCountsProc,
  vecPseudotime,
  vecSizeFactors,
  matDispersionsH1,
  matMuH1,
  matDropoutH1,
  scaKbyGeneH1,
  matDispersionsH0,
  matMuH0,
  matDropoutH0,
  scaKbyGeneH0,
  scaWindowRadius=NULL ){
  
  # Compute log likelihoods
  matboolNotZeroObserved <- matCountsProc >0 & !is.na(matCountsProc)
  matboolZero <- matCountsProc==0
  
  vecLogLikFull <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
      vecMu=matMuH1[i,]*vecSizeFactors,
      vecDispEst=matDispersionsH1[i,], 
      vecDropoutRateEst=matDropoutH1[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  vecLogLikRed <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
      vecMu=matMuH0[i,]*vecSizeFactors,
      vecDispEst=matDispersionsH0[i,], 
      vecDropoutRateEst=matDropoutH0[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  
  # (II) Differential expression analysis
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDeltaDegFreedom <- scaKbyGeneH1 - scaKbyGeneH0
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH = p.adjust(vecPvalue, method = "BH")
  
  dfModelFreeDEAnalysis =   as.data.frame(cbind(
    "Gene" = row.names(matCountsProc),
    "p"=as.numeric(vecPvalue),
    "adj.p"=as.numeric(vecPvalueBH),
    "loglik_full"=vecLogLikFull,
    "loglik_red"=vecLogLikRed,
    "deviance"=vecDeviance,
    "mean"=vecMuClusterH0,
    "dispersion_H0"=vecDispersionsH0,
    "converged_H0"=rep(boolConvergenceZINBH0!=0,dim(matCountsProc)[1]),
    "converged_full"=rep(all(boolConvergenceZINBH1),dim(matCountsProc)[1]),
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}