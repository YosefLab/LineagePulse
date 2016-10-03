#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++   Differential expression analysis  +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Differential expression analysis
#' 
#' Performs differential expression analysis based on previously estimated 
#' null and alternative models.
#' (I) Compute loglikelihood of data under null H0 and alternative H1 model.
#' (II) Differential expression analysis as loglikelihood ratio test.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param matMuH1: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial mean parameters
#'    under alternative model H1.
#' @param matDispersionsH1: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial dispersion parameters
#'    under alternative model H1.
#' @param matDropoutH1: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#'    These are the observation-wise point estimates, not the
#'    logistic functions. Fit under alternative model H1.
#' @param scaKbyGeneH1: (scalar) Degrees of freedom
#'    by gene used in alternative model H1. The logistic
#'    dropout model is ignored as it is shared between
#'    alternative and null model.
#' @param matMuH0: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial mean parameters
#'    under null model H0.
#' @param matDispersionsH0: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial dispersion parameters
#'    under null model H0.
#' @param matDropoutH0: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#'    These are the observation-wise point estimates, not the
#'    logistic functions. Fit under null model H0.
#' @param scaKbyGeneH0: (scalar) Degrees of freedom
#'    by gene used in null model H0. The logistic
#'    dropout model is ignored as it is shared between
#'    alternative and null model.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes x reported variables) 
#'    Summary of differential expression analysis, sorted by adj.p:
#'    {Gene: gene ID,
#'    p: raw p-value, 
#'    adj.p: BH corrected p-value, 
#'    loglik_full: loglikelihood of alternative model H1,
#'    loglik_red: loglikelihood of null model H0,
#'    deviance: loglikelihood ratio test statistic (the deviance),
#'    mean_H0: inferred gene-wise mean parameter (constant null model),
#'    dispersion_H0: inferred gene-wise dispersion parameter (constant null model)}
#' @export

runDEAnalysis <- function(matCountsProc,
  vecSizeFactors,
  lsMuModelH1,
  lsDispModelH1,
  lsMuModelH0,
  lsDispModelH0,
  lsDropModel,
  scaKbyGeneH1,
  scaKbyGeneH0,
  scaWindowRadius=NULL ){
  
  # (I) Compute log likelihoods
  lsLL <- bplapply( seq(1,dim(matCountsProc)[1]), function(i){
    vecCounts <- matCountsProc[i,]
    vecboolNotZeroObserved <- vecCounts >0 & !is.na(vecCounts)
    vecboolZero <- vecCounts==0
    
    # Decompress parameters by gene
    vecMuParamH0 <- decompressMeans( vecMuModel=lsMuModelH0$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH0$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH0 <- decompressDispersions( vecDispModel=lsDispModelH0$matDispModel[i,],
      lsDispModel=lsDispersionModelH0$lsDispModelGlobal,
      vecInterval=NULL )
    vecDropoutParamH0 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH0,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
    
    # Decompress parameters by gene
    vecMuParamH1 <- decompressMeans( vecMuModel=lsMuModelH1$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH1 <- decompressDispersions( vecDispModel=lsDispModelH1$matDispModel[i,],
      lsDispModel=lsDispersionModelH1$lsDispModelGlobal,
      vecInterval=NULL )
    vecDropoutParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH1,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
    
    # Compute loglikelihoods
    scaLogLikFull <- evalLogLikGene(vecCounts=vecCounts,
      vecSizeFactors=vecSizeFactors,
      vecMu=vecMuParamH1,
      vecDispEst=vecDispParamH1, 
      vecDropoutRateEst=vecDropoutParamH1,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
    
    scaLogLikRed <- evalLogLikGene(vecCounts=vecCounts,
      vecMu=vecMuParamH0,
      vecSizeFactors=vecSizeFactors,
      vecDispEst=vecDispParamH0, 
      vecDropoutRateEst=vecDropoutParamH0,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
    
    return(list( scaLogLikFull=scaLogLikFull,
      scaLogLikRed=scaLogLikRed ))
  })
  vecLogLikFull <- sapply(lsLL, function(LL_i) LL_i[["scaLogLikFull"]])
  vecLogLikRed <- sapply(lsLL, function(LL_i) LL_i[["scaLogLikRed"]])
  
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
    "mean_H0"=lsMuModelH0$matMuModel,
    "dispersion_H0"=lsDispModelH0$matDispModel,
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}