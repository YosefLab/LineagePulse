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
#' @param: objectLineagePulse: (LineagePulseObject)
#' LineagePulseObject with fitted null and alternative models.
#' 
#' @return objectLineagePulse: (LineagePulseObject)
#' LineagePulseObject with analysis summary table (dfResults)
#' added.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
runDEAnalysis <- function(objectLineagePulse){
  
  # (I) Compute log likelihoods
  vecLogLikFull <- evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                lsMuModel=objectLineagePulse@lsMuModelH1,
                                lsDispModel=objectLineagePulse@lsDispModelH1, 
                                lsDropModel=objectLineagePulse@lsDropModel,
                                matWeights=objectLineagePulse@matWeights)
  vecLogLikRed <- evalLogLikMatrix(matCounts=objectLineagePulse@matCountsProc,
                                   lsMuModel=objectLineagePulse@lsMuModelH0,
                                   lsDispModel=objectLineagePulse@lsDispModelH0, 
                                   lsDropModel=objectLineagePulse@lsDropModel,
                                   matWeights=objectLineagePulse@matWeights)

  # (II) Differential expression analysis
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDFbyGeneH1 <- objectLineagePulse@lsMuModelH1$lsMuModelGlobal$scaDegFreedom + 
    objectLineagePulse@lsDispModelH1$lsDispModelGlobal$scaDegFreedom
  scaDFbyGeneH0 <- objectLineagePulse@lsMuModelH0$lsMuModelGlobal$scaDegFreedom + 
    objectLineagePulse@lsDispModelH0$lsDispModelGlobal$scaDegFreedom
  scaDeltaDegFreedom <- scaDFbyGeneH1 - scaDFbyGeneH0
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH <- p.adjust(vecPvalue, method = "BH")
  
  if(objectLineagePulse@lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
    vecMu <- as.vector(objectLineagePulse@lsMuModelH0$matMuModel)
  } else {
    vecMu <- as.vector(objectLineagePulse@lsMuModelConst$matMuModel)
  }
  dfResults <- data.frame(
    gene=rownames(objectLineagePulse@matCountsProc),
    p=as.numeric(vecPvalue),
    adj.p=as.numeric(vecPvalueBH),
    mean_H0=vecMu,
    dispersion_H0=as.vector(objectLineagePulse@lsDispModelH0$matDispModel),
    loglik_full=vecLogLikFull,
    loglik_red=vecLogLikRed,
    deviance=vecDeviance,
    df_full=scaDFbyGeneH1,
    df_red=scaDFbyGeneH0,
    stringsAsFactors=FALSE)
  rownames(dfResults) <- dfResults$gene
  
  dfResults <- dfResults[match(objectLineagePulse@vecAllGenes,dfResults$gene),]
  rownames(dfResults) <- objectLineagePulse@vecAllGenes
  vecboolAllZero <- !(objectLineagePulse@vecAllGenes %in% rownames(objectLineagePulse@matCountsProc))
  dfResults$allZero <- vecboolAllZero
  
  objectLineagePulse@dfResults <- dfResults
  return(objectLineagePulse)
}