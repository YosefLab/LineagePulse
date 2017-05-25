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
#' @param: objLP: (LineagePulseObject)
#' LineagePulseObject with fitted null and alternative models.
#' 
#' @return objLP: (LineagePulseObject)
#' LineagePulseObject with analysis summary table (dfResults)
#' added.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
runDEAnalysis <- function(objLP){
  
  # (I) Compute log likelihoods
  vecLogLikFull <- evalLogLikMatrix(matCounts=objLP@matCountsProc,
                                lsMuModel=objLP@lsMuModelH1,
                                lsDispModel=objLP@lsDispModelH1, 
                                lsDropModel=objLP@lsDropModel,
                                matWeights=objLP@matWeights)
  vecLogLikRed <- evalLogLikMatrix(matCounts=objLP@matCountsProc,
                                   lsMuModel=objLP@lsMuModelH0,
                                   lsDispModel=objLP@lsDispModelH0, 
                                   lsDropModel=objLP@lsDropModel,
                                   matWeights=objLP@matWeights)

  # (II) Differential expression analysis
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDFbyGeneH1 <- objLP@lsMuModelH1$lsMuModelGlobal$scaDegFreedom + 
    objLP@lsDispModelH1$lsDispModelGlobal$scaDegFreedom
  scaDFbyGeneH0 <- objLP@lsMuModelH0$lsMuModelGlobal$scaDegFreedom + 
    objLP@lsDispModelH0$lsDispModelGlobal$scaDegFreedom
  scaDeltaDegFreedom <- scaDFbyGeneH1 - scaDFbyGeneH0
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH <- p.adjust(vecPvalue, method = "BH")
  
  if(objLP@lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
    vecMu <- as.vector(objLP@lsMuModelH0$matMuModel)
  } else {
    vecMu <- as.vector(objLP@lsMuModelConst$matMuModel)
  }
  dfResults <- data.frame(
    gene=rownames(objLP@matCountsProc),
    p=as.numeric(vecPvalue),
    adj.p=as.numeric(vecPvalueBH),
    mean_H0=vecMu,
    loglik_full=vecLogLikFull,
    loglik_red=vecLogLikRed,
    deviance=vecDeviance,
    df_full=scaDFbyGeneH1,
    df_red=scaDFbyGeneH0,
    stringsAsFactors=FALSE)
  rownames(dfResults) <- dfResults$gene
  
  dfResults <- dfResults[match(objLP@vecAllGenes,dfResults$gene),]
  rownames(dfResults) <- objLP@vecAllGenes
  vecboolAllZero <- !(objLP@vecAllGenes %in% rownames(objLP@matCountsProc))
  dfResults$allZero <- vecboolAllZero
  
  objLP@dfResults <- dfResults
  return(objLP)
}