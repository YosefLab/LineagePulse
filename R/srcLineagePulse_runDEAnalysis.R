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
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for differential expression.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis.
#' \item loglik_full: Loglikelihood of full model.
#' \item loglik_red: Loglikelihood of reduced model.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item mean: Inferred mean parameter of constant model of first batch.
#' From combined samples in case-ctrl. 
#' \item allZero (bool) Whether there were no observed non-zero observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#' Entries only present in case-only DE analysis:
#' \itemize{
#' \item converge_impulse: Convergence status of optim for 
#' impulse model fit (full model).
#' \item converge_const: Convergence status of optim for 
#' constant model fit (reduced model).
#' }
#' Entries only present in mixture model DE analysis:
#' \itemize{
#' \item converge_mixture: Convergence status of optim for 
#' mixture model fit (full model).
#' \item converge_const: Convergence status of optim for 
#' constant model fit (reduced model).
#' }
#'    
#' @author David Sebastian Fischer
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