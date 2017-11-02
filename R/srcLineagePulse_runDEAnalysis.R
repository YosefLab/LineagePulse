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
#' @param objLP (LineagePulseObject)
#' LineagePulseObject with fitted null and alternative models.
#' 
#' @return objLP (LineagePulseObject)
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
    vecLogLikFull <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH1(objLP),
        lsDispModel=lsDispModelH1(objLP), 
        lsDropModel=lsDropModel(objLP),
        matWeights=matWeights(objLP))
    vecLogLikRed <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH0(objLP),
        lsDispModel=lsDispModelH0(objLP), 
        lsDropModel=lsDropModel(objLP),
        matWeights=matWeights(objLP))
    
    # (II) Differential expression analysis
    # Compute difference in degrees of freedom 
    # between null model and alternative model.
    scaDFbyGeneH1 <- lsMuModelH1(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH1(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDFbyGeneH0 <- lsMuModelH0(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH0(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDeltaDegFreedom <- scaDFbyGeneH1 - scaDFbyGeneH0
    # Compute test statistic: Deviance
    vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
    # Get p-values from Chi-square distribution (assumption about null model)
    vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
    # Multiple testing correction (Benjamini-Hochberg)
    vecPvalueBH <- p.adjust(vecPvalue, method = "BH")
    
    if(lsMuModelH0(objLP)$lsMuModelGlobal$strMuModel=="constant"){
        vecMu <- as.vector(lsMuModelH0(objLP)$matMuModel)
    } else {
        vecMu <- as.vector(lsMuModelConst(objLP)$matMuModel)
    }
    dfResults <- data.frame(
        gene=rownames(matCountsProc(objLP)),
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
    
    dfResults <- dfResults[match(vecAllGenes(objLP), dfResults$gene),]
    rownames(dfResults) <- vecAllGenes(objLP)
    vecboolAllZero <- !(vecAllGenes(objLP) %in% rownames(matCountsProc(objLP)))
    dfResults$allZero <- vecboolAllZero
    
    dfResults(objLP) <- dfResults
    return(objLP)
}