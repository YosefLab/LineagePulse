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
#' \item p: P-value for differential expression with ZINB noise.
#' \item mean: Inferred mean parameter of constant model of first batch.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis with NB noise.
#' \item p_nb: P-value for differential expression with ZINB noise.
#' \item padj_nb: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis with NB noise.
#' \item loglik_full_zinb: Loglikelihood of full model with ZINB noise.
#' \item loglik_red_zinb: Loglikelihood of reduced model with ZINB noise.
#' \item loglik_full_nb: Loglikelihood of full model with NB noise.
#' \item loglik_red_nb: Loglikelihood of reduced model with NB noise.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item allZero (bool) Whether there were no observed non-zero observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#'    
#' @author David Sebastian Fischer
runDEAnalysis <- function(objLP){
    
    ### (I) Compute log likelihoods
    ## ZINB model
    vecLogLikFull_ZINB <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH1(objLP),
        lsDispModel=lsDispModelH1(objLP), 
        lsDropModel=lsDropModel(objLP),
        matWeights=matWeights(objLP))
    vecLogLikRed_ZINB <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH0(objLP),
        lsDispModel=lsDispModelH0(objLP), 
        lsDropModel=lsDropModel(objLP),
        matWeights=matWeights(objLP))
    
    ## NB model
    vecLogLikFull_NB <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH1_NB(objLP),
        lsDropModel=list(lsDropModelGlobal=list(strDropModel="none")),
        lsDispModel=lsDispModelH1_NB(objLP), 
        matWeights=matWeights(objLP))
    vecLogLikRed_NB <- evalLogLikMatrix(
        matCounts=matCountsProc(objLP),
        lsMuModel=lsMuModelH0_NB(objLP),
        lsDispModel=lsDispModelH0_NB(objLP),
        lsDropModel=list(lsDropModelGlobal=list(strDropModel="none")),
        matWeights=matWeights(objLP))
    
    ### (II) Differential expression analysis
    ## ZINB model
    # Compute difference in degrees of freedom 
    # between null model and alternative model.
    scaDFbyGeneH1_ZINB <- lsMuModelH1(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH1(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDFbyGeneH0_ZINB <- lsMuModelH0(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH0(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDeltaDegFreedom_ZINB <- scaDFbyGeneH1_ZINB - scaDFbyGeneH0_ZINB
    # Compute test statistic: Deviance
    vecDeviance_ZINB <- 2*(vecLogLikFull_ZINB - vecLogLikRed_ZINB)
    # Get p-values from Chi-square distribution (assumption about null model)
    vecPvalue_ZINB <- pchisq(vecDeviance_ZINB, 
                             df=scaDeltaDegFreedom_ZINB, lower.tail=FALSE)
    # Multiple testing correction (Benjamini-Hochberg)
    vecPvalueBH_ZINB <- p.adjust(vecPvalue_ZINB, method = "BH")
    
    ## NB model
    scaDFbyGeneH1_NB <- lsMuModelH1(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH1(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDFbyGeneH0_NB <- lsMuModelH0(objLP)$lsMuModelGlobal$scaDegFreedom + 
        lsDispModelH0(objLP)$lsDispModelGlobal$scaDegFreedom
    scaDeltaDegFreedom_NB <- scaDFbyGeneH1_NB - scaDFbyGeneH0_NB
    vecDeviance_NB <- 2*(vecLogLikFull_NB - vecLogLikRed_NB)
    vecPvalue_NB <- pchisq(vecDeviance_NB, 
                           df=scaDeltaDegFreedom_NB, lower.tail=FALSE)
    vecPvalueBH_NB <- p.adjust(vecPvalue_NB, method = "BH")
    
    ### (III) Summarise results
    if(lsMuModelH0(objLP)$lsMuModelGlobal$strMuModel=="constant"){
        vecMu <- as.vector(lsMuModelH0_NB(objLP)$matMuModel)
    } else {
        vecMu <- NA
    }
    dfResults <- data.frame(
        gene=rownames(matCountsProc(objLP)),
        p=vecPvalue_ZINB,
        padj=vecPvalueBH_ZINB,
        mean_H0=vecMu,
        p_nb=vecPvalue_ZINB,
        padj_nb=vecPvalueBH_NB,
        df_full_zinb=scaDFbyGeneH1_ZINB,
        df_red_zinb=scaDFbyGeneH0_ZINB,
        df_full_nb=scaDFbyGeneH1_NB,
        df_red_nb=scaDFbyGeneH0_NB,
        loglik_full_zinb=vecLogLikFull_ZINB,
        loglik_red_zinb=vecLogLikRed_ZINB,
        loglik_full_nb=vecLogLikFull_NB,
        loglik_red_nb=vecLogLikRed_NB,
        stringsAsFactors=FALSE)
    rownames(dfResults) <- dfResults$gene
    
    dfResults <- dfResults[match(vecAllGenes(objLP), dfResults$gene),]
    rownames(dfResults) <- vecAllGenes(objLP)
    vecboolAllZero <- !(vecAllGenes(objLP) %in% rownames(matCountsProc(objLP)))
    dfResults$allZero <- vecboolAllZero
    
    dfResults(objLP) <- dfResults
    return(objLP)
}

#' Test for existance of drop-out with log-likelihood ratio test
#' 
#' Performs one test for entire data set.
#' 
#' @seealso Called by user.
#' 
#' @param objLP (LineagePulseObject)
#' LineagePulseObject with fitted null and alternative models.
#' 
#' @return  (data frame)
#' Summary of hypothesis test
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for existance of drop-out in data set. If high, the
#' data can be explained with models based on NB noise and zero-inflation
#' is not necessary: The null hypothesis of no zero-inflation cannot be
#' rejected.
#' \item loglik_zinb: Loglikelihood of full model with ZINB noise (all genes).
#' \item loglik_nb: Loglikelihood of reduced model with NB noise (all genes).
#' \item df_full: Degrees of freedom of full model with ZINB noise (all genes).
#' \item df_red: Degrees of freedom of reduced model with NB noise (all genes).
#' }
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 100,
#'     scaNConst = 10,
#'     scaNLin = 10,
#'     scaNImp = 10,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 30),
#'     vecGeneWiseDropoutRates = rep(0.1, 30))
#' matDropoutPredictors <- as.matrix(data.frame(
#'     log_means = log(rowMeans(lsSimulatedData$counts)+1) ))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "splines", scaDFSplinesMu = 6,
#'     strDropModel="logistic", 
#'     matPiConstPredictors = matDropoutPredictors)
#' testDropout(objLP)$p
#' 
#' @author David Sebastian Fischer
#' 
#' @export
testDropout <- function(objLP){
    
    ## Degrees of freedom used in mean, dispersion and drop-out models:
    # mean and dispersion models:
    scaDF_ZINB <- sum(objLP$dfResults$df_full_zinb)
    scaDF_NB <- sum(objLP$dfResults$df_full_nb)
    # drop-out models:
    if(lsDropModel(objLP)$lsDropModelGlobal$strDropFitGroup == "PerCell") {
        # one model per cell
        scaDF_ZINB <- scaDF_ZINB +
            dim(lsDropModel(objLP)$matDropoutLinModel)[1] * 
            dim(lsDropModel(objLP)$matDropoutLinModel)[2]
    } else if(lsDropModel(objLP)$lsDropModelGlobal$strDropFitGroup == "AllCells"){
        # one model shared across all cells
        scaDF_ZINB <- scaDF_ZINB +
            dim(lsDropModel(objLP)$matDropoutLinModel)[2]
    }
    scaDeltaDegFreedom_dropout <- scaDF_ZINB - scaDF_NB
    ## The likelihood considered is the likelihood across the entire data set
    scaLogLikAllGenes_ZINB <- sum(objLP$dfResults$loglik_full_zinb)
    scaLogLikAllGenes_NB <- sum(objLP$dfResults$loglik_full_nb)
    scaDeviance_dropout <- 2*(scaLogLikAllGenes_ZINB - scaLogLikAllGenes_NB)
    scaPvalue_dropout <- pchisq(scaDeviance_dropout, 
                                df=scaDeltaDegFreedom_dropout, lower.tail=FALSE)
    
    return(data.frame(
        p=scaPvalue_dropout,
        df_delta=scaDeltaDegFreedom_dropout,
        deviance=scaDeviance_dropout,
        df_full=scaDF_ZINB,
        df_red=scaDF_NB,
        loglik_full=scaLogLikAllGenes_ZINB,
        loglik_red=scaLogLikAllGenes_NB,
        stringsAsFactors=FALSE)
    )
}