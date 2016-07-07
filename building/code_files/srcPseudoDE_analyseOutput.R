#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Analyse LineagePulse output    ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(gplots)
library(ggplot2)
library(reshape2)

#' Compare ImpulseDE2 output against other differential expression method
#' 
#' The comparison is performed based on the adjusted p-values of differential
#' expression. This methods allows the user to explore visually how 
#' ImpulseDE2 output differs from similar methods.
#' 
#' @seealso Auxillary method not called during ImpulseDE2 running.
#' Called separately by user.
#' 
#' @param matQval: matrix with p-values by methods
#' 
#' @return NULL
#' 
#' @export

anlayseOuput <- function(
  folderLineagePulseOutput,
  folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20,
  dfGeneAnnotation=NULL){
  
  graphics.off()
  
  # Load data
  setwd(folderLineagePulseOutput)
  load("PseudoDE_matCountsProc.RData")
  load("PseudoDE_matDropout.RData")
  load("PseudoDE_matProbNB.RData")
  load("PseudoDE_matMu.RData")
  matCountsProcLP <- matCountsProc
  
  load("ImpulseDE2_matCountDataProc.RData")
  load("ImpulseDE2_dfAnnotationProc.RData")
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_vecDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_matSizeFactors.RData")
  
  if(!is.null(dfGeneAnnotation)){
    rownames(matCountsProcLP) <- as.vector(dfGeneAnnotation[match(rownames(matCountsProcLP), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(matDropout) <- as.vector(dfGeneAnnotation[match(rownames(matDropout), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(matProbNB) <- as.vector(dfGeneAnnotation[match(rownames(matProbNB), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(matMu) <- as.vector(dfGeneAnnotation[match(rownames(matMu), rownames(dfGeneAnnotation)),"gene_short_name"])
    
    rownames(matCountDataProc) <- as.vector(dfGeneAnnotation[match(rownames(matCountDataProc), rownames(dfGeneAnnotation)),"gene_short_name"])
    dfImpulseResults$Gene <- as.vector(dfGeneAnnotation[match(rownames(dfImpulseResults), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(dfImpulseResults) <- as.vector(dfGeneAnnotation[match(rownames(dfImpulseResults), rownames(dfGeneAnnotation)),"gene_short_name"])
    names(vecDEGenes) <- as.vector(dfGeneAnnotation[match(names(vecDEGenes), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(lsImpulseFits$parameters_case) <- as.vector(dfGeneAnnotation[match(rownames(lsImpulseFits$parameters_case), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(lsImpulseFits$values_case) <- as.vector(dfGeneAnnotation[match(rownames(lsImpulseFits$values_case), rownames(dfGeneAnnotation)),"gene_short_name"])
    rownames(matSizeFactors) <- as.vector(dfGeneAnnotation[match(rownames(matSizeFactors), rownames(dfGeneAnnotation)),"gene_short_name"])
  }
  
  # Plot
  setwd(folderPDFs)
  # 1. CDF q-values
  vecX <- seq(-100,-1,by=0.5)
  vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResults$p)))/log(10) <= thres, na.rm=TRUE)})
  pdf(paste0(folderPDFs,"/LineagePulse_ECDF-pvalues.pdf"),width=7,height=7)
  plot(vecX,vecCDF1,
    col="black",pch=4,type="l",
    ylim=c(0,max(vecCDF1,na.rm=TRUE)),
    xlab="log_10 p-value", 
    ylab=paste0("Cumulative p-value distribution"),
    main=paste0("Cumulative p-values distribution LineagePulse"),
    cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  legend(x="topleft",
    legend=c("LineagePulse"),
    fill=c("black"))
  dev.off()
  graphics.off()
  
  # Additional: histogram p-value distribution
  pdf(paste0(folderPDFs,"/LineagePulse_Histogram-logpvalues.pdf"),width=7,height=7)
  hist(log(as.numeric(as.vector(dfImpulseResults$p)))/log(10),
    breaks=20,
    xlab="log_10 adjusted p-value)", 
    ylab=paste0("Number of observations"),
    main=paste0("Histogram of adjusted p-value"),
    cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  dev.off()
  graphics.off()
  pdf(paste0(folderPDFs,"/LineagePulse_Histogram-pvalues.pdf"),width=7,height=7)
  hist(as.numeric(as.vector(dfImpulseResults$p)),
    breaks=20,
    xlab="log_10 adjusted p-value)", 
    ylab=paste0("Number of observations"),
    main=paste0("Histogram of adjusted p-value"),
    cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  dev.off()
  graphics.off()
  
  scaNumGenes <- dim(matCountDataProc)[1]
  scaNumCells <- dim(matCountDataProc)[2]
  scaNIDs <- 200
  scaFracExprIsHigh <- 0.5
  # 2. Plot highly expressed genes
  vecIDsHighExpr <- apply(matCountDataProc, 1, function(gene){sum(gene > 10) >= scaFracExprIsHigh*scaNumCells})
  vecIDsTopQvalHighExpre <- as.vector(dfImpulseResults[match(rownames(dfImpulseResults), rownames(matCountDataProc[vecIDsHighExpr,])),]$Gene)
  vecIDsTopQvalHighExpre <- (vecIDsTopQvalHighExpre[!is.na(vecIDsTopQvalHighExpre)])[1:scaNIDs]
  plotDEGenes( vecGeneIDs=vecIDsTopQvalHighExpre,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_ImpulseTraces_HighlyExpressedLowQval"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    boolLogPlot=TRUE,
    NPARAM=6)
  graphics.off()
  
  # 3. Plot top q-val
  vecIDsTopQval <- as.vector(dfImpulseResults[1:scaNIDs,]$Gene)
  plotDEGenes( vecGeneIDs=vecIDsTopQval,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_LineagePulse_ImpulseTraces_LowQval"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    boolLogPlot=TRUE,
    NPARAM=6)
  graphics.off()
  
  # 4. Plot worst significant q-val
  scaThres <- 10^(-3)
  idxIDatThres <- min(which(dfImpulseResults$adj.p > scaThres))-1
  vecIDsWorstQval <- as.vector(dfImpulseResults[max(1,idxIDatThres-scaNIDs):idxIDatThres,]$Gene)
  plotDEGenes( vecGeneIDs=vecIDsWorstQval,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_LineagePulse_ImpulseTraces_HighButSignificantQval"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    boolLogPlot=TRUE,
    NPARAM=6)
  graphics.off()
  
  # 5. Look at drop out rate vs mean fitting
  vecSumNotDropout <- apply(matDropout, 2, function(cell){sum(1-cell, na.rm=TRUE)})
  vecSumProbNBRates <- apply(matProbNB, 2, function(cell){sum(cell, na.rm=TRUE)})
  vecDepth <- apply(matCountsProcLP, 2, function(cell){sum(cell, na.rm=TRUE)})
  dfScatter <- data.frame(
    x=log(vecDepth),
    y=log(vecSumNotDropout))
  g <- ggplot(dfScatter, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot cumulative NB mixture probability versus sequencing depth by cell") +
    xlab(paste0("log cumulative NB mixture probability")) +
    ylab(paste0("log sequencing depth")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_Scatter_CumulativeNBProbvsDepth.pdf"),width=7,height=7)
  print(g)
  dev.off()
  graphics.off()
  
  return(NULL)
}