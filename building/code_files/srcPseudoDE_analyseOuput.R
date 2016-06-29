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
  strSCMode="batch",
  scaWindowRadis=20){
  
  graphics.off()
  
  # Load data
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
  
  # 1. CDF p-values
  vecX <- seq(-100,-1,by=0.5)
  vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResults$adj.p))/log(10) <= thres, na.rm=TRUE)})
  pdf(paste0(folderPDFs,"/LineagePulse_ECDF-pvalues.pdf"),width=7,height=7)
  plot(vecX,vecCDF1,
    col="black",pch=4,type="l",
    ylim=c(0,max(vecCDF1,na.rm=TRUE)),
    xlab="log_10(p-value)", 
    ylab=paste0("Cumulative p-value distribution"),
    main=paste0("Cumulative p-values distribution LineagePulse"),
    cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  legend(x="topleft",
    legend=c("LineagePulse"),
    fill=c("black"))
  dev.off()
  graphics.off()
  
  scaNumGenes <- dim(matCountDataProc)[1]
  scaNumCells <- dim(matCountDataProc)[2]
  scaNIDs <- 200
  scaFracExprIsHigh <- 0.2
  # 2. Plot highly expressed genes
  vecIDsHighExpr <- apply(matCountDataProc, 1, function(gene){sum(gene > 10) >= scaFracExprIsHigh*scaNumCells})
  vecIDsTopQvalHighExpre <- (dfImpulseResults[rownames(matCountDataProc[vecIDsHighExpr,]),]$Gene)[1:scaNIDs]
  plotDEGenes( vecGeneIDs=vecIDsTopQvalHighExpre,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_LineagePulseFits_HighExpr"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    NPARAM=6)
  
  # 3. Plot top q-val
  vecIDsTopQval <- dfImpulseResults[1:scaNIDs,]$Gene
  plotDEGenes( vecGeneIDs=vecIDsTopQval,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_LineagePulseFits_TopQval"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    NPARAM=6)
  
  # 4. Plot worst significant q-val
  scaThres <- 10^(-3)
  idxIDatThres <- min(which(dfImpulseResults$adj.p > scaThres))-1
  vecIDsWorstQval <- dfImpulseResults[max(1,idxIDatThres-scaNIDs):idxIDatThres,]$Gene
  plotDEGenes( vecGeneIDs=vecIDsWorstQval,
    matCountDataProc=matCountDataProc,
    matTranslationFactors=NULL,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix=paste0("_LineagePulseFits_WorstQval"), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=NULL,
    strNameMethod2=NULL,
    strMode="singlecell",
    strSCMode="continuous",
    NPARAM=6)
  
  # 5. Look at drop out rate vs mean fitting
  vecSumNotDropout <- apply(matDropout, 2, function(cell){sum(1-cell, na.rm=TRUE)})
  vecSumProbNBRates <- apply(matProbNB, 2, function(cell){sum(cell, na.rm=TRUE)})
  vecDepth <- apply(matCountsProcLP, 2, function(cell){sum(cell, na.rm=TRUE)})
  plot(log(vecDepth), vecSumNotDropout)
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
  pdf(paste0("CumulNBProb-Depth_scatter.pdf"),width=7,height=7)
  print(g)
  dev.off()
  graphics.off()
  
  # Plot full lines coloured by depth
  g <- ggplot() + geom_line()
  for(cell in 1:scNumCells){
    dfLogisticScatter <- data.frame(
      x=log(matMu[,cell]),
      y=log(matDropout[,cell]))
    g <- g + geom_point(data=dfLogisticScatter, aes(x=x, y=y),  size=4)
  }
  print(g)
  
  return(NULL)
}