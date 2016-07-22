#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Analyse LineagePulse output    ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(gplots)
library(ggplot2)
library(reshape2)

#' Compute AUC of logistic function
#' 
#' This is a measure for the intensity of the predicted drop-out effect in the
#' given cell. Assumes that drop-out rate(mean) relationship in 
#' monotonously decreasing.
#' 
#' @seealso Called by axuillary \code{anlayseOuput}.
#' 
#' @param vecLinearModel: (numeric vector length 2)
#'    Two parameters of linear model of logistic function.
#' 
#' @return scaAUC: (scalar) Area under the logistic curve
#'    from 0 to Inf.
#' 
#' @export

computeAUCLogistic <- function(vecLinearModel){
  
  # The logistic function f is
  # f(x) = 1/(1+e^(-(a1+a2*x))
  # where a1,a2 is the linear model.
  # Note that in our case a1 and a2 are always negative
  # because the probability of drop-out decreases with
  # increasing mean. 
  # f*(x) = 1/(1+e^(+(-a1-a2*x))
  # where -a1-a2*x >== 0 forall x >= 0
  # The integral over f* is
  # F*(x) = -ln(1+e^(-(-a1-a2*x)))+C 
  # Accordingly, we compute the AUC = F*(x=Inf) - F(x=0)
  # where F*(x=Inf) = 0 and F(x=0) is to be evaluted:
  scaAUC <- log(1+exp(vecLinearModel[1]))
  
  # Catch exception in which a1,a2 are not <0:
  if(vecLinearModel[1]>0 | vecLinearModel[2]>0){
    scaAUC <- NA
  }
  
  return(scaAUC)
}

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
  load("PseudoDE_matDropoutLinModel.RData")
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
  # a) Sum of drop out rates of a cell
  vecSumNotDropout <- apply(matDropout, 2, function(cell){sum(1-cell, na.rm=TRUE)})
  # b) AUC of logistic curve as measure for drop-out intensity
  vecAUCLogistic <- apply(matDropoutLinModel, 1, function(cellmodel){
    computeAUCLogistic(cellmodel)
  })
  # c) Inflexion point of logistic
  vecInflexLogistic <- apply(matDropoutLinModel, 1, function(cellmodel){
    return(-cellmodel[1]/cellmodel[2])
  })
  # Sequencing depth as a comparative experimental measure for drop-out
  vecDepth <- apply(matCountsProcLP, 2, function(cell){sum(cell, na.rm=TRUE)})
  
  dfScatter1 <- data.frame(
    x=log(vecDepth),
    y=log(vecSumNotDropout))
  g1 <- ggplot(dfScatter1, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot cumulative NB mixture probability versus sequencing depth by cell") +
    xlab(paste0("log sequencing depth")) +
    ylab(paste0("log cumulative NB mixture probability")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_Scatter_CumulativeNBProbvsDepth.pdf"),width=7,height=7)
  print(g1)
  dev.off()
  
  dfScatter2 <- data.frame(
    x=log(vecDepth),
    y=log(vecAUCLogistic))
  g2 <- ggplot(dfScatter2, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot AUC logistic drop-out model versus sequencing depth by cell") +
    xlab(paste0("log sequencing depth")) +
    ylab(paste0("AUC logistic drop-out model")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_Scatter_AUCDropoutvsDepth.pdf"),width=7,height=7)
  print(g2)
  dev.off()
  graphics.off()
  
  dfScatter3 <- data.frame(
    x=vecDepth,
    y=vecInflexLogistic)
  g3 <- ggplot(dfScatter3, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot inflexion point logistic drop-out model versus sequencing depth by cell") +
    xlab(paste0("Sequencing depth")) +
    ylab(paste0("Inflexion point logistic drop-out model")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_Scatter_InflexionDropoutvsDepth.pdf"),width=7,height=7)
  print(g3)
  dev.off()
  graphics.off()
  
  dfScatter3b <- data.frame(
    x=log(vecDepth),
    y=log(vecInflexLogistic-min(vecInflexLogistic)+1))
  g3b <- ggplot(dfScatter3b, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot inflexion point logistic drop-out model versus sequencing depth by cell") +
    xlab(paste0("log Sequencing depth")) +
    ylab(paste0("log Inflexion point logistic drop-out model (adjusted)")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_Scatter_logInflexionDropoutvslogDepth.pdf"),width=7,height=7)
  print(g3b)
  dev.off()
  graphics.off()
  
  # Plot drop-out rate as function of cell
  pdf(paste0(folderPDFs,"/LineagePulse_LogisticScatter_DropoutvsMeanbyCell.pdf"),width=7,height=7)
  for(cell in seq(1,dim(matDropout)[2])){
    plot( log(matMu[,cell])/log(10), 
      matDropout[,cell], 
      xlim=c(-4,max(log(matMu[,cell])/log(10))),
      ylab="Dropout rate",
      xlab="log_10 mu",
      main=paste0("cell_",cell," (Sequencing depth: ",round(vecDepth[cell]/10^5,1),"e5)"))
  }
  dev.off()
  graphics.off()
  # Summary plot with all logistic curves in one plot
  # Colour by sequencing depth ordering
  # Draw samples from dropout logistic model to save memory
  vecMu=sapply(seq(1,log(max(matMu)/10^(-4))/log(2)), function(x) 10^(-4)*2^x)
  matDropoutSampled <- matrix(NA, nrow=length(vecMu), ncol=dim(matDropout)[2])
  for(cell in seq(1,dim(matDropout)[2])){
    vecCellModel <- matDropoutLinModel[cell,]
    matDropoutSampled[,cell] <- 1/(1+exp(-(vecCellModel[1]+vecCellModel[2]*vecMu)))
  }
  colnames(matDropoutSampled) <- colnames(matDropout)
  # Reshape matrices into single column first:
  # Look at first scaNSummaryCells cells
  scaNSummaryCells <- dim(matMu)[2]
  dfMuMolten <- melt(matMu[,1:scaNSummaryCells])
  dfDropoutMolten <- melt(matDropout[,1:scaNSummaryCells])
  dfLines1 <- data.frame(
    mu=log(dfMuMolten$value),
    dropout=dfDropoutMolten$value,
    depth=vecDepth[match(dfDropoutMolten$Var2,names(vecDepth))])
  
  dfDropoutSampledMolten <- melt(matDropoutSampled[,1:scaNSummaryCells])
  dfLines1 <- data.frame(
    mu=log(vecMu)/log(10),
    dropout=dfDropoutSampledMolten$value,
    depth=vecDepth[match(dfDropoutSampledMolten$Var2,names(vecDepth))])
  g4 <- ggplot(dfLines1, aes(x=mu, y=dropout, group=depth)) +
    geom_line(aes(colour=depth)) +
    scale_colour_gradient(low="red",high="green") +
    labs(title="Logistic drop-out model versus sequencing depth by cell") +
    xlab(paste0("log_10 mean parameter")) +
    ylab(paste0("Drop-out rate")) +
    xlim(c(-1,2))
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf(paste0(folderPDFs,"/LineagePulse_AllCurves_LogisticDropoutvsMubyDepth.pdf"),width=7,height=7)
  print(g4)
  dev.off()
  graphics.off()
  
  return(NULL)
}