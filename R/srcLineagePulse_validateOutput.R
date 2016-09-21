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
#' @seealso Called by \code{anlayseOuput}.
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++  Validation metrics for ZINB fits  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Generate metrics for validation of ZINB fit
#' 
#' Compute and plot summary metrics that allow the user to judge
#' whether the ZINB fit is sensible.
#' 
#' @seealso Called by \code{runLineagePulse} or separately by user.
#' 
#' @param matQval: matrix with p-values by methods
#' 
#' @return NULL
#' 
#' @export

validateOuput <- function(
  dirLineagePulseTempFiles,
  dirValidationOut,
  strSCMode="continuous",
  scaWindowRadis=20,
  dfGeneAnnotation=NULL){
  
  graphics.off()
  
  # Load data
  setwd(folderLineagePulseOutput)
  load("LineagePulse_matCountsProc.RData")
  load("LineagePulse_matMuH0.RData")
  load("LineagePulse_matDispersionsH0.RData")
  load("LineagePulse_matDropoutH0.RData")
  load("LineagePulse_matMuH1.RData")
  load("LineagePulse_matDispersionsH1.RData")
  load("LineagePulse_matDropoutH1.RData")
  load("LineagePulse_matZH1.RData")
  load("LineagePulse_matDropoutLinModel.RData")
  load("LineagePulse_dfDEAnalysis.RData")
  
  # Initialise
  setwd(dirValidationOut)
  scaNumGenes <- dim(matCountDataProc)[1]
  scaNumCells <- dim(matCountDataProc)[2]
  
  # 1. ECDF q-values
  vecX <- seq(max(-100,min(dfDEAnalysis$adj.p)),0,by=0.5)
  vecCDF1 <- sapply(vecX, function(thres){
    sum( log(as.numeric(as.vector(dfDEAnalysis$adj.p)))/log(10) <= thres, na.rm=TRUE)})
  pdf("LineagePulse_ECDF-qvalues.pdf",width=7,height=7)
  plot(vecX,vecCDF1,
    col="black",pch=4,type="l",
    ylim=c(0,max(vecCDF1,na.rm=TRUE)),
    xlab="log_10 q-value", 
    ylab=paste0("Cumulative q-value distribution"),
    main=paste0("Cumulative q-values distribution LineagePulse"),
    cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  legend(x="topleft",
    legend=c("LineagePulse"),
    fill=c("black"))
  dev.off()
  graphics.off()

  scaNIDs <- 200
  scaFracExprIsHigh <- 0.5
  
  # 2. Plot Fits
  # a). Plot highly expressed genes
  vecIDsHighExpr <- apply(matCountDataProc, 1, function(gene){
    sum(gene > 10) >= scaFracExprIsHigh*scaNumCells})
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
  
  # b) Plot top q-val
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
  
  # c) Plot worst significant q-val
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
  
  # 3. Look at drop out rate vs mean fitting
  # Sequencing depth as a comparative experimental measure for drop-out
  vecDepth <- apply(matCountsProcLP, 2, function(cell){sum(cell, na.rm=TRUE)})
  
  # a) Sum of drop out rates of a cell
  vecSumNotDropout <- apply(matDropout, 2, function(cell){sum(1-cell, na.rm=TRUE)})
  
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
  pdf("LineagePulse_Scatter_CumulativeNBProbvsDepth.pdf",width=7,height=7)
  print(g1)
  dev.off()
  # b) AUC of logistic curve as measure for drop-out intensity
  vecAUCLogistic <- apply(matDropoutLinModel, 1, function(cellmodel){
    computeAUCLogistic(cellmodel)
  })
  
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
  pdf("LineagePulse_Scatter_AUCDropoutvsDepth.pdf",width=7,height=7)
  print(g2)
  dev.off()
  graphics.off()
  
  # c) Inflexion point of logistic
  vecInflexLogistic <- apply(matDropoutLinModel, 1, function(cellmodel){
    return(-cellmodel[1]/cellmodel[2])
  })
  
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
  pdf("LineagePulse_Scatter_InflexionDropoutvsDepth.pdf",width=7,height=7)
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
  pdf("LineagePulse_Scatter_logInflexionDropoutvslogDepth.pdf",width=7,height=7)
  print(g3b)
  dev.off()
  graphics.off()
  
  # 4. Plot drop-out rate as function of cell
  pdf("LineagePulse_LogisticScatter_DropoutvsMeanbyCell.pdf",width=7,height=7)
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
  
  # 5. Logistic dropout rate fit as function of mean coloured by
  # sequencing depth.
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
  pdf("LineagePulse_AllCurves_LogisticDropoutvsMubyDepth.pdf",width=7,height=7)
  print(g4)
  dev.off()
  graphics.off()
  
  # 6. Plot variance measure as function of mean
  # a) Plot dispersion parameter fit as function of mean parameter
  # From H0 to have one parameter each per gene.
  dfScatterDispvsMean <- data.frame(
    x=log(matMuH0[,1])/log(10),
    y=matDispersions[,1]) )
  gScatterDispvsMean <- ggplot(dfScatterDispvsMean, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Dispersion parameter as function of mean parameter under H0.") +
    xlab(paste0("log10 H0 mean parameter")) +
    ylab(paste0("H0 dispersion parameter")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_Scatter_DispvsMean.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  graphics.off()
  
  # b) Plot coefficient of variation as function of mean parameter
  # Compute coefficient of variation under H0
  vecVar <- matMuH0[,1] + matMuH0[,1]^2 / matDispersions[,1]
  vecCV <- vecVar / matMuH0[,1]^2
  # From H0 to have one parameter each per gene.
  dfScatterDispvsMean <- data.frame(
    x=log(matMuH0[,1])/log(10),
    y=vecCV )
  gScatterDispvsMean <- ggplot(dfScatterDispvsMean, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Dispersion parameter as function of mean parameter under H0.") +
    xlab(paste0("log10 H0 mean parameter")) +
    ylab(paste0("coefficient of variation")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_Scatter_CVvsMean.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  graphics.off()
  
  return(NULL)
}