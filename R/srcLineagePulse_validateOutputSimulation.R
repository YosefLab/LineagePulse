#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++     Analyse LineagePulse output on simulated data given model    ++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(gplots)
library(ggplot2)
library(reshape2)

#' Generate metrics for comparison of ZINB fit against underlying model
#' 
#' Compare inferred parameters with underlying model parameters.
#' 
#' @seealso Auxillary method not called during ImpulseDE2 running.
#' Called separately by user.
#' 
#' @param matQval: matrix with p-values by methods
#' 
#' @return NULL
#' 
#' @export

validateOuput <- function(
  dirLineagePulseTempFiles,
  dirValidationOut,
  matMuHidden,
  matDispHidden,
  matDropoutHidden,
  matDropoutLinModelHidden,
  vecSizeFactorsHidden,
  vecConstIDs,
  vecImpulseIDs,
  scaWindowRadis=NULL ){
  
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
  
  vecAnalysedGenes <- rownames(matCountsProc)
  vecAnalysedCells <- colnames(matCountsProc)
  matMuHidden <- matMuHidden[vecAnalysedGenes,vecAnalysedCells]
  matDispHidden <- matDispHidden[vecAnalysedGenes,vecAnalysedCells]
  matDropoutHidden <- matDropoutHidden[vecAnalysedGenes,vecAnalysedCells]
  matDropoutLinModelHidden <- matDropoutLinModelHidden[vecAnalysedCells,]
  vecSizeFactorsHidden <- vecSizeFactorsHidden[vecAnalysedCells]
  
  # Initialise
  setwd(dirValidationOut)
  scaNumGenes <- dim(matCountDataProc)[1]
  vecConstIDs <- vecConstIDs[vecConstIDs %in% rownames(matCountDataProc)]
  vecImpulseIDs <- vecImpulseIDs[vecImpulseIDs %in% rownames(matCountDataProc)]
  scaNumGenesConst <- length(vecConstIDs)
  scaNumGenesImpulse <- length(vecImpulseIDs)
  scaNumCells <- dim(matCountDataProc)[2]
  
  # 1.) Scatter plots of parameter sets
  # a) Mean parameters H0
  dfScatterMuModelvsMuInferredH0 <- data.frame(
    x=log(matMuH0[,1])/log(10),
    y=log(apply(matMuHidden,1, function(gene) mean(gene, na.rm=TRUE)) )/log(10) )
  gScatterDispvsMean <- ggplot(dfScatterMuModelvsMuInferredH0, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Inferred mean (H0) versus model mean parameter.") +
    xlab(paste0("log10 H0 inferred mean parameter")) +
    ylab(paste0("log10 gene-wise average of model mean parameter")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_MuModelvsMuInferredH0.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  graphics.off()
  
  # b) Mean parameters H1: For one cell only
  dfScatterMuModelvsMuInferredH1 <- data.frame(
    x=log(matMuH1[,1])/log(10),
    y=log(matMuHidden[,1])/log(10) )
  gScatterDispvsMean <- ggplot(dfScatterMuModelvsMuInferredH1, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Inferred mean (H1) of first cell versus model mean parameter.") +
    xlab(paste0("log10 H0 inferred mean parameter")) +
    ylab(paste0("log10 model mean parameter")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_MuModelvsMuInferredH1.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  graphics.off()
  
  # c) Dispersion parameters H0 (using dispersion of first cell from model)
  dfScatterDispModelvsDispInferredH0 <- data.frame(
    x=matDispersionH0[,1],
    y=matDispersionHidden[,1] )
  gScatterDispvsMean <- ggplot(dfScatterDispModelvsDispInferredH0, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Inferred dispersion (H0) versus model dispersion parameter.") +
    xlab(paste0("H0 inferred dispersion parameter")) +
    ylab(paste0("model dispersion parameter of first cell")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_DispModelvsDispInferredH0.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  
  # d) Mean parameters H1: For one cell only
  dfScatterDispModelvsDispInferredH1 <- data.frame(
    x=matDispersionH1[,1],
    y=matDispersionHidden[,1] )
  gScatterDispvsMean <- ggplot(dfScatterDispModelvsDispInferredH1, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Inferred dispersion (H0) versus model dispersion parameter.") +
    xlab(paste0("H1 inferred dispersion parameter of first cell")) +
    ylab(paste0("model dispersion parameter of first cell")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_DispModelvsDispInferredH1.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  
  # e) Dropout model
  # PDF with one page per cell: Plot data and inferred and true logistic model
  
  # 2. Mean model by gene
  # Plot data, true and inferred model by gene. qvalue
  
  # 3. LRT hist
  # Overall deviation comparison: LRT
  # Compare under negative binomial and under ZINB model
  # Compute loglikelihood of true underlying model
  # Recover zero predictions:
  matMuHiddenTemp <- matMuHidden
  matMuHiddenTemp[matMuHiddenTemp==0] <- 10^(-4)
  vecNBLLHiddenModel <- sapply(seq(1,scaNumGenes), function(gene){
    sum(dnbinom(x=matCountsProc[gene,],
      mu=matMuHiddenTemp[gene,],
      size=matDispHidden[gene,],
      log=TRUE), na.rm=TRUE)
  })
  vecZINBLLHiddenModel <- unlist(lapply( seq(1,scaNumGenes), function(i){
    evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
      vecMu=matMuHiddenTemp[i,],
      vecSizeFactors=vecSizeFactorsHidden,
      vecDispEst=matDispHidden[i,], 
      vecDropoutRateEst=matDropoutHidden[i,],
      vecboolNotZeroObserved=matData[i,]>0 & !is.na(matData[i,]), 
      vecboolZero=matData[i,]==0,
      scaWindowRadius=scaWindowRadius)
  }))
  # Compute loglikelihood of inferred model
  # Recover zero predictions:
  matMuTemp <- matMuH1
  matMuTemp[matMuTemp==0] <- 10^(-4)
  vecNBLLH1Model <- sapply(seq(1,scaNumGenes), function(gene){
    sum(dnbinom(x=matCountsProc[gene,],
      mu=matMuTemp[gene,],
      size=matDispersionsH1[gene,],
      log=TRUE), na.rm=TRUE)
  })
  vecZINBLLH1Model <- unlist(lapply( seq(1,scaNumGenes), function(i){
    evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
      vecMu=matMuTemp[i,],
      vecSizeFactors=vecSizeFactorsHidden,
      vecDispEst=matDispersionsH1[i,], 
      vecDropoutRateEst=matDropoutH1[i,],
      vecboolNotZeroObserved=matData[i,]>0 & !is.na(matData[i,]), 
      vecboolZero=matData[i,]==0,
      scaWindowRadius=scaWindowRadius)
  }))
  # LRT
  scaQThres <- 10^(-3)
  vecLRT_NB <- vecNBLLH1Model-vecNBLLHiddenModel
  vecLRT_ZINB <- vecZINBLLH1Model-vecZINBLLHiddenModel
  dfLRT_NB <- data.frame(value=vecLRT_NB)
  dfLRT_ZINB <- data.frame(value=vecLRT_ZINB)
  gHistLRT_NB <- ggplot( dfLRT_NB, aes(value)) +
    geom_histogram(alpha = 1, position = 'identity') +
    ggtitle(paste0("log LRT value under NB model")) +
    xlab("Difference in loglikelihood") +
    ylab("Frequency")
  graphics.off()
  pdf("LineagePulse_AnalyseSimulated_Hist_LRT_InferredH1vsHidden_NB.pdf")
  print(gHistLRT)
  dev.off()
  gHistLRT_ZINB <- ggplot( dfLRT_ZINB, aes(value)) +
    geom_histogram(alpha = 1, position = 'identity') +
    ggtitle(paste0("log LRT value under ZINB model")) +
    xlab("Difference in loglikelihood") +
    ylab("Frequency")
  graphics.off()
  pdf("LineagePulse_AnalyseSimulated_Hist_LRT_InferredH1vsHidden_ZINB.pdf")
  print(gHistLRT_ZINB)
  dev.off() 
  
  # 4. ECDF q-values divided by const and impulse true model
  vecX <- seq(max(-100,min(dfDEAnalysis$adj.p)),0,by=0.5)
  vecCDFConst <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecConstIDs,]$adj.p)))/log(10) <= thres, na.rm=TRUE)})
  vecCDFImpulse <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecImpulseIDs,]$adj.p)))/log(10) <= thres, na.rm=TRUE)})
  dfECDFQvalByModel <- data.frame( thres=vecX,
    qval=c(vecCDFConst,vecCDFImpulse),
    model=c(rep("const",length(vecCDFConst)), rep("impulse",length(vecCDFImpulse))) )
  gECDFQvalByModel <- ggplot(dfECDFQvalByModel, aes(x=thres, y=qval, group=model)) +
    geom_line(aes(colour = group))
  pdf("LineagePulse_AnalyseSimulated_ECDF_QvalByModelClass.pdf",width=7,height=7)
  print(gECDFQvalByModel)
  dev.off() 
  
  # 5. Q-value as function of model
  vecHiddenModelType <- rep("const", length(vecAnalysedGenes))
  names(vecHiddenModelType) <- vecAnalysedGenes
  vecHiddenModelType[vecImpulseIDs] <- "impulse"
  # a) Q-value as function of average true mean parameter
  # divided by const and impulse true model
  dfScatterQvalvsMuModelByModel <- data.frame(
    mu=apply(matMuHidden, 1, function(gene) mean(gene, na.rm=TRUE)),
    qval=dfDEAnalysis[vecAnalysedGenes,]$adj.p,
    model=vecHiddenModelType )
  gScatterQvalvsMean <- ggplot(dfScatterQvalvsMuModelByModel, aes(x=mu, y=qval, group=model)) +
    geom_point(aes(colour = group)) + 
    labs(title="Q-value versus true by underlying model type") +
    xlab(paste0("average true mean")) +
    ylab(paste0("log10 Q-value")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_QvalvsMuModelByModel.pdf",width=7,height=7)
  print(gScatterQvalvsMean)
  dev.off()
  
  # b) Q-value as function of true dispersion parameter
  # divided by const and impulse true model
  dfScatterQvalvsDispModelByModel <- data.frame(
    disp=apply(matDispHidden, 1, function(gene) mean(gene, na.rm=TRUE)),
    qval=dfDEAnalysis[vecAnalysedGenes,]$adj.p,
    model=vecHiddenModelType )
  gScatterQvalvsDisp <- ggplot(dfScatterQvalvsDispModelByModel, aes(x=disp, y=qval, group=model)) +
    geom_point(aes(colour = group)) + 
    labs(title="Q-value versus true by underlying model type") +
    xlab(paste0("average true dispersion")) +
    ylab(paste0("log10 Q-value")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_QvalvsDispModelByModel.pdf",width=7,height=7)
  print(gScatterQvalvsDisp)
  dev.off()
  
  # c) Q-value as function of number of drop-outs
  # divided by const and impulse true model
  dfScatterQvalvsDispModelByModel <- data.frame(
    sumdroprate=apply(matDropoutHidden, 1, function(gene) sum(gene, na.rm=TRUE)),
    qval=dfDEAnalysis[vecAnalysedGenes,]$adj.p,
    model=vecHiddenModelType )
  gScatterQvalvsDisp <- ggplot(dfScatterQvalvsDispModelByModel, aes(x=sumdroprate, y=qval, group=model)) +
    geom_point(aes(colour = group)) + 
    labs(title="Q-value versus sum of true dropout rates by underlying model type") +
    xlab(paste0("sum of true dropout rates")) +
    ylab(paste0("log10 Q-value")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_Scatter_QvalvsPiModelByModel.pdf",width=7,height=7)
  print(gScatterQvalvsDisp)
  dev.off()
    
  return(NULL)
}