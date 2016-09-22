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

validateOuputSimulation <- function(
  dirLineagePulseTempFiles,
  dirValidationOut,
  matMuHidden,
  matDispHidden,
  matDropoutHidden,
  matDropoutLinModelHidden,
  vecSizeFactorsHidden,
  vecConstIDs,
  vecImpulseIDs,
  scaWindowRadius=NULL ){
  
  dirCurrent <- getwd()
  
  # TODO: add correlation R2 to 1. and null pval distri to 4.
  # Find way to compare logistic noise models
  
  # Load data
  print(paste0("Load data from LineagePulse output directory: ",dirLineagePulseTempFiles))
  setwd(dirLineagePulseTempFiles)
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
  scaNumGenes <- dim(matCountsProc)[1]
  vecConstIDs <- vecConstIDs[vecConstIDs %in% rownames(matCountsProc)]
  vecImpulseIDs <- vecImpulseIDs[vecImpulseIDs %in% rownames(matCountsProc)]
  scaNumGenesConst <- length(vecConstIDs)
  scaNumGenesImpulse <- length(vecImpulseIDs)
  scaNumCells <- dim(matCountsProc)[2]
  
  vecHiddenModelType <- rep("const", length(vecAnalysedGenes))
  names(vecHiddenModelType) <- vecAnalysedGenes
  vecHiddenModelType[vecImpulseIDs] <- "impulse"
  
  # 1.) Scatter plots of parameter sets
  print("# 1. Scatter plots of inferred versus true parameter sets")
  # a) Mean parameters H0
  print("# a) Mean parameters H0")
  dfScatterMuModelvsMuInferredH0 <- data.frame(
    x=log(matMuH0[,1])/log(10),
    y=log(apply(matMuHidden,1, function(gene) mean(gene, na.rm=TRUE)) )/log(10),
    model=vecHiddenModelType )
  gScatterDispvsMean <- ggplot(dfScatterMuModelvsMuInferredH0, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Inferred mean (H0) versus\n underlying model mean parameter.") +
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
  print("# b) Mean parameters H1: For one cell only")
  dfScatterMuModelvsMuInferredH1 <- data.frame(
    x=log(matMuH1[,1])/log(10),
    y=log(matMuHidden[,1])/log(10),
    model=vecHiddenModelType )
  gScatterDispvsMean <- ggplot(dfScatterMuModelvsMuInferredH1, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Inferred mean (H1) of first cell versus\n underlying model mean parameter.") +
    xlab(paste0("log10 H1 inferred mean parameter")) +
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
  print("# c) Dispersion parameters H0 (using dispersion of first cell from model)")
  dfScatterDispModelvsDispInferredH0 <- data.frame(
    x=matDispersionsH0[,1],
    y=matDispHidden[,1],
    model=vecHiddenModelType )
  gScatterDispvsMean <- ggplot(dfScatterDispModelvsDispInferredH0, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Inferred dispersion (H0) versus\n underlying model dispersion parameter.") +
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
  
  # d) Dispersion parameters H1: For one cell only
  print("# d) Dispersion parameters H1: For one cell only")
  dfScatterDispModelvsDispInferredH1 <- data.frame(
    x=matDispersionsH1[,1],
    y=matDispHidden[,1],
    model=vecHiddenModelType )
  gScatterDispvsMean <- ggplot(dfScatterDispModelvsDispInferredH1, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Inferred dispersion (H0) versus\n underlying model dispersion parameter.") +
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
  print("# 3. LRT histogram")
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
    ylab("Frequency") +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  graphics.off()
  pdf("LineagePulse_AnalyseSimulated_Hist_LRT_InferredH1vsHidden_NB.pdf")
  print(gHistLRT_NB)
  dev.off()
  gHistLRT_ZINB <- ggplot( dfLRT_ZINB, aes(value)) +
    geom_histogram(alpha = 1, position = 'identity') +
    ggtitle(paste0("log LRT value under ZINB model")) +
    xlab("Difference in loglikelihood") +
    ylab("Frequency") +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  graphics.off()
  pdf("LineagePulse_AnalyseSimulated_Hist_LRT_InferredH1vsHidden_ZINB.pdf")
  print(gHistLRT_ZINB)
  dev.off() 
  
  # 4. ECDF q-values divided by const and impulse true model
  print("# 4. ECDF q-values divided by const and impulse true model")
  vecX <- seq(max(-100,log(min(dfDEAnalysis$adj.p))/log(10)),0,by=0.5)
  vecCDFConst <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecConstIDs,]$p)))/log(10) <= thres, na.rm=TRUE)})
  vecCDFImpulse <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecImpulseIDs,]$p)))/log(10) <= thres, na.rm=TRUE)})
  dfECDFPvalByModel <- data.frame( thres=vecX,
    pval=c(vecCDFConst,vecCDFImpulse),
    model=c(rep("const",length(vecCDFConst)), rep("impulse",length(vecCDFImpulse))) )
  gECDFPvalByModel <- ggplot(dfECDFPvalByModel, aes(x=thres, y=pval, group=model)) +
    geom_line(aes(colour = model)) +
    labs(title="ECDF log10 p-values by underlying model type") +
    xlab(paste0("log10 p-value threshold")) +
    ylab(paste0("empricial cumulative density")) + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AnalyseSimulated_ECDF_PvalByModelClass.pdf",width=7,height=7)
  print(gECDFPvalByModel)
  dev.off() 
  
  # 5. Q-value as function of model parameters
  print("# 5. Q-value as function of model parameters")
  # a) Q-value as function of average true mean parameter
  print("# a) Q-value as function of average true mean parameter")
  # divided by const and impulse true model
  dfScatterQvalvsMuModelByModel <- data.frame(
    mu=apply(matMuHidden, 1, function(gene) mean(gene, na.rm=TRUE)),
    qval=log(dfDEAnalysis[vecAnalysedGenes,]$adj.p)/log(10),
    model=vecHiddenModelType )
  gScatterQvalvsMean <- ggplot(dfScatterQvalvsMuModelByModel, aes(x=mu, y=qval, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Q-value versus true mean\n by underlying model type") +
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
  print("# b) Q-value as function of true dispersion parameter")
  # divided by const and impulse true model
  dfScatterQvalvsDispModelByModel <- data.frame(
    disp=apply(matDispHidden, 1, function(gene) mean(gene, na.rm=TRUE)),
    qval=log(dfDEAnalysis[vecAnalysedGenes,]$adj.p)/log(10),
    model=vecHiddenModelType )
  gScatterQvalvsDisp <- ggplot(dfScatterQvalvsDispModelByModel, aes(x=disp, y=qval, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Q-value versus true dispersion\n by underlying model type") +
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
  print("# c) Q-value as function of number of drop-outs")
  # divided by const and impulse true model
  dfScatterQvalvsDispModelByModel <- data.frame(
    sumdroprate=apply(matDropoutHidden, 1, function(gene) sum(gene, na.rm=TRUE)),
    qval=log(dfDEAnalysis[vecAnalysedGenes,]$adj.p)/log(10),
    model=vecHiddenModelType )
  gScatterQvalvsDisp <- ggplot(dfScatterQvalvsDispModelByModel, aes(x=sumdroprate, y=qval, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Q-value versus sum of true dropout rates\n by underlying model type") +
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
    
  print("Done generating validation plots with models underlying simulation as reference.")
  # Switch back to original directory
  setwd(dirCurrent)
}