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
#' @seealso Auxillary method not called by LineagePulse wrapper.
#' Called separately by user.
#' 
#' @param matQval: matrix with p-values by methods
#' 
#' @return NULL
#' 
#' @export

validateOuputSimulation <- function(
  dirOutLineagePulse,
  dirOutSimulation,
  dirOutValidation,
  scaWindowRadius=NULL ){
  
  dirCurrent <- getwd()
  
  # Load data
  print(paste0("Load data from Simulation output directory: ",dirOutSimulation))
  setwd(dirOutSimulation)
  load(file=file.path(dirOutSimulation,"Simulation_vecPT.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_vecConstIDs.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_vecImpulseIDs.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matImpulseModelHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matMuHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_vecSizeFactorsHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matDispHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matSampledDataHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matDropoutLinModelHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matDropoutRatesHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matDropoutsHidden.RData"))
  load(file=file.path(dirOutSimulation,"Simulation_matSampledCountsObserved.RData"))

  print(paste0("Load data from LineagePulse output directory: ",dirOutLineagePulse))
  setwd(dirOutLineagePulse)
  # Load data
  load("LineagePulse_matCountsProc.RData")
  load("LineagePulse_vecPseudotimeProc.RData")
  load("LineagePulse_vecSizeFactors.RData")
  load("LineagePulse_lsMuModelH0.RData")
  load("LineagePulse_lsDispModelH0.RData")
  load("LineagePulse_lsMuModelH1.RData")
  load("LineagePulse_lsDispModelH1.RData")
  load("LineagePulse_lsDropModel.RData")
  #load("LineagePulse_matZH1.RData")
  load("LineagePulse_dfDEAnalysis.RData")
  
  vecAnalysedGenes <- rownames(matCountsProc)
  vecAnalysedCells <- colnames(matCountsProc)
  matMuHidden <- matMuHidden[vecAnalysedGenes,vecAnalysedCells]
  matDispHidden <- matDispHidden[vecAnalysedGenes,vecAnalysedCells]
  matDropoutRatesHidden <- matDropoutRatesHidden[vecAnalysedGenes,vecAnalysedCells]
  matDropoutLinModelHidden <- matDropoutLinModelHidden[vecAnalysedCells,]
  vecSizeFactorsHidden <- vecSizeFactorsHidden[vecAnalysedCells]
  
  # Initialise
  setwd(dirOutValidation)
  scaNumGenes <- dim(matCountsProc)[1]
  vecConstIDs <- vecConstIDs[vecConstIDs %in% rownames(matCountsProc)]
  vecImpulseIDs <- vecImpulseIDs[vecImpulseIDs %in% rownames(matCountsProc)]
  scaNumGenesConst <- length(vecConstIDs)
  scaNumGenesImpulse <- length(vecImpulseIDs)
  scaNumCells <- dim(matCountsProc)[2]
  
  scaEps <- 10^(-5) # Constant to be added to paramters before taking log 
  # to avoid errors at zero values.
  
  vecHiddenModelType <- rep("const", length(vecAnalysedGenes))
  names(vecHiddenModelType) <- vecAnalysedGenes
  vecHiddenModelType[vecImpulseIDs] <- "impulse"
  
  # 1.) Scatter plots of parameter sets
  print("# 1. Scatter plots of inferred versus true parameter sets")
  # a) Mean parameters H0
  print("# a) Mean parameters H0")
  # Generate data frame for ggplot
  dfScatterMuModelvsMuInferredH0 <- data.frame(
    x=log(lsMuModelH0$matMuModel)/log(10),
    y=log(apply(matMuHidden,1, function(gene) mean(gene, na.rm=TRUE)) )/log(10),
    model=vecHiddenModelType )
  # Get correlation
  scaCorrMuModelvsMuInferredH0 <- round(cor(dfScatterMuModelvsMuInferredH0$x,
    dfScatterMuModelvsMuInferredH0$y), 2)
  # Generate plots
  gScatterMuModelvsMuInferredH0 <- ggplot(dfScatterMuModelvsMuInferredH0, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Inferred mean parameter (H0) versus\n underlying model mean parameter.") +
    xlab(paste0("log10 H0 inferred mean parameter")) +
    ylab(paste0("log10 gene-wise average of model mean parameter")) + 
    scale_fill_continuous(name = "Count") + 
    geom_text(x=4/5*(max(dfScatterMuModelvsMuInferredH0$x)-min(dfScatterMuModelvsMuInferredH0$x)),
      y=min(dfScatterMuModelvsMuInferredH0$x),
      label=paste0("R^2=",scaCorrMuModelvsMuInferredH0),
      size=5 ) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  # Plot
  pdf("LineagePulseSim_Scatter_MuModelvsMuInferredH0.pdf",width=7,height=7)
  print(gScatterMuModelvsMuInferredH0)
  dev.off()
  graphics.off()
  
  # b) Mean parameters H1: For one cell only
  print("# b) Mean parameters H1: Average over all cells")
  vecMuParamAveH1 <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
    mean( decompressMeansByGene(vecMuModel=lsMuModelH1$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL) )
  }))
  # Generate data frame for ggplot
  dfScatterMuModelvsMuInferredH1 <- data.frame(
    x=log(vecMuParamAveH1)/log(10),
    y=log(apply(matMuHidden, 1, function(i) mean(i, na.rm=TRUE)))/log(10),
    model=vecHiddenModelType )
  # Get correlation
  scaCorrMuModelvsMuInferredH1 <- round(cor(dfScatterMuModelvsMuInferredH1$x,
    dfScatterMuModelvsMuInferredH1$y), 2)
  # Generate plots
  gScatterMuModelvsMuInferredH1 <- ggplot(dfScatterMuModelvsMuInferredH1, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title="Average inferred mean parameter (H1) versus\n average underlying model mean parameter.") +
    xlab(paste0("log10 H1 inferred mean parameter")) +
    ylab(paste0("log10 model mean parameter")) + 
    scale_fill_continuous(name = "Count") + 
    geom_text(x=4/5*(max(dfScatterMuModelvsMuInferredH1$x)-min(dfScatterMuModelvsMuInferredH1$x)),
      y=min(dfScatterMuModelvsMuInferredH1$x),
      label=paste0("R^2=",scaCorrMuModelvsMuInferredH1),
      size=5 ) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  # Plot
  pdf("LineagePulseSim_Scatter_MuModelvsMuInferredH1.pdf",width=7,height=7)
  print(gScatterMuModelvsMuInferredH1)
  dev.off()
  graphics.off()
  
  # c) Dispersion parameters H0 (using dispersion of first cell from model)
  print("# c) Dispersion parameters H0")
  # Generate data frame for ggplot
  dfScatterDispModelvsDispInferredH0 <- data.frame(
    x=lsDispModelH0$matDispModel,
    y=apply(matDispHidden, 1, function(i) median(i, na.rm=TRUE)),
    model=vecHiddenModelType )
  # Outlier detection
  boolOutlier <- dfScatterDispModelvsDispInferredH0$x < 0.01 |
    dfScatterDispModelvsDispInferredH0$x > 5
  dfScatterDispModelvsDispInferredH0 <- dfScatterDispModelvsDispInferredH0[!boolOutlier,]
  # Get correlation
  scaCorrDispModelvsDispInferredH0 <- round(cor(dfScatterDispModelvsDispInferredH0$x,
    dfScatterDispModelvsDispInferredH0$y), 2)
  # Generate plots
  gDispModelvsDispInferredH0 <- ggplot(
    dfScatterDispModelvsDispInferredH0, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title=paste0("Inferred dispersion parameter (H0) versus\n",
      " median underlying model dispersion parameter")) +
    xlab(paste0("inferred dispersion parameter (H0)")) +
    ylab(paste0("median underlying model dispersion parameter")) + 
    scale_fill_continuous(name = "Count") + 
    geom_text(x=4/5*(max(dfScatterDispModelvsDispInferredH0$x)-min(dfScatterDispModelvsDispInferredH0$x)),
      y=min(dfScatterDispModelvsDispInferredH0$x),
      label=paste0("R^2=",scaCorrDispModelvsDispInferredH0),
      size=5 ) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  # Plot
  pdf("LineagePulseSim_Scatter_DispModelMedianvsDispInferredH0.pdf",width=7,height=7)
  print(gDispModelvsDispInferredH0)
  dev.off()
  
  # d) Dispersion parameters H1: Median
  print("# d) Dispersion parameters H1: Median")
  vecDispParamH1 <- sapply(seq(1,scaNumGenes), function(i){
    vecDispParamH1Gene_i <- decompressDispByGene(vecDispModel=lsDispModelH1$matDispModel[i,],
      lsDispModelGlobal=lsDispModelH1$lsDispModelGlobal,
      vecInterval=NULL)
    return(median(vecDispParamH1Gene_i, na.rm=TRUE))
  })
  # Generate data frame for ggplot
  dfScatterDispModelvsDispInferredMedianH1 <- data.frame(
    x=vecDispParamH1,
    y=matDispHidden[,1],
    model=vecHiddenModelType )
  # Outlier detection
  boolOutlier <- dfScatterDispModelvsDispInferredMedianH1$x < 0.01 |
    dfScatterDispModelvsDispInferredMedianH1$x > 5
  dfScatterDispModelvsDispInferredMedianH1 <- dfScatterDispModelvsDispInferredMedianH1[!boolOutlier,]
  # Get correlation
  scaCorrDispModelvsDispInferredMedianH1 <- round(cor(dfScatterDispModelvsDispInferredMedianH1$x,
    dfScatterDispModelvsDispInferredMedianH1$y), 2)
  # Generate plots
  gScatterDispModelvsDispInferredMedianH1 <- ggplot(
    dfScatterDispModelvsDispInferredMedianH1, aes(x=x, y=y, group=model)) +
    geom_point(aes(colour = model)) + 
    labs(title=paste0("median inferred dispersion parameter (H1) versus\n",
      " median underlying model dispersion parameter.")) +
    xlab(paste0("median inferred dispersion parameter (H1)")) +
    ylab(paste0("median model dispersion parameter")) + 
    scale_fill_continuous(name = "Count") +
    geom_text(x=4/5*(max(dfScatterDispModelvsDispInferredMedianH1$x)-min(dfScatterDispModelvsDispInferredMedianH1$x)),
      y=min(dfScatterDispModelvsDispInferredMedianH1$x),
      label=paste0("R^2=",scaCorrDispModelvsDispInferredMedianH1),
      size=5 ) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  # Plot
  pdf("LineagePulseSim_Scatter_DispModelMedianvsDispInferredMedianH1.pdf",width=7,height=7)
  print(gScatterDispModelvsDispInferredMedianH1)
  dev.off()
  
  # 2. Dropout model
  NCellsToPlot <- 100
  NCellsToPlot <- min(scaNumCells,NCellsToPlot)
  print("# 2. Plot true and inferred dropout models by cell")
  # PDF with one page per cell: Plot data and inferred and true logistic model
  # Reshape matrices into single column first:
  lsGplotsDropoutModelsByCell <- list()
  for(cell in seq(1,NCellsToPlot)){
    lsMuHiddenSort <- sort(log(matMuHidden[,cell]+scaEps)/log(10), index.return=TRUE)
    vecMuHiddenSort <- lsMuHiddenSort$x
    vecDropoutHiddenSort <- (matDropoutRatesHidden[,cell])[lsMuHiddenSort$ix]
    vecMuParamH1 <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
      decompressMeansByGene(vecMuModel=lsMuModelH1$matMuModel[i,],
        lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
        vecInterval=cell)
    }))
    vecDropParamH1 <- decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[cell,],
      vecMu=vecMuParamH1,
      matPiConstPredictors=lsDropModel$matPiConstPredictors )
    lsMuH1Sort <- sort(log(vecMuParamH1+scaEps)/log(10), index.return=TRUE)
    vecMuH1Sort <- lsMuH1Sort$x
    vecDropParamH1Sort <- vecDropParamH1[lsMuH1Sort$ix] 
    dfDropoutModelsByCell <- data.frame(
      mu=c(vecMuHiddenSort, vecMuH1Sort),
      dropout=c(vecDropoutHiddenSort, vecDropParamH1Sort),
      model=c(rep("true", scaNumGenes), rep("H1", scaNumGenes)) )
    lsGplotsDropoutModelsByCell[[cell]] <- ggplot(dfDropoutModelsByCell) +
      geom_line(aes(x=mu, y=dropout, colour=model)) +
      labs(title=paste0("Cell ", cell, ": Inferred and true logistic drop-out model")) +
      xlab(paste0("log_10 true mean parameter")) +
      ylab(paste0("dropout rate")) +
      geom_density(aes(mu, colour=model)) +          
      theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))
  }
  pdf("LineagePulseSim_DropoutModelByCell.pdf",width=7,height=7)
  for(gplot in lsGplotsDropoutModelsByCell){
    print(gplot)
  }
  dev.off()
  graphics.off()
  
  # 3. Mean model by gene
  scaNGenesToPlot <- 50
  print("# 3. Plot inferred pseudotime expression models by gene.")
  # Plot data, true and inferred model by gene.
  if(!is.null(vecConstIDs) & !is.null(vecImpulseIDs)){
    # Plot constant genes
    print("# a) Plot constant genes.")
    vecIDsToPlot <- vecConstIDs[seq(1,min(scaNGenesToPlot,length(vecConstIDs)))]
    lsGplotsConstantIDs <- list()
    for(id in vecIDsToPlot){
      scaConstModelRefParam <- matMuHidden[id,1]
      vecImpulseModelRefParam <- NULL
      
      lsGplotsConstantIDs[[match(id, vecIDsToPlot)]] <- plotGene(vecCounts=matCountsProc[id,],
        vecPseudotime=vecPseudotimeProc,
        vecDropoutRates=matDropoutRatesHidden[id,],
        vecImpulseModelParam=lsMuModelH1$matMuModel[id,],
        vecImpulseModelRefParam=vecImpulseModelRefParam,
        scaConstModelParam=lsMuModelH0$matMuModel[id,],
        scaConstModelRefParam=scaConstModelRefParam,
        strGeneID=id,
        strTitleSuffix=paste0("Q-value ", round(dfDEAnalysis[id,"adj.p"],3)))
    }
    pdf("LineagePulseSim_ExpressionTracesConstant.pdf")
    for(gplot in lsGplotsConstantIDs){
      print(gplot)
    }
    dev.off()
    graphics.off()
    
    # Plot impulse genes
    print("# b) Plot impulse genes.")
    vecIDsToPlot <- vecImpulseIDs[seq(1,min(scaNGenesToPlot,length(vecImpulseIDs)))]
    lsGplotsImpulseIDs <- list()
    for(id in vecIDsToPlot){
      scaConstModelRefParam <- NULL
      vecImpulseModelRefParam <- matImpulseModelHidden[id,]
      vecImpulseModelRefParam[2] <- log(vecImpulseModelRefParam[2])
      vecImpulseModelRefParam[3] <- log(vecImpulseModelRefParam[3])
      vecImpulseModelRefParam[4] <- log(vecImpulseModelRefParam[4])
      
      lsGplotsImpulseIDs[[match(id, vecIDsToPlot)]] <- plotGene(vecCounts=matCountsProc[id,],
        vecPseudotime=vecPseudotimeProc,
        vecDropoutRates=matDropoutRatesHidden[id,],
        vecImpulseModelParam=lsMuModelH1$matMuModel[id,],
        vecImpulseModelRefParam=vecImpulseModelRefParam,
        scaConstModelParam=lsMuModelH0$matMuModel[id,],
        scaConstModelRefParam=scaConstModelRefParam,
        strGeneID=id,
        strTitleSuffix=paste0("Q-value ", round(dfDEAnalysis[id,"adj.p"],3)))
    }
    pdf("LineagePulseSim_ExpressionTracesImpulse.pdf")
    for(gplot in lsGplotsImpulseIDs){
      print(gplot)
    }
    dev.off()
    graphics.off()
    
  } else {
    vecIDsToPlot <- rownames(matCountsProc)[seq(1,min(scaNGenesToPlot,scaNumGenes))]
    lsGplotsGeneIDs <- list()
    for(id in vecIDsToPlot){
      if(id %in% vecConstIDs){
        scaConstModelRefParam <- matMuHidden[id,1]
      } else {
        scaConstModelRefParam <- NULL
      }
      if(id %in% vecImpulseIDs){
        vecImpulseModelRefParam <- matImpulseModelHidden[id,]
        vecImpulseModelRefParam[2] <- log(vecImpulseModelRefParam[2])
        vecImpulseModelRefParam[3] <- log(vecImpulseModelRefParam[3])
        vecImpulseModelRefParam[4] <- log(vecImpulseModelRefParam[4])
      } else {
        vecImpulseModelRefParam <- NULL
      }
      lsGplotsGeneIDs[[match(id, vecIDsToPlot)]] <- plotGene(vecCounts=matCountsProc[id,],
        vecPseudotime=vecPseudotimeProc,
        vecDropoutRates=matDropoutRatesHidden[id,],
        vecImpulseModelParam=lsMuModelH1$matMuModel[id,],
        vecImpulseModelRefParam=vecImpulseModelRefParam,
        scaConstModelParam=lsMuModelH0$matMuModel[id,],
        scaConstModelRefParam=scaConstModelRefParam,
        strGeneID=id,
        strTitleSuffix=paste0("Q-value ", round(dfDEAnalysis[id,"adj.p"],3)))
    }
    pdf("LineagePulseSim_ExpressionTraces.pdf")
    for(gplot in lsGplotsGeneIDs){
      print(gplot)
    }
    dev.off()
    graphics.off()
  }
  
  # 4. LRT hist
  print("# 4. LRT histogram")
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
    evalLogLikGene(vecCounts=matCountsProc[i,],
      vecMu=matMuHiddenTemp[i,],
      vecSizeFactors=vecSizeFactorsHidden,
      vecDispEst=matDispHidden[i,], 
      vecDropoutRateEst=matDropoutRatesHidden[i,],
      vecboolNotZeroObserved=matCountsProc[i,]>0 & !is.na(matCountsProc[i,]), 
      vecboolZero=matCountsProc[i,]==0,
      scaWindowRadius=scaWindowRadius)
  }))
  # Compute loglikelihood of inferred model
  vecNBLLH1Model <- sapply(seq(1,scaNumGenes), function(i){
    vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH1 <- decompressDispByGene( vecDispModel=lsDispModelH1$matDispModel[i,],
      lsDispModel=lsDispModelH1$lsDispModelGlobal,
      vecInterval=NULL )
    scaLL <- sum(dnbinom(x=matCountsProc[i,],
      mu=vecMuParamH1,
      size=vecDispParamH1,
      log=TRUE), na.rm=TRUE)
    return(scaLL)
  })
  vecZINBLLH1Model <- unlist(lapply( seq(1,scaNumGenes), function(i){
    vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH1 <- decompressDispByGene( vecDispModel=lsDispModelH1$matDispModel[i,],
      lsDispModel=lsDispModelH1$lsDispModelGlobal,
      vecInterval=NULL )
    vecDropoutParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH1,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
    scaLL <- evalLogLikGene(vecCounts=matCountsProc[i,],
      vecMu=vecMuParamH1,
      vecSizeFactors=vecSizeFactors,
      vecDispEst=vecDispParamH1, 
      vecDropoutRateEst=vecDropoutParamH1,
      vecboolNotZeroObserved=matCountsProc[i,]>0 & !is.na(matCountsProc[i,]), 
      vecboolZero=matCountsProc[i,]==0,
      scaWindowRadius=scaWindowRadius)
    return(scaLL)
  }))
  # LRT
  scaQThres <- 10^(-3)
  vecLRT_NB <- vecNBLLH1Model-vecNBLLHiddenModel
  vecLRT_ZINB <- vecZINBLLH1Model-vecZINBLLHiddenModel
  dfLRT_NB <- data.frame(value=vecLRT_NB)
  dfLRT_ZINB <- data.frame(value=vecLRT_ZINB)
  gHistLRT_NB <- ggplot( dfLRT_NB, aes(value)) +
    geom_histogram(alpha = 1, position = 'identity') +
    ggtitle(paste0("Log likelihoodratio under NB model\nLL(Inferred model H1) - LL(underlying model)")) +
    xlab("Difference in loglikelihood") +
    ylab("Frequency") +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  graphics.off()
  pdf("LineagePulseSim_Hist_LRT_InferredH1vsHidden_NB.pdf")
  print(gHistLRT_NB)
  dev.off()
  gHistLRT_ZINB <- ggplot( dfLRT_ZINB, aes(value)) +
    geom_histogram(alpha = 1, position = 'identity') +
    ggtitle(paste0("Log likelihoodratio under ZINB model\nLL(Inferred model H1) - LL(underlying model)")) +
    xlab("Difference in loglikelihood") +
    ylab("Frequency") +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  graphics.off()
  pdf("LineagePulseSim_Hist_LRT_InferredH1vsHidden_ZINB.pdf")
  print(gHistLRT_ZINB)
  dev.off() 
  
  # 5. ECDF q-values divided by const and impulse true model
  print("# 5. ECDF q-values divided by const and impulse true model")
  vecX <- seq(max(-100,round(log(min(dfDEAnalysis$adj.p))/log(10))),0,by=0.5)
  # Compute observed ECDF for constant and for impulse distributed genes
  vecCDFConst <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecConstIDs,]$p)))/log(10) <= thres, na.rm=TRUE)})
  vecCDFImpulse <- sapply(vecX, function(thres){
    sum(log(as.numeric(as.vector(dfDEAnalysis[vecImpulseIDs,]$p)))/log(10) <= thres, na.rm=TRUE)})
  # Compute CDF under null p-value distribution: uniform
  # Scale by number of observations for constant model, against
  # which the null is compared
  vecCDFBackgroundNull <- 10^(vecX)*scaNumGenesConst
  dfECDFPvalByModel <- data.frame( thres=vecX,
    pval=c(vecCDFConst,vecCDFImpulse,vecCDFBackgroundNull),
    model=c(rep("const",length(vecCDFConst)), 
      rep("impulse",length(vecCDFImpulse)),
      rep("null_distribution",length(vecCDFBackgroundNull))) )
  gECDFPvalByModel <- ggplot(dfECDFPvalByModel, aes(x=thres, y=pval, group=model)) +
    geom_line(aes(colour = model)) +
    labs(title="ECDF log10 p-values by underlying model type") +
    xlab(paste0("log10 p-value threshold")) +
    ylab(paste0("empricial cumulative density")) + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulseSim_ECDF_PvalByModelClass.pdf",width=7,height=7)
  print(gECDFPvalByModel)
  dev.off() 
  
  # 6. Q-value as function of model parameters
  print("# 6. Q-value as function of model parameters")
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
  pdf("LineagePulseSim_Scatter_QvalvsMuModelByModel.pdf",width=7,height=7)
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
  pdf("LineagePulseSim_Scatter_QvalvsDispModelByModel.pdf",width=7,height=7)
  print(gScatterQvalvsDisp)
  dev.off()
  
  # c) Q-value as function of number of drop-outs
  print("# c) Q-value as function of number of drop-outs")
  # divided by const and impulse true model
  dfScatterQvalvsDispModelByModel <- data.frame(
    sumdroprate=apply(matDropoutRatesHidden, 1, function(gene) sum(gene, na.rm=TRUE)),
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
  pdf("LineagePulseSim_Scatter_QvalvsPiModelByModel.pdf",width=7,height=7)
  print(gScatterQvalvsDisp)
  dev.off()
    
  print("Done generating validation plots with models underlying simulation as reference.")
  # Switch back to original directory
  setwd(dirCurrent)
}