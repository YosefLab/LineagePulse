#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++  Validation metrics for ZINB fits  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Generate metrics for validation of ZINB fit
#' 
#' Compute and plot summary metrics that allow the user to judge
#' whether the ZINB fit is sensible. All metrics are saved as plots
#' to pdfs.
#' 
#' @seealso Called by \code{runLineagePulse} or separately by user.
#' 
#' @param dirOutLineagePulse: (str directory)
#'    Directory to which detailed LineagePulse output is saved to.
#'    Model matrices are loaded from here.
#' @param dirOutValidation: (str directory)
#'    Directory to which validation metrics/plots are saved to.
#'    This is the output directory for this function.
#' 
#' @return NULL: Output saved to dirOutValidation.
#'    
#' @author David Sebastian Fischer
#' 
#' @export

validateOutput <- function(dirOutLineagePulse,
  dirOutValidation ){
  
  # Parameters
  # Point 2:
  scaNIDs <- 100  # Genes/Cells to plot for each group
  scaThres <- 10^(-3) # Expression threshold for highly expressed group
  scaFracExprIsHigh <- 0.5 # Fraction of cells above 
  #expression threshold to be in the highly expressed group
  
  scaEps <- 10^(-5) # Constant to be added to paramters before taking log 
  # to avoid errors at zero values.
  
  # Get current wd to change back to at the end
  dirCurrent <- getwd()
  
  # Load data
  print(paste0("Load data from LP output directory: ",dirOutLineagePulse))
  setwd(dirOutLineagePulse)
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsInputParam.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_matCountsProc.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_vecPseudotimeProc.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsMuModelH0.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsDispModelH0.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsMuModelH1.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsDispModelH1.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_lsDropModel.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_dfDEAnalysis.RData"))
  load(file=file.path(dirOutLineagePulse,"LineagePulse_vecNormConst.RData"))
  
  # Initialise
  setwd(dirOutValidation)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  scaWindowRadius <- lsInputParam$scaWindowRadius
  scaNIDs <- min(scaNIDs, scaNumGenes)
  
  # 1. ECDF q-values
  print("# 1. Print ECDF of q-values.")
  vecX <- seq(max(-30,min(log(dfDEAnalysis$adj.p)/log(10))),0,by=0.5)
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
  
  # 2. Plot Fits
  print("# 2. Print model fits sorted into groups.")
  # a) Plot highly expressed genes.
  print("# a) Plot top significant of highly expressed genes.")
  vecboolIDsHighExpr <- apply(matCountsProc, 1, function(gene){
    sum(gene > 50) >= scaFracExprIsHigh*scaNumCells})
  vecIDsHighExpr <- rownames(matCountsProc)[vecboolIDsHighExpr]
  vecIDsTopQvalHighExpre <- as.vector(dfDEAnalysis[dfDEAnalysis$Gene %in% vecIDsHighExpr,]$Gene)
  vecIDsTopQvalHighExpre <- (vecIDsTopQvalHighExpre[!is.na(vecIDsTopQvalHighExpre)])[
    1:min(scaNIDs,sum(!is.na(vecIDsTopQvalHighExpre)))]
  if(!is.na(vecIDsTopQvalHighExpre)){
    lsGplotsHighExpr <- list()
    for(id in vecIDsTopQvalHighExpre){
      # Decompress parameters by gene
      vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[id,],
        lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
        vecInterval=NULL )
      vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
        vecMu=vecMuParamH1,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[id,] )
      if(lsMuModelH1$lsMuModelGlobal$strMuModel=="impulse"){
        vecImpulseModelParam=lsMuModelH1$matMuModel[id,]
      } else {
        vecImpulseModelParam <- NULL
      }
      if(lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
        scaConstModelParam=lsMuModelH0$matMuModel[id,]
      } else {
        scaConstModelParam <- NULL
      }
      
      lsGplotsHighExpr[[match(id, vecIDsTopQvalHighExpre)]] <- plotGene(vecCounts=matCountsProc[id,],
        vecPseudotime=vecPseudotimeProc,
        vecPiParamH1=vecPiParamH1,
        vecImpulseModelParam=vecImpulseModelParam,
        scaConstModelParam=scaConstModelParam,
        strGeneID=id,
        strTitleSuffix=paste0("Q-value ", dfDEAnalysis[id,"adj.p"]))
    }
  } else {
    lsGplotsHighExpr <- NULL
  }
  pdf("LineagePulse_ImpulseTraces_HighlyExpressedLowQval.pdf")
  for(gplot in lsGplotsHighExpr){
    print(gplot)
  }
  dev.off()
  graphics.off()
  
  # b) Plot top q-val
  print("# b) Plot top significant genes.")
  vecIDsTopQval <- as.vector(dfDEAnalysis[1:min(scaNIDs,dim(dfDEAnalysis)[1]),]$Gene)
  lsGplotsTopQval <- list()
  for(id in vecIDsTopQval){
    # Decompress parameters by gene
    vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[id,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH1,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[id,] )
    if(lsMuModelH1$lsMuModelGlobal$strMuModel=="impulse"){
      vecImpulseModelParam=lsMuModelH1$matMuModel[id,]
    } else {
      vecImpulseModelParam <- NULL
    }
    if(lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
      scaConstModelParam=lsMuModelH0$matMuModel[id,]
    } else {
      scaConstModelParam <- NULL
    }
    
    lsGplotsTopQval[[match(id, vecIDsTopQval)]] <- plotGene(vecCounts=matCountsProc[id,],
      vecPseudotime=vecPseudotimeProc,
      vecPiParamH1=vecPiParamH1,
      vecImpulseModelParam=vecImpulseModelParam,
      scaConstModelParam=scaConstModelParam,
      strGeneID=id,
      strTitleSuffix=paste0("Q-value ", dfDEAnalysis[id,"adj.p"]))
  }
  pdf("LineagePulse_ImpulseTraces_LowQval.pdf")
  for(gplot in lsGplotsTopQval){
    print(gplot)
  }
  dev.off()
  graphics.off()
  
  # c) Plot top q-val with point size adjusted for LL ratio of models
  print("# c) Plot top significant genes with point size adjusted for LL ratio of models.")
  vecIDsTopQval <- as.vector(dfDEAnalysis[1:min(scaNIDs,dim(dfDEAnalysis)[1]),]$Gene)
  lsGplotsTopQval <- list()
  for(id in vecIDsTopQval){
    # Decompress parameters by gene
    vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[id,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH1 <- decompressDispByGene( vecDispModel=lsDispModelH1$matDispModel[id,],
      lsDispModelGlobal=lsDispModelH1$lsDispModelGlobal,
      vecInterval=NULL )
    vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH1,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[id,] )
    
    vecMuParamH0 <- decompressMeansByGene( vecMuModel=lsMuModelH0$matMuModel[id,],
      lsMuModelGlobal=lsMuModelH0$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH0 <- decompressDispByGene( vecDispModel=lsDispModelH0$matDispModel[id,],
      lsDispModelGlobal=lsDispModelH0$lsDispModelGlobal,
      vecInterval=NULL )
    vecPiParamH0 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH0,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[id,] )
    
    if(lsMuModelH1$lsMuModelGlobal$strMuModel=="impulse"){
      vecImpulseModelParam=lsMuModelH1$matMuModel[id,]
    } else {
      vecImpulseModelParam <- NULL
    }
    if(lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
      scaConstModelParam=lsMuModelH0$matMuModel[id,]
    } else {
      scaConstModelParam <- NULL
    }
    
    lsGplotsTopQval[[match(id, vecIDsTopQval)]] <- plotGene(vecCounts=matCountsProc[id,],
      vecPseudotime=vecPseudotimeProc,
      vecPiParamH1=vecPiParamH1,
      vecImpulseModelParam=vecImpulseModelParam,
      scaConstModelParam=scaConstModelParam,
      strGeneID=id,
      strTitleSuffix=paste0("Q-value ", dfDEAnalysis[id,"adj.p"]),
      boolScaleByLL=TRUE,
      vecMuParamH0=vecMuParamH0,
      vecMuParamH1=vecMuParamH1,
      vecDispParamH0=vecDispParamH0,
      vecDispParamH1=vecDispParamH1,
      vecPiParamH0=vecPiParamH0,
      vecNormConst=vecNormConst,
      scaWindowRadius=lsInputParam$scaWindowRadius)
  }
  pdf("LineagePulse_ImpulseTraces_LowQval_LLRatioByPoint.pdf")
  for(gplot in lsGplotsTopQval){
    print(gplot)
  }
  dev.off()
  graphics.off()
  
  # d) Plot worst significant q-val
  print("# d) Plot least significant genes.")
  vecboolSigIDs <- which(as.vector(dfDEAnalysis$adj.p) > scaThres)
  if(length(vecboolSigIDs)>0){
    vecIDsWorstQval <- as.vector(dfDEAnalysis[vecboolSigIDs[length(vecboolSigIDs):
        max(1,length(vecboolSigIDs)-scaNIDs)],]$Gene)
    lsGplotsHighQval <- list()
    for(id in vecIDsWorstQval){
      # Decompress parameters by gene
      vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[id,],
        lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
        vecInterval=NULL )
      vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
        vecMu=vecMuParamH1,
        vecPiConstPredictors=lsDropModel$matPiConstPredictors[id,] )
      if(lsMuModelH1$lsMuModelGlobal$strMuModel=="impulse"){
        vecImpulseModelParam=lsMuModelH1$matMuModel[id,]
      } else {
        vecImpulseModelParam <- NULL
      }
      if(lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
        scaConstModelParam=lsMuModelH0$matMuModel[id,]
      } else {
        scaConstModelParam <- NULL
      }
      
      lsGplotsHighQval[[match(id, vecIDsWorstQval)]] <- plotGene(vecCounts=matCountsProc[id,],
        vecPseudotime=vecPseudotimeProc,
        vecPiParamH1=vecPiParamH1,
        vecImpulseModelParam=vecImpulseModelParam,
        scaConstModelParam=scaConstModelParam,
        strGeneID=id,
        strTitleSuffix=paste0("Q-value ", round(dfDEAnalysis[id,"adj.p"],3)))
    }
  } else {
    lsGplotsHighQval <- NULL
  }
  pdf("LineagePulse_ImpulseTraces_HighQval.pdf")
  for(gplot in lsGplotsHighQval){
    print(gplot)
  }
  dev.off()
  graphics.off()
  
  # 3. Look at dropout function
  print("# 3. Dropout function")
  # Sequencing depth as a comparative experimental measure for drop-out
  vecDepth <- apply(matCountsProc, 2, function(cell){sum(cell, na.rm=TRUE)})
  names(vecDepth) <- colnames(matCountsProc)
  
  # a) Inflexion point of logistic
  print("# b) Scatter plot inflexion point logistic drop-out model versus sequencing depth by cell")
  vecInflexLogistic <- apply(lsDropModel$matDropoutLinModel, 1, function(cellmodel){
    return(-cellmodel[1]/cellmodel[2])
  })
  
  dfScatter3 <- data.frame(
    x=vecDepth,
    y=vecInflexLogistic)
  g3 <- ggplot(dfScatter3, aes(x=x, y=y)) +
    geom_point() + 
    labs(title="Scatter plot inflexion point logistic drop-out model\n versus sequencing depth by cell") +
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
  
  # 4. Plot drop-out rate as function of cell
  print("# 4. Plot drop-out rate as function of mean parameter by cell.")
  pdf("LineagePulse_LogisticScatter_DropoutvsMeanbyCell.pdf",width=7,height=7)
  for(cell in seq(1,min(scaNumCells,scaNIDs))){
    # Decompress parameters by gene
    vecMuParamH1 <- do.call(rbind, lapply(seq(1,scaNumGenes), function(i){
      decompressMeansByGene(vecMuModel=lsMuModelH1$matMuModel[i,],
        lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
        vecInterval=cell)
    }))
    vecDropParamH1 <- decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[cell,],
      vecMu=vecMuParamH1,
      matPiConstPredictors=lsDropModel$matPiConstPredictors )
    
    plot( log(vecMuParamH1)/log(10), 
      vecDropParamH1, 
      xlim=c(-4,max(log(vecMuParamH1)/log(10))),
      ylab="dropout rate H1",
      xlab="log_10 mu H1",
      main=paste0("cell_",cell," (Sequencing depth: ",round(vecDepth[cell]/10^5,1),"e5)"))
  }
  dev.off()
  graphics.off()
  
  # 5. Logistic dropout rate fit as function of mean coloured by
  # sequencing depth.
  scaNCells5 <- min(scaNumCells,1000)
  print(paste0("# 5. Plot logistic dropout rate fit as function of mean",
    " coloured by sequencing depth."))
  # Define mu values at which logistic function is sampled
  scaMaxMu <- 10^4
  vecMu=sapply(seq( 1, round(log(scaMaxMu)/log(2)), by=0.2 ), function(x) 2^x)
  matDropoutSampled <- matrix(NA, nrow=length(vecMu), ncol=scaNCells5)
  for(cell in seq(1,scaNCells5)){
    matDropoutSampled[,cell] <- decompressDropoutRateByCell(vecDropModel=lsDropModel$matDropoutLinModel[cell,],
      vecMu=vecMu,
      matPiConstPredictors=lsDropModel$matPiConstPredictors )
  }
  colnames(matDropoutSampled) <- rownames(lsDropModel$matDropoutLinModel[seq(1,scaNCells5),])
  
  dfDropoutSampledMolten <- melt(matDropoutSampled)
  # Reset non-uniue depth values
  vecDepthUnique <- vecDepth
  while(any(duplicated(vecDepthUnique))){
    vecDepthUnique[duplicated(vecDepthUnique)] <- vecDepthUnique[duplicated(vecDepthUnique)]+0.001
  }
  dfLines1 <- data.frame(
    mu=rep(log(vecMu)/log(10), scaNCells5),
    dropout=dfDropoutSampledMolten$value,
    log_depth=log(vecDepthUnique[match(dfDropoutSampledMolten$Var2,names(vecDepthUnique))]))
  
  gLogisticDropoutvsMubyDepth <- ggplot(dfLines1, aes(x=mu, y=dropout, group=log_depth)) +
    geom_line(aes(colour=log_depth)) +
    scale_colour_gradient(low="red",high="green") +
    labs(title="Logistic drop-out model versus sequencing depth by cell") +
    xlab(paste0("log_10 mean parameter H1")) +
    ylab(paste0("dropout rate")) +
    #    xlim(c(-1,2)) + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_AllCurves_LogisticDropoutvsMubyDepth.pdf",width=7,height=7)
  print(gLogisticDropoutvsMubyDepth)
  dev.off()
  graphics.off()
  
  # 6. Plot variance measure as function of mean
  print("# 6. Plot variance measure as function of mean")
  if(lsMuModelH0$lsMuModelGlobal$strMuModel=="constant"){
    vecMu <- lsMuModelH0$matMuModel
    # Transform to vector without col name
    vecMu <- array(vecMu)
  } else {
    print(paste0("H0 is not constant model, set gene-wise mean to ",
      "mean of all mean estimates of a gene"))
    vecMu <- sapply(seq(1,scaNumGenes), function(i){
      vecMuParamH0 <- decompressMeansByGene( vecMuModel=lsMuModelH0$matMuModel[i,],
        lsMuModelGlobal=lsMuModelH0$lsMuModelGlobal,
        vecInterval=NULL )
      return(mean(vecMuParamH0, na.rm=TRUE))
    })
  }
  if(lsDispModelH0$lsDispModelGlobal$strDispModel=="constant"){
    vecDisp <- lsDispModelH0$matDispModel
    # Transform to vector without col name
    vecDisp <- array(vecDisp)
  } else {
    print("WARNING: H0 dispersion model is not constant.",
      " Dispersion parameter in following plots chosen",
      " as median parameter estimated.")
    vecDisp <- sapply(seq(1,scaNumGenes), function(i){
      vecDispParamH0 <- decompressDispByGene( vecDispModel=lsDispModelH0$matDispModel[i,],
        lsDispModelGlobal=lsDispModelH0$lsDispModelGlobal,
        vecInterval=NULL )
      return(median(vecDispParamH0, na.rm=TRUE))
    })
  }
  # Do Outlier detection
  print("Remove dispersion outliers.")
  vecDisp[vecDisp < 0.001] <- NA
  vecDisp[vecDisp > 1000] <- NA
  # a) Plot dispersion parameter fit as function of mean parameter
  print("# a) Plot dispersion parameter fit as function of mean parameter")
  # From H0 to have one parameter each per gene.
  dfScatterDispvsMean <- data.frame(
    mu=log(vecMu)/log(10),
    disp=vecDisp )
  gScatterDispvsMean <- ggplot(dfScatterDispvsMean, aes(x=mu, y=disp)) +
    geom_point() + 
    labs(title="Dispersion parameter as function of\n mean parameter under H0.") +
    xlab(paste0("log10 mean parameter H0")) +
    ylab(paste0("dispersion parameter H0")) + 
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
  print("# b) Plot coefficient of variation as function of mean parameter")
  # Compute coefficient of variation under H0
  vecSD <- sqrt( vecMu + vecMu^2 / vecDisp )
  vecCV <- vecSD / vecMu
  # From H0 to have one parameter each per gene.
  dfScatterCVvsMean <- data.frame(
    mu=log(vecMu)/log(10),
    cv=vecCV )
  gScatterDispvsMean <- ggplot(dfScatterCVvsMean, aes(x=mu, y=cv)) +
    geom_point() + 
    labs(title="Coefficient of variation as function of\n mean parameter under H0.") +
    xlab(paste0("log10 mean parameter H0")) +
    ylab(paste0("coefficient of variation H0")) + 
    scale_fill_continuous(name = "Count") + 
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  pdf("LineagePulse_Scatter_CVvsMean.pdf",width=7,height=7)
  print(gScatterDispvsMean)
  dev.off()
  graphics.off()
  
  print("Done generating validation plots.")
  setwd(dirCurrent)
}