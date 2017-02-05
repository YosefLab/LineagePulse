#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  Plot counts and model for one gene  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot counts and model for one gene
#' 
#' Plot counts and model for one gene as a scatter plot with model trajectories:
#' Both as counts or log10 counts vs pseudotime.
#' Dropout rates can be given and are visualised as the colour of the observation
#' points.
#' A reference model trajectory can be added.
#' 
#' @param objectLineagePulse: (LineagePulseObject) LineagePulseObject
#'    base plot on.
#' @param strGeneID: (str) Name of gene, used for title of plot.
#' @param vecReferenceMuParam: (numeric vector length number of cells)
#'    [Default NULL] Reference mean trajectory which can be plotted
#' @param strTitleSuffix: (str) String to be added to title.
#' @param boolColourByDropout: (bool) Whether to colour scatter
#'    plot by posterior of drop-out.
#' 
#' @return gGenePlot: (ggplot object)
#'    Model rajectories and scatter plot for given gene.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
plotGene <- function(objectLineagePulse,
                     strGeneID,
                     vecReferenceMuParam=NULL,
                     strTitleSuffix=NULL,
                     boolLogPlot=TRUE,
                     boolColourByDropout=TRUE){
  
  ### 1. Extract data and models
  vecCounts <- objectLineagePulse@matCountsProc[strGeneID,]
  vecCountsNorm <- vecCounts
  for(confounder in seq(1, length(objectLineagePulse@lsMuModelH0$lsMuModelGlobal$vecConfounders))){
    vecCountsNorm <- vecCountsNorm/(objectLineagePulse@lsMuModelH1$lsmatBatchModel[[confounder]][strGeneID,][
      objectLineagePulse@lsMuModelH1$lsMuModelGlobal$lsvecidxBatchAssign[[confounder]]])
  }
  vecMuParamH0 <- decompressMeansByGene( vecMuModel=objectLineagePulse@lsMuModelH0$matMuModel[strGeneID,],
                                         lsvecBatchModel=NULL,
                                         lsMuModelGlobal=objectLineagePulse@lsMuModelH0$lsMuModelGlobal,
                                         vecInterval=NULL )
  vecMuParamH1 <- decompressMeansByGene( vecMuModel=objectLineagePulse@lsMuModelH1$matMuModel[strGeneID,],
                                         lsvecBatchModel=NULL,
                                         lsMuModelGlobal=objectLineagePulse@lsMuModelH1$lsMuModelGlobal,
                                         vecInterval=NULL )
  vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=objectLineagePulse@lsDropModel$matDropoutLinModel,
                                               vecMu=vecMuParamH1,
                                               vecPiConstPredictors=objectLineagePulse@lsDropModel$matPiConstPredictors[strGeneID,] )
  if(boolLogPlot){
    vecCountsNorm <- log(vecCountsNorm)/log(10)
    vecMuParamH0 <- log(vecMuParamH0)/log(10)
    vecMuParamH1 <- log(vecMuParamH1)/log(10)
  }
  
  ### 2. Data scatter plot
  # Set drop-out rates as constant for visualistion if not given.
  if(boolColourByDropout){
    dfScatterCounts <- data.frame(
      pseudotime=objectLineagePulse@lsMuModelH1$lsMuModelGlobal$vecPseudotime,
      counts=vecCountsNorm,
      dropout_rate=vecPiParamH1)
    gGenePlot <- ggplot() +
      geom_point(data=dfScatterCounts, aes(x=pseudotime, y=counts, colour=dropout_rate), show.legend=TRUE)
  } else {
    dfScatterCounts <- data.frame(
      pseudotime=objectLineagePulse@lsMuModelH1$lsMuModelGlobal$vecPseudotime,
      counts=vecCountsNorm)
    gGenePlot <- ggplot() +
      geom_point(data=dfScatterCounts, aes(x=pseudotime, y=counts), show.legend=TRUE)
  }
  
  gGenePlot <- gGenePlot + 
    labs(title=paste0(strGeneID, "\nlog10 q-value=", round(objectLineagePulse@dfResults[strGeneID,]$adj.p,2) )) +
    xlab(paste0("pseudotime")) +
    scale_colour_gradient(high="red",low="green",limits=c(0, 1)) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14))
  if(boolLogPlot) gGenePlot <- gGenePlot + ylab(paste0("log10 counts"))
  else gGenePlot <- gGenePlot + ylab(paste0("counts"))
  
  # Set plotting threshold based on observed data
  scaMaxPlot <- 2 * max(vecCounts)
  
  ### 3. Add models to plot
  vecMuParamH0[vecMuParamH0 > scaMaxPlot] <- NA
  vecMuParamH1[vecMuParamH1 > scaMaxPlot] <- NA
  if(is.null(vecReferenceMuParam)){
    dfLineImpulse <- data.frame( pseudotime=rep(objectLineagePulse@lsMuModelH1$lsMuModelGlobal$vecPseudotime, 2),
                                 counts=c(vecMuParamH0,
                                          vecMuParamH1),
                                 model=c(rep("H0", length(vecMuParamH0)), 
                                         rep("H1", length(vecMuParamH1))) )
  } else {
    dfLineImpulse <- data.frame( pseudotime=rep(objectLineagePulse@lsMuModelH1$lsMuModelGlobal$vecPseudotime, 3),
                                 counts=c(vecMuParamH0,
                                          vecMuParamH1,
                                          vecReferenceMuParam),
                                 model=c(rep("H0", length(vecMuParamH0)), 
                                         rep("H1", length(vecMuParamH1)),
                                         rep("Reference", length(vecReferenceMuParam))) )
  }
  gGenePlot <- gGenePlot + geom_line(data=dfLineImpulse, 
                                     aes(x=pseudotime, y=counts, linetype=model),
                                     show.legend=TRUE)
  
  return(gGenePlot)
}