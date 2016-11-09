#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  Plot counts and model for one gene  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot counts and model for one gene
#' 
#' Plot counts and model for one gene as a scatter plot with line plots. 
#' Dropout rates can be given and are visualised as the colour of the observation
#' points. Inferred and underlying reference models can be plotted as lines.
#' 
#' @seealso Called by \code{runLineagePulse} or separately by user.
#' 
#' @param vecCounts: (numeric vector length cells) Observed counts for given
#'    gene.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#' @param vecPiParamH1: (numeric vector length cells) [Default NULL]
#'    Dropout rates for observations under alternative model.
#' @param vecImpulseModelParam: (numeric vector length impulse parameters=6)
#'    [Default NULL] Impulse model parameters for given gene. 
#'    Amplitudes in log format, like LineagePulse output!
#' @param vecImpulseModelRefParam: (numeric vector length impulse parameters=6)
#'    [Default NULL] Reference impulse model parameters for given gene.
#'    Amplitudes in log format, like LineagePulse output!
#' @param strNameImpulseModelRef: (str) [Default NULL] 
#'    Name of reference impulse model for given gene.
#' @param scaConstModelParam: (scalar) [Default NULL]
#'    Constant model parameter for given gene: expression mean.
#' @param scaConstModelRefParam: (scalar) [Default NULL]
#'    Reference constant model parameter for given gene: expression mean.
#' @param strGeneID: (str) Name of gene, used for title of plot.
#' @param strTitleSuffix: (str) String to be added to title.
#' 
#' @return gGenePlot: (ggplot object) Plot which can be printed or saved
#'    to pdf. Inferred impulse and constant model are plotted as 
#'    a constant line, the underlying models as dashed.
#' 
#' @export

plotGene <- function(vecCounts,
  vecPseudotime,
  vecPiParamH1=NULL,
  vecImpulseModelParam=NULL,
  vecImpulseModelRefParam=NULL,
  strNameImpulseModelRef=NULL,
  scaConstModelParam=NULL,
  scaConstModelRefParam=NULL,
  strGeneID,
  strTitleSuffix=NULL,
  boolScaleByLL=FALSE,
  vecMuParamH0=NULL,
  vecMuParamH1=NULL,
  vecDispParamH0=NULL,
  vecDispParamH1=NULL,
  vecPiParamH0=NULL,
  vecNormConst=NULL,
  scaWindowRadius=NULL){
  
  scaNumCells <- length(vecCounts)
  # Set drop-out rates as constant for visualistion if not given.
  if(is.null(vecPiParamH1)){ vecPiParamH1 <- array(0, scaNumCells) }
  
  if(boolScaleByLL){
    # Compute observation-wise loglikelihoods under null model
    vecLogLikH0 <- sapply(seq(1,scaNumCells), 
      function(j){
        if(is.null(scaWindowRadius)){
          vecInterval <- j
        } else {
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        }
        
        scaLogLikCell <- evalLogLikGene(
          vecCounts=vecCounts[vecInterval],
          vecMu=vecMuParamH0[j]*vecNormConst[vecInterval],
          vecNormConst=vecNormConst,
          vecDisp=rep(vecDispParamH0[j], length(vecInterval)), 
          vecPi=vecPiParamH0[vecInterval], 
          vecboolNotZero= !is.na(vecCounts[vecInterval]) & vecCounts[vecInterval]>0, 
          vecboolZero= !is.na(vecCounts[vecInterval]) & vecCounts[vecInterval]==0,
          scaWindowRadius=scaWindowRadius )
        return(scaLogLikCell)
      }) 
    
    # Compute observation-wise loglikelihoods under alternative model
    vecLogLikH1 <- sapply(seq(1,scaNumCells), 
      function(j){
        if(is.null(scaWindowRadius)){
          vecInterval <- j
        } else {
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        }
        
        scaLogLikCell <- evalLogLikGene(
          vecCounts=vecCounts[vecInterval],
          vecMu=vecMuParamH1[j]*vecNormConst[vecInterval],
          vecNormConst=vecNormConst,
          vecDisp=rep(vecDispParamH1[j], length(vecInterval)), 
          vecPi=vecPiParamH1[vecInterval], 
          vecboolNotZero= !is.na(vecCounts[vecInterval]) & vecCounts[vecInterval]>0, 
          vecboolZero= !is.na(vecCounts[vecInterval]) & vecCounts[vecInterval]==0,
          scaWindowRadius=scaWindowRadius )
        return(scaLogLikCell)
      })
    vecLogLikRatio <- vecLogLikH1-vecLogLikH0

    dfScatterCounts <- data.frame(
      pt=vecPseudotime,
      counts=vecCounts,
      dropout_rate=vecPiParamH1,
      loglik_ratio=vecLogLikRatio)
    gGenePlot <- ggplot(dfScatterCounts, aes(x=pt, y=counts)) +
      geom_point(aes(colour=dropout_rate, size=loglik_ratio), show.legend=TRUE)
  } else {
    dfScatterCounts <- data.frame(
      pt=vecPseudotime,
      counts=vecCounts,
      dropout_rate=vecPiParamH1)
    gGenePlot <- ggplot(dfScatterCounts, aes(x=pt, y=counts)) +
      geom_point(aes(colour=dropout_rate), show.legend=TRUE)
  }
  
  gGenePlot <- gGenePlot + 
    labs(title=paste0(strGeneID," ",strTitleSuffix)) +
    xlab(paste0("pseudotime")) +
    ylab(paste0("counts")) +
    scale_colour_gradient(high="red",low="green",limits=c(0, 1)) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14)) +
    annotate("text", x=(max(vecPseudotime)-min(vecPseudotime))*3/4,
      y=max(vecCounts)*9/10, label="[Line legend]\nUnderlying model: dashed,\nInferred model: solid")
  # Set plotting threshold based on observed data
  scaMaxPlot <- 2 * max(vecCounts)
  
  # Add models to plot
  scaNPoints <- 100 # Points to plot for each model
  vecPTCoord <- seq(min(vecPseudotime), max(vecPseudotime),
    by=(max(vecPseudotime)-min(vecPseudotime))/(scaNPoints-1) )
  # Impulse model
  if(!is.null(vecImpulseModelParam)){
    vecImpulseValues <- evalImpulseModel_comp(vecTheta=vecImpulseModelParam,
      vecTimepoints=vecPTCoord)
    vecImpulseValues[vecImpulseValues > scaMaxPlot] <- NA
    dfLineImpulse <- data.frame( pt=vecPTCoord,
      impulse=vecImpulseValues)
    gGenePlot <- gGenePlot + 
      geom_line(data=dfLineImpulse, 
        aes(x=pt, y=impulse),
        linetype="solid",
        show.legend=TRUE)
  }
  # Reference impulse model (e.g. true model if handling simulated data)
  if(!is.null(vecImpulseModelRefParam)){
    vecImpulseValuesRef <- evalImpulseModel_comp(vecTheta=vecImpulseModelRefParam,
      vecTimepoints=vecPTCoord)
    vecImpulseValuesRef[vecImpulseValuesRef > scaMaxPlot] <- NA
    dfLineImpulse <- data.frame( pt=vecPTCoord,
      impulse=vecImpulseValuesRef)
    gGenePlot <- gGenePlot + geom_line(data=dfLineImpulse, 
      aes(x=pt, y=impulse),
      linetype="dashed",
      show.legend=TRUE)
  }
  # Constant model
  if(!is.null(scaConstModelParam)){
    vecConstValues <- c(scaConstModelParam,scaConstModelParam)
    vecConstValues[vecConstValues > scaMaxPlot] <- NA
    vecPTCoordConst <- c(min(vecPseudotime), max(vecPseudotime))
    dfLineConst <- data.frame( pt=vecPTCoordConst,
      const=vecConstValues)
    gGenePlot <- gGenePlot + geom_line(data=dfLineConst,
      aes(x=pt, y=const),
      linetype="solid",
      show.legend=TRUE)
  }
  # Reference constant model
  if(!is.null(scaConstModelRefParam)){
    vecConstRefValues <- c(scaConstModelRefParam,scaConstModelRefParam)
    vecConstRefValues[vecConstRefValues > scaMaxPlot] <- NA
    vecPTCoordConst <- c(min(vecPseudotime), max(vecPseudotime))
    dfLineConstRef <- data.frame( pt=vecPTCoordConst,
      const=vecConstRefValues)
    gGenePlot <- gGenePlot + geom_line(data=dfLineConstRef, 
      aes(x=pt, y=const),
      linetype="dashed",
      show.legend=TRUE)
  }
  
  return(gGenePlot)
}