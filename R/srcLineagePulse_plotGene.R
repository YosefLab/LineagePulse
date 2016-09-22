#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  Plot counts and model for one gene  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot counts and model for one gene
#' 
#' Plot counts and model for one gene as a scatter plot with line plots. 
#' Dropout rates can be given and are visualised as the colour of the observation
#' points. 
#' 
#' @seealso Called by \code{runLineagePulse} or separately by user.
#' 
#' @param vecCounts: (numeric vector length cells) Observed counts for given
#'    gene.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#' @param vecDropoutRates: (numeric vector length cells) [Default NULL]
#'    Dropout rates for observations.
#' @param vecImpulseModelParam: (numeric vector length impulse parameters=6)
#'    [Default NULL] Impulse model parameters for given gene.
#' @param vecImpulseModelRefParam: (numeric vector length impulse parameters=6)
#'    [Default NULL] Reference impulse model parameters for given gene.
#' @param strNameImpulseModelRef: (str) [Default NULL] 
#'    Name of reference impulse model for given gene.
#' @param scaConstModelParam: (scalar) [Default NULL]
#'    Constant model parameter for given gene: expression mean.
#' @param strGeneID: (str) Name of gene, used for title of plot.
#' @param strTitleSuffix: (str) String to be added to title.
#' 
#' @return gGenePlot: (ggplot object) Plot which can be printed or saved
#'    to pdf.
#' 
#' @export

plotGene <- function(vecCounts,
  vecPseudotime,
  vecDropoutRates=NULL,
  vecImpulseModelParam=NULL,
  vecImpulseModelRefParam=NULL,
  strNameImpulseModelRef=NULL,
  scaConstModelParam=NULL,
  strGeneID,
  strTitleSuffix=NULL){
  
  scaNumCells <- length(vecCounts)
  # Set drop-out rates as constant for visualistion if not given.
  if(is.null(vecDropoutRates)){ vecDropoutRates <- array(0, scaNumCells) }
  dfScatterCounts <- data.frame(
    pt=vecPseudotime,
    counts=vecCounts,
    dropout_rate=vecDropoutRates )
  gGenePlot <- ggplot(dfScatterCounts, aes(x=pt, y=counts)) +
    geom_point(aes(colour=dropout_rate)) + 
    labs(title=paste0(strGeneID," ",strTitleSuffix)) +
    xlab(paste0("pseudotime")) +
    ylab(paste0("counts")) +
    scale_colour_gradient(high="red",low="green",limits=c(0, 1)) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  
  # Add models to plot
  scaNPoints <- 100 # Points to plot for each model
  vecPTCoord <- seq(min(vecPseudotime), max(vecPseudotime),
    by=(max(vecPseudotime)-min(vecPseudotime))/(scaNPoints-1) )
  # Impulse model
  if(!is.null(vecImpulseModelParam)){
    vecImpulseValues <- calcImpulse_comp(vecTheta=vecImpulseModelParam,
      vecTimepoints=vecPTCoord)
    dfLineImpulse <- data.frame( pt=vecPTCoord,
      impulse=vecImpulseValues)
    gGenePlot <- gGenePlot + geom_line(data=dfLineImpulse, 
      aes(x=pt, y=impulse))
  }
  # Reference impulse model (e.g. true model if handling simulated data)
  if(!is.null(vecImpulseModelRefParam)){
    vecImpulseValuesRef <- calcImpulse_comp(vecTheta=vecImpulseModelRefParam,
      vecTimepoints=vecPTCoord)
    dfLineImpulse <- data.frame( pt=vecPTCoord,
      impulse=vecImpulseValuesRef)
    gGenePlot <- gGenePlot + geom_line(data=dfLineImpulse, 
      aes(x=pt, y=impulse))
  }
  # Constant model
  if(!is.null(scaConstModelParam)){
    vecConstValues <- c(scaConstModelParam,scaConstModelParam)
    vecPTCoordConst <- c(min(vecPseudotime), max(vecPseudotime))
    dfLineConst <- data.frame( pt=vecPTCoordConst,
      const=vecConstValues)
    gGenePlot <- gGenePlot + geom_line(data=dfLineConst, 
      aes(x=pt, y=const))
  }
  
  return(gGenePlot)
}