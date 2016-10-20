#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++  Plot comparative q-value ECDF  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot comparative q-value ECDF
#' 
#' Loads DE results data frame from list of given LineagePulse ouput folders
#' and plots ECDF for each output into single plot.
#' 
#' @seealso Auxillary method not called by LineagePulse wrapper.
#' Called separately by user.
#' 
#' @param lsDirOutLineagePulse: (list) [Default NULL]
#'    List of LineagePulse output directories from which
#'    dfDEResults is loaded.
#' @param filePDF: (directory/filename)
#'    File to which plot is saved (pdf).
#' 
#' @return NULL
#' 
#' @export

plotComparativeECDF <- function(
  lsDirOutLineagePulse,
  filePDF=NULL,
  strRunLabels,
  strTitle="Comparative Q-value ECDF"){
  
  lsVecECDF <- list()
  vecX <- seq(-100,-1,by=0.1)
  for(dirLP in lsDirOutLineagePulse){
    load(paste0(dirLP,"/LineagePulse_dfDEAnalysis.RData"))
    lsVecECDF[length(lsVecECDF)+1] <- list(sapply(vecX, function(thres){
      sum(log(as.numeric(as.vector(dfDEAnalysis$adj.p)))/log(10) <= thres, na.rm=TRUE)}))
  }
  
  # Generate data frame for ggplot
  dfECDF <- data.frame(
    logq=rep(vecX, length(lsVecECDF)),
    counts=unlist(lsVecECDF),
    run=unlist(lapply(strRunLabels, function(run) rep(run, length(vecX)))) )
  # ggplot
  gECDF <- ggplot(dfECDF, aes(x=logq, y=counts, group=run)) +
    geom_line(aes(colour=run)) +
    labs(title=strTitle) +
    xlab(paste0("log_10(q-value)")) +
    ylab(paste0("Cumulative q-value distribution")) +
    xlim(-15,-1) + 
    ylim(0,max(unlist(lsVecECDF))) +
    theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
  # lot
  if(!is.null(filePDF)){
    pdf(filePDF,width=7,height=7)
    print(gECDF)
    dev.off()
  }
  print(gECDF)
  
  return(NULL)
}