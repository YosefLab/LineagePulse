#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++     Plot pseudotime clustering of cells    +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot pseudotime clustering of cells
#' 
#' Summarises the 1D clustering of cells in pseudotime by two plots:
#' The empirical distribution function of cells in pseudotime and
#' the empricial cumulative distribution function of cells in 
#' pseudotime with cluster borders indicated in each.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param objectLineagePulse: (LineagePulseObject)
#' LineagePulse output object which was run with 
#' strMuModel==clusters.
#' 
#' @return (ggplot object)
#' Cell clustering in pseudotime visualisation.
#'    
#' @author David Sebastian Fischer
#' 
#' @export
plotPseudotimeClustering <- function(objectLineagePulse){
  
  dfClusterBorders <- data.frame( borders=sapply( seq(1,max(objectLineagePulse@lsMuModelH1$lsMuModelGlobal$vecClusterAssign)-1), 
    function(centroid){(lsResultsClustering$Centroids[centroid]+lsResultsClustering$Centroids[centroid+1])/2} ))
  dfPseudotime <- data.frame( pseudotime=as.vector(vecPseudotime) )
  
  plotEDF <- ggplot() +
    geom_density(data=dfPseudotime, aes(x=pseudotime), colour="black", bw=1) +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Density estimation of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical probability density")
  
  plotECDF <- ggplot() +
    stat_ecdf(data=dfPseudotime, aes(x=pseudotime), colour="red") +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Empirical cumulative density of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical cumulative probability density")
  
  pdf(strPDFname)
  print(plotEDF)
  print(plotECDF)
  dev.off()
  
  return(NULL)
}