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
#' @param lsObjectLineagePulse: (list of LineagePulseObjects)
#' List of LineagePulseObjects with DE results table to 
#' base ECDF on.
#' @param filePDF: (directory/filename)
#' File to which plot is saved (pdf).
#' @param strRunLabels: (string vector arbitrary length)
#' Names of runs to be plotted. Correspond
#' to directories supplied in lsDirOutLineagePulse.
#' @param strTitle: (string) Title of plot.
#' 
#' @return (ggplot2 object) ECDF plot
#' 
#' @author David Sebastian Fischer
#' 
#' @export
plotComparativeECDF <- function(
    lsObjectLineagePulse,
    filePDF=NULL,
    strRunLabels,
    strTitle="Comparative Q-value ECDF"){
    
    if(length(lsDirOutLineagePulse) != length(strRunLabels)){
        stop("Number of supplied lsObjectLineagePulse and run names dont match.")
    }
    lsVecECDF <- list()
    vecNAnalysed <- array(NA, length(lsObjectLineagePulse))
    vecX <- seq(-100,-1,by=0.1)
    for(object in lsObjectLineagePulse){
        lsVecECDF[length(lsVecECDF)+1] <- list(sapply(vecX, function(thres){
            sum(log(as.numeric(as.vector(object@dfDEAnalysis$adj.p)))/log(10) <= thres, na.rm=TRUE)}))
        vecNAnalysed[length(lsVecECDF)] <- sum(!is.na(as.numeric(as.vector(object@dfDEAnalysis$adj.p))))
    }
    
    # Generate data frame for ggplot
    dfECDF <- data.frame(
        logq=rep(vecX, length(lsVecECDF)),
        counts=unlist(lsVecECDF),
        data_set=unlist(lapply( seq(1, length(strRunLabels)), function(i){
            rep( paste(strRunLabels[i], "[", vecNAnalysed[i],"]"),  
                 length(vecX))
        })) )
    # ggplot
    gECDF <- ggplot(dfECDF, aes(x=logq, y=counts, group=data_set)) +
        geom_line(aes(colour=data_set)) +
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
    return(gECDF)
}