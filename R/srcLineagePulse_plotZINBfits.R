#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Plot ZINB fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param vecGeneIDs: (string vector) Gene names to be plotted. 
#'    Elements must be rownames of matCounts.
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecMu: (numeric vector number of genes)
#'    Constant mean parameter estimate (H0) of genes.
#'    Must be contain names of vecGeneIDs as labels.
#' @param vecDispersions: (numeric vector number of genes)
#'    Constant dispersion parameter estimate (H0) of genes.
#'    Must be contain names of vecGeneIDs as labels.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matProbNB: (probability matrix genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param strPDFname: (str) Name of .pdf with plots.
#' 
#' @return NULL
#' @export

plotZINBfits <- function(vecGeneIDs, 
  matCounts,
  vecMu, 
  vecDispersions,
  matProbNB,
  scaWindowRadius,
  strPDFname="LineagePulse_ZINBfits.pdf"){
  
  scaNumCells <- dim(matCounts)[2]
    
  # Only keep target genes in input parameters
  matCounts <- matCounts[vecGeneIDs,]
  vecMu <- vecMu[vecGeneIDs]
  vecDispersions <- vecDispersions[vecGeneIDs]
  matProbNB <- matProbNB[vecGeneIDs,]
  
  # Only for batch/singlecell plotting:
  # Width of negative binomial pdf (ylim) in time units for plotting
  # E.g. PDF_WIDTH <- 1 means that the inferred negative binomial
  # density of a sample will be scaled into the interval [0,1]
  # and plotted at sample x time units vertically in the interval
  # [x,x+1] time units, with the y axis (the counts) being the
  # horizontal axis of the pdf.
  #PDF_WIDTH <- 1
  
  # Scale count data by size factors for plotting:
  # The impulse models are fit based on normalised means.
  # Therefore, the model curves follow the normalised
  # count data and not the raw count data. However, 
  # fitting was still performed based on raw count data.
  #matCounts <- matCounts / matrix(vecNormConst, 
  #  nrow=dim(matCounts)[1], 
  #  ncol=dim(matCounts)[2], 
  #  byrow=TRUE)
  
  # Define grid for printing plots
  par(mfrow=c(3,3), xpd=TRUE)
  scaLegendInset <- -0.65
  
  # Create list of plots
  lsPlots <- list()
  for (geneID in vecGeneIDs){
    vecXCoordPDF <- seq(0, max(matCounts[geneID,],na.rm=TRUE))
    vecYCoordPDF <- dnbinom(vecXCoordPDF,
      mu=rep(vecMu[geneID],scaNumCells),
      size=rep(vecDispersions[geneID],scaNumCells))
    # Plot both
    vecMixtures <- c("Negative binomial","Dropout")
    dfData <- data.frame( counts=as.vector(matCounts[geneID,] ),
      mixture=vecMixtures[as.numeric(matProbNB[geneID,] < 0.5)+1]  )
    dfDensity <- data.frame( counts=vecXCoordPDF, 
      density=vecYCoordPDF*(dim(dfData)[1]))
    lsPlots[[length(lsPlots)+1]] <- suppressMessages( 
      ggplot( ) +
        geom_histogram( data=dfData, 
          aes(x=counts, y=..count.., fill=mixture), 
          alpha=0.5 ) +
        geom_line( data=dfDensity, 
          aes(x=counts, y=density), 
          col = "blue") +
        labs(title=paste0(geneID," Observed cells: ",length(dfData$counts),
          "\n Non-zero counts: ", sum(dfData$counts!=0))) +
        xlab("Counts") +
        ylab("Frequency") 
    )
  }
  
  # Write plots into .pdf
  pdf(strPDFname,height=6.0,width=9.0)
  for(plotGene in lsPlots){
    suppressMessages( print(plotGene) )
  }
  dev.off()
}