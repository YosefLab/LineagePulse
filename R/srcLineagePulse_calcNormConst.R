#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++     Compute Size factors    ++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute size factors for a scRNA-seq dataset
#' 
#' Either use externally supplied normalisation constants or set these to one.
#' 
#' Size factors scale the negative binomial likelihood
#' model of a gene to the sequencing depth of each sample.
#' As of now, no stable normalisation constants are known for
#' scRNAseq data and they are set to 1 here.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCountsProc: (matrix genes x samples)
#'    Count data.
#' @param vecNormConstExternal: (numeric vector number of cells) 
#'    Model scaling factors, one per cell. These factors will linearly 
#'    scale the mean model for evaluation of the loglikelihood. 
#'    Must be named according to the column names of matCounts.
#'    Supplied by user.
#' 
#' @return vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors, one per cell.
#'    
#' @author David Sebastian Fischer
#' 
#' @export

calcNormConst <- function(objectLineagePulse,
                          vecNormConstExternal){
  
  if(!is.null(vecNormConstExternal)){
    vecNormConst <- vecNormConstExternal
  } else {
    print("All size factors are set to one.")
    vecNormConst <- array(1, dim(objectLineagePulse@matCountsProc)[2])
    names(vecNormConst) <- colnames(objectLineagePulse@matCountsProc)
  }
  
  if(any(vecNormConst==0)){
    warning("WARNING IN LINEAGEPULSE: Found size factors==0, setting these to 1.")
    vecNormConst[vecNormConst==0] <- 1
  }
  
  objectLineagePulse@vecNormConst <- as.vector(vecNormConst)
  return(objectLineagePulse)
}