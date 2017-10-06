#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++     Compute Size factors    ++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute size factors for a LineagePulse-object
#' 
#' Either use externally supplied normalisation constants or set these to one.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param objLP: (LineagePulse-object)
#' Object to fit normalization constants on.
#' @param vecNormConstExternal: (numeric vector number of cells) 
#' Model scaling factors supplied by user, one per cell. 
#' 
#' @return vecNormConst: (numeric vector number of cells) 
#' Model scaling factors to be used, one per cell.
#' 
#' @author David Sebastian Fischer
calcNormConst <- function(objLP,
                          vecNormConstExternal){
  
  if(!is.null(vecNormConstExternal)){
    vecNormConst <- vecNormConstExternal
  } else {
    print("# All size factors are set to one.")
    vecNormConst <- array(1, dim(objLP@matCountsProc)[2])
    names(vecNormConst) <- colnames(objLP@matCountsProc)
  }
  
  if(any(vecNormConst==0)){
    warning("WARNING IN LINEAGEPULSE: Found size factors==0, setting these to 1.")
    vecNormConst[vecNormConst==0] <- 1
  }
  
  objLP@vecNormConst <- as.vector(vecNormConst)
  return(objLP)
}