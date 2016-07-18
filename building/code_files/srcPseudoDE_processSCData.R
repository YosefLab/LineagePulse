#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Process single cell data    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Prepare single cell data for analysis
#' 
#' Check that input is correctly supplied and formatted. Then process input
#' data for analysis.
#' (I) Helper functions:
#'    checkNull() Check whether object was supplied (is not NULL).
#'    checkCounts() Checks whether elements are count data.
#'    checkNumeric() Checks whether elements are numeric.
#'    checkCounts() Checks whether elements are count data.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCounts: (matrix genes x samples)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @return matCountsProc: (matrix genes x samples)
#'    Processed count data of all cells, unobserved entries are NA.
#' @return vecPseudotimeProc: (numerical vector length number of cells)
#'    Processed pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Names of elements are cell names.
#' @export

processSCData <- function(matCounts,
  vecPseudotime,
  scaSmallRun=scaSmallRun,
  boolDEAnalysisImpulseModel,
  boolDEAnalysisModelFree,
  verbose ){
  
  # Check whether object was supplied (is not NULL).
  checkNull <- function(objectInput,strObjectInput){
    if(is.null(objectInput)){
      stop(paste0( "ERROR: ", strObjectInput," was not given as input." ))
    }
  }
  # Checks whether elements are numeric
  checkNumeric <- function(matInput, strMatInput){
    if(any(!is.numeric(matInput))){
      stop(paste0( "ERROR: ", strMatInput, " contains non-numeric elements. Requires count data." ))
    }
  }
  # Checks whether elements are count data: non-negative integer finite numeric elements.
  # Note that NA are allowed.
  checkCounts <- function(matInput, strMatInput){
    checkNumeric(matInput, strMatInput)
    if(any(matInput[!is.na(matInput)] %% 1 != 0)){
      stop(paste0( "ERROR: ", strMatInput, " contains non-integer elements. Requires count data." ))
    }
    if(any(!is.finite(matInput[!is.na(matInput)]))){
      stop(paste0( "ERROR: ", strMatInput, " contains infinite elements. Requires count data." ))
    }
    if(any(matInput[!is.na(matInput)]<0)){
      stop(paste0( "ERROR: ", strMatInput, " contains negative elements. Requires count data." ))
    }
  }

  # Name genes if names are not given.
  nameGenes <- function(matInput){
    if(is.null(rownames(matInput))){
      rownames(matInput) <- paste0("Region_", seq(1,nrow(matInput)))
    }
    if(any(is.na(rownames(matInput)))){
      scaNAs <- sum(is.na(rownames(matInput)))
      idxNAs <- is.na(rownames(matInput))
      (rownames(matInput))[idxNAs] <- paste0("Gene_", seq(1,scaNAs))
    }
    return(matInput)
  }
  
  # (I) Check input
  # 1. matCounts
  checkNull(matCounts,"matCounts")
  checkCounts(matCounts,"matCounts")
  
  # 2. vecPseudotime
  checkNull(vecPseudotime,"vecPseudotime")
  checkNumeric(vecPseudotime,"vecPseudotime")
  if(!all(names(vecPseudotime) %in% colnames(matCounts))){
    stop("ERROR: Not all cells in vecPseudotime are given matCounts.")
  }
  
  # 3. Check scaSmallRun
  if(!is.null(scaSmallRun)){
    checkCounts(scaSmallRun, "scaSmallRun")
    if(scaSmallRun > dim(matCounts)[1]){
      stop(paste0( "ERROR: scaSmallRun (",scaSmallRun,
        ") larger then data set (",dim(matCounts)[1]," genes)."))
    }
  }
  
  # 4. boolDEAnalysisImpulseModel and boolDEAnalysisModelFree
  if(!boolDEAnalysisImpulseModel & !boolDEAnalysisModelFree){
    stop(paste0("ERROR: No differential analysis mode selected.",
      "Set either boolDEAnalysisImpulseModel or boolDEAnalysisModelFree as TRUE."))
  }
  
  # (II) Process data
  # Convert from data frame to matrix for ziber()
  if(is.list(matCounts)){
    matCounts <- data.matrix(matCounts)
  }
  # Take out cells with NA pseudotime coordinate
  vecidxPT <- !is.na(vecPseudotime)
  vecPseudotimeProc <- vecPseudotime[vecidxPT]
  matCountsProc <- matCounts[,names(vecPseudotimeProc)]
  # Remove all zero or NA genes/cells
  vecidxGenes <- apply(matCountsProc, 1, function(gene){any(gene>0 & is.finite(gene) & !is.na(gene))})
  vecidxCells <- apply(matCountsProc, 2, function(cell){any(cell>0 & is.finite(cell) & !is.na(cell))})
  vecPseudotimeProc <- vecPseudotimeProc[vecidxCells]
  matCountsProc <- matCountsProc[vecidxGenes,vecidxCells]
  # Reduce data set to small run size if required.
  # Keep full data set for size factor estimation.
  matCountsProcFull <- matCountsProc
  if(!is.null(scaSmallRun)){
    scaNRows <- min(scaSmallRun,dim(matCountsProc)[1])
    matCountsProc <- matCountsProc[1:scaNRows,]
  }
  
  # Print summary of processing
  print(paste0(sum(!vecidxPT), " out of ", length(vecidxPT), " cells did not have a pseudotime coordinate and were excluded."))
  print(paste0(sum(!vecidxGenes), " out of ", length(vecidxGenes), " genes did not contain non-zero observations and are excluded from analysis."))
  print(paste0(sum(!vecidxCells), " out of ", length(vecidxCells), " cells did not contain non-zero observations and are excluded from analysis."))
  if(!is.null(scaSmallRun)){
    print(paste0("Operating on subset of data, set by scaSmallRun: ",scaSmallRun," out of ",dim(matCountsProcFull)[1]," genes."))
  }
  
  # Name nameless gene:
  matCountsProc <- nameGenes(matCountsProc)
  
  return(list(matCountsProc=matCountsProc,
    matCountsProcFull=matCountsProcFull,
    vecPseudotimeProc=vecPseudotimeProc))
}