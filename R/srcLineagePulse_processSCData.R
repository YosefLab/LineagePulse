#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Process single cell data    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Prepare single cell data for analysis
#' 
#' Check that input is correctly supplied and formatted. Then process input
#' data for analysis.
#' Helper functions:
#'    checkNull() Check whether object was supplied (is not NULL).
#'    checkNumeric() Checks whether elements are numeric.
#'    checkCounts() Checks whether elements are count data.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCounts: (matrix genes x samples)
#'    Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param strMuModel: (str) {"constant"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param boolCoEstDispMean: (bool) [Default TRUE]
#'    Whether mean and dispersion parameters are to be co-estimated
#'    (simulatneous optimisation). Only available for certain 
#'    dispersion and mean models:
#'    dispersion models: constant.
#'    mean models: constant, cluster, sliding window vector, impulse.
#'    Note that co-estimation in model estimation B (without drop-
#'    out model estimation) leads to a single step estimation as 
#'    mean and dispersion parameter don't have to be iterated over.
#'    This makes estimation of large data sets with complex H1 mean
#'    model (e.g. impulse) possible, as the drop-out model can be 
#'    estimated based on H0 (boolEstimateNoiseBasedOnH0) so that
#'    the complex model only has to be estimated once (simultaneous
#'    with the dispersion parameters). This may generally lead to better
#'    convergence as the steps in coordinate-ascent are in a larger
#'    space, closer to full gradient ascent. Setting to TRUE
#'    is encouraged. Optimisation routines for individual mean 
#'    and dispersion fitting (if FALSE) exist, but these may be viewed
#'    as non-deprecated parts of an earlier implementation of the
#'    alorithm.
#' @param boolVecWindowsAsBFGS: (bool) [Default FALSE] Whether
#'    mean parameters of a gene are co-estimated in "windows"
#'    mode with BFGS algorithm (optimisation with dimensionality
#'    of number of cells) or estimated one by one, conditioned
#'    one the latest estimates of neighbours. The latter case
#'    (boolVecWindowsAsBFGS=FALSE) is coordinate ascent within the gene
#'    and each mean parameter is optimised once only.
#' @param dirOut: (str directory) [Default NULL]
#'    Directory to which detailed output is saved to.
#'    Defaults to current working directory if NULL.
#'
#' @return list: (length 4)
#'    \itemize{
#'        \item matCountsProc: (matrix genes x samples)
#'    Processed count data of all cells, unobserved entries are NA,
#'    trimmed by scaSmallRun.
#'        \item matCountsProcFull: (matrix genes x samples)
#'    Processed count data of all cells, unobserved entries are NA.
#'        \item vecPseudotimeProc: (numerical vector length number of cells)
#'    Processed pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Names of elements are cell names.
#'        \item dirOut: (str directory)
#'    Directory to which detailed output is saved to.
#'    Defaults to current working directory if NULL.
#'    }

#' @export

processSCData <- function(matCounts,
  matPiConstPredictors,
  vecPseudotime,
  scaSmallRun,
  strMuModel,
  strDispModel,
  scaWindowRadius,
  boolCoEstDispMean,
  boolVecWindowsAsBFGS,
  dirOut ){
  
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
  
  # (I) Check input data
  # 1. matCounts
  checkNull(matCounts,"matCounts")
  checkCounts(matCounts,"matCounts")
  
  # 2. matPiConstPredictors
  if(!is.null(matPiConstPredictors)){
    checkNumeric(matPiConstPredictors,"matPiConstPredictors")
    if(!is.null(rownames(matCounts))){
      if(!rownames(matCounts) %in% rownames(matPiConstPredictors)){
        stop(paste0("ERROR: Some genes named in rows of matCounts do not ",
          "occur in rows of matPiConstPredictors."))   
      }
    } else {
      if(!is.null(rownames(matPiConstPredictors))){
        stop(paste0("ERROR: Named genes in matPiConstPredictors",
          " but not in matCounts."))
      }
    }
  }
  
  # 3. vecPseudotime
  checkNull(vecPseudotime,"vecPseudotime")
  checkNumeric(vecPseudotime,"vecPseudotime")
  if(!all(names(vecPseudotime) %in% colnames(matCounts))){
    stop("ERROR: Not all cells in vecPseudotime are given matCounts.")
  }
  
  # (II) Check settings
  # 1. Check scaSmallRun
  if(!is.null(scaSmallRun)){
    checkCounts(scaSmallRun, "scaSmallRun")
    if(scaSmallRun > dim(matCounts)[1]){
      stop(paste0( "ERROR: scaSmallRun (",scaSmallRun,
        ") larger then data set (",dim(matCounts)[1]," genes)."))
    }
  }
  
  # 2. Check general mean estimation settings
  if(!is.null(scaWindowRadius)){
    if(scaWindowRadius==0){
      stop(paste0("LineagePulse received scaWindowRadius=0.", 
        " Set to NULL if smoothing is not desired."))
    }
  }
  if(!is.null(scaWindowRadius) & !(strMuModel %in% c("windows","impulse","constant"))){
    stop(paste0("Smooting via scaWindowRadius can only be applied in strMuModel=",
      "constant, windows and impulse. Set scaWindowRadius=NULL or adjust strMuModel.",
      " Given: strMuModel=", strMuModel, " scaWindowRadius=", scaWindowRadius))
  }
  if(boolVecWindowsAsBFGS & !strMuModel=="windows"){
    stop(paste0("boolVecWindowsAsBFGS set to TRUE but strMuModel=",
      strMuModel, ". boolVecWindowsAsBFGS is only used",
      " with strMuModel=windows."))
  }
  if(strMuModel=="windows" & is.null(scaWindowRadius)){
    stop(paste0("Cannot use strMuModel=windows with scaWindowRadius=NULL.",
      " Means have to be regularised by smoothing in windows mode."))
  }
  
  # 3. Check single mean-/dispersion estimation models
  if(!boolCoEstDispMean){
    # a) Mean models
    if(!(strMuModel %in% c("windows","clusters","impulse","constant"))){
      stop(paste0("strMuModel not recognised: ", strMuModel, 
        " Must be one of: windows, clusters, impulse, constant."))
    }
    # b) Dispersion models
    if(!(strDispModel %in% c("constant"))){
      stop(paste0("strDispModel not recognised: ", strDispModel, 
        " Must be one of: constant."))
    }
  }
  
  # 4. Check co-estimation models
  # Note: These functions are implemented separately from
  # single mean and dispersion estimation.
  if(boolCoEstDispMean){
    if(!(strMuModel %in% c("windows","clusters","impulse","constant"))){
      stop(paste0("strMuModel not recognised: ", strMuModel, 
        " Must be one of: windows, clusters, impulse, constant if co-estimation",
        " of mean and dispersion is used (boolCoEstDispMean=TRUE)."))
    }
    if(strMuModel=="windows" & !boolVecWindowsAsBFGS){
      stop(paste0("The dispersion parameter co-estimation couples all LL",
        " terms of a gene so that there is no computational",
        " advantage in not estimating all mean parameters",
        " simultaneously. Set boolVecWindowsAsBFGS=TRUE",
        " if (strMuModel==windows and boolCoEstDispMean=TRUE)."))
    }
    if(strMuModel=="windows")
      if(!(strDispModel %in% c("constant"))){
        stop(paste0("strDispModel not recognised: ", strDispModel, 
          " Must be one of: constant if co-estimation",
          " of mean and dispersion is used (boolCoEstDispMean=TRUE)."))
      }
  }
  
  # (III) Process data
  # Convert from data frame to matrix
  if(is.list(matCounts)){
    matCounts <- data.matrix(matCounts)
  }
  # Name genes if names not given
  if(is.null(rownames(matCounts))){
    rownames(matCounts) <- paste0("Gene_", seq(1,nrow(matCounts)))
    rownames(matPiConstPredictors) <- rownames(matCounts)
  }
  # Take out cells with NA pseudotime coordinate
  vecidxPT <- !is.na(vecPseudotime)
  vecPseudotimeProc <- sort(vecPseudotime[vecidxPT])
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
  
  # Reduce matPiConstPredictors to genes in matCountsProc
  matPiConstPredictorsProc <- matPiConstPredictors[rownames(matCountsProc)]
  
  # Process output directory
  if(is.null(dirOut)){
    dirOut <- getwd()
  }
  
  return(list(matCountsProc=matCountsProc,
    matCountsProcFull=matCountsProcFull,
    matPiConstPredictorsProc=matPiConstPredictorsProc,
    vecPseudotimeProc=vecPseudotimeProc,
    dirOut=dirOut ))
}