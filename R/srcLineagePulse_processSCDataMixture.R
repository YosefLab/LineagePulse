#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Process single cell data    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Prepare single cell data for analysis
#' 
#' Check that input is correctly supplied and formatted. Then process input
#' data for analysis. This function catches downstream errors arising from 
#' incorrect input parameter combinations and explains the input error
#' to the user. Input passing this function should not lead to errors
#' arising from incorrect input.
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
#' @param vecNormConstExternal: (numeric vector number of cells) 
#'    Model scaling factors, one per cell. These factors will linearly 
#'    scale the mean model for evaluation of the loglikelihood. 
#'    Must be named according to the column names of matCounts.
#'    Supplied by user.
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
#'    
#' @author David Sebastian Fischer
#' 
#' @export
processSCDataMixture <- function(matCounts,
                                 vecFixedAssignments,
                                 matPiConstPredictors,
                                 vecNormConstExternal,
                                 strDispModel,
                                 scaMaxEstimationCyclesEMlike,
                                 scaMaxEstimationCyclesDropModel){
  
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
  
  # 3. vecNormConstExternal
  if(!is.null(vecNormConstExternal)){
    checkNumeric(vecNormConstExternal,"vecNormConstExternal")
    if(!all(colnames(matCounts) %in% names(vecNormConstExternal))){
      stop("ERROR: Not all cells in matCounts are given in vecNormConstExternal")
    }
  }
  
  # 4. scaMaxEstimationCyclesEMlike
  checkNumeric(scaMaxEstimationCyclesEMlike, "scaMaxEstimationCyclesEMlike")
  
  # 5. scaMaxEstimationCyclesDropModel
  checkNumeric(scaMaxEstimationCyclesDropModel, "scaMaxEstimationCyclesDropModel")
  
  # 6. vecFixedAssignments
  if(!is.null(vecNormConstExternal)){
    checkNumeric(vecFixedAssignments, "vecFixedAssignments")
    if(!all(colnames(matCounts) %in% names(vecFixedAssignments))){
      stop("ERROR: Not all cells in matCounts are given in vecFixedAssignments")
    }
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
  
  # 2. Check single dispersion estimation models
  if(!(strDispModel %in% c("constant"))){
    stop(paste0("strDispModel not recognised: ", strDispModel, 
                " Must be one of: constant."))
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
  # Remove all zero or NA genes/cells
  vecidxGenes <- apply(matCounts, 1, function(gene){any(gene>0 & is.finite(gene) & !is.na(gene))})
  vecidxCells <- apply(matCounts, 2, function(cell){any(cell>0 & is.finite(cell) & !is.na(cell))})
  matCountsProc <- matCounts[vecidxGenes,vecidxCells]
  # Keep target normalisation constants
  if(!is.null(vecNormConstExternal)){
    vecNormConstExternalProc <- vecNormConstExternal[colnames(matCountsProc)]
  } else {
    vecNormConstExternalProc <- NULL
  }
  if(!is.null(vecFixedAssignments)){
    vecFixedAssignmentsProc <- vecFixedAssignments[colnames(matCountsProc)]
    # Change to numbers from 1 to k given clusters
    vecGivenMixtures <- unique(vecFixedAssignmentsProc[!is.na(vecFixedAssignmentsProc)])
    vecGivenMixturesProc <- seq(1, length(vecGivenMixtures))
    vecFixedAssignmentsProc[!is.na(vecFixedAssignmentsProc)] <- vecGivenMixturesProc[
      match(vecFixedAssignmentsProc[!is.na(vecFixedAssignmentsProc)],
            vecGivenMixtures)]
  } else {
    vecFixedAssignmentsProc <- NULL
  }
  
  
  # Print summary of processing
  print(paste0(sum(!vecidxGenes), " out of ", length(vecidxGenes), " genes did not contain non-zero observations and are excluded from analysis."))
  print(paste0(sum(!vecidxCells), " out of ", length(vecidxCells), " cells did not contain non-zero observations and are excluded from analysis."))
  
  # Reduce matPiConstPredictors to genes in matCountsProc
  matPiConstPredictorsProc <- matPiConstPredictors[rownames(matCountsProc)]
  
  return(list(matCountsProc=matCountsProc,
              vecFixedAssignmentsProc=vecFixedAssignmentsProc,
              vecNormConstExternalProc=vecNormConstExternalProc,
              matPiConstPredictorsProc=matPiConstPredictorsProc ))
}