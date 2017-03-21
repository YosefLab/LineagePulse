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
																 dfAnnotation,
																 vecConfounders,
																 boolFixedPopulations,
																 vecNCentroidsPerPop,
																 vecH0Pop,
                                 scaNMixtures,
                                 matPiConstPredictors,
                                 vecNormConstExternal,
                                 strDispModel,
                                 scaMaxEstimationCyclesEMlike,
                                 scaMaxEstimationCyclesDropModel,
																 boolVerbose,
																 boolSuperVerbose,
																 STR_VERSION){
  
  # Check whether object was supplied (is not NULL).
  checkNull <- function(objectInput,strObjectInput){
    if(is.null(objectInput)){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strObjectInput," was not given as input." ))
    }
  }
  # Checks whether elements are numeric
  checkNumeric <- function(matInput, strMatInput){
    if(any(!is.numeric(matInput))){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strMatInput, " contains non-numeric elements. Requires count data." ))
    }
  }
  # Checks whether elements are count data: non-negative integer finite numeric elements.
  # Note that NA are allowed.
  checkCounts <- function(matInput, strMatInput){
    checkNumeric(matInput, strMatInput)
    if(any(matInput[!is.na(matInput)] %% 1 != 0)){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strMatInput, " contains non-integer elements. Requires count data." ))
    }
    if(any(!is.finite(matInput[!is.na(matInput)]))){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strMatInput, " contains infinite elements. Requires count data." ))
    }
    if(any(matInput[!is.na(matInput)]<0)){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strMatInput, " contains negative elements. Requires count data." ))
    }
  }
  # Checks whether elements are logical
  checkLogical <- function(boolElement, strBoolElement){
    if(!is.logical(boolElement)){
      stop(paste0( "ERROR IN INPUT DATA CHECK: ", strBoolElement, " must be logical (TRUE or FALSE)." ))
    }
  }
  
  # (I) Check input data
  # 1. matCounts
  checkNull(matCounts,"matCounts")
  checkCounts(matCounts,"matCounts")
  
  # 2. dfAnnotation, vecConfounders, boolFixedPopulations, vecNCentroidsPerPop
  if(!is.null(dfAnnotation)){
  	# Check that all cells are mentioned in dfAnnotation
  	if(!all(colnames(matCounts) %in% rownames(dfAnnotation))){
  		stop(paste0("ERROR IN INPUT DATA CHECK: Not all cells given in matCounts (colnames) are given in dfAnnotation (rownames)."))
  	}
  	# Check structure
  	if(any(rownames(dfAnnotation)!=dfAnnotation$cell)){
  		stop(paste0("ERROR IN INPUT DATA CHECK: ",
  								"Cell IDs in rownames(dfAnnotation) are not the same as cell IDs in dfAnnotation$Samples."))
  	}
  	if(!is.null(vecConfounders)){
  		if(!all(vecConfounders %in% colnames(dfAnnotation))){
  			stop(paste0("ERROR IN INPUT DATA CHECK: ",
  									"Not all confounders given in vecConfounders given in dfAnnotation (columns)."))
  		}
  		if(any(is.null(dfAnnotation[,vecConfounders]) | is.na(dfAnnotation[,vecConfounders]))){
  			stop(paste0("ERROR IN INPUT DATA CHECK: ",
  									"Supply batch assignments for all cells and all confounders given in vecConfounders"))
  		}
  	}
  	if(boolFixedPopulations){
  		if(!"populations" %in% colnames(dfAnnotation)){
  			stop(paste0("ERROR IN INPUT DATA CHECK: ",
  									"population does not occur in dfAnnotation (columns)."))
  		} else {
  			if(!is.null(scaNMixtures)){
  				if(scaNMixtures < length(unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)]))){
  					stop(paste0("ERROR IN INPUT DATA CHECK: ",
  											"Number of populations in dfAnnotation$populations (",
  											length(unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)])),
  											"):", unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)]),
  											" must not exceed maximum number of populations given by scaNMixtures=",
  											scaNMixtures))
  				}
  			}
  		}
  	  if(!is.null(vecNCentroidsPerPop)){
  	    vecFixedPop <- unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)])
  	    # Check names
  	    if(!all(names(vecNCentroidsPerPop) %in% vecFixedPop)){
  	      stop(paste0("ERROR IN INPUT DATA CHECK: ",
  	                  "Not all names(vecNCentroidsPerPop) are in dfAnnotation$populations."))
  	    }
  	    if(!all(vecFixedPop %in% names(vecNCentroidsPerPop))){
  	      stop(paste0("ERROR IN INPUT DATA CHECK: ",
  	                  "Not all dfAnnotation$populations are in names(vecNCentroidsPerPop)."))
  	    }
  	    # Check content
  	    checkNumeric(vecNCentroidsPerPop)
  	    if(sum(vecNCentroidsPerPop)==scaNMixtures){
  	      warning(paste0("WARNING IN INPUT DATA CHECK: ",
  	                     "Sum of centroids per fixed populations",
  	                     "(vecNCentroidsPerPop) is equal ",
  	                     "to total number of centroids (scaNMixtures). ",
  	                     "There are no free centroids."))
  	    }
  	    if(sum(vecNCentroidsPerPop)>scaNMixtures){
  	      stop(paste0("WARNING IN INPUT DATA CHECK: ",
  	                  "Sum of centroids per fixed populations",
  	                  "(vecNCentroidsPerPop) is larger ",
  	                  "than total number of centroids (scaNMixtures)."))
  	    }
  	    if(!is.null(vecH0Pop)){
  	      if(!all(vecH0Pop %in% vecFixedPop)){
  	        stop(paste0("ERROR IN INPUT DATA CHECK: ",
  	                    "Not all vecH0Pop are in dfAnnotation$populations."))
  	      }
  	    }
  	  }
  	} else {
  	  if(!is.null(vecNCentroidsPerPop)){
  	    stop(paste0("ERROR IN INPUT DATA CHECK: ",
  	                "Do not supply vecNCentroidsPerPop if boolFixedPopulations=FALSE."))
  	  } 
  	}
  } else {
  	if(!is.null(vecConfounders)){
  		stop(paste0("ERROR IN INPUT DATA CHECK: ",
  		            "vecConfounders supplied but no dfAnnotation which carries batch structure."))
  	}
  	if(boolFixedPopulations){
  		stop(paste0("ERROR IN INPUT DATA CHECK: ",
  		            "boolFixedPopulations=TRUE but no dfAnnotation which carries population assignments."))
  	}
    if(!is.null(vecH0Pop)){
      stop(paste0("ERROR IN INPUT DATA CHECK: ",
                  "vecH0Pop supplied but no dfAnnotation which carries population assignments."))
    }
  }
  
  # 3. matPiConstPredictors
  if(!is.null(matPiConstPredictors)){
    checkNumeric(as.matrix(matPiConstPredictors),"matPiConstPredictors")
    if(!is.null(rownames(matCounts))){
      if(any(!rownames(matCounts) %in% rownames(matPiConstPredictors))){
        stop(paste0("ERROR IN INPUT DATA CHECK: Some genes named in rows of matCounts do not ",
                    "occur in rows of matPiConstPredictors."))   
      }
    } else {
      if(!is.null(rownames(matPiConstPredictors))){
        stop(paste0("ERROR IN INPUT DATA CHECK: Named genes in matPiConstPredictors",
                    " but not in matCounts."))
      }
    }
  }
  
  # 4. vecNormConstExternal
  if(!is.null(vecNormConstExternal)){
    checkNumeric(vecNormConstExternal,"vecNormConstExternal")
    if(!all(colnames(matCounts) %in% names(vecNormConstExternal))){
      stop("ERROR IN INPUT DATA CHECK: ",
      		 "Not all cells in matCounts are given in vecNormConstExternal")
    }
  }
  
  # 5. scaMaxEstimationCyclesEMlike
  checkNumeric(scaMaxEstimationCyclesEMlike, "scaMaxEstimationCyclesEMlike")
  
  # 6. scaMaxEstimationCyclesDropModel
  checkNumeric(scaMaxEstimationCyclesDropModel, "scaMaxEstimationCyclesDropModel")
  
  # 7. scaNMixtures
  checkNull(scaNMixtures, "scaNMixtures")
  checkNumeric(scaNMixtures, "scaNMixtures")
  
  # 8. booleans
  checkLogical(boolFixedPopulations, "boolFixedPopulations")
  checkLogical(boolVerbose, "boolVerbose")
  checkLogical(boolSuperVerbose, "boolSuperVerbose")
  
  # (II) Check settings
  # Check single dispersion estimation models
  if(!(strDispModel %in% c("constant"))){
    stop(paste0("ERROR IN INPUT DATA CHECK: ",
    						"strDispModel not recognised: ", strDispModel, 
                " Must be one of: constant."))
  }
  
  # (III) Process data
  strReport <- ""
  # Convert from data frame to matrix
  if(is.list(matCounts)){
    matCounts <- data.matrix(matCounts)
  }
  # Name genes if names not given
  if(is.null(rownames(matCounts))){
    rownames(matCounts) <- paste0("Gene_", seq(1,nrow(matCounts)))
    rownames(matPiConstPredictors) <- rownames(matCounts)
  }
  matCountsProc <- matCounts[,dfAnnotation$cell]
  # Remove all zero or NA genes/cells
  vecidxGenes <- apply(matCountsProc, 1, function(gene) any(gene>0 & is.finite(gene) & !is.na(gene)) )
  vecidxCells <- apply(matCountsProc, 2, function(cell) any(cell>0 & is.finite(cell) & !is.na(cell)) )
  dfAnnotationProc <- dfAnnotation[vecidxCells,]
  matCountsProc <- matCountsProc[vecidxGenes,vecidxCells]
  # Keep target normalisation constants
  if(!is.null(vecNormConstExternal)){
    vecNormConstExternalProc <- vecNormConstExternal[colnames(matCountsProc)]
  } else {
    vecNormConstExternalProc <- NULL
  }
  if(boolFixedPopulations){
    vecFixedAssignmentsUnique <-  unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)])
    dfAnnotation$populations <- match(dfAnnotation$populations, vecFixedAssignmentsUnique)
    if(is.null(vecNCentroidsPerPop)){
      # Create vecNCentroidsPerPop if boolFixedPopulations and this was not given
      # with one centroid per population.
      vecFixedPop <- unique(dfAnnotation$populations[!is.na(dfAnnotation$populations)])
      vecNCentroidsPerPop <- rep(1, length(vecFixedPop))
      names(vecNCentroidsPerPop) <- vecFixedPop
    } else {
      # Make sure order is right
      vecNCentroidsPerPop <- vecNCentroidsPerPop[vecFixedAssignmentsUnique]
    }
  }
  
  # Print summary of processing
  strMessage <- paste0("LineagePulse ",STR_VERSION)
  strReport <- paste0(strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  strMessage <- paste0("--- Data preprocessing")
  strReport <- paste0(strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  strMessage <- paste0("# ", sum(!vecidxGenes), " out of ", length(vecidxGenes), 
                       " genes did not contain non-zero observations and are excluded from analysis.")
  strReport <- paste0(strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  strMessage <- paste0("# ", sum(!vecidxCells), " out of ", length(vecidxCells), 
                       " cells did not contain non-zero observations and are excluded from analysis.")
  strReport <- paste0(strReport, strMessage, "\n")
  if(boolVerbose) print(strMessage)
  
  # Reduce matPiConstPredictors to genes in matCountsProc
  matPiConstPredictorsProc <- matPiConstPredictors[rownames(matCountsProc),]
  
  objectLineagePulse <- new('LineagePulseObject',
                            dfAnnotationProc    = dfAnnotationProc,
                            dfResults           = NULL,
                            lsDispModelH0       = NULL,
                            lsDispModelH1       = NULL,
                            lsDispModelConst    = NULL,
                            lsDropModel         = NULL,
                            lsFitZINBReporters  = NULL,
                            lsMuModelH1         = NULL,
                            lsMuModelH0         = NULL,
                            lsMuModelConst      = NULL,
                            matCountsProc       = matCountsProc,
                            matWeights          = NULL,
                            strReport           = strReport,
                            vecAllGenes         = rownames(matCounts),
  													vecConfounders      = vecConfounders,
  													boolFixedPopulations= boolFixedPopulations,
  													vecNCentroidsPerPop = vecNCentroidsPerPop,
  													vecH0Pop            = vecH0Pop,
                            vecNormConst        = NULL,
  													strVersion          = STR_VERSION)
  
  return(list(objectLineagePulse=objectLineagePulse,
              vecNormConstExternalProc=vecNormConstExternalProc,
              matPiConstPredictorsProc=matPiConstPredictorsProc ))
}