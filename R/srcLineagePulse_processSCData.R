#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Process single cell data    ++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Prepare single cell data for analysis
#' 
#' Check that input is correctly supplied and formatted. Then process input
#' data for analysis. This function catches downstream errors arising from 
#' incorrect input parameter combinations and explains the input error
#' to the user. Input passing this function should not lead to errors
#' arising from incorrect input.
#' Helper functions:
#' checkNull() Check whether object was supplied (is not NULL).
#' checkNumeric() Checks whether elements are numeric.
#' checkCounts() Checks whether elements are count data.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param counts (matrix genes x samples)
#' (matrix genes x cells (sparseMatrix or standard) or file)
#' Matrix: Count data of all cells, unobserved entries are NA.
#' file: .mtx file from which count matrix is to be read.
#' @param dfAnnotation (data frame cells x meta characteristics)
#' Annotation table which contains meta data on cells.
#' @param vecConfoundersMu (vector of strings number of confounders on  mean)
#' [Default NULL] Confounders to correct for in mu batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param vecConfoundersDisp 
#' (vector of strings number of confounders on  dispersion)
#' [Default NULL] Confounders to correct for in dispersion batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param matPiConstPredictors (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param vecNormConstExternal (numeric vector number of cells) 
#' Model scaling factors, one per cell. These factors will linearly 
#' scale the mean model for evaluation of the loglikelihood. 
#' Must be named according to the column names of matCounts.
#' Supplied by user.
#' @param strDispModelFull (str) {"constant"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' population structure in the alternative model (H1).
#' @param strDispModelRed (str) {"constant"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' population structure in the null model (H0).
#' @param strMuModel (str) {"constant"}
#' [Default "impulse"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' population structure in the alternative model (H1).
#' @param scaDFSplinesDisp (sca) [Default 3] 
#' If strDispModelFull=="splines" or strDispModelRed=="splines", 
#' the degrees of freedom of the natural
#' cubic spline to be used as a dispersion parameter model.
#' @param scaDFSplinesMu (sca) [Default 3] 
#' If strMuModel=="splines", the degrees of freedom of the natural
#' cubic spline to be used as a mean parameter model.
#' @param scaMaxEstimationCycles (integer) [Default 20] Maximum number 
#' of estimation cycles performed in fitZINB(). One cycle
#' contain one estimation of of each parameter of the 
#' zero-inflated negative binomial model as coordinate ascent.
#' @param boolVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation in high detail with local 
#' convergence flags and step-by-step loglikelihood computation.
#'
#' @return list (length 3)
#' \itemize{
#'     \item objLP (LineagePulseObject)
#'     Initialisation of LineagePulseObject.
#'     \item matCountsProcFull (matrix genes x samples)
#' Processed count data of all cells, unobserved entries are NA.
#'     \item matPiConstPredictorsProc 
#' (numeric matrix genes x number of constant gene-wise 
#' drop-out predictors) Processed version of matPiConstPredictors.
#' }
#' 
#' @author David Sebastian Fischer
processSCData <- function(
    counts,
    dfAnnotation,
    vecConfoundersDisp,
    vecConfoundersMu,
    matPiConstPredictors,
    vecNormConstExternal,
    strDispModelFull,
    strDispModelRed,
    strMuModel,
    scaDFSplinesDisp,
    scaDFSplinesMu,
    scaMaxEstimationCycles,
    boolVerbose,
    boolSuperVerbose){
    
    # read count object from file if supplied as string
    if(is(counts, "character")) {
        # check that format is .mtx
        strFormat <- unlist(strsplit(x = counts, split = "[.]"))
        strFormat <- strFormat[length(strFormat)]
        if(strFormat!="mtx"){
            stop("ERROR IN INPUT DATA: ",
                 "Input matrix file is not .mtx. ",
                 "Supply counts as .mtx file or as R matrix object.")
        }
        counts <- readMM(counts)
    }
    
    # Check whether object was supplied (is not NULL).
    checkNull <- function(objectInput,strObjectInput){
        if(is.null(objectInput)){
            stop( "ERROR: ", strObjectInput,
                  " was not given as input." )
        }
    }
    # Checks whether elements are numeric
    checkNumeric <- function(matInput, strMatInput){
        if(any(!is.numeric(matInput))){
            stop( "ERROR: ", strMatInput, 
                  " contains non-numeric elements." )
        }
    }
    # Checks whether elements are count data: 
    # non-negative integer finite numeric elements.
    # Note that NA are allowed.
    checkCounts <- function(matInput, strMatInput){
        checkNumeric(matInput, strMatInput)
        if(any(matInput[!is.na(matInput)] %% 1 != 0)){
            stop( "ERROR: ", strMatInput, 
                  " contains non-integer elements.",
                  " Requires count data." )
        }
        if(any(!is.finite(matInput[!is.na(matInput)]))){
            stop( "ERROR: ", strMatInput, 
                  " contains infinite elements.",
                  " Requires count data." )
        }
        if(any(matInput[!is.na(matInput)]<0)){
            stop( "ERROR: ", strMatInput, 
                  " contains negative elements.",
                  " Requires count data." )
        }
    }
    # Checks whether elements are logical
    checkLogical <- function(boolElement, strBoolElement){
        if(!is.logical(boolElement)){
            stop( "ERROR IN INPUT DATA CHECK: ", 
                  strBoolElement, 
                  " must be logical (TRUE or FALSE)." )
        }
    }
    
    # (I) Check input data
    # 1. counts
    checkNull(counts,"counts")
    checkCounts(counts,"counts")
    
    # 2. dfAnnotation, vecConfounders
    if(!is.null(dfAnnotation)){
        # Check that all cells are mentioned in dfAnnotation
        if(!all(colnames(counts) %in% rownames(dfAnnotation))){
            stop("Not all cells given in counts (colnames)",
                 " are given in dfAnnotation (rownames).")
        }
        # Check structure
        if(strMuModel %in% c("splines", "impulse")){
            checkNull(dfAnnotation$continuous,"dfAnnotation$continuous")
            checkNumeric(dfAnnotation$continuous,"dfAnnotation$continuous")
        }
        if(any(rownames(dfAnnotation)!=dfAnnotation$cell)){
            stop("Cell IDs in rownames(dfAnnotation)",
                 " are not the same as cell IDs",
                 " in dfAnnotation$Samples.")
        }
        if(!is.null(vecConfoundersMu)){
            if(!all(vecConfoundersMu %in% colnames(dfAnnotation))){
                stop("Not all confounders given in",
                     " vecConfoundersMu given in",
                     " dfAnnotation (columns).")
            }
            if(any(is.null(dfAnnotation[,vecConfoundersMu]) | 
                   is.na(dfAnnotation[,vecConfoundersMu]))){
                stop("Supply batch assignments for",
                     " all cells and all confounders given",
                     " in vecConfoundersMu")
            }
        }
    } else {
        stop("Supply dfAnnotation.")
    }
    
    # 3. matPiConstPredictors
    if(!is.null(matPiConstPredictors)){
        checkNumeric(as.matrix(matPiConstPredictors),"matPiConstPredictors")
        if(!is.null(rownames(counts))){
            if(any(!rownames(counts) %in% rownames(matPiConstPredictors))){
                stop("ERROR: Some genes named in rows of counts do not ",
                     "occur in rows of matPiConstPredictors.")
            }
        } else {
            if(!is.null(rownames(matPiConstPredictors))){
                stop("ERROR: Named genes in matPiConstPredictors",
                     " but not in counts.")
            }
        }
    }
    
    # 4. vecNormConstExternal
    if(!is.null(vecNormConstExternal)){
        checkNumeric(vecNormConstExternal,"vecNormConstExternal")
        if(!all(colnames(counts) %in% names(vecNormConstExternal))){
            stop("ERROR IN INPUT DATA CHECK: ",
                 "Not all cells in counts are given in vecNormConstExternal")
        }
    }
    
    # 5. scaMaxEstimationCycles
    checkNumeric(scaMaxEstimationCycles, "scaMaxEstimationCycles")
    
    # 6. booleans
    checkLogical(boolVerbose, "boolVerbose")
    checkLogical(boolSuperVerbose, "boolSuperVerbose")
    
    # (II) Check settings
    # 1. Check gene-wise models
    if(!(strMuModel %in% c("groups","impulse","splines","constant"))){
        stop("strMuModel not recognised: ", strMuModel, 
             " Must be one of: groups, splines, impulse, constant.")
    }
    if(!(strDispModelFull %in% c("groups","splines","constant"))){
        stop("strDispModelFull not recognised: ", strDispModelFull, 
             " Must be one of: groups, splines, constant.")
    }
    if(!(strDispModelRed %in% c("groups","splines","constant"))){
        stop("strDispModelRed not recognised: ", strDispModelRed, 
             " Must be one of: groups, splines, constant.")
    }
    if(!(strDispModelRed %in% c("groups","splines","constant"))){
        stop("strDispModelRed not recognised: ", strDispModelRed, 
             " Must be one of: groups, splines, constant.")
    }
    
    # (III) Process data
    strReport <- ""
    # Convert from data frame to matrix
    if(is.list(counts)){
        counts <- data.matrix(counts)
    }
    # sort count matrix according to dfAnnotation
    if(!all(colnames(counts) == dfAnnotation$cell)) {
        counts <- counts[,as.double(match(colnames(counts), 
                                          dfAnnotation$cell))]
    }
    # make sparse if is not sparse
    if(!is(object = counts, class2 = "sparseMatrix")){
        counts <- sparseMatrix(
            i = (which(counts>0)-1) %% dim(counts)[1],
            j = (which(counts>0)-1) %/% dim(counts)[1], 
            x = counts[which(counts>0)], 
            dims = dim(counts),
            dimnames = dimnames(counts),
            symmetric = FALSE,
            index1 = FALSE)
    }
    # Name genes if names not given
    if(is.null(rownames(counts))){
        rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
        rownames(matPiConstPredictors) <- rownames(counts)
    }
    # Take out cells with NA continuous covariate
    if(strMuModel %in% c("splines", "impulse")){
        vecboolPTnotNA <- !is.na(dfAnnotation$continuous)
        dfAnnotationProc <- dfAnnotation[vecboolPTnotNA,]
        counts <- counts[,as.double(which(vecboolPTnotNA))]
    } else {
        dfAnnotationProc <- dfAnnotation
    }
    # Remove all zero or NA genes/cells
    vecAllGenes <- rownames(counts)
    vecboolGenes <- Matrix::rowSums(counts, na.rm = TRUE) > 0
    vecboolCells <- Matrix::colSums(counts, na.rm = TRUE) > 0
    dfAnnotationProc <- dfAnnotationProc[which(vecboolCells),]
    counts <- counts[as.double(which(vecboolGenes)), 
                     as.double(which(vecboolCells))]
    # Keep target normalisation constants
    if(!is.null(vecNormConstExternal)) {
        vecNormConstExternalProc <- vecNormConstExternal[colnames(counts)]
    } else {
        vecNormConstExternalProc <- NULL
    }
    
    # Print summary of processing
    strMessage <- paste0("LineagePulse for count data: v", 
                         packageDescription("LineagePulse", 
                                            fields = "Version"))
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    strMessage <- paste0("--- Data preprocessing")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    if(strMuModel %in% c("splines","impulse")){
        strMessage <- paste0("# ", sum(!vecboolPTnotNA), " out of ", 
                             length(vecboolPTnotNA), 
                             " cells did not have a continuous covariate",
                             " and were excluded.")
        strReport <- paste0(strReport, strMessage, "\n")
        if(boolVerbose) message(strMessage)
    }
    
    strMessage <- paste0("# ", sum(!vecboolGenes), " out of ", 
                         length(vecboolGenes), 
                         " genes did not contain non-zero observations",
                         " and are excluded from analysis.")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    strMessage <- paste0("# ", sum(!vecboolCells), " out of ", 
                         length(vecboolCells), 
                         " cells did not contain non-zero observations",
                         " and are excluded from analysis.")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    # Reduce matPiConstPredictors to genes in counts
    matPiConstPredictorsProc <- 
        matPiConstPredictors[rownames(counts),,drop=FALSE]
    
    objLP <- new(
        'LineagePulseObject',
        dfAnnotationProc    = dfAnnotationProc,
        dfResults           = NULL,
        lsDispModelH0       = NULL,
        lsDispModelH1       = NULL,
        lsDispModelConst    = NULL,
        lsDropModel         = NULL,
        lsFitConvergence    = NULL,
        lsMuModelH1         = NULL,
        lsMuModelH0         = NULL,
        lsMuModelConst      = NULL,
        matCountsProc       = counts,
        matWeights          = NULL,
        scaDFSplinesDisp    = scaDFSplinesDisp,
        scaDFSplinesMu      = scaDFSplinesMu,
        strReport           = strReport,
        vecAllGenes         = vecAllGenes,
        vecConfoundersMu    = vecConfoundersMu,
        vecConfoundersDisp  = vecConfoundersDisp,
        boolFixedPopulations= FALSE,
        vecNCentroidsPerPop = NULL,
        vecH0Pop            = NULL,
        vecNormConst        = NULL,
        strVersion          = packageDescription("LineagePulse", 
                                                 fields = "Version"))
    
    return(list(objLP=objLP,
                vecNormConstExternalProc=vecNormConstExternalProc,
                matPiConstPredictorsProc=matPiConstPredictorsProc ))
}