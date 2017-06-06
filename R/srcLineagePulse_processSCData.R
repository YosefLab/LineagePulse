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
#' checkNull() Check whether object was supplied (is not NULL).
#' checkNumeric() Checks whether elements are numeric.
#' checkCounts() Checks whether elements are count data.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCounts: (matrix genes x samples)
#' Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param vecNormConstExternal: (numeric vector number of cells) 
#' Model scaling factors, one per cell. These factors will linearly 
#' scale the mean model for evaluation of the loglikelihood. 
#' Must be named according to the column names of matCounts.
#' Supplied by user.
#' @param vecPseudotime: (numerical vector length number of cells)
#' Pseudotime coordinates (1D) of cells: One scalar per cell.
#' Has to be named: Names of elements are cell names.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#' on which ImpulseDE2 is supposed to be run, the full
#' data set is only used for size factor estimation.
#' @param strMuModel: (str) {"constant"}
#' [Default "impulse"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
#' @param strDispModelFull: (str) {"constant"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
#' @param strDispModelRed: (str) {"constant"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' pseudotime in the null model (H0).
#'
#' @return list: (length 3)
#' \itemize{
#'     \item objLP: (LineagePulseObject)
#'     Initialisation of LineagePulseObject.
#'     \item matCountsProcFull: (matrix genes x samples)
#' Processed count data of all cells, unobserved entries are NA.
#'     \item vecPseudotimeProc: (numerical vector length number of cells)
#' Names of elements are cell names.
#' }
#' 
#' @author David Sebastian Fischer
#' 
#' @export

processSCData <- function(
    matCounts,
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
    
    # Check whether object was supplied (is not NULL).
    checkNull <- function(objectInput,strObjectInput){
        if(is.null(objectInput)){
            stop(paste0( "ERROR: ", strObjectInput," was not given as input." ))
        }
    }
    # Checks whether elements are numeric
    checkNumeric <- function(matInput, strMatInput){
        if(any(!is.numeric(matInput))){
            stop(paste0( "ERROR: ", strMatInput, " contains non-numeric elements." ))
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
    
    # 2. dfAnnotation, vecConfounders
    if(!is.null(dfAnnotation)){
        # Check that all cells are mentioned in dfAnnotation
        if(!all(colnames(matCounts) %in% rownames(dfAnnotation))){
            stop(paste0("Not all cells given in matCounts (colnames) are given in dfAnnotation (rownames)."))
        }
        # Check structure
        if(strMuModel=="impulse"){
            checkNull(dfAnnotation$pseudotime,"dfAnnotation$pseudotime")
            checkNumeric(dfAnnotation$pseudotime,"dfAnnotation$pseudotime")
        }
        if(any(rownames(dfAnnotation)!=dfAnnotation$cell)){
            stop(paste0("Cell IDs in rownames(dfAnnotation) are not the same as cell IDs in dfAnnotation$Samples."))
        }
        if(!is.null(vecConfoundersMu)){
            if(!all(vecConfoundersMu %in% colnames(dfAnnotation))){
                stop(paste0("Not all confounders given in vecConfoundersMu given in dfAnnotation (columns)."))
            }
            if(any(is.null(dfAnnotation[,vecConfoundersMu]) | is.na(dfAnnotation[,vecConfoundersMu]))){
                stop(paste0("Supply batch assignments for all cells and all confounders given in vecConfoundersMu"))
            }
        }
    } else {
        stop(paste0("Supply dfAnnotation."))
    }
    
    # 3. matPiConstPredictors
    if(!is.null(matPiConstPredictors)){
        checkNumeric(as.matrix(matPiConstPredictors),"matPiConstPredictors")
        if(!is.null(rownames(matCounts))){
            if(any(!rownames(matCounts) %in% rownames(matPiConstPredictors))){
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
    
    # 4. vecNormConstExternal
    if(!is.null(vecNormConstExternal)){
        checkNumeric(vecNormConstExternal,"vecNormConstExternal")
        if(!all(colnames(matCounts) %in% names(vecNormConstExternal))){
            stop("ERROR IN INPUT DATA CHECK: ",
                 "Not all cells in matCounts are given in vecNormConstExternal")
        }
    }
    
    # 5. scaMaxEstimationCycles
    checkNumeric(scaMaxEstimationCycles, "scaMaxEstimationCycles")
    
    # 6. booleans
    checkLogical(boolVerbose, "boolVerbose")
    checkLogical(boolSuperVerbose, "boolSuperVerbose")
    
    # (II) Check settings
    # 1. Check general mean estimation settings
    
    # 2. Check single mean-/dispersion estimation models
    
    # 3. Check co-estimation models
    # Note: These functions are implemented separately from
    # single mean and dispersion estimation.
    if(!(strMuModel %in% c("groups","impulse","splines","constant"))){
        stop(paste0("strMuModel not recognised: ", strMuModel, 
                    " Must be one of: groups, splines, impulse, constant."))
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
    # Take out cells with NA pseudotime coordinate
    if(strMuModel=="impulse"){
        vecidxPTnotNA <- !is.na(dfAnnotation$pseudotime)
        dfAnnotationProc <- dfAnnotation[vecidxPTnotNA,]
        vecidxPTsort <- sort(dfAnnotationProc$pseudotime,
                             decreasing=FALSE, index.return=TRUE)$ix
        dfAnnotationProc <- dfAnnotation[vecidxPTsort,]
        matCountsProc <- matCounts[,dfAnnotationProc$cell]
    } else {
        dfAnnotationProc <- dfAnnotation
        matCountsProc <- matCounts
    }
    # Remove all zero or NA genes/cells
    vecidxGenes <- apply(matCountsProc, 1, function(gene){any(gene>0 & is.finite(gene) & !is.na(gene))})
    vecidxCells <- apply(matCountsProc, 2, function(cell){any(cell>0 & is.finite(cell) & !is.na(cell))})
    dfAnnotationProc <- dfAnnotationProc[vecidxCells,]
    matCountsProc <- matCountsProc[vecidxGenes,vecidxCells]
    # Keep target normalisation constants
    if(!is.null(vecNormConstExternal)) vecNormConstExternalProc <- vecNormConstExternal[colnames(matCountsProc)]
    else vecNormConstExternalProc <- NULL
    
    # Print summary of processing
    strMessage <- paste0("LineagePulse for count data: v0.99")#, packageDescription("LineagePulse", fields = "Version"))
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    strMessage <- paste0("--- Data preprocessing")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    if(strMuModel=="impulse"){
        strMessage <- paste0("# ", sum(!vecidxPTnotNA), " out of ", length(vecidxPTnotNA), 
                             " cells did not have a pseudotime coordinate and were excluded.")
        strReport <- paste0(strReport, strMessage, "\n")
        if(boolVerbose) print(strMessage)
    }
    
    strMessage <- paste0("# ", sum(!vecidxGenes), " out of ", length(vecidxGenes), 
                         " genes did not contain non-zero observations and are excluded from analysis.")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    strMessage <- paste0("# ", sum(!vecidxCells), " out of ", length(vecidxCells), 
                         " cells did not contain non-zero observations and are excluded from analysis.")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    # Reduce matPiConstPredictors to genes in matCountsProc
    matPiConstPredictorsProc <- matPiConstPredictors[rownames(matCountsProc),,drop=FALSE]
    
    objLP <- new(
        'LineagePulseObject',
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
        matCountsProc       = sparseMatrix(
            i = (which(matCountsProc>0)-1) %% dim(matCountsProc)[1],
            j = (which(matCountsProc>0)-1) %/% dim(matCountsProc)[1], 
            x = matCountsProc[which(matCountsProc>0)], 
            dims = dim(matCountsProc),
            dimnames = dimnames(matCountsProc),
            symmetric = FALSE,
            index1 = FALSE),
        matWeights          = NULL,
        scaDFSplinesDisp    = scaDFSplinesDisp,
        scaDFSplinesMu      = scaDFSplinesMu,
        strReport           = strReport,
        vecAllGenes         = rownames(matCounts),
        vecConfoundersMu    = vecConfoundersMu,
        vecConfoundersDisp  = vecConfoundersDisp,
        boolFixedPopulations= FALSE,
        vecNCentroidsPerPop = NULL,
        vecH0Pop            = NULL,
        vecNormConst        = NULL,
        strVersion          = "0.99")#packageDescription("LineagePulse", fields = "Version"))
    
    return(list(objLP=objLP,
                vecNormConstExternalProc=vecNormConstExternalProc,
                matPiConstPredictorsProc=matPiConstPredictorsProc ))
}