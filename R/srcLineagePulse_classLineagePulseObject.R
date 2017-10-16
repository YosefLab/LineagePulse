###############################################################################
#################     LineagePulse output container class     #################
###############################################################################

### 1. Define output container class

# Define class unions for slots
setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('matrixORNULL', members=c('matrix', 'NULL'))
setClassUnion('characterORNULL', members=c('character', 'NULL'))
setClassUnion('listORNULL', members=c('list', 'NULL'))
setClassUnion('data.frameORNULL', members=c('data.frame', 'NULL'))

#' Container class for LineagePulse output
#' 
#' LineagePulse output and intermediate results such as model fits.
#' 
#' @slot dfAnnotationProc  (data frame cells x meta characteristics)
#' Annotation table which contains meta data on cells.
#' @slot dfResults (data frame samples x reported characteristics) 
#' Summary of fitting procedure and 
#' differential expression results for each gene.
#' @slot lsDispModelH1 (list)
#' Object containing description of gene-wise dispersion parameter models of H1.
#' @slot lsDispModelH0 (list)
#' Object containing description of gene-wise dispersion parameter models of H0.
#' @slot lsDispModelConst (list)
#' Object containing description of gene-wise dispersion parameter models of
#' constant model (necessary if H0 is not constant (mixture models)).
#' @slot lsDropModel (list)
#' @slot lsMuModelH0 (list)
#' Object containing description of gene-wise mean parameter models of H0.
#' @slot lsMuModelH1 (list)
#' Object containing description of gene-wise mean parameter models of H1.
#' @slot lsMuModelConst (list)
#' Object containing description of gene-wise means parameter models of
#' constant model (necessary if H0 is not constant (mixture models)).
#' @slot lsFitZINBReporters (list)
#' Fitting reporter summary.
#' @slot matCountsProc (count matrix genes x cells)
#' Sparse matrix of counts, non-observed values are NA.
#' @slot matWeights (numeric matrix cells x mixtures)
#' Assignments of cells to mixtures.
#' @slot scaDFSplinesDisp (scalar) 
#' If strDispModelFull=="splines" or strDispModelRed=="splines", 
#' the degrees of freedom of the natural
#' cubic spline to be used as a dispersion parameter model.
#' @slot scaDFSplinesMu (scalar) 
#' If strMuModel=="splines", the degrees of freedom of the natural
#' cubic spline to be used as a mean parameter model.
#' @slot strReport (str) LineagePulse stdout report (log).
#' @slot vecAllGenes  (vector of strings) 
#' All genes originally supplied, including all zero genes.
#' @slot vecConfoundersDisp (vector of strings number of confounders on  dispersion)
#' Confounders to correct for in dispersion batch
#' correction model, are subset of column names of
#' dfAnnotation which describe condounding variables.
#' @slot vecConfoundersMu (vector of strings number of confounders on  mean)
#' Confounders to correct for in mu batch
#' correction model, are subset of column names of
#' dfAnnotation which describe condounding variables.
#' @slot scaOmega (scalar)
#' Regularisation parameter on entropy penalty on mixture weights.
#' Needed for mixture model mode only.
#' @slot boolFixedPopulations
#' Whether a population of cells is to be fixed to be 
#' assigned to a subset of centroids (RSA scenario).
#' Needed for mixture model mode only.
#' @slot vecNCentroidsPerPop
#' Number of centroids per population (RSA scenario).
#' Needed for mixture model mode only.
#' @slot vecH0Pop (integer vector)
#' Needed for mixture model mode only.
#' @slot vecNormConst (numeric vector length number of cells)
#' Normalisation constants (library size factors) for each cell.
#' @slot strVersion (str) Version of LineagePulse that was run.
#' 
#' @name LineagePulseObject-class
#'   
#' @author David Sebastian Fischer
setClass(
  'LineagePulseObject',
  slots = c(
    dfAnnotationProc    = "data.frameORNULL",
    dfResults           = "data.frameORNULL",
    lsDispModelH0       = "listORNULL",
    lsDispModelH1       = "listORNULL",
    lsDispModelConst    = "listORNULL",
    lsDropModel         = "listORNULL",
    lsMuModelH0         = "listORNULL",
    lsMuModelH1         = "listORNULL",
    lsMuModelConst      = "listORNULL",
    lsFitZINBReporters  = "listORNULL",
    matCountsProc       = "dgCMatrix",
    matWeights          = "matrixORNULL",
    scaDFSplinesDisp    = "numericORNULL",
    scaDFSplinesMu      = "numericORNULL",
    strReport           = "characterORNULL",
    vecAllGenes         = "characterORNULL",
    vecConfoundersDisp  = "characterORNULL",
    vecConfoundersMu    = "characterORNULL",
    scaOmega            = "numericORNULL",
    boolFixedPopulations= "logical",
    vecNCentroidsPerPop = "numericORNULL",
    vecH0Pop            = "characterORNULL",
    vecNormConst        = "numericORNULL",
    strVersion          = "character")
)

### 2. Enable accession of private elements via functions
### which carry the same name as the element.

#' LineagePulseObject accession methods
#' 
#' Get internal data of LineagePulse output object.
#' 
#' @param objLP (LineagePulse-Object)  A LineagePulse output object to extract
#' object from.
#' 
#' @return The internal data object specified by the function.
#' 
#' @name get_accessors
#' @rdname get_accessors
#' @aliases
#' get_lsMuModelH0
#' get_lsMuModelH1
#' get_lsMuModelConst
#' get_lsDispModelH0
#' get_lsDispModelH1
#' get_lsDispModelConst
#' get_lsDropModel
#' get_matCountsProc
#' get_matWeights
#' get_scaDFSplinesDisp
#' get_scaDFSplinesMu
#' get_strReport
#' get_vecAllGenes
#' get_vecConfoundersDisp
#' get_vecConfoundersMu
#' get_scaOmega
#' get_boolFixedPopulation
#' get_vecH0Pop
#' get_vecNormConst
#' get_strVersion
#' 
#' @examples    
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 10,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' # get hidden objects within LineagePulse object
#' dfAnnotationProc <- get_dfAnnotationProc(objLP)
#' lsMuModelH0 <- get_lsMuModelH0(objLP)
#' lsMuModelH1 <- get_lsMuModelH1(objLP)
#' lsMuModelConst <- get_lsMuModelConst(objLP)
#' lsDispModelH0 <- get_lsDispModelH0(objLP)
#' lsDispModelH1 <- get_lsDispModelH1(objLP)
#' lsDispModelConst <- get_lsDispModelConst(objLP)
#' lsDropModel <- get_lsDropModel(objLP)
#' matCountDataProc <- get_matCountsProc(objLP)
#' matWeights <- get_matWeights(objLP)
#' scaDFSplinesDisp <- get_scaDFSplinesDisp(objLP) 
#' scaDFSplinesMu <- get_scaDFSplinesMu(objLP) 
#' strReport <- get_strReport(objLP) 
#' vecAllGenes <- get_vecAllGenes(objLP) 
#' vecConfoundersDisp <- get_vecConfoundersDisp(objLP) 
#' vecConfoundersMu <- get_vecConfoundersMu(objLP) 
#' scaOmega <- get_scaOmega(objLP) 
#' boolFixedPopulations <- get_boolFixedPopulations(objLP)
#' vecH0Pop <- get_vecH0Pop(objLP) 
#' vecNormConst <- get_vecNormConst(objLP) 
#' strVersion <- get_strVersion(objLP)
#' 
#' @author David Sebastian Fischer
NULL

#' @rdname get_accessors
#' @export
get_dfAnnotationProc <- function(objLP) 
    return(objLP@dfAnnotationProc)

#' @rdname get_accessors
#' @export
get_lsMuModelH0 <- function(objLP) 
    return(objLP@lsMuModelH0)

#' @rdname get_accessors
#' @export
get_lsMuModelH1 <- function(objLP) 
    return(objLP@lsMuModelH1)

#' @rdname get_accessors
#' @export
get_lsMuModelConst <- function(objLP) 
    return(objLP@lsMuModelConst)

#' @rdname get_accessors
#' @export
get_lsDispModelH0 <- function(objLP) 
    return(objLP@lsDispModelH0)

#' @rdname get_accessors
#' @export
get_lsDispModelH1 <- function(objLP) 
    return(objLP@lsDispModelH1)

#' @rdname get_accessors
#' @export
get_lsDispModelConst <- function(objLP) 
    return(objLP@lsDispModelConst)

#' @rdname get_accessors
#' @export
get_lsDropModel <- function(objLP) 
    return(objLP@lsDropModel)

#' @rdname get_accessors
#' @export
get_matCountsProc <- function(objLP) 
    return(objLP@matCountsProc)

#' @rdname get_accessors
#' @export
get_matWeights <- function(objLP) 
    return(objLP@matWeights)

#' @rdname get_accessors
#' @export
get_scaDFSplinesDisp <- function(objLP) 
    return(objLP@scaDFSplinesDisp)

#' @rdname get_accessors
#' @export
get_scaDFSplinesMu <- function(objLP) 
    return(objLP@scaDFSplinesMu)

#' @rdname get_accessors
#' @export
get_strReport <- function(objLP) 
    return(objLP@strReport)

#' @rdname get_accessors
#' @export
get_vecAllGenes <- function(objLP) 
    return(objLP@vecAllGenes)

#' @rdname get_accessors
#' @export
get_vecConfoundersDisp <- function(objLP) 
    return(objLP@vecConfoundersDisp)

#' @rdname get_accessors
#' @export
get_vecConfoundersMu <- function(objLP) 
    return(objLP@vecConfoundersMu)

#' @rdname get_accessors
#' @export
get_scaOmega <- function(objLP) 
    return(objLP@scaOmega)

#' @rdname get_accessors
#' @export
get_boolFixedPopulations <- function(objLP) 
    return(objLP@boolFixedPopulations)

#' @rdname get_accessors
#' @export
get_vecH0Pop <- function(objLP) 
    return(objLP@vecH0Pop)

#' @rdname get_accessors
#' @export
get_vecNormConst <- function(objLP) 
    return(objLP@vecNormConst)

#' @rdname get_accessors
#' @export
get_strVersion <- function(objLP) 
    return(objLP@strVersion)

### 2. Enable accession of public elements via list-like
### properties of LineagePulseObject.

# a) Enable names()
#' List-like accessor methods for LineagePulseObject: names()
#' 
#' names() function for LineagePulseObject.
#' Allow usage of LineagePulse ouput object like a list with
#' respect to the most relevant output:
#' dfLineagePulseResults and vecDEGenes.
#' List of all available list-object like accessors:
#' \link{names,LineagePulseObject-method},
#' \link{[[,LineagePulseObject,character,missing-method},
#' \link{$,LineagePulseObject-method}.
#' 
#' @param x (LineagePulseObject) LineagePulse output object.
#' 
#' @return Names of elements in x available via list-like accessors.
#' 
#' @name names,LineagePulseObject-method
#' 
#' @examples    
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 10,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' names(objLP)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
setMethod('names', 'LineagePulseObject', function(x) {
  return( c("dfResults") )
})

# b) Enable object[[ element ]] operator
#' List-like accessor methods for LineagePulseObject: names()
#' 
#' names() function for LineagePulseObject.
#' Allow usage of LineagePulse ouput object like a list with
#' respect to the most relevant output:
#' dfLineagePulseResults and vecDEGenes.
#' List of all available list-object like accessors:
#' \link{names,LineagePulseObject-method},
#' \link{[[,LineagePulseObject,character,missing-method},
#' \link{$,LineagePulseObject-method}.
#' 
#' @param x (LineagePulseObject) LineagePulse output object.
#' @param i (str) Element from x list to be retrieved.
#' @param j () Ignored argument to generic.
#' @param ...  () Ignored argument to generic.
#' 
#' @return Target element from x.
#' 
#' @name [[,LineagePulseObject,character,missing-method
#' 
#' @examples    
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 10,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' head(objLP[["dfResults"]])
#' 
#' @author David Sebastian Fischer
#' 
#' @export
setMethod('[[', c('LineagePulseObject', 'character', 'missing'), function(x, i, j, ...){
  if(identical(i, "dfResults")){ return(x@dfResults)
  } else { return(NULL) }
})

# c) Enable object$element operator, which relies on [[ ]]
#' List-like accessor methods for LineagePulseObject: $
#' 
#' $ accessor function for LineagePulseObject, relies on [[ ]].
#' Allow usage of LineagePulse ouput object like a list with
#' respect to the most relevant output:
#' dfLineagePulseResults and vecDEGenes.
#' List of all available list-object like accessors:
#' \link{names,LineagePulseObject-method},
#' \link{[[,LineagePulseObject,character,missing-method},
#' \link{$,LineagePulseObject-method}.
#' 
#' @param x (LineagePulseObject) LineagePulse output object.
#' @param name (str) Element from x list to be retrieved.
#' 
#' @return Target element from x.
#' 
#' @name $,LineagePulseObject-method
#' 
#' @examples    
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 10,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' head(objLP$dfResults)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
setMethod('$', 'LineagePulseObject', function(x, name) x[[name]] )

### 3. Functions on LineagePulseObject that perform specific tasks

# a) Print report

#' Print LineagePulse report
#' 
#' Print LineagePulse report string to .txt file or stdout.
#'
#' @param object (LineagePulseObject) Output object of LineagePulse.
#' @param file (file) [DEFAULT ""] File to print report to. Default is stdout.
#' 
#' @return nothing to return
#'  
#' @examples    
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 10,
#'     scaNConst = 2,
#'     scaNLin = 2,
#'     scaNImp = 2,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 6),
#'     vecGeneWiseDropoutRates = rep(0.1, 6))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' writeReport(objLP, file="")
#'   
#' @author David Sebastian Fischer  
#' 
#' @export
writeReport <- function(object, file="") {
    write(object@strReport, file=file, ncolumns=1)
    return(NULL)
}
