################################################################################
#################     LineagePulse output container class     ##################
################################################################################

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
#' @slot dfAnnotationProc
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
#' @param object (objectLineagePulse)  A LineagePulse output object.
#' 
#' @return The internal data object specified by the function.
#' 
#' @aliases get_lsMuModelH1,LineagePulseObject-method
#'    get_lsDispModelH1,LineagePulseObject-method
#'    get_lsMuModelH0,LineagePulseObject-method
#'    get_lsDispModelH0,LineagePulseObject-method
#'    get_DropModel,LineagePulseObject-method
#'    get_lsFitZINBReporters
#'    get_dfAnnotationProc,LineagePulseObject-method
#'    get_vecNormConst,LineagePulseObject-method
#'    get_scaNProc,LineagePulseObject-method
#'    get_scaQThres,LineagePulseObject-method
#'    get_strReport,LineagePulseObject-method
#' 
#' @name LineagePulseObject_Accessors
NULL

### I. Set generic function which defines string as a function:
### setGeneric('funName', function(object) standardGeneric('funName'), valueClass = 'funOutputClass')
### II. Define function on LineagePulseObject:
### setMethod('funName', 'LineagePulseObject', function(object) object@funName)

#' @return (list) lsMuModelH1
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_lsMuModelH1', function(object) standardGeneric('get_lsMuModelH1'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_lsMuModelH1', 'LineagePulseObject', function(object) object@lsMuModelH1)

#' @return (list) lsDispModelH1
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_lsDispModelH1', function(object) standardGeneric('get_lsDispModelH1'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_lsDispModelH1', 'LineagePulseObject', function(object) object@lsDispModelH1)

#' @return (list) lsMuModelH0
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_lsMuModelH0', function(object) standardGeneric('get_lsMuModelH0'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_lsMuModelH0', 'LineagePulseObject', function(object) object@lsMuModelH0)

#' @return (list) lsDispModelH0
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_lsDispModelH0', function(object) standardGeneric('get_lsDispModelH0'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_lsDispModelH0', 'LineagePulseObject', function(object) object@lsDispModelH0)

#' @return (list) lsDropModel
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_DropModel', function(object) standardGeneric('get_DropModel'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_DropModel', 'LineagePulseObject', function(object) object@lsDropModel)

#' @return (list) lsFitZINBReporters
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_lsFitZINBReporters', function(object) standardGeneric('get_lsFitZINBReporters'), valueClass = 'listORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_lsFitZINBReporters', 'LineagePulseObject', function(object) object@lsFitZINBReporters)

#' @return (data frame size genes x reported characteristics) dfAnnotationProc
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_dfAnnotationProc', function(object) standardGeneric('get_dfAnnotationProc'), valueClass = 'data.frame')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_dfAnnotationProc', 'LineagePulseObject', function(object) object@dfAnnotationProc)

#' @return (numeric vector length number of samples) vecNormConst
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_vecNormConst', function(object) standardGeneric('get_vecNormConst'), valueClass = 'numeric')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_vecNormConst', 'LineagePulseObject', function(object) object@vecNormConst)

#' @return (scalar) scaNProc
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_scaNProc', function(object) standardGeneric('get_scaNProc'), valueClass = 'numeric')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_scaNProc', 'LineagePulseObject', function(object) object@scaNProc)

#' @return (scalar) scaQThres
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_scaQThres', function(object) standardGeneric('get_scaQThres'), valueClass = 'numericOrNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_scaQThres', 'LineagePulseObject', function(object) object@scaQThres)

#' @return (str) strReport
#' @name LineagePulseObject_Generics_Accessors
#' @export
setGeneric('get_strReport', function(object) standardGeneric('get_strReport'), valueClass = 'characterORNULL')
#' @name LineagePulseObject_Accessors
#' @export
setMethod('get_strReport', 'LineagePulseObject', function(object) object@strReport)

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
#' @export
writeReport <- function(object, file="") write(object@strReport, file=file, ncolumns=1)
