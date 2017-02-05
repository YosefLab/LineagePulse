################################################################################
#################     LineagePulse output container class     ##################
################################################################################

### 1. Define output container class

# Define class unions for slots
setClassUnion('numericORNULL', members = c('numeric', 'NULL'))
setClassUnion('matrixORNULL', members = c('matrix', 'NULL'))
setClassUnion('characterORNULL', members = c('character', 'NULL'))
setClassUnion('listORNULL', members = c('list', 'NULL'))
setClassUnion('data.frameORNULL', members = c('data.frame', 'NULL'))

#' Container class for LineagePulse output
#' 
#' LineagePulse output and intermediate results such as model fits.
#' 
#' @slot dfDEAnalysis (data frame samples x reported characteristics) 
#'    Summary of fitting procedure and 
#'    differential expression results for each gene.
#'    \itemize{
#'      \item Gene: Gene ID.
#'      \item p: P-value for differential expression.
#'      \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#'      for differential expression analysis.
#'      \item loglik_full: Loglikelihood of full model.
#'      \item loglik_red: Loglikelihood of reduced model.
#'      \item df_full: Degrees of freedom of full model.
#'      \item df_red: Degrees of freedom of reduced model
#'      \item mean: Inferred mean parameter of constant model of first batch.
#'      From combined samples in case-ctrl. 
#'      \item allZero (bool) Whether there were no observed non-zero observations of this gene.
#'      If TRUE, fitting and DE analsysis were skipped and entry is NA.
#'    }
#'    Entries only present in case-only DE analysis:
#'    \itemize{
#'      \item converge_impulse: Convergence status of optim for 
#'      impulse model fit (full model).
#'      \item converge_const: Convergence status of optim for 
#'      constant model fit (reduced model).
#'    }
#'    Entries only present in mixture model DE analysis:
#'    \itemize{
#'      \item converge_mixture: Convergence status of optim for 
#'      mixture model fit (full model).
#'      \item converge_const: Convergence status of optim for 
#'      constant model fit (reduced model).
#'    }
#'    Entries only present if boolIdentifyTransients is TRUE:
#'    \itemize{
#'      \item converge_sigmoid: Convergence status of optim for 
#'      sigmoid model fit to samples of case condition.
#'      \item impulseTOsigmoid_p: P-value of loglikelihood ratio test
#'      impulse model fit versus sigmoidal model on samples of case condition.
#'      \item impulseTOsigmoid_padj: Benjamini-Hochberg 
#'      false-discovery rate corrected p-value of loglikelihood ratio test
#'      impulse model fit versus sigmoid model on samples of case condition.
#'      \item sigmoidTOconst_p: P-value of loglikelihood ratio test
#'      sigmoidal model fit versus constant model on samples of case condition.
#'      \item sigmoidTOconst_padj: Benjamini-Hochberg 
#'      false-discovery rate corrected p-value of loglikelihood ratio test
#'      sigmoidal model fit versus constant model on samples of case condition.
#'      \item isTransient (bool) Whether gene is transiently
#'      activated or deactivated and differentially expressed.
#'      \item isMonotonous (bool) Whether gene is not transiently
#'      activated or deactivated and differentially expressed. This scenario
#'      corresponds to a montonous expression level increase or decrease.
#'    }
#' @slot lsMuModelH1
#' @slot lsDispModelH1
#' @slot lsMuModelH0
#' @slot lsDispModelH0
#' @slot lsDropModel
#' @slot dfAnnotationProc
#' @slot strReport (str) LineagePulse stdout report.
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
    lsDropModel         = "listORNULL",
    lsMuModelH0         = "listORNULL",
    lsMuModelH1         = "listORNULL",
    lsFitZINBReporters  = "listORNULL",
    matCountsProc       = "matrixORNULL",
    matWeights          = "matrixORNULL",
    scaWindowRadius     = "numericORNULL",
    strReport           = "characterORNULL",
    vecAllGenes         = "characterORNULL",
    vecConfounders      = "characterORNULL",
    boolFixedPopulations= "logical",
    vecNormConst        = "numericORNULL",
    strVersion          = "character")
)

### 2. Enable accession of private elements via functions
### which carry the same name as the element.

#' LineagePulseObject accessor method generics
#' 
#' Generics for methods which operate on LineagePulseObject.
#'  
#' @param object (object) Object from which to retrieve data.
#' 
#' @aliases  get_lsMuModelH1
#'    get_lsDispModelH1
#'    get_lsMuModelH0
#'    get_lsDispModelH0
#'    get_DropModel
#'    get_lsFitZINBReporters
#'    get_dfAnnotationProc 
#'    get_vecNormConst
#'    get_scaNProc
#'    get_scaQThres
#'    get_strReport
#' 
#' @name LineagePulseObject_Generics_Accessors
NULL

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

# a) Enable printing of report to .txt file

#' Print LineagePulse report string to .txt file
#'
#' @param object (LineagePulseObject) Output object of LineagePulse.
#' @param fileReport (file) File to print report to.
#'  
#' @name writeReportToFile
#' @export
setGeneric('writeReportToFile', function(object, fileReport) standardGeneric('writeReportToFile'))

#' Print LineagePulse report string to .txt file
#' 
#' @param object (LineagePulseObject) Output object of LineagePulse.
#' @param fileReport (file) File to print report to.
#' 
#' @return NULL
#' 
#' @author David Sebastian Fischer
#' 
#' @name writeReportToFile,LineagePulseObject,character-method
#' @export
setMethod('writeReportToFile', signature(object='LineagePulseObject', fileReport='character'), 
          function(object, fileReport) write(object@strReport, file=fileReport, ncolumns=1) 
)

# b) Enable printing of report to stdout

#' Print LineagePulse report string to stdout
#'
#' @param object (LineagePulseObject) Output object of LineagePulse.
#'  
#' @name writeReportToStdout
#' @export
setGeneric('writeReportToStdout', function(object) standardGeneric('writeReportToStdout'))

#' Print LineagePulse report string to stdout
#' 
#' @param object (LineagePulseObject) Output object of LineagePulse.
#' 
#' @author David Sebastian Fischer
#' 
#' @name writeReportToStdout,LineagePulseObject-method
#' @export
setMethod('writeReportToStdout', signature(object='LineagePulseObject'), 
          function(object) write(object@strReport, file="", ncolumns=1) 
)
