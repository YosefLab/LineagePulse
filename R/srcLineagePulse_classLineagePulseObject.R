###############################################################################
#################     LineagePulse output container class     #################
###############################################################################

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
#' Object containing description of 
#' gene-wise dispersion parameter models of H1.
#' @slot lsDispModelH0 (list)
#' Object containing description of 
#' gene-wise dispersion parameter models of H0.
#' @slot lsDispModelConst (list)
#' Object containing description of gene-wise dispersion parameter models of
#' constant model (necessary if H0 is not constant (mixture models)).
#' @slot lsDispModelH1_NB (list)
#' Object containing description of 
#' gene-wise dispersion parameter models of H1 with NB noise model.
#' @slot lsDispModelH0_NB (list)
#' Object containing description of 
#' gene-wise dispersion parameter models of H0 with NB noise model.
#' @slot lsDropModel (list)
#' @slot lsMuModelH0 (list)
#' Object containing description of gene-wise mean parameter models of H0.
#' @slot lsMuModelH1 (list)
#' Object containing description of gene-wise mean parameter models of H1.
#' @slot lsMuModelConst (list)
#' Object containing description of gene-wise means parameter models of
#' constant model (necessary if H0 is not constant (mixture models)).
#' @slot lsMuModelH0_NB (list)
#' Object containing description of gene-wise mean parameter models of H0
#' with NB noise model.
#' @slot lsMuModelH1_NB (list)
#' Object containing description of gene-wise mean parameter models of H1.
#' with NB noise model.
#' @slot lsFitConvergence (list)
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
#' @slot vecConfoundersDisp 
#' (vector of strings number of confounders on  dispersion)
#' Confounders to correct for in dispersion batch
#' correction model, are subset of column names of
#' dfAnnotation which describe condounding variables.
#' @slot vecConfoundersMu 
#' (vector of strings number of confounders on  mean)
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
        lsDispModelH0_NB    = "listORNULL",
        lsDispModelH1_NB    = "listORNULL",
        lsDropModel         = "listORNULL",
        lsMuModelH0         = "listORNULL",
        lsMuModelH1         = "listORNULL",
        lsMuModelConst      = "listORNULL",
        lsMuModelH0_NB      = "listORNULL",
        lsMuModelH1_NB      = "listORNULL",
        lsFitConvergence    = "listORNULL",
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
#' @param objLP (LineagePulse-Object)  
#' A LineagePulse output object to extract from.
#' 
#' @return The internal data object specified by the function.
#' 
#' @name accessors
#' @rdname accessors
#' @aliases
#' dfAnnotationProc
#' dfResults
#' lsMuModelH0
#' lsMuModelH1
#' lsMuModelConst
#' lsMuModelH0_NB
#' lsMuModelH1_NB
#' lsDispModelH0
#' lsDispModelH1
#' lsDispModelConst
#' lsDropModel
#' lsFitConvergence
#' matCountsProc
#' matWeights
#' scaDFSplinesDisp
#' scaDFSplinesMu
#' strReport
#' vecAllGenes
#' vecConfoundersDisp
#' vecConfoundersMu
#' scaOmega
#' boolFixedPopulation
#' vecH0Pop
#' vecNormConst
#' strVersion
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
#' dfAnnotationProc <- dfAnnotationProc(objLP)
#' dfResults <- dfResults(objLP)
#' lsMuModelH0 <- lsMuModelH0(objLP)
#' lsMuModelH1 <- lsMuModelH1(objLP)
#' lsMuModelConst <- lsMuModelConst(objLP)
#' lsMuModelH0_NB <- lsMuModelH0_NB(objLP)
#' lsMuModelH1_NB <- lsMuModelH1_NB(objLP)
#' lsDispModelH0 <- lsDispModelH0(objLP)
#' lsDispModelH1 <- lsDispModelH1(objLP)
#' lsDispModelConst <- lsDispModelConst(objLP)
#' lsDropModel <- lsDropModel(objLP)
#' lsFitConvergence <- lsFitConvergence(objLP)
#' matCountDataProc <- matCountsProc(objLP)
#' matWeights <- matWeights(objLP)
#' scaDFSplinesDisp <- scaDFSplinesDisp(objLP) 
#' scaDFSplinesMu <- scaDFSplinesMu(objLP) 
#' strReport <- strReport(objLP) 
#' vecAllGenes <- vecAllGenes(objLP) 
#' vecConfoundersDisp <- vecConfoundersDisp(objLP) 
#' vecConfoundersMu <- vecConfoundersMu(objLP) 
#' scaOmega <- scaOmega(objLP) 
#' boolFixedPopulations <- boolFixedPopulations(objLP)
#' vecH0Pop <- vecH0Pop(objLP) 
#' vecNormConst <- vecNormConst(objLP) 
#' strVersion <- strVersion(objLP)
#' 
#' @author David Sebastian Fischer
NULL

#' @rdname accessors
#' @export
dfAnnotationProc <- function(objLP) 
    return(objLP@dfAnnotationProc)

#' @rdname accessors
#' @export
dfResults <- function(objLP) 
    return(objLP@dfResults)

#' @rdname accessors
#' @export
lsMuModelH0 <- function(objLP) 
    return(objLP@lsMuModelH0)

#' @rdname accessors
#' @export
lsMuModelH1 <- function(objLP) 
    return(objLP@lsMuModelH1)

#' @rdname accessors
#' @export
lsMuModelConst <- function(objLP) 
    return(objLP@lsMuModelConst)

#' @rdname accessors
#' @export
lsMuModelH0_NB <- function(objLP) 
    return(objLP@lsMuModelH0_NB)

#' @rdname accessors
#' @export
lsMuModelH1_NB <- function(objLP) 
    return(objLP@lsMuModelH1_NB)

#' @rdname accessors
#' @export
lsDispModelH0 <- function(objLP) 
    return(objLP@lsDispModelH0)

#' @rdname accessors
#' @export
lsDispModelH1 <- function(objLP) 
    return(objLP@lsDispModelH1)

#' @rdname accessors
#' @export
lsDispModelConst <- function(objLP) 
    return(objLP@lsDispModelConst)

#' @rdname accessors
#' @export
lsDispModelH0_NB <- function(objLP) 
    return(objLP@lsDispModelH0_NB)

#' @rdname accessors
#' @export
lsDispModelH1_NB <- function(objLP) 
    return(objLP@lsDispModelH1_NB)

#' @rdname accessors
#' @export
lsDropModel <- function(objLP) 
    return(objLP@lsDropModel)

#' @rdname accessors
#' @export
lsFitConvergence <- function(objLP) 
    return(objLP@lsFitConvergence)

#' @rdname accessors
#' @export
matCountsProc <- function(objLP) 
    return(objLP@matCountsProc)

#' @rdname accessors
#' @export
matWeights <- function(objLP) 
    return(objLP@matWeights)

#' @rdname accessors
#' @export
scaDFSplinesDisp <- function(objLP) 
    return(objLP@scaDFSplinesDisp)

#' @rdname accessors
#' @export
scaDFSplinesMu <- function(objLP) 
    return(objLP@scaDFSplinesMu)

#' @rdname accessors
#' @export
strReport <- function(objLP) 
    return(objLP@strReport)

#' @rdname accessors
#' @export
vecAllGenes <- function(objLP) 
    return(objLP@vecAllGenes)

#' @rdname accessors
#' @export
vecConfoundersDisp <- function(objLP) 
    return(objLP@vecConfoundersDisp)

#' @rdname accessors
#' @export
vecConfoundersMu <- function(objLP) 
    return(objLP@vecConfoundersMu)

#' @rdname accessors
#' @export
scaOmega <- function(objLP) 
    return(objLP@scaOmega)

#' @rdname accessors
#' @export
boolFixedPopulations <- function(objLP) 
    return(objLP@boolFixedPopulations)

#' @rdname accessors
#' @export
vecH0Pop <- function(objLP) 
    return(objLP@vecH0Pop)

#' @rdname accessors
#' @export
vecNormConst <- function(objLP) 
    return(objLP@vecNormConst)

#' @rdname accessors
#' @export
strVersion <- function(objLP) 
    return(objLP@strVersion)

#' @rdname accessors
#' @export
strReport <- function(objLP) 
    return(objLP@strReport)

#' LineagePulseObject setters
#' 
#' Setters for LineagePulse-Object. Set the entry defined by the function
#' name to in objLP to value.
#' 
#' @param objLP (LineagePulse-Object)  
#' A LineagePulse output object to write into.
#' @param value (str) The value to be set.
#' 
#' @return The LineagePulse output object with updated entry.
#' 
#' @name LPsetters
#' @rdname LPsetters
#' @aliases
#' `dfAnnotationProc<-`
#' `dfResults<-`
#' `lsMuModelH0<-`
#' `lsMuModelH1<-`
#' `lsMuModelConst<-`
#' `lsMuModelH0_NB<-`
#' `lsMuModelH1_NB<-`
#' `lsDispModelH0<-`
#' `lsDispModelH1<-`
#' `lsDispModelConst<-`
#' `lsDispModelH0_NB<-`
#' `lsDispModelH1_NB<-`
#' `lsDropModel<-`
#' `lsFitConvergence<-`
#' `matCountsProc<-`
#' `matWeights<-`
#' `scaDFSplinesDisp<-`
#' `scaDFSplinesMu<-`
#' `strReport<-`
#' `vecAllGenes<-`
#' `vecConfoundersDisp<-`
#' `vecConfoundersMu<-`
#' `scaOmega<-`
#' `boolFixedPopulation<-`
#' `vecH0Pop<-`
#' `vecNormConst<-`
#' `strVersion<-`
#' 
#' @author David Sebastian Fischer
NULL

#' @rdname LPsetters
#' @export
`dfAnnotationProc<-` <- function(objLP, value) {
    objLP@dfAnnotationProc <- value
    objLP
}

#' @rdname LPsetters
#' @export
`dfResults<-` <- function(objLP, value) {
    objLP@dfResults <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsMuModelH0<-` <- function(objLP, value) {
    objLP@lsMuModelH0 <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsMuModelH1<-` <- function(objLP, value) {
    objLP@lsMuModelH1 <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsMuModelConst<-` <- function(objLP, value) {
    objLP@lsMuModelConst <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsMuModelH0_NB<-` <- function(objLP, value) {
    objLP@lsMuModelH0_NB <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsMuModelH1_NB<-` <- function(objLP, value) {
    objLP@lsMuModelH1_NB <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDispModelH0<-` <- function(objLP, value) {
    objLP@lsDispModelH0 <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDispModelH1<-` <- function(objLP, value) {
    objLP@lsDispModelH1 <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDispModelConst<-` <- function(objLP, value) {
    objLP@lsDispModelConst <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDispModelH0_NB<-` <- function(objLP, value) {
    objLP@lsDispModelH0_NB <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDispModelH1_NB<-` <- function(objLP, value) {
    objLP@lsDispModelH1_NB <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsDropModel<-` <- function(objLP, value) {
    objLP@lsDropModel <- value
    objLP
}

#' @rdname LPsetters
#' @export
`lsFitConvergence<-` <- function(objLP, value) {
    objLP@lsFitConvergence <- value
    objLP
}

#' @rdname LPsetters
#' @export
`matCountsProc<-` <- function(objLP, value) {
    objLP@matCountsProc <- value
    objLP
}

#' @rdname LPsetters
#' @export
`matWeights<-` <- function(objLP, value) {
    objLP@matWeights <- value
    objLP
}

#' @rdname LPsetters
#' @export
`scaDFSplinesDisp<-` <- function(objLP, value) {
    objLP@scaDFSplinesDisp <- value
    objLP
}

#' @rdname LPsetters
#' @export
`scaDFSplinesMu<-` <- function(objLP, value) {
    objLP@scaDFSplinesMu <- value
    objLP
}

#' @rdname LPsetters
#' @export
`strReport<-` <- function(objLP, value) {
    objLP@strReport <- value
    objLP
}

#' @rdname LPsetters
#' @export
`vecAllGenes<-` <- function(objLP, value) {
    objLP@vecAllGenes <- value
    objLP
}

#' @rdname LPsetters
#' @export
`vecConfoundersDisp<-` <- function(objLP, value) {
    objLP@vecConfoundersDisp <- value
    objLP
}

#' @rdname LPsetters
#' @export
`vecConfoundersMu<-` <- function(objLP, value) {
    objLP@vecConfoundersMu <- value
    objLP
}

#' @rdname LPsetters
#' @export
`scaOmega<-` <- function(objLP, value) {
    objLP@scaOmega <- value
    objLP
}

#' @rdname LPsetters
#' @export
`boolFixedPopulation<-` <- function(objLP, value) {
    objLP@boolFixedPopulation <- value
    objLP
}

#' @rdname LPsetters
#' @export
`vecH0Pop<-` <- function(objLP, value) {
    objLP@vecH0Pop <- value
    objLP
}

#' @rdname LPsetters
#' @export
`vecNormConst<-` <- function(objLP, value) {
    objLP@vecNormConst <- value
    objLP
}

#' @rdname LPsetters
#' @export
`strVersion<-` <- function(objLP, value) {
    objLP@strVersion <- value
    objLP
}

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
setMethod('[[', c('LineagePulseObject', 'character', 'missing'), 
          function(x, i, j, ...){
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
#' @param file (file) [DEFAULT ""] 
#' File to print report to. Default is stdout.
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
