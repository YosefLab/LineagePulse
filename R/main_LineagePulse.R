#' @import BiocParallel
#' @import circlize
#' @importFrom compiler cmpfun
#' @import ComplexHeatmap
#' @import gplots
#' @import ggplot2
#' @importFrom grDevices dev.off graphics.off pdf
#' @importFrom grid gpar
#' @import knitr
#' @import Matrix
#' @import methods
#' @import RColorBrewer
#' @import SingleCellExperiment
#' @import splines
#' @importFrom stats density dnbinom median optim p.adjust pchisq
#' quantile rnbinom rnorm runif sd lm qnbinom rbinom
#' @import SummarizedExperiment
#' @importFrom utils packageDescription
NULL

###############################################################################
### Libraries and source code
###############################################################################
# this section is for building if the code is not used as package but just as a
# collection of functions. Used for development.

#library(BiocParallel)
#library(compiler)
#library(ggplot2)
#library(Matrix)
#library(splines)
#library(SummarizedExperiment)

### CHANGE THIS PATH TO THE PATH IN WHICH YOU HAVE LineagePulse/R ###
#setwd("~/gitDevelopment/LineagePulse/R")

#source("srcLineagePulse_evalLogLikZINB.R")
#source("srcLineagePulse_calcPostDrop.R")
#source("srcLineagePulse_calcNormConst.R")
#source("srcLineagePulse_classLineagePulseObject.R")
#source("srcLineagePulse_decompressParameters.R")
#source("srcLineagePulse_evalDropoutModel.R")
#source("srcLineagePulse_evalImpulseModel.R")
#source("srcLineagePulse_evalLogLikNB.R")
#source("srcLineagePulse_evalLogLikZINB.R")
#source("srcLineagePulse_fitNB_fitMeanDispersion.R")
#source("srcLineagePulse_fitNB.R")
#source("srcLineagePulse_fitZINB_fitMeanDispersion.R")
#source("srcLineagePulse_fitZINB_fitDropout.R")
#source("srcLineagePulse_fitZINB.R")
#source("srcLineagePulse_fitZINB_WrapperLP.R")
#source("srcLineagePulse_getFits.R")
#source("srcLineagePulse_initialiseImpulseParameters.R")
#source("srcLineagePulse_initialiseImpulseParametersNB.R")
#source("srcLineagePulse_plotGene.R")
#source("srcLineagePulse_processSCData.R")
#source("srcLineagePulse_runDEAnalysis.R")
#source("srcLineagePulse_simulateDataSet.R")
#source("srcLineagePulse_sortGeneTrajectories.R")

###############################################################################
### Main function
###############################################################################

#' LineagePulse wrapper: Differential expression analysis on scRNA-seq
#' 
#' This function performs all steps of longitudinal or discrete differential
#' expression analysis in a continuous covariate (such as pseudotime) 
#' or according to a grouping (such as clusters or dconditions).
#' 
#' This function is the wrapper function for the LineagePulse algorithm
#' which performs differential expression analysis in pseudotime. 
#' Note that LineagePulse has many input parameters but
#' only few will be relevant for you and you will be able to leave the 
#' remaining ones as their defaults. Read up on specific input parameters
#' in the input parameter annotation of this function in the vignette.
#' 
#' @aliases LineagePulse wrapper, main function
#' 
#' @param counts (matrix genes x cells (sparseMatrix or standard), 
#' SummarizedExperiment or file)
#' Matrix: Count data of all cells, unobserved entries are NA.
#' SummarizedExperiment or SingleCellExperiment: 
#' Count data of all cells in assay(counts)
#' and annotation data can be supplied as colData(counts) or separately
#' via dfAnnotation.
#' file: .mtx file from which count matrix is to be read.
#' @param dfAnnotation (data frame cells x meta characteristics)
#' [Default NULL]
#' Annotation table which contains meta data on cells.
#' This data frame may be supplied as colData(counts) if
#' counts is a SummerizedExperiment or SingleCellExperiment object.
#' May contain the following columns
#' cell: Cell IDs.
#' continuous: Pseudotemporal coordinates of cells.
#' Confounder1: Batch labels of cells with respect 
#' to first confounder. Name is arbitrary: Could
#' for example be "patient" with batch labels
#' patientA, patientB, patientC.
#' Confounder2: As Confounder1 for another confounding
#' variable.
#' ... ConfounderX.
#' population: Fixed population assignments (for
#' strMuModel="MM"). Cells not assigned have to be NA.
#' groups: Discrete grouping of cells (e.g. clusters or experimental
#' conditions which are to be used as popuation structure if 
#' strMuModel or strDispModel are "groups").
#' rownames: Must be IDs from column cell.
#' Remaining entries in table are ignored.
#' @param vecConfoundersMu 
#' (vector of strings number of confounders on  mean)
#' [Default NULL] Confounders to correct for in mu batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param vecConfoundersDisp 
#' (vector of strings number of confounders on  dispersion)
#' [Default NULL] Confounders to correct for in dispersion batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param strMuModel (str) {"constant", "groups", "MM",
#' "splines","impulse"}
#' [Default "splines"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' population structure in the alternative model (H1).
#' @param strDispModelRed (str) {"constant", "groups", "splines"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' population structure in the null model (H0).
#' @param strDispModelFull (str) {"constant", "groups", "splines"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' population structure in the alternative model (H1).
#' @param strDropModel (str) {"logistic_ofMu", "logistic"}
#' [Default "logistic"] Definition of drop-out model.
#' "logistic_ofMu" - include the fitted mean in the linear model
#' of the drop-out rate and use offset and matPiConstPredictors.
#' "logistic" - only use offset and matPiConstPredictors.
#' @param strDropFitGroup (str) {"PerCell", "AllCells"}
#' [Defaul "PerCell"] Definition of groups on cells on which
#' separate drop-out model parameterisations are fit.
#' "PerCell" - one parametersiation (fit) per cell
#' "ForAllCells" - one parametersiation (fit) for all cells
#' @param scaDFSplinesMu (sca) [Default 6] 
#' If strMuModel=="splines", the degrees of freedom of the natural
#' cubic spline to be used as a mean parameter model.
#' @param scaDFSplinesDisp (sca) [Default 3] 
#' If strDispModelFull=="splines" or strDispModelRed=="splines", 
#' the degrees of freedom of the natural
#' cubic spline to be used as a dispersion parameter model.
#' @param matPiConstPredictors (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied
#' @param vecNormConstExternal (numeric vector number of cells) 
#' Model scaling factors, one per cell. These factors will linearly 
#' scale the mean model for evaluation of the loglikelihood. 
#' Must be named according to the column names of matCounts.
#' @param boolEstimateNoiseBasedOnH0 (bool) [Default TRUE]
#' Whether to co-estimate logistic drop-out model with the 
#' constant null model or with the alternative model. The
#' co-estimation with the noise model typically extends the
#' run-time of this model-estimation step strongly. While
#' the drop-out model is more accurate if estimated based on
#' a more realistic model expression model (the alternative
#' model), a trade-off for speed over accuracy can be taken
#' and the dropout model can be chosen to be estimated based
#' on the constant null expression model (set to TRUE).
#' @param scaMaxEstimationCycles (integer) [Default 20] Maximum number 
#' of estimation cycles performed in fitZINB(). One cycle
#' contain one estimation of of each parameter of the 
#' zero-inflated negative binomial model as coordinate ascent.
#' @param scaNProc (scalar) [Default 1] Number of processes for 
#' parallelisation.
#' @param boolVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation in high detail with local 
#' convergence flags and step-by-step loglikelihood computation.
#' 
#' @return dfDEAnalysis (data frame genes x reported variables) 
#' Summary of differential expression analysis:
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for differential expression with ZINB noise.
#' \item mean: Inferred mean parameter of constant model of first batch.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis with NB noise.
#' \item p_nb: P-value for differential expression with ZINB noise.
#' \item padj_nb: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis with NB noise.
#' \item loglik_full_zinb: Loglikelihood of full model with ZINB noise.
#' \item loglik_red_zinb: Loglikelihood of reduced model with ZINB noise.
#' \item loglik_full_nb: Loglikelihood of full model with NB noise.
#' \item loglik_red_nb: Loglikelihood of reduced model with NB noise.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item allZero (bool) Whether there were no observed non-zero observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 100,
#'     scaNConst = 10,
#'     scaNLin = 10,
#'     scaNImp = 10,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 30),
#'     vecGeneWiseDropoutRates = rep(0.1, 30))
#' matDropoutPredictors <- as.matrix(data.frame(
#'     log_means = log(rowMeans(lsSimulatedData$counts)+1) ))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "splines", scaDFSplinesMu = 6,
#'     strDropModel = "logistic", 
#'     matPiConstPredictors = matDropoutPredictors)
#' tail(objLP$dfResults)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runLineagePulse <- function(
    counts,
    dfAnnotation=NULL,
    vecConfoundersMu=NULL,
    vecConfoundersDisp=NULL,
    strMuModel="splines",
    strDispModelFull="constant",
    strDispModelRed="constant",
    strDropModel="logistic",
    strDropFitGroup="PerCell",
    scaDFSplinesMu=6,
    scaDFSplinesDisp=3,
    matPiConstPredictors=NULL,
    vecNormConstExternal=NULL,
    boolEstimateNoiseBasedOnH0=TRUE,
    scaMaxEstimationCycles=20,
    scaNProc=1,
    boolVerbose=TRUE,
    boolSuperVerbose=FALSE ){
    
    # 1. Data preprocessing
    # Extract count matrix if handed SummarizedExperiment
    # or SingleCellExperiment which extends SummarizedExperiment
    if (is(counts, "SummarizedExperiment")){ 
        counts <- assay(counts)
        if(is.null(dfAnnotation)){
            dfAnnotation <- colData(counts)
        }
    }
    
    vecAllGenes <- rownames(counts)
    lsProcessedSCData <- processSCData(
        counts=counts,
        dfAnnotation=dfAnnotation,
        vecConfoundersMu=vecConfoundersMu,
        vecConfoundersDisp=vecConfoundersDisp,
        matPiConstPredictors=matPiConstPredictors,
        vecNormConstExternal=vecNormConstExternal,
        strMuModel=strMuModel,
        strDispModelFull=strDispModelFull,
        strDispModelRed=strDispModelRed,
        scaDFSplinesMu=scaDFSplinesMu,
        scaDFSplinesDisp=scaDFSplinesDisp,
        scaMaxEstimationCycles=scaMaxEstimationCycles,
        boolVerbose=boolVerbose,
        boolSuperVerbose=boolSuperVerbose)
    objLP <- lsProcessedSCData$objLP
    vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
    matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
    
    # Inialise parallelisation
    # Set the parallelisation environment in BiocParallel:
    if(scaNProc > 1){
        register(MulticoreParam(workers=scaNProc)) 
    } else if(scaNProc == 1) {
        # For debugging in serial mode
        register(SerialParam())
    }
    
    # 2. Compute normalisation constants
    strMessage <- paste0("--- Compute normalisation constants:")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    vecNormConst(objLP) <- calcNormConst(objLP=objLP,
                           vecNormConstExternal=vecNormConstExternalProc)
    
    # 3. Fit ZINB and NB models for both H1 and H0.
    strMessage <- paste0("--- Fit ZINB model for both H1 and H0.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_fitmm <- system.time({
        objLP <- fitLPModels(
            objLP=objLP,
            matPiConstPredictors=matPiConstPredictorsProc,
            boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
            strMuModel=strMuModel,
            strDispModelFull=strDispModelFull,
            strDispModelRed=strDispModelRed,
            strDropModel=strDropModel,
            strDropFitGroup=strDropFitGroup,
            scaMaxEstimationCycles=scaMaxEstimationCycles,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose )
    })
    strMessage <- paste0("Time elapsed during ZINB fitting: ",
                         round(tm_fitmm["elapsed"]/60,2)," min")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    # 4. Differential expression analysis:
    strMessage <- paste0("--- Run differential expression analysis.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_deanalysis_mf <- system.time({
        objLP <- runDEAnalysis( objLP=objLP )
    })
    
    strMessage <- paste0("Finished runLineagePulse().")
    strReport(objLP) <- paste0(strReport(objLP), strMessage)
    if(boolVerbose) message(strMessage)
    
    return(objLP)
}