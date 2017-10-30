#' @import BiocParallel
#' @importFrom compiler cmpfun
#' @import gplots
#' @import ggplot2
#' @importFrom grDevices dev.off graphics.off pdf
#' @import knitr
#' @import Matrix
#' @import methods
#' @import splines
#' @importFrom stats dnbinom median optim p.adjust pchisq 
#' rnbinom rnorm runif sd lm qnbinom rbinom
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
#source("srcLineagePulse_fitZINB_fitMeanDispersion.R")
#source("srcLineagePulse_fitZINB_fitDropout.R")
#source("srcLineagePulse_fitZINB.R")
#source("srcLineagePulse_fitZINB_WrapperLP.R")
#source("srcLineagePulse_getFits.R")
#source("srcLineagePulse_initialiseImpulseParameters.R")
#source("srcLineagePulse_plotGene.R")
#source("srcLineagePulse_processSCData.R")
#source("srcLineagePulse_runDEAnalysis.R")
#source("srcLineagePulse_simulateDataSet.R")
#source("srcLineagePulse_sortGeneTrajectories.R")

###############################################################################
### Main function
###############################################################################

#' LineagePulse wrapper: Differential expression analysis in pseudotime
#' 
#' This function performs all steps of longitudinal differential
#' expression analysis in pseudotime for you.
#' 
#' This function is the wrapper function for the LineagePulse algorithm
#' which performs differential expression analysis in pseudotime. 
#' Note that LineagePulse has many input parameters but
#' only few will be relevant for you and you will be able to leave the 
#' remaining ones as their defaults. Read up on specific input parameters
#' in the input parameter annotation of this function or follow this short
#' guide:
#' 
#' MINIMAL INPUT
#' 1. Supply data: Count data (matCounts) and pseudotime coordinates
#' of the cells (vecPseudotime). You may decide to also provide cell-wise
#' normalisation factors (such as factors accounting for sequencing depth)
#' (vecNormConstExternal).
#' 2. Supply cell meta data: Data frame (dfAnnotation) which contains
#' all cell ids (the colnames of matCounts) in a column "cell",
#' and cell-wise pseudotime coordinates (scalar) in a a column "pseudpotime".
#' Rownames must be the ids in column "cell".
#' 3. Chose the model constraining mean (strMuModel) and dispersion 
#' parameters (strDispModel) for each gene. 
#' 
#' ADDITIONAL FACULTATIVE SETTINGS
#' 5. Supply gene-specific drop-out predictors if wanted 
#' (matPiConstPredictors).
#' 6. Set optimisation parameters (boolEstimateNoiseBasedOnH0,
#' scaMaxEstimationCycles).
#' 7. Chose the number of processes you want to use (scaNProc), LineagePulse
#' is parallelised on all computation intensive steps. Note that
#' the current parallelisation scheme runs on Unix (MacOS) and Linux but
#' not on Windows. Set scaNProc to NA and start a BiocParallel environment
#' before the call to this function to use parallelisation on Windows.
#' 8. Set the level of detail with which you want to follow
#' progress through text printed on the console during a run
#' (boolVerbose, boolSuperVerbose).
#' 
#' Finally, after running LineagePulse, you may continue to work
#' on your data set by:
#' A) Regenerating observation-wise parameters,
#' such as the mean parameter matrix which represents the hidden
#' expression states, with the functions
#' \link{getFitsMean}, \link{getFitsDispersion}, \link{getFitsDropout}.
#' B) You can also compute the observation-wise probability of 
#' dropout with \link{getPostDrop}.
#' C) You can extract LineagePulse model normalised and/or
#' batch corrected data with \link{getNormData}, 
#' D) You can have a closer look at the global expression
#' trajecotries of the genes with \link{sortGeneTrajectories}.
#' E) You can look at gene-wise model fits with \link{plotGene}.
#' 
#' @aliases LineagePulse wrapper, main function
#' 
#' @param counts (matrix genes x cells (sparseMatrix or standard), 
#' SummarizedExperiment or file)
#' Matrix: Count data of all cells, unobserved entries are NA.
#' SummarizedExperiment: Count data of all cells in assay(counts).
#' file: .mtx file from which count matrix is to be read.
#' @param dfAnnotation (data frame cells x meta characteristics)
#' Annotation table which contains meta data on cells.
#' May contain the following columns
#' cell: Cell IDs.
#' pseudotime: Pseudotemporal coordinates of cells.
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
#' [Default "impulse"] Model according to which the mean
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
#' [Default "logistic_ofMu"] Definition of drop-out model.
#' "logistic_ofMu" - include the fitted mean in the linear model
#' of the drop-out rate and use offset and matPiConstPredictors.
#' "logistic" - only use offset and matPiConstPredictors.
#' @param strDropFitGroup (str) {"PerCell", "AllCells"}
#' [Defaul "PerCell"] Definition of groups on cells on which
#' separate drop-out model parameterisations are fit.
#' "PerCell" - one parametersiation (fit) per cell
#' "ForAllCells" - one parametersiation (fit) for all cells
#' @param scaDFSplinesMu (sca) [Default 3] 
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
#' gene: gene ID,
#' p: raw p-value, 
#' adj.p: BH corrected p-value, 
#' loglik_full: loglikelihood of alternative model H1,
#' loglik_red: loglikelihood of null model H0,
#' deviance: loglikelihood ratio test statistic (the deviance),
#' mean_H0: inferred gene-wise mean parameter (constant null model),
#' dispersion_H0: inferred gene-wise dispersion parameter 
#' (constant null model)
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
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' tail(objLP$dfResults)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runLineagePulse <- function(
    counts,
    dfAnnotation,
    vecConfoundersMu=NULL,
    vecConfoundersDisp=NULL,
    strMuModel="impulse",
    strDispModelFull="constant",
    strDispModelRed="constant",
    strDropModel="logistic_ofMu",
    strDropFitGroup="PerCell",
    scaDFSplinesMu=3,
    scaDFSplinesDisp=3,
    matPiConstPredictors=NULL,
    vecNormConstExternal=NULL,
    boolEstimateNoiseBasedOnH0=TRUE,
    scaMaxEstimationCycles=20,
    scaNProc=1,
    boolVerbose=TRUE,
    boolSuperVerbose=FALSE ){
    
    STRVERSION <- packageDescription("LineagePulse", fields = "Version")
    
    # 1. Data preprocessing
    # Extract count matrix if handed SummarizedExperiment
    if (class(counts) == "SummarizedExperiment"){ 
        counts <- assay(counts)
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
        boolSuperVerbose=boolSuperVerbose,
        STRVERSION=STRVERSION)
    objLP <- lsProcessedSCData$objLP
    vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
    matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
    
    # Clear memory
    rm(counts)
    rm(matPiConstPredictors)
    rm(dfAnnotation)
    rm(lsProcessedSCData)
    gc()
    
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
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    objLP <- calcNormConst(objLP=objLP,
                           vecNormConstExternal=vecNormConstExternalProc)
    
    # 3. Fit ZINB model for both H1 and H0.
    strMessage <- paste0("--- Fit ZINB model for both H1 and H0.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_fitmm <- system.time({
        objLP <- fitContinuousModels(
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
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    
    # 4. Differential expression analysis:
    strMessage <- paste0("--- Run differential expression analysis.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_deanalysis_mf <- system.time({
        objLP <- runDEAnalysis( objLP=objLP )
    })
    strMessage <- paste0(
        "Time elapsed during differential expression analysis: ",
        round(tm_deanalysis_mf["elapsed"]/60,2)," min")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    strMessage <- paste0("Finished runLineagePulse().")
    objLP@strReport <- paste0(objLP@strReport, strMessage)
    if(boolVerbose) message(strMessage)
    
    return(objLP)
}