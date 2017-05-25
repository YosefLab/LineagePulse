#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit zero-inflated negative binomial model to data
#' 
#' This is the algorithmic core wrapper of LineagePulse that carries out
#' the entire parameter estimation for both H0 and H1.
#'
#' Fit alternative H1 and null H0 zero-inflated negative binomial model 
#' to a data set using cycles of coordinate ascent. The algorithm first
#' fits the either H1 or H0 together with the logistic dropout model by
#' iterating over cell-wise (dropout models) and gene-wise (negative 
#' binomial models) parameters. Subsequently, the remaining model (H0 
#' or H1) is estimated by iterating over zero-inflated negative binomial
#' mean and dispersion parameter estimation condition on the previously
#' estimated logistic drop-out model.
#' 
#' Estimation of H0 and H1 are therefore separate coordinate ascent 
#' procedures yielding different local optima on the overall zero-inflated
#' negative binomial loglikelihood function with different gene-wise constraints.
#' 
#' Convergence is tracked with the the loglikelihood of the entire 
#' data matrix. Every step is a maximum likelihood estimation of the 
#' target parameters conditioned on the remaining parameter estimates. 
#' Therefore, convergence to a local optimum is guaranteed if the algorithm
#' is run until convergence. Parallelisation of each estimation step 
#' is implemented where conditional independences of parameter estimations
#' allow so. 
#' 
#' Convergence can be followed with verbose=TRUE (at each 
#' iteration) or at each step (boolSuperVerbose=TRUE). Variables for the
#' logistic drop-out model are a constant and the estimated mean parameter
#' and other constant gene-specific variables (such as GC-conten) in 
#' matPiConstPredictors. Three modes are available for modelling the mean
#' parameter: As a gene-wise constant (the default null model), by cluster 
#' (this is fast as neighbourhoods don't have to be evaluated), 
#' sliding windows (using neighbourhood smoothing), and as an impulse model.
#' 
#' To save memory, not the entire parameter matrix (genes x cells) but
#' the parmater models are stored in the objects lsMuModel, lsDispModel
#' and lsDropModel. These objects are described in detail in the annotation
#' of the return values of the function. In short, these object contain
#' the gene/cell-wise parameters of the model used to constrain the parameter
#' in question and the predictors necessary to evaluate the parameter model
#' to receive the observation-wise paramter values. Example: Impulse model
#' for the mean parameter: lsMuModel contains the parameter estimates for an
#' impulse model for each gene and pseudotime coordinates. Therefore, the
#' mean parameter for each observation can be computed as the value of the
#' impulse model evaluated at the pseudotime points for each gene.
#' 
#' @seealso Called by \code{runLineagePulse}. Calls parameter estimation
#' wrappers:
#' \code{fitPiZINB}, \code{fitZINBMu}, \code{fitZINBDisp} and
#' \code{fitZINBMuDisp}.
#' Calls \code{evalLogLikMatrix} to follow convergence.
#' 
#' @param objLP: (LineagePulseObject)
#' LineagePulseObject for which null and alternative model are to be fitted.
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param strMuModel: (str) {"constant"}
#' [Default "impulse"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
#' @param boolEstimateNoiseBasedOnH0: (bool) [Default: FALSE]
#' Whether to co-estimate logistic drop-out model with the 
#' constant null model or with the alternative model. The
#' co-estimation with the noise model typically extends the
#' run-time of this model-estimation step strongly. While
#' the drop-out model is more accurate if estimated based on
#' a more realistic model expression model (the alternative
#' model), a trade-off for speed over accuracy can be taken
#' and the dropout model can be chosen to be estimated based
#' on the constant null expression model (set to TRUE).
#' @param scaMaxEstimationCycles: (integer) [Default 20] Maximium number 
#' of estimation cycles performed in fitZINB(). One cycle
#' contain one estimation of of each parameter of the 
#' zero-inflated negative binomial model as coordinate ascent.
#' @param verbose: (bool) Whether to follow convergence of the 
#' iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose: (bool) Whether to follow convergence of the 
#' iterative parameter estimation in high detail with local 
#' convergence flags and step-by-step loglikelihood computation.
#' 
#' @return objLP: (LineagePulseObject)
#' LineagePulseObject with models with and fitting reporters added.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
fitContinuousModels <- function(
    objLP,
    matPiConstPredictors,
    strMuModel="windows",
    strDispModelFull="constant",
    strDispModelRed="constant",
    strDropModel="logistic_ofMu",
    strDropFitGroup="PerCell",
    boolEstimateNoiseBasedOnH0=TRUE,
    scaMaxEstimationCycles=20,
    boolVerbose=FALSE,
    boolSuperVerbose=FALSE ){
    
    ####################################################
    # Set sequence of model estimation: Which mean model
    # (H0 or H1) is co-estimated with noise model and which
    # is estimated conditioned on the noise model of this estimation.
    if(boolEstimateNoiseBasedOnH0){
        strMuModelA <- "constant"
        strMuModelB <- strMuModel
        strDispModelA <- strDispModelRed
        strDispModelB <- strDispModelFull
        strNameModelA <- paste0("H0: mu=",strMuModelA, " disp=", strDispModelA)
        strNameModelB <- paste0("H1: mu=",strMuModelB, " disp=", strDispModelB)
    } else {
        strMuModelA <- strMuModel
        strMuModelB <- "constant"
        strDispModelA <- strDispModelFull
        strDispModelB <- strDispModelRed
        strNameModelA <- paste0("H1: mu=",strMuModelA, " disp=", strDispModelA)
        strNameModelB <- paste0("H0: mu=",strMuModelB, " disp=", strDispModelB)
    }
    
    ####################################################
    # Fit model A
    strMessage <- paste0("### a) Fit negative binomial model A (",
                         strNameModelA,") with noise model.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    tm_cycle <- system.time({
        lsFitsModelA <- fitZINB(
            matCounts=objLP@matCountsProc,
            dfAnnotation=objLP@dfAnnotationProc,
            vecConfoundersMu=objLP@vecConfoundersMu,
            vecConfoundersDisp=objLP@vecConfoundersDisp,
            vecNormConst=objLP@vecNormConst,
            lsDropModel=NULL,
            strMuModel=strMuModelA,
            strDispModel=strDispModelA,
            strDropModel=strDropModel,
            strDropFitGroup=strDropFitGroup,
            scaDFSplinesDisp=objLP@scaDFSplinesDisp,
            scaDFSplinesMu=objLP@scaDFSplinesMu,
            matPiConstPredictors=matPiConstPredictors,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelA <- lsFitsModelA$lsMuModel
    lsDispModelA <- lsFitsModelA$lsDispModel
    lsDropModel <- lsFitsModelA$lsDropModel
    boolConvergenceModelA <- lsFitsModelA$boolConvergenceModel
    vecEMLogLikModelA <- lsFitsModelA$vecEMLogLikModel
    objLP@strReport <- paste0(objLP@strReport,
                              lsFitsModelA$strReport)
    rm(lsFitsModelA)
    
    strMessage <- paste0("Finished fitting zero-inflated negative binomial ",
                         "model A with noise model in ", round(tm_cycle["elapsed"]/60,2)," min.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    ####################################################
    # Fit model B
    strMessage <- paste0("### b) Fit negative binomial model B (",
                         strNameModelB,").")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    tm_cycleB <- system.time({
        lsFitsModelB <- fitZINB(
            matCounts=objLP@matCountsProc,
            dfAnnotation=objLP@dfAnnotationProc,
            vecConfoundersMu=objLP@vecConfoundersMu,
            vecConfoundersDisp=objLP@vecConfoundersDisp,
            vecNormConst=objLP@vecNormConst,
            lsDropModel=lsDropModel,
            strMuModel=strMuModelB,
            strDispModel=strDispModelB,
            strDropFitGroup=strDropFitGroup,
            scaDFSplinesDisp=objLP@scaDFSplinesDisp,
            scaDFSplinesMu=objLP@scaDFSplinesMu,
            matPiConstPredictors=matPiConstPredictors,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelB <- lsFitsModelB$lsMuModel
    lsDispModelB <- lsFitsModelB$lsDispModel
    boolConvergenceModelB <- lsFitsModelB$boolConvergenceModel
    vecEMLogLikModelB <- lsFitsModelB$vecEMLogLikModel
    objLP@strReport <- paste0(objLP@strReport,
                              lsFitsModelB$strReport)
    rm(lsFitsModelB)
    
    strMessage <- paste0("Finished fitting zero-inflated negative binomial ",
                         "model B in ", round(tm_cycleB["elapsed"]/60,2)," min.")
    objLP@strReport <- paste0(objLP@strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    if(boolEstimateNoiseBasedOnH0){
        lsFitZINBReporters <- list(
            boolConvergenceH1=boolConvergenceModelB,
            boolConvergenceH0=boolConvergenceModelA,
            vecEMLogLikH1=vecEMLogLikModelB,
            vecEMLogLikH0=vecEMLogLikModelA)
        objLP@lsMuModelH1 <- lsMuModelB
        objLP@lsDispModelH1 <- lsDispModelB
        objLP@lsMuModelH0 <- lsMuModelA
        objLP@lsDispModelH0 <- lsDispModelA
    } else {
        lsFitZINBReporters <- list(
            boolConvergenceH1=boolConvergenceModelA,
            boolConvergenceH0=boolConvergenceModelB,
            vecEMLogLikH1=vecEMLogLikModelA,
            vecEMLogLikH0=vecEMLogLikModelB)
        objLP@lsMuModelH1 <- lsMuModelA
        objLP@lsDispModelH1 <- lsDispModelA
        objLP@lsMuModelH0 <- lsMuModelB
        objLP@lsDispModelH0 <- lsDispModelB
    }
    objLP@lsDropModel <- lsDropModel
    objLP@lsFitZINBReporters <- lsFitZINBReporters
    
    return(objLP)
}