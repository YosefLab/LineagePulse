#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++     Fit full and alternative model    ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit zero-inflated negative binomial models necessary for LineagePulse
#' hypothesis test to data
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
#' @seealso Called by \code{runLineagePulse}. 
#' Calls model estimation wrappers: \code{fitContinuousModels}.
#' 
#' @param objLP (LineagePulseObject)
#' LineagePulseObject to which null and alternative model are to be fitted.
#' @param matPiConstPredictors (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param strMuModel (str) {"constant", "groups", "MM",
#' "splines","impulse"}
#' [Default "impulse"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
#' @param strDispModelRed (str) {"constant", "groups", "splines"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' pseudotime in the null model (H0).
#' @param strDispModelFull (str) {"constant", "groups", "splines"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' pseudotime in the alternative model (H1).
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
#' @param boolEstimateNoiseBasedOnH0 (bool) [Default FALSE]
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
#' @param boolVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose (bool) Whether to follow convergence of the 
#' iterative parameter estimation in high detail with local 
#' convergence flags and step-by-step loglikelihood computation.
#' 
#' @return objLP (LineagePulseObject)
#' LineagePulseObject with models with and fitting reporters added.
#' 
#' @author David Sebastian Fischer
fitContinuousModels <- function(
    objLP,
    matPiConstPredictors,
    strMuModel="constant",
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