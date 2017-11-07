#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++     Fit full and alternative model    +++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit all models necessary for LineagePulse
#' 
#' Fit alternative H1 and null H0 on both ZINB and NB noise model
#' to a data set using cycles of coordinate ascent. The algorithm first
#' fits the either H1 or H0 together with the dropout model by
#' iterating over cell-wise (dropout models) and gene-wise (negative 
#' binomial models) parameters. Subsequently, the remaining model (H0 
#' or H1) is estimated by iterating over zero-inflated negative binomial
#' mean and dispersion parameter estimation condition on the previously
#' estimated logistic drop-out model. The NB noise model based models are
#' estimated in parallel across genes.
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
fitLPModels <- function(
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
    strMessage <- paste0("### a) Fit ZINB model A (",
                         strNameModelA,") with noise model.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_cycle <- system.time({
        lsFitsModelA <- fitModel(
            matCounts=matCountsProc(objLP),
            dfAnnotation=dfAnnotationProc(objLP),
            vecConfoundersMu=vecConfoundersMu(objLP),
            vecConfoundersDisp=vecConfoundersDisp(objLP),
            vecNormConst=vecNormConst(objLP),
            lsDropModel=NULL,
            strMuModel=strMuModelA,
            strDispModel=strDispModelA,
            strDropModel=strDropModel,
            strDropFitGroup=strDropFitGroup,
            scaDFSplinesDisp=scaDFSplinesDisp(objLP),
            scaDFSplinesMu=scaDFSplinesMu(objLP),
            matPiConstPredictors=matPiConstPredictors,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelA <- lsFitsModelA$lsMuModel
    lsDispModelA <- lsFitsModelA$lsDispModel
    lsDropModel <- lsFitsModelA$lsDropModel
    boolConvergenceModelA <- lsFitsModelA$boolConvergenceModel
    vecEMLogLikModelA <- lsFitsModelA$vecEMLogLikModel
    strReport(objLP) <- paste0(strReport(objLP),
                              lsFitsModelA$strReport)
    
    strMessage <- paste0(
        "Finished fitting zero-inflated negative binomial ",
        "model A with noise model in ", 
        round(tm_cycle["elapsed"]/60,2)," min.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    ####################################################
    # Fit model B
    strMessage <- paste0("### b) Fit ZINB model B (",
                         strNameModelB,").")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_cycleB <- system.time({
        lsFitsModelB <- fitModel(
            matCounts=matCountsProc(objLP),
            dfAnnotation=dfAnnotationProc(objLP),
            vecConfoundersMu=vecConfoundersMu(objLP),
            vecConfoundersDisp=vecConfoundersDisp(objLP),
            vecNormConst=vecNormConst(objLP),
            lsDropModel=lsDropModel,
            strMuModel=strMuModelB,
            strDispModel=strDispModelB,
            strDropModel=strDropModel,
            strDropFitGroup=strDropFitGroup,
            scaDFSplinesDisp=scaDFSplinesDisp(objLP),
            scaDFSplinesMu=scaDFSplinesMu(objLP),
            matPiConstPredictors=matPiConstPredictors,
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelB <- lsFitsModelB$lsMuModel
    lsDispModelB <- lsFitsModelB$lsDispModel
    boolConvergenceModelB <- lsFitsModelB$boolConvergenceModel
    vecEMLogLikModelB <- lsFitsModelB$vecEMLogLikModel
    strReport(objLP) <- paste0(strReport(objLP),
                              lsFitsModelB$strReport)
    
    strMessage <- paste0(
        "Finished fitting zero-inflated negative binomial ",
        "model B in ", round(tm_cycleB["elapsed"]/60,2)," min.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)

    ####################################################
    # Fit model C
    strMessage <- paste0("### c) Fit NB model A (",
                         strNameModelA,").")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_cycleA_NB <- system.time({
        lsFitsModelA_NB <- fitModel(
            matCounts=matCountsProc(objLP),
            dfAnnotation=dfAnnotationProc(objLP),
            vecConfoundersMu=vecConfoundersMu(objLP),
            vecConfoundersDisp=vecConfoundersDisp(objLP),
            vecNormConst=vecNormConst(objLP),
            lsDropModel=NULL,
            strMuModel=strMuModelA,
            strDispModel=strDispModelA,
            strDropModel="none",
            scaDFSplinesDisp=scaDFSplinesDisp(objLP),
            scaDFSplinesMu=scaDFSplinesMu(objLP),
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelA_NB <- lsFitsModelA_NB$lsMuModel
    lsDispModelA_NB <- lsFitsModelA_NB$lsDispModel
    boolConvergenceModelA_NB <- lsFitsModelA_NB$boolConvergenceModel
    strReport(objLP) <- paste0(strReport(objLP),
                               lsFitsModelA_NB$strReport)
    
    strMessage <- paste0(
        "Finished fitting NB ",
        "model B in ", round(tm_cycleA_NB["elapsed"]/60,2)," min.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    ####################################################
    # Fit model D
    strMessage <- paste0("### d) Fit NB model B (",
                         strNameModelB,").")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    tm_cycleB_NB <- system.time({
        lsFitsModelB_NB <- fitModel(
            matCounts=matCountsProc(objLP),
            dfAnnotation=dfAnnotationProc(objLP),
            vecConfoundersMu=vecConfoundersMu(objLP),
            vecConfoundersDisp=vecConfoundersDisp(objLP),
            vecNormConst=vecNormConst(objLP),
            lsDropModel=NULL,
            strMuModel=strMuModelB,
            strDispModel=strDispModelB,
            strDropModel="none",
            scaDFSplinesDisp=scaDFSplinesDisp(objLP),
            scaDFSplinesMu=scaDFSplinesMu(objLP),
            boolVerbose=boolVerbose,
            boolSuperVerbose=boolSuperVerbose)
    })
    lsMuModelB_NB <- lsFitsModelB_NB$lsMuModel
    lsDispModelB_NB <- lsFitsModelB_NB$lsDispModel
    boolConvergenceModelB_NB <- lsFitsModelB_NB$boolConvergenceModel
    strReport(objLP) <- paste0(strReport(objLP),
                              lsFitsModelB_NB$strReport)
    
    strMessage <- paste0(
        "Finished fitting NB ",
        "model B in ", round(tm_cycleB_NB["elapsed"]/60,2)," min.")
    strReport(objLP) <- paste0(strReport(objLP), strMessage, "\n")
    if(boolVerbose) message(strMessage)
        
    if(boolEstimateNoiseBasedOnH0){
        lsFitConvergence <- list(
            boolConvergenceH1=boolConvergenceModelB,
            boolConvergenceH0=boolConvergenceModelA,
            vecEMLogLikH1=vecEMLogLikModelB,
            vecEMLogLikH0=vecEMLogLikModelA,
            boolConvergenceH1_NB=boolConvergenceModelB_NB,
            boolConvergenceH0_NB=boolConvergenceModelA_NB)
        lsMuModelH1(objLP) <- lsMuModelB
        lsDispModelH1(objLP) <- lsDispModelB
        lsMuModelH0(objLP) <- lsMuModelA
        lsDispModelH0(objLP) <- lsDispModelA
        lsMuModelH1_NB(objLP) <- lsMuModelB_NB
        lsDispModelH1_NB(objLP) <- lsDispModelB_NB
        lsMuModelH0_NB(objLP) <- lsMuModelA_NB
        lsDispModelH0_NB(objLP) <- lsDispModelA_NB
    } else {
        lsFitConvergence <- list(
            boolConvergenceH1=boolConvergenceModelA,
            boolConvergenceH0=boolConvergenceModelB,
            vecEMLogLikH1=vecEMLogLikModelA,
            vecEMLogLikH0=vecEMLogLikModelB,
            boolConvergenceH1_NB=boolConvergenceModelA_NB,
            boolConvergenceH0_NB=boolConvergenceModelB_NB)
        lsMuModelH1(objLP) <- lsMuModelA
        lsDispModelH1(objLP) <- lsDispModelA
        lsMuModelH0(objLP) <- lsMuModelB
        lsDispModelH0(objLP) <- lsDispModelB
        lsMuModelH1_NB(objLP) <- lsMuModelA_NB
        lsDispModelH1_NB(objLP) <- lsDispModelA_NB
        lsMuModelH0_NB(objLP) <- lsMuModelB_NB
        lsDispModelH0_NB(objLP) <- lsDispModelB_NB
    }
    lsDropModel(objLP) <- lsDropModel
    lsFitConvergence(objLP) <- lsFitConvergence
    
    return(objLP)
}