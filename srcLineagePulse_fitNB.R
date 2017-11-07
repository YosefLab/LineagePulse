#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++     Fit ZINB model to a data set    +++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit negative binomial model to data
#' 
#' This function fits a NB model with variable input. The NB model can be used
#' to perform model selection to test for drop-out.
#' 
#' To save memory, not the entire parameter matrix (genes x cells) but
#' the parmater models are stored in the objects lsMuModel and lsDispModel. 
#' In short, these object contain
#' the gene/cell-wise parameters of the model used to constrain the parameter
#' in question and the predictors necessary to evaluate the parameter model
#' to receive the observation-wise paramter values.
#' 
#' @seealso Called by \code{fitContinuousModels}. 
#' Calls parameter estimation wrappers:
#' \code{fitPiZINB}, \code{fitZINBMuDisp}.
#' Calls \code{evalLogLikMatrix} to follow convergence.
#' 
#' @param matCounts (matrix genes x cells)
#' Count data of all cells, unobserved entries are NA.
#' @param dfAnnotation (data frame cells x meta characteristics)
#' Annotation table which contains meta data on cells.
#' @param vecConfoundersMu (vector of strings number of confounders on  mean)
#' [Default NULL] Confounders to correct for in mu batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param vecConfoundersDisp 
#' (vector of strings number of confounders on dispersion)
#' [Default NULL] Confounders to correct for in dispersion batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param vecNormConst (numeric vector number of cells) 
#' Model scaling factors, one per cell. These factors linearly 
#' scale the mean model for evaluation of the loglikelihood.
#' @param scaDFSplinesMu (sca) [Default NULL] 
#' If strMuModel=="splines", the degrees of freedom of the natural
#' cubic spline to be used as a mean parameter model.
#' @param scaDFSplinesDisp (sca) [Default NULL] 
#' If strDispModelFull=="splines" or strDispModelRed=="splines", 
#' the degrees of freedom of the natural
#' cubic spline to be used as a dispersion parameter model.
#' @param matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' @param matMuModelInit (numeric matrix genes x mu model parameters)
#' [Default NULL]
#' Contains initialisation of mean model parameters 
#' according to the used model.
#' @param lsmatBatchModelInitMu (list) [Default NULL]
#' Initialisation of batch correction models for mean parameter.
#' @param matDispModelInit (numeric matrix genes x disp model parameters)
#' [Default NULL]
#' Contains initialisation of dispersion model parameters 
#' according to the used model.
#' @param lsmatBatchModelInitDisp (list) [Default NULL]
#' Initialisation of batch correction models for dispersion parameter.
#' @param strMuModel (str) {"constant", "groups", "MM",
#' "splines","impulse"}
#' [Default "impulse"] Model according to which the mean
#' parameter is fit to each gene as a function of 
#' population structure in the alternative model (H1).
#' @param strDispModel (str) {"constant", "groups", "splines"}
#' [Default "constant"] Model according to which dispersion
#' parameter is fit to each gene as a function of 
#' population structure in the given model.
#' @param boolVerbose (bool) [Default TRUE]
#' Whether to follow convergence of the 
#' iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose (bool) [Default TRUE]
#' Whether to follow convergence of the 
#' iterative parameter estimation in high detail with local 
#' convergence flags and step-by-step loglikelihood computation.
#' 
#' @return list
#' \itemize{
#' \item lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' \item lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' \item matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' \item boolConvergenceModel: (bool) 
#' Convergence status of model estimation.
#' \item strReport: (str) Log of model estimation to be added to 
#' overall log.
#' }
fitNB <- function(
    matCounts,
    dfAnnotation,
    vecConfoundersMu=NULL,
    vecConfoundersDisp=NULL,
    vecNormConst,
    scaDFSplinesMu=NULL,
    scaDFSplinesDisp=NULL,
    matWeights=NULL,
    matMuModelInit=NULL,
    lsmatBatchModelInitMu=NULL,
    matDispModelInit=NULL,
    lsmatBatchModelInitDisp=NULL,
    strMuModel,
    strDispModel,
    boolVerbose=TRUE,
    boolSuperVerbose=TRUE){
    
    ####################################################
    # Internal Numerical Estimation Parameters:
    # Numerical optmisation of impulse model hyperparameters
    MAXIT_BFGS_MuDisp <- 1000 # optim default is 1000
    RELTOL_BFGS_MuDisp <- 10^(-8) 
    # optim default is sqrt(.Machine$double.eps)=1e-8
    # Lowering RELTOL_BFGS_MuDisp gives drastic run time improvements.
    # Set to 10^(-4) to maintain sensible fits 
    # without running far into saturation
    # in the objective (loglikelihood).
    ####################################################
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    strReport <- ""
    
    ### Initialise
    # a) Mu model
    lsMuModel <- initMuModel(
        matCounts=matCounts,
        dfAnnotation=dfAnnotation,
        vecConfoundersMu=vecConfoundersMu,
        vecNormConst=vecNormConst,
        scaDFSplinesMu=scaDFSplinesMu,
        matWeights=matWeights,
        matMuModelInit=matMuModelInit,
        lsmatBatchModelInitMu=lsmatBatchModelInitMu,
        strMuModel=strMuModel,
        MAXIT_BFGS_MuDisp=MAXIT_BFGS_MuDisp,
        RELTOL_BFGS_MuDisp=RELTOL_BFGS_MuDisp)
    
    # b) Dispersion model
    lsDispModel <- initDispModel(
        matCounts=matCounts,
        dfAnnotation=dfAnnotation,
        vecConfoundersDisp=vecConfoundersDisp,
        scaDFSplinesDisp=scaDFSplinesDisp,
        matWeights=matWeights,
        matDispModelInit=matDispModelInit,
        lsmatBatchModelInitDisp=lsmatBatchModelInitDisp,
        strDispModel=strDispModel,
        strMuModel=strMuModel,
        MAXIT_BFGS_MuDisp=MAXIT_BFGS_MuDisp,
        RELTOL_BFGS_MuDisp=RELTOL_BFGS_MuDisp)
    
    # Evaluate initialisation loglikelihood
    scaLogLikNew <- sum(evalLogLikMatrixNB(
        matCounts=matCounts,
        lsMuModel=lsMuModel,
        lsDispModel=lsDispModel, 
        matWeights=matWeights ))
    strMessage <- paste0("#  .   Initialisation: ",
                         "ll ", scaLogLikNew)
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    ### Gene-wise parameter estimation:
    # Estimate mean and dispersion parameters simultaneously.
    # a/b) Negative binomial mean AND dispersion parameter.
    tm_mudisp <- system.time({
        lsFitMuDisp <- fitNBMuDisp(matCountsProc=matCounts,
                                     lsMuModel=lsMuModel,
                                     lsDispModel=lsDispModel,
                                     matWeights=matWeights)
    })
    lsDispModel$matDispModel <- lsFitMuDisp$matDispModel
    colnames(lsDispModel$matDispModel) <- NULL 
    # Need this so that column names dont grow to par.par.par...
    lsDispModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelDisp
    if(!is.null(lsDispModel$lsmatBatchModel)) {
        for(i in seq_along(lsDispModel$lsmatBatchModel)) {
            colnames(lsDispModel$lsmatBatchModel[[i]]) <- NULL 
        }
        # Need this so that column names dont grow to par.par.par...
    }
    lsMuModel$matMuModel <- lsFitMuDisp$matMuModel
    colnames(lsMuModel$matMuModel) <- NULL 
    # Need this so that column names dont grow to par.par.par...
    lsMuModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelMu
    if(!is.null(lsMuModel$lsmatBatchModel)) {
        for(i in seq_along(lsMuModel$lsmatBatchModel)) {
            colnames(lsMuModel$lsmatBatchModel[[i]]) <- NULL
        } 
        # Need this so that column names dont grow to par.par.par...
    }
    vecboolMuDispEstConverged <- lsFitMuDisp$vecConvergence
    
    # Report convergence problems
    if(any(vecboolMuDispEstConverged != 0)){
        strMessage <- paste0(
            "Mean/dispersion estimation did not converge in ", 
            sum(vecboolMuDispEstConverged[vecboolMuDispEstConverged>0]), 
            " cases [codes: ",
            paste(unique(vecboolMuDispEstConverged[
                vecboolMuDispEstConverged!=0]), collapse=","), "].")
        strReport <- paste0(strReport, strMessage, "\n")
        if(boolSuperVerbose) message(strMessage)
    }
    
    scaLogLikNew <- sum(lsFitMuDisp$vecLL)
    
    strMessage <- paste0(
        "# Mean+Disp co-estimation: ",
        "ll ", scaLogLikNew, " in ",
        round(tm_mudisp["elapsed"]/60,2)," min.")
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    # Name model matrix rows
    rownames(lsMuModel$matMuModel) <- rownames(matCounts)
    if(!is.null(lsMuModel$lsmatBatchModel)) {
        for(i in seq_along(lsMuModel$lsmatBatchModel)) {
            rownames(lsMuModel$lsmatBatchModel[[i]]) <- rownames(matCounts)
        }
    }
    rownames(lsDispModel$matDispModel) <- rownames(matCounts)
    if(!is.null(lsDispModel$lsmatBatchModel)) {
        for(i in seq_along(lsDispModel$lsmatBatchModel)) {
            rownames(lsDispModel$lsmatBatchModel[[i]]) <- rownames(matCounts)
        }
    }
    # Name model matrix columns
    if(strMuModel=="groups") {
        colnames(lsMuModel$matMuModel) <- unique(dfAnnotation$groups)
    } else if(strMuModel=="impulse") {
        colnames(lsMuModel$matMuModel) <- 
            c("beta1", "beta2", "h0", "h1", "h2", "t1", "t2")
    }

    return(list(
        lsMuModel=lsMuModel,
        lsDispModel=lsDispModel,
        strReport=strReport ))
}