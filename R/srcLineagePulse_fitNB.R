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
    RELTOL_BFGS_MuDisp <- 10^(-4) 
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
    # Mean parameters (mu): Gene-wise mean of non-zero observations.
    # Impulse model: Initialised to constant (mean).
    lsMuModel <- list(
        matMuModel=NULL,
        lsmatBatchModel=NULL,
        lsMuModelGlobal=list(
            scaDegFreedom=NULL, # general parameters
            strMuModel=strMuModel,
            vecNormConst=vecNormConst,
            scaNumCells=scaNumCells,
            vecPseudotime=dfAnnotation$pseudotime, # continuos models
            scaNSplines=NULL, # spline models
            matSplineBasis=NULL,
            vecidxGroups=NULL, # group models
            vecGroups=NULL,
            scaNGroups=NULL,
            scaNMix=NULL, # mixture models
            vecConfounders=vecConfoundersMu, # batch models
            scaNConfounders=NULL,
            vecNBatches=NULL, 
            # auxillary, can be computed as sapply(lsvecidxBatchAssign, max)
            lsvecidxBatchAssign=NULL, 
            # object called to determine assignment of cell to batch
            lsvecBatchUnique=NULL, 
            # Kept so that idx vector can be associated with batch names
            MAXIT_BFGS_MuDisp=MAXIT_BFGS_MuDisp,
            RELTOL_BFGS_MuDisp=RELTOL_BFGS_MuDisp) )
    # Initialise mean model parameters
    if(is.null(matMuModelInit)){
        vecMuModelInit <- Matrix::rowMeans(matCounts, na.rm = TRUE) 
        vecMuModelInit[vecMuModelInit < 10^(-10)] <- 10^(-10)
        if(strMuModel=="constant"){
            lsMuModel$matMuModel <- matrix(
                vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
        } else if(strMuModel=="impulse"){
            lsMuModel$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=7)
            lsMuModel$matMuModel[,c(3:5)] <- matrix(
                vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE)
        } else if(strMuModel=="splines"){
            # Note : counts ~ exp(W*C) where W is spline basis 
            # and C are coefficients
            # Therefore: log(counts) ~ W*C is a very similar problem.
            # We ignore ZINB noise and initialise the coefficients based 
            # on this problem with gaussian noise 
            # (linear model optimised with RSS and without intercept: 
            # notation ylm(~0+x)).
            vecPTSpline <- ns(x = dfAnnotation$pseudotime, 
                              df = scaDFSplinesMu, intercept = TRUE)
            lsMuModel$matMuModel <- do.call(
                rbind, lapply(seq_len(scaNumGenes) ,function(i){
                    vecCounts <- as.vector(log(matCounts[as.double(i),]+1))
                    lmFit <- lm(vecCounts ~ 0+vecPTSpline)
                    return(as.vector(lmFit$coef))
                }))
        } else if(strMuModel=="groups"){
            lsMuModel$matMuModel <- matrix(
                vecMuModelInit, nrow=scaNumGenes, 
                ncol=length(unique(dfAnnotation$groups)), byrow=FALSE)
        } else if(strMuModel=="MM"){
            lsMuModel$matMuModel <- matrix(
                vecMuModelInit, nrow=scaNumGenes, 
                ncol=dim(matWeights)[2], byrow=FALSE)
        } else {
            stop("ERROR fitZINB(): strMuModel=", strMuModel, 
                 " not recognised.")
        }
    } else {
        lsMuModel$matMuModel <- matMuModelInit
    }
    if(strMuModel=="MM"){
        lsMuModel$lsMuModelGlobal$scaNMix <- dim(matWeights)[2]
    }
    if(strMuModel=="groups"){
        vecGroupsMu <- as.vector(dfAnnotation$groups)
        lsMuModel$lsMuModelGlobal$scaNGroups <- length(unique(vecGroupsMu))
        lsMuModel$lsMuModelGlobal$vecGroups <- unique(vecGroupsMu)
        lsMuModel$lsMuModelGlobal$vecidxGroups <- match(vecGroupsMu, 
                                                        unique(vecGroupsMu))
    }
    if(strMuModel=="splines"){
        lsMuModel$lsMuModelGlobal$vecTimepoints <- 
            sort(unique(dfAnnotation$pseudotime ))
        lsMuModel$lsMuModelGlobal$vecindTimepointAssign <- match(
            dfAnnotation$pseudotime, lsMuModel$lsMuModelGlobal$vecTimepoints)
        lsMuModel$lsMuModelGlobal$scaNSplines <- scaDFSplinesMu
        lsMuModel$lsMuModelGlobal$matSplineBasis <- 
            ns(x = dfAnnotation$pseudotime, 
               df = scaDFSplinesMu, intercept = TRUE)[,,drop=FALSE]
    }
    if(strMuModel=="impulse"){
        lsMuModel$lsMuModelGlobal$vecTimepoints <- 
            sort(unique(dfAnnotation$pseudotime ))
        lsMuModel$lsMuModelGlobal$vecindTimepointAssign <- match(
            dfAnnotation$pseudotime, lsMuModel$lsMuModelGlobal$vecTimepoints)
    }
    if(is.null(lsmatBatchModelInitMu)){
        # Initialise batch model parameters
        if(!is.null(vecConfoundersMu)){
            lsvecBatchAssign <- 
                lapply(vecConfoundersMu, function(confounder) {
                    dfAnnotation[[confounder]] 
                })
            # Initialise batch correction factors to 1
            lsMuModel$lsmatBatchModel <- 
                lapply(lsvecBatchAssign, function(vecBatchAssign){
                    matrix(1, nrow=scaNumGenes,
                           ncol=length(unique(vecBatchAssign)) )
                })
        } else {
            lsMuModel$lsmatBatchModel <- list(matrix(NA, nrow=scaNumGenes, 
                                                     ncol=1 ))
        }
    } else {
        lsMuModel$lsmatBatchModel <- lsmatBatchModelInitMu
    }
    # Add global batch parameters
    if(!is.null(vecConfoundersMu)){
        lsMuModel$lsMuModelGlobal$scaNConfounders <- length(vecConfoundersMu)
        lsvecBatchesMu <- 
            lapply(vecConfoundersMu, function(confounder) { 
                dfAnnotation[[confounder]] 
            })
        lsMuModel$lsMuModelGlobal$lsvecBatchUnique <- 
            lapply(lsvecBatchesMu, function(vecBatchAssign) { 
                unique(vecBatchAssign) 
            })
        lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign <- 
            lapply(lsvecBatchesMu, function(vecBatchAssign) { 
                match(vecBatchAssign, unique(vecBatchAssign)) 
            })
        lsMuModel$lsMuModelGlobal$vecNBatches <- 
            sapply(lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign, max)
    } 
    # Mu model + batch correction model
    lsMuModel$lsMuModelGlobal$scaDegFreedom <- 
        dim(lsMuModel$matMuModel)[2] + 
        sum(sapply(lsMuModel$lsmatBatchModel, function(mat) dim(mat)[2]-1 ))
    
    # b) Dispersion model
    # Dispersions: Low dispersion factor yielding high variance which makes
    # cost function screening easy in the first iteration.
    lsDispModel <- list(matDispModel=NULL,
                        lsmatBatchModel=NULL,
                        lsDispModelGlobal=list(
                            scaDegFreedom=NULL, # general parameters
                            strDispModel=strDispModel,
                            scaNumCells=scaNumCells,
                            vecPseudotime=dfAnnotation$pseudotime,
                            scaNSplines=NULL, # spline models
                            matSplineBasis=NULL,
                            vecidxGroups=NULL, # group models
                            vecGroups=NULL,
                            scaNGroups=NULL,
                            scaNMix=NULL, # mixture model
                            vecConfounders=vecConfoundersDisp, # batch model
                            scaNConfounders=NULL,
                            vecNBatches=NULL,
                            lsvecidxBatchAssign=NULL,
                            lsvecBatchUnique=NULL ) )
    if(is.null(matDispModelInit)){
        if(strDispModel=="constant"){
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, 
                                               ncol=1, byrow=FALSE)
        } else if(strDispModel=="MM"){
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, 
                                               ncol=dim(matWeights)[2])
        } else if(strDispModel=="splines"){
            lsDispModel$matDispModel <- matrix(0, nrow=scaNumGenes, 
                                               ncol=scaDFSplinesDisp)
            lsDispModel$matDispModel[,1] <- 1
        } else if(strDispModel=="groups"){
            lsDispModel$matDispModel <- 
                matrix(1, nrow=scaNumGenes, 
                       ncol=length(unique(dfAnnotation$groups)))
        } else {
            stop("ERROR fitZINB(): strDispModel=", strDispModel, 
                 " not recognised.")
        } 
    } else {
        lsDispModel$matDispModel <- matDispModelInit
    }
    if(is.null(lsmatBatchModelInitDisp)){
        # Initialise batch model parameters
        if(!is.null(vecConfoundersDisp)){
            lsvecBatchAssign <- 
                lapply(vecConfoundersDisp, function(confounder) {
                    dfAnnotation[[confounder]] 
                })
            # Initialise batch correction factors to 1
            lsDispModel$lsmatBatchModel <- 
                lapply(lsvecBatchAssign, function(vecBatchAssign){
                    matrix(1, nrow=scaNumGenes,
                           ncol=length(unique(vecBatchAssign)) )
                })
        } else {
            lsDispModel$lsmatBatchModel <- 
                list(matrix(NA, nrow=scaNumGenes, ncol=1 ))
        }
    } else {
        lsDispModel$lsmatBatchModel <- lsmatBatchModelInitDisp
    }
    if(strDispModel=="groups"){
        vecGroupsDisp <- as.vector(dfAnnotation$groups)
        lsDispModel$lsDispModelGlobal$scaNGroups <- 
            length(unique(vecGroupsDisp))
        lsDispModel$lsDispModelGlobal$vecGroups <- 
            vecGroupsDisp
        lsDispModel$lsDispModelGlobal$vecidxGroups <- 
            match(vecGroupsDisp, unique(vecGroupsDisp))
    }
    if(strDispModel=="splines"){
        lsDispModel$lsDispModelGlobal$vecTimepoints <- 
            sort(unique(dfAnnotation$pseudotime ))
        lsDispModel$lsDispModelGlobal$vecindTimepointAssign <- match(
            dfAnnotation$pseudotime, 
            lsDispModel$lsDispModelGlobal$vecTimepoints)
        lsDispModel$lsDispModelGlobal$scaNSplines <- scaDFSplinesDisp
        lsDispModel$lsDispModelGlobal$matSplineBasis <- 
            ns(x = dfAnnotation$pseudotime, 
               df = scaDFSplinesDisp, intercept = TRUE)[,,drop=FALSE]
    }
    if(strMuModel=="MM"){ # Also have to fill this in if only strMuModel is MM
        lsDispModel$lsDispModelGlobal$scaNMix <- dim(matWeights)[2]
    }
    if(!is.null(vecConfoundersDisp)) {
        lsDispModel$lsDispModelGlobal$scaNConfounders <- 
            length(vecConfoundersDisp)
        lsvecBatchesDisp <- 
            lapply(vecConfoundersDisp, function(confounder) {
                as.vector(dfAnnotation[[confounder]]) 
            })
        lsDispModel$lsDispModelGlobal$lsvecBatchUnique <- 
            lapply(lsvecBatchesDisp, function(vecBatchAssign) {
                unique(vecBatchAssign) 
            })
        lsDispModel$lsDispModelGlobal$lsvecidxBatchAssign <- 
            lapply(lsvecBatchesDisp, function(vecBatchAssign) {
                match(vecBatchAssign, unique(vecBatchAssign)) 
            })
        lsDispModel$lsDispModelGlobal$vecNBatches <- 
            sapply(lsDispModel$lsDispModelGlobal$lsvecidxBatchAssign, max)
    }
    # disp model + batch correction model
    lsDispModel$lsDispModelGlobal$scaDegFreedom <- 
        dim(lsDispModel$matDispModel)[2] + 
        sum(lsDispModel$lsDispModelGlobal$vecNBatches) - 
        length(lsDispModel$lsDispModelGlobal$vecNBatches) 
    
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