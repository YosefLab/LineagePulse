#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++     Fit ZINB model to a data set    +++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit zero-inflated negative binomial model to data
#' 
#' This function fits a ZINB model with variable input. It is the wrapper for
#' the individual fits of the full and alternative models in LineagePulse.
#' 
#' The estimation is iterative coordinate ascent over gene-wise and cell-wise
#' model if the drop-out model is not set a priori. If the drop-out model
#' is given, the estimation is a single M-like step of the iterative 
#' coordinate ascent. 
#' 
#' Convergence of iterative coordinate ascent is tracked with the the 
#' loglikelihood of the entire data matrix. 
#' Every step is a maximum likelihood estimation of the 
#' target parameters conditioned on the remaining parameter estimates. 
#' Therefore, convergence to a local optimum is guaranteed if the algorithm
#' is run until convergence. Parallelisation of each estimation step 
#' is implemented where conditional independences of parameter estimations
#' allow so. 
#' 
#' Convergence can be followed with verbose=TRUE (at each 
#' iteration) or at each step (boolSuperVerbose=TRUE).
#' 
#' To save memory, not the entire parameter matrix (genes x cells) but
#' the parmater models are stored in the objects lsMuModel, lsDispModel
#' and lsDropModel. In short, these object contain
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
#' @param matPiConstPredictors (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) [Default NULL]
#' Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param lsDropModel (list) [Default NULL]
#' Object containing description of cell-wise drop-out parameter models.
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
#' @param scaMaxEstimationCycles (integer) [Default 20] Maximum number 
#' of estimation cycles performed in fitZINB(). One cycle
#' contain one estimation of of each parameter of the 
#' zero-inflated negative binomial model as coordinate ascent.
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
#' \item lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' \item matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' \item boolConvergenceModel: (bool) 
#' Convergence status of model estimation.
#' \item vecEMLogLikModel: (numeric vector number of genes) 
#' Likelihood of model fits by iterative coordinate ascent iteration.
#' \item strReport: (str) Log of model estimation to be added to 
#' overall log.
#' }
fitZINB <- function(
    matCounts,
    dfAnnotation,
    vecConfoundersMu=NULL,
    vecConfoundersDisp=NULL,
    vecNormConst,
    scaDFSplinesMu=NULL,
    scaDFSplinesDisp=NULL,
    matWeights=NULL,
    matPiConstPredictors=NULL,
    lsDropModel=NULL,
    matMuModelInit=NULL,
    lsmatBatchModelInitMu=NULL,
    matDispModelInit=NULL,
    lsmatBatchModelInitDisp=NULL,
    strMuModel,
    strDispModel,
    strDropModel="logistic_ofMu",
    strDropFitGroup="PerCell",
    scaMaxEstimationCycles=20,
    boolVerbose=TRUE,
    boolSuperVerbose=TRUE){
    
    ####################################################
    # Internal Numerical Estimation Parameters:
    # Minimim fractional liklihood increment necessary to
    # continue EM-iterations:
    scaPrecEM <- 1-10^(-4)
    # Numerical optmisation of impulse model hyperparameters
    MAXIT_BFGS_MuDisp <- 1000 # optim default is 1000
    RELTOL_BFGS_MuDisp <- 10^(-4) 
    # optim default is sqrt(.Machine$double.eps)=1e-8
    # Lowering RELTOL_BFGS_MuDisp gives drastic run time improvements.
    # Set to 10^(-4) to maintain sensible fits 
    # without running far into saturation
    # in the objective (loglikelihood).
    # Numerical optmisation of dropout model hyperparameters
    MAXIT_BFGS_Pi <- 10000
    RELTOL_BFGS_Pi <- 10^(-4)
    ####################################################
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    boolFitDrop <- is.null(lsDropModel)
    if(!boolFitDrop) scaMaxEstimationCycles <- 1 
    # Do not need estimation cycle if drop-out model is fixed
    vecEMLogLikModel <- array(NA, scaMaxEstimationCycles)
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
    
    # c) Drop-out model: Only if this is ot given.
    # Dropout model: Initialise as offset=0 and log(mu)  parameter which
    # is forced to be negative during fitting, as -1. 
    # The parameter corresponding to log(mu) may not be initialised too 
    # close to zero, as the cost function cannot always pick up the signal 
    # in such cases, leading to an MLE with this parameter untouched.
    boolExternalDropModel <- TRUE
    if(is.null(lsDropModel)){
        boolExternalDropModel <- FALSE
        
        lsDropModel <- list(
            matDropoutLinModel=NULL,
            matPiConstPredictors=matPiConstPredictors,
            lsDropModelGlobal=list(
                strDropModel=strDropModel,
                strDropFitGroup=strDropFitGroup,
                scaNumGenes=scaNumGenes,
                scaNumCells=scaNumCells,
                MAXIT_BFGS_Pi=MAXIT_BFGS_Pi,
                RELTOL_BFGS_Pi=RELTOL_BFGS_Pi))
        # Target initialisation drop-out rate: 0.99, linear model mu
        # parameter = -1 -> solve for offset of linear model:
        scaPiTarget <- 0.99
        if(!is.null(matPiConstPredictors)){
            scaConstPredictors <- dim(matPiConstPredictors)[2]
        } else { scaConstPredictors <- 0 }
        if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic_ofMu"){
            scaPiLinModelMuParam <- -1
            scaPiLinModelOffset <- log(scaPiTarget) - log(1-scaPiTarget) - 
                scaPiLinModelMuParam*log(min(vecMuModelInit, na.rm=TRUE))
            lsDropModel$matDropoutLinModel <- cbind(
                rep(scaPiLinModelOffset, scaNumCells), 
                rep(scaPiLinModelMuParam, scaNumCells),
                matrix(0, nrow=scaNumCells, ncol=scaConstPredictors))
        } else if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic") {
            scaPiLinModelOffset <- log(scaPiTarget) - log(1-scaPiTarget)
            lsDropModel$matDropoutLinModel <- cbind(
                rep(scaPiLinModelOffset, scaNumCells), 
                matrix(0, nrow=scaNumCells, ncol=scaConstPredictors))
        }
    }
    
    # Evaluate initialisation loglikelihood
    scaLogLikNew <- sum(evalLogLikMatrix(
        matCounts=matCounts,
        lsMuModel=lsMuModel,
        lsDispModel=lsDispModel, 
        lsDropModel=lsDropModel,
        matWeights=matWeights ))
    strMessage <- paste0("#  .   Initialisation: ",
                         "ll ", scaLogLikNew)
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) message(strMessage)
    
    ### Iteration
    scaIter <- 1
    scaLogLikOld <- NA
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & 
                          scaIter <= scaMaxEstimationCycles)){
        # If drop-out model was supplied (boolFitDrop==FALSE), 
        # drop-out model estimation is skipped and 
        # the mean-dispersion model is estimated 
        # conditioned on the supplied drpo-out model. 
        # No iteration is required: 
        # the loop is exited via tha scaIter <= scaMaxEstimationCycles
        # condition (scaMaxEstimationCycles is set to 1 if boolFitDrop==FALSE).
        tm_iter <- system.time({
            if(boolFitDrop){
                #####  1. Cell-wise parameter estimation
                # Dropout rate
                tm_pi <- system.time({
                    lsFitPi <- fitZINBPi(
                        matCounts=matCounts,
                        lsMuModel=lsMuModel,
                        lsDispModel=lsDispModel,
                        lsDropModel=lsDropModel)
                    lsDropModel$matDropoutLinModel <- 
                        lsFitPi$matDropoutLinModel
                    vecboolPiEstConverged <- lsFitPi$vecboolConverged
                    vecLL <- lsFitPi$vecLL
                })
                colnames(lsDropModel$matDropoutLinModel) <- NULL 
                # Want this so that column names dont grow to par.par.par...
                
                if(any(vecboolPiEstConverged != 0)){
                    strMessage <- paste0(
                        "Dropout estimation did not converge in ", 
                        sum(vecboolPiEstConverged), " cases [codes: ",
                        paste(unique(vecboolPiEstConverged[
                            vecboolPiEstConverged!=0])), "].")
                    strReport <- paste0(strReport, strMessage, "\n")
                    if(boolSuperVerbose) message(strMessage)
                }
                if(any(vecboolPiEstConverged==1001)){
                    strMessage <- paste0(
                        "Fatal dropout estimation error in ", 
                        sum(vecboolPiEstConverged==1001), " cases.")
                    strReport <- paste0(strReport, strMessage, "\n")
                    if(boolSuperVerbose) message(strMessage)
                }
                strMessage <- paste0(
                    "# ",scaIter,".   Drop-out estimation: ",
                    "ll     ", sum(vecLL), " in ",
                    round(tm_pi["elapsed"]/60,2)," min.")
                strReport <- paste0(strReport, strMessage, "\n")
                if(boolSuperVerbose) message(strMessage)
            }
            
            ##### 2. Gene-wise parameter estimation:
            # Estimate mean and dispersion parameters simultaneously.
            # a/b) Negative binomial mean AND dispersion parameter.
            tm_mudisp <- system.time({
                lsFitMuDisp <- fitZINBMuDisp(matCountsProc=matCounts,
                                             lsMuModel=lsMuModel,
                                             lsDispModel=lsDispModel,
                                             lsDropModel=lsDropModel,
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
            vecboolMuEstConverged <- lsFitMuDisp$vecConvergence
            vecboolDispEstConverged <- lsFitMuDisp$vecConvergence
            
            # Evaluate Likelihood
            scaLogLikOld <- scaLogLikNew
            scaLogLikNew <- sum(lsFitMuDisp$vecLL)
        })
        
        # Iteration complete
        if(any(vecboolDispEstConverged != 0)){
            strMessage <- paste0(
                "(Mean-) Dispersion estimation did not converge in ", 
                sum(vecboolDispEstConverged[vecboolDispEstConverged>0]), 
                " cases [codes: ",
                paste(unique(vecboolDispEstConverged[
                    vecboolDispEstConverged!=0]), collapse=","), "].")
            strReport <- paste0(strReport, strMessage, "\n")
            if(boolSuperVerbose) message(strMessage)
        }
        
        strMessage <- paste0(
            "# ",scaIter, ".   Mean+Disp co-estimation: ",
            "ll ", scaLogLikNew, " in ",
            round(tm_mudisp["elapsed"]/60,2)," min.")
        strReport <- paste0(strReport, strMessage, "\n")
        if(boolSuperVerbose) message(strMessage)
        
        strMessage <- paste0(
            "# ",scaIter, ".  Iteration with ",
            "ll   ", scaLogLikNew, " in ",
            round(tm_iter["elapsed"]/60,2)," min.")
        if(boolVerbose & !boolSuperVerbose) message(strMessage)
        
        vecEMLogLikModel[scaIter] <- scaLogLikNew
        scaIter <- scaIter+1
    }
    
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
    if(!boolExternalDropModel) {
        rownames(lsDropModel$matDropoutLinModel) <- colnames(matCounts)
    }
    # Name model matrix columns
    if(strMuModel=="groups") {
        colnames(lsMuModel$matMuModel) <- unique(dfAnnotation$groups)
    } else if(strMuModel=="impulse") {
        colnames(lsMuModel$matMuModel) <- 
            c("beta1", "beta2", "h0", "h1", "h2", "t1", "t2")
    }
    
    # Evaluate convergence
    if(all(as.logical(vecboolDispEstConverged)) &
       all(as.logical(vecboolMuEstConverged)) &
       scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
        boolConvergenceModel <- TRUE
    } else { boolConvergenceModel <- FALSE }
    
    return(list(
        lsMuModel=lsMuModel,
        lsDispModel=lsDispModel,
        lsDropModel=lsDropModel,
        boolConvergenceModel=boolConvergenceModel,
        vecEMLogLikModel=vecEMLogLikModel,
        strReport=strReport ))
}