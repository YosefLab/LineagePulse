#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++     Fit ZINB model to a data set    +++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit (zero-inflated) negative binomial model to data
#' 
#' This function fits a ZINB or NB model with variable input. 
#' It is the wrapper for the individual fits of the full and alternative 
#' models in LineagePulse.
#' 
#' For ZINB models with drop-out model estimation:
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
#' @param strDropModel (str) {"logistic_ofMu", "logistic", "none"}
#' [Default "logistic_ofMu"] Definition of drop-out model.
#' "logistic_ofMu" - include the fitted mean in the linear model
#' of the drop-out rate and use offset and matPiConstPredictors.
#' "logistic" - only use offset and matPiConstPredictors.
#' "none" - negative binomial noise model without zero-inflation.
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
#' 
#' @author David Sebastian Fischer
fitModel <- function(
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
    RELTOL_BFGS_MuDisp <- 10^(-8) 
    # optim default is sqrt(.Machine$double.eps)=1e-8
    # Lowering RELTOL_BFGS_MuDisp gives drastic run time improvements.
    # Set to 10^(-4) to maintain sensible fits 
    # without running far into saturation
    # in the objective (loglikelihood).
    # Numerical optmisation of dropout model hyperparameters
    MAXIT_BFGS_Pi <- 10000
    RELTOL_BFGS_Pi <- 10^(-8)
    ####################################################
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    
    boolFitDrop <- is.null(lsDropModel) & strDropModel != "none"
    if(!boolFitDrop) scaMaxEstimationCycles <- 1 
    # Do not need estimation cycle if drop-out model is fixed
    vecEMLogLikModel <- array(NA, scaMaxEstimationCycles)
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
    
    # c) Drop-out model: This object is also initialised with
    # NB noise but will not receive fits at any point.
    lsDropModel <- initDropModel(
        matCounts=matCounts,
        matPiConstPredictors=matPiConstPredictors,
        lsDropModel=lsDropModel,
        strDropModel=strDropModel,
        strDropFitGroup=strDropFitGroup,
        MAXIT_BFGS_Pi=MAXIT_BFGS_Pi,
        RELTOL_BFGS_Pi=RELTOL_BFGS_Pi)
    
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
                    lsFitPi <- fitPi(
                        matCounts=matCounts,
                        lsMuModel=lsMuModel,
                        lsDispModel=lsDispModel,
                        lsDropModel=lsDropModel)
                    lsDropModel$matDropoutLinModel <- 
                        lsFitPi$matDropoutLinModel
                    vecPiEstConverged <- lsFitPi$vecConverged
                    vecLL <- lsFitPi$vecLL
                })
                
                # Summary and warnings of drop-out model estimation step
                if(any(vecPiEstConverged != 0)){
                    strMessage <- paste0(
                        "Dropout estimation did not converge in ", 
                        sum(vecPiEstConverged != 0), " cases [codes: ",
                        paste(unique(vecPiEstConverged[
                            vecPiEstConverged!=0])), "].")
                    strReport <- paste0(strReport, strMessage, "\n")
                    if(boolSuperVerbose) message(strMessage)
                }
                if(any(vecPiEstConverged==1001)){
                    strMessage <- paste0(
                        "Fatal dropout estimation error in ", 
                        sum(vecPiEstConverged==1001), " cases.")
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
                lsFitMuDisp <- fitMuDisp(matCountsProc=matCounts,
                                         lsMuModel=lsMuModel,
                                         lsDispModel=lsDispModel,
                                         lsDropModel=lsDropModel,
                                         matWeights=matWeights)
            })
            lsDispModel$matDispModel <- lsFitMuDisp$matDispModel
            lsDispModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelDisp
            lsMuModel$matMuModel <- lsFitMuDisp$matMuModel
            lsMuModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelMu
            vecMuDispEstConverged <- lsFitMuDisp$vecConvergence
            
            # Summary and warnings of mu-disp model estimation step
            if(any(vecMuDispEstConverged != 0)){
                strMessage <- paste0(
                    "(Mean-) Dispersion estimation did not converge in ", 
                    sum(vecMuDispEstConverged != 0), " cases [codes: ",
                    paste(unique(vecMuDispEstConverged[
                        vecMuDispEstConverged!=0]), collapse=","), "].")
                strReport <- paste0(strReport, strMessage, "\n")
                if(boolSuperVerbose) message(strMessage)
            }
            strMessage <- paste0(
                "# ",scaIter, ".   Mean+Disp co-estimation: ",
                "ll ", scaLogLikNew, " in ",
                round(tm_mudisp["elapsed"]/60,2)," min.")
            strReport <- paste0(strReport, strMessage, "\n")
            if(boolSuperVerbose) message(strMessage)
        })
        
        # Complete iteration.
        scaLogLikOld <- scaLogLikNew
        scaLogLikNew <- sum(lsFitMuDisp$vecLL)
        vecEMLogLikModel[scaIter] <- scaLogLikNew
        scaIter <- scaIter+1
        
        # Only print summary of iteration if not already printing
        # step-wise summaries.
        strMessage <- paste0(
            "# ",scaIter-1, ".  Iteration with ",
            "ll   ", scaLogLikNew, " in ",
            round(tm_iter["elapsed"]/60,2)," min.")
        if(boolVerbose & !boolSuperVerbose) message(strMessage)
    }
    
    # Evaluate convergence
    if(all(vecMuDispEstConverged == 0) &
       scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
        boolConvergenceModel <- TRUE
    } else { 
        boolConvergenceModel <- FALSE 
    }
    
    return(list(
        lsMuModel=lsMuModel,
        lsDispModel=lsDispModel,
        lsDropModel=lsDropModel,
        boolConvergenceModel=boolConvergenceModel,
        vecEMLogLikModel=vecEMLogLikModel,
        strReport=strReport ))
}