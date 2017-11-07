#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++     Dropout model container object    +++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Initialise drop-out model container object
#' 
#' Either use supplied fits from previous fitting or initialise 
#' from count data.
#' 
#' @seealso Called by \code{fitModel}. 
#' 
#' @param matCounts (matrix genes x cells)
#' Count data of all cells, unobserved entries are NA.
#' @param matPiConstPredictors (numeric matrix genes x number of constant
#' gene-wise drop-out predictors) [Default NULL]
#' Predictors for logistic drop-out 
#' fit other than offset and mean parameter (i.e. parameters which
#' are constant for all observations in a gene and externally supplied.)
#' Is null if no constant predictors are supplied.
#' @param lsDropModel (list) [Default NULL]
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
#' @param MAXIT_BFGS_Pi (sca)
#' Maximum number of iterations in BFGS estimation of dropout models.
#' This is a control parameter to optim().
#' @param RELTOL_BFGS_Pi (sca) 
#' Relative tolerance of BFGS estimation of dropout models.
#' This is a control parameter to optim().
#' 
#' @return lsDropModel (list)
#' Initialisation of drop-out model object.
#' 
#' @author David Sebastian Fischer
initDropModel <- function(
    matCounts,
    matPiConstPredictors,
    lsDropModel,
    strDropModel,
    strDropFitGroup,
    MAXIT_BFGS_Pi,
    RELTOL_BFGS_Pi){
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    
    # Dropout model: Initialise as offset=0 and log(mu)  parameter which
    # is forced to be negative during fitting, as -1. 
    # The parameter corresponding to log(mu) may not be initialised too 
    # close to zero, as the cost function cannot always pick up the signal 
    # in such cases, leading to an MLE with this parameter untouched.
    if(is.null(lsDropModel) & strDropModel != "none"){
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
    } else if(is.null(lsDropModel) & strDropModel == "none"){
        # This corresponds to a NB noise model, do not need initialisation
        # of model matrices only set global parameters which are used to
        # evaluate type of noise model.
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
    }
    # Note that if: !is.null(lsDropModel) & strDropModel != "none",
    # The prefit lsDropModel is returned.
    
    return(lsDropModel)
}