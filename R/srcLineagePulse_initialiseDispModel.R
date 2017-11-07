#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++     Dispersion model container object    ++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Initialise dispersion model container object
#' 
#' Either use supplied fits from previous fitting or initialise 
#' from count data.
#' 
#' @seealso Called by \code{fitModel}. 
#' 
#' @param matCounts (matrix genes x cells)
#' Count data of all cells, unobserved entries are NA.
#' @param dfAnnotation (data frame cells x meta characteristics)
#' Annotation table which contains meta data on cells.
#' @param vecConfoundersDisp 
#' (vector of strings number of confounders on dispersion)
#' [Default NULL] Confounders to correct for in dispersion batch
#' correction model, must be subset of column names of
#' dfAnnotation which describe condounding variables.
#' @param scaDFSplinesDisp (sca) [Default NULL] 
#' If strDispModelFull=="splines" or strDispModelRed=="splines", 
#' the degrees of freedom of the natural
#' cubic spline to be used as a dispersion parameter model.
#' @param matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
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
#' @param MAXIT_BFGS_MuDisp (sca)
#' Maximum number of iterations in BFGS estimation of Mu/Disp models.
#' This is a control parameter to optim().
#' @param RELTOL_BFGS_MuDisp (sca) 
#' Relative tolerance of BFGS estimation of Mu/Disp models.
#' This is a control parameter to optim().
#' 
#' @return lsDispModel (list)
#' Initialisation of dispersion model object.
#' 
#' @author David Sebastian Fischer
initDispModel <- function(
    matCounts,
    dfAnnotation,
    vecConfoundersDisp,
    scaDFSplinesDisp,
    matWeights,
    matDispModelInit,
    lsmatBatchModelInitDisp,
    strDispModel,
    strMuModel,
    MAXIT_BFGS_MuDisp,
    RELTOL_BFGS_MuDisp) {
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    
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
    
    return(lsDispModel)
}