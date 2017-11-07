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