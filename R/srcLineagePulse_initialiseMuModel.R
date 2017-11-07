initMuModel <- function(
    matCounts,
    dfAnnotation,
    vecConfoundersMu,
    vecNormConst,
    scaDFSplinesMu,
    matWeights,
    matMuModelInit,
    lsmatBatchModelInitMu,
    strMuModel,
    MAXIT_BFGS_MuDisp,
    RELTOL_BFGS_MuDisp) {
    
    scaNumGenes <- nrow(matCounts)
    scaNumCells <- ncol(matCounts)
    
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
    
    return(lsMuModel)
}