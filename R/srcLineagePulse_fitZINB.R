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
    strDropModel,
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
    RELTOL_BFGS_MuDisp <- 10^(-4) # optim default is sqrt(.Machine$double.eps)=1e-8
    # Lowering RELTOL_BFGS_MuDisp gives drastic run time improvements.
    # Set to 10^(-4) to maintain sensible fits without running far into saturation
    # in the objective (loglikelihood).
    # Numerical optmisation of dropout model hyperparameters
    MAXIT_BFGS_Pi <- 10000
    RELTOL_BFGS_Pi <- 10^(-4)
    ####################################################
    
    scaNumGenes <- dim(matCounts)[1]
    scaNumCells <- dim(matCounts)[2]  
    boolFitDrop <- is.null(lsDropModel)
    if(!boolFitDrop) scaMaxEstimationCycles <- 1 # Do not need estimation cycle if drop-out model is fixed
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
            vecNBatches=NULL, # can be computed as sapply(lsvecidxBatchAssign, max), keep for speed
            lsvecidxBatchAssign=NULL, # object called to determine assignment of cell to batch
            lsvecBatchUnique=NULL, # Kept so that idx vector can be associated with batch names
            MAXIT_BFGS_MuDisp=MAXIT_BFGS_MuDisp,
            RELTOL_BFGS_MuDisp=RELTOL_BFGS_MuDisp) )
    # Initialise mean model parameters
    if(is.null(matMuModelInit)){
        vecMuModelInit <- apply(matCounts, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
        vecMuModelInit[vecMuModelInit < 10^(-10)] <- 10^(-10)
        if(strMuModel=="constant"){
            lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=1, byrow=FALSE)
        } else if(strMuModel=="impulse"){
            lsMuModel$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=7)
            lsMuModel$matMuModel[,c(3:5)] <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=3, byrow=FALSE)
        } else if(strMuModel=="splines"){
            lsMuModel$matMuModel <- matrix(1, nrow=scaNumGenes, ncol=scaDFSplinesMu)
        } else if(strMuModel=="groups"){
            lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=length(unique(dfAnnotation$groups)), byrow=FALSE)
        } else if(strMuModel=="MM"){
            lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=dim(matWeights)[2], byrow=FALSE)
        } else  if(strMuModel=="windows"){
            lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
        } else {
            stop(paste0("ERROR fitZINB(): strMuModel=", strMuModel, " not recognised."))
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
        lsMuModel$lsMuModelGlobal$vecidxGroups <- match(vecGroupsMu, unique(vecGroupsMu))
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
            lsvecBatchAssign <- lapply(vecConfoundersMu, function(confounder) dfAnnotation[[confounder]] )
            lsMuModel$lsmatBatchModel <- lapply(lsvecBatchAssign, function(vecBatchAssign){
                matrix(1, nrow=scaNumGenes,
                       ncol=length(unique(vecBatchAssign)) ) # Initialise batch correction factors to 1
            })
        } else {
            lsMuModel$lsmatBatchModel <- list(matrix(NA, nrow=scaNumGenes, ncol=1 ))
        }
    } else {
        lsMuModel$lsmatBatchModel <- lsmatBatchModelInitMu
    }
    # Add global batch parameters
    if(!is.null(vecConfoundersMu)){
        lsMuModel$lsMuModelGlobal$scaNConfounders <- length(vecConfoundersMu)
        lsvecBatchesMu <- lapply(vecConfoundersMu, function(confounder) dfAnnotation[[confounder]] )
        lsMuModel$lsMuModelGlobal$lsvecBatchUnique <-lapply(lsvecBatchesMu, function(vecBatchAssign) unique(vecBatchAssign) )
        lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign <- lapply(lsvecBatchesMu, function(vecBatchAssign) match(vecBatchAssign, unique(vecBatchAssign)) )
        lsMuModel$lsMuModelGlobal$vecNBatches <- sapply(lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign, max)
    } 
    lsMuModel$lsMuModelGlobal$scaDegFreedom <- dim(lsMuModel$matMuModel)[2] + # Mu model
        sum(sapply(lsMuModel$lsmatBatchModel, function(mat) dim(mat)[2]-1 )) # Batch correction model
    
    # b) Dispersion model
    # Dispersions: Low dispersion factor yielding high variance which makes
    # cost function screening easy in the first iteration.
    lsDispModel <- list(matDispModel=NULL,
                        lsmatBatchModel=NULL,
                        lsDispModelGlobal=list(
                            scaDegFreedom=NULL, # general parameters
                            strDispModel=strDispModel,
                            scaNumCells=scaNumCells,
                            vecPseudotime=dfAnnotation$pseudotime, # continuous models
                            scaNSplines=NULL, ## spline models
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
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=1, byrow=FALSE)
        } else if(strDispModel=="MM"){
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=dim(matWeights)[2])
        } else if(strDispModel=="splines"){
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=scaDFSplinesDisp)
        } else if(strDispModel=="groups"){
            lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=length(unique(dfAnnotation$groups)))
        } else {
            stop(paste0("ERROR fitZINB(): strDispModel=", strDispModel, " not recognised."))
        } 
    } else {
        lsDispModel$matDispModel <- matDispModelInit
    }
    if(is.null(lsmatBatchModelInitDisp)){
        # Initialise batch model parameters
        if(!is.null(vecConfoundersDisp)){
            lsvecBatchAssign <- lapply(vecConfoundersDisp, function(confounder) dfAnnotation[[confounder]] )
            lsDispModel$lsmatBatchModel <- lapply(lsvecBatchAssign, function(vecBatchAssign){
                matrix(1, nrow=scaNumGenes,
                       ncol=length(unique(vecBatchAssign)) ) # Initialise batch correction factors to 1
            })
        } else {
            lsDispModel$lsmatBatchModel <- list(matrix(NA, nrow=scaNumGenes, ncol=1 ))
        }
    } else {
        lsDispModel$lsmatBatchModel <- lsmatBatchModelInitDisp
    }
    if(strDispModel=="groups"){
        vecGroupsDisp <- as.vector(dfAnnotation$groups)
        lsDispModel$lsDispModelGlobal$scaNGroups <- length(unique(vecGroupsDisp))
        lsDispModel$lsDispModelGlobal$vecGroups <- vecGroupsDisp
        lsDispModel$lsDispModelGlobal$vecidxGroups <- match(vecGroupsDisp, unique(vecGroupsDisp))
    }
    if(strDispModel=="splines"){
        lsDispModel$lsDispModelGlobal$vecTimepoints <- 
            sort(unique(dfAnnotation$pseudotime ))
        lsDispModel$lsDispModelGlobal$vecindTimepointAssign <- match(
            dfAnnotation$pseudotime, lsDispModel$lsDispModelGlobal$vecTimepoints)
        lsDispModel$lsDispModelGlobal$scaNSplines <- scaDFSplinesDisp
        lsDispModel$lsDispModelGlobal$matSplineBasis <- 
            ns(x = dfAnnotation$pseudotime, 
               df = scaDFSplinesDisp, intercept = TRUE)[,,drop=FALSE]
    }
    if(strMuModel=="MM"){ # Not a typo! Also have to fill this in if only strMuModel is MM
        lsDispModel$lsDispModelGlobal$scaNMix <- dim(matWeights)[2]
    }
    if(!is.null(vecConfoundersDisp)) {
        lsDispModel$lsDispModelGlobal$scaNConfounders <- length(vecConfoundersDisp)
        lsvecBatchesDisp <- lapply(vecConfoundersDisp, function(confounder) as.vector(dfAnnotation[[confounder]]) )
        lsDispModel$lsDispModelGlobal$lsvecBatchUnique <-lapply(lsvecBatchesDisp, function(vecBatchAssign) unique(vecBatchAssign) )
        lsDispModel$lsDispModelGlobal$lsvecidxBatchAssign <- lapply(lsvecBatchesDisp, function(vecBatchAssign) match(vecBatchAssign, unique(vecBatchAssign)) )
        lsDispModel$lsDispModelGlobal$vecNBatches <- sapply(lsDispModel$lsDispModelGlobal$lsvecidxBatchAssign, max)
    }
    lsDispModel$lsDispModelGlobal$scaDegFreedom <- dim(lsDispModel$matDispModel)[2] + # Disp model
        sum(lsDispModel$lsDispModelGlobal$vecNBatches) - 
        length(lsDispModel$lsDispModelGlobal$vecNBatches) # Batch correction model
    
    # c) Drop-out model: Only if this is ot given.
    # Dropout model: Initialise as offset=0 and log(mu)  parameter which
    # is forced to be negative during fitting, as -1. The parameter corresponding
    # to log(mu) may not be initialised too close to zero, as the cost function 
    # cannot always pick up the signal in such cases, leading to an MLE with this 
    # parameter untouched.
    boolExternalDropModel <- TRUE
    if(is.null(lsDropModel)){
        boolExternalDropModel <- FALSE
        
        lsDropModel <- list(matDropoutLinModel=NULL,
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
        } else if(lsDropModel$lsDropModelGlobal$strDropModel=="logistic"){
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
                         "ll          ", scaLogLikNew)
    strReport <- paste0(strReport, strMessage, "\n")
    if(boolVerbose) print(strMessage)
    
    ### Iteration
    scaIter <- 1
    scaLogLikOld <- NA
    while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxEstimationCycles)){
        # If drop-out model was supplied (boolFitDrop==FALSE), drop-out model
        # estimation is skipped and the mean-dispersion model is estimated 
        # conditioned on the supplied drpo-out model. No iteration is required
        # and the while loop is excited via tha scaIter <= scaMaxEstimationCycles
        # condition (scaMaxEstimationCycles is set to 1 if boolFitDrop==FALSE).
        tm_iter <- system.time({
            if(boolFitDrop){
                #####  1. Cell-wise parameter estimation
                # Dropout rate
                tm_pi <- system.time({
                    lsFitPi <- fitZINBPi(matCounts=matCounts,
                                         lsMuModel=lsMuModel,
                                         lsDispModel=lsDispModel,
                                         lsDropModel=lsDropModel)
                    lsDropModel$matDropoutLinModel <- lsFitPi$matDropoutLinModel
                    vecboolPiEstConverged <- lsFitPi$vecboolConverged
                    vecLL <- lsFitPi$vecLL
                })
                colnames(lsDropModel$matDropoutLinModel) <- NULL # Want this so that column names dont grow to par.par.par...
                
                if(any(vecboolPiEstConverged != 0)){
                    strMessage <- paste0("Dropout estimation did not converge in ", 
                                         sum(vecboolPiEstConverged), " cases [codes: ",
                                         paste(unique(vecboolPiEstConverged[vecboolPiEstConverged!=0])), "].")
                    strReport <- paste0(strReport, strMessage, "\n")
                    if(boolSuperVerbose) print(strMessage)
                }
                if(any(vecboolPiEstConverged==1001)){
                    strMessage <- paste0("Fatal dropout estimation error in ", 
                                         sum(vecboolPiEstConverged==1001), " cases.")
                    strReport <- paste0(strReport, strMessage, "\n")
                    if(boolSuperVerbose) print(strMessage)
                }
                strMessage <- paste0("# ",scaIter,".   Drop-out estimation: ",
                                     "ll     ", sum(vecLL), " in ",
                                     round(tm_pi["elapsed"]/60,2)," min.")
                strReport <- paste0(strReport, strMessage, "\n")
                if(boolSuperVerbose) print(strMessage)
            }
            
            # reference for values from optim
            #vecLogLikRef <- evalLogLikMatrix(matCounts=matCounts,
            #                                 lsMuModel=lsMuModel,
            #                                 lsDispModel=lsDispModel, 
            #                                 lsDropModel=lsDropModel,
            #                                 matWeights=matWeights )
            #strMessage <- paste0("# ",scaIter, ".   Dropout co-estimation complete: ",
            #                     "reference ll of    ", sum(vecLogLikRef), ".")
            #print(strMessage)
            
            ##### 2. Gene-wise parameter estimation:
            # Estimate mean and dispersion parameters simultaneously.
            # a/b) Negative binomial mean AND dispersion parameter.
            tm_mudisp <- system.time({
                lsFitMuDisp <- fitZINBMuDisp(matCounts=matCounts,
                                             lsMuModel=lsMuModel,
                                             lsDispModel=lsDispModel,
                                             lsDropModel=lsDropModel,
                                             matWeights=matWeights)
            })
            lsDispModel$matDispModel <- lsFitMuDisp$matDispModel
            colnames(lsDispModel$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
            lsDispModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelDisp
            if(!is.null(lsDispModel$lsmatBatchModel)) {
                for(i in seq(1, length(lsDispModel$lsmatBatchModel))) colnames(lsDispModel$lsmatBatchModel[[i]]) <- NULL # Need this so that column names dont grow to par.par.par...
            }
            lsMuModel$matMuModel <- lsFitMuDisp$matMuModel
            colnames(lsMuModel$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
            lsMuModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModelMu
            if(!is.null(lsMuModel$lsmatBatchModel)) {
                for(i in seq(1, length(lsMuModel$lsmatBatchModel))) colnames(lsMuModel$lsmatBatchModel[[i]]) <- NULL # Need this so that column names dont grow to par.par.par...
            }
            vecboolMuEstConverged <- lsFitMuDisp$vecConvergence
            vecboolDispEstConverged <- lsFitMuDisp$vecConvergence
            
            # reference for values from optim
            #vecLogLikRef <- evalLogLikMatrix(matCounts=matCounts,
            #                                 lsMuModel=lsMuModel,
            #                                 lsDispModel=lsDispModel, 
            #                                 lsDropModel=lsDropModel,
            #                                 matWeights=matWeights )
            #strMessage <- paste0("# ",scaIter, ".   Mean+Disp co-estimation: ",
            #                     "reference ll  ", sum(vecLogLikRef), ".")
            #print(strMessage)
            
            # Evaluate Likelihood
            scaLogLikOld <- scaLogLikNew
            scaLogLikNew <- sum(lsFitMuDisp$vecLL)
        })
        
        # Iteration complete
        if(any(vecboolDispEstConverged != 0)){
            strMessage <- paste0("(Mean-) Dispersion estimation did not converge in ", 
                                 sum(vecboolDispEstConverged[vecboolDispEstConverged>0]), " cases [codes: ",
                                 paste(unique(vecboolDispEstConverged[vecboolDispEstConverged!=0]), collapse=","), "].")
            strReport <- paste0(strReport, strMessage, "\n")
            if(boolSuperVerbose) print(strMessage)
        }
        
        strMessage <- paste0("# ",scaIter, ".   Mean+Disp co-estimation: ",
                             "ll ", scaLogLikNew, " in ",
                             round(tm_mudisp["elapsed"]/60,2)," min.")
        strReport <- paste0(strReport, strMessage, "\n")
        if(boolSuperVerbose) print(strMessage)
        
        strMessage <- paste0("# ",scaIter, ".  Iteration with ",
                             "ll   ", scaLogLikNew, " in ",
                             round(tm_iter["elapsed"]/60,2)," min.")
        if(boolVerbose & !boolSuperVerbose) print(strMessage)
        
        vecEMLogLikModel[scaIter] <- scaLogLikNew
        scaIter <- scaIter+1
    }
    
    # Name model matrix rows
    rownames(lsMuModel$matMuModel) <- rownames(matCounts)
    if(!is.null(lsMuModel$lsmatBatchModel)) {
        for(i in seq(1, length(lsMuModel$lsmatBatchModel))) rownames(lsMuModel$lsmatBatchModel[[i]]) <- rownames(matCounts)
    }
    rownames(lsDispModel$matDispModel) <- rownames(matCounts)
    if(!is.null(lsDispModel$lsmatBatchModel)) {
        for(i in seq(1, length(lsDispModel$lsmatBatchModel))) rownames(lsDispModel$lsmatBatchModel[[i]]) <- rownames(matCounts)
    }
    if(!boolExternalDropModel) rownames(lsDropModel$matDropoutLinModel) <- colnames(matCounts)
    # Name model matrix columns
    if(strMuModel=="clusters") colnames(lsMuModel$matMuModel) <- unique(dfAnnotation$clusters)
    else if(strMuModel=="windows") colnames(lsMuModel$matMuModel) <- colnames(matCounts)
    else if(strMuModel=="impulse") colnames(lsMuModel$matMuModel) <- c("beta1", "beta2", "h0", "h1", "h2", "t1", "t2")
    
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