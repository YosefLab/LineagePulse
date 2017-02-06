fitZINB <- function(matCounts,
										dfAnnotation,
										vecConfounders,
										vecNormConst,
										matWeights=NULL,
										matPiConstPredictors=NULL,
										scaWindowRadius=NULL,
										boolVecWindowsAsBFGS=FALSE,
										lsDropModel=NULL,
										matMuModelInit=NULL,
										lsmatBatchModelInit=NULL,
										matDispModelInit=NULL,
										strMuModel,
										strDispModel,
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
	lsMuModel <- list(matMuModel=NA,
										lsmatBatchModel=NA,
										lsMuModelGlobal=list(scaDegFreedom=NA,
										                     strMuModel=strMuModel,
																				 vecNormConst=vecNormConst,
																				 scaNumCells=scaNumCells,
																				 vecPseudotime=dfAnnotation$pseudotime,
																				 vecidxClusterAssign=NA,
																				 scaClusterK=NA,
																				 vecConfounders=vecConfounders,
																				 lsvecidxBatchAssign=NULL, # object called to determine assignment of cell to batch
																				 lsvecidxBatchUnique=NULL, # Can be computer from lsvecidxBatchAssign but kept for speed
																				 lsvecBatchUnique=NULL, # Kept so that idx vector can be associated with batch names
																				 boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
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
		} else if(strMuModel=="clusters"){
			lsMuModel$matMuModel <- matrix(vecMuModelInit, nrow=scaNumGenes, ncol=length(unique(dfAnnotation$clusters)), byrow=FALSE)
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
	if(strMuModel=="clusters"){
	  lsMuModel$lsMuModelGlobal$vecidxClusterAssign <- match(dfAnnotation$clusters, unique(dfAnnotation$clusters))
	  lsMuModel$lsMuModelGlobal$scaClusterK <- length(unique(dfAnnotation$clusters))
	}
	if(is.null(lsmatBatchModelInit)){
	  # Initialise batch model parameters
	  if(!is.null(vecConfounders)){
	    lsvecBatchAssign <- lapply(vecConfounders, function(confounder) dfAnnotation[[confounder]] )
	    lsMuModel$lsmatBatchModel <- lapply(lsvecBatchAssign, function(vecBatchAssign){
	      matrix(1, nrow=scaNumGenes,
	             ncol=length(unique(vecBatchAssign)) ) # Initialise batch correction factors to 1
	    })
	  } else {
	    lsvecBatchAssign <- list(rep(1,scaNumCells))
	    lsMuModel$lsmatBatchModel <- lapply(lsvecBatchAssign, function(vecBatchAssign){
	      matrix(NA, nrow=scaNumGenes, ncol=1 )
	    })
	  }
	} else {
	  lsMuModel$lsmatBatchModel <- lsmatBatchModelInit
	}
	# Add global batch parameters
	if(!is.null(vecConfounders)){
	  lsvecBatchAssign <- lapply(vecConfounders, function(confounder) dfAnnotation[[confounder]] )
	  lsMuModel$lsMuModelGlobal$lsvecBatchUnique <-lapply(lsvecBatchAssign, function(vecBatchAssign) unique(vecBatchAssign) )
	  lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign <- lapply(lsvecBatchAssign, function(vecBatchAssign) match(vecBatchAssign, unique(vecBatchAssign)) )
	  lsMuModel$lsMuModelGlobal$lsvecidxBatchUnique <- lapply(lsvecBatchAssign, function(vecBatchAssign) seq(1,length(unique(vecBatchAssign))) )
	} else {
	  lsvecBatchAssign <- list(rep(1,scaNumCells))
	  lsMuModel$lsMuModelGlobal$lsvecBatchUnique <- NA
	  lsMuModel$lsMuModelGlobal$lsvecidxBatchAssign <- lapply(lsvecBatchAssign, function(vecBatchAssign) match(vecBatchAssign, unique(vecBatchAssign)) )
	  lsMuModel$lsMuModelGlobal$lsvecidxBatchUnique <- lapply(lsvecBatchAssign, function(vecBatchAssign) seq(1,length(unique(vecBatchAssign))) )
	}
	lsMuModel$lsMuModelGlobal$scaDegFreedom <- dim(lsMuModel$matMuModel)[2] + # Mu model
	  sum(sapply(lsMuModel$lsmatBatchModel, function(mat) dim(mat)[2]-1 )) # Batch correction model
	
	# b) Dispersion model
	# Dispersions: Low dispersion factor yielding high variance which makes
	# cost function screening easy in the first iteration.
	lsDispModel <- list(matDispModel=NA,
											lsDispModelGlobal=list(scaDegFreedom=NA,
											                       strDispModel=strDispModel,
																						 scaNumCells=scaNumCells,
																						 vecPseudotime=dfAnnotation$pseudotime,
																						 vecClusterAssign=dfAnnotation$clusters) )
	if(is.null(matDispModelInit)){
		if(strDispModel=="constant"){
			lsDispModel$matDispModel <- matrix(1, nrow=scaNumGenes, ncol=1, byrow=FALSE)
		} else {
			stop(paste0("ERROR fitZINB(): strDispModel=", strDispModel, " not recognised."))
		} 
	} else {
		lsDispModel$matDispModel <- matDispModelInit
	}
	lsDispModel$lsDispModelGlobal$scaDegFreedom <- dim(lsDispModel$matDispModel)[2] # Disp model
	
	
	# c) Drop-out model: Only if this is ot given.
	# Dropout model: Initialise as offset=0 and log(mu)  parameter which
	# is forced to be negative during fitting, as -1. The parameter corresponding
	# to log(mu) may not be initialised too close to zero, as the cost function 
	# cannot always pick up the signal in such cases, leading to an MLE with this 
	# parameter untouched.
	boolExternalDropModel <- TRUE
	if(is.null(lsDropModel)){
	  boolExternalDropModel <- FALSE
		lsDropModel <- list(matDropoutLinModel=NA,
												matPiConstPredictors=matPiConstPredictors,
												lsPiOptimHyperparam=list(
													MAXIT_BFGS_Pi=MAXIT_BFGS_Pi,
													RELTOL_BFGS_Pi=RELTOL_BFGS_Pi))
		# Target initialisation drop-out rate: 0.99, linear model mu
		# parameter = -1 -> solve for offset of linear model:
		scaPiTarget <- 0.99
		scaPiLinModelMuParam <- -1
		scaPiLinModelOffset <- log(scaPiTarget) - log(1-scaPiTarget) - 
			scaPiLinModelMuParam*log(min(vecMuModelInit, na.rm=TRUE))
		scaPredictors <- 2
		if(!is.null(matPiConstPredictors)){
			scaPredictors <- scaPredictors + dim(matPiConstPredictors)[2]
		}
		lsDropModel$matDropoutLinModel <- cbind(
			rep(scaPiLinModelOffset, scaNumCells), 
			rep(scaPiLinModelMuParam, scaNumCells),
			matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
	}
	
	# Evaluate initialisation loglikelihood
	scaLogLikNew <- sum(evalLogLikMatrix(matCounts=matCounts,
	                                       lsMuModel=lsMuModel,
	                                       lsDispModel=lsDispModel, 
	                                       lsDropModel=lsDropModel,
	                                       scaWindowRadius=scaWindowRadius,
	                                       matWeights=matWeights ))
	strMessage <- paste0("#  .   Initialisation complete: ",
	                     "log likelihood of         ", scaLogLikNew)
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
				# Drop-out estimation is independent between cells and can be parallelised.
				tm_pi <- system.time({
					lsFitsPi <- bplapply(seq(1, scaNumCells), function(cell){
						if(!is.null(scaWindowRadius)){
							scaindIntervalStart <- max(1,cell-scaWindowRadius)
							scaindIntervalEnd <- min(scaNumCells,cell+scaWindowRadius)
							vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
						} else {
							vecInterval <- cell
						}
						
						lsFitPi <- fitPiZINB(
							vecCounts=matCounts[,cell],
							lsMuModel=lsMuModel,
							lsDispModel=lsDispModel,
							lsDropModel=lsDropModel,
							vecidxInterval=vecInterval,
							idxTarget=cell)
						return(lsFitPi)
					})
					lsDropModel$matDropoutLinModel <-  do.call(rbind, lapply(lsFitsPi, function(cell) cell$vecLinModel))
					vecboolPiEstConverged <- sapply(lsFitsPi, function(cell) cell$scaConvergence)
					vecLL <- sapply(lsFitsPi, function(cell) cell$scaLL)  
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
				strMessage <- paste0("# ",scaIter,".   Drop-out estimation complete: ",
				                     "loglikelihood of     ", sum(vecLL), " in ",
				                     round(tm_pi["elapsed"]/60,2)," min.")
				strReport <- paste0(strReport, strMessage, "\n")
				if(boolSuperVerbose) print(strMessage)
			}
			
			##### 2. Gene-wise parameter estimation:
			# Estimate mean and dispersion parameters simultaneously.
			# a/b) Negative binomial mean AND dispersion parameter.
			tm_mudisp <- system.time({
				lsFitMuDisp <- fitZINBMuDisp(matCounts=matCounts,
																		 lsMuModel=lsMuModel,
																		 lsDispModel=lsDispModel,
																		 lsDropModel=lsDropModel,
																		 scaWindowRadius=scaWindowRadius,
																		 matWeights=matWeights)
			})
			lsDispModel$matDispModel <- lsFitMuDisp$matDispModel
			colnames(lsDispModel$matDispModel) <- NULL # Need this so that column names dont grow to par.par.par...
			lsMuModel$matMuModel <- lsFitMuDisp$matMuModel
			colnames(lsMuModel$matMuModel) <- NULL # Need this so that column names dont grow to par.par.par...
			lsMuModel$lsmatBatchModel <- lsFitMuDisp$lsmatBatchModel
			for(i in seq(1, length(lsMuModel$lsmatBatchModel))) colnames(lsMuModel$lsmatBatchModel[[i]]) <- NULL # Need this so that column names dont grow to par.par.par...
			
			vecboolMuEstConverged <- lsFitMuDisp$vecConvergence
			vecboolDispEstConverged <- lsFitMuDisp$vecConvergence
			
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
		
		strMessage <- paste0("# ",scaIter, ".   Mean+Disp co-estimation complete: ",
		             "loglikelihood of ", scaLogLikNew, " in ",
		             round(tm_mudisp["elapsed"]/60,2)," min.")
		strReport <- paste0(strReport, strMessage, "\n")
		if(boolSuperVerbose) print(strMessage)
		
		strMessage <- paste0("# ",scaIter, ". complete with ",
		                     "log likelihood of   ", scaLogLikNew, " in ",
		                     round(tm_iter["elapsed"]/60,2)," min.")
		if(boolVerbose & !boolSuperVerbose) print(strMessage)
		
		vecEMLogLikModel[scaIter] <- scaLogLikNew
		scaIter <- scaIter+1
	}
	
	# Name model matrix rows
	rownames(lsMuModel$matMuModel) <- rownames(matCounts)
	for(i in seq(1, length(lsMuModel$lsmatBatchModel))) rownames(lsMuModel$lsmatBatchModel[[i]]) <- rownames(matCounts)
	rownames(lsDispModel$matDispModel) <- rownames(matCounts)
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