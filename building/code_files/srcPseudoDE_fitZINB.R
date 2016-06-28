#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for dispersion fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial overdispersion on single gene given
#' the drop-out rate and negative binomial mean parameter. The
#' dispersion parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' 
#' @aliases evalLogLikHurdleNB_com
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecY: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuEst: (vector number of cells) Negative binomial
#'    mean parameter estimate of clusters to which cells
#'    belong.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell. 
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikDispNB <- function(scaTheta,
  vecY,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved, 
  vecboolZero){ 
  
  # Log linker function to fit positive dispersions
  scaDispEst <- exp(scaTheta)
  
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  scaDispEst[scaDispEst < .Machine$double.eps] <- .Machine$double.eps
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  scaDispEst[scaDispEst > 1/.Machine$double.eps] <- 1/.Machine$double.eps
  
  # Note on handling very low probabilities: vecLikZeros
  # typically does not have zero elements as it has the 
  # the summand drop-out rate. Also the log cannot be
  # broken up over the sum to dnbinom. In contrast to that,
  # the log is taken in dnbinom for vecLikNonzeros to avoid 
  # zero probabilities. Zero probabilities are handled
  # through substitution of the minimal probability under
  # machine precision. The model is scaled according to the
  # size factors which are used for mean parameter inference.
  # Likelihood of zero counts:
  vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
    dnbinom(
      vecY[vecboolZero], 
      mu=vecMuEst[vecboolZero]*vecSizeFactors[vecboolZero], 
      size=scaDispEst, 
      log=FALSE) +
    vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
      sum(vecLikZeros==0)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLikNonzeros <- (1-vecDropoutRateEst[vecboolNotZeroObserved])*
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMuEst[vecboolNotZeroObserved]*vecSizeFactors[vecboolNotZeroObserved], 
      size=scaDispEst, 
      log=FALSE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( log(vecLikNonzeros[vecLikNonzeros!=0]) +
      sum(vecLikNonzeros==0)*log(.Machine$double.eps) )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Fit zero-inflated negative binomial model to data
#' 
#' Fit zero-inflated negative binomial model to data: One mean per cluster 
#' and either one dispersion parameter across all observations in all cluster
#' for a gene or one dispersion parameter per cluster per gene. Dropout rate
#' and dispersion factor inferred here are used as hyperparamters in the
#' impulse fitting stage of PseudoDE.
#' To parallelise this code, replace lapply by bplapply in dispersion
#' factor and drop-out rate estimation and uncomment the BiocParallel
#' register() command at the beginning of this function.
#' Counts can be imputed under the ZINB model as:
#' matCountsProcImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecSpikeInGenes: (string vector) Names of genes
#'    which correspond to external RNA spike-ins. Currently
#'    not used.
#' @param scaWindowRadius: (integer) 
#'    Interval length for mean paramter estimation 
#' @param boolOneDispPerGene: (bool) [Default TRUE]
#'    Whether one negative binomial dispersion factor is fitted
#'    per gene or per gene for each cluster.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' @param boolSuperVerbose: (bool) Whether to follow EM-algorithm
#'    progress in high detail with local convergence flags. 
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' 
#' @return (list length 6)
#'    \itemize{
#'      \item vecDispersions: (numeric matrix genes x clusters)
#'        Inferred negative binomial dispersions.
#'      \item matDropout: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial drop out rates.
#'      \item matProbNB: (numeric matrix genes x cells)
#'        Inferred probabilities of every observation to be generated
#'        from the negative binomial component of the zero-inflated
#'        negative binomial mixture model.
#'      \item matCountsProcImputed: (numeric matrix genes x cells)
#'        Data predicted with inferred zero-inflated negative binomial
#'        model.
#'      \item matMuCluster: (numeric matrix genes x clusters)
#'        Inferred negative binomial cluster means.
#'      \item boolConvergence: (bool) Convergence of EM algorithm.
#'    }
#' @export

fitZINB <- function(matCountsProc, 
  lsResultsClustering,
  vecSizeFactors,
  vecSpikeInGenes=NULL,
  boolOneDispPerGene=TRUE,
  scaWindowRadius=NULL,
  scaMaxiterEM=100,
  verbose=FALSE,
  boolSuperVerbose=FALSE,
  nProc=1){
  
  # Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-5)
  
  # Store EM convergence
  vecEMLogLik <- array(NA,scaMaxiterEM)
  
  # Set number of processes to be used for parallelisation
  # This function is currently not parallelised to reduce memory usage.
  # Read function description for further information.
  register(MulticoreParam(nProc))
  
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]  
  
  # (I) Initialisation
  # The initial model is estimated under the assumption that all zero-counts
  # are drop-outs.
  matMuCluster <- array(NA,c(scaNumGenes,lsResultsClustering$K))
  matMu <- array(NA,c(scaNumGenes,scaNumCells))
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  
  # E-step:
  # Posterior of dropout: matZ
  if(boolSuperVerbose){print("Initialisation E-step: Estimtate posterior of dropout")}
  matZ <- t(apply(matCountsProc==0, 1, as.numeric))
  
  # (II) EM itertion
  scaIter <- 1
  scaLogLikNew <- 0
  scaLogLikOld <- 0
  while(scaIter == 1 | scaIter == 2 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){
    # M-step:
    # a) Negative binomial mean parameter
    # Use MLE of mean parameter of negative binomial distribution: weighted average.
    # Data are scaled by size factors.
    if(boolSuperVerbose){print("M-step: Estimtate negative binomial mean parameters")}
    
    if(!is.null(scaWindowRadius)){
      # Estimate mean parameter for each cell as ZINB model for cells within pseudotime
      # interval with cell density centred at target cell.
      print("only valid for size factors==1")
      #for(j in seq(1,scaNumCells)){
      do.call(cbind, bplapply(seq(1,scaNumCells), 
        function(j){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          #matMu[,j] <- rowSums(
          return( rowSums(
            matCountsProc[,vecInterval]*
              (1-matZ)[,vecInterval],
            na.rm=TRUE) / rowSums((1-matZ)[,vecInterval])
          )
        })
      )
      matMu[matMu==0 | is.na(matMu)] <- 1/scaNumCells
    } else {
      # Estimate mean parameter by cluster.
      print("only valid for size factors==1")
      for(k in seq(1,max(vecindClusterAssign))){
        matMuCluster[,k] <- rowSums(
          matCountsProc[,vecindClusterAssign==k]*
            (1-matZ)[,vecindClusterAssign==k],
          na.rm=TRUE) / rowSums((1-matZ)[,vecindClusterAssign==k])
      }
      # Add pseudocounts to zeros
      matMuCluster[matMuCluster==0 | is.na(matMuCluster)] <- 1/scaNumCells
      matMu <- matMuCluster[,vecindClusterAssign]
    }
    
    # b) Dropout rate
    # Fit dropout rate with GLM
    if(boolSuperVerbose){print("M-step: Estimtate dropout rate")}
    vecPiFit <- bplapply(seq(1,scaNumCells), function(j) {
      glm( matZ[,j] ~ log(matMu[,j]*vecSizeFactors[j]),
        family=binomial(link=logit),
        control=list(maxit=1000)
      )[c("converged","fitted.values")]
    })
    vecboolConvergedGLMpi <- sapply(vecPiFit, function(x) x[[1]])
    matDropout <- sapply(vecPiFit, function(x) x[[2]])
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(!vecboolConvergedGLMpi), " cases."))
    }
    
    # c) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood
    if(boolSuperVerbose){print("M-step: Estimtate negative binomial dispersion parameters")}
    matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
    matboolZero <- matCountsProc == 0
    scaDispGues <- 1
    if(boolOneDispPerGene){  
      vecDispFit <- bplapply(seq(1,scaNumGenes), function(i){
        tryCatch({ unlist(optim(
          par=log(scaDispGues),
          fn=evalLogLikDispNB,
          vecY=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          method="BFGS",
          control=list(maxit=1000, fnscale=-1)
        )[c("par","convergence")] )
        }, error=function(strErrorMsg){
          print(paste0("ERROR: Fitting negative binomial dispersion parameter: fitZINB().",
            " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
          print(paste0("vecCounts ", paste(matCountsProc[i,],sep=" ")))
          print(paste0("vecMeans ", paste(matMu[i,],sep=" ")))
          print(paste0("scaDispEst ", paste(scaDispGues,sep=" ")))
          print(paste0("vecDropout ", paste(matDropout[i,],sep=" ")))
          print(paste0("vecSizeFactors ", paste(vecSizeFactors,sep=" ")))
          lsErrorCausingGene <- list(matCountsProc[i,], matMu[i,], log(scaDispGues), vecSizeFactors, matDropout[i,])
          names(lsErrorCausingGene) <- c("vecCounts", "veMu", "logscaDispEst", "vecSizeFactors", "vecDropout")
          save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
          stop(strErrorMsg)
        })
      })
    } else {
      #  not coded
      stop("not coded")
    }
    vecboolConvergedGLMdisp <- sapply(vecDispFit, function(x) x["convergence"])
    vecDispersions <- sapply(vecDispFit, function(fit){ exp(fit["par"]) })
    matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    matDispersions[matDispersions<=0] <- min(matDispersions[matDispersions>0])
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(vecboolConvergedGLMdisp), " cases."))
    }
    
    # E-step:
    if(boolSuperVerbose){print("E-step: Estimtate posterior of dropout")}
    matZ <- matDropout/(matDropout + (1-matDropout)*
        dnbinom(0, mu = matMu, size = matDispersions) )
    matZ[matCountsProc > 0] <- 0
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    matboolNotZeroObserved <- matCountsProc >0 & !is.na(matCountsProc)
    matboolZero <- matCountsProc==0
    scaLogLikNew <- sum(sapply( seq(1,scaNumGenes), function(i){
      evalLogLikZINB_PseudoDE_comp(vecY=matCountsProc[i,],
        vecMu=matMu[i,],
        vecDispEst=matDispersions[i,], 
        vecDropoutRateEst=matDropout[i,],
        vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
        vecboolZero=matboolZero[i,])
    }))
    
    # EM-iteration complete
    if(verbose){print(paste0("Completed Iteration ", scaIter, " with data log likelihood of ", scaLogLikNew))}
    vecEMLogLik[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  # Evaluate convergence
  if(all(as.logical(vecboolConvergedGLMdisp)) & all(as.logical(vecboolConvergedGLMpi)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergence <- TRUE
  } else {
    boolConvergence <- FALSE
  }
  
  # Compute mixture probabilities and imputed counts
  matProbNB <- 1 - matZ
  
  # Name rows and columns of output
  rownames(matDropout) <- rownames(matCountsProc)
  colnames(matDropout) <- colnames(matCountsProc)
  rownames(matMuCluster) <- rownames(matCountsProc)
  rownames(matMu) <- rownames(matCountsProc)
  colnames(matMu) <- colnames(matCountsProc)
  names(vecDispersions) <- rownames(matCountsProc)
  rownames(matDispersions) <- rownames(matCountsProc)
  colnames(matDispersions) <- colnames(matCountsProc)
  rownames(matProbNB) <- rownames(matCountsProc)
  colnames(matProbNB) <- colnames(matCountsProc)
  
  # Check dispersions
  if(any(is.na(matDispersions) | !is.finite(matDispersions))){
    matDispersions[is.na(matDispersions) | !is.finite(matDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions. Set to 1.")
  }
  
  return(list( vecDispersions=vecDispersions,
    matDropout=matDropout, 
    matProbNB=matProbNB,
    matMuCluster=matMuCluster,
    matMu=matMu,
    boolConvergence=boolConvergence,
    vecEMLogLik=vecEMLogLik ))
}