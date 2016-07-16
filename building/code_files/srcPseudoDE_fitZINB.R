#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of mean parameter estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param matDropoutLinMod: (matrix number of cells x 2) Logistic linear
#'    model parameters of the dropout rate as a function of the mean.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecDisp,
  vecSizeFactors,
  vecDropoutRateEst,
  matDropoutLinMod=NULL,
  vecboolNotZeroObserved,
  vecboolZero){ 
  
  # Log linker function to fit positive means
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  if(TRUE){
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=scaMu*vecSizeFactors,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    # Correct drop-out rate estimate for mean estimate
    # This is not used in the current framework
    vecDropoutRateEst <- 1/(1+exp(-matDropoutLinMod[,1]-log(scaMu)*matDropoutLinMod[,2]))
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=scaMu*vecSizeFactors,
      vecDispEst=vecDisp, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' Note that this fitting routine is meant for fitting a single mean
#' to a set of observations and can be used in wrappers which compute
#' every mean separately (assuming independence).
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaDisp: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitMuZINB_LinPulse <- function(vecCounts,
  vecDisp,
  vecNormConst,
  vecDropoutRateEst,
  vecProbNB){ 
  
  if(all(vecNormConst==1)){
    # Closed form maximum likelihood estimator
    scaMu <- sum(vecCounts*vecProbNB, na.rm=TRUE)/sum(vecProbNB, na.rm=TRUE)
  } else {
    # Numerical maximum likelihood estimator
    scaMu <- tryCatch({
      exp(unlist(optimise(
        evalLogLikMuZINB_LinPulse_comp,
        vecCounts=vecCounts,
        vecDisp=vecDisp,
        vecDropoutRateEst=vecDropoutRateEst,
        vecNormConst=vecNormConst,
        vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
        vecboolObserved=!is.na(vecCounts),
        lower = log(.Machine$double.eps),
        upper = log(max(vecNormConst*vecCounts, na.rm=TRUE)+1),
        maximum = TRUE)["maximum"]))
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuZINB_LinPulse().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
      print(paste0("vecDisp ", paste(vecDisp,collapse=" ")))
      print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
      lsErrorCausingGene <- list(vecCounts, vecDisp, vecDropoutRateEst, vecNormConst)
      names(lsErrorCausingGene) <- c("vecCounts", "vecDisp", "vecDropoutRateEst","vecNormConst")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      stop(strErrorMsg)
    })
  }
  
  # Catch boundary of likelihood domain on mu space
  if(is.na(scaMu) | scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
}

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
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
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

evalLogLikDispZINB_LinPulse <- function(scaTheta,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL){ 
  
  # Log linker function to fit positive dispersions
  scaDispEst <- exp(scaTheta)
  
  # Prevent dispersion estimate from shrinking to zero
  # to avoid numerical errors:
  # Could also write as if statement, this accomodates for vectors later.
  scaDispEst[scaDispEst < .Machine$double.eps] <- .Machine$double.eps
  # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  scaDispEst[scaDispEst > 1/.Machine$double.eps] <- 1/.Machine$double.eps
  
  vecDispersions <- rep(scaDispEst, length(vecCounts))
  
  if(is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikZINB_LinPulse_comp( vecCounts=vecCounts,
      vecMu=vecMuEst*vecSizeFactors,
      vecDispEst=vecDispersions, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    scaLogLik <- evalLogLikSmoothZINB_LinPulse( vecCounts=vecCounts,
      vecMu=vecMuEst*vecSizeFactors,
      vecDispEst=vecDispersions, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius )
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuEst: (vector number of cells) Negative binomial
#'    mean parameter estimate of clusters to which cells
#'    belong.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell. 
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitDispZINB_LinPulse <- function( scaDispGuess,
  vecCounts,
  vecMuEst,
  vecSizeFactors,
  vecDropoutRateEst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL){ 
  
  fitDisp <- tryCatch({ unlist(optim(
    par=log(scaDispGuess),
    fn=evalLogLikDispZINB_LinPulse_comp,
    vecCounts=vecCounts,
    vecMuEst=vecMuEst,
    vecSizeFactors=vecSizeFactors,
    vecDropoutRateEst=vecDropoutRateEst,
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero,
    scaWindowRadius=scaWindowRadius,
    method="BFGS",
    control=list(maxit=1000, fnscale=-1)
  )[c("par","convergence")] )
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting zero-inflated negative binomial dispersion parameter: fitDispZINB_LinPulse().",
      " Wrote report into LineagePulse_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("vecMuEst ", paste(vecMuEst,collapse=" ")))
    print(paste0("vecDropout ", paste(vecDropout,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    lsErrorCausingGene <- list(vecCounts, vecMuEst, log(scaDispGuess), vecSizeFactors, vecDropout)
    names(lsErrorCausingGene) <- c("vecCounts", "vecMuEst", "logscaDispEst", "vecSizeFactors", "vecDropout")
    save(lsErrorCausingGene,file=file.path(getwd(),"LineagePulse_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  fitDisp["par"] <- exp(fitDisp["par"])
  
  # Catch boundary of likelihood domain on mu space
  if(fitDisp["par"] < .Machine$double.eps){fitDisp["par"] <- .Machine$double.eps}
    # Prevent dispersion estimate from growing to infinity
  # to avoid numerical errors:
  if(fitDisp["par"] > 1/.Machine$double.eps){fitDisp["par"] <- 1/.Machine$double.eps}
  
  return(fitDisp)
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
#'    Smoothing interval length.
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
  scaPrecEM <- 1-10^(-6)
  
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
  
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  # E-step:
  # Posterior of dropout: matZ
  if(boolSuperVerbose){print("Initialisation E-step: Estimtate posterior of dropout")}
  matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    as.numeric(matCountsProc[i,]==0)
  }))
  matDropout <- matZ
  
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
      # Note that this corresponds to maximising smoothed log likelihood but instead
      # of using the implemented cost function evalLogLikSmoothZINB_LinPulse for an entire
      # gene, the optimisation problem is broken up into 1D problems for each mean.
      matMu <- do.call(cbind, bplapply(seq(1,scaNumCells), function(j){
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        vecMu <- sapply(seq(1,scaNumGenes), function(i){
          scaMu <- fitMuZINB_LinPulse(vecCounts=matCountsProc[i,vecInterval],
            vecDisp=matDispersions[i,vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            vecDropoutRateEst=matDropout[i,vecInterval],
            vecProbNB=1-matZ[i,vecInterval] )
          return(scaMu)
        })
        return(vecMu)
      }))
    } else {
      # Estimate mean parameter by cluster. No smoothing is used.
      matMuCluster <- do.call(cbind, lapply(seq(1,max(vecindClusterAssign)), function(k){
        vecMu <- unlist(bplapply(seq(1,scaNumGenes), function(i){
          scaMu <- fitMuZINB_LinPulse(vecCounts=matCountsProc[i,vecindClusterAssign==k],
            vecDisp=matDispersions[i,vecindClusterAssign==k],
            vecNormConst=vecSizeFactors[vecindClusterAssign==k],
            vecDropoutRateEst=matDropout[i,vecindClusterAssign==k],
            vecProbNB=1-matZ[i,vecindClusterAssign==k])
          return(scaMu)
        }))
        return(vecMu)
      }))
      matMu <- matMuCluster[,vecindClusterAssign]
    }
    
    # b) Dropout rate
    # Fit dropout rate with GLM
    if(boolSuperVerbose){print("Constrained E-step: Estimtate dropout rate")}
    vecPiFit <- bplapply(seq(1,scaNumCells), function(j) {
      glm( matZ[,j] ~ log(matMu[,j]*vecSizeFactors[j]),
        family=binomial(link=logit),
        control=list(maxit=1000)
      )[c("converged","fitted.values","coefficients")]
    })
    vecboolConvergedGLMpi <- sapply(vecPiFit, function(x) x[[1]])
    matDropout <- sapply(vecPiFit, function(x) x[[2]])
    matDropoutLinModel <- do.call(rbind, lapply(vecPiFit, function(x) x[[3]]))
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(!vecboolConvergedGLMpi), " cases."))
    }
    
    # c) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood.
    if(boolSuperVerbose){print("M-step: Estimtate negative binomial dispersion parameters")}
    matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
    matboolZero <- matCountsProc == 0
    scaDispGuess <- 1
    if(boolOneDispPerGene){  
      vecDispFit <- bplapply(seq(1,scaNumGenes), function(i){
        fitDispZINB_LinPulse(scaDispGuess=scaDispGuess,
          vecCounts=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius )
      })
    } else {
      #  not coded
      stop("not coded")
    }
    vecboolConvergedGLMdisp <- sapply(vecDispFit, function(fit) fit["convergence"])
    vecDispersions <- sapply(vecDispFit, function(fit){ fit["par"] })
    matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(vecboolConvergedGLMdisp), " cases."))
    }
    
    # E-step:
    if(boolSuperVerbose){print("E-step: Estimtate posterior of dropout")}
    # Use closed form expression of negative binomial density at x=0:
    matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
    }))
    matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      matDropout[i,]/(matDropout[i,] + (1-matDropout[i,])*matNBZero[i,])
    }))
    #matZ <- matDropout/(matDropout + (1-matDropout)*
    #    dnbinom(0, mu = matMu, size = matDispersions) )
    matZ[matCountsProc > 0] <- 0
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    matboolNotZeroObserved <- matCountsProc>0 & !is.na(matCountsProc)
    matboolZero <- matCountsProc==0
    if(is.null(scaWindowRadius)){
      scaLogLikNew <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,])
        })
      ))
    } else {
      scaLogLikNew <- sum(unlist(
        bplapply( seq(1,scaNumGenes), function(i){
          evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
            vecMu=matMu[i,],
            vecDispEst=matDispersions[i,], 
            vecDropoutRateEst=matDropout[i,],
            vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
            vecboolZero=matboolZero[i,],
            scaWindowRadius=scaWindowRadius)
        })
      ))
    }
    
    # EM-iteration complete
    if(verbose){print(paste0("Completed iteration ", scaIter, " with log likelihood of ", scaLogLikNew))}
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
  rownames(matMu) <- rownames(matCountsProc)
  colnames(matMu) <- colnames(matCountsProc)
  names(vecDispersions) <- rownames(matCountsProc)
  rownames(matDispersions) <- rownames(matCountsProc)
  colnames(matDispersions) <- colnames(matCountsProc)
  rownames(matProbNB) <- rownames(matCountsProc)
  colnames(matProbNB) <- colnames(matCountsProc)
  if(is.null(scaWindowRadius)){
    rownames(matMuCluster) <- rownames(matCountsProc)
  } else {
    matMuCluster <- NA
  }
  
  # Check dispersions
  if(any(is.na(matDispersions) | !is.finite(matDispersions))){
    matDispersions[is.na(matDispersions) | !is.finite(matDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions. Set to 1.")
  }
  
  return(list( vecDispersions=vecDispersions,
    matDropout=matDropout,
    matDropoutLinModel=matDropoutLinModel,
    matProbNB=matProbNB,
    matMuCluster=matMuCluster,
    matMu=matMu,
    boolConvergence=boolConvergence,
    vecEMLogLik=vecEMLogLik ))
}