#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Note: Most problems are numerical and not conceptual problems.

# Previously encountered problems:
# Genes with many zeros and few very high counts: If dispersion is initialised
# too high (=1), the increase in LL cannot always be picked up by dnbinom which
# disables convergence.
# Still didnt see convergence, this might reoccur during fitting. Removed -Inf
# masking from likelihood function for NB part. Seeing convergence of means now.
# Still problems with drop-out rate -> Use more genes to get more reliable 
# traces, still problems. Take out masking of dropout contribution of nonzeros
# in likelihood and force logistic to be decreasing in means - good cnvergence
# until 3rd/4th iter.
# Try replace mean optimise with optim and put masking in likelihood again: Note
# that optim becomes unstable if Inf is returned by objective, while optimise
# seems to be forced out of the "bad" parameter region in which dnbinom yields
# LL=-Inf, optim needs the low value masking to not break down. Optim can be
# initialised outside of the bad parameter region to avoid it. Now only 
# problems in mean estimation (3rd iteration). Remove bounds from BFGS
# optimisation. Cleared problem by setting mean initialisation to log(max+1)
# instead of log(mean+1), was it stuck in bad parameter region again? maybe.
# Initialise mean and drop-out rate to prior value, converging fully!
# Compared BFGS against Nelder-Mead in optim for mean estimation on 200x400
# points: Completed 4 iter with BFGS in 36.03 min on 2 cores 
# (LL=24536057.4015036). Completed 4 iter with Nelder-Mead in 35.22 min on 
# 2 cores (LL=-24536058.6125272). Exact same data set. Both yield 172 DE genes.
# Problem with impulse mode: LL returned from optim does not match LL
# computed from parameters returned by optim. Effect not so heavy in first 
# iteration but very heavy after that. Not due to matLinModel, is it dispersion?
# Dispersion is initialised low, if initialised higher, effects are stronger
# in first iteration already. Is this inaccuracy in reporting parameters? Wouldnt
# be so bad in first iteration with low dispersion param!?
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

fitZINBMu <- function( matCountsProc,
  vecPseudotime=NULL,
  vecindClusterAssign=NULL,
  matMu,
  matDispersions,
  vecSizeFactors,
  matDropout,
  matConstPredictorsPi=NULL,
  matLinModelPi=NULL,
  matZ,
  matImpulseParam=NULL,
  scaNumGenes,
  scaNumCells,
  scaWindowRadius,
  boolDynamicPi,
  boolBFGSEstOfMu,
  strMuModel ){

  if(strMuModel=="windows"){
    # Estimate mean parameter for each cell as ZINB model for cells within pseudotime
    # interval with cell density centred at target cell.
    # Note that this corresponds to maximising smoothed log likelihood but instead
    # of using the implemented cost function evalLogLikSmoothZINB_LinPulse for an entire
    # gene, the optimisation problem is broken up into 1D problems for each mean.
    matMu <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      # Note: Mean parameter estimates of gene i depend on each other:
      # Either estimate all parameters for gene i together 
      # (quasi-Newton estimation with BFGS: boolBFGSEstOfMu=TRUE)
      # or use the latest updates of the remaining parameters during 
      # one-by-one estimation (coordinate ascent, boolBFGSEstOfMu=TRUE).
      if(boolBFGSEstOfMu){
        vecMu <- fitMuVecZINB_LinPulse(
          vecCounts=matCountsProc[i,],
          vecMu=matMu[i,],
          vecDisp=matDispersions[i,],
          vecNormConst=vecSizeFactors,
          vecDropoutRateEst=matDropout[i,],
          vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,]),
          matLinModelPi=matLinModelPi,
          scaWindowRadius=scaWindowRadius,
          boolDynamicPi=boolDynamicPi )
      } else {
        vecMu <- matMu[i,]
        for(j in seq(1,scaNumCells)){
          scaindIntervalStart <- max(1,j-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
          vecMu[j] <- fitMuZINB_LinPulse(vecCounts=matCountsProc[i,vecInterval],
            vecMu=vecMu[vecInterval],
            vecDisp=matDispersions[i,vecInterval],
            vecNormConst=vecSizeFactors[vecInterval],
            vecDropoutRateEst=matDropout[i,vecInterval],
            vecProbNB=1-matZ[i,vecInterval],
            vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,]),
            matLinModelPi=matLinModelPi[vecInterval,],
            scaTarget=match(j,vecInterval),
            scaWindowRadius=scaWindowRadius,
            strMuModel=strMuModel,
            boolDynamicPi=boolDynamicPi )
          # vecDropout is re-estimated in mean estimation based on the new mean.
        }
      }
      return(vecMu)
    }))
  } else if(strMuModel=="clusters"){
    # Estimate mean parameter by cluster. No smoothing is used.
    matMuCluster <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
      vecMu <- sapply(seq(1,max(vecindClusterAssign)), function(k){
        vecInterval <- vecindClusterAssign==k
        scaMu <- fitMuZINB_LinPulse(
          vecCounts=matCountsProc[i,vecInterval],
          vecMu=matMu[i,vecInterval],
          vecDisp=matDispersions[i,vecInterval],
          vecNormConst=vecSizeFactors[vecInterval],
          vecDropoutRateEst=matDropout[i,vecInterval],
          vecProbNB=1-matZ[i,vecInterval],
          vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,]),
          matLinModelPi=matLinModelPi[vecInterval,],
          strMuModel=strMuModel,
          boolDynamicPi=boolDynamicPi )
        return(scaMu)
      })
      return(vecMu)
    }))
    matMu <- matMuCluster[,vecindClusterAssign]
  } else if(strMuModel=="impulse"){
    lsImpulseFits <- bplapply(seq(1,scaNumGenes), function(i){
      lsImpulseFit <- fitMuImpulseZINB_LinPulse(
        vecCounts=matCountsProc[i,],
        vecDisp=matDispersions[i,],
        vecNormConst=vecSizeFactors,
        vecDropoutRateEst=matDropout[i,],
        vecProbNB=1-matZ[i,],
        vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,]),
        matLinModelPi=matLinModelPi,
        vecPseudotime=vecPseudotime,
        vecImpulseParam=matImpulseParam[i,],
        boolDynamicPi=boolDynamicPi )
      return(lsImpulseFit)
    })
    matMu <- do.call(rbind, lapply(lsImpulseFits, function(x) x$vecImpulseValue))
    matImpulseParam <- do.call(rbind, lapply(lsImpulseFits, function(x) x$vecBestFitParam))
  } else if(strMuModel=="constant"){
    vecMu <- unlist(bplapply( seq(1,scaNumGenes), function(i){
      fitMuConstZINB( vecCounts=matCountsProc[i,],
        vecDisp=matDispersions[i,],
        vecNormConst=vecSizeFactors,
        vecDropoutRateEst=matDropout[i,],
        matLinModelPi=matLinModelPi,
        scaWindowRadius=scaWindowRadius )
    }))
    matMu <- matrix(vecMu, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  }
  
  return(list( matMu=matMu,
    matImpulseParam=matImpulseParam ))
}

#' Fit zero-inflated negative binomial model to data
#' 
#' Fit zero-inflated negative binomial model to data by composite 
#' coordinate ascent, every step is a partial MLE step. The scheme
#' iterates over cell-wise (drop-out models) and gene-wise (negative 
#' binomial models) parameters. Convergence is tracked with the the
#' loglikelihood of the entire data matrix. Every step is a strict MLE
#' of the target parameters given the remaining parameter estimates,
#' therefore, convergence to a local optimum is guaranteed. Parallelisation
#' can be performed as independences of parameters are exploited where
#' possible. Convergence can be followed with verbose=TRUE (at each 
#' iteration) or at each step (boolSuperVerbose=TRUE). Variables for the
#' logistic drop-out model are a constant and the estimated mean parameter
#' and other constant gene-specific variables (such as GC-conten) in 
#' matConstPredictorsPi. Three modes are available for modelling the mean
#' parameter: By cluster (this is fast as neighbourhoods don't have to
#' be evaluated), sliding windows (recommended) and as an impulse model 
#' (under construction).
#' Counts can be imputed under the ZINB model as:
#' matCountsProcImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
#' 
#' @seealso Called by \code{runLineagePulse}.
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
#'    Smoothing interval radius.
#' @param boolOneDispPerGene: (bool) [Default TRUE]
#'    Whether one negative binomial dispersion factor is fitted
#'    per gene or per gene for each cluster.
#' @param vecPseudotime: (numerical vector number of cells)
#'    Pseudotime coordinates of cells.
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
  strMuModel="windows",
  vecPseudotime=NULL,
  scaMaxiterEM=100,
  verbose=FALSE,
  boolSuperVerbose=FALSE,
  nProc=1){
  
  # Compute degrees of freedom of model for each gene
  # Drop-out model is ignored, would be counted at each gene.
  # Mean model:
  if(strMuModel=="windows"){
    # One mean parameter per cell
    scaKbyGene <- dim(matCountsProc)[2]
  } else if(strMuModel=="clusters"){
    # One mean parameter per cluster
    scaKbyGene <- lsResultsClustering$K
  } else if(strMuModel=="impulse"){
    # Six impulse model parameter to model means
    scaKbyGene <- 6
  } else if(strMuModel=="constant"){
    # One constant mean
    scaKbyGene <- 1
  } else {
    stop(paste0("ERROR in fitZINB(): strMuModel not recognised: ", strMuModel))
  }
  # Dispersion model
  if(boolOneDispPerGene){
    # One dispersion factor per gene.
    scaKbyGene <- scaKbyGene + 1
  } else {
    stop(paste0("ERROR in fitZINB(): boolOneDispPerGene set to FALSE."))
  }
  
  # Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-6)
  
  # Store EM convergence
  vecEMLogLik <- array(NA,scaMaxiterEM)
  scaPredictors <- 2
  boolDynamicPi <- TRUE
  boolBFGSEstOfMu <- TRUE
  matConstPredictorsPi <- NULL
  
  # Set smoothing mode
  boolSmoothed <- FALSE
  if(!is.null(scaWindowRadius)){
    if(scaWindowRadius>0){
      boolSmoothed <- TRUE
    } else {
      scaWindowRadius <- NULL
    }
  }
  
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
  # Pre-compute constant matrices
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
  matboolZero <- matCountsProc == 0
  
  # Initialise parameters:
  # Mean parameters (mu): Gene-wise mean of non-zero observations
  # Dispersions: Low dispersion factor yielding high variance which makes
  # cost function screening easy in the first iteration.
  # Drop-out rate: Set to 0.5, ie a constant model with all parameters zero.
  # As the parameter corresponding to log(mu) is forced to be negative later,
  # this parameter is initialised as very close to zero but negative.  
  vecMu <- apply(matCountsProc, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
  vecMu[vecMu < .Machine$double.eps] <- .Machine$double.eps
  matMu <- matrix(vecMu, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  matDispersions <- matrix(0.001, nrow=scaNumGenes, ncol=scaNumCells)
  matDropout <- matrix(0.5, nrow=scaNumGenes, ncol=scaNumCells)
  matLinModelPi <- cbind(rep(0, scaNumCells), rep(-10^(-10), scaNumCells),
    matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
  if(strMuModel=="impulse"){ 
    matImpulseParam <- matrix(1, nrow=scaNumGenes, ncol=6)
    matImpulseParam[,1:3] <- log(matrix(vecMu, nrow=scaNumGenes, ncol=3, byrow=FALSE))
  } else { matImpulseParam <- matrix(NA, nrow=scaNumGenes, ncol=6) }
  if(boolSuperVerbose){
    scaLogLikInit <- evalLogLikMatrix(matCounts=matCountsProc,
      matMu=matMu,
      vecSizeFactors=vecSizeFactors,
      matDispersions=matDispersions, 
      matDropout=matDropout, 
      matboolNotZeroObserved=matboolNotZeroObserved, 
      matboolZero=matboolZero,
      scaWindowRadius=scaWindowRadius,
      boolSmoothed=boolSmoothed )
    print(paste0("Completed initialisation with ",
      "log likelihood of         ", scaLogLikInit))
  }
  
  # (II) EM itertion
  scaIter <- 1
  scaLogLikNew <- 0
  scaLogLikOld <- 0
  
  while(scaIter == 1 | scaIter == 2 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){
    ##### 1. Gene-wise parameter estimation: 
    # a) Negative binomial mean parameter
    # Only compute posterior if using closed form estimator for mean:
    # Posterior is not necessary in all other cases, expect for parameter
    # initialisation of impulse model.
    if((all(vecSizeFactors==1) & !boolDynamicPi) | strMuModel=="impulse"){    
      matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
        (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
      }))
      matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
        vecZ <- sapply(seq(1,scaNumCells), function(j){
          if(matboolNotZeroObserved[i,j]){
            scaZ <- 0
          } else {
            if(boolSmoothed){
              scaindIntervalStart <- max(1,j-scaWindowRadius)
              scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
              vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
            } else {
              vecInterval <- j
            }
            scaZ <- sum(matDropout[i,j]/(matDropout[i,j] + 
                (1-matDropout[i,j])*matNBZero[i,vecInterval])) *
              1/length(vecInterval)
          }
          return(scaZ)
        })
        return(vecZ)
      }))
    }
    lsMuFits <- fitZINBMu( matCountsProc=matCountsProc,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign,
      matMu=matMu,
      matDispersions=matDispersions,
      vecSizeFactors=vecSizeFactors,
      matDropout=matDropout,
      matConstPredictorsPi=matConstPredictorsPi,
      matLinModelPi=matLinModelPi,
      matZ=matZ,
      matImpulseParam=matImpulseParam,
      scaNumGenes=scaNumGenes,
      scaNumCells=scaNumGenes,
      scaWindowRadius=scaWindowRadius,
      boolDynamicPi=boolDynamicPi,
      boolBFGSEstOfMu=boolBFGSEstOfMu,
      strMuModel=strMuModel )
    matMu <- lsMuFits$matMu
    matImpulseParam <- lsMuFits$matImpulseParam
    # These udates are done during mean estimation too
    # Reestimate drop-out rates based on new means
    matDropout <- do.call(cbind, bplapply(seq(1, scaNumCells), function(cell){
      vecLinModelOut <- cbind(1,log(matMu[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
      vecDropout <- 1/(1+exp(-vecLinModelOut))
      return(vecDropout)
    }))
    if(boolSuperVerbose){
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMu,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersions, 
        matDropout=matDropout,
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius,
        boolSmoothed=boolSmoothed )
      print(paste0("# ",scaIter, ".1) Mean estimation complete: ",
        "loglikelihood of       ", scaLogLikTemp))
    }
    
    # b) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood.
    if(boolOneDispPerGene){  
      vecDispFit <- bplapply(seq(1,scaNumGenes), function(i){
        fitDispZINB_LinPulse(scaDispGuess=matDispersions[i,1],
          vecCounts=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius,
          boolSmoothed=boolSmoothed )
      })
    } else {
      #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
      stop("Non constant-dispersion factors are not yet implemented. david.seb.fischer@gmail.com")
    }
    vecboolConvergedGLMdisp <- sapply(vecDispFit, function(fit) fit["convergence"])
    vecDispersions <- sapply(vecDispFit, function(fit){ fit["par"] })
    matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    if(boolSuperVerbose){
      if(any(vecboolConvergedGLMdisp != 0)){
        print(paste0("Dispersion estimation did not converge in ", 
          sum(vecboolConvergedGLMdisp), " cases."))
      }
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMu,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersions, 
        matDropout=matDropout, 
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius,
        boolSmoothed=boolSmoothed )
      print(paste0("# ",scaIter, ".2) Dispersion estimation complete: ",
        "loglikelihood of ", scaLogLikTemp))
    }  
    
    #####  2. Cell-wise parameter estimation
    # Dropout rate
    # Drop-out estimation is independent between cells and can be parallelised.
    matLinModelPi <- do.call(rbind, 
      bplapply(seq(1, scaNumCells), function(cell){
        if(boolSmoothed){
          scaindIntervalStart <- max(1,cell-scaWindowRadius)
          scaindIntervalEnd <- min(scaNumCells,cell+scaWindowRadius)
          vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        } else {
          vecInterval <- cell
        }
        
        vecLinModelPi <- fitPiZINB_LinPulse(
          vecLinModelPi=matLinModelPi[cell,],
          matPredictorsPi=cbind(1,log(matMu[,cell]),matConstPredictorsPi),
          vecCounts=matCountsProc[,cell],
          matMu=matMu[,vecInterval],
          matDisp=matDispersions[,vecInterval],
          scaNormConst=vecSizeFactors[cell],
          boolSmoothed=boolSmoothed )
        return(vecLinModelPi)
      }))
    matDropout <- do.call(cbind, 
      lapply(seq(1, scaNumCells), function(cell){
        vecLinModelOut <- cbind(1,log(matMu[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
        vecDropout <- 1/(1+exp(-vecLinModelOut))
        return(vecDropout)
      }))
    if(boolSuperVerbose){
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMu,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersions, 
        matDropout=matDropout, 
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius,
        boolSmoothed=boolSmoothed )
      print(paste0("# ",scaIter,".3) Drop-out estimation complete: ",
        "loglikelihood of   ", scaLogLikTemp))
    }  
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    scaLogLikNew <- evalLogLikMatrix(matCounts=matCountsProc,
      matMu=matMu,
      vecSizeFactors=vecSizeFactors,
      matDispersions=matDispersions, 
      matDropout=matDropout, 
      matboolNotZeroObserved=matboolNotZeroObserved, 
      matboolZero=matboolZero,
      scaWindowRadius=scaWindowRadius,
      boolSmoothed=boolSmoothed )
    
    # EM-iteration complete
    if(!boolSuperVerbose){
      if(verbose){print(paste0("# ",scaIter, ".) complete with ",
        "log likelihood of ", scaLogLikNew))}
    }
    vecEMLogLik[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  print("Exited EM-iteration.")
  
  # Evaluate convergenc
  if(all(as.logical(vecboolConvergedGLMdisp)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergence <- TRUE
  } else {
    boolConvergence <- FALSE
  }
  
  # Compute mixture probabilities and imputed counts  
  matNBZero <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    (matDispersions[i,]/(matDispersions[i,]+matMu[i,]))^matDispersions[i,]
  }))
  matZ <- do.call(rbind, bplapply(seq(1,scaNumGenes), function(i){
    vecZ <- sapply(seq(1,scaNumCells), function(j){
      if(matboolNotZeroObserved[i,j]){
        scaZ <- 0
      } else {
        scaindIntervalStart <- max(1,j-scaWindowRadius)
        scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
        vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
        scaZ <- sum(matDropout[i,j]/(matDropout[i,j] + 
            (1-matDropout[i,j])*matNBZero[i,vecInterval])) *
          1/length(vecInterval)
      }
      return(scaZ)
    })
    return(vecZ)
  }))
  matProbNB <- 1 - matZ
  
  # Name rows and columns of output
  rownames(matMu) <- rownames(matCountsProc)
  colnames(matMu) <- colnames(matCountsProc)
  rownames(matDispersions) <- rownames(matCountsProc)
  colnames(matDispersions) <- colnames(matCountsProc)
  rownames(matDropout) <- rownames(matCountsProc)
  colnames(matDropout) <- colnames(matCountsProc)
  rownames(matLinModelPi) <- colnames(matCountsProc)
  rownames(matProbNB) <- rownames(matCountsProc)
  colnames(matProbNB) <- colnames(matCountsProc)
  
  return(list( matMu=matMu,
    matDispersions=matDispersions,
    matDropout=matDropout,
    matDropoutLinModel=matLinModelPi,
    matProbNB=matProbNB,
    matImpulseParam=matImpulseParam,
    boolConvergence=boolConvergence,
    vecEMLogLik=vecEMLogLik,
    scaKbyGene=scaKbyGene ))
}