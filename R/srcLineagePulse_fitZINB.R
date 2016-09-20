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
  matMu=NULL,
  matDispersions,
  vecSizeFactors,
  matDropout,
  matZ=NULL,
  matConstPredictorsPi=NULL,
  matLinModelPi=NULL,
  matImpulseParam=NULL,
  scaWindowRadius,
  boolDynamicPi,
  boolBFGSEstOfMu=NULL,
  strMuModel ){

  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
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
    scaKbyGeneH1 <- dim(matCountsProc)[2]
  } else if(strMuModel=="clusters"){
    # One mean parameter per cluster
    scaKbyGeneH1 <- lsResultsClustering$K
  } else if(strMuModel=="impulse"){
    # Six impulse model parameter to model means
    scaKbyGeneH1 <- 6
  } else if(strMuModel=="constant"){
    # One constant mean
    scaKbyGeneH1 <- 1
  } else {
    stop(paste0("ERROR in fitZINB(): strMuModel not recognised: ", strMuModel))
  }
  # Dispersion model
  if(boolOneDispPerGene){
    # One dispersion factor per gene.
    scaKbyGeneH1 <- scaKbyGeneH1 + 1
  } else {
    stop(paste0("ERROR in fitZINB(): boolOneDispPerGene set to FALSE."))
  }
  
  # Parameters:
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations:
  scaPrecEM <- 1-10^(-4)
  
  # Store EM convergence
  vecEMLogLikH1 <- array(NA,scaMaxiterEM)
  vecEMLogLikH0 <- array(NA,scaMaxiterEM)
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
  
  # Pre-compute constant matrices
  matSizeFactors <- matrix(vecSizeFactors, nrow=scaNumGenes, ncol=scaNumCells, byrow=TRUE)
  matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
  matboolZero <- matCountsProc == 0
  
  print("### a) Fit alternative negative binomial model and noise model.")
  # (I) Initialisation
  # Initialise parameters:
  # Mean parameters (mu): Gene-wise mean of non-zero observations
  # Dispersions: Low dispersion factor yielding high variance which makes
  # cost function screening easy in the first iteration.
  # Drop-out rate: Set to 0.5, ie a constant model with all parameters zero.
  # As the parameter corresponding to log(mu) is forced to be negative later,
  # this parameter is initialised as very close to zero but negative.  
  vecMuH1 <- apply(matCountsProc, 1, function(gene) mean(gene[gene>0], na.rm=TRUE))
  vecMuH1[vecMuH1 < .Machine$double.eps] <- .Machine$double.eps
  matMuH1 <- matrix(vecMuH1, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
  matDispersionsH1 <- matrix(0.001, nrow=scaNumGenes, ncol=scaNumCells)
  matDropoutH1 <- matrix(0.5, nrow=scaNumGenes, ncol=scaNumCells)
  matLinModelPi <- cbind(rep(0, scaNumCells), rep(-10^(-10), scaNumCells),
    matrix(0, nrow=scaNumCells, ncol=scaPredictors-2))
  if(strMuModel=="impulse"){ 
    matImpulseParam <- matrix(1, nrow=scaNumGenes, ncol=6)
    matImpulseParam[,1:3] <- log(matrix(vecMuH1, nrow=scaNumGenes, ncol=3, byrow=FALSE))
  } else { matImpulseParam <- matrix(NA, nrow=scaNumGenes, ncol=6) }
  scaLogLikInit <- evalLogLikMatrix(matCounts=matCountsProc,
    matMu=matMuH1,
    vecSizeFactors=vecSizeFactors,
    matDispersions=matDispersionsH1, 
    matDropout=matDropoutH1, 
    matboolNotZeroObserved=matboolNotZeroObserved, 
    matboolZero=matboolZero,
    scaWindowRadius=scaWindowRadius )
  if(verbose){
    print(paste0("Completed initialisation with ",
      "log likelihood of         ", scaLogLikInit))
  }
  
  # (II) Estimation itertion: 
  # Estimate H1 negative binomial model and logistic noise model
  scaIter <- 1
  scaLogLikNew <- scaLogLikInit
  scaLogLikOld <- NA
  
  while(scaIter == 1 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){    
    #####  1. Cell-wise parameter estimation
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
          matPredictorsPi=cbind(1,log(matMuH1[,cell]),matConstPredictorsPi),
          vecCounts=matCountsProc[,cell],
          matMu=matMuH1[,vecInterval],
          matDisp=matDispersionsH1[,vecInterval],
          scaNormConst=vecSizeFactors[cell],
          boolSmoothed=boolSmoothed )
        return(vecLinModelPi)
      }))
    matDropoutH1 <- do.call(cbind, 
      lapply(seq(1, scaNumCells), function(cell){
        vecLinModelOut <- cbind(1,log(matMuH1[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
        vecDropout <- 1/(1+exp(-vecLinModelOut))
        return(vecDropout)
      }))
    if(boolSuperVerbose){
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMuH1,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersionsH1, 
        matDropout=matDropoutH1, 
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius )
      print(paste0("# ",scaIter,".1) Drop-out estimation complete: ",
        "loglikelihood of   ", scaLogLikTemp))
    }  
    
    ##### 2. Gene-wise parameter estimation: 
    # a) Negative binomial mean parameter
    # Only compute posterior if using closed form estimator for mean:
    # Posterior is not necessary in all other cases, expect for parameter
    # initialisation of impulse model.
    matZH1 <- calcProbNB( matMu=matMuH1,
      matDispersions=matDispersionsH1,
      matDropout=matDropoutH1,
      matboolZero=matboolZero,
      matboolNotZeroObserved=matboolNotZeroObserved,
      scaWindowRadius=scaWindowRadius )
    
    lsMuFitsH1 <- fitZINBMu( matCountsProc=matCountsProc,
      vecPseudotime=vecPseudotime,
      vecindClusterAssign=vecindClusterAssign,
      matMu=matMuH1,
      matDispersions=matDispersionsH1,
      vecSizeFactors=vecSizeFactors,
      matDropout=matDropoutH1,
      matZ=matZH1,
      matConstPredictorsPi=matConstPredictorsPi,
      matLinModelPi=matLinModelPi,
      matImpulseParam=matImpulseParam,
      scaWindowRadius=scaWindowRadius,
      boolDynamicPi=boolDynamicPi,
      boolBFGSEstOfMu=boolBFGSEstOfMu,
      strMuModel=strMuModel )
    matMuH1 <- lsMuFitsH1$matMu
    matImpulseParam <- lsMuFitsH1$matImpulseParam
    # These udates are done during mean estimation too
    # Reestimate drop-out rates based on new means
    matDropoutH1 <- do.call(cbind, bplapply(seq(1, scaNumCells), function(cell){
      vecLinModelOut <- cbind(1,log(matMuH1[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
      vecDropout <- 1/(1+exp(-vecLinModelOut))
      return(vecDropout)
    }))
    if(boolSuperVerbose){
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMuH1,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersionsH1, 
        matDropout=matDropoutH1,
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius )
      print(paste0("# ",scaIter, ".2) Mean estimation complete: ",
        "loglikelihood of       ", scaLogLikTemp))
    }
    
    # b) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood.
    if(boolOneDispPerGene){  
      vecDispFitH1 <- bplapply(seq(1,scaNumGenes), function(i){
        fitDispZINB_LinPulse(scaDispGuess=matDispersionsH1[i,1],
          vecCounts=matCountsProc[i,],
          vecMuEst=matMuH1[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropoutH1[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius,
          boolSmoothed=boolSmoothed )
      })
    } else {
      #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
      stop("Non constant-dispersion factors are not yet implemented. david.seb.fischer@gmail.com")
    }
    vecboolConvergedGLMdispH1 <- sapply(vecDispFitH1, function(fit) fit["convergence"])
    vecDispersionsH1 <- sapply(vecDispFitH1, function(fit){ fit["par"] })
    matDispersionsH1 <- matrix(vecDispersionsH1, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE) 
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    scaLogLikNew <- evalLogLikMatrix(matCounts=matCountsProc,
      matMu=matMuH1,
      vecSizeFactors=vecSizeFactors,
      matDispersions=matDispersionsH1, 
      matDropout=matDropoutH1, 
      matboolNotZeroObserved=matboolNotZeroObserved, 
      matboolZero=matboolZero,
      scaWindowRadius=scaWindowRadius )
    
    # Iteration complete
    if(boolSuperVerbose){
      if(any(vecboolConvergedGLMdispH1 != 0)){
        print(paste0("Dispersion estimation did not converge in ", 
          sum(vecboolConvergedGLMdispH1), " cases."))
      }
      print(paste0("# ",scaIter, ".3) Dispersion estimation complete: ",
        "loglikelihood of ", scaLogLikNew))
    } else {
      if(verbose){print(paste0("# ",scaIter, ".) complete with ",
        "log likelihood of ", scaLogLikNew))}
    }
    vecEMLogLikH1[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  print("Finished fitting alternative negative binomial model and noise model.")

  # Evaluate convergence
  if(all(as.logical(vecboolConvergedGLMdispH1)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergenceH1 <- TRUE
  } else { boolConvergenceH1 <- FALSE }

  print("### b) Fit null negative binomial model.")
  # Initialise dispersion factors to those estimated
  # for alternative model.
  matDispersionsH0 <- matDispersionsH1
  matDropoutH0 <- matDropoutH1
  
  # One mean per gene and one dispersion factor
  scaKbyGeneH0 <- 1+1
  
  scaIter <- 1
  scaLogLikNew <- NA
  scaLogLikOld <- NA
  
  while(scaIter == 1 | scaIter == 2 | 
      (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){
    ##### Gene-wise parameter estimation: 
    # a) Negative binomial mean parameter
    # Only compute posterior if using closed form estimator for mean:
    # Posterior is not necessary in all other cases, expect for parameter
    # initialisation of impulse model.
    lsMuFitsH0 <- fitZINBMu( matCountsProc=matCountsProc,
      matDispersions=matDispersionsH0,
      vecSizeFactors=vecSizeFactors,
      matDropout=matDropoutH0,
      matConstPredictorsPi=matConstPredictorsPi,
      matLinModelPi=matLinModelPi,
      scaWindowRadius=scaWindowRadius,
      boolDynamicPi=boolDynamicPi,
      strMuModel="constant" )
    matMuH0 <- lsMuFitsH0$matMu
    # These udates are done during mean estimation too
    # Reestimate drop-out rates based on new means
    matDropoutH0 <- do.call(cbind, bplapply(seq(1, scaNumCells), function(cell){
      vecLinModelOut <- cbind(1,log(matMuH0[,cell]),matConstPredictorsPi) %*% matLinModelPi[cell,]
      vecDropout <- 1/(1+exp(-vecLinModelOut))
      return(vecDropout)
    }))
    if(boolSuperVerbose){
      scaLogLikTemp <- evalLogLikMatrix(matCounts=matCountsProc,
        matMu=matMuH0,
        vecSizeFactors=vecSizeFactors,
        matDispersions=matDispersionsH0, 
        matDropout=matDropoutH0,
        matboolNotZeroObserved=matboolNotZeroObserved, 
        matboolZero=matboolZero,
        scaWindowRadius=scaWindowRadius )
      print(paste0("# ",scaIter, ".1) Mean estimation complete: ",
        "loglikelihood of       ", scaLogLikTemp))
    }
    
    # b) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood.
    if(boolOneDispPerGene){  
      vecDispFitH0 <- bplapply(seq(1,scaNumGenes), function(i){
        fitDispZINB_LinPulse(scaDispGuess=matDispersionsH0[i,1],
          vecCounts=matCountsProc[i,],
          vecMuEst=matMuH0[i,],
          vecSizeFactors=vecSizeFactors,
          vecDropoutRateEst=matDropoutH0[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          scaWindowRadius=scaWindowRadius,
          boolSmoothed=boolSmoothed )
      })
    } else {
      #  Not coded yet. Contact david.seb.fischer@gmail.com if desired.
      stop("Non constant-dispersion factors are not yet implemented. david.seb.fischer@gmail.com")
    }
    vecboolConvergedGLMdispH0 <- sapply(vecDispFitH0, function(fit) fit["convergence"])
    vecDispersionsH0 <- sapply(vecDispFitH0, function(fit){ fit["par"] })
    matDispersionsH0 <- matrix(vecDispersionsH0, nrow=scaNumGenes, ncol=scaNumCells, byrow=FALSE)
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    scaLogLikNew <- evalLogLikMatrix(matCounts=matCountsProc,
      matMu=matMuH0,
      vecSizeFactors=vecSizeFactors,
      matDispersions=matDispersionsH0, 
      matDropout=matDropoutH0, 
      matboolNotZeroObserved=matboolNotZeroObserved, 
      matboolZero=matboolZero,
      scaWindowRadius=scaWindowRadius )
    
    # EM-iteration complete
    if(boolSuperVerbose){
      if(any(vecboolConvergedGLMdispH0 != 0)){
        print(paste0("Dispersion estimation did not converge in ", 
          sum(vecboolConvergedGLMdispH0), " cases."))
      }
      print(paste0("# ",scaIter, ".2) Dispersion estimation complete: ",
        "loglikelihood of ", scaLogLikNew))
    } else {
      if(verbose){print(paste0("# ",scaIter, ".) complete with ",
        "log likelihood of ", scaLogLikNew))}
    }
    vecEMLogLikH0[scaIter] <- scaLogLikNew
    scaIter <- scaIter+1
  }
  print("Finished fitting alternative mean model.")
  
  # Evaluate convergence
  if(all(as.logical(vecboolConvergedGLMdispH0)) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergenceH0 <- TRUE
  } else { boolConvergenceH0 <- FALSE }
  
  # Name rows and columns of output
  rownames(matMuH1) <- rownames(matCountsProc)
  colnames(matMuH1) <- colnames(matCountsProc)
  rownames(matDispersionsH1) <- rownames(matCountsProc)
  colnames(matDispersionsH1) <- colnames(matCountsProc)
  rownames(matDropoutH1) <- rownames(matCountsProc)
  colnames(matDropoutH1) <- colnames(matCountsProc)
  rownames(matMuH0) <- rownames(matCountsProc)
  colnames(matMuH0) <- colnames(matCountsProc)
  rownames(matDispersionsH0) <- rownames(matCountsProc)
  colnames(matDispersionsH0) <- colnames(matCountsProc)
  rownames(matDropoutH0) <- rownames(matCountsProc)
  colnames(matDropoutH0) <- colnames(matCountsProc)
  rownames(matLinModelPi) <- colnames(matCountsProc)
  
  lsFitZINBReporters <- list( boolConvergenceH1=boolConvergenceH1,
    boolConvergenceH0=boolConvergenceH0,
    vecEMLogLikH1=vecEMLogLikH1,
    vecEMLogLikH0=vecEMLogLikH0,
    scaKbyGeneH1=scaKbyGeneH1,
    scaKbyGeneH0=scaKbyGeneH0 )
  
  return(list( matMuH1=matMuH1,
    matDispersionsH1=matDispersionsH1,
    matDropoutH1=matDropoutH1,
    matMuH0=matMuH0,
    matDispersionsH0=matDispersionsH0,
    matDropoutH0=matDropoutH0,
    matDropoutLinModel=matLinModelPi,
    matImpulseParam=matImpulseParam,
    lsFitZINBReporters=lsFitZINBReporters ))
}