#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++    Model-free DE analysis  ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Model-free differential expression analysis
#' 
#' Performs model-free differential expression analysis. This functiion is
#' broken up into:
#' (I) Fit null model.
#' (II) Compute likelihood of full and reduced model.
#' (II) Differential expression analysis.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param vecDispersions: (numeric matrix genes x clusters)
#'    Inferred negative binomial dispersions.
#' @param matMuCluster: (numeric matrix genes x clusters)
#'    Inferred negative binomial cluster means.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param boolConvergenceZINBH1: (bool) Convergence of EM algorithm for
#'    full model.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes) 
#'    Summary of model-free differential expression analysis.
#' @export

runModelFreeDEAnalysis <- function(matCountsProc,
  vecPseudotime,
  vecSizeFactors,
  lsResultsClusteringH1,
  vecDispersionsH1,
  matMuClusterH1,
  matDropoutH1,
  boolConvergenceZINBH1,
  scaWindowRadius=NULL,
  boolOneDispPerGene=TRUE,
  nProc=1,
  scaMaxiterEM = 100,
  verbose=TRUE ){
  
  # (I) Fit null model
  # Null model is a single cluster: Create Clustering.
  lsH0Clustering <- list()
  lsH0Clustering$Assignments <- rep(1,dim(matCountsProc)[2])
  lsH0Clustering$Centroids <- mean(vecPseudotime, na.rm=TRUE)
  lsH0Clustering$K <- 1
  
  # Fit zero-inflated negative binomial null model
  lsResZINBFitsH0 <- fitZINB( matCountsProc=matCountsProc, 
    lsResultsClustering=lsH0Clustering,
    vecSizeFactors=vecSizeFactors,
    vecSpikeInGenes=NULL,
    boolOneDispPerGene=boolOneDispPerGene,
    scaWindowRadius=NULL,
    nProc=nProc,
    scaMaxiterEM=scaMaxiterEM,
    verbose=verbose )
  vecMuClusterH0  <- lsResZINBFitsH0$matMuCluster
  vecDispersionsH0 <- lsResZINBFitsH0$vecDispersions
  matDropoutH0 <- lsResZINBFitsH0$matDropout
  boolConvergenceZINBH0 <- lsResZINBFitsH0$boolConvergence
  matMuH0 <- matrix(vecMuClusterH0, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  matDispersionsH0 <- matrix(vecDispersionsH0, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  
  # (II) Fit alternative model
  # Alternative model was fit furing hyperparameter estimation
  # if clusters were used for mean estimation during hyperparameter
  # estimation.
  if(is.null(scaWindowRadius)){
    matMuH1 <- matMuClusterH1[,lsResultsClusteringH1$Assignments]
    matDispersionsH1 <- matrix(vecDispersionsH1, nrow=dim(matCountsProc)[1], 
      ncol=dim(matCountsProc)[2], byrow=FALSE)
  } else {
    # Fit zero-inflated negative binomial alternative model
    # without sliding mean estimation but with mean estimation
    # by cluster.
    lsResZINBFitsH1 <- fitZINB( matCountsProc=matCountsProc, 
      lsResultsClustering=lsResultsClusteringH1,
      vecSizeFactors=vecSizeFactors,
      vecSpikeInGenes=NULL,
      boolOneDispPerGene=boolOneDispPerGene,
      scaWindowRadius=NULL,
      nProc=nProc,
      scaMaxiterEM=scaMaxiterEM,
      verbose=verbose )
    vecMuClusterH1  <- lsResZINBFitsH1$matMuCluster
    vecDispersionsH1 <- lsResZINBFitsH1$vecDispersions
    matDropoutH1 <- lsResZINBFitsH1$matDropout
    boolConvergenceZINBH1 <- lsResZINBFitsH1$boolConvergence
    matMuH1 <- matrix(vecMuClusterH1, nrow=dim(matCountsProc)[1], 
      ncol=dim(matCountsProc)[2], byrow=FALSE)
    matDispersionsH1 <- matrix(vecDispersionsH1, nrow=dim(matCountsProc)[1], 
      ncol=dim(matCountsProc)[2], byrow=FALSE)
  }
  
  # (III) Compute log likelihoods
  matboolNotZeroObserved <- matCountsProc >0 & !is.na(matCountsProc)
  matboolZero <- matCountsProc==0
  vecLogLikFull <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_PseudoDE_comp(vecY=matCountsProc[i,],
      vecMu=matMuH1[i,]*vecSizeFactors,
      vecDispEst=matDispersionsH1[i,], 
      vecDropoutRateEst=matDropoutH1[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  vecLogLikRed <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_PseudoDE_comp(vecY=matCountsProc[i,],
      vecMu=matMuH0[i,]*vecSizeFactors,
      vecDispEst=matDispersionsH0[i,], 
      vecDropoutRateEst=matDropoutH0[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  
  # (IV) Differential expression analysis
  # scaK: Number of clusters used in full model
  scaK <- lsResultsClusteringH1$K
  # Dropout rate models are shared and do not influence
  # the difference in degrees of freedom of the model.
  # Full model: K x mean and 1 or K x dispersion parameters.
  if(boolOneDispPerGene){ scaDegFreedomFull <- scaK + 1
  }else{ scaDegFreedomFull <- scaK*2 }
  # Null model: One dispersion estimate and one overall mean estimate
  scaDegFreedomRed <- 2
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDeltaDegFreedom <- scaDegFreedomFull - scaDegFreedomRed
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH = p.adjust(vecPvalue, method = "BH")
  
  dfModelFreeDEAnalysis =   as.data.frame(cbind(
    "Gene" = row.names(matCountsProc),
    "p"=as.numeric(vecPvalue),
    "adj.p"=as.numeric(vecPvalueBH),
    "loglik_full"=vecLogLikFull,
    "loglik_red"=vecLogLikRed,
    "deviance"=vecDeviance,
    "mean"=vecMuClusterH0,
    "dispersion_H0"=vecDispersionsH0,
    "converged_H0"=rep(boolConvergenceZINBH0!=0,dim(matCountsProc)[1]),
    "converged_full"=rep(all(boolConvergenceZINBH1),dim(matCountsProc)[1]),
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}