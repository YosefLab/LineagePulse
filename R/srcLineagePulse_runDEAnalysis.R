#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++   Differential expression analysis  +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Differential expression analysis
#' 
#' Performs differential expression analysis based on previously estimated 
#' null and alternative models.
#' (I) Compute loglikelihood of data under null H0 and alternative H1 model.
#' (II) Differential expression analysis as loglikelihood ratio test.
#' 
#' @seealso Called by \code{runLineagePulse}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param lsMuModelH1: (list length 2)
#'    All objects necessary to compute H1 mean parameters for all
#'    observations.
#'    \itemize{
#'      \item matMuModel: (numerical matrix genes x number of model parameters)
#'      Parameters of mean model for each gene.
#'      \item lsMuModelGlobal: (list) Global variables for mean model,
#'      common to all genes.
#'      \itemize{
#'        \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'      "windows"} Name of the mean model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'        \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'      of a gene are simultaneously estiamted as a vector with BFGS
#'      in windows mode.
#'        \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'      for BFGS estimation of impulse model with optim (termination criterium).
#'        \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'      change in objective function for BFGS estimation of impulse 
#'      model with optim (termination criterium).
#'      }
#'    }
#' @param lsDispModelH1: (list length 2)
#'    All objects necessary to compute H1 dispersion parameters for all
#'    observations.
#'    \itemize{
#'      \item matDispModel: (numerical matrix genes x number of model parameters)
#'    Parameters of dispersion model for each gene.
#'      \item lsDispModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'      \itemize{
#'        \item strDispModel: (str) {"constant"} 
#'      Name of the dispersion model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'      }
#'    }
#' @param lsMuModelH0: (list length 2)
#'    All objects necessary to compute H0 mean parameters for all
#'    observations.
#'    \itemize{
#'      \item matMuModel: (numerical matrix genes x number of model parameters)
#'      Parameters of mean model for each gene.
#'      \item lsMuModelGlobal: (list) Global variables for mean model,
#'      common to all genes.
#'      \itemize{
#'        \item strMuModel: (str) {"constant", "impulse", "clusters", 
#'      "windows"} Name of the mean model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'        \item boolVecWindowsAsBFGS: (bool) Whether mean parameters
#'      of a gene are simultaneously estiamted as a vector with BFGS
#'      in windows mode.
#'        \item MAXIT_BFGS_Impulse: (int) Maximum number of iterations
#'      for BFGS estimation of impulse model with optim (termination criterium).
#'        \item RELTOL_BFGS_Impulse: (scalar) Relative tolerance of
#'      change in objective function for BFGS estimation of impulse 
#'      model with optim (termination criterium).
#'      }
#'    }
#' @param lsDispModelH0: (list length 2)
#'    All objects necessary to compute H0 dispersion parameters for all
#'    observations.
#'    \itemize{
#'      \item matDispModel: (numerical matrix genes x number of model parameters)
#'    Parameters of dispersion model for each gene.
#'      \item lsDispModelGlobal: (list) Global variables for mean model,
#'    common to all genes.
#'      \itemize{
#'        \item strDispModel: (str) {"constant"} 
#'      Name of the dispersion model.
#'        \item scaNumCells: (scalar) [Default NA] Number of cells
#'      for which model is evaluated. Used for constant model.
#'        \item vecPseudotime: (numerical vector number of cells)
#'      [Default NA] Pseudotime coordinates of cells. Used for
#'      impulse model.
#'        \item vecindClusterAssign: (integer vector length number of
#'      cells) [Default NA] Index of cluster assigned to each cell.
#'      Used for clusters model.
#'      }
#'    }
#' @param lsDropModel: (list length 2)
#'    All objects necessary to compute drop-out parameters for all
#'    observations, omitting mean parameters (which are stored in lsMeanModel).
#'    \itemize{
#'      \item matDropoutLinModel: (numeric matrix cells x number of model parameters)
#'      {offset parameter, log(mu) parameter, parameters belonging to
#'      constant predictors}
#'      Parameters of drop-out model for each cell
#'      \item matPiConstPredictors: (numeric matrix genes x number of constant
#'      gene-wise drop-out predictors) Predictors for logistic drop-out 
#'      fit other than offset and mean parameter (i.e. parameters which
#'      are constant for all observations in a gene and externally supplied.)
#'      Is null if no constant predictors are supplied.
#'    }
#' @param scaKbyGeneH1: (scalar) Degrees of freedom
#'    by gene used in null model H1. The logistic
#'    dropout model is ignored as it is shared between
#'    alternative and null model.
#' @param scaKbyGeneH0: (scalar) Degrees of freedom
#'    by gene used in null model H0. The logistic
#'    dropout model is ignored as it is shared between
#'    alternative and null model.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes x reported variables) 
#'    Summary of differential expression analysis, sorted by adj.p:
#'    Gene: gene ID,
#'    p: raw p-value, 
#'    adj.p: FDR corrected p-value, 
#'    loglik_full: loglikelihood of alternative model H1,
#'    loglik_red: loglikelihood of null model H0,
#'    deviance: loglikelihood ratio test statistic (the deviance),
#'    mean_H0: inferred gene-wise mean parameter (constant null model),
#'    dispersion_H0: inferred gene-wise dispersion parameter (constant null model)
#'    
#' @author David Sebastian Fischer
#' 
#' @export

runDEAnalysis <- function(matCountsProc,
  vecNormConst,
  lsMuModelH1,
  lsDispModelH1,
  lsMuModelH0,
  lsDispModelH0,
  lsDropModel,
  scaKbyGeneH1,
  scaKbyGeneH0,
  scaWindowRadius=NULL ){
  
  # (I) Compute log likelihoods
  lsLL <- bplapply( seq(1,dim(matCountsProc)[1]), function(i){
    vecCounts <- matCountsProc[i,]
    vecboolNotZero <- !is.na(vecCounts) & vecCounts>0
    vecboolZero <- !is.na(vecCounts) & vecCounts==0
    
    # Decompress parameters by gene H0
    vecMuParamH0 <- decompressMeansByGene( vecMuModel=lsMuModelH0$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH0$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH0 <- decompressDispByGene( vecDispModel=lsDispModelH0$matDispModel[i,],
      lsDispModelGlobal=lsDispModelH0$lsDispModelGlobal,
      vecInterval=NULL )
    vecPiParamH0 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH0,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
    
    # Decompress parameters by gene H1
    vecMuParamH1 <- decompressMeansByGene( vecMuModel=lsMuModelH1$matMuModel[i,],
      lsMuModelGlobal=lsMuModelH1$lsMuModelGlobal,
      vecInterval=NULL )
    vecDispParamH1 <- decompressDispByGene( vecDispModel=lsDispModelH1$matDispModel[i,],
      lsDispModelGlobal=lsDispModelH1$lsDispModelGlobal,
      vecInterval=NULL )
    vecPiParamH1 <- decompressDropoutRateByGene( matDropModel=lsDropModel$matDropoutLinModel,
      vecMu=vecMuParamH1,
      vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,] )
    
    # Compute loglikelihoods
    scaLogLikFull <- evalLogLikGene(vecCounts=vecCounts,
      vecMu=vecMuParamH1,
      vecNormConst=vecNormConst,
      vecDisp=vecDispParamH1, 
      vecPi=vecPiParamH1,
      vecboolNotZero=vecboolNotZero, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
    
    scaLogLikRed <- evalLogLikGene(vecCounts=vecCounts,
      vecMu=vecMuParamH0,
      vecNormConst=vecNormConst,
      vecDisp=vecDispParamH0, 
      vecPi=vecPiParamH0,
      vecboolNotZero=vecboolNotZero, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
    
    return(list( scaLogLikFull=scaLogLikFull,
      scaLogLikRed=scaLogLikRed ))
  })
  vecLogLikFull <- sapply(lsLL, function(LL_i) LL_i[["scaLogLikFull"]])
  vecLogLikRed <- sapply(lsLL, function(LL_i) LL_i[["scaLogLikRed"]])
  
  # (II) Differential expression analysis
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDeltaDegFreedom <- scaKbyGeneH1 - scaKbyGeneH0
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
    "mean_H0"=array(lsMuModelH0$matMuModel),
    "dispersion_H0"=array(lsDispModelH0$matDispModel),
    stringsAsFactors = FALSE))
  rownames(dfModelFreeDEAnalysis) <- dfModelFreeDEAnalysis$Gene
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}