################################################################################
#######################     LineagePulse package     ###########################
################################################################################

################################################################################
### Libraries and source code
################################################################################

library(BiocParallel)
library(compiler)
library(ggplot2)
setwd("~/gitDevelopment/LineagePulse/R")

source("srcLineagePulse_evalLogLikZINB.R")

source("srcLineagePulse_clusterCellsInPseudotime.R")
source("srcLineagePulse_calcPostDrop.R")
source("srcLineagePulse_calcNormConst.R")
source("srcLineagePulse_classLineagePulseObject.R")
source("srcLineagePulse_decompressParameters.R")
source("srcLineagePulse_estimateMMAssignments.R")
source("srcLineagePulse_evalDropoutModel.R")
source("srcLineagePulse_evalImpulseModel.R")
source("srcLineagePulse_fitZINB_cofitMeanDispersion.R")
source("srcLineagePulse_fitZINB_fitDropout.R")
source("srcLineagePulse_fitZINB.R")
source("srcLineagePulse_fitZINB_WrapperMixture.R")
source("srcLineagePulse_fitZINB_WrapperLP.R")
source("srcLineagePulse_initialiseImpulseParameters.R")
source("srcLineagePulse_initialiseCentroids.R")
source("srcLineagePulse_plotComparativeECDF.R")
source("srcLineagePulse_plotGene.R")
source("srcLineagePulse_plotPseudotimeClustering.R")
source("srcLineagePulse_processSCData.R")
source("srcLineagePulse_processSCDataMixture.R")
source("srcLineagePulse_runDEAnalysis.R")
source("srcLineagePulse_simulateDataSet.R")
source("srcLineagePulse_sortGeneTrajectories.R")
source("srcLineagePulse_validateOutput.R")
source("srcLineagePulse_validateOutputSimulation.R")

################################################################################
### Main function
################################################################################

#' LineagePulse wrapper: Differential expression analysis in pseudotime
#' 
#' This function performs all steps of longitudinal differential
#' expression analysis in pseudotime for you.
#' 
#' This function is the wrapper function for the LineagePulse algorithm
#' which performs differential expression analysis in pseudotime. 
#' Note that LineagePulse has many input parameters but
#' only few will be relevant for you and you will be able to leave the 
#' remaining ones as their defaults. Read up on specific input parameters
#' in the input parameter annotation of this function or follow this short
#' guide:
#' 
#' MINIMAL INPUT
#' 1. Supply data: Count data (matCounts) and pseudotime coordinates
#' of the cells (vecPseudotime). You may decide to also provide cell-wise
#' normalisation factors (such as factors accounting for sequencing depth)
#' (vecNormConstExternal).
#' 2. Supply cell meta data: Data frame (dfAnnotation) which contains
#' all cell ids (the colnames of matCounts) in a column "cell",
#' and cell-wise pseudotime coordinates (scalar) in a a column "pseudpotime".
#' Rownames must be the ids in column "cell".
#' 3. Chose the model constraining mean (strMuModel) and dispersion 
#' parameters (strDispModel) for each gene. If you run clusters,
#' you may decide to force the number of clusters (scaKCluster)
#' rather than using internal model selection or use clusters
#' based on true time (time of sampling) (boolClusterInPseudotime).
#' 
#' ADDITIONAL FACULTATIVE SETTINGS
#' 4. Decide whether you want to use local negative binomial model
#' smooting (scaWindowRadius). 
#' 5. Supply gene-specific drop-out predictors if wanted 
#' (matPiConstPredictors).
#' 6. Set optimisation parameters (boolEstimateNoiseBasedOnH0,
#' boolVecWindowsAsBFGS, boolCoEstDispMean, scaMaxEstimationCycles).
#' 7. Chose the number of processes you want to use (scaNProc), LineagePulse
#' is parallelised on all computation intensive steps. Note that
#' the current parallelisation scheme runs on Unix (MacOS) and Linux but
#' not on windows. Adjust the code section within this wrapper to
#' parallelise on windows.
#' 8. Set the level of detail with which you want to follow
#' progress through text printed on the console during a run
#' (verbose, boolSuperVerbose).
#' 
#' Finally, after running LineagePulse, you may continue to work
#' on your data set by:
#' A) Regenerating observation-wise parameters,
#' such as the mean parameter matrix which represents the hidden
#' expression states, with the functions in srcLineagePulse_decompressParameters.R.
#' B) You can also compute the observation-wise probability of 
#' dropout with calcPostDrop.
#' C) You can have a closer look at the global expression
#' trajecotries of the genes with sortGeneTrajectories.
#' D) You can look at gene-wise model fits with plotGene().
#' E) If you use strMuModel=="clusters", you can plot the
#' observed density of cells across pseudotime and the cluster 
#' boundaries in this 1D space with plotPseudotimeClustering.
#' 
#' @aliases LineagePulse wrapper, main function
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param dfAnnotation: (data frame cells x meta characteristics)
#'    Annotation table which contains meta data on cells.
#'    May contain the following columns
#'    cell: Cell IDs.
#'    pseudotime: Pseudotemporal coordinates of cells.
#'    Confounder1: Batch labels of cells with respect 
#'    to first confounder. Name is arbitrary: Could
#'    for example be "patient" with batch labels
#'    patientA, patientB, patientC.
#'    Confounder2: As Confounder1 for another confounding
#'    variable.
#'    ... ConfounderX.
#'    population: Fixed population assignments (for
#'    strMuModel="MM"). Cells not assigned have to be NA.
#'    clusters: External clustering results assigning each cell
#'    to one cluster ID. Used if strMuModel=="clusters" and this
#'    column is given, otherwise clustering is internally generated.
#'    rownames: Must be IDs from column cell.
#'    Remaining entries in table are ignored.
#' @param vecConfounders: Confounders to correct for in batch
#'    correction modeling, must be subset of column names of
#'    dfAnnotation which describe condounding variables.
#' @param strMuModel: (str) {"constant", "cluster", "MM",
#'    "windows","impulse"}
#'    [Default "impulse"] Model according to which the mean
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param strDispModel: (str) {"constant"}
#'    [Default "constant"] Model according to which dispersion
#'    parameter is fit to each gene as a function of 
#'    pseudotime in the alternative model (H1).
#' @param matPiConstPredictors: (numeric matrix genes x number of constant
#'    gene-wise drop-out predictors) Predictors for logistic drop-out 
#'    fit other than offset and mean parameter (i.e. parameters which
#'    are constant for all observations in a gene and externally supplied.)
#'    Is null if no constant predictors are supplied
#' @param vecNormConstExternal: (numeric vector number of cells) 
#'    Model scaling factors, one per cell. These factors will linearly 
#'    scale the mean model for evaluation of the loglikelihood. 
#'    Must be named according to the column names of matCounts.
#' @param scaKCluster: (integer) [Default NULL] Forces number of centroids
#'    in K-means to be K: setting this to an integer (not NULL) skips model
#'    selection in clusterting.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius of cells within pseudotemporal
#'    ordering. Each negative binomial model inferred on
#'    observation [gene i, cell j] is fit and evaluated on 
#'    the observations [gene i, cells in neighbourhood of j],
#'    the model is locally smoothed in pseudotime.
#' @param boolEstimateNoiseBasedOnH0: (bool) [Default: FALSE]
#'    Whether to co-estimate logistic drop-out model with the 
#'    constant null model or with the alternative model. The
#'    co-estimation with the noise model typically extends the
#'    run-time of this model-estimation step strongly. While
#'    the drop-out model is more accurate if estimated based on
#'    a more realistic model expression model (the alternative
#'    model), a trade-off for speed over accuracy can be taken
#'    and the dropout model can be chosen to be estimated based
#'    on the constant null expression model (set to TRUE).
#' @param boolVecWindowsAsBFGS: (bool) [Default FALSE] Whether
#'    mean parameters of a gene are co-estimated in "windows"
#'    mode with BFGS algorithm (optimisation with dimensionality
#'    of number of cells) or estimated one by one, conditioned
#'    one the latest estimates of neighbours. The latter case
#'    (boolVecWindowsAsBFGS=FALSE) is coordinate ascent within the gene
#'    and each mean parameter is optimised once only.
#' @param boolCoEstDispMean: (bool) [Default TRUE]
#'    Whether mean and dispersion parameters are to be co-estimated
#'    (simulatneous optimisation). Only available for certain 
#'    dispersion and mean models.
#'    Note that co-estimation in model estimation B (without drop-
#'    out model estimation) leads to a single step estimation as 
#'    mean and dispersion parameter don't have to be iterated over.
#'    This makes estimation of large data sets with complex H1 mean
#'    model (e.g. impulse) possible, as the drop-out model can be 
#'    estimated based on H0 (boolEstimateNoiseBasedOnH0) so that
#'    the complex model only has to be estimated once (simultaneous
#'    with the dispersion parameters). This may generally lead to better
#'    convergence as the steps in coordinate-ascent are in a larger
#'    space, closer to full gradient ascent. Setting to TRUE
#'    is encouraged. Optimisation routines for individual mean 
#'    and dispersion fitting (if FALSE) exist, but these may be viewed
#'    as non-deprecated parts of an earlier implementation of the
#'    alorithm.
#' @param scaMaxEstimationCycles: (integer) [Default 20] Maximum number 
#'    of estimation cycles performed in fitZINB(). One cycle
#'    contain one estimation of of each parameter of the 
#'    zero-inflated negative binomial model as coordinate ascent.
#' @param scaNProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param verbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation with one report per cycle.
#' @param boolSuperVerbose: (bool) Whether to follow convergence of the 
#'    iterative parameter estimation in high detail with local 
#'    convergence flags and step-by-step loglikelihood computation.
#' 
#' @return dfDEAnalysis: (data frame genes x reported variables) 
#'    Summary of differential expression analysis, sorted by adj.p:
#'    {Gene: gene ID,
#'    p: raw p-value, 
#'    adj.p: BH corrected p-value, 
#'    loglik_full: loglikelihood of alternative model H1,
#'    loglik_red: loglikelihood of null model H0,
#'    deviance: loglikelihood ratio test statistic (the deviance),
#'    mean_H0: inferred gene-wise mean parameter (constant null model),
#'    dispersion_H0: inferred gene-wise dispersion parameter (constant null model)}
#'    
#' @author David Sebastian Fischer
#' 
#' @export
runLineagePulse <- function(matCounts,
														dfAnnotation,
														vecConfounders=NULL,
														strMuModel="impulse",
														strDispModel="constant",
														matPiConstPredictors=NULL,
														vecNormConstExternal=NULL,
														scaKClusters=NULL,
														scaWindowRadius=NULL,
														boolEstimateNoiseBasedOnH0=FALSE,
														boolVecWindowsAsBFGS=FALSE,
														boolCoEstDispMean=TRUE,
														scaMaxEstimationCycles=20,
														scaNProc=1,
														boolVerbose=TRUE,
														boolSuperVerbose=FALSE ){
	
  STR_VERSION <- "v0.99"
  
  # 1. Data preprocessing
	vecAllGenes <- rownames(matCounts)
	lsProcessedSCData <- processSCData( matCounts=matCounts,
																			dfAnnotation=dfAnnotation,
																			vecConfounders=vecConfounders,
																			matPiConstPredictors=matPiConstPredictors,
																			vecNormConstExternal=vecNormConstExternal,
																			strMuModel=strMuModel,
																			strDispModel=strDispModel,
																			scaWindowRadius=scaWindowRadius,
																			boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
																			boolCoEstDispMean=boolCoEstDispMean,
																			scaMaxEstimationCycles=scaMaxEstimationCycles,
																			boolVerbose=boolVerbose,
																			boolSuperVerbose=boolSuperVerbose,
																			STR_VERSION=STR_VERSION)
	objectLineagePulse <- lsProcessedSCData$objectLineagePulse
	vecNormConstExternalProc <- lsProcessedSCData$vecNormConstExternalProc
	matPiConstPredictorsProc <- lsProcessedSCData$matPiConstPredictorsProc
	
	# Clear memory
	rm(matCounts)
	rm(matPiConstPredictors)
	rm(dfAnnotation)
	rm(lsProcessedSCData)
	
	# Inialise parallelisation
	# Set the parallelisation environment in BiocParallel:
	if(scaNProc > 1){
		register(MulticoreParam(workers=scaNProc)) 
		#register(SnowParam(workers=scaNProc, timeout=60*60*24*7))
	} else {
		# For debugging in serial mode
		register(SerialParam())
	}
	
	if(strMuModel=="clusters"){
	  # 2. Cluster cells in pseudo-time
	  strMessage <- paste0("--- Run/read clustering.")
	  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	  if(boolVerbose) print(strMessage)
	  
	  tm_clustering <- system.time({
	    if(is.null(objectLineagePulse@dfAnnotationProc$clusters)){
	      # Internal clustering
	      if(!is.null(scaKClusters)){
	        strMessage <- paste0("Use internally created clustering (K-means) with ",
	                                                     scaKClusters, " clusters.")
	      } else {
	        strMessage <- paste0("Use internally created clustering (K-means) with ",
	                             " internal model selection to select number of clusters.")
	      }
	      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	      if(boolSuperVerbose) print(strMessage)
	      
	      objectLineagePulse@dfAnnotationProc$clusters <- clusterCellsInPseudotime(vecPseudotime=objectLineagePulse@dfAnnotationProc$pseudotime,
	                                                                               scaKexternal=scaKClusters)
	      
	      strMessage <- paste0("Chose the number of clusters K as ", 
	                           max(objectLineagePulse@dfAnnotationProc$clusters),
	                           " based on the Gap-statistic.")
	      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	      if(boolSuperVerbose) print(strMessage)
	    } else {
	      # External clustering
	      strMessage <- paste0("Use externally created clustering (K-means) supplied through ",
	                           "dfAnnotation$clusters with ",
	                           length(unique(objectLineagePulse@dfAnnotationProc$clusters)),
	                           " clusters.")
	      objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	      if(boolSuperVerbose) print(strMessage)
	    }
	  })
	  strMessage <- paste0("Time elapsed during clustering: ",round(tm_clustering["elapsed"]/60,2), " min")
	  objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	  if(boolVerbose) print(strMessage)
	}
	
	# 3. Compute normalisation constants
	strMessage <- paste0("--- Compute normalisation constants:")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	if(boolVerbose) print(strMessage)
	
	objectLineagePulse <- calcNormConst(objectLineagePulse=objectLineagePulse,
																			vecNormConstExternal=vecNormConstExternalProc)
	
	# 5. Fit ZINB model for both H1 and H0.
	strMessage <- paste0("--- Fit ZINB model for both H1 and H0.")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	if(boolVerbose) print(strMessage)
	
	tm_fitmm <- system.time({
		objectLineagePulse <- fitNullAlternative( objectLineagePulse=objectLineagePulse,
																							matPiConstPredictors=matPiConstPredictorsProc,
																							boolEstimateNoiseBasedOnH0=boolEstimateNoiseBasedOnH0,
																							boolVecWindowsAsBFGS=boolVecWindowsAsBFGS,
																							strMuModel=strMuModel,
																							strDispModel=strDispModel,
																							scaMaxEstimationCycles=scaMaxEstimationCycles,
																							boolVerbose=boolVerbose,
																							boolSuperVerbose=boolSuperVerbose )
	})
	strMessage <- paste0("Time elapsed during ZINB fitting: ",
	                     round(tm_fitmm["elapsed"]/60,2)," min")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	if(boolVerbose) print(strMessage)
	
	
	# 6. Differential expression analysis:
	strMessage <- paste0("--- Run differential expression analysis.")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	if(boolVerbose) print(strMessage)
	
	tm_deanalysis_mf <- system.time({
		objectLineagePulse <- runDEAnalysis( objectLineagePulse=objectLineagePulse )
	})
	strMessage <- paste0("Time elapsed during differential expression analysis: ",
	                     round(tm_deanalysis_mf["elapsed"]/60,2)," min")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage, "\n")
	if(boolVerbose) print(strMessage)
	
	strMessage <- paste0("Finished runLineagePulse().")
	objectLineagePulse@strReport <- paste0(objectLineagePulse@strReport, strMessage)
  if(boolVerbose) print(strMessage)
  
	return(objectLineagePulse)
}