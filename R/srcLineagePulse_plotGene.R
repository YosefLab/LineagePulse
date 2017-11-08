#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  Plot counts and model for one gene  ++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot counts and model for one gene
#' 
#' Plot counts and model for one gene as a scatter plot with model trajectories:
#' Both as counts or log10 counts vs continuous covariate .
#' Dropout rates can be given and are visualised as the colour of the observation
#' points.
#' A reference model trajectory can be added.
#' Contour plots: 
#' A contour plot aids the interpretation of the model trajectory
#' in context of the observed data if local overplotting occurs 
#' because of a non-uniform distribution of cells in the
#' continuous covariate (a very common scenario in pseudotime for example).
#' In such a case, one might mistake cells with high counts from an
#' interval with high cellular intensity as a indication of increased
#' expression in this interval. 
#' Instead, one can explain such local maxima based on the 
#' simple phenomenom that more samples were taken from the underyling 
#' distribution in that continuous covariate interval 
#  which increases the probability of sampling far in the tail of the 
#' distribution. 
#' To migitate this visually, we added the option to add "contours": 
#' A contour is an alternative to a confidence interval
#' which includes information on the local sampling density. 
#' The contour is specific to a number of cells n. 
#' It is an expression trajectory in continuous covariate which lies
#' at the line above which n cells are expected within 
#' a local continuous covariate bin given its sampling density and 
#' the H1 non-constant expression model. 
#' Contours therefore quantify the increase in maximum observed 
#' counts in continuous covariate intervals due to sampling density, 
#' which is not due to the expression model. 
#' Note that this effect is corrected for in the model fitting and 
#' that the contour plots just aids visual interpretation of 
#' the fit to the data!
#' 
#' @param objLP (LineagePulseObject) LineagePulseObject
#' base plot on.
#' @param strGeneID (str) Name of gene, used for title of plot.
#' @param vecReferenceMuParam (numeric vector length number of cells)
#' [Default NULL] Reference mean trajectory which can be plotted
#' @param strTitleSuffix (str) String to be added to title.
#' @param boolLogPlot (bool) [Default TRUE]
#' Whether to log transform y-axis.
#' @param boolColourByDropout (bool) [Default TRUE]
#' Whether to colour scatter plot by posterior of drop-out.
#' @param boolH1NormCounts (bool) [Default FALSE] 
#' Whether to show normalised counts
#' (size factors and H1 batch factor estimates) as oppose to raw counts.
#' @param boolLineageContour (bool) [Default FALSE]
#' Whether to the "lineage contour" lines to the scatter plot.
#' @param boolTime (bool) [Default TRUE]
#' Whether continuous covariate is time, this simplifies the scatter
#' plot strongly. Show mean expression per time point as orange point
#' and 25\% and 75\% quantile of observations per time point as error bars.
#' @param bwDensity (bandwith numeric or string) [Default NULL]
#' Bandwith to be used to kernel density smooting
#' of cell density in continuous covariate 
#' (used if boolLineageContour=TRUE).
#' If not set, defaults to stats:density() default.
#' @param scaGgplot2Size (scalar) size in ggplot2 scatter
#' @param scaGgplot2Alpha (scalar) alpha in ggplot2 scatter
#' 
#' @return gplotGene (ggplot object)
#' Model rajectories and scatter plot for given gene.
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 100,
#'     scaNConst = 10,
#'     scaNLin = 10,
#'     scaNImp = 10,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 30),
#'     vecGeneWiseDropoutRates = rep(0.1, 30))
#' matDropoutPredictors <- as.matrix(data.frame(
#'     log_means = log(rowMeans(lsSimulatedData$counts)+1) ))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "splines", scaDFSplinesMu = 6,
#'     strDropModel="logistic", 
#'     matPiConstPredictors = matDropoutPredictors)
#' gplotExprProfile <- plotGene(
#'     objLP = objLP,
#'     strGeneID = rownames(lsSimulatedData$counts)[1],
#'     boolLineageContour = FALSE)
#' #print(gplotExprProfile)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
plotGene <- function(
    objLP,
    strGeneID,
    vecReferenceMuParam=NULL,
    strTitleSuffix=NULL,
    boolLogPlot=TRUE,
    boolColourByDropout=TRUE,
    boolH1NormCounts=FALSE,
    boolLineageContour=FALSE,
    boolTime=FALSE,
    bwDensity = NULL,
    scaGgplot2Size = 0.5,
    scaGgplot2Alpha = 0.5){
    
    cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                    "#D55E00", "#CC79A7", "grey", "red", "pink")
    ### 0. Check input
    ## Can only use one colour scale!
    if(boolLineageContour & boolColourByDropout){
        warning("Only one colour scale can be used. ",
                "Set either boolLineageContour or",
                " boolColourByDropout to FALSE")
        boolColourByDropout <- FALSE
    }
    if(!lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
       c("splines", "impulse") & boolTime) {
        stop("Selecting the time simplification for",
             "continuous models (boolTime=TRUE) does",
             "not make sense if strMuModel is not either",
             "splines or impulse.")
    }
    ### 1. Extract data and models
    vecCounts <- matCountsProc(objLP)[strGeneID,,sparse=FALSE]
    matCountsGene <- t(as.matrix(vecCounts))
    rownames(matCountsGene) <- strGeneID
    if(boolH1NormCounts){
        vecCountsToPlot <- as.vector(getNormData(
            matCounts=matCountsGene,
            lsMuModel=lsMuModelH1(objLP)))
    } else {
        vecCountsToPlot <- vecCounts
    }
    vecMuParamH0 <- decompressMeansByGene(
        vecMuModel=lsMuModelH0(objLP)$matMuModel[strGeneID,],
        lsvecBatchModel=NULL,
        lsMuModelGlobal=lsMuModelH0(objLP)$lsMuModelGlobal,
        vecInterval=NULL )
    vecMuParamH1 <- decompressMeansByGene(
        vecMuModel=lsMuModelH1(objLP)$matMuModel[strGeneID,],
        lsvecBatchModel=NULL,
        lsMuModelGlobal=lsMuModelH1(objLP)$lsMuModelGlobal,
        vecInterval=NULL )
    vecMuParamH1_NB <- decompressMeansByGene(
        vecMuModel=lsMuModelH1_NB(objLP)$matMuModel[strGeneID,],
        lsvecBatchModel=NULL,
        lsMuModelGlobal=lsMuModelH1_NB(objLP)$lsMuModelGlobal,
        vecInterval=NULL )
    vecDropPosterior <- as.vector(calcPostDrop_Matrix(
        matCounts=matCountsGene,
        lsMuModel=lsMuModelH1(objLP),
        lsDispModel=lsDispModelH1(objLP), 
        lsDropModel=lsDropModel(objLP),
        vecIDs=strGeneID))
    if(boolLogPlot){
        vecCountsToPlotForTime <- vecCountsToPlot 
        # Keep as log these as mean has to be taken before log transform
        # if boolTime=TRUE is used.
        vecCountsToPlot <- log(vecCountsToPlot+1)/log(10)
        vecMuParamH0 <- log(vecMuParamH0+1)/log(10)
        vecMuParamH1 <- log(vecMuParamH1+1)/log(10)
        vecMuParamH1_NB <- log(vecMuParamH1_NB+1)/log(10)
    }
    
    ### 2. Data scatter plot
    # Set drop-out rates as constant for visualistion if not given.
    if(!boolTime) {
        dfScatterCounts <- data.frame( counts=vecCountsToPlot )
        if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
           c("splines", "impulse")){
            dfScatterCounts$x <- lsMuModelH1(objLP)$lsMuModelGlobal$vecContinuousCovar
            
            if(boolColourByDropout){
                dfScatterCounts$dropout_posterior <- vecDropPosterior
                gplotGene <- ggplot() +
                    geom_point(data=dfScatterCounts, aes(
                        x=x, y=counts, colour=dropout_posterior), 
                        size = scaGgplot2Size, alpha = scaGgplot2Alpha, 
                        show.legend=TRUE)
            } else {
                gplotGene <- ggplot() +
                    geom_point(data=dfScatterCounts, aes(
                        x=x, y=counts), size = scaGgplot2Size,  
                        alpha = scaGgplot2Alpha, show.legend=TRUE)
            }
            
        } else if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% c("groups")){
            dfScatterCounts$x <- seq_len(lsMuModelH1(objLP)$lsMuModelGlobal$scaNumCells)
            dfScatterCounts$groups <- dfAnnotationProc(objLP)$groups
            
            if(boolColourByDropout){
                dfScatterCounts$dropout_posterior <- vecDropPosterior
                gplotGene <- ggplot() +
                    geom_point(data=dfScatterCounts, aes(
                        x=x, y=counts, colour=dropout_posterior, shape = groups), 
                        size = scaGgplot2Size,  alpha = scaGgplot2Alpha, 
                        show.legend=TRUE) +
                    scale_colour_gradient(high="red",low="green",limits=c(0, 1)) 
            } else {
                gplotGene <- ggplot() +
                    geom_point(data=dfScatterCounts, aes(
                        x=x, y=counts, shape = groups), 
                        size = scaGgplot2Size,  alpha = scaGgplot2Alpha, 
                        show.legend=TRUE)
            }
        }
    } else {
        # simplify scatter plot if time is continuous covariate
        # because many observatinos falls on the same covariate 
        # value.
        matQuantilesByTp <- do.call(rbind, tapply(
            vecCountsToPlotForTime, dfAnnot$continuous, quantile))
        vecMeanCount <- tapply(vecCountsToPlotForTime, dfAnnot$continuous, 
                               mean, na.rm=TRUE)
        if(boolLogPlot) {
            vecMeanCount <- log(vecMeanCount + 1)
            matQuantilesByTp <- log(matQuantilesByTp + 1)
        }
        vecTime <- tapply(dfAnnot$continuous, dfAnnot$continuous, 
                               unique)
        dfScatterCounts <- data.frame( 
            mean_count=vecMeanCount,
            time=vecTime,
            quantile_0=matQuantilesByTp[,1],
            quantile_25=matQuantilesByTp[,2],
            quantile_50=matQuantilesByTp[,3],
            quantile_75=matQuantilesByTp[,4],
            quantile_100=matQuantilesByTp[,5])
        gplotGene <- ggplot() + geom_point(data=dfScatterCounts, aes(
            x=time, y=mean_count), colour="orange") +
            geom_errorbar(data = dfScatterCounts, aes(
                   x=time, ymin=quantile_25, ymax=quantile_75))
    }
    
    # Set plotting threshold based on observed data
    scaMaxPlot <- 2 * max(vecCounts)
    
    ### 3. Add models to plot
    vecMuParamH0[vecMuParamH0 > scaMaxPlot] <- NA
    vecMuParamH1[vecMuParamH1 > scaMaxPlot] <- NA
    vecMuParamH1_NB[vecMuParamH1_NB > scaMaxPlot] <- NA
    if(is.null(vecReferenceMuParam)){
        if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
           c("splines", "impulse")){
            dfLineImpulse <- data.frame(
                x=rep(lsMuModelH1(objLP)$lsMuModelGlobal$vecContinuousCovar, 3),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecMuParamH1_NB),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("H1(NB)", length(vecMuParamH1_NB))) )
        } else if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
                  c("groups")){
            dfLineImpulse <- data.frame(
                x=rep(seq_len(lsMuModelH1(objLP)$lsMuModelGlobal$scaNumCells), 3),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecMuParamH1_NB),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("H1(NB)", length(vecMuParamH1_NB))) )
        }
    } else {
        if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
           c("splines", "impulse")){
            dfLineImpulse <- data.frame(
                x=rep(lsMuModelH1(objLP)$lsMuModelGlobal$vecContinuousCovar, 4),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecMuParamH1_NB,
                         vecReferenceMuParam),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("H1(NB)", length(vecMuParamH1_NB)),
                        rep("Reference", length(vecReferenceMuParam))) )
        } else if(lsMuModelH1(objLP)$lsMuModelGlobal$strMuModel %in% 
                  c("groups")){
            dfLineImpulse <- data.frame(
                x=rep(seq_len(lsMuModelH1(objLP)$lsMuModelGlobal$scaNumCells), 4),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecMuParamH1_NB,
                         vecReferenceMuParam),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("H1(NB)", length(vecMuParamH1)),
                        rep("Reference", length(vecReferenceMuParam))) )
        }
    }
    gplotGene <- gplotGene + geom_line(data=dfLineImpulse, aes(
        x=x, y=counts, linetype=model), size = 1, show.legend=TRUE)
    if(boolH1NormCounts) {
        gplotGene <- gplotGene + ylab("normalised counts")
    } else {
        gplotGene <- gplotGene + ylab("counts")
    }
    
    if(boolLineageContour) {
        # decompress log free
        if(boolLogPlot) {
            vecMuParamH1_nonlog <- 10^(vecMuParamH1)
        } else {
            vecMuParamH1_nonlog <- vecMuParamH1
        }
        vecDispParamH1 <- decompressDispByGene(
            vecDispModel=lsDispModelH1(objLP)$matDispModel[strGeneID,],
            lsvecBatchModel=NULL,
            lsDispModelGlobal=lsDispModelH1(objLP)$lsDispModelGlobal,
            vecInterval=NULL )
        # copmute density estimate of cells in continuous covariate
        if(is.null(bwDensity)) {
            bwDensity <- "nrd0"
        }
        lsDensity <- density(x = dfAnnotationProc(objLP)$continuous, bw=bwDensity, 
                             n = length(dfAnnotationProc(objLP)$continuous))
        scaNKDEbin <- 10
        scaWindowRadPT <- (max(dfAnnotationProc(objLP)$continuous, 
                               na.rm = TRUE) - 
                               min(dfAnnotationProc(objLP)$continuous, 
                                   na.rm = TRUE)) / (2*scaNKDEbin)
        vecObsInBin <- sapply(dfAnnotationProc(objLP)$continuous, function(pt) {
            sum(abs(pt-dfAnnotationProc(objLP)$continuous) <= scaWindowRadPT)
        })
        dfExpectObsCI <- data.frame(
            continuous = rep(dfAnnotationProc(objLP)$continuous,3), 
            trajectory_contour = c(
                qnbinom(p = 1 - sapply(
                    vecObsInBin, function(x) min(x,1) ) / vecObsInBin,
                    size=vecDispParamH1, mu = vecMuParamH1_nonlog),
                qnbinom(p = 1 - sapply(
                    vecObsInBin, function(x) min(x,5) ) / vecObsInBin, 
                    size=vecDispParamH1, mu = vecMuParamH1_nonlog),
                qnbinom(p = 1 - sapply(
                    vecObsInBin, function(x) min(x,10) ) / vecObsInBin, 
                    size=vecDispParamH1, mu = vecMuParamH1_nonlog)),
            ncells = factor(c(rep("1 cell", length(vecMuParamH1_nonlog)),
                              rep("5 cells", length(vecMuParamH1_nonlog)),
                              rep("10 cells", length(vecMuParamH1_nonlog))),
                            levels = c("1 cell", "5 cells", "10 cells")),
            stringsAsFactors = FALSE
        )
        if(boolLogPlot){
            dfExpectObsCI$trajectory_contour <- 
                log(dfExpectObsCI$trajectory_contour)/log(10)
        }
        gplotGene <- gplotGene + geom_line(data = dfExpectObsCI, aes(
            x = continuous, y = trajectory_contour, colour = ncells)) +
            scale_colour_manual(values = cbbPalette, name = "contours")
        scale_alpha_continuous(guide=FALSE)
    } 
    
    gplotGene <- gplotGene + 
        labs(title=paste0(
            strGeneID, "\nlog10 q-value=", 
            round(log(dfResults(objLP)[strGeneID,]$padj,2)/log(10)) )) +
        xlab(paste0("continuous covariate")) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14,face="bold"),
              title=element_text(size=14,face="bold"),
              legend.text=element_text(size=14))
    if(boolLogPlot) {
        gplotGene <- gplotGene + ylab(paste0("log10(counts+10)"))
    } else { 
        gplotGene <- gplotGene + ylab(paste0("counts"))
    }
    
    return(gplotGene)
}

#' Plot density of cells in continuous covariate
#' 
#' Uses kernel density estimate.
#' 
#' @param objLP (LineagePulseObject) LineagePulseObject to base plot on.
#' 
#' @return gplotKDE (ggplot object)
#' ggplot2 kernel density estimator plot. 
#' 
#' @examples
#' lsSimulatedData <- simulateContinuousDataSet(
#'     scaNCells = 100,
#'     scaNConst = 10,
#'     scaNLin = 10,
#'     scaNImp = 10,
#'     scaMumax = 100,
#'     scaSDMuAmplitude = 3,
#'     vecNormConstExternal=NULL,
#'     vecDispExternal=rep(20, 30),
#'     vecGeneWiseDropoutRates = rep(0.1, 30))
#' objLP <- runLineagePulse(
#'     counts = lsSimulatedData$counts,
#'     dfAnnotation = lsSimulatedData$annot,
#'     strMuModel = "impulse")
#' gplotCellDensity <- plotCellDensity(objLP = objLP)
#' #print(gplotCellDensity)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
plotCellDensity <- function(objLP){
    
    gplotKDE <- ggplot() + geom_density(data = dfAnnotationProc(objLP), aes(
        x = continuous))
    
    return(gplotKDE)
}