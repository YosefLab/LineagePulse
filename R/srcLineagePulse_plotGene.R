#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++  Plot counts and model for one gene  +++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot counts and model for one gene
#' 
#' Plot counts and model for one gene as a scatter plot with model trajectories:
#' Both as counts or log10 counts vs pseudotime.
#' Dropout rates can be given and are visualised as the colour of the observation
#' points.
#' A reference model trajectory can be added.
#' Contour plots: A contour plot aids the interpretation of the model trajectory
#' in context of the observed data if local overplotting occurs because of a
#' non-uniform distribution of cells in pseudotime (a very common scenario).
#' In such a case, one might mistake cells with high counts from a pseudotime
#' interval with high cellular intensity as a indication of increased expression
#' in this interval. Instead, one can explain such local maxima based on the 
#' simple phenomenom that more samples were taken from the underyling distribution
#' in that pseudotime interval which increases the probability of sampling
#' far in the tail of the distribution. To migitate this visually, we added the
#' option to add "contours": A contour is an alternative to a confidence interval
#' which includes information on the local sampling density. The contour is specific
#' to a number of cells n. It is an expression trajectory in pseudotime which lies
#' at the line above which n cells are expected within a local pseudotime bin given
#' its sampling density and the H1 non-constant expression model. Contours therefore
#' quantify the increase in maximum observed counts in pseudotime intervals due to
#' sampling density, which is not due to the expression model. Note that this effect
#' is corrected for in the model fitting and that the contour plots just aids visual
#' interpretation of the fit to the data!
#' 
#' @param objLP (LineagePulseObject) LineagePulseObject
#'    base plot on.
#' @param strGeneID (str) Name of gene, used for title of plot.
#' @param vecReferenceMuParam (numeric vector length number of cells)
#'    [Default NULL] Reference mean trajectory which can be plotted
#' @param strTitleSuffix (str) String to be added to title.
#' @param boolColourByDropout (bool) Whether to colour scatter
#' plot by posterior of drop-out.
#' @param boolLineageContour (bool) [Default FALSE]
#' Whether to the "lineage contour" lines to the scatter plot.
#' @param bwDensity (bandwith: numeric or string) [Default NULL]
#' Bandwith to be used to kernel density smooting
#' of cell density in pseudotime (used if boolExpectedDistrPanel=TRUE).
#' If not set, defaults to stats:density() default.
#' 
#' @return gGenePlot: (ggplot object)
#'    Model rajectories and scatter plot for given gene.
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
    bwDensity = NULL){
    
    cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "grey", "red", "pink")
    ### 1. Extract data and models
    vecCounts <- objLP@matCountsProc[strGeneID,,sparse=FALSE]
    matCountsGene <- t(as.matrix(vecCounts))
    rownames(matCountsGene) <- strGeneID
    if(boolH1NormCounts){
        vecCountsToPlot <- as.vector(getNormData(
            matCounts=matCountsGene,
            lsMuModel=objLP@lsMuModelH1))
    } else {
        vecCountsToPlot <- vecCounts
    }
    vecMuParamH0 <- decompressMeansByGene(
        vecMuModel=objLP@lsMuModelH0$matMuModel[strGeneID,],
        lsvecBatchModel=NULL,
        lsMuModelGlobal=objLP@lsMuModelH0$lsMuModelGlobal,
        vecInterval=NULL )
    vecMuParamH1 <- decompressMeansByGene(
        vecMuModel=objLP@lsMuModelH1$matMuModel[strGeneID,],
        lsvecBatchModel=NULL,
        lsMuModelGlobal=objLP@lsMuModelH1$lsMuModelGlobal,
        vecInterval=NULL )
    vecDropPosterior <- as.vector(calcPostDrop_Matrix(
        matCounts=matCountsGene,
        lsMuModel=objLP@lsMuModelH1,
        lsDispModel=objLP@lsDispModelH1, 
        lsDropModel=objLP@lsDropModel,
        vecIDs=strGeneID))
    if(boolLogPlot){
        vecCountsToPlot <- log(vecCountsToPlot)/log(10)
        vecMuParamH0 <- log(vecMuParamH0)/log(10)
        vecMuParamH1 <- log(vecMuParamH1)/log(10)
    }
    
    ### 2. Data scatter plot
    # Set drop-out rates as constant for visualistion if not given.
    dfScatterCounts <- data.frame( counts=vecCountsToPlot )
    if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("splines", "impulse")){
        dfScatterCounts$x <- objLP@lsMuModelH1$lsMuModelGlobal$vecPseudotime
        
        if(boolColourByDropout){
            dfScatterCounts$dropout_posterior <- vecDropPosterior
            gGenePlot <- ggplot() +
                geom_point(data=dfScatterCounts, aes(x=x, y=counts, colour=dropout_posterior), show.legend=TRUE)
        } else {
            gGenePlot <- ggplot() +
                geom_point(data=dfScatterCounts, aes(x=x, y=counts), show.legend=TRUE)
        }
        
    } else if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("groups")){
        dfScatterCounts$x <- seq(1, objLP@lsMuModelH1$lsMuModelGlobal$scaNumCells, by=1)
        dfScatterCounts$groups <- objLP@dfAnnotationProc$groups
        
        if(boolColourByDropout){
            dfScatterCounts$dropout_posterior <- vecDropPosterior
            gGenePlot <- ggplot() +
                geom_point(data=dfScatterCounts, aes(
                    x=x, y=counts, colour=dropout_posterior, shape = groups), show.legend=TRUE) +
                scale_colour_gradient(high="red",low="green",limits=c(0, 1)) 
        } else {
            gGenePlot <- ggplot() +
                geom_point(data=dfScatterCounts, aes(
                    x=x, y=counts, shape = groups), show.legend=TRUE)
        }
    }
    
    # Set plotting threshold based on observed data
    scaMaxPlot <- 2 * max(vecCounts)
    
    ### 3. Add models to plot
    vecMuParamH0[vecMuParamH0 > scaMaxPlot] <- NA
    vecMuParamH1[vecMuParamH1 > scaMaxPlot] <- NA
    if(is.null(vecReferenceMuParam)){
        if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("splines", "impulse")){
            dfLineImpulse <- data.frame(
                x=rep(objLP@lsMuModelH1$lsMuModelGlobal$vecPseudotime, 2),
                counts=c(vecMuParamH0,
                         vecMuParamH1),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1))) )
        } else if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("groups")){
            dfLineImpulse <- data.frame(
                x=rep(seq(1,objLP@lsMuModelH1$lsMuModelGlobal$scaNumCells,by=1), 2),
                counts=c(vecMuParamH0,
                         vecMuParamH1),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1))) )
        }
    } else {
        if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("splines", "impulse")){
            dfLineImpulse <- data.frame(
                x=rep(objLP@lsMuModelH1$lsMuModelGlobal$vecPseudotime, 3),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecReferenceMuParam),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("Reference", length(vecReferenceMuParam))) )
        } else if(objLP@lsMuModelH1$lsMuModelGlobal$strMuModel %in% c("groups")){
            dfLineImpulse <- data.frame(
                x=rep(seq(1,objLP@lsMuModelH1$lsMuModelGlobal$scaNumCells,by=1), 3),
                counts=c(vecMuParamH0,
                         vecMuParamH1,
                         vecReferenceMuParam),
                model=c(rep("H0", length(vecMuParamH0)), 
                        rep("H1", length(vecMuParamH1)),
                        rep("Reference", length(vecReferenceMuParam))) )
        }
    }
    gGenePlot <- gGenePlot + geom_line(
        data=dfLineImpulse, 
        aes(x=x, y=counts, linetype=model),
        show.legend=TRUE)
    if(boolH1NormCounts) {
        gGenePlot <- gGenePlot + ylab("normalised counts")
    } else {
        gGenePlot <- gGenePlot + ylab("counts")
    }
    
    if(boolExpectedDistrPanel) {
        # decompress log free
        if(boolLogPlot) {
            vecMuParamH1_nonlog <- 10^(vecMuParamH1)
        } else {
            vecMuParamH1_nonlog <- vecMuParamH1
        }
        vecDispParamH1 <- decompressDispByGene(
            vecDispModel=objLP@lsDispModelH1$matDispModel[strGeneID,],
            lsvecBatchModel=NULL,
            lsDispModelGlobal=objLP@lsDispModelH1$lsDispModelGlobal,
            vecInterval=NULL )
        # copmute density estimate of cells in pseudotime
        if(is.null(bwDensity)) {
            bwDensity <- "nrd0"
        }
        lsDensity <- density(x = objLP@dfAnnotationProc$pseudotime, bw=bwDensity, 
                             n = length(objLP@dfAnnotationProc$pseudotime))
        scaNKDEbin <- 10
        scaWindowRadPT <- (max(objLP@dfAnnotationProc$pseudotime, na.rm = TRUE) - 
                               min(objLP@dfAnnotationProc$pseudotime, na.rm = TRUE))/(2*scaNKDEbin)
        vecObsInBin <- sapply(objLP@dfAnnotationProc$pseudotime, function(pt) {
            sum(abs(pt-objLP@dfAnnotationProc$pseudotime) <= scaWindowRadPT)
        })
        dfExpectObsCI <- data.frame(
            pseudotime = rep(objLP@dfAnnotationProc$pseudotime,3), 
            contour = c(qnbinom(p = 1-sapply(vecObsInBin, function(x) min(x,1) )/vecObsInBin, size=vecDispParamH1, mu = vecMuParamH1_nonlog),
                        qnbinom(p = 1-sapply(vecObsInBin, function(x) min(x,5) )/vecObsInBin, size=vecDispParamH1, mu = vecMuParamH1_nonlog),
                        qnbinom(p = 1-sapply(vecObsInBin, function(x) min(x,10) )/vecObsInBin, size=vecDispParamH1, mu = vecMuParamH1_nonlog)),
            ncells = factor(c(rep("1 cell", length(vecMuParamH1_nonlog)),
                              rep("5 cells", length(vecMuParamH1_nonlog)),
                              rep("10 cells", length(vecMuParamH1_nonlog))),
                            levels = c("1 cell", "5 cells", "10 cells")),
            stringsAsFactors = FALSE
        )
        if(boolLogPlot){
            dfExpectObsCI$contour <- log(dfExpectObsCI$contour)/log(10)
        }
        gGenePlot <- gGenePlot + geom_line(data = dfExpectObsCI, aes(
            x = pseudotime, y = contour, colour = ncells)) +
            scale_colour_manual(values = cbbPalette, name = "contours")
        scale_alpha_continuous(guide=FALSE)
    } 
    
    gGenePlot <- gGenePlot + 
        labs(title=paste0(
            strGeneID, "\nlog10 q-value=", 
            round(log(objLP@dfResults[strGeneID,]$adj.p,2)/log(10)) )) +
        xlab(paste0("x")) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14,face="bold"),
              title=element_text(size=14,face="bold"),
              legend.text=element_text(size=14))
    if(boolLogPlot) {
        gGenePlot <- gGenePlot + ylab(paste0("log10 counts"))
    } else { 
        gGenePlot <- gGenePlot + ylab(paste0("counts"))
    }
    
    return(gGenePlot)
}

#' Plot density of cells in pseudotime
#' 
#' @param objLP (LineagePulseObject) LineagePulseObject to base plot on.
#'    
#' @return gGenePlot: (ggplot object)
#' ggplot2 kernel density estimator plot. 
#'    
#' @author David Sebastian Fischer
#' 
#' @export
plotCellDensity <- function(
    objLP,
    bwDensity = NULL){
    
    gplotKSD <- ggplot() + geom_density(data = objLP@dfAnnotationProc, aes(
        x = pseudotime))
    
    return(gplotKSD)
}