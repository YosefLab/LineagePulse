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
#' 
#' @param objLP: (LineagePulseObject) LineagePulseObject
#'    base plot on.
#' @param strGeneID: (str) Name of gene, used for title of plot.
#' @param vecReferenceMuParam: (numeric vector length number of cells)
#'    [Default NULL] Reference mean trajectory which can be plotted
#' @param strTitleSuffix: (str) String to be added to title.
#' @param boolColourByDropout: (bool) Whether to colour scatter
#'    plot by posterior of drop-out.
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
    boolH1NormCounts=FALSE ){
    
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
                    x=x, y=counts, colour=dropout_posterior, shape = groups), show.legend=TRUE)
        } else {
            gGenePlot <- ggplot() +
                geom_point(data=dfScatterCounts, aes(x=x, y=counts, shape = groups), show.legend=TRUE)
        }
    }
    
    gGenePlot <- gGenePlot + 
        labs(title=paste0(
            strGeneID, "\nlog10 q-value=", 
            round(log(objLP@dfResults[strGeneID,]$adj.p,2)/log(10)) )) +
        xlab(paste0("x")) +
        scale_colour_gradient(high="red",low="green",limits=c(0, 1)) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14,face="bold"),
              title=element_text(size=14,face="bold"),
              legend.text=element_text(size=14))
    if(boolLogPlot) gGenePlot <- gGenePlot + ylab(paste0("log10 counts"))
    else gGenePlot <- gGenePlot + ylab(paste0("counts"))
    
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
    
    return(gGenePlot)
}