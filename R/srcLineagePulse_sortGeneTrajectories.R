#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++    Cluster expression mean trajectories  +++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cluster expression mean trajectories
#' 
#' Sorts inferred gene trajectories by peak time in pseudotime. Optional: Can
#' create a heatmap of the gene trajectories sorted according to peak time.
#' The heatmap is based on z-scores.
#' 
#' @seealso Called by \code{LineagePulse_main} or by user.
#'
#' @param vecIDs (vector of strings)
#' Names of genes to cluster.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param dirHeatmap (str directory) [Default NULL]
#' Directory to which heatmap is saved to. Not 
#' created if NULL. Need to have lsMuModel$lsMuModelGlobal$vecPseudotime
#' defined.
#' 
#' @return vecSortedGenes (string vector number of IDs)
#' IDs sorted by peak time in pseudotime.
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
#' #sortGeneTrajectories(
#' #    vecIDs = rownames(objLP$objLP@lsMuModelH1$matMuModel),
#' #    lsMuModel = objLP@lsMuModelH1,
#' #    dirHeatmap = "your_target_directory"))
#' 
#' @author David Sebastian Fischer
#' 
#' @export
sortGeneTrajectories <- function(
    vecIDs,
    lsMuModel,
    dirHeatmap){
    
    # Check IDs are in model matrix
    if(any(!(vecIDs %in% rownames(lsMuModel$matMuModel)))){
        message(paste0("ERROR: Some IDs in vecIDs are not given",
                     " in model lsMuModel$matMuModel."))
        return(NULL)
    }
    
    # (I) Decompress mean parameters in a memory saving way:
    # Round to integer counts.
    matMuParam <- do.call(rbind, bplapply(vecIDs, function(id){
        round( decompressMeansByGene(
            vecMuModel=lsMuModel$matMuModel[id,],
            lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
            vecInterval=NULL ) )
    }))
    rownames(matMuParam) <- vecIDs
    
    # (II) Sort genes based on peak time.
    # Find peak time (cell index) of each gene.
    vecindPeak <- apply(matMuParam, 1, function(gene){
        (sort(gene, index.return=TRUE)$ix)[1]
    })
    # Sort genes based on peak time
    vecSortedGenes <- vecIDs[sort(vecindPeak, index.return=TRUE)$ix]
    
    # (III) Create heatmap of sorted genes.
    if(!is.null(dirHeatmap)){
        # Row normalise expression values
        matMuNorm <- do.call(rbind, lapply(vecSortedGenes, function(gene){
            scaSD <- sd(matMuParam[gene,], na.rm=TRUE)
            if(scaSD==0){ scaSD <- 1 }
            (matMuParam[gene,]-mean(matMuParam[gene,], na.rm=TRUE))/scaSD
        }))
        # Set column names: Hack tick labeling of heatmap.2: Only
        # shows columns as lables which are not names NA
        scaCells <- length(lsMuModel$lsMuModelGlobal$vecPseudotime)
        vecTicks <- array(NA, scaCells)
        scaDistBetweenCells <- round(scaCells/10)
        vecindTicks <- sapply(seq(0,9), function(i) 1+i*scaDistBetweenCells)
        vecindTicks[10] <- scaCells
        vecTicks[vecindTicks] <- round(
            lsMuModel$lsMuModelGlobal$vecPseudotime[vecindTicks], 0)
        colnames(matMuNorm) <- vecTicks
        # Set upper bounds of z-scores to visualise
        scaMuNormLowBound <- max(min(matMuNorm, na.rm=TRUE), -5)
        scaMuNormUpperBound <- min(max(matMuNorm, na.rm=TRUE), 5)
        scaBreaks <- 50
        breaks <- seq(scaMuNormLowBound,scaMuNormUpperBound,
                      by=(scaMuNormUpperBound - scaMuNormLowBound) /
                          (scaBreaks-1))
        hm.colors <- colorpanel( length(breaks)-1, "red", "blue" )
        
        # Plot genes sorted by peak time
        pdf(paste0(dirHeatmap, "/LineagePulse_GenesSortedByPeakTime.pdf"))
        heatmap.2(matMuNorm, 
                  dendrogram="none", Rowv=FALSE, Colv=FALSE, 
                  xlab = "pseudotime", ylab =  "genes",
                  labRow=NA,# Supress gene names
                  trace="none",
                  density.info="none",
                  lmat=rbind( c(3,4), c(2,1) ), lhei=c(1,3), lwid=c(1,3),
                  key.title = "", key.xlab = "z-score", key.ylab = NULL,
                  symkey=FALSE, 
                  breaks=breaks, 
                  col=hm.colors, 
                  scale="none")
        dev.off()
        graphics.off()
        
        # Plot clusters by gene
        graphics.off()
        pdf(paste0(dirHeatmap, 
                   "/LineagePulse_GenesClusteredByTrajectory.pdf"))
        heatmap.2(matMuNorm, 
                  dendrogram="row", Rowv=TRUE, Colv=FALSE, 
                  xlab = "pseudotime", ylab =  "genes",
                  labRow=NA,# Supress gene names
                  trace="none",
                  density.info="none",
                  lmat=rbind( c(3,4), c(2,1) ), lhei=c(1,3), lwid=c(1,3),
                  key.title = "", key.xlab = "z-score", key.ylab = NULL,
                  symkey=FALSE, 
                  breaks=breaks, 
                  col=hm.colors, 
                  scale="none")
        dev.off()
        graphics.off()
    }
    
    return(vecSortedGenes)
}