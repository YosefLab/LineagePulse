initialiseCentroidsFromCells <- function(matCounts,
                                         lsvecFixedCentrByPop=NULL,
                                         vecAssignPop=NULL,
                                         scaN){
  
  scaNCells <- dim(matCounts)[2]
  # Centroid array
  vecidxCentroidCells <- array(NA, scaN)
  # Work on correlation-based distance
  matDist <- 1-cor(matCounts)
  
  # Assign centroids for populations set a priori
  scaNCentrAssign <- 0
  for(pop in lsvecFixedCentrByPop){
    scaNCentrPop <- length(pop)
    # Chose cell with minimum cumulative distance to other cells
    # -> cells which lies in centre of population.
    vecidxCurrentCells <- which(vecAssignPop==names(pop))
    vecidxCumulDist <- sort(apply(matDist[vecidxCurrentCells,vecidxCurrentCells], 1, 
                                   function(cell) sum(cell, na.rm=TRUE) ),
                             decreasing=FALSE, index.return=TRUE)$ix
    vecidxCentroidCells[(scaNCentrAssign+1):(scaNCentrAssign+scaNCentrPop)] <- 
      vecidxCumulDist[1:scaNCentrPop]
    scaNCentrAssign <- scaNCentrAssign+scaNCentrPop
  }
  
  # Greedy: Iteratively select cell with high cumulative distance
  # Not maximal to avoid outlier cells as centroids.
  # to all other previously selected cells.
  # Scales linear in scaN.
  for(centroid in which(is.na(vecidxCentroidCells))){
    vecidxCellsCovered <- c(vecidxCentroidCells[!is.na(vecidxCentroidCells)], which(!is.na(vecAssignPop)))
    vecidxCellsLeft <- setdiff(seq(1,scaNCells), vecidxCellsCovered)
    if(length(vecidxCellsCovered)>0){
      vecidxCumulDistSort <- sort(apply(matDist[vecidxCellsLeft,vecidxCellsCovered,drop=FALSE], 1, function(cell) sum(cell, na.rm=TRUE)), index.return=TRUE)$ix
      vecidxCentroidCells[centroid] <- vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
    } else {
      vecidxCumulDistSort <- sort(apply(matDist, 1, function(cell) sum(cell, na.rm=TRUE)), index.return=TRUE)$ix
      vecidxCentroidCells[centroid] <- vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
    }
  }
  
  return(vecidxCentroidCells)
}