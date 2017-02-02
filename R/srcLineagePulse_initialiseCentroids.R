initialiseCentroidsFromCells <- function(matCounts,
                                         vecFixedCells=NULL,
                                         scaN){
  
  scaNCells <- dim(matCounts)[2]
  # Centroid array
  vecidxCentroidCells <- array(NA, scaN)
  # Work on correlation-based distance
  matDist <- 1-cor(matCounts)
  
  # Assign centroids for populations set a priori
  vecPopulationsSet <- unique(vecFixedCells[!is.na(vecFixedCells)])
  for(pop in vecPopulationsSet){
    # Chose cell with minimum cumulative distance to other cells
    # -> cells which lies in centre of population.
    vecidxCurrentCells <- which(vecFixedCells==pop)
    vecCumulDist <- apply(matDist[vecidxCurrentCells,vecidxCurrentCells], 1, function(cell) sum(cell, na.rm=TRUE) )
    vecidxCentroidCells[pop] <- which.min(vecCumulDist)
  }
  
  # Greedy: Iteratively select cell with high cumulative distance
  # Not maximal to avoid outlier cells as centroids.
  # to all other previously selected cells.
  # Scales linear in scaN.
  for(centroid in which(is.na(vecidxCentroidCells))){
    vecidxCellsCovered <- c(vecidxCentroidCells[!is.na(vecidxCentroidCells)], which(!is.na(vecFixedCells)))
    vecidxCellsLeft <- setdiff(seq(1,scaNCells), vecidxCellsCovered)
    vecidxCumulDistSort <- sort(apply(matDist[vecidxCellsLeft,vecidxCellsCovered], 1, function(cell) sum(cell, na.rm=TRUE)), index.return=TRUE)$ix
    vecidxCentroidCells[centroid] <- vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
  }
  
  return(vecidxCentroidCells)
}