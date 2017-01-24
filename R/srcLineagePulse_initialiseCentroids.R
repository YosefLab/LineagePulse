initialiseCentroidsFromCells <- function(matCounts, 
                                         scaN){
  scaNCells <- dim(matCounts)[2]
  # Work on correlation-based distance
  matDist <- 1-cor(matCounts)
  
  # Greedy: Iteratively select cell with maximal cumulative distance
  # to all other previously selected cells.
  # Scales linear in scaN.
  vecidxCells <- array(NA, scaN)
  vecDistSort <- sort(matDist)
  vecidxCells[c(1,2)] <- which(matDist == vecDistSort[round(length(vecDistSort)/2)], arr.ind = TRUE)[1,]
  if(scaN>2){
    for(i in seq(3,scaN)){
      vecidxCurrentCells <- vecidxCells[!is.na(vecidxCells)]
      vecidxCellsLeft <- setdiff(seq(1,scaNCells), vecidxCurrentCells)
      vecidxCumulDistSort <- sort(apply(matDist[vecidxCellsLeft,vecidxCurrentCells], 1, function(cell) sum(cell, na.rm=TRUE)), index.return=TRUE)$ix
      vecidxCells[i] <- vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
    }
  }
  return(vecidxCells)
}