# implement column correlation for sparseMatrix
compColCorSparse <- function(matSparse){
    scaNRow <- nrow(matSparse)
    vecColMeans <- colMeans(matSparse)
    vecColSums <- colSums(matSparse)
    
    # cov_j,k=sum_i^n( (x_ij-mu_j)*(x_ik-mu_k) )
    matCov <- (# sum_i^n x_ij*x_ik
        crossprod(matSparse) -
            # sum_i^n(x_ij*mu_k+x_ik*mu_j)=
            # mu_k*sum_i^n(x_ij)+mu_j*sum_i^n(x_ik)=
            # mu_k*sum_j+mu_j*sum_k=
            # mu_k*mu_j*n+mu_j*mu_k*n=
            # 2*mu_k*mu_j*n=
            # 2*mu_j*sum_k
            2*(vecColMeans %o% vecColSums) + 
            # sum_i^n mu_j*mu_k
            scaNRow*vecColMeans %o% vecColMeans)/ 
        (scaNRow-1)
    
    vecSD <- sqrt(diag(matCov)) 
    matCor <- matCov/(vecSD %o% vecSD)
    matCor <- as.matrix(matCor[,,sparse=FALSE]) # turn into standard matrix
    return(matCor)
}

initialiseCentroidsFromCells <- function(
    matCounts,
    lsvecFixedCentrByPop=NULL,
    vecAssignPop=NULL,
    scaN){
    
    scaNCells <- dim(matCounts)[2]
    # Centroid array
    vecidxCentroidCells <- array(NA, scaN)
    # Work on correlation-based distance
    matDist <- 1-compColCorSparse(matSparse=matCounts)
    
    # Assign centroids for populations set a priori
    scaNCentrAssign <- 0
    for(pop in lsvecFixedCentrByPop){
        scaNCentrPop <- length(pop)
        # Chose cell with minimum cumulative distance to other cells
        # -> cells which lies in centre of population.
        vecidxCurrentCells <- which(vecAssignPop==names(pop))
        vecidxCumulDist <- order(apply(
            matDist[vecidxCurrentCells, vecidxCurrentCells], 1, 
            function(cell) sum(cell, na.rm=TRUE) ),
            decreasing=FALSE)
        vecidxCentroidCells[(scaNCentrAssign+1):(scaNCentrAssign+scaNCentrPop)] <- 
            vecidxCumulDist[1:scaNCentrPop]
        scaNCentrAssign <- scaNCentrAssign+scaNCentrPop
    }
    
    # Greedy: Iteratively select cell with high cumulative distance
    # Not maximal to avoid outlier cells as centroids.
    # to all other previously selected cells.
    # Scales linear in scaN.
    for(centroid in which(is.na(vecidxCentroidCells))){
        vecidxCellsCovered <- c(vecidxCentroidCells[!is.na(vecidxCentroidCells)], 
                                which(!is.na(vecAssignPop)))
        vecidxCellsLeft <- setdiff(seq(1,scaNCells), vecidxCellsCovered)
        if(length(vecidxCellsCovered)>0){
            vecidxCumulDistSort <- order(apply(
                matDist[vecidxCellsLeft,vecidxCellsCovered,drop=FALSE], 1, 
                function(cell) sum(cell, na.rm=TRUE)))
            vecidxCentroidCells[centroid] <- 
                vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
        } else {
            vecidxCumulDistSort <- order(apply(
                matDist, 1, function(cell) sum(cell, na.rm=TRUE)))
            vecidxCentroidCells[centroid] <- 
                vecidxCumulDistSort[round(length(vecidxCumulDistSort)/2)]
        }
    }
    
    return(vecidxCentroidCells)
}