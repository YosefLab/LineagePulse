#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++     Estimate impulse parameters for initialisation    ++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Estimate impulse model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model and are based on proximity binned data.
#' 
#' @seealso Called by impulse model fitting wrappers:
#' \code{fitImpulseZINB}.
#' 
#' @param vecCounts: (count vector number of samples) Count data.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
#' @param vecMu: (vector number of cells) Negative binomial
#' mean parameter estimates.
#' @param vecDisp: (vector number of cells) Negative binomial
#' dispersion parameter estimates.
#' @param vecDrop: (vector number of cells) Dropout parameter
#' estimates. 
#' 
#' @return list (length 2)
#' \itemize{
#'     \item peak: (numeric vector number of impulse
#' model parameters) Impulse model parameter initialisation 
#' corresponding to a peak.
#'     \item valley: (numeric vector number of impulse
#' model parameters) Impulse model parameter initialisation 
#' corresponding to a valley.
#' }
#' 
#' @author David Sebastian Fischer
initialiseImpulseParameters <- function(
    vecCounts,
    lsMuModelGlobal,
    vecMu,
    vecDisp,
    vecDrop){
    
    # Compute posterior of drop-out
    vecPostNB <- 1 - calcPostDrop_Vector(
        vecMu=vecMu,
        vecDisp=vecDisp,
        vecDrop=vecDrop,
        vecboolZero= !is.na(vecCounts) & vecCounts==0,
        vecboolNotZero= !is.na(vecCounts) & vecCounts>0 )
    
    # Observations are pooled to give rough estimates of expression
    # levels in local environment. The pooling are either the clusters
    # of cells (if the mean model is clusters) or cells grouped
    # by proximity in pseudotime.
    if(lsMuModelGlobal$strMuModel=="groups"){
        vecidxGroups <- lsMuModelGlobal$vecindClusterAssign
    } else {
        scaGroups <- 10
        scaCellsPerGroup <- round(length(vecCounts)/scaGroups)
        vecidxGroups <- array(NA, length(vecCounts))
        scaidxNew <- 0
        for(k in seq(1,scaGroups)){
            # Define clusters as groups of cells of uniform size
            scaidxLast <- scaidxNew + 1
            scaidxNew <- scaidxLast + scaCellsPerGroup - 1
            # Pick up remaining cells in last cluster
            if(k==scaGroups){scaidxNew=length(vecCounts)}
            vecidxGroups[seq(scaidxLast, scaidxNew)] <- k
        }
    }
    
    # Compute characteristics of the groups: Expression level
    # and time coordinate.
    vecGroupExprMean <- array(NA, scaGroups)
    vecGroupTimeCoord <- array(NA, scaGroups)
    for(k in unique(vecidxGroups)){
        vecidxK <- match(k, vecidxGroups) 
        vecGroupExprMean[k] <- sum(vecCounts[vecidxK] / 
                                       lsMuModelGlobal$vecNormConst[vecidxK]*
                                       vecPostNB[vecidxK], na.rm=TRUE)/
            sum(vecPostNB[vecidxK], na.rm=TRUE)
        vecGroupTimeCoord[k] <- mean(lsMuModelGlobal$vecPseudotime[vecidxK], na.rm=TRUE)  
    }
    # Catch exception: sum(vecPostNB[vecTimepointAssign==tp])==0
    vecGroupExprMean[is.na(vecGroupExprMean)] <- 0
    
    # Comute characteristics of the expression trajectory in groups.
    nTimepts <- length(vecGroupTimeCoord)
    scaMaxMiddleMean <- max(vecGroupExprMean[2:(nTimepts-1)], na.rm=TRUE)
    scaMinMiddleMean <- min(vecGroupExprMean[2:(nTimepts-1)], na.rm=TRUE)
    # +1 to push indicices up from middle stretch to entire window (first is omitted here)
    indMaxMiddleMean <- match(scaMaxMiddleMean,vecGroupExprMean[2:(nTimepts-1)]) + 1
    indMinMiddleMean <- match(scaMinMiddleMean,vecGroupExprMean[2:(nTimepts-1)]) + 1
    # Gradients between neighbouring points
    vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
        (vecGroupExprMean[x+1]-vecGroupExprMean[x]) / 
            (vecGroupTimeCoord[x+1]-vecGroupTimeCoord[x])}) )
    vecGradients[is.na(vecGradients) | !is.finite(vecGradients)] <- 0
    
    # Compute peak initialisation
    # Beta: Has to be negative, Theta1: Low, Theta2: High, Theta3: Low
    # t1: Around first observed inflexion point, t2: Around second observed inflexion point
    indLowerInflexionPoint <- match(
        max(vecGradients[1:(indMaxMiddleMean-1)], na.rm=TRUE), 
        vecGradients[1:(indMaxMiddleMean-1)])
    indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(
        min(vecGradients[indMaxMiddleMean:length(vecGradients)], na.rm=TRUE), 
        vecGradients[indMaxMiddleMean:length(vecGradients)])
    vecParamGuessPeak <- c(
        1,
        1,
        log(vecGroupExprMean[1]+1),
        log(scaMaxMiddleMean+1),
        log(vecGroupExprMean[nTimepts]+1),
        (vecGroupTimeCoord[indLowerInflexionPoint] + 
             vecGroupTimeCoord[indLowerInflexionPoint+1])/2,
        (vecGroupTimeCoord[indUpperInflexionPoint] + 
             vecGroupTimeCoord[indUpperInflexionPoint+1])/2)
    
    # Compute valley initialisation
    # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
    # t1: Around first observed inflexion point, t2: Around second observed inflexion point
    indLowerInflexionPoint <- match(
        min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), 
        vecGradients[1:(indMinMiddleMean-1)])
    indUpperInflexionPoint <- indMinMiddleMean - 1 + match(
        max(vecGradients[indMinMiddleMean:(nTimepts-1)], na.rm=TRUE), 
        vecGradients[indMinMiddleMean:(nTimepts-1)])
    vecParamGuessValley <- c(
        1,
        1,
        log(vecGroupExprMean[1]+1),
        log(scaMinMiddleMean+1),
        log(vecGroupExprMean[nTimepts]+1),
        (vecGroupTimeCoord[indLowerInflexionPoint] + 
             vecGroupTimeCoord[indLowerInflexionPoint+1])/2,
        (vecGroupTimeCoord[indUpperInflexionPoint] + 
             vecGroupTimeCoord[indUpperInflexionPoint+1])/2 )
    
    return( list(peak=vecParamGuessPeak,
                 valley=vecParamGuessValley) )
}