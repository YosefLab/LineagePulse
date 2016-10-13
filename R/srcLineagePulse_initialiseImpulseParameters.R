#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++     Estimate impulse parameters for initialisation    ++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Estimate impulse model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model.
#' 
#' @seealso Called by \code{fitImpulse_gene}.
#' 
#' @param vecTimepoints: (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecCounts: (count vector number of samples) Count data.
#' #' @param vecDropoutRate: (probability vector number of samples) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel for each gene and cell.
#' @param vecProbNB: (probability vector number of samples) [Default NULL]
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecTimepointAssign: (numeric vector number samples) 
#'    Timepoints assigned to samples.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#'    
#' @return vecParamGuessPeak: (numeric vector number of impulse
#'    model parameters) Impulse model parameter initialisation 
#'    corresponding to a peak.
#' @export

initialiseImpulseParametes <- function(vecTimepoints, 
  vecCounts,
  vecDropoutRate=NULL,
  vecProbNB=NULL,
  strSCMode="clustered",
  vecTimepointAssign, 
  vecNormConst,
  strMode){
  
  # Compute general statistics for initialisation:
  # Expression means by timepoint
  if(strMode=="batch" | strMode=="longitudinal"){
    vecExpressionMeans <- sapply(vecTimepoints,
      function(tp){mean(vecCounts[vecTimepointAssign==tp], na.rm=TRUE)})
  } else if(strMode=="singlecell"){
    if(strSCMode=="clustered"){
      vecExpressionMeans <- unlist(sapply( vecTimepoints, function(tp){
        sum((vecCounts/vecNormConst*vecProbNB)[vecTimepointAssign==tp], na.rm=TRUE)/
          sum(vecProbNB[vecTimepointAssign==tp], na.rm=TRUE)
      } ))
      # Catch exception: sum(vecProbNB[vecTimepointAssign==tp])==0
      vecExpressionMeans[is.na(vecExpressionMeans)] <- 0
    } else if(strSCMode=="continuous"){
      vecTimepointAssignSort <- sort(vecTimepointAssign, index.return=TRUE)
      vecCountsSort <- vecCounts[vecTimepointAssignSort$ix]
      vecProbNBSort <- vecProbNB[vecTimepointAssignSort$ix]
      scaWindows <- 10
      scaCellsPerClus <- round(length(vecCountsSort)/scaWindows)
      vecExpressionMeans <- array(NA, scaWindows)
      vecTimepoints <- array(NA, scaWindows)
      scaidxNew <- 0
      for(k in seq(1,scaWindows)){
        # Define clusters as groups of cells of uniform size
        scaidxLast <- scaidxNew + 1
        scaidxNew <- scaidxLast + scaCellsPerClus - 1
        # Pick up remaining cells in last cluster
        if(k==scaWindows){scaidxNew=length(vecCountsSort)}
        vecidxK <- seq(scaidxLast, scaidxNew)
        # Infer negative binomial mean parameter
        # Mean estimation here has to be replaced if size factors
        # are used.
        vecExpressionMeans[k] <- sum((vecCountsSort*vecProbNBSort)[vecidxK], na.rm=TRUE)/
          sum(vecProbNBSort[vecidxK], na.rm=TRUE)
        vecTimepoints[k] <- mean((vecTimepointAssignSort$x)[vecidxK], na.rm=TRUE)
      }
      # Catch exception: sum(vecProbNB[vecTimepointAssign==tp])==0
      vecExpressionMeans[is.na(vecExpressionMeans)] <- 0
    } else {
      stop(paste0("ERROR: Unrecognised strSCMode in fitImpulse(): ",strSCMode))
    }
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  nTimepts <- length(vecTimepoints)
  scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  # +1 to push indicices up from middle stretch to entire window (first is omitted here)
  indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  # Gradients between neighbouring points
  vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
    (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
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
  vecParamGuessPeak <- c(1,log(vecExpressionMeans[1]+1),
    log(scaMaxMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2)
  
  # Compute valley initialisation
  # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(
    min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), 
    vecGradients[1:(indMinMiddleMean-1)])
  indUpperInflexionPoint <- indMinMiddleMean - 1 + match(
    max(vecGradients[indMinMiddleMean:(nTimepts-1)], na.rm=TRUE), 
    vecGradients[indMinMiddleMean:(nTimepts-1)])
  vecParamGuessValley <- c(1,log(vecExpressionMeans[1]+1),
    log(scaMinMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2 )
  
  lsParamGuesses <- list(peak=vecParamGuessPeak, valley=vecParamGuessValley)
  return(lsParamGuesses)
}