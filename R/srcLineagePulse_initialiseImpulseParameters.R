#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++     Estimate impulse parameters for initialisation    +++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Estimate impulse model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model and are based natural cubic spline fits 
#' in log space with a gaussian noise model.
#' 
#' @seealso Called by impulse model fitting wrappers:
#' \code{fitImpulseZINB}.
#' 
#' @param vecCounts (count vector number of samples) Count data.
#' @param lsMuModelGlobal (list)
#' Object containing meta-data of gene-wise mean parameter models.
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
    lsMuModelGlobal){
    
    # Estimate via natural cubic spline fit and initialise based on that fit
    vecContinuousCovarSpline <- ns(x = lsMuModelGlobal$vecContinuousCovar, 
                              df = 4, intercept = TRUE)
    vecSplineFit <- exp(lm( log(vecCounts/lsMuModelGlobal$vecNormConst+1) ~ 
                                0+vecContinuousCovarSpline)$fitted.values)
    
    # Compute peak initialisation
    # Beta: Has to be negative, Theta1: Low, Theta2: High, Theta3: Low
    # t1: Around first observed inflexion point, t2: 
    # Around second observed inflexion point
    # beta: set so that curve is roughly linear (ie very smooth):
    # reach half of difference between expr at start and inflexion point
    # half way between time of start and inflexsion point.
    # h0+1/4*(h0-h1)=h0+(h0-h1)/(1+exp(-b1*((tmin+(t1-tmin)/2)-t1)))
    # => 1/4 = 1/(1+exp(-b1*(tmin-t1)/2))
    # => 4 = 1+exp(-b1*(tmin-t1)/2)
    # => log(3) = -b1*(tmin-t1)/2
    # => b1 = -(2*log(3))/(tmin-t1)
    # => b1 = (2*log(3))/(t1-tmin)
    # and similar for t2
    scaT1Peak <- mean(lsMuModelGlobal$vecContinuousCovar[
        c(1, which.max(vecSplineFit[seq(2,length(vecSplineFit))])+1)])
    scaT2Peak <- mean(lsMuModelGlobal$vecContinuousCovar[
        c(which.max(vecSplineFit[seq(1,length(vecSplineFit)-1)]), 
          length(vecSplineFit))])
    # Catch exception that continuous covariate values are not unique and 
    # inflexion point falls on boundary of space:
    # Set inflexion points on first unique value next to boundary.
    if(scaT1Peak == lsMuModelGlobal$vecContinuousCovar[1]) {
        scaT1Peak <- min(lsMuModelGlobal$vecContinuousCovar[
            lsMuModelGlobal$vecContinuousCovar != scaT1Peak])
    }
    if(scaT2Peak == lsMuModelGlobal$vecContinuousCovar[
        length(lsMuModelGlobal$vecContinuousCovar)]) {
        scaT1Peak <- max(lsMuModelGlobal$vecContinuousCovar[
            lsMuModelGlobal$vecContinuousCovar != scaT2Peak])
    }
    vecParamGuessPeak <- c(
        2*log(3)/(scaT1Peak-min(lsMuModelGlobal$vecContinuousCovar)), # beta 1
        2*log(3)/(max(lsMuModelGlobal$vecContinuousCovar)-scaT2Peak), # beta 2
        log(vecSplineFit[1]+1), # h0
        log(max(vecSplineFit)+1), # h1
        log(vecSplineFit[length(vecSplineFit)]+1), # h2
        scaT1Peak, # t1
        scaT2Peak ) # t2
    
    # Compute valley initialisation
    # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
    # t1: Around first observed inflexion point, t2: 
    # Around second observed inflexion point
    scaT1Valley <- mean(lsMuModelGlobal$vecContinuousCovar[
            c(1, which.min(vecSplineFit[seq(2,length(vecSplineFit))])+1)])
    scaT2Valley <- mean(lsMuModelGlobal$vecContinuousCovar[
            c(which.min(vecSplineFit[seq(1,length(vecSplineFit)-1)]), 
              length(vecSplineFit))])
    # Catch exception that continuous covariate values are not unique and 
    # inflexion point falls on boundary of space:
    # Set inflexion points on first unique value next to boundary.
    if(scaT1Valley == lsMuModelGlobal$vecContinuousCovar[1]) {
        scaT1Valley <- min(lsMuModelGlobal$vecContinuousCovar[
            lsMuModelGlobal$vecContinuousCovar != scaT1Valley])
    }
    if(scaT2Valley == lsMuModelGlobal$vecContinuousCovar[
        length(lsMuModelGlobal$vecContinuousCovar)]) {
        scaT2Valley <- max(lsMuModelGlobal$vecContinuousCovar[
            lsMuModelGlobal$vecContinuousCovar != scaT2Valley])
    }
    vecParamGuessValley <- c(
        2*log(3)/(scaT1Valley-min(lsMuModelGlobal$vecContinuousCovar)), # beta 1
        2*log(3)/(max(lsMuModelGlobal$vecContinuousCovar)-scaT2Valley), # beta 2
        log(vecSplineFit[1]+1), # h0
        log(min(vecSplineFit)+1), # h1
        log(vecSplineFit[length(vecSplineFit)]+1), # h2
        scaT1Valley, # t1
        scaT2Valley ) # t2
    
    return( list(peak=vecParamGuessPeak,
                 valley=vecParamGuessValley) )
}