#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++    Compute value of impulse model   +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of impulse function given parameters.
#' 
#' Compute value of impulse function given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
#' 
#' @seealso Compiled version: \link{evalImpulseModel_comp}
#' 
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#'    \{beta, h0, h1, h2, t1, t2\}
#'    Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#'    Time points to be evaluated.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#'     Model values for given time points.
#'     
#' @author David Sebastian Fischer
#' 
#' @export
evalImpulseModel <- function(vecImpulseParam,
                             vecTimepoints){
  
  # beta1 is vecImpulseParam[1]
  # beta2 is vecImpulseParam[2]
  # h0 is vecImpulseParam[3]
  # h1 is vecImpulseParam[4]
  # h2 is vecImpulseParam[5]
  # t1 is vecImpulseParam[6]
  # t2 is vecImpulseParam[7]
  
  vecImpulseValue <- sapply(vecTimepoints, function(t){
    (1/vecImpulseParam[4]) * 
      (vecImpulseParam[3] + (vecImpulseParam[4]-vecImpulseParam[3])*
         (1/(1+exp(-vecImpulseParam[1]*(t-vecImpulseParam[6]))))) *
      (vecImpulseParam[5] + (vecImpulseParam[4]-vecImpulseParam[5])*
         (1/(1+exp(vecImpulseParam[2]*(t-vecImpulseParam[7])))))
  })
  
  # Catch lower bound on mu space
  #vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}

#' Compiled function: evalImpulseModel
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalImpulseModel}.
#' 
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#'    \{beta, h0, h1, h2, t1, t2\}
#'    Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#'    Time points to be evaluated.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#'     Model values for given time points.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
evalImpulseModel_comp <- cmpfun(evalImpulseModel)
