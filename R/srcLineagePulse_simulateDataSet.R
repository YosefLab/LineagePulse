#' Simulate a data set for LinagePulse

#' Simulates a data set with genes with constant and impulse
#' expression traces. Expression strength and variation in impulse
#' like traces are parameterised and random. All temporary files
#' are saved into dirOutSimulation and only the objects necessary
#' for running LineagePulse (the count matrix and the pseudotime
#' vector are returned). The remaining objects representing hidden
#' parameters can be used to evaluate parameter estimates. Cells
#' are distributed uniformly in pseudotime.
#' 
#' @seealso Called by separately by user.
#' 
#' @param scaNCells: (scalar) Number of cells in data set.
#' @param scaNConst: (scalar) Number of constant genes in data set.
#' @param scaNImp: (scalar) Number of impulse distributed genes in data set.
#' @param scaPTmax: (scalar) [Default 100]
#'    Maximum pseudotime coordinate of all cells. This doesnt
#'    really matter as cells are uniformly distributed in the 
#'    current implementation.
#' @param scaMumax: (scalar) [Default 1000]
#'    Maximum expression mean parameter to be used.
#' @param scaSDImpulseAmplitude: (scalar) [Default 1]
#'    Standard deviation of normal distribution form which the 
#'    amplitude change within an impulse trace is drawn.
#' @param vecSizeFactorsExternal: (numeric vector number of cells)
#'    [Default NULL]
#'    Size factors for data set. Size factors are set to 1 if this is
#'    not specified (NULL).
#' @param matDropoutModelExternal: (numeric matrix cells x 2) 
#'    [Default NULL] External drop-out model, has to have
#'    one row for each simulated cell.
#' @param dirOutSimulation: (str directory)
#'    Directory to which simulated parameter objects are 
#'    saved to.
#' 
#' @return vecPT: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#' @return matSampledCountsObserved: (matrix genes x cells)
#'    Sampled count data of all cells after drop-out.
#' 
#' @export

simulateDataSet <- function(scaNCells,
  scaNConst,
  scaNImp,
  scaPTmax=100,
  scaMumax=1000,
  scaSDImpulseAmplitude=1,
  vecSizeFactorsExternal=NULL,
  matDropoutModelExternal=NULL,
  dirOutSimulation){
  
  ####
  # Internal functions
  # Evalute impulse model at time points
  evalImpulse <- function(t,beta,t1,t2,h0,h1,h2){
    return(1/h1* (h0+(h1-h0)*1/(1+exp(-beta*(t-t1))))*
        (h2+(h1-h2)*1/(1+exp(beta*(t-t2)))))
  }
  # Evaluate logisitic model at mean parameter values
  evalLogistic <- function(mu, a1, a2){
    return(1/(1+exp(-(a1+a2*mu))))
  }
  
  ####
  # Simulate data
  # 1. Distribute cells across pseudotime
  vecPT <- seq(0, scaPTmax, by=scaPTmax/(scaNCells-1))
  names(vecPT) <- paste0("_",seq(1,scaNCells))
  
  # 2. Create hidden data set
  vecConstIDs <- paste0(rep("_",scaNConst),c(1:scaNConst))
  vecImpulseIDs <- paste0(rep("_",scaNImp),c((scaNConst+1):(scaNConst+scaNImp)))
  # a. Draw means from uniform (first half of genes): one mean per gene
  vecMuConstHidden <- runif(scaNConst)*scaMumax
  matMuConstHidden <- matrix(vecMuConstHidden,
    nrow=scaNConst,
    ncol=scaNCells,
    byrow=FALSE )
  rownames(matMuConstHidden) <- vecConstIDs
  colnames(matMuConstHidden) <- names(vecPT)
  
  # b. Draw means from impulse model
  beta <- runif(scaNImp)*2+0.5
  ta <- runif(scaNImp)*scaPTmax
  tb <- runif(scaNImp)*scaPTmax
  t1 <- ta
  t1[tb < ta] <- tb[tb < ta]
  t2 <- tb
  t2[tb < ta] <- ta[tb < ta]
  h0 <- runif(scaNImp)*scaMumax
  h1 <- h0*abs(rnorm(n=scaNImp, mean=0,sd=scaSDImpulseAmplitude))
  h2 <- h0*abs(rnorm(n=scaNImp, mean=0,sd=scaSDImpulseAmplitude))
  h0[h0<0.00001] <-0.00001
  h1[h1<0.00001] <-0.00001
  h2[h2<0.00001] <-0.00001
  lsMuImpulseHidden <- lapply(seq(1,scaNImp), function(gene){
    evalImpulse(t=vecPT,
      beta=beta[gene],
      t1=t1[gene],
      t2=t2[gene],
      h0=h0[gene],
      h1=h1[gene],
      h2=h2[gene])
  })
  matImpulseModelHidden <- cbind(beta, h0, h1, h2, t1, t2)
  matMuImpulseHidden <- do.call(rbind, lsMuImpulseHidden)
  rownames(matImpulseModelHidden) <- vecImpulseIDs
  rownames(matMuImpulseHidden) <- vecImpulseIDs
  colnames(matMuImpulseHidden) <- names(vecPT)
  
  # c. Periodic flutuations
  # Not coded
  
  # d. Merge data
  matMuHidden <- do.call(rbind, list(matMuConstHidden, matMuImpulseHidden))
  
  # Add size factors
  if(is.null(vecSizeFactorsExternal)){
    print("Setting size factors uniformly =1")
    vecSizeFactorsHidden <- array(1, dim(matMuHidden)[2])
  } else {
    print("Use externally supplied size factors.")
    if(length(vecSizeFactorsExternal) != scaNCells){
      stop("vecSizeFactorsExternal has to be of the length of the number of cells scaNCells.")
    }
    vecSizeFactorsHidden <- vecSizeFactorsExternal
  }
  names(vecSizeFactorsHidden) <- colnames(matMuImpulseHidden)
  
  matMuHidden <- matMuHidden*matrix(vecSizeFactorsHidden,
    nrow=dim(matMuHidden)[1],
    ncol=dim(matMuHidden)[2], byrow=TRUE)
  
  # e. draw dispersions by gene
  vecDispHidden <- runif(dim(matMuHidden)[1])*2+0.05
  matDispHidden <- matrix(vecDispHidden,nrow=dim(matMuHidden)[1],
    ncol=dim(matMuHidden)[2], byrow=FALSE)
  rownames(matDispHidden) <- rownames(matMuHidden)
  colnames(matDispHidden) <- names(vecPT)
  
  # f. add noise - draw from negative binomial
  matSampledDataHidden <- do.call(rbind, lapply(seq(1,dim(matMuHidden)[1]), function(gene){
    sapply(seq(1,scaNCells), function(cell){
      rnbinom(n=1, mu=matMuHidden[gene,cell], size=vecDispHidden[gene])
    })
  }))
  rownames(matSampledDataHidden) <- rownames(matMuHidden)
  colnames(matSampledDataHidden) <- names(vecPT)
  
  # 3. Apply drop out
  # a. Set drop out models
  if(is.null(matDropoutModelExternal)){
    a1 <- c(-2,-2,-4)
    a2 <- c(0.01,0.001,0.01)
    a1 <- array(a1, scaNCells)
    a2 <- array(a2, scaNCells)
    matDropoutLinModelHidden <- cbind(a1,a2)
  } else {
    if(scaNCells!=dim(matDropoutModelExternal)[1]){
      stop("Size of matDropoutModelExternal does not correspond to size of simulated data set.")
    } else {
      matDropoutLinModelHidden <- matDropoutModelExternal
    }
  }
  rownames(matDropoutLinModelHidden) <- names(vecPT)
  
  # b. Draw drop-out rates
  lsDropoutRatesHidden <- lapply(seq(1,scaNCells), function(cell){
    return(1 - evalLogistic(mu=matMuHidden[,cell], 
      a1=matDropoutLinModelHidden[cell,1], 
      a2=matDropoutLinModelHidden[cell,2]))
  })
  matDropoutRatesHidden <- do.call(cbind, lsDropoutRatesHidden)
  rownames(matDropoutRatesHidden) <- rownames(matMuHidden)
  colnames(matDropoutRatesHidden) <- names(vecPT)
  
  # c. Draw drop-outs from rates: Bernoulli experiments
  lsDropoutsHidden <- lapply(seq(1,dim(matDropoutRatesHidden)[1]), function(gene){
    rbinom(n=rep(1, dim(matDropoutRatesHidden)[2]), 
      size=rep(1, dim(matDropoutRatesHidden)[2]),
      prob=matDropoutRatesHidden[gene,])
  })
  matDropoutsHidden <- do.call(rbind, lsDropoutsHidden)
  rownames(matDropoutsHidden) <- rownames(matMuHidden)
  colnames(matDropoutsHidden) <- names(vecPT)
  
  # d. Create observed data: merge hidden data with drop-outs
  matSampledDataObserved <- matSampledDataHidden
  matSampledDataObserved[matDropoutsHidden==1] <- 0
  
  #  Counts
  matSampledCountsObserved <- round(matSampledDataObserved)
  
  # Save simulation
  save(vecPT,file=file.path(dirOutSimulation,"Simulation_vecPT.RData"))
  save(vecConstIDs,file=file.path(dirOutSimulation,"Simulation_vecConstIDs.RData"))
  save(vecImpulseIDs,file=file.path(dirOutSimulation,"Simulation_vecImpulseIDs.RData"))
  save(matImpulseModelHidden,file=file.path(dirOutSimulation,"Simulation_matImpulseModelHidden.RData"))
  save(matMuHidden,file=file.path(dirOutSimulation,"Simulation_matMuHidden.RData"))
  save(vecSizeFactorsHidden,file=file.path(dirOutSimulation,"Simulation_vecSizeFactorsHidden.RData"))
  save(matDispHidden,file=file.path(dirOutSimulation,"Simulation_matDispHidden.RData"))
  save(matSampledDataHidden,file=file.path(dirOutSimulation,"Simulation_matSampledDataHidden.RData"))
  save(matDropoutLinModelHidden,file=file.path(dirOutSimulation,"Simulation_matDropoutLinModelHidden.RData"))
  save(matDropoutRatesHidden,file=file.path(dirOutSimulation,"Simulation_matDropoutRatesHidden.RData"))
  save(matDropoutsHidden,file=file.path(dirOutSimulation,"Simulation_matDropoutsHidden.RData"))
  save(matSampledCountsObserved,file=file.path(dirOutSimulation,"Simulation_matSampledCountsObserved.RData"))
  
  return(list( vecPT=vecPT,
    matSampledCountsObserved=matSampledCountsObserved ))
}