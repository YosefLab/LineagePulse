#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++   Simulate a data set  +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
#' Set either vecGeneWiseDropoutRates or matDropoutModelExternal.
#' 
#' @seealso Called by separately by user.
#' 
#' @param scaNCells (scalar) Number of cells in data set.
#' @param scaNConst (scalar) Number of constant expression profiles (genes) in data set.
#' @param scaNLin (scalar) Number of linear expression profiles (genes) in data set.
#' @param scaNImp (scalar) Number of impulse model expression profiles (genes) in data set.
#' @param scaMumax (scalar) [Default 1000]
#' Maximum expression mean parameter to be used.
#' @param scaSDMuAmplitude (scalar) [Default 1]
#' Standard deviation of normal distribution form which the 
#' amplitude change within an impulse trace is drawn.
#' @param vecNormConstExternal (numeric vector number of cells)
#' [Default NULL]
#' Size factors for data set. Size factors are set to 1 if this is
#' not specified (NULL).
#' @param vecDispExternal (numeric vector number of genes)
#' [Default NULL]
#' Dispersion parameters per gene supplied by user.s
#' @param vecGeneWiseDropoutRates (numeric vector number of cells) 
#' [Default NULL] One drop-out rate per gene.
#' @param matDropoutModelExternal (numeric matrix cells x 2) 
#' [Default NULL] External drop-out model, has to have
#' one row for each simulated cell.
#' 
#' @return list (length 2)
#' \itemize{
#' \item pseudotime (numerical vector length number of cells)
#' Pseudotime coordinates (1D) of cells. One scalar per cell.
#' \item counts (matrix genes x cells)
#' Sampled count data of all cells after drop-out.
#' }
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
#' plot(lsSimulatedData$annot$pseudotime, lsSimulatedData$counts[1,])
#' 
#' @author David Sebastian Fischer
#' 
#' @export
simulateContinuousDataSet <- function(
    scaNCells,
    scaNConst,
    scaNLin,
    scaNImp,
    scaMumax=100,
    scaSDMuAmplitude=1,
    vecNormConstExternal=NULL,
    vecDispExternal=NULL,
    vecGeneWiseDropoutRates=NULL,
    matDropoutModelExternal=NULL){
    
    ####
    # Internal functions
    # Evalute impulse model at time points
    evalImpulse <- function(t,beta,t1,t2,h0,h1,h2){
        return(1/h1* (h0+(h1-h0)*1/(1+exp(-beta*(t-t1))))*
                   (h2+(h1-h2)*1/(1+exp(beta*(t-t2)))))
    }
    
    SCA_MIN_MU <- 10^(-5)
    
    ### 0. Check input
    if(!is.null(vecNormConstExternal)) {
        if(scaNConst + scaNLin + scaNImp != length(vecNormConstExternal)){
            stop("scaNConst + scaNLin + scaNImp has to be length(vecNormConstExternal)")
        }
    }
    if(!is.null(vecDispExternal)) {
        if(scaNConst + scaNLin + scaNImp != length(vecDispExternal)){
            stop("scaNConst + scaNLin + scaNImp has to be length(vecDispExternal)")
        }
    }
    if(!is.null(vecGeneWiseDropoutRates)) {
        if(scaNConst + scaNLin + scaNImp != length(vecGeneWiseDropoutRates)){
            stop("scaNConst + scaNLin + scaNImp has to be length(vecGeneWiseDropoutRates)")
        }
    }
    if(!is.null(matDropoutModelExternal)) {
        if(scaNCells != dim(matDropoutModelExternal)[1]){
            stop("scaNCells has to be dim(matDropoutModelExternal)[1]")
        }
        if(dim(matDropoutModelExternal)[1] != 2){
            stop("dim(matDropoutModelExternal)[1] has to be 2 (offsets and mean).")
        }
    }
    
    ### 1. Distribute cells across pseudotime
    vecPT <- seq(0, 1, by=1/(scaNCells-1))
    names(vecPT) <- paste0("cell_",seq(1,scaNCells))
    dfAnnot <- data.frame(
        cell = names(vecPT),
        pseudotime = vecPT,
        row.names = names(vecPT),
        stringsAsFactors = FALSE
    )
    
    ### 2. Create underlying count matrix
    print("Draw mean trajectories")
    if(scaNConst>0) {
        vecConstIDs <- paste0(rep("gene_", scaNConst),
                              c(1:scaNConst))
    } else {
        vecConstIDs <- NULL
    }
    
    if(scaNLin>0) {
        vecLinIDs <- paste0(rep("gene_",scaNLin), 
                            c((scaNConst+1):(scaNConst+scaNLin)))
    } else {
        vecLinIDs <- NULL
    }
    
    if(scaNImp>0) {
        vecImpulseIDs <- paste0(rep("gene_",scaNImp), 
                                c((scaNConst+scaNLin+1):(scaNConst+scaNLin+scaNImp)))
    } else {
        vecImpulseIDs <- NULL
    }
    
    ## Draw means from uniform (first half of genes): one mean per gene
    vecMuConstHidden <- runif(scaNConst)*scaMumax
    vecMuConstHidden[vecMuConstHidden < SCA_MIN_MU] <- SCA_MIN_MU
    matMuConstHidden <- matrix(vecMuConstHidden,
                               nrow=scaNConst,
                               ncol=scaNCells,
                               byrow=FALSE )
    rownames(matMuConstHidden) <- vecConstIDs
    colnames(matMuConstHidden) <- names(vecPT)
    
    ## Draw means from linear model
    vecMuLinStart <- runif(scaNLin)*scaMumax
    vecMuLinStart[vecMuLinStart < SCA_MIN_MU] <- SCA_MIN_MU
    vecMuLinEnd <- vecMuLinStart * abs(1+rnorm(n=scaNLin, mean=1,sd=scaSDMuAmplitude))
    vecMuLinEnd[vecMuLinEnd < SCA_MIN_MU] <- SCA_MIN_MU
    matMuLinHidden <- do.call(rbind, lapply(seq(1, scaNLin), function(i) {
        vecMuLinStart[i] + vecPT/max(vecPT)*(vecMuLinEnd[i] - vecMuLinStart[i])
    }))
    rownames(matMuLinHidden) <- vecLinIDs
    colnames(matMuLinHidden) <- names(vecPT)
    
    ## Draw means from impulse model
    beta <- runif(scaNImp)*2+0.5
    ta <- runif(scaNImp)
    tb <- runif(scaNImp)
    t1 <- ta
    t1[tb < ta] <- tb[tb < ta]
    t2 <- tb
    t2[tb < ta] <- ta[tb < ta]
    h0 <- runif(scaNImp)*scaMumax
    h1 <- h0 * abs(1+rnorm(n=scaNImp, mean=0,sd=scaSDMuAmplitude))
    h2 <- h0 * abs(1+rnorm(n=scaNImp, mean=0,sd=scaSDMuAmplitude))
    h0[h0 < SCA_MIN_MU] <- SCA_MIN_MU
    h1[h1 < SCA_MIN_MU] <- SCA_MIN_MU
    h2[h2 < SCA_MIN_MU] <- SCA_MIN_MU
    if(scaNImp>0){
        lsMuImpulseHidden <- lapply(seq(1,scaNImp), function(gene){
            evalImpulse(t=vecPT,
                        beta=beta[gene],
                        t1=t1[gene],
                        t2=t2[gene],
                        h0=h0[gene],
                        h1=h1[gene],
                        h2=h2[gene])
        })
    } else { lsMuImpulseHidden <- list(matrix(0, nrow=0, ncol=scaNCells)) }
    matImpulseModelHidden <- cbind(beta, h0, h1, h2, t1, t2)
    matMuImpulseHidden <- do.call(rbind, lsMuImpulseHidden)
    rownames(matImpulseModelHidden) <- vecImpulseIDs
    rownames(matMuImpulseHidden) <- vecImpulseIDs
    colnames(matMuImpulseHidden) <- names(vecPT)
    
    ## Merge data
    matMuHidden <- do.call(rbind, list(matMuConstHidden, 
                                       matMuLinHidden,
                                       matMuImpulseHidden))
    
    ## Add size factors
    if(is.null(vecNormConstExternal)){
        print("Setting size factors uniformly =1")
        vecNormConstHidden <- array(1, dim(matMuHidden)[2])
    } else {
        print("Use externally supplied size factors.")
        if(length(vecNormConstExternal) != scaNCells){
            stop("vecNormConstExternal has to be of the length of the number of cells scaNCells.")
        }
        vecNormConstHidden <- vecNormConstExternal
    }
    names(vecNormConstHidden) <- colnames(matMuImpulseHidden)
    
    matMuHidden <- matMuHidden*matrix(vecNormConstHidden,
                                      nrow=dim(matMuHidden)[1],
                                      ncol=dim(matMuHidden)[2], byrow=TRUE)
    
    ## draw dispersions by gene
    print("Draw dispersion")
    if(is.null(vecDispExternal)) {
        vecDispHidden <- 3 + rnorm(n = dim(matMuHidden)[1], mean = 0, sd = 0.5)
        vecDispHidden[vecDispHidden < 0.05] <- 0.05
    } else {
        vecDispHidden <- vecDispExternal
    }
    matDispHidden <- matrix(vecDispHidden,nrow=dim(matMuHidden)[1],
                            ncol=dim(matMuHidden)[2], byrow=FALSE)
    rownames(matDispHidden) <- rownames(matMuHidden)
    colnames(matDispHidden) <- names(vecPT)
    
    ## add noise - draw from negative binomial
    print("Simulate negative binomial noise")
    matSampledDataHidden <- do.call(rbind, lapply(seq(1,dim(matMuHidden)[1]), function(gene){
        sapply(seq(1,scaNCells), function(cell){
            rnbinom(n=1, mu=matMuHidden[gene,cell], size=vecDispHidden[gene])
        })
    }))
    rownames(matSampledDataHidden) <- rownames(matMuHidden)
    colnames(matSampledDataHidden) <- names(vecPT)
    
    ### 3. Apply drop out
    print("Simulate drop-out")
    ## Generate underlying drop-out rate matrix
    if(!is.null(matDropoutModelExternal) & !is.null(vecGeneWiseDropoutRates)) {
        stop("Supply either matDropoutModelExternal or vecGeneWiseDropoutRates.")
    }
    if(!is.null(matDropoutModelExternal)){
        lsDropoutRatesHidden <- lapply(seq(1,scaNCells), function(cell){
            decompressDropoutRateByCell(
                vecDropModel=matDropoutModelExternal[cell,],
                vecMu=matMuHidden[,cell],
                matPiConstPredictors=NULL,
                lsDropModelGlobal=list(strDropModel = "logistic_ofMu",
                                       scaNumCells = scaNCells))
        })
        matDropoutRatesHidden <- do.call(cbind, lsDropoutRatesHidden)
        rownames(matDropoutRatesHidden) <- rownames(matMuHidden)
        colnames(matDropoutRatesHidden) <- names(vecPT)
    } else if(!is.null(vecGeneWiseDropoutRates)) {
        matDropoutRatesHidden <- matrix(
            vecGeneWiseDropoutRates, 
            nrow = length(vecGeneWiseDropoutRates), ncol = dim(matMuHidden)[2], byrow = FALSE)
    } else {
        stop("Supply either matDropoutModelExternal or vecGeneWiseDropoutRates.")
    } 
    ## Draw drop-outs from rates: Bernoulli experiments
    lsDropoutsHidden <- lapply(seq(1,dim(matDropoutRatesHidden)[1]), function(i){
        rbinom(n=rep(1, dim(matDropoutRatesHidden)[2]), 
               size=rep(1, dim(matDropoutRatesHidden)[2]),
               prob=matDropoutRatesHidden[i,])
    })
    matDropoutsHidden <- do.call(rbind, lsDropoutsHidden)
    rownames(matDropoutsHidden) <- rownames(matMuHidden)
    colnames(matDropoutsHidden) <- names(vecPT)
    
    ## Create observed data: merge hidden data with drop-outs
    matSampledDataObserved <- matSampledDataHidden
    matSampledDataObserved[matDropoutsHidden==1] <- 0
    
    ## Round to counts
    matSampledCountsObserved <- round(matSampledDataObserved)
    
    return(list( annot=dfAnnot,
                 counts=matSampledCountsObserved ))
}