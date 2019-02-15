#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++    Likelihood functions    ++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute loglikelihood of negative binomial model
#' for a vector of counts.
#' 
#' This liklihood function is appropriate for sequencing data without drop 
#' outs. This is the core function used for every likelihood evaluation 
#' in LineagePulse, such as maximum likelihood-based estimation. 
#' It operates on a vector of counts, such as observations of a gene. 
#' Note that for the sake of numerical stability, 
#' lower bounds on loglikelihood terms are implemented.
#' 
#' @param vecCounts (count vector number of samples)
#' Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#' Negative binomial mean parameter.
#' @param vecDisp (scalar vector number of samples) 
#' Negative binomial dispersion parameters.
#' @param vecidxNotZero (bool vector number of samples)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of samples)
#' Whether observation is zero.
#' 
#' @return scaLogLik (scalar) Likelihood under zero-inflated
#' negative binomial model.
#' 
#' @author David Sebastian Fischer
evalLogLikNB <- function(
    vecCounts,
    vecMu,
    vecDisp, 
    vecidxNotZero, 
    vecidxZero){  
    
    scaLogLikPrecLim <- -700
    
    # Likelihood of zero counts:
    vecLogLikZeros <- vecDisp[vecidxZero] * 
        (log(vecDisp[vecidxZero]) - log(vecDisp[vecidxZero] + vecMu[vecidxZero])) 
    vecLogLikZeros[vecLogLikZeros < scaLogLikPrecLim] <- scaLogLikPrecLim
    # Likelihood of non-zero counts:
    vecLogLikNonzeros <- dnbinom(
            vecCounts[vecidxNotZero], 
            mu=vecMu[vecidxNotZero], 
            size=vecDisp[vecidxNotZero], 
            log=TRUE)
    vecLogLikNonzeros[vecLogLikNonzeros < scaLogLikPrecLim] <- 
        scaLogLikPrecLim
    
    return(sum(vecLogLikZeros) + sum(vecLogLikNonzeros))
}


#' Compute loglikelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). This
#' is the core function used for every likelihood evaluation 
#' in LineagePulse, such as maximum likelihood-based estimation. 
#' It operates on a vector of counts, such as observations of a gene. 
#' Note that for the sake of numerical stability, 
#' lower bounds on loglikelihood terms are implemented.
#' 
#' @param vecCounts (count vector number of samples)
#' Observed read counts, not observed are NA.
#' @param vecMu (vector number of samples) 
#' Negative binomial mean parameter.
#' @param vecDisp (scalar vector number of samples) 
#' Negative binomial dispersion parameters.
#' @param vecPi (probability vector number of samples) 
#' Drop-out rate estimates.
#' @param vecidxNotZero (bool vector number of samples)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of samples)
#' Whether observation is zero.
#' 
#' @return scaLogLik (scalar) Likelihood under zero-inflated
#' negative binomial model.
#' 
#' @author David Sebastian Fischer
evalLogLikZINB <- function(
    vecCounts,
    vecMu,
    vecDisp, 
    vecPi, 
    vecidxNotZero, 
    vecidxZero){  
    
    # Note on handling very low probabilities: vecLikZeros
    # typically does not have zero elements as it has the 
    # the summand drop-out rate. Also the log cannot be
    # broken up over the sum to dnbinom. In contrast to that,
    # the log is taken in dnbinom for vecLikNonzeros to avoid 
    # zero probabilities. Zero probabilities are handled
    # through substitution of the minimal probability under
    # machine precision.
    scaLogLikPrecLim <- -700
    
    # Likelihood of zero counts:
    vecLogLikZeros <- 
        log((1-vecPi[vecidxZero])*
                (vecDisp[vecidxZero] / 
                     (vecDisp[vecidxZero] + 
                          vecMu[vecidxZero]))^vecDisp[vecidxZero] +
                vecPi[vecidxZero])
    vecLogLikZeros[vecLogLikZeros < scaLogLikPrecLim] <- scaLogLikPrecLim
    scaLogLikZeros <- sum(vecLogLikZeros)
    # Likelihood of non-zero counts:
    vecLogLikNonzeros <- log(1-vecPi[vecidxNotZero]) +
        dnbinom(
            vecCounts[vecidxNotZero], 
            mu=vecMu[vecidxNotZero], 
            size=vecDisp[vecidxNotZero], 
            log=TRUE)
    vecLogLikNonzeros[vecLogLikNonzeros < scaLogLikPrecLim] <- 
        scaLogLikPrecLim
    scaLogLikNonzeros <- sum(vecLogLikNonzeros)
    # Compute likelihood of all data:
    scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
    # Maximise log likelihood: Return likelihood as value 
    # to optimisation routine
    return(scaLogLik)
}

#' Wrapper for log likelihood of (zero-inflated) negative binomial model
#' for a vector of counts.
#'
#' Chooses whether to evaluate ZINB or NB model based on input.
#'
#' @param vecCounts (count vector number of cells)
#' Observed read counts, not observed are NA.
#' @param vecMu (numeric vector number of cells) 
#' Negative binomial mean parameter.
#' @param vecNormConst (numeric vector number of cells) 
#' Model scaling factors, one per cell.
#' @param vecDisp (numeric vector number of cells) 
#' Negative binomial dispersion parameters.
#' @param vecPi (probability vector number of cells) 
#' Drop-out rate estimates. Set to NULL to evaluate NB
#' likelihood.
#' @param vecidxNotZero (bool vector number of cells)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of cells)
#' Whether observation is zero.
#' 
#' @return scaLogLik (scalar) Likelihood under zero-inflated
#' negative binomial model.
#' 
#' @author David Sebastian Fischer
evalLogLikGene <- function(
    vecCounts,
    vecMu,
    vecNormConst,
    vecDisp, 
    vecPi,
    vecidxNotZero, 
    vecidxZero ){
    
    if(!is.null(vecPi)){
        scaLogLik <- evalLogLikZINB(
            vecCounts=vecCounts,
            vecMu=vecMu*vecNormConst,
            vecDisp=vecDisp, 
            vecPi=vecPi,
            vecidxNotZero=vecidxNotZero, 
            vecidxZero=vecidxZero )
    } else {
        scaLogLik <- evalLogLikNB(
            vecCounts=vecCounts,
            vecMu=vecMu*vecNormConst,
            vecDisp=vecDisp, 
            vecidxNotZero=vecidxNotZero, 
            vecidxZero=vecidxZero )
    }
    return(scaLogLik)
}

#' Wrapper for log likelihood of zero-inflated negative binomial model
#' for a vector of counts.
#' 
#' NOT YET SUPPORTED.
#' LINEAGEPULSE CODE WILL BE EXTENDED BY MODULAR FUNCTIONALITIES
#' AND THIS IS ONE INSTANCE OF A PLACEHOLDER USED FOR DEVELOPING.
#' Chooses whether to evaluate ZINB or NB model based on input.
#'
#' @param vecCounts (count vector number of cells)
#' Observed read counts, not observed are NA.
#' @param matMuParam (numeric matrix number of cells x number of mixtures) 
#' Negative binomial mean parameter matrix with one mean per
#' cell and per mixture.
#' @param vecNormConst (numeric vector number of cells) 
#' Model scaling factors, one per cell.
#' @param vecDisp (numeric vector number of cells) 
#' Negative binomial dispersion parameters.
#' @param vecPi (probability vector number of cells) 
#' Drop-out rate estimates. Set to NULL to evaluate NB model.
#' @param vecidxNotZero (bool vector number of cells)
#' Whether observation is larger than zero.
#' @param vecidxZero (bool vector number of cells)
#' Whether observation is zero.
#' @param scaNCells (scalar)
#' Number of cells (auxillary parameter, this is length of vecCounts).
#' 
#' @return NULL
#' This will be: scaLogLik (scalar) Likelihood under zero-inflated
#' negative binomial model.
#' 
#' @author David Sebastian Fischer
evalLogLikGeneMM <- function(
    vecCounts,
    matMuParam,
    vecNormConst,
    vecDisp, 
    vecPi,
    vecidxNotZero, 
    vecidxZero,
    scaNCells){
    
    stop("LineagePulse ERROR: evalLogLikGeneMM() not yet supported. ",
         "Contact developer.")
    return(NULL)
}

#' Wrapper for log likelihood of (zero-inflated) negative binomial model
#' for a matrix of counts (parallelised).
#' 
#' This likelihood function is a wrapper computes loglikelihood
#' of entire data set by parallelising loglikelihood computation
#' over genes.
#' 
#' @seealso Called directly by \code{fitZINB} to track
#' convergence of estimation iteration on entire data set.
#'
#' @param matCounts (count matrix genes x cells)
#' Observed read counts, not observed are NA.
#' @param lsMuModel (list)
#' Object containing description of gene-wise mean parameter models.
#' @param lsDispModel (list)
#' Object containing description of gene-wise dispersion parameter models.
#' @param lsDropModel (list)
#' Object containing description of cell-wise drop-out parameter models.
#' Set to NULL to evaluate NB likelihood as oppose to ZINB.
#' @param matWeights (numeric matrix cells x mixtures) [Default NULL]
#' Assignments of cells to mixtures (for strMuModel="MM").
#' @param boolConstModelOnMMLL (bool) [Default FALSE]
#' Whether to evaluate constant gene-wise mean model 
#' on mixture model likelihood.
#' This is necessary for numeric stability when testing whether mixture model
#' fits are worse than constant fits.
#'
#' @return vecLogLik (vector length number of genes) 
#' Loglikelihood of each gene under zero-inflated negative binomial model.
#' 
#' @author David Sebastian Fischer
evalLogLikMatrix <- function(
    matCounts,
    lsMuModel,
    lsDispModel, 
    lsDropModel,
    matWeights=NULL,
    boolConstModelOnMMLL=FALSE){
    
    # Make sure usage of boolConstModelOnMMLL is correct
    if(boolConstModelOnMMLL &
       (is.null(matWeights) | 
        lsMuModel$lsMuModelGlobal$strMuModel!="constant")){
        stop("boolConstModelOnMMLL usage wrong in evalLogLikMatrix().")
    }
    
    # Parallelize likelihood computation across genes.
    vecLogLik <- unlist(
        bplapply( seq_len(nrow(matCounts)), function(i){
            vecCounts <- matCounts[as.double(i),] 
            # Iterator index on sparse matrix must be double!
            if(lsMuModel$lsMuModelGlobal$strMuModel=="MM"){
                # Expand cells x mixtures parameter matrices for given gene.
                matMuParam <- decompressMuByGeneMM(
                    vecMuModel=lsMuModel$matMuModel[i,],
                    lsvecBatchModel=lapply(
                        lsMuModel$lsmatBatchModel, 
                        function(mat) mat[i,] ),
                    lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                    vecInterval=NULL)
                
                matDispParam <- decompressDispByGeneMM(
                    lsDispModel$matDispModel[i,],
                    lsvecBatchModel=lapply(
                        lsDispModel$lsmatBatchModel, 
                        function(mat) mat[i,] ),
                    lsDispModelGlobal=lsMuModel$lsDispModelGlobal,
                    vecInterval=NULL )
                
                if(lsDropModel$lsDropModelGlobal$strDropModel != "none") {
                    # ZINB noise model
                    matDropParam <- do.call(
                        cbind, lapply(seq_len(ncol(matWeights)), function(m){
                            decompressDropoutRateByGene(
                                matDropModel=lsDropModel$matDropoutLinModel,
                                vecMu=matMuParam[,m],
                                vecPiConstPredictors=
                                    lsDropModel$matPiConstPredictors[i,],
                                lsDropModelGlobal=lsDropModel$lsDropModelGlobal )
                        }))
                } else {
                    # NB noise model
                    matDropParam <- NULL
                }
                # Evalute mixture model likelihood given parameter matrices.
                scaLL <- evalLogLikGeneMM(
                    vecCounts=vecCounts,
                    matMuParam=matMuParam,
                    vecNormConst=lsMuModel$lsMuModelGlobal$vecNormConst,
                    matDispParam=matDispParam,
                    matDropParam=matDropParam,
                    matWeights=matWeights,
                    vecidxNotZero= which(!is.na(vecCounts) & vecCounts>0), 
                    vecidxZero= which(!is.na(vecCounts) & vecCounts==0),
                    scaNCells=length(vecCounts) )
                
            } else {
                # Decompress parameter vectors for given gene.
                vecMuParam <- decompressMeansByGene(
                    vecMuModel=lsMuModel$matMuModel[i,],
                    lsvecBatchModel=lapply(lsMuModel$lsmatBatchModel, 
                                           function(mat) mat[i,] ),
                    lsMuModelGlobal=lsMuModel$lsMuModelGlobal,
                    vecInterval=NULL )
                vecDispParam <- decompressDispByGene(
                    vecDispModel=lsDispModel$matDispModel[i,],
                    lsvecBatchModel=lapply(lsDispModel$lsmatBatchModel, 
                                           function(mat) mat[i,] ),
                    lsDispModelGlobal=lsDispModel$lsDispModelGlobal,
                    vecInterval=NULL )
                if(lsDropModel$lsDropModelGlobal$strDropModel != "none") {
                    # ZINB noise model
                    vecPiParam <- decompressDropoutRateByGene(
                        matDropModel=lsDropModel$matDropoutLinModel,
                        vecMu=vecMuParam,
                        vecPiConstPredictors=lsDropModel$matPiConstPredictors[i,],
                        lsDropModelGlobal=lsDropModel$lsDropModelGlobal)
                } else {
                    # NB noise model
                    vecPiParam <- NULL
                }
                
                # Evaluate loglikelihood of given gene.
                if(boolConstModelOnMMLL) {
                    scaLL <- evalLogLikGeneMM(
                        vecCounts=vecCounts,
                        matMuParam=do.call(
                            cbind, lapply(seq_len(ncol(matWeights)), 
                                          function(m) vecMuParam )),
                        vecNormConst=lsMuModel$lsMuModelGlobal$vecNormConst,
                        matDispParam=do.call(
                            cbind, lapply(seq_len(ncol(matWeights)), 
                                          function(m) vecDispParam )),
                        matDropParam=do.call(
                            cbind, lapply(seq_len(ncol(matWeights)), 
                                          function(m) vecPiParam )),
                        matWeights=matWeights,
                        vecidxNotZero= which(!is.na(vecCounts) & vecCounts>0), 
                        vecidxZero= which(!is.na(vecCounts) & vecCounts==0),
                        scaNCells=length(vecCounts) )
                } else {
                    scaLL <- evalLogLikGene(
                        vecCounts=vecCounts,
                        vecMu=vecMuParam,
                        vecNormConst=lsMuModel$lsMuModelGlobal$vecNormConst,
                        vecDisp=vecDispParam, 
                        vecPi=vecPiParam,
                        vecidxNotZero= which(!is.na(vecCounts) & vecCounts>0), 
                        vecidxZero= which(!is.na(vecCounts) & vecCounts==0) )
                }
            }
            return(scaLL)
        })
    )
    return(vecLogLik)
}