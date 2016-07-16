setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files")
source("srcPseudoDE_evalLogLikZINB.R")
source("srcPseudoDE_fitZINB.R")
# Compile function
evalLogLikZINB_LinPulse_comp <- cmpfun(evalLogLikZINB_LinPulse)
evalLogLikSmoothZINB_LinPulse_comp <- cmpfun(evalLogLikSmoothZINB_LinPulse)
evalLogLikMuZINB_LinPulse_comp <- cmpfun(evalLogLikMuZINB_LinPulse)
evalLogLikDispZINB_LinPulse_comp <- cmpfun(evalLogLikDispZINB_LinPulse)
scaLogLikNewSmooth <- sum(sapply( seq(1,scaNumGenes), function(i){
  sum(sapply(seq(1,scaNumCells), 
    function(j){
      scaindIntervalStart <- max(1,j-scaWindowRadius)
      scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
      vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
      scaLogLikCell <- sum(log(matDropout[i,vecInterval]*(matCountsProc[i,vecInterval]==0) +
          (1-matDropout[i,vecInterval])*
          dnbinom(matCountsProc[i,vecInterval],
            mu=rep(matMu[i,j], length(vecInterval)), 
            size=rep(matDispersions[i,j], length(vecInterval)), 
            log=FALSE)))
      return(scaLogLikCell)
    }))
}))
scaLogLikNewSmoothFun <- sum(unlist(
  bplapply( seq(1,scaNumGenes), function(i){
    evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matCountsProc[i,],
      vecMu=matMu[i,],
      vecDispEst=matDispersions[i,], 
      vecDropoutRateEst=matDropout[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,],
      scaWindowRadius=scaWindowRadius)
  })
))
scaLogLikNewSmooth
scaLogLikNewSmoothFun