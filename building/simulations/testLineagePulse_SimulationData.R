# Simulate ZINB data

# Set scenario:
# Pseudotime interval
PTmax <- 100
Ncells <- 100
Mumax <- 1000
evalImpulse <- function(t,beta,t1,t2,h0,h1,h2){
  return(1/h1*(h0+(h1-h0)*1/(1+exp(-beta*(t-t1))))*(h2+(h1-h2)*1/(1+exp(beta*(t-t2)))))
}
evalLogistic <- function(mu, a1, a2){
  return(1/(1+exp(-(a1+a2*mu))))
}

# 1. Distribute cells across pseudotime
vecPT <- seq(0, PTmax, by=PTmax/(Ncells-1))

# 2. Create hidden data set
# a. Draw means from uniform (first half of genes): one mean per gene
Nconst <-10
vecMuConst <- runif(Nconst)*Mumax
matHiddenDataConst <- matrix(vecMuConst,
  nrow=Nconst,
  ncol=Ncells,
  byrow=FALSE )

# b. Draw means from impulse model
NImp <- 10
beta <- runif(NImp)*2+0.5
#t1 <- seq(0, PTmax, by=PTmax/(NImp-1))
#t2 <- seq(1, 1+PTmax*2, by=2*PTmax/(NImp-1))
t1 <- runif(NImp)*PTmax
t2 <- runif(NImp)*PTmax
h0 <- runif(NImp)*Mumax
h1 <- runif(NImp)*Mumax
h2 <- runif(NImp)*Mumax
lsHiddenDataImp <- lapply(seq(1,NImp), function(gene){
  evalImpulse(t=vecPT,
    beta=beta[gene],
    t1=t1[gene],
    t2=t2[gene],
    h0=h0[gene],
    h1=h1[gene],
    h2=h2[gene])
})
matHiddenDataImp <- do.call(rbind, lsHiddenDataImp)

# c. Periodic flutuations

# d. Merge data
matHiddenData <- do.call(rbind, list(matHiddenDataConst, matHiddenDataImp))

# e. draw dispersions by gene
vecPhi <- runif(dim(matHiddenData)[1]*5)+0.0001

# f. add noise - draw from negative binomial
matSampledData <- do.call(rbind, lapply(seq(1,dim(matHiddenData)[1]), function(gene){
  sapply(seq(1,dim(matHiddenData)[2]), function(cell){
    rnbinom(n=1, mu=matHiddenData[gene,cell], size=vecPhi[gene])
  })
}))

# 3. Apply drop out
# a. Set drop out models
a1 <- c(-2,-2,-4)
a2 <- c(0.01,0.001,0.01)
a1 <- array(a1, Ncells)
a2 <- array(a2, Ncells)

# b. Draw drop-out rates
lsDropoutRates <- lapply(seq(1,Ncells), function(cell){
  return(1 - evalLogistic(mu=matHiddenData[,cell], 
    a1=a1[cell], 
    a2=a2[cell]))
})
matDropoutRates <- do.call(cbind, lsDropoutRates)

# c. Draw drop-outs from rates: Bernoulli experiments
lsDropouts <- lapply(seq(1,dim(matDropoutRates)[1]), function(gene){
  rbinom(n=rep(1, dim(matDropoutRates)[2]), 
  size=rep(1, dim(matDropoutRates)[2]),
  prob=matDropoutRates[gene,])
})
matDropouts <- do.call(rbind, lsDropouts)

# d. Create observed data: merge hidden data with drop-outs
matObsData <- matSampledData
matObsData[matDropouts==1] <- 0

plot(vecPT,matObsData[12,])

# 4. Run LineagePulse
NCORES = 2

# 1. Counts
matData <- round(matObsData)
colnames(matData) <- paste0("_",seq(1,dim(matData)[2]))
rownames(matData) <- paste0("_",seq(1,dim(matData)[1]))
names(vecPT) <- colnames(matData)

source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
evalLogLikZINB_LinPulse_comp <- cmpfun(evalLogLikZINB_LinPulse)
evalLogLikSmoothZINB_LinPulse_comp <- cmpfun(evalLogLikSmoothZINB_LinPulse)
evalLogLikMuZINB_LinPulse_comp <- cmpfun(evalLogLikMuZINB_LinPulse)
evalLogLikDispZINB_LinPulse_comp <- cmpfun(evalLogLikDispZINB_LinPulse)
evalLogLikPiZINB_LinPulse_comp <- cmpfun(evalLogLikPiZINB_LinPulse)
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
lsDEresults <- runPseudoDE(
  matCounts=matData,
  vecPseudotime=vecPT,
  K=6,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit=TRUE,
  boolOneDispPerGene = TRUE,
  scaWindowRadius=20,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  boolPlotZINBfits=FALSE,
  scaMaxiterEM=100,
  nProc=NCORES,
  verbose=TRUE )
load("PseudoDE_vecDispersions.RData")
vecDispersionsInferred <- vecDispersions
load("PseudoDE_matDropout.RData")
load("PseudoDE_matDropoutLinModel.RData")
matDropoutInferred <- matDropout
load("PseudoDE_matMu.RData")

#---
matCountsProc=matData
vecPseudotime=vecPT
K=6
scaSmallRun=NULL
boolPseudotime = TRUE
boolContPseudotimeFit=TRUE
boolOneDispPerGene = TRUE
scaWindowRadius=20
boolDEAnalysisImpulseModel = TRUE
boolDEAnalysisModelFree = FALSE
boolPlotZINBfits=FALSE
scaMaxiterEM=100
nProc=NCORES
lsResultsClustering=NULL
boolOneDispPerGene=TRUE
verbose=TRUE
matConstPredictorsPi=NULL
scaMaxiterEM=100
boolSuperVerbose=FALSE
vecSizeFactors <- rep(1,dim(matCountsProc)[2])
#--
#mean estimation
vecCounts=matCountsProc[i,vecInterval]
vecMu=vecMu[vecInterval]
vecDisp=matDispersions[i,vecInterval]
vecNormConst=vecSizeFactors[vecInterval]
vecDropoutRateEst=vecDropout[vecInterval]
vecProbNB=1-matZ[i,]
vecPredictorsPi=cbind(1,NA,matConstPredictorsPi[i,])
matLinModelPi=matLinModelPi
scaTarget=match(j,vecInterval)
vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0
vecboolZero=vecCounts==0
scaTheta <- 1
#----

# Open .pdf
matInferredData <- lsDEresults$lsImpulseDE2results$
  lsImpulseFits$values_case
vecGeneIDs <- rownames(matData)
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
pdf(paste("LineagePulse_SimulatedGenes.pdf",sep=""),height=6.0,width=9.0)

# Define grid for printing plots
if (length(vecGeneIDs) == 1){
  par(mfrow=c(1,1), xpd=TRUE)
  scaLegendInset <- -0.15
} else if (length(vecGeneIDs) <= 4){
  par(mfrow=c(2,2), xpd=TRUE)
  scaLegendInset <- -0.35
} else if (length(vecGeneIDs) <= 6){
  par(mfrow=c(2,3), xpd=TRUE)
  scaLegendInset <- -0.3
} else {
  par(mfrow=c(3,3), xpd=TRUE)
  scaLegendInset <- -0.65
}
# Time points for plotting of impulse model
vecX <- seq(min(vecPT),max(vecPT),min(0.1,(max(vecPT)-min(vecPT))/100))
# Find elements in vecX corresponding (closest) to observed time point
indVecXObs <- unlist(lapply( 
  vecPT, function(t){match( min(abs(vecX-t)), abs(vecX-t) )} 
))
for (scaGeneID in seq(1,dim(matData)[1])){
  plot(1:Ncells, log(matData[scaGeneID,]+1)/log(2), 
    col="black", 
    pch=3,
    xlab="Pseudotime",
    ylab=paste0("log_2 Impulse fit and count data"),
    main=rownames(matData)[scaGeneID])
  points(1:Ncells, log(matHiddenData[scaGeneID,]+1)/log(2), 
    col="red",
    type="l")
  points(1:Ncells, log(matInferredData[scaGeneID,]+1)/log(2),
    col="blue",
    type="l")
  #legend(x="topright",
  #  legend=c("True model","Sampled data after drop-out","Inferred model"),
  #  fill=c("red","black","blue"))
}
# Close .pdf
dev.off()

# Overall deviation comparison: LRT
# Compute loglikelihood of impulse imputed
# Recover zero predictions:
matHiddenDataTemp <- matHiddenData
matHiddenDataTemp[matHiddenDataTemp==0] <- 10^(-4)
vecLLImpulseImputedByGene <- sapply(seq(1,dim(matData)[1]), function(gene){
  sum(dnbinom(x=matData[gene,],
    mu=matHiddenDataTemp[gene,],
    size=vecPhi[gene],
    log=TRUE), na.rm=TRUE)
})
vecLLImpulseImputedByGene <- unlist(lapply( seq(1,dim(matData)[1]), function(i){
  evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matData[i,],
    vecMu=matHiddenDataTemp[i,],
    vecSizeFactors=array(1,dim(matData)[2]),
    vecDispEst=array(vecPhi[i],dim(matData)[2]), 
    vecDropoutRateEst=matDropoutRates[i,],
    vecboolNotZeroObserved=matData[i,]>0 & !is.na(matData[i,]), 
    vecboolZero=matData[i,]==0,
    scaWindowRadius=20)
}))
# Compute loglikelihood of baseline imputed
# Recover zero predictions:
matInferredDataTemp <- matInferredData
matInferredDataTemp[matInferredDataTemp==0] <- 10^(-4)
vecLLBaselineImputedByGene <- sapply(seq(1,dim(matData)[1]), function(gene){
  sum(dnbinom(x=matData[gene,],
    mu=matInferredDataTemp[gene,],
    size=vecPhi[gene],
    log=TRUE), na.rm=TRUE)
})
vecLLBaselineImputedByGene <- unlist(lapply( seq(1,dim(matData)[1]), function(i){
  evalLogLikSmoothZINB_LinPulse_comp(vecCounts=matData[i,],
    vecMu=matInferredDataTemp[i,],
    vecSizeFactors=array(1,dim(matData)[2]),
    vecDispEst=array(vecDispersionsInferred[i],dim(matData)[2]), 
    vecDropoutRateEst=matDropoutInferred[i,],
    vecboolNotZeroObserved=(matData[i,]>0 & !is.na(matData[i,])), 
    vecboolZero=(matData[i,]==0),
    scaWindowRadius=20)
}))
# LRT
scaQThres <- 10^(-3)
vecIDs <- rownames(matData)
vecboolIDs <- rownames(matData) %in% vecIDs
vecLRT <- vecLLImpulseImputedByGene[vecboolIDs]-vecLLBaselineImputedByGene[vecboolIDs]
dfLRT <- data.frame(value=vecLRT)
gHistLRT <- ggplot( dfLRT, aes(value)) +
  geom_histogram(alpha = 1, position = 'identity') +
  ggtitle(paste0("log LRT value")) +
  xlab("Difference in loglikelihood") +
  ylab("Frequency")
graphics.off()
pdf("Hist_LRTvalues.pdf")
print(gHistLRT)
dev.off()
print(paste0("Mean logLRT value: ", mean(vecLRT)))
print(paste0("Stdv logLRT value: ", sd(vecLRT)))

# Analyse overfitting
# View as function of drop-out
vecNumDropouts <- apply(matDropouts, 1, sum)
plot(vecNumDropouts,vecLRT)
# View as function of mean
vecSampleMean <- apply(matSampledData, 1, mean)
plot(vecSampleMean,vecLRT)
