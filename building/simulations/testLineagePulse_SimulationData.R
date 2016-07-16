# Simulate ZINB data

# Set scenario:
# Pseudotime interval
PTmax <- 100
Ncells <- 400
Mumax <- 10000
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
Nconst <- 100
vecMuConst <- runif(Nconst)*Mumax
matHiddenDataConst <- matrix(vecMuConst, nrow=Nconst, ncol=Ncells, byrow=FALSE)

# b. Draw means from impulse model
NImp <- 100
beta <- rep(1,NImp)
t1 <- seq(0, PTmax, by=PTmax/(NImp-1))
t2 <- seq(1, 1+PTmax*2, by=2*PTmax/(NImp-1))
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
vecPhi <- runif(dim(matHiddenData)[1])*10+0.1

# f. add noise - draw from negative binomial
matSampledData <- do.call(rbind, lapply(seq(1,dim(matHiddenData)[1]), function(gene){
  rnbinom(n=dim(matHiddenData)[2],
    size=vecPhi[gene],
    mu=matHiddenData[gene,])
}))

# 3. Apply drop out
# a. Set drop out models
a1 <- c(-2,-2,-20)
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
matDropoutRates <- matrix(0.5, nrow=dim(matHiddenData)[1], ncol=Ncells)

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

plot(vecPT,matObsData[128,])

# 4. Run LineagePulse
NCORES = 2

# 1. Counts
matData <- round(matObsData)
colnames(matData) <- paste0("_",seq(1,dim(matData)[2]))
rownames(matData) <- paste0("_",seq(1,dim(matData)[1]))
names(vecPT) <- colnames(matData)

source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
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
  boolDEAnalysisModelFree = TRUE,
  boolPlotZINBfits=FALSE,
  scaMaxiterEM=100,
  nProc=NCORES,
  verbose=TRUE )

scaGeneID <- 130
plot(1:400, matCountsProc[scaGeneID,])
points(1:400, matHiddenData[scaGeneID,], col="green")
points(1:400, matMu[scaGeneID,], col="red")
points(1:400, matMu[scaGeneID,], col="yellow")

scaGeneID<- 150
vecDispersions[scaGeneID]-vecDispersions1[scaGeneID]
max(matZ1[scaGeneID,]-matZ[scaGeneID,])
max(matMu1[scaGeneID,]-matMu[scaGeneID,])
plot(1:400, matMu[scaGeneID,], col="red")
points(1:400, matMu1[scaGeneID,], col="yellow")
plot(matMu[scaGeneID,], matZ[scaGeneID,], col="green")
points(matMu[scaGeneID,], matDropout[scaGeneID,])
