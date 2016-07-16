setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
load("PseudoDE_matCountsProc.RData")
load("PseudoDE_matCountsProcFull.RData")
load("PseudoDE_vecPseudotimeProc.RData")  
load("PseudoDE_lsResultsClustering.RData")  
load("PseudoDE_dfAnnotation.RData")
load("PseudoDE_vecSizeFactors.RData")  
load("PseudoDE_vecDispersions.RData")
load("PseudoDE_matDropout.RData")
load("PseudoDE_matProbNB.RData")
load("PseudoDE_matMuCluster.RData")
load("PseudoDE_matMu.RData")
load("PseudoDE_vecEMLogLik.RData")
vecClusterAssignments <- lsResultsClustering$Assignments
names(vecClusterAssignments) <- colnames(matCountsProc)
lsInputToImpulseDE2 <- list(matDropout, matProbNB, matMuCluster, 
  vecClusterAssignments, lsResultsClustering$Centroids )
names(lsInputToImpulseDE2) <- c("matDropout", "matProbNB", "matMuCluster", 
  "vecClusterAssignments", "vecCentroids")

nProc <- 2
strSCMode <- "continuous"
scaWindowRadius <- 20

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/R/ImpulseDE2_main.R")
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
lsImpulseDE2results <- runImpulseDE2(
  matCountData = matCountsProc, 
  dfAnnotation = dfAnnotation,
  strCaseName = "case", 
  strControlName = NULL, 
  strMode = "singlecell",
  strSCMode = strSCMode,
  scaWindowRadius=scaWindowRadius,
  nProc = nProc, 
  Q_value = 10^(-5),
  boolPlotting = TRUE,
  lsPseudo = lsInputToImpulseDE2,
  vecDispersionsExternal = vecDispersions,
  vecSizeFactorsExternal = vecSizeFactors,
  boolRunDESeq2 = FALSE,
  boolSimplePlot = FALSE, 
  boolLogPlot = TRUE,
  scaSmallRun=20)