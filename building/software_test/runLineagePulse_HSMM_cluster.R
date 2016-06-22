rm(list = ls())
NCORES = 16

setwd("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input")
load("PseudoDE_matCountsRd.RData")
load("PseudoDE_vecPT.RData")

source("/data/yosef2/users/fischerd/code/LineagePulse/building/code_files/PseudoDE_main.R")
setwd("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/output")
lsDEresults <- runPseudoDE(matCounts=matCountsRd,
  vecPseudotime=vecPT,
  K=6,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit=FALSE,
  boolPlotZINBfits=FALSE,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  nProc=NCORES,
  scaMaxiterEM=100)
