source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/NoState3/PseudoDE_dfAnnotation.RData")
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/NoState3/PseudoDE_matCountsProcFull.RData")
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/NoState3/PseudoDE_vecPseudotimeProc.RData")
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
lsDEresults <- runPseudoDE(
  matCounts=matCountsProcFull[1:200,],
  vecPseudotime=vecPseudotimeProc,
  K=6,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit=TRUE,
  boolOneDispPerGene = TRUE,
  scaWindowRadius=20,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  boolPlotZINBfits=TRUE,
  scaMaxiterEM=4,
  nProc=2,
  verbose=TRUE )

load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_lsImpulseDE2results.RData")
load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matDropout.RData")
names(lsImpulseDE2results)
head(lsImpulseDE2results$dfImpulseResults)
head(lsImpulseDE2results$lsImpulseFits$parameters_case)
