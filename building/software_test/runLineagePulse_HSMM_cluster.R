rm(list = ls())
NCORES = 16

setwd("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input")
load("PseudoDE_HSMM.RData")
load("PseudoDE_HSMM_sample_sheetRAW.RData")
load("PseudoDE_HSMM_gene_annotationRAW.RData")
load("PseudoDE_dfCountsHSMM_SC.RData")
load("PseudoDE_dfFpkmHSMM_SC.RData")


if(TRUE){
  # Investigate distribution of cells over pseudotime
  vecPTpointsAll_HSMM <- as.vector(pData(HSMM)$Pseudotime)
  names(vecPTpointsAll_HSMM) <- as.vector(rownames(pData(HSMM)))
  vecPTpoints_HSMM <- vecPTpointsAll_HSMM[!is.na(vecPTpointsAll_HSMM)]
  print(paste0("HSMM: Total cells: ",length(vecPTpointsAll_HSMM),", Non NA: ",length(vecPTpoints_HSMM)))
}else{
  # Use real observation time
  vecPTpointsAll_HSMM <- HSMM_sample_sheetRAW$Hours
  names(vecPTpointsAll_HSMM) <- rownames(HSMM_sample_sheetRAW)
  vecPTpoints_HSMM <- vecPTpointsAll_HSMM[!is.na(vecPTpointsAll_HSMM)]
}

# Run PseudoDE
matCounts <- data.matrix(dfCountsHSMM_SC)
vecPT <- vecPTpoints_HSMM
matCounts <- matCounts[,colnames(matCounts) %in% names(vecPTpoints_HSMM)]
vecPT <- vecPT[colnames(matCounts)]

matCountsRd <- matCounts[apply(matCounts,1,function(gene){any(gene>10) & mean(gene)>=1}),]
matCountsRd <- round(matCountsRd)

source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
setwd("data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/output")
lsDEresults <- runPseudoDE(matCounts=matCountsRd,
  vecPseudotime=vecPT,
  K=NULL,
  scaSmallRun=NULL,
  boolPseudotime = TRUE,
  boolContPseudotimeFit=TRUE,
  boolPlotZINBfits=FALSE,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  nProc=NCORES,
  scaMaxiterEM=100)
