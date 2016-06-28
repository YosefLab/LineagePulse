rm(list = ls())
NCORES = 16

# 1. Counts
dfscRNA <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input/HSMM_matCountsNoState3.tab", sep="\t", header=F)
dfGeneIDs <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input/HSMM_matCountsNoState3_genes.tab",sep="\t",header=F)
dfCellIDs <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input/HSMM_matCountsNoState3_samples.tab",sep="\t",header=F)
matData <- round(t(apply(dfscRNA,1,as.numeric)))
colnames(matData) <- dfCellIDs$V1
rownames(matData) <- dfGeneIDs$V1

# 2. Annotation
dfAnnotation <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/input/HSMM_dfAnnotationNoState3.tab",header=T, stringsAsFactors = F)
vecPT <- dfAnnotation$Pseudotime
names(vecPT) <- rownames(dfAnnotation)

source("/data/yosef2/users/fischerd/code/LineagePulse/building/code_files/PseudoDE_main.R")
setwd("/data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/output/NoState3")
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
  boolPlotZINBfits=TRUE,
  scaMaxiterEM=100,
  nProc=NCORES,
  verbose=TRUE )