rm(list = ls())
NCORES = 16

# 1. Counts
dfscRNA <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/input/Lung_matCountsNoState1Cont.tab", sep="\t", header=F)
dfGeneIDs <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/input/Lung_matCountsNoState1Cont_genes.tab",sep="\t",header=F)
dfCellIDs <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/input/Lung_matCountsNoState1Cont_samples.tab",sep="\t",header=F)
matData <- round(t(apply(dfscRNA,1,as.numeric)))
colnames(matData) <- dfCellIDs$V1
rownames(matData) <- dfGeneIDs$V1

# 2. Annotation
dfAnnotation <- read.table("/data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/input/Lung_dfAnnotationNoState1Cont.tab",header=T, stringsAsFactors = F)
vecPT <- dfAnnotation$Pseudotime
names(vecPT) <- rownames(dfAnnotation)

source("/data/yosef2/users/fischerd/code/LineagePulse/building/code_files/PseudoDE_main.R")
setwd("/data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/output/NoState1Cont")
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