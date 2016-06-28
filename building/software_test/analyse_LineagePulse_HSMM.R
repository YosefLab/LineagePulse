rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/NoState3/PseudoDE_lsImpulseDE2results.RData")
dfImpulseResults <- lsImpulseDE2results$dfImpulseResults
hist(log(dfImpulseResults$adj.p)/log(10))

# Enrichment analysis
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
load("PseudoDE_HSMM.RData")
load("PseudoDE_HSMM_sample_sheetRAW.RData")
load("PseudoDE_HSMM_gene_annotationRAW.RData")
load("PseudoDE_dfCountsHSMM_SC.RData")
load("PseudoDE_dfFpkmHSMM_SC.RData")

N<-500
vecIDsTopHitsrsem <- as.vector(dfImpulseResults[1:N,"Gene"])
vecIDsTopHitsname <- as.vector(HSMM_gene_annotationRAW[rownames(HSMM_gene_annotationRAW) %in% vecIDsTopHitsrsem,"gene_short_name"])
write(vecIDsTopHitsname, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/clusterruns/output/vecIDsTopHitsname.tab")
# to GO enrichment

# manual inspection
N<-5
vecIDsTopHitsrsem <- as.vector(dfImpulseResults[1:N,"Gene"])
vecIDsTopHitsname <- as.vector(HSMM_gene_annotationRAW[rownames(HSMM_gene_annotationRAW) %in% vecIDsTopHitsrsem,"gene_short_name"])
write(vecIDsTopHitsname, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/clusterruns/output/vecIDsTopHitsnameManual.tab")
