rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/notsmoothed/NoState3/PseudoDE_lsImpulseDE2results.RData")
dfImpulseResults <- lsImpulseDE2results$dfImpulseResults

# Enrichment analysis
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
load("PseudoDE_HSMM.RData")
load("PseudoDE_HSMM_sample_sheetRAW.RData")
load("PseudoDE_HSMM_gene_annotationRAW.RData")
load("PseudoDE_dfCountsHSMM_SC.RData")
load("PseudoDE_dfFpkmHSMM_SC.RData")

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/R/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/srcPseudoDE_analyseOutput.R")
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/outputLineagePulse/notsmoothed/NoState3"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/pdfs/NoState3"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20,
  dfGeneAnnotation=HSMM_gene_annotationRAW)

# For enrichment analysis of top N hits
N<-500
vecIDsTopHitsrsem <- as.vector(dfImpulseResults[1:N,"Gene"])
vecIDsTopHitsname <- as.vector(HSMM_gene_annotationRAW[match(vecIDsTopHitsrsem, rownames(HSMM_gene_annotationRAW)),"gene_short_name"])
write(vecIDsTopHitsname, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/enrichment/vecIDsTopHitsname.tab")
# to GO enrichment

# manual inspection
N<-5
vecIDsTopHitsrsem <- as.vector(dfImpulseResults[1:N,"Gene"])
vecIDsTopHitsname <- as.vector(HSMM_gene_annotationRAW[match(vecIDsTopHitsrsem, rownames(HSMM_gene_annotationRAW)),"gene_short_name"])
write(vecIDsTopHitsname, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/enrichment/vecIDsTopHitsnameManual.tab")

# Fot enrichment analysis of all significant hits at threshold:
qthres <- 10^-5
vecIDsTopHitsrsem <- as.vector(dfImpulseResults[dfImpulseResults$adj.p<=qthres,"Gene"])
vecIDsTopHitsname <- as.vector(HSMM_gene_annotationRAW[rownames(HSMM_gene_annotationRAW) %in% vecIDsTopHitsrsem,"gene_short_name"])
write(vecIDsTopHitsname, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/enrichment/vecIDs_SignificantHits.tab")
# to GO enrichment

# Myogenic factors
strMYOD1 <- rownames(HSMM_gene_annotationRAW)[which(HSMM_gene_annotationRAW$gene_short_name=="MYOD1")]
dfImpulseResults[as.vector(dfImpulseResults$Gene)==strMYOD1,]
strMEF2C <- rownames(HSMM_gene_annotationRAW)[which(HSMM_gene_annotationRAW$gene_short_name=="MEF2C")]
dfImpulseResults[as.vector(dfImpulseResults$Gene)==strMEF2C,]
strMYF5 <- rownames(HSMM_gene_annotationRAW)[which(HSMM_gene_annotationRAW$gene_short_name=="MYF5")]
dfImpulseResults[as.vector(dfImpulseResults$Gene)==strMYF5,]
strMYF6 <- rownames(HSMM_gene_annotationRAW)[which(HSMM_gene_annotationRAW$gene_short_name=="MYF6")]
dfImpulseResults[as.vector(dfImpulseResults$Gene)==strMYF6,]
