rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState1/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo1 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState2/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo2 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState3/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo3 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState3Rev/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo3Rev <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState1Cont/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo1Cont <- dfImpulseResults

# Plot highly expressed genes
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/R/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/srcPseudoDE_analyseOutput.R")
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState1"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState1"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState2"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState2"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState3"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState3"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState1Cont"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState1Cont"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/NoState3Rev"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState3Rev"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
  
scaThres <- 10^(-10)
sum(dfImpulseResultsNo1$adj.p < scaThres)
sum(dfImpulseResultsNo2$adj.p < scaThres)
sum(dfImpulseResultsNo3$adj.p < scaThres)

# CDF q-values
vecX <- seq(-300,-1,by=1)
vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResultsNo1$adj.p))/log(10) <= thres, na.rm=TRUE)})
vecCDF2 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResultsNo1Cont$adj.p))/log(10) <= thres, na.rm=TRUE)})
vecCDF3 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResultsNo2$adj.p))/log(10) <= thres, na.rm=TRUE)})
vecCDF4 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResultsNo3$adj.p))/log(10) <= thres, na.rm=TRUE)})
vecCDF5 <- sapply(vecX, function(thres){sum(log(as.numeric(dfImpulseResultsNo3Rev$adj.p))/log(10) <= thres, na.rm=TRUE)})
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/ECDF-pvalues_DropStates.pdf",width=7,height=7)
plot(vecX,vecCDF1,
  col="red",pch=4,type="l",
  xlim=c(-30,-1),
  ylim=c(0,max(max(vecCDF1,na.rm=TRUE),max(vecCDF2,na.rm=TRUE))),
  xlab="-log_10(p-value)", 
  ylab=paste0("Cumulative p-value distribution"),
  main=paste0("Lungeputhelium scRNA-seq"))
points(vecX,vecCDF2,
  col="green",pch=4,type="l")
points(vecX,vecCDF3,
  col="blue",pch=4,type="l")
#points(vecX,vecCDF4,
#  col="black",pch=4,type="l")
points(vecX,vecCDF5,
  col="orange",pch=4,type="l")
legend(x="topleft",
  legend=c("Drop state 1",
    "Drop state 1 (continuous branch)",
    "Drop state 2",
    "Drop state 3", 
    "Drop state 3 (reversed pseudotime)"),
  fill=c("red",
    "green",
    "blue",
    "black",
    "orange"))
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/ECDF-pvalues_DropStates_IntervallowQval.pdf",width=7,height=7)
plot(vecX,vecCDF1,
  col="red",pch=4,type="l",
  xlim=c(-20,-1),
  ylim=c(0,max(max(vecCDF1,na.rm=TRUE),max(vecCDF2,na.rm=TRUE))),
  xlab="-log_10(p-value)",
  ylab=paste0("Cumulative p-value distribution"),
  main=paste0("Lungepithelium scRNA-seq"))
points(vecX,vecCDF2,
  col="green",pch=4,type="l")
points(vecX,vecCDF3,
  col="blue",pch=4,type="l")
points(vecX,vecCDF4,
  col="black",pch=4,type="l")
points(vecX,vecCDF5,
  col="orange",pch=4,type="l")
legend(x="topleft",
  legend=c("Drop state 1",
    "Drop state 1 (continuous branch)",
    "Drop state 2",
    "Drop state 3", 
    "Drop state 3 (reversed pseudotime)"),
  fill=c("red",
    "green",
    "blue",
    "black",
    "orange"))
dev.off()
graphics.off()

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
