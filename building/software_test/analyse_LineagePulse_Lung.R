rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState1/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo1 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState2/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo2 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState3/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo3 <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState3Rev/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo3Rev <- dfImpulseResults
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState1Cont/ImpulseDE2_dfImpulseResults.RData")
dfImpulseResultsNo1Cont <- dfImpulseResults

# Plot highly expressed genes
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/R/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/srcPseudoDE_analyseOutput.R")
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState1"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState1"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState2"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState2"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState3"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState3"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState1Cont"
folderPDFs <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/NoState1Cont"
anlayseOuput(
  folderLineagePulseOutput=folderLineagePulseOutput,
  folderPDFs=folderPDFs,
  strSCMode="continuous",
  scaWindowRadis=20)
folderLineagePulseOutput <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/outputLineagePulse/smoothed1InImpulse/NoState3Rev"
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

# CDF p-values
vecX <- seq(-300,-1,by=1)
vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResultsNo1$p)))/log(10) <= thres, na.rm=TRUE)})
vecCDF1Cont <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResultsNo1Cont$p)))/log(10) <= thres, na.rm=TRUE)})
vecCDF2 <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResultsNo2$p)))/log(10) <= thres, na.rm=TRUE)})
vecCDF3 <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResultsNo3$p)))/log(10) <= thres, na.rm=TRUE)})
vecCDF3Rev <- sapply(vecX, function(thres){sum(log(as.numeric(as.vector(dfImpulseResultsNo3Rev$p)))/log(10) <= thres, na.rm=TRUE)})
scaMaxCDF <- max(c(vecCDF1,vecCDF1Cont,vecCDF2,vecCDF3,vecCDF3Rev), na.rm=TRUE)
# Get number of successfully fitted genes
scaFitted1 <- sum(!is.na(as.numeric(as.vector(dfImpulseResultsNo1$p))))
scaFitted1Cont <- sum(!is.na(as.numeric(as.vector(dfImpulseResultsNo1Cont$p))))
scaFitted2 <- sum(!is.na(as.numeric(as.vector(dfImpulseResultsNo2$p))))
scaFitted3 <- sum(!is.na(as.numeric(as.vector(dfImpulseResultsNo3$p))))
scaFitted3Rev <- sum(!is.na(as.numeric(as.vector(dfImpulseResultsNo3Rev$p))))
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/pdfs/ECDF-pvalues_DropStates.pdf",width=7,height=7)
plot(vecX,vecCDF1,
  col="red",pch=4,type="l",
  xlim=c(-30,-1),
  ylim=c(0,scaMaxCDF),
  xlab="log_10(p-value)", 
  ylab=paste0("Cumulative p-value distribution"),
  main=paste0("Lungeputhelium scRNA-seq"))
points(vecX,vecCDF1Cont,
  col="orange",pch=4,type="l")
points(vecX,vecCDF2,
  col="green",pch=4,type="l")
points(vecX,vecCDF3,
  col="blue",pch=4,type="l")
points(vecX,vecCDF3Rev,
  col="black",pch=4,type="l")
legend(x="topleft",
  legend=c(
    paste0("Drop state 1 [", scaFitted1, "]"),
    paste0("Drop state 1 (continuous branch) [", scaFitted1Cont, "]"),
    paste0("Drop state 2 [", scaFitted2, "]"),
    paste0("Drop state 3 [", scaFitted3, "]"),
    paste0("Drop state 3 (reversed pseudotime) [", scaFitted3Rev, "]")),
  fill=c(
    "red",
    "orange",
    "green",
    "blue",
    "black")
)
dev.off()
graphics.off()

hist(log(dfImpulseResults$adj.p)/log(10))