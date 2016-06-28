rm(list = ls())
NCORES = 3

# Load data
dfCountsHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/rsem_readCountsTable.txt", header=FALSE, colClasses="numeric")
dfFpkmHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/rsem_fpkmTable.txt", header=FALSE, colClasses="numeric")
dfCellsHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/cell_list.txt", header=FALSE, colClasses="character")
dfGenesHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/gene_list.txt", header=FALSE, colClasses="character")

vecSampleNames <- apply(dfCellsHSMM, 1, function(sample){ unlist(strsplit(sample,split="/"))[2] })
vecSampleNames <- unlist(lapply(vecSampleNames, function(sample){ unlist(strsplit(sample,split="_1"))[1] }))
colnames(dfCountsHSMM) <- vecSampleNames
rownames(dfCountsHSMM) <- dfGenesHSMM[,1]
colnames(dfFpkmHSMM) <- vecSampleNames
rownames(dfFpkmHSMM) <- dfGenesHSMM[,1]

# get single cell samples
fileSampleAnnot <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/SraRunTable_scRNAMonocle.txt"
dfSampleAnnot <- read.table(fileSampleAnnot, header=TRUE, sep="\t", colClasses="character")
vecSCSamples <- dfSampleAnnot[dfSampleAnnot$library_protocol_s == "Single-cell",]$Run_s
dfCountsHSMM_SC <- dfCountsHSMM[,vecSCSamples]
dfFpkmHSMM_SC <- dfFpkmHSMM[,vecSCSamples]

# Make monolcle annotation
HSMM_gene_annotationRAW <- data.frame( gene_short_name=dfGenesHSMM[,2], biotype="N/A" )
rownames(HSMM_gene_annotationRAW) <- dfGenesHSMM[,1]
HSMM_sample_sheetRAW <- data.frame( Library=dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$SRA_Study_s, 
  Media="GM",
  Hours=0,
  stringsAsFactors=FALSE)
rownames(HSMM_sample_sheetRAW) <- dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$Run_s
# Differentiation medium DM from growth medium GM after t=0h
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T0",]$Media <- "DM"
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T24",]$Hours <- 24
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T48",]$Hours <- 48
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T72",]$Hours <- 72

# Get pseudotime
#source("http://bioconductor.org/biocLite.R")
#biocLite("monocle")
#biocLite(c("Biobase"))
#devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
#source("https://bioconductor.org/biocLite.R")

detach("package:monocle", unload=TRUE)
library(monocle)
library(Biobase)

pd <- new("AnnotatedDataFrame", data=HSMM_sample_sheetRAW)
fd <- new("AnnotatedDataFrame", data=HSMM_gene_annotationRAW)
HSMM <- newCellDataSet( cellData=as.matrix(dfFpkmHSMM_SC), phenoData=pd, featureData=fd )
print(head(fData(HSMM)))

# filter
HSMM <- detectGenes(HSMM, min_expr=0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr="expression~Media",cores=NCORES)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes <- intersect(ordering_genes, expressed_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, use_irlba=F)
# num_paths is number of leaves of tree
HSMM <- orderCells(HSMM, num_paths=2, reverse=T)
print(pData(HSMM))
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/MST_state.pdf")
plot_spanning_tree(HSMM)
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/MST_time.pdf")
plot_spanning_tree(HSMM, color_by ="Hours")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/MST_pseudotime.pdf")
plot_spanning_tree(HSMM, color_by ="Pseudotime")
dev.off()
plot_spanning_tree(HSMM, color_by ="State")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/data")
save(HSMM,file=file.path(getwd(),"PseudoDE_HSMM.RData"))
save(HSMM_sample_sheetRAW,file=file.path(getwd(),"PseudoDE_HSMM_sample_sheetRAW.RData"))
save(HSMM_gene_annotationRAW,file=file.path(getwd(),"PseudoDE_HSMM_gene_annotationRAW.RData"))
save(dfCountsHSMM_SC,file=file.path(getwd(),"PseudoDE_dfCountsHSMM_SC.RData"))
save(dfFpkmHSMM_SC,file=file.path(getwd(),"PseudoDE_dfFpkmHSMM_SC.RData"))

matCountsAll <- dfCountsHSMM_SC
dfAnnotationAllMonocle <- pData(HSMM)
HSMM_filtered <- HSMM[expressed_genes, pData(HSMM)$State != 3]
matCountsNoState3 <- dfCountsHSMM_SC[,rownames(pData(HSMM[,pData(HSMM)$State != 3]))]
dfAnnotationNoState3 <- pData(HSMM_filtered)

write.table(matCountsAll, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsAll.tab")
write.table(as.vector(rownames(matCountsAll)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsAll_genes.tab")
write.table(as.vector(colnames(matCountsAll)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsAll_samples.tab")
write.table(dfAnnotationAllMonocle, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_dfAnnotationAllMonocle.tab")

write.table(matCountsNoState3, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsNoState3.tab")
write.table(as.vector(rownames(matCountsNoState3)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsNoState3_genes.tab")
write.table(as.vector(colnames(matCountsNoState3)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_matCountsNoState3_samples.tab")
write.table(dfAnnotationNoState3, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/inputLineagePulse/HSMM_dfAnnotationNoState3.tab")

if(FALSE){
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/data")
  load("PseudoDE_HSMM.RData")
  load("PseudoDE_HSMM_sample_sheetRAW.RData")
  load("PseudoDE_HSMM_gene_annotationRAW.RData")
  load("PseudoDE_dfCountsHSMM_SC.RData")
  load("PseudoDE_dfFpkmHSMM_SC.RData")
}

# Find contaminating
HSMM_filtered <- HSMM[expressed_genes,]
cds_subset <- HSMM_filtered[row.names(subset(fData(HSMM_filtered), gene_short_name %in% c("PDGFRA", "SPHK1", "CDK1", "MEF2C", "MYH3"))),]
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/MarkerGenesInPT_pseudotime.pdf")
plot_genes_in_pseudotime(cds_subset, color_by="Hours")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/MarkerGenesInPT_state.pdf")
plot_genes_in_pseudotime(cds_subset, color_by="State")
dev.off()
plot_genes_in_pseudotime(cds_subset, color_by="State")
# No clear separation into contamined and not contaminated?

# DE analysis all cells
HSMM_filtered <- HSMM[expressed_genes,]

full_model_fits <- fitModel(HSMM_filtered,  modelFormulaStr="expression~VGAM::bs(Pseudotime)", cores=3)
reduced_model_fits <- fitModel(HSMM_filtered, modelFormulaStr="expression~1", cores=3)
pseudotime_test_res <- compareModels(full_model_fits, reduced_model_fits)
# Merge the test results with the metadata for each gene, including HUGO symbol, etc.
pseudotime_test_res <- merge(fData(HSMM), pseudotime_test_res, by="row.names")

lsDEAnalysisAllCells <- list(HSMM_filtered=HSMM_filtered,
  full_model_fits=full_model_fits,
  reduced_model_fits=reduced_model_fits,
  pseudotime_test_res=pseudotime_test_res)
save(lsDEAnalysisAllCells, "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/data/lsDEAnalysisAllCells.RData")

# DE analysis without state 3
HSMM_filtered <- HSMM[expressed_genes, pData(HSMM)$State != 3]

full_model_fits <- fitModel(HSMM_filtered,  modelFormulaStr="expression~VGAM::bs(Pseudotime)", cores=3)
reduced_model_fits <- fitModel(HSMM_filtered, modelFormulaStr="expression~1", cores=3)
pseudotime_test_res <- compareModels(full_model_fits, reduced_model_fits)
# Merge the test results with the metadata for each gene, including HUGO symbol, etc.
pseudotime_test_res <- merge(fData(HSMM), pseudotime_test_res, by="row.names")

lsDEAnalysisNoState3 <- list(HSMM_filtered=HSMM_filtered,
  full_model_fits=full_model_fits,
  reduced_model_fits=reduced_model_fits,
  pseudotime_test_res=pseudotime_test_res)
save(lsDEAnalysisNoState3, "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/monocle/data/lsDEAnalysisNoState3.RData")

if(FALSE){
  setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
  load("PseudoDE_HSMM.RData")
  load("PseudoDE_HSMM_sample_sheetRAW.RData")
  load("PseudoDE_HSMM_gene_annotationRAW.RData")
  load("PseudoDE_dfCountsHSMM_SC.RData")
  load("PseudoDE_dfFpkmHSMM_SC.RData")
  
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_dfAnnotation.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsProc.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matProbNB.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matMuCluster.RData")
}

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

dfPseudotime <- data.frame( pseudotime=as.vector(vecPT) )
plotEDF <- ggplot() +
  geom_density(data=dfPseudotime, aes(x=pseudotime), colour="black", bw=1) +
  labs(title="Density estimation of cells in pseudotime") +
  xlab("pseudotime") +
  ylab("empirical probability density")
print(plotEDF)

matCountsRd <- matCounts[apply(matCounts,1,function(gene){any(gene>10) & mean(gene)>=1}),]
matCountsRd <- round(matCountsRd)
save(matCountsRd,file=file.path(getwd(),"PseudoDE_matCountsRd.RData"))
save(vecPT,file=file.path(getwd(),"PseudoDE_vecPT.RData"))

nProc=3
source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
lsDEresults <- runPseudoDE(matCounts=matCountsRd,
  vecPseudotime=vecPT,
  K=6,
  scaSmallRun=20,
  boolPseudotime = TRUE,
  boolContPseudotimeFit=TRUE,
  boolPlotZINBfits=FALSE,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  nProc=nProc,
  scaMaxiterEM=20)

dfDEImpulse <- data.frame( lsDEresults$lsImpulseDE2results$dfImpulseResults[c("Gene","adj.p")], stringsAsFactors = FALSE)
dfDEModelfree <- data.frame( lsDEresults$dfModelFreeDEAnalysis[c("Gene","adj.p")], stringsAsFactors = FALSE)
dfDEModelfree <- dfDEModelfree[match(dfDEModelfree$Gene,dfDEImpulse$Gene),]
dfComparison <- cbind(dfDEImpulse,dfDEModelfree)
dfComparison

if(FALSE){
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_dfAnnotation.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matCountsProc.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matCountsProcFull.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matProbNB.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matDropout.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_vecDispersions.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_vecSizeFactors.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matMuCluster.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out/PseudoDE_matMu.RData")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
  load("ImpulseDE2_matCountDataProc.RData")
  load("ImpulseDE2_dfAnnotationProc.RData")
  # Load Impulse output
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_vecDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_dfDESeq2Results.RData")
  load("ImpulseDE2_lsMatTranslationFactors.RData")
  load("ImpulseDE2_matSizeFactors.RData")
  source("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/building/code_files/PseudoDE_main.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
  graphics.off()
  plotDEGenes(vecGeneIDs=vecDEGenes, 
    matCountDataProc=matCountDataProc, 
    lsMatTranslationFactors=lsMatTranslationFactors, 
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits, 
    matMuCluster=matMuCluster,
    vecClusterAssignments=lsResultsClustering$Assignments,
    dfImpulseResults=dfImpulseResults, 
    vecRefPval=NULL, 
    strCaseName="case", 
    strControlName=NULL, 
    strMode="singlecell",
    strSCMode="continuous",
    strFileNameSuffix = "DEgenes", 
    strPlotTitleSuffix = "", 
    strPlotSubtitle = "",
    strNameMethod1="ImpulseDE2", 
    strNameMethod2=NULL,
    boolSimplePlot=FALSE, 
    boolLogPlot=TRUE,
    NPARAM=6)
}

if(FALSE){
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_dfAnnotation.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_vecDispersions.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsInputToImpulseDE2.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsProc.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matProbNB.RData")
  matDropout <- lsInputToImpulseDE2$matDropout
  matProbNB <- lsInputToImpulseDE2$matProbNB
  matCountsImputed <- lsInputToImpulseDE2$matCountsImputed
}
if(FALSE){
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_coefs_mu.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_fit_mu.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_Y.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_zhat.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_thetahat.RData")
  #print(coefs_mu)
  print(fit_mu)
  print(thetahat)
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_fit_pi.RData")
  print(fit_pi)
}
