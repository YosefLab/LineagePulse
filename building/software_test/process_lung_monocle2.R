detach("package:monocle", unload=TRUE)
library(Biobase)
devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
library(monocle)

rm(list = ls())
NCORES = 3

# Load data
dfCountsLung <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/rsem/rsem_readCountsTable.txt", header=FALSE, colClasses="numeric")
dfFpkmLung <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/rsem/rsem_fpkmTable.txt", header=FALSE, colClasses="numeric")
dfCellsLung <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/rsem/cell_list.txt", header=FALSE, colClasses="character")
dfGenesLung <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/rsem/gene_list.txt", header=FALSE, colClasses="character")

vecSampleNames <- apply(dfCellsLung, 1, function(sample){ unlist(strsplit(sample,split="/"))[2] })
vecSampleNames <- unlist(lapply(vecSampleNames, function(sample){ unlist(strsplit(sample,split="_1"))[1] }))
colnames(dfCountsLung) <- vecSampleNames
rownames(dfCountsLung) <- dfGenesLung[,1]
colnames(dfFpkmLung) <- vecSampleNames
rownames(dfFpkmLung) <- dfGenesLung[,1]

# get single cell samples - all are sc!
fileSampleAnnot <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/SraRunTable_Lung.txt"
dfSampleAnnot <- read.table(fileSampleAnnot, header=TRUE, sep="\t", colClasses="character")
vecSCSamples <- dfSampleAnnot$Run_s
dfCountsLung_SC <- dfCountsLung[,vecSCSamples]
dfFpkmLung_SC <- dfFpkmLung[,vecSCSamples]

# Make monocle annotation
Lung_gene_annotationRAW <- data.frame( gene_short_name=dfGenesLung[,2], Media="M" )
rownames(Lung_gene_annotationRAW) <- dfGenesLung[,1]
Lung_sample_sheetRAW <- data.frame( Library=dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$SRA_Study_s, 
  Hours=0,
  stringsAsFactors=FALSE)
rownames(Lung_sample_sheetRAW) <- dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$Run_s
# Differentiation medium DM from growth medium GM after t=0h
Lung_sample_sheetRAW[,"Hours"] <- sapply(as.vector(dfSampleAnnot$age_s),function(cell){
  age <- unlist(strsplit(cell, split=" "))
  return(as.numeric(age[length(age)]))
})

# Filter cells and genes -> otherwise estimateDispersions fails!
matCountsLung <- as.matrix(round(dfCountsLung_SC))
vecTotalCounts <- apply(matCountsLung,2,function(cell){sum(cell,na.rm=TRUE)})
vecboolEmptyCells <- apply(matCountsLung,2,function(cell){!all(cell[!is.na(cell)]==0)})
vecboolEmptyGenes <- apply(matCountsLung,1,function(gene){!all(gene[!is.na(gene)]==0)})
matCountsLung <- matCountsLung[vecboolEmptyGenes,vecboolEmptyCells]
Lung_gene_annotation <- Lung_gene_annotationRAW[rownames(matCountsLung),]
Lung_sample_sheet <- Lung_sample_sheetRAW[colnames(matCountsLung),]
pd <- new("AnnotatedDataFrame", data=Lung_sample_sheet)
fd <- new("AnnotatedDataFrame", data=Lung_gene_annotation)
Lung <- newCellDataSet(as.matrix(matCountsLung),
  phenoData = pd,
  featureData = fd,
  expressionFamily=negbinomial() )
Lung <- estimateSizeFactors(Lung)
Lung <- estimateDispersions(Lung)

# This example performs a greatly simplified version of the single-cell RNA-Seq
# analysis of skeletal myoblast differentiation
# described in Trapnell, Cacchiarelli et al (Nature Biotechnology, 2014).

# Count how many cells each gene is expressed in, and how many
# genes are expressed in each cell
Lung <- detectGenes(Lung, min_expr = 0.1)

# Get a list of genes expressed in at least 50 cells
expressed_genes <- row.names(subset(fData(Lung), num_cells_expressed >= 50))

# Test the above genes for differential expression in response from switch from GM to DM
# Note: this step can take several hours on a single core, so you might want to parallelize it
# with the 'cores' argument 
diff_test_res <- differentialGeneTest(Lung[expressed_genes,], fullModelFormulaStr="~Hours", cores=3)

# Use the differentially expressed genes as the basis for ordering the cells
# by progress through differentiation

# First: collet the list of genes that are significantly differentially expressed (at FDR < 1%)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

#ordering_genes <- intersect(ordering_genes, row.names(subset(fData(Lung), biotype %in% c("protein_coding"))))
ordering_genes <- intersect(ordering_genes, row.names(subset(fData(Lung))))

# Second: mark those genes as the ones used for ordering
Lung <- setOrderingFilter(Lung, ordering_genes)

# Third: perform dimensionality reduction using ICA
Lung <- reduceDimension(Lung, reduction_method="ICA")

# Fourth: compute the minimum spanning tree in the reduced space and use it to order the cells.
# Note that we're allowing a branch with two outcomes in the biological process
Lung <- orderCells(Lung, num_paths=2)

lsMonocole2input <- list(Lung, expressed_genes, diff_test_res, ordering_genes)
save(lsMonocole2input, file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/data/lsMonocole2input.RData")

### SAVE PROGRESS UP TO HERE
#load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/data/lsMonocole2input_v1.RData")
Lung <- lsMonocole2input[[1]]
expressed_genes <- lsMonocole2input[[2]]
diff_test_res <- lsMonocole2input[[3]]
ordering_genes <- lsMonocole2input[[4]]

LungRev <- orderCells(Lung, num_paths=2, reverse=T)

plot_cell_trajectory(Lung, color_by="Hours")
plot_cell_trajectory(Lung, color_by="Pseudotime")
plot_cell_trajectory(Lung, color_by="State")
plot_cell_trajectory(LungRev, color_by="Pseudotime")
plot_cell_trajectory(LungRev, color_by="State")
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_pseudotime.pdf")
plot_cell_trajectory(Lung, color_by="Pseudotime")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_pseudotimeReverse.pdf")
plot_cell_trajectory(LungRev, color_by="Pseudotime")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_time.pdf")
plot_cell_trajectory(Lung, color_by="Hours")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_state.pdf")
plot_cell_trajectory(Lung, color_by="State")
dev.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_stateReverse.pdf")
plot_cell_trajectory(LungRev, color_by="State")
dev.off()

# State 1 in Lung is state 2 in LungRev
Lung_filtered1 <- LungRev[expressed_genes, pData(LungRev)$State != 2]
matCountsNoState1 <- matCountsLung[,rownames(pData(LungRev[,pData(LungRev)$State != 2]))]
dfAnnotationNoState1 <- pData(Lung_filtered1)
Lung_filtered2 <- Lung[expressed_genes, pData(Lung)$State != 2]
matCountsNoState2 <- matCountsLung[,rownames(pData(Lung[,pData(Lung)$State != 2]))]
dfAnnotationNoState2 <- pData(Lung_filtered2)
Lung_filtered3 <- Lung[expressed_genes, pData(Lung)$State != 3]
matCountsNoState3 <- matCountsLung[,rownames(pData(Lung[,pData(Lung)$State != 3]))]
dfAnnotationNoState3 <- pData(Lung_filtered3)

plot_cell_trajectory(Lung, color_by="Pseudotime") 
# Find pseudotime coordinate of branch opint visually
Lung_filtered1Cont <- LungRev[expressed_genes, pData(Lung)$Pseudotime > 13.7]
matCountsNoState1Cont <- matCountsLung[,rownames(pData(LungRev[,pData(Lung)$Pseudotime > 13.7]))]
dfAnnotationNoState1Cont <- pData(Lung_filtered1Cont)
plot_cell_trajectory(Lung_filtered1Cont, color_by="Pseudotime")

Lung_filtered3Rev <- LungRev[expressed_genes, pData(LungRev)$State != 3]
matCountsNoState3Rev <- matCountsLung[,rownames(pData(LungRev[,pData(LungRev)$State != 3]))]
dfAnnotationNoState3Rev <- pData(Lung_filtered3Rev)
plot_cell_trajectory(Lung_filtered3Rev, color_by="Pseudotime")

write.table(matCountsNoState1, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1.tab")
write.table(as.vector(rownames(matCountsNoState1)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1_genes.tab")
write.table(as.vector(colnames(matCountsNoState1)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1_samples.tab")
write.table(dfAnnotationNoState1, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_dfAnnotationNoState1.tab")

write.table(matCountsNoState2, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState2.tab")
write.table(as.vector(rownames(matCountsNoState2)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState2_genes.tab")
write.table(as.vector(colnames(matCountsNoState2)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState2_samples.tab")
write.table(dfAnnotationNoState2, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_dfAnnotationNoState2.tab")

write.table(matCountsNoState3, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3.tab")
write.table(as.vector(rownames(matCountsNoState3)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3_genes.tab")
write.table(as.vector(colnames(matCountsNoState3)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3_samples.tab")
write.table(dfAnnotationNoState3, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_dfAnnotationNoState3.tab")

write.table(matCountsNoState1Cont, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1Cont.tab")
write.table(as.vector(rownames(matCountsNoState1Cont)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1Cont_genes.tab")
write.table(as.vector(colnames(matCountsNoState1Cont)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState1Cont_samples.tab")
write.table(dfAnnotationNoState1Cont, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_dfAnnotationNoState1Cont.tab")

write.table(matCountsNoState3Rev, row.names = FALSE, col.names = FALSE, sep="\t",
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3Rev.tab")
write.table(as.vector(rownames(matCountsNoState3Rev)), sep="\t",  row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3Rev_genes.tab")
write.table(as.vector(colnames(matCountsNoState3Rev)), sep="\t", row.names = FALSE, col.names = FALSE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_matCountsNoState3Rev_samples.tab")
write.table(dfAnnotationNoState3Rev, sep="\t", row.names = TRUE, col.names = TRUE,
  file="/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/inputLineagePulse/Lung_dfAnnotationNoState3Rev.tab")

lung_genes <- row.names(subset(fData(Lung), gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
plot_genes_jitter(Lung[lung_genes,], grouping="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn","Cdk1"))),]
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Scgb1a1", "Foxj1", "Pdpn","Ager","Sftpc","Sftpb"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/Lungepithelium/monocle2/MST_monocle2_MarkerGenesInPT.pdf")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Scgb1a1"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Foxj1"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Pdpn"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Ager"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Sftpc"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
cds_subset <- Lung[row.names(subset(fData(Lung), gene_short_name %in% c("Sftpb"))),]
plot_genes_in_pseudotime(cds_subset, color_by="State")
dev.off()

BEAM_res <- BEAM(Lung, branch_point=1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4,
  cores = 1,
  use_gene_short_name = T,
  show_rownames = T)

# Let's remove all the cells that Monocle tagged as being interstitial 
# mesenchymal cells (e.g. fibroblast contamination) so they don't interfere 
# with downstream analysis
Lung_filtered <- Lung[expressed_genes, pData(Lung)$State != 3]
#Lung_filtered <- Lung

# Grab a few specific genes so we can plot their kinetics in pseudotime
cds_subset <- Lung_filtered[row.names(subset(fData(Lung_filtered), gene_short_name %in% c("Cdk1"))),]
pdf("muscle_marker_pseudotime_plot.pdf")
plot_genes_in_pseudotime(cds_subset, color_by="State")
dev.off()