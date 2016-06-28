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

# Make monolcle annotation
Lung_gene_annotationRAW <- data.frame( gene_short_name=dfGenesLung[,2], biotype="N/A" )
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

# Get pseudotime
# monocle
#biocLite("monocle")
#source("http://bioconductor.org/biocLite.R")
# monocle2
biocLite(c("Biobase"))
devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
source("https://bioconductor.org/biocLite.R")

#detach("package:monocle", unload=TRUE)
library(monocle)
#library(data.table)
library(Biobase)


pd <- new("AnnotatedDataFrame", data=Lung_sample_sheetRAW)
fd <- new("AnnotatedDataFrame", data=Lung_gene_annotationRAW)
Lung <- newCellDataSet( cellData=as(as.matrix(dfFpkmLung_SC), "sparseMatrix"), 
  phenoData=pd, 
  featureData=fd,
  expressionFamily=negbinomial.size())
Lung <- estimateSizeFactors(Lung,method="blind")

# filter
Lung <- detectGenes(Lung, min_expr=0.1)
print(head(fData(Lung)))
print(head(pData(Lung)))
expressed_genes <- row.names(subset(fData(Lung), num_cells_expressed >= 10))

pData(Lung)$Total_mRNAs <- Matrix::colSums(exprs(Lung))
Lung <- Lung[,row.names(subset(pData(Lung), Total_mRNAs >= 10000 & Total_mRNAs <= 40000))]
Lung <- detectGenes(Lung, min_expr = 0.1)

diff_test_res <- differentialGeneTest(Lung[expressed_genes,], fullModelFormulaStr="expression~Hours",cores=NCORES)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes <- intersect(ordering_genes, expressed_genes)
Lung <- setOrderingFilter(Lung, ordering_genes)
Lung <- reduceDimension(Lung)
# num_paths is number of leaves of tree
Lung <- orderCells(Lung, num_paths=2, reverse=F)

plot_cell_trajectory(Lung, color_by="Hours")



plot_spanning_tree(Lung)
print(pData(Lung))

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
save(Lung,file=file.path(getwd(),"PseudoDE_Lung.RData"))
save(Lung_sample_sheetRAW,file=file.path(getwd(),"PseudoDE_Lung_sample_sheetRAW.RData"))
save(Lung_gene_annotationRAW,file=file.path(getwd(),"PseudoDE_Lung_gene_annotationRAW.RData"))
save(dfCountsLung_SC,file=file.path(getwd(),"PseudoDE_dfCountsLung_SC.RData"))
save(dfFpkmLung_SC,file=file.path(getwd(),"PseudoDE_dfFpkmLung_SC.RData"))
if(FALSE){
  setwd("/Users/davidsebastianfischer/MasterThesis/code/LineagePulse/software_test_out")
  load("PseudoDE_Lung.RData")
  load("PseudoDE_Lung_sample_sheetRAW.RData")
  load("PseudoDE_Lung_gene_annotationRAW.RData")
  load("PseudoDE_dfCountsLung_SC.RData")
  load("PseudoDE_dfFpkmLung_SC.RData")
  
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_dfAnnotation.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsProc.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matProbNB.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matMuCluster.RData")
}

if(TRUE){
  # Investigate distribution of cells over pseudotime
  vecPTpointsAll_Lung <- as.vector(pData(Lung)$Pseudotime)
  names(vecPTpointsAll_Lung) <- as.vector(rownames(pData(Lung)))
  vecPTpoints_Lung <- vecPTpointsAll_Lung[!is.na(vecPTpointsAll_Lung)]
  print(paste0("Lung: Total cells: ",length(vecPTpointsAll_Lung),", Non NA: ",length(vecPTpoints_Lung)))
}else{
  # Use real observation time
  vecPTpointsAll_Lung <- Lung_sample_sheetRAW$Hours
  names(vecPTpointsAll_Lung) <- rownames(Lung_sample_sheetRAW)
  vecPTpoints_Lung <- vecPTpointsAll_Lung[!is.na(vecPTpointsAll_Lung)]
}

# Run PseudoDE
matCounts <- data.matrix(dfCountsLung_SC)
vecPT <- vecPTpoints_Lung
matCounts <- matCounts[,colnames(matCounts) %in% names(vecPTpoints_Lung)]
vecPT <- vecPT[colnames(matCounts)]
