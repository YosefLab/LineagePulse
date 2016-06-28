#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
### redirect stdout/stderr 
#PBS -e localhost:/data/yosef2/users/fischerd/code/scriptreports/160627_scRNAseqLungDEAnalysis_NoState2.err
#PBS -o localhost:/data/yosef2/users/fischerd/code/scriptreports/160627_scRNAseqLungDEAnalysis_NoState2.out

#PBS -m ae
#PBS -M fischerd@berkeley.edu

### ############################################################################
Rscript /data/yosef2/users/fischerd/data/scRNAseq_Lung/LineagePulse/runLineagePulse_Lung_NoState2_cluster.R