#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
### redirect stdout/stderr 
#PBS -e localhost:/data/yosef2/users/fischerd/code/scriptreports/160622_scRNAseqHSMMComAnalysis.err
#PBS -o localhost:/data/yosef2/users/fischerd/code/scriptreports/160622_scRNAseqHSMMComAnalysis.out

#PBS -m ae
#PBS -M fischerd@berkeley.edu

### ############################################################################
Rscript /data/yosef2/users/fischerd/data/scRNAseq_Monocle/LineagePulse/runLineagePulse_HSMM_contPT_Impulse.R