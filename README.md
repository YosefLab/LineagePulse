## LineagePulse ##
### Differential expression analysis and model fitting for single-cell RNA-seq data in pseudotime ###

LineagePulse is a differential expression algorithm for single-cell RNA-seq (scRNA-seq) data.
Its main application is differential expression analysis across pseudotemporal orderings (developmental trajectories) of scRNA-seq data.
LineagePulse can also be used for group-wise differential expression tests.
As LineagePulse is based on model fitting, it can also be used for batch correction and gene expression trajectory estimation only.

### INSTALLATION

Clone the GitHub LineagePulse repository first into you local target directory
and then install the package via "R CMD INSTALL .".

### VIGNETTE
View the vignette (.html) in your browser by double clicking on "vignettes/LineagePulse_Tutorial.Rmd".

### Problems with BiocParallel:
LineagePulse is parallelised with BiocParallel. 
If BiocParallel does not start any processes during the model estimation, add 
options(bphost="localhost")
into your R script before the first call to any LineagePulse function.