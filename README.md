# LineagePulse

# Problems with BiocParallel:
LineagePulse is parallelised with BiocParallel. 
If BiocParallel does not start any processes during the model estimation, add 
options(bphost="localhost")
into your R script before the first call to any LineagePulse function.