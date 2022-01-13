# ImageMonitoring_SparseLearning
Here are the computer codes for my manuscript ``Statistical Quality Control Using Image Intelligence: A Sparse Learning Approach". 

Most heavy lifting is done in imspc.c, which is called in main.R and corr.R. To use the code, do the following first on the command line (then you can run R files):
R CMD SHLIB imspc.c

main.R: includes the numerical examples where image noise are i.i.d.

corr.R: includes the numerical examples where image noise are correlated.
