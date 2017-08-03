# Code Supplement

This repository contains the code supplement for the paper

### The Impact of Model Assumptions in Scalar-on-Image Regression     
*Clara Happ, Sonja Greven, Volker Schmid (Department of Statistics, LMU Munich, Munich, Germany)*     
*for the Alzheimer's Disease Neuroimaging Initiative.*  
A preprint is available on [arXiv](https://arxiv.org/abs/1707.02233).

It provides:  
 *  Usage examples for all models used in the paper       
 * `R` implementations of methods, if not already available  
 * `R` functions for all measures developed in the paper  
  
## Usage ##

`R` functions are directly applicable. The `C` implementation of the Bayesian GMRF models requires compilation. Change to the `C` subdirectory and run the following code in the command line (tested under Linux only)

1. `R CMD SHLIB utilities/*.c` (compiles all utility functions)  
2. `R CMD SHLIB mainGibbs_GMRF.c utilities/*.o` (compiles main for GMRF)  
3. `R CMD SHLIB mainGibbs_HyperparamsFixed.c utilities/*.o` (compiles main for SparseGMRF) 

Make sure that the `Makevars` file is in the same directory as the main files. 


## Bug reports ##

Please use [GitHub issues](https://github.com/ClaraHapp/SOIR/issues) for reporting bugs or issues.

