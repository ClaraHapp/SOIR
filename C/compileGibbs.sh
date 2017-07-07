R CMD SHLIB utilities/*.c #compile all utility functions
R CMD SHLIB  mainGibbs.c utilities/*.o # compile main
