/* General utility functions */

// print timestamp
void timestamp();

// calculate scalar product of arrays a, b, both with length n
double scalarProduct(double* a, double* b, int n);

// calculate Matrix-Vector multiplication A*b (must be both double!)
void multMatrixVector(double* A, double* b, double* res, int nrow, int ncol);

// sampling without replacement for random update
void sample(int n, int *res);

// check upper quantile
void checkUpper(double val, double *upper, int len);

// check lower quantile
void checkLower(double val, double *lower, int len);