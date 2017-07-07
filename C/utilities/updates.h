/* Update functions used in Gibbs Sampler */

// update variance parameter sigma_epsilon
void updateSigEpsInformative(double* sigEps, double shapeEps, double scaleEps, int nPat, double* Y, double* FixEfAlpha, double* xbeta);


// update variance parameter sigma_beta without Ising field
void updateSigBetaInformative_betaOnly(double* sigBeta, double shapeBeta, double scaleBeta, int nVox, double* beta, int* numNeigh, double* sumPredNeigh);


// update coefficients for intercept and fixed effects
void updateAlpha(double* alpha, int nPat, double* Y, double* xbeta, int nFixEf, double* FixEf, double* FixEfInvRoot, double sigEps);


// update beta and gamma via single-site Gibbs sampler; based on Goldsmith et al., 2014
void updateBetaGamma(int nPat, int nVox, int nDelta, //array sizes
                 int l, // actual voxel index
                 double* Xl, double* Y, // data
                 double* beta, int* gamma, // parameters
                 double* InProd, int* delta, double* xbeta0, double* FixEfAlpha, int* numNeigh, double* sumPredNeigh, //auxiliary variables 
                 double sigEps, double sigBeta, double a, double b); //hyperparameters


// update beta via single-site Gibbs sampler
void updateBetaOnly(int nPat, int nVox, int nDelta, //array sizes
                 int l, // actual voxel index
                 double* Xl, double* Y, // data
                 double* beta, // parameters
                 double* InProd, int* delta, double* xbeta0, double* FixEfAlpha, int* numNeigh, double* sumPredNeigh, //auxiliary variables 
                 double sigEps, double sigBeta); //hyperparameters