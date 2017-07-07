/* main Gibbs sampling algorithm for hierarchical scalar-on-image regression model, based on Goldsmith et al, 2014 */

// include libraries
#include <stdio.h>
#include <time.h> // for timestamp
#include <float.h> // for DBL_MAX
#include <R.h> 
#include <Rmath.h>
#include <math.h> // for pi = M_PI

// include utility functions
#include "updates.h" // Update functions used in Gibbs Sampler
#include "utilities.h" // General utility functions

//macros
#define matrix(array, nRow,i,j) array[j*(nRow)+i] // for simpler matrix access


//main Gibbs sampling algorithm for hierarchical scalar-on-image regression model, to be called from R
/*
 * @args:	K: number of iterations
 * 			nPat: number of patients
 * 			nVox: numer of voxels in (vectorized) images
 *			nFixEf: number of fixed effects (incl. intercept)
 *			nDelta: maximal number of neighbours in neighbourhood structure
 * 			nMCMC: thinning interval
 *			burnin: the length of the burnin period
 * 			X: matrix of (vectorized) images, dim = nVox*nPat
 * 			FixEf: matrix of fixed effects (incl. intercept), dim = nPat * nFixEf
 *			Y: response vector of length nPat
 * 			delta: matrix of neighbourhood structure, dim = nVox * nDelta
 * 			FixEfInvRoot: auxiliary variable, matrix root of inverse of t(FixEf)%*%FixEf, dim = nFixEf * nFixEf
 *			sigEps, sigBeta: variance hyperparameters
 *			a, b: Ising parameters a and b
 * 			alphaRes, betaRes, gammaRes: empty verctors for results
 * 			betaLower, betaUpper: empty vectors for pointwise quantile-based CI
 */ 
void allGibbs(int* K, int* nPat, int* nVox, int* nFixEf, int* nDelta, int* nMCMC, int* burnin, // array sizes
              double* X, double* FixEf, double* Y, // data
              int* delta, // neighbourhood structure
              double* FixEfInvRoot, // auxiliary variable, needed for alpha generation
              double* sigEps, double* sigBeta, // variance hyperparameters
              double* a, double* b, //hyperparameters for Ising field
              int* calcCI,  // should pointwise quantile CI be calculated (2.5 and 97.5 quantiles)
              double* alphaRes, double* betaRes, double* gammaRes, 
              double* betaLower, double* betaUpper, double* alphaLower, double* alphaUpper) //files for results 
{
  
  /* auxiliary variables */
  
  // for basic code (iteration variables etc)
  unsigned int i, j, l, k, writeEach;
  double sum;
  
  // model parameters
  double alpha[*nFixEf];
  double beta[*nVox];
  int gamma[*nVox];
  
  // auxiliary arrays
  double Xl[*nPat];
  double InProd[*nVox];
  double xbeta[*nPat];
  double xbeta0[*nPat];
  double FixEfAlpha[*nPat];
  int numNeigh[*nVox];
  double sumPredNeigh[*nVox];
  int updateOrder[*nVox];
  
  // arrays for quantile calculation
  int nQuant;
  double** lowerTailBeta;
  double** upperTailBeta;
  double** lowerTailAlpha;
  double** upperTailAlpha;
  
  // allocate memory if CIs are really needed
  if(*calcCI != 0)
  {
    nQuant = (int)(0.025* (*K - *burnin) / *nMCMC);
    
    if(nQuant == 0)
      error("Number of iterations is too low for quantile calculation!");
    
    lowerTailBeta = calloc(*nVox, sizeof(double *));
    for (i=0; i<*nVox; i++)
      lowerTailBeta[i] = calloc(nQuant, sizeof(double));
    
    upperTailBeta = calloc(*nVox, sizeof(double *));
    for (i=0; i<*nVox; i++)
      upperTailBeta[i] = calloc(nQuant, sizeof(double));
    
    lowerTailAlpha = calloc(*nFixEf, sizeof(double *));
    for (i=0; i<*nFixEf; i++)
      lowerTailAlpha[i] = calloc(nQuant, sizeof(double));
    
    upperTailAlpha = calloc(*nFixEf, sizeof(double *));
    for (i=0; i<*nFixEf; i++)
      upperTailAlpha[i] = calloc(nQuant, sizeof(double));
  }
  
  
  // results
  int countRes;
  
  /* for random number generation */
  GetRNGstate(); // for random number generation
  
  
  /* initialize variables */
  
  countRes = 0;
  
  
  // parameters: alpha
  for(i = 0; i < *nFixEf; i++)
  {
    alpha[i] = 0;
    alphaRes[i] = 0;
  }
  
  // parameters: beta & gamma
  for(i = 0; i < *nVox; i++)
  {
    if(runif(0,1) < 0.5) // with probability 0.5: initialize voxel with 1 / N(0, sigBeta)
    {
      gamma[i] = 1;
      beta[i] = rnorm(0, sqrt(*sigBeta));
    }
    else // with probability 0.5: initialize voxel with 0 / 0
    {
      gamma[i] = 0;
      beta[i] = 0;
    }
    
    betaRes[i] = 0;
    gammaRes[i] = 0;
    
  }
  
  // auxiliary arrays
  multMatrixVector(X, beta, xbeta, *nVox, *nPat); // xBeta = X * beta (should be = 0 if beta = 0)
  multMatrixVector(FixEf, alpha, FixEfAlpha, *nPat, *nFixEf); // FixEfAlpha = FixEf * alpha
  
  for(i = 0; i < *nVox; i++) // vector of inner products, sum(X[i,]^2) 
  {
    InProd[i] = 0;
    
    for(l = 0; l < *nPat; l++)
      InProd[i] = InProd[i] + matrix(X, *nVox,i,l)* matrix(X, *nVox,i,l);
  }
  
  
  for(l = 0; l < *nVox; l++)
  {
    numNeigh[l] =0;
    sumPredNeigh[l] = 0;
    
    for(i = 0; i < *nDelta; i++)
    {
      if(matrix(delta, *nVox, l, i) !=  NA_INTEGER) // matrix(delta, nVox, l, i)-1 is neighbour of l
      {
        numNeigh[l]++;
        
        sumPredNeigh[l] = sumPredNeigh[l] + beta[matrix(delta, *nVox, l, i)-1] * gamma[matrix(delta, *nVox, l, i)-1];
      }
    }	
  }	
  
  // quantile auxiliary arrays
  if(*calcCI != 0)
  {
    for(i = 0; i < *nVox; i++)
    {
      betaUpper[i] = 0;
      betaLower[i] = 0;
      
      for(j = 0; j < nQuant; j++)
      {
        lowerTailBeta[i][j] = DBL_MAX;
        upperTailBeta[i][j] = -DBL_MAX;
      } 
    }
    
    for(i = 0; i < *nFixEf; i++)
    {
      alphaUpper[i] = 0;
      alphaLower[i] = 0;
      
      for(j = 0; j < nQuant; j++)
      {
        lowerTailAlpha[i][j] = DBL_MAX;
        upperTailAlpha[i][j] = -DBL_MAX;
      } 
    }
  }
  
  /* Gibbs algorithm */
  
  for(k = 0; k < (*K); k++) // K Gibbs steps
  {
    /*if(k%1000 == 0)
     {
     timestamp();
     printf("k = %04d (%.2f %% of total)\n", k, (double)k/(double)*K);
     }*/
    
    sample(*nVox, updateOrder); // sample update order randomly
    
    // update beta & gamma: sweep over all voxels
    for(j = 0; j < (*nVox); j++) 
    { 		
      /* compute auxiliary variables for local voxel l */
      
      // get index of next voxel
      l = updateOrder[j];
      
      //Xl = X[l,], xbeta0 = xbeta - beta[l]*Xl 
      for(i = 0; i < (*nPat); i++)
      {
        Xl[i] = matrix(X, *nVox,l,i);
        
        xbeta0[i] = xbeta[i] - beta[l]*Xl[i];
      }
      
      /* do local iteration */
      updateBetaGamma(*nPat, *nVox, *nDelta, l, Xl, Y, beta, gamma, InProd, delta, xbeta0, FixEfAlpha, numNeigh, sumPredNeigh, *sigEps, *sigBeta, *a, *b);
      
      /* update (global!) X * Beta based on current iteration */
      
      for(i = 0; i < (*nPat); i++)
        xbeta[i] = xbeta0[i] + beta[l]*Xl[i];
    }
    
    
    /* update alpha (impose flat prior) and FixEfAlpha */
    
    updateAlpha(alpha, *nPat, Y, xbeta, *nFixEf, FixEf, FixEfInvRoot, *sigEps);
    
    // FixEfAlpha = FixEf * alpha
    multMatrixVector(FixEf, alpha, FixEfAlpha, *nPat, *nFixEf);
    
    
    /* update variables for convergence analysis */
    if(k >= *burnin)
    {
      // norm of gamma
      sum = 0;
      for(i = 0; i < *nVox; i++) 
        sum = sum + (double)gamma[i]; // gamma[i]^2 = gamma[i]
      
      
      if(k % (*nMCMC) == 0)
      {
        countRes++;
        
        
        for(i = 0; i < *nFixEf; i++)
          alphaRes[i] = alphaRes[i] + alpha[i];
        
        for(i = 0; i < *nVox; i++)
        {
          betaRes[i] = betaRes[i] + beta[i];
          gammaRes[i] = gammaRes[i] + gamma[i];
        }
        
        if(*calcCI != 0) // quantiles
        {
          for(i = 0; i < *nVox; i++)
          {
            checkUpper(beta[i], upperTailBeta[i], nQuant);
            checkLower(beta[i], lowerTailBeta[i], nQuant);
          }
          
          for(i = 0; i < *nFixEf; i++)
          {
            checkUpper(alpha[i], upperTailAlpha[i], nQuant);
            checkLower(alpha[i], lowerTailAlpha[i], nQuant);
          }
        }
      }
    }
  } 
  
  
  for(i = 0; i < *nFixEf; i++)
    alphaRes[i] = alphaRes[i]/countRes;
  
  for(i = 0; i < *nVox; i++)
  {
    betaRes[i] = betaRes[i] /countRes;
    gammaRes[i] = (double) gammaRes[i] /countRes;
  }
  
  
  if(*calcCI != 0)
  {
    
    for(i = 0; i < *nVox; i++)
    {
      betaUpper[i] = upperTailBeta[i][0];
      betaLower[i] = lowerTailBeta[i][nQuant - 1];
    }
    
    free(upperTailBeta);
    free(lowerTailBeta);
    
    
    for(i = 0; i < *nFixEf; i++)
    {
      alphaUpper[i] = upperTailAlpha[i][0];
      alphaLower[i] = lowerTailAlpha[i][nQuant - 1];
    }
    
    free(upperTailAlpha);
    free(lowerTailAlpha);
  }
  
  PutRNGstate(); // finish random number generation
  
  return;
}
