/* Update functions used in Gibbs Sampler */
#include <stdio.h>
#include <time.h> // for timestamp
#include <R.h> 
#include <Rmath.h>
#include <math.h> // for pi = M_PI

#include "utilities.h"

//macros
#define matrix(array, nRow,i,j) array[j*(nRow)+i] // for simpler matrix access


// update variance parameter sigma_epsilon
void updateSigEpsInformative(double* sigEps, double shapeEps, double scaleEps, int nPat, double* Y, double* FixEfAlpha, double* xbeta)
{
	double sum;
	unsigned int i;

	sum = 0;

	for(i = 0; i < nPat; i++)
		sum = sum + pow(Y[i] - FixEfAlpha[i] - xbeta[i],2);

	*sigEps =  1/rgamma(shapeEps + nPat/2, 1/(scaleEps + sum/2)); 
  
	return;
}


// update variance parameter sigma_beta
void updateSigBetaInformative_betaOnly(double* sigBeta, double shapeBeta, double scaleBeta, int nVox, double* beta, int* numNeigh, double* sumPredNeigh)
{
	double sum;
	unsigned int l;

	sum = 0;

	for(l = 0; l < nVox; l++)
	{
		// non-empty neighbourhood
		if(numNeigh[l] != 0) 
			sum = sum + pow(beta[l] - sumPredNeigh[l]/(double)numNeigh[l],2)*numNeigh[l]/2;
	}

	*sigBeta = 1/rgamma(shapeBeta + nVox/2, 1/(scaleBeta + sum)); // was ...rgamma(nZ, sum) before 04.06.2014
 
	return;
}


// update coefficients for intercept and fixed effects
void updateAlpha(double* alpha, int nPat, double* Y, double* xbeta, int nFixEf, double* FixEf, double* FixEfInvRoot, double sigEps)
{
	unsigned int i, m, l;
	double meanAlpha, sum;

	double alphaTmp[nFixEf];

	// alphaTmp ~ N(1/sqrt(sigEps)*FixEfInvRoot*t(FixEf)*tmp, I)
	for(i = 0; i < nFixEf; i++) // incl. intercept!
	{
		meanAlpha = 0;

		for(m = 0; m < nFixEf; m++)
		{
			sum = 0;
			for(l = 0; l < nPat; l++)
				sum =  sum + matrix(FixEf, nPat, l,m)*(Y[l] - xbeta[l]);

			meanAlpha = meanAlpha + sum*matrix(FixEfInvRoot,nFixEf, i,m);
		}

		meanAlpha = 1/sqrt(sigEps)*meanAlpha;

		alphaTmp[i] = rnorm(meanAlpha, 1);  
	}


	for(i = 0; i < nFixEf; i++)
	{
		alpha[i] = 0;
    
		for(l = 0; l < nFixEf; l++)
			alpha[i] = alpha[i] + matrix(FixEfInvRoot, nFixEf, l, i)*alphaTmp[l];

		alpha[i] = sqrt(sigEps)*alpha[i];
	}

	return;
}


// update beta and gamma via single-site Gibbs sampler; based on Goldsmith et al., 2014
void updateBetaGamma(int nPat, int nVox, int nDelta, //array sizes
                 int l, // actual voxel index
                 double* Xl, double* Y, // data
                 double* beta, int* gamma, // parameters
                 double* InProd, int* delta, double* xbeta0, double* FixEfAlpha, int* numNeigh, double* sumPredNeigh, //auxiliary variables 
                 double sigEps, double sigBeta, double a, double b) //hyperparameters
{
	/* auxiliary variables */
	unsigned int i;
	int indNeigh;
	double sigEpsInv, sigBetaInv, invNeighSize, betaBar, betaStar, Sigl, Mul, tmp, allNeigh, g, ppos;


	double resid[nPat]; // = Y- FixEfAlpha- xbeta0

	/* set auxiliary variables */
	sigEpsInv = 1/sigEps;
	sigBetaInv = 1/sigBeta;

	/* check neighbourhood */

	// empty neighbourhood -> continue with next voxel		
	if(numNeigh[l] == 0) 
	{
		gamma[l] = 0;
		beta[l] = 0;

		return;
	} 


	/* generate nonzero coefficient */

	betaBar = sumPredNeigh[l]/ (double) numNeigh[l]; //mean of neighbouring voxels

	//check neighbourhood of neighbours of l
	invNeighSize = 0; // sum of inverse neighbourhood sizes
	tmp = 0;

	// generate proposal sample from normal distribution with mean=Mul, sigma = sqrt(Sigl) 	
	Sigl = 1/(sigEpsInv * InProd[l] + sigBetaInv * numNeigh[l] );

	for(i = 0; i < nPat; i++)
		resid[i] = Y[i] - FixEfAlpha[i] - xbeta0[i];

	Mul = Sigl * ( sigEpsInv * scalarProduct(resid, Xl, nPat) + sigBetaInv * numNeigh[l] *  betaBar);

	betaStar = rnorm(Mul, sqrt(Sigl));

	/* compute g_l and posterior probability */

	tmp = 0;
	for(i = 0; i < nDelta; i++)
	{
		if( matrix(delta, nVox, l, i) !=  NA_INTEGER) //if delta[l,i] != NA
			tmp = tmp + (1-2* gamma[ matrix(delta, nVox, l, i) - 1]); // correct for indices starting with 0 (C) vs. indices starting with 1 (R)
	}

	g = sqrt(2*M_PI* sigBeta/numNeigh[l])*exp((-.5 * sigEpsInv)*(2*betaStar* scalarProduct(resid,Xl, nPat) - pow(betaStar,2) * scalarProduct(Xl, Xl, nPat))
	                                         + .5 * sigBetaInv*numNeigh[l]*pow(betaStar-betaBar,2)- a + b*tmp); 

	//unscaled:
	ppos = 1/(1+g);

	/* Bernoulli choice based on posterior probability */

	if(rbinom(1, ppos)==1)
	{
		//update sumPredNeigh for all neighbours of l
		for(i = 0; i < nDelta; i++)
		{
			if( matrix(delta, nVox, l, i) !=  NA_INTEGER) // if delta[l,i] != NA: delta[l,i] -1 is neighbour of l
				sumPredNeigh[matrix(delta, nVox, l, i)-1] = sumPredNeigh[matrix(delta, nVox, l, i)-1] - gamma[l]*beta[l] + betaStar; // correct for indices starting with 0 (C) vs. indices starting with 1 (R)
		}

		// update gamma and beta
		gamma[l] = 1;
		beta[l] = betaStar;
	}
	else 
	{
		//update sumPredNeigh for all neighbours of l
		for(i = 0; i < nDelta; i++)
		{
			if( matrix(delta, nVox, l, i) !=  NA_INTEGER) // if delta[l,i] != NA
				sumPredNeigh[matrix(delta, nVox, l, i)-1] = sumPredNeigh[matrix(delta, nVox, l, i)-1] - gamma[l]*beta[l]; // correct for indices starting with 0 (C) vs. indices starting with 1 (R)
		}

		// update gamma and beta
		gamma[l] = 0;
		beta[l] = 0;
	}
	
	return;
}


// update beta via single-site Gibbs sampler; based on Goldsmith et al., 2014
// ignore Ising field gamma (=set each voxel to 1)
void updateBetaOnly(int nPat, int nVox, int nDelta, //array sizes
                 int l, // actual voxel index
                 double* Xl, double* Y, // data
                 double* beta, // parameters
                 double* InProd, int* delta, double* xbeta0, double* FixEfAlpha, int* numNeigh, double* sumPredNeigh, //auxiliary variables 
                 double sigEps, double sigBeta) //hyperparameters
{
	/* auxiliary variables */
	unsigned int i;
	int indNeigh;
	double sigEpsInv, sigBetaInv, invNeighSize, betaBar, betaStar, Sigl, Mul, tmp, allNeigh, g, ppos;


	double resid[nPat]; // = Y- FixEfAlpha- xbeta0

	/* set auxiliary variables */
	sigEpsInv = 1/sigEps;
	sigBetaInv = 1/sigBeta;

	/* check neighbourhood */

	// empty neighbourhood -> continue with next voxel		
	if(numNeigh[l] == 0) 
	{
		beta[l] = 0;

		return;
	} 

	/* generate nonzero coefficient */

	betaBar = sumPredNeigh[l]/ (double) numNeigh[l]; //mean of neighbouring voxels

	//check neighbourhood of neighbours of l
	invNeighSize = 0; // sum of inverse neighbourhood sizes
	tmp = 0;

	// generate proposal sample from normal distribution with mean=Mul, sigma = sqrt(Sigl) 	
	Sigl = 1/(sigEpsInv * InProd[l] + sigBetaInv * numNeigh[l] );

	for(i = 0; i < nPat; i++)
		resid[i] = Y[i] - FixEfAlpha[i] - xbeta0[i];

	Mul = Sigl * ( sigEpsInv * scalarProduct(resid, Xl, nPat) + sigBetaInv * numNeigh[l] *  betaBar);

	betaStar = rnorm(Mul, sqrt(Sigl));


	//update sumPredNeigh for all neighbours of l
		for(i = 0; i < nDelta; i++)
		{
			if( matrix(delta, nVox, l, i) !=  NA_INTEGER) // if delta[l,i] != NA: delta[l,i] -1 is neighbour of l
				sumPredNeigh[matrix(delta, nVox, l, i)-1] = sumPredNeigh[matrix(delta, nVox, l, i)-1] - beta[l] + betaStar; // correct for indices starting with 0 (C) vs. indices starting with 1 (R)
		}

	
	beta[l] = betaStar;
	
	return;
}