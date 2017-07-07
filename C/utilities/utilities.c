/* General utility functions */
#include <stdio.h>
#include <time.h> // for timestamp
#include <R.h> 
#include <Rmath.h>
#include <math.h> // for pi = M_PI

//macros
#define matrix(array, nRow,i,j) array[j*(nRow)+i] // for simpler matrix access


// print timestamp (from http://cc.byexamples.com)
void timestamp()
{
	time_t ltime; /* calendar time */
	struct tm *Tm;

	ltime=time(NULL); /* get current cal time */
	Tm=localtime(&ltime);

	printf("%02d.%02d.%04d, %02d:%02d:%02d ",
	       Tm->tm_mday,
	       Tm->tm_mon+1,
	       Tm->tm_year+1900,
	       Tm->tm_hour,
	       Tm->tm_min,
	       Tm->tm_sec);
}


// calculate scalar product of arrays a, b, both with length n
double scalarProduct(double* a, double* b, int n)
{
	double sum = 0;
	unsigned int i;

	for(i = 0; i < n; i++)
		sum = sum + a[i]*b[i];

	return(sum);	
}


// calculate Matrix-Vector multiplication A*b (must be both double!)
void multMatrixVector(double* A, double* b, double* res, int nrow, int ncol)
{
	unsigned int i, j;

	for(i = 0; i < nrow; i++)
	{
		res[i] = 0;

		for(j = 0; j < ncol; j++)
			res[i] = res[i] + matrix(A, nrow, i,j)*b[j];
	}

	return;
}


// sampling without replacement for random update; adapted from https://stat.ethz.ch/pipermail/r-help/2009-April/194193.html
void sample(int n, int *res)
{
	int i, j, N;
	int tmp[n];

	for (i = 0; i < n; i++)
		tmp[i] = i;

	N = n;

	for (i = 0; i < N; i++)
	{
		j = n * unif_rand();
		res[i] = tmp[j]; 
		tmp[j] = tmp[--n]; //i.e. n = n-1; tmp[j] = tmp[n];
	}

	return;
}


// check upper quantile
// params: val: a critical value to be tested
//         upper: an array containing the upper quantiles. Must be sorted
//        len: an integer, the length of upper
void checkUpper(double val, double *upper, int len)
{
  int i,j;
  
  if(val < upper[0]) // val smaller than smallest value
    return;
  
  i = len - 1;
  while(val <= upper[i] & i > 0)
    i--;
  
  for(j = 0; j < i; j++)
    upper[j] = upper[j+1];
  
  upper[i] = val;
  
  return;
}


// check lower quantile
// params: val: a critical value to be tested
//         lower: an array containing the lower quantiles. Must be sorted
//        len: an integer, the length of lower
void checkLower(double val, double *lower, int len)
{
  int i,j;
  
  if(val > lower[len -1]) // val bigger than biggest value
    return;
  
  i = 0;
  while(val >= lower[i] & i < len-1)
    i++;
  
  for(j = len-1; j > i; j--)
    lower[j] = lower[j-1];
  
  lower[i] = val;
  
  return;
}