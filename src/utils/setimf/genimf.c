#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "taus113-v2.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_double.h>
#include "genimf.h"

//----------------------------------------------------------------
// Generates an array of masses of length N based on the string
// called parse.  mmin, mmax describe the possible range of masses.
// Algorithm based on Pflamm-Altenburg 2008
// 
// parse should be of the form
//       mass1(alpha1)mass2(alpha2)mass3(alpha3)...(alphaN)
//
// example (Kroupa IMF)
//       0.01(0.3)0.08(1.3))0.5(2.3)1.0(2.3)
//
// Note: upper mass limit set by mmax, not parse.
//----------------------------------------------------------------
double* generateIMF(char *parse, long N, double mmin, double mmax)
{
  //-----------------------------------------------
  // Parse the input imf for mass, alphas
  //-----------------------------------------------
  int n=0;
  int i=0, j=0;
  double integ=0.0, mtemp=0.0;

  while (parse[i] != '\0'){
    if (parse[i] == '(' || parse[i] == ')')
      n++;
    i++;
  }
  
  char **sects = (char**)malloc(sizeof(char*)*n);
  for (i=0; i<n; i++)
    sects[i] = (char*)malloc(sizeof(char)*1024);

  i=0;
  char delims[] = "()";
  char *result = NULL;
  result = strtok( parse, delims );
  while( result != NULL ) {
    strncpy(sects[i] , result, 1024);
    result = strtok( NULL, delims );
    i++;
  } 

  n = n/2;
  double *alpha = (double*)malloc(sizeof(double)*n);
  double *mass  = (double*)malloc(sizeof(double)*(n+1));
  double *psi   = (double*)malloc(sizeof(double)*n);
  double *lims  = (double*)malloc(sizeof(double)*(n+1));
  double *m     = (double*)malloc(sizeof(double)*N);
  
  for (i=0; i<2*n; i+=2)
    {
      mass[i/2]  = strtod(sects[i], NULL);
      alpha[i/2] = strtod(sects[i+1], NULL);
    }
  
  mass[n] = POS_INFINITY;  
  lims[0] = 0.0;
  lims[n] = 1.0;
  psi[0] = 1.0;
  
  //-----------------------------------------------
  // Determine psi[i] as well as lims[i]
  //-----------------------------------------------
  for (i=1; i<n; i++)
    psi[i] = psi[i-1]*pow(mass[i], alpha[i]-alpha[i-1]);
  
  for (i=0; i<n; i++)
    integ += (psi[i]/(1-alpha[i]))*(pow(mass[i+1], 1-alpha[i]) - pow(mass[i], 1-alpha[i]));

  for (i=0; i<n; i++)
    psi[i] = psi[i]/integ;

  for (i=1; i<n; i++)
    {
      double val = 0.0;
      for (j=0; j<i; j++)
	val += psi[j]/(1-alpha[j])*(pow(mass[j+1], 1-alpha[j])-pow(mass[j], 1-alpha[j]));
      lims[i] = val;
    }

  //-----------------------------------------------
  // Generate IMF
  //-----------------------------------------------
  for(i=0; i<N; i++){
    do {
      double X = rng_t113_dbl();
      int j = 0;
      while (X > lims[j+1]) {
	j++;
      }
      mtemp = pow((1.0-alpha[j])/psi[j]*(X-lims[j])+pow(mass[j],1.0-alpha[j]), 1.0/(1.0-alpha[j]));
      if (isnan(mtemp)) {
	fprintf(stderr, "Oops!  m=NaN.  Please make coefficients more precise.\n");
	exit(-127);
      }
    } while (mtemp < mmin || mtemp > mmax) ;

    m[i] = mtemp;
  }

  return m;
}
