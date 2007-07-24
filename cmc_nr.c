/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

long FindZero_r(long kmin, long kmax, double r){
	/* binary search:
	 * given the array star[].r and the two indices kmin and kmax,
	 * with the conditions 
	 * 1) array is monotonic in its indices,
	 * 2) kmin<kmax,
	 * find and return the index k, such that star[k].r<r<star[k+1].r */
	long ktry;

        if ((star[kmin].r>r && kmin>1) || star[kmax].r<r) {
          dprintf("r is outside kmin kmax!!\n");
          dprintf("star[kmin].r= %lf, star[kmax].r= %lf, kmin= %li, kmax= %li, r=%lf\n", 
                star[kmin].r, star[kmax].r, kmin, kmax, r);
        };

	do {
		ktry = (kmin+kmax+1)/2;
		if (star[ktry].r<r){
			kmin = ktry;
		} else {
			kmax = ktry-1;
		}
	} while (kmax!=kmin);

	return kmin;
}

double sigma_r(double r){
	/* binary search:
	 * given the array sigma_array.r[] and the two indices kmin and kmax,
	 * with the conditions 
	 * 1) array is monotonic in its indices,
	 * 2) kmin<kmax,
	 * find the index k, such that sigma_array.r[k]<r<sigma_array.r[k+1] */
	long ktry, kmin=1, kmax=sigma_array.n;
	do {
		ktry = (kmin+kmax+1)/2;
		if (sigma_array.r[ktry]<r){
			kmin = ktry;
		} else {
			kmax = ktry-1;
		}
	} while (kmax!=kmin);
	
	/* don't even bother with interpolating */
	return(sigma_array.sigma[kmin]);
}

#if 0

#define JMAX 500
#define FUNC(k, r) (star[(k)].r - (r))

long FindZero_r(long x1, long x2, double r)
{
	int j;
	long xmid, rtb, dx, fdx, rdx;
	double f, fmid;

	f = FUNC(x1, r);
	fmid = FUNC(x2, r);
	if (f * fmid >= 0.0) {
		eprintf("x1=%ld x2=%ld r=%g f=%g fmid=%g\n", x1, x2, r, f, fmid);
		nrerror("Root must be bracketed for bisection in FindZero_r");
	}
	if (f<0.0){
		rdx = -1;
		fdx = 1;
		dx = x2 - x1;
		rtb = x1;
	} else {
		rdx = 0;
		fdx = -1;
		dx = x1 - x2;
		rtb = x2;
	}
	for (j = 1; j <= JMAX; j++) {
		xmid = rtb + (dx /= 2);
		fmid = FUNC(xmid, r);
		if (fmid <= 0.0)
			rtb = xmid;
		if (dx == 0) {
			if (fmid > 0.0)
				nrerror("f > 0 in FindZero_r");
			/* look for change of sign */
			while (FUNC(rtb, r) < 0.0)
				rtb += fdx;
			return(rtb + rdx);
		} else if (fmid == 0.0) {
			return(rtb + rdx);
		}
	}
	
	/* the routine only gets here if the root bisection fails */
	nrerror("Too many bisections in FindZero_r");
	return(0); /* never get here */
}

#undef FUNC
#undef JMAX
#endif

/*
#define FUNC(k, E, J) \
       (2.0 * (E - star[k].phi) - SQR(J / star[k].r))
*/
/* #define FUNC(j, k, E, J) (2.0 * (E - (star[k].phi)) - SQR(J / star[k].r)) */
#define FUNC(j, k, E, J) (2.0 * SQR(star[(k)].r) * ((E) - (star[(k)].phi + PHI_S(star[k].r, j))) - SQR(J))
//#define FUNC(j, k, E, J) (2.0 * ((E) - (star[(k)].phi + PHI_S(star[k].r, j))) - SQR(J/star[(k)].r))
long FindZero_Q(long j, long kmin, long kmax, double E, double J){
	/* another binary search:
	 * anologous to above, except FUNC(k) may be decreasing 
	 * rather than increasing */
	long ktry;

	if(FUNC(j, kmin, E, J)<FUNC(j, kmax, E, J)){
		do {
			ktry = (kmin+kmax+1)/2;
			if (FUNC(j, ktry, E,J)<0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	} else {
		do {
			ktry = (kmin+kmax+1)/2;
			if (FUNC(j, ktry, E,J)>0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	}


	return kmin;
}

/* calculates the radial velocity in the interval [star[k].r, star[k+1].r] */
struct calc_vr_params {
  long index, k;
  double E, J;
};

double calc_vr_within_interval(double r, void *p) {
  double pot;
  struct calc_vr_params *params= (struct calc_vr_params *) p;
  long index, k;
  double E, J;

  index= params->index; k= params->k;
  E= params->E; J= params->J;
  pot= (star[k].phi + (star[k + 1].phi - star[k].phi) 
			 * (1.0/star[k].r - 1.0/r) /
			 (1.0/star[k].r - 1.0/star[k + 1].r));

  return(2.0 * E - fb_sqr(J/r) - 2.0 * (pot + PHI_S(r, index)));
};

double calc_vr(double r, long index, double E, double J) {
  long k;
  struct Interval inter;
  struct calc_vr_params params;
  double vr;

  if (SEARCH_GRID) {
    inter= search_grid_get_interval(r_grid, r);
  } else {
    inter.min= 1;
    inter.max= clus.N_MAX+1;
  };
  if (inter.min== inter.max-1) {
    k= inter.min;
  } else {
    k= FindZero_r(inter.min, inter.max, r);
  };

  params.k= k; params.index= index;
  params.E= E; params.J= J;
  /*dprintf("Interval in calc_vr for index %li is [%li,%li]\n", index, k, k+1);*/
  vr= calc_vr_within_interval(r, (void *) &params);

  return(vr);
};

/* Calculates the root of vr with the additional constraint that when 
 * substituting it back into calc_vr_within_interval leads to vr>= 0.0 .
 */
double find_root_vr(long index, long k, double E, double J) {
  struct calc_vr_params p;
  int status, not_converged;
  double r_low=-1.0, r_high=-1.0, apsis;
  long iter;
  gsl_function F;

  not_converged= 1;
  iter= APSIDES_MAX_ITER;
  p.index= index; p.k= k;
  p.E= E; p.J= J;
  F.function= &calc_vr_within_interval;
  F.params= &p;

  status= gsl_root_fsolver_set(q_root, &F, star[k].r, star[k+1].r);
  if (status) {
    eprintf("Initialization of root solver failed! Error Code: %i\n", status);
  exit(1);
  };

  while(not_converged && iter) {
    status= gsl_root_fsolver_iterate(q_root);
    if (!status) {
      r_low= gsl_root_fsolver_x_lower(q_root);
      r_high= gsl_root_fsolver_x_upper(q_root);
      not_converged= (gsl_root_test_interval(r_low, r_high, APSIDES_PRECISION, 0.)==GSL_CONTINUE);
    } else {
      if (status== GSL_EBADFUNC) {
        eprintf("Iteration encountered a singular point, i.e. vr= Inf or NaN!\n");
        exit(1);
      } else {
        eprintf("Iteration stopped due to an unknown error! Error Code: %i\n", status);
        exit(1);
      };
    };
    iter++;
  };

  if (!not_converged) {
    /*
     *dprintf("Found zero in interval [%.12g,%.12g]\n", r_low, r_high);
     *dprintf("Interval width is %g\n", r_high-r_low);
     *dprintf("Star interval width is %g\n", star[k+1].r-star[k].r);
     *dprintf("Values of vr range from %g to %g\n", GSL_FN_EVAL(&F, r_low), GSL_FN_EVAL(&F, r_high));
     *dprintf("Interval in find_root_vr for index %li is [%li, %li]\n", index, k, k+1);
     */
  } else {
    dprintf("Found NO zero in interval [%.12g,%.12g]\n", r_low, r_high);
    dprintf("Interval width is %g\n", r_high-r_low);
    dprintf("Values of vr range from %g to %g\n", GSL_FN_EVAL(&F, r_low), GSL_FN_EVAL(&F, r_high));
    exit(1);
  };


  if (r_high-r_low>APSIDES_PRECISION) dprintf("Wrong assumption!!!\n");

  apsis= gsl_root_fsolver_root(q_root);
  if (GSL_FN_EVAL(&F,apsis)< 0.) {
    if (GSL_FN_EVAL(&F, r_low)<0.) {
      apsis= r_high;
    };
    if (GSL_FN_EVAL(&F, r_high)< 0.) {
      apsis= r_low;
    };
  };

  if (GSL_FN_EVAL(&F, apsis)< 0.) {
    eprintf("Residual is negative! Apsis= %g, vr(apsis)= %g\n", apsis, GSL_FN_EVAL(&F, apsis));
    exit(1);
  }

  return(apsis);
};

#if 0

/* Find Zero OF Q */
/* Finding the root of the function defined in 'func' by bisection.  */
/* Requires the root to be bracketed between x1 and x2. */
#define JMAX 500
#define FUNC(k, E, J) \
	(2.0 * (E - star[k].phi) - SQR(J / star[k].r))

long FindZero_Q_old(long x1, long x2, double E, double J)
{
	void nrerror(char error_text[]);
	int j;
	long xmid, rtb, dx, fdx, rdx;
	double f, fmid;

	f = FUNC(x1, E, J);
	fmid = FUNC(x2, E, J);
	if (f * fmid >= 0.0) {
		eprintf("x1=%ld x2=%ld E=%g J=%g f=%g fmid=%g\n", x1, x2, E, J, f, fmid);
		nrerror("Root must be bracketed for bisection in FindZero_Q");
	}
	rtb = f < 0.0 ? (rdx = -1, fdx = 1,  dx = x2 - x1, x1) 
		      : (rdx =  0, fdx = -1, dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		xmid = rtb + (dx /= 2);
		fmid = FUNC(xmid, E, J);
		if (fmid <= 0.0)
			rtb = xmid;
		if (dx == 0) {
			if (fmid > 0.0)
				nrerror("f > 0 in FindZero_Q");
			/* look for change of sign */
			while ((FUNC(rtb, E, J)) < 0.0)
				rtb += fdx;
			return(rtb + rdx);
		} else if (fmid == 0.0) {
			return(rtb + rdx);
		}
	}

	/* the routine only gets here if the root bisection fails */
	nrerror("Too many bisections in FindZero_Q");
	return(0); /* never get here */
}

#undef FUNC
#undef JMAX

/* this implemantation breaks the code when -ffloat-store option is used 
 * so I renamed it.*/
/* find zero of Q */
/* Requires the root to be bracketed between x1 and x2. */
#define JMAX 500
#define FUNC(k, E, J) (2.0 * SQR(star[(k)].r) * ((E) - star[(k)].phi) - SQR(J))
long new_FindZero_Q(long x1, long x2, double E, double J)
{
	int j;
	long xmid, rtb, dx, fdx, rdx;
	double f, fmid;

	f = FUNC(x1, E, J);
	fmid = FUNC(x2, E, J);
	if (f * fmid >= 0.0) {
		eprintf("x1=%ld x2=%ld E=%g J=%g f=%g fmid=%g\n", x1, x2, E, J, f, fmid);
		nrerror("Root must be bracketed for bisection in FindZero_Q");
	}
	rtb = f < 0.0 ? (rdx = -1, fdx = 1, dx = x2 - x1, x1) : (rdx = 0, fdx = -1, dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		xmid = rtb + (dx /= 2);
		fmid = FUNC(xmid, E, J);
		if (fmid <= 0.0)
			rtb = xmid;
		if (dx == 0) {
			if (fmid > 0.0)
				nrerror("f > 0 in FindZero_Q");
			/* look for change of sign */
			while (FUNC(rtb, E, J) < 0.0)
				rtb += fdx;
			return(rtb + rdx);
		} else if (fmid == 0.0) {
			return(rtb + rdx);
		}
	}
	
	/* the routine only gets here if the root bisection fails */
	nrerror("Too many bisections in FindZero_Q");
	return(0); /* never get here */
}

#undef FUNC
#undef JMAX

#endif
