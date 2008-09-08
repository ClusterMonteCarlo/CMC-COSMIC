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
#ifndef EXPERIMENTAL
#define FUNC(j, k, E, J) (2.0 * SQR(star[(k)].r) * ((E) - (star[(k)].phi + PHI_S(star[k].r, j))) - SQR(J))
#else
#define FUNC(j, k, E, J) (2.0 * ((E) - (star[(k)].phi + PHI_S(star[k].r, j))) - SQR((J)/star[(k)].r))
#endif
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
double calc_pot_within_interval(double r, void *p) {
  double pot;
  struct calc_vr_params *params= (struct calc_vr_params *) p;
  long index, k;
  double E, J;

  index= params->index; k= params->k;
  E= params->E; J= params->J;
  
  if (r< (star[k].r-DBL_EPSILON) || r> star[k+1].r+DBL_EPSILON) {
    eprintf("r= %g is not in [%g,%g]! r-r_low= %g, r-r_high= %g\n", r, 
        star[k].r, star[k+1].r, r-star[k].r, r-star[k+1].r);
    exit_cleanly(-1);
  };
  if (fabs(r-star[k].r)< DBL_EPSILON) {
    pot= star[k].phi;
  } else if (fabs(r-star[k+1].r)< DBL_EPSILON) {
    pot= star[k+1].phi;
  } else {
    if (r< star[1].r) {
      pot= star[0].phi-cenma.m*madhoc/r;
    } else {
      pot= (star[k].phi + (star[k + 1].phi - star[k].phi) 
                           * (1.0/star[k].r - 1.0/r) /
                           (1.0/star[k].r - 1.0/star[k + 1].r));
    }
  };

  return(pot);
};

double calc_pot_in_interval(double r, long k) {
  double pot;
  if (r< (star[k].r-DBL_EPSILON) || r> star[k+1].r+DBL_EPSILON) {
    eprintf("r= %g is not in [%g,%g]! r-r_low= %g, r-r_high= %g\n", r, 
        star[k].r, star[k+1].r, r-star[k].r, r-star[k+1].r);
    exit_cleanly(-1);
  };
  if (fabs(r-star[k].r)< DBL_EPSILON) {
    pot= star[k].phi;
  } else if (fabs(r-star[k+1].r)< DBL_EPSILON) {
    pot= star[k+1].phi;
  } else {
    if (r< star[1].r) {
      pot= star[0].phi-cenma.m*madhoc/r;
    } else {
      pot= (star[k].phi + (star[k + 1].phi - star[k].phi) 
                           * (1.0/star[k].r - 1.0/r) /
                           (1.0/star[k].r - 1.0/star[k + 1].r));
    } 
  };

  return(pot);
};


long find_zero_Q(long j, long kmin, long kmax, long double E, long double J){
  /* another binary search:
   * anologous to above, except FUNC(k) may be decreasing 
   * rather than increasing */
  long ktry, kmax1, fevals;
  long double rmax, pot_max, rmin, pot_min;
  long double rtry, pot_try;
  long double qmin, qmax, qtry, func;

  fevals= 0;
  kmax1= kmax;
  rmax= star[kmax].r; pot_max= star[kmax].phi;
  rmin= star[kmin].r; pot_min= star[kmin].phi;
  qmin= function_q(j, rmin, pot_min, E, J);
  qmax= function_q(j, rmax, pot_max, E, J);
  fevals+= 2;
  /*
   *if (j==3265) {
   *  printf("search starts at kmin=%li, rmin=%Lg, q[kmin]=%Lg\n", kmin, rmin, qmin);
   *  printf("search starts at kmax=%li, rmax=%Lg, q[kmax]=%Lg pot[kmax]=%Lg\n", kmax, rmax, qmax, pot_max);
   *} else {
   *  dprintf("search starts at kmin=%li, rmin=%Lg, q[kmin]=%Lg\n", kmin, rmin, qmin);
   *  dprintf("search starts at kmax=%li, rmax=%Lg, q[kmax]=%Lg\n", kmax, rmax, qmax);
   *  dprintf("E=%Lg, J=%Lg, pot_min=%Lg, pot_max=%Lg\n", E, J, pot_min, pot_max);
   *}
   */
  if(qmin< qmax){
    //dprintf("increasing\n");
    do {
      ktry = (kmin+kmax+1)/2;
      rtry= star[ktry].r; pot_try= star[ktry].phi;
      qtry= function_q(j, rtry, pot_try, E, J);
      fevals++;
      //dprintf("ktry=%li, q[ktry]=%g\n", ktry, (double) qtry);
      if (qtry<0.e0){
        kmin = ktry;
      } else {
        kmax = ktry-1;
      }
    } while (kmax!=kmin);
  } else {
    //dprintf("decreasing\n");
    do {
      ktry = (kmin+kmax+1)/2;
      rtry= star[ktry].r; pot_try= star[ktry].phi;
      qtry= function_q(j, rtry, pot_try, E, J);
      fevals++;
      if (qtry>0.e0){
        kmin = ktry;
      } else {
        kmax = ktry-1;
      }
    } while (kmax!=kmin);
  };

  /*
   *if (j==3265) {
   *  qtry= function_q(j, star[kmin].r, star[kmin].phi, E, J);
   *  printf("search ends at kmin=%li, q[kmin]=%Lg\n", kmin, qtry);
   *  qtry= function_q(j, star[kmin+1].r, star[kmin+1].phi, E, J);
   *  printf("search ends at rmin=%g, q[kmin+1]=%Lg\n", star[kmin].r, qtry);
   *  func= FUNC(j, kmax1, E, J);
   *  printf("E=%Lg, J=%Lg, FUNC[max]=%Lg, N_MAX+1=%li\n", E, J, func, clus.N_MAX +1);
   *}
   */
  //printf("f_z_Q fevals=%li\n", fevals);
  return (kmin);
} 

/* Calculate the square of vr ! */
double calc_vr_within_interval(double r, void *p) {
  struct calc_vr_params *params= (struct calc_vr_params *) p;
  long index, k;
  double pot, E, J, vr;

  index= params->index; k= params->k;
  E= params->E; J= params->J;
  pot= calc_pot_within_interval(r, p);
  vr= function_q(index, r, pot, E, J);
  return(vr);
};

/* Calculates the square of vr */
double calc_vr_in_interval(double r, long index, long k, double E, double J) {
  struct calc_vr_params params;
  double vr;

  params. index = index;
  params. k     = k;
  params. E     = E;
  params. J     = J;

  vr= calc_vr_within_interval(r, &params);
  
  return(vr);
};

double calc_Q_within_interval(double r, void *p) {
  struct calc_vr_params *params= (struct calc_vr_params *) p;
  double pot, E, J;
  long index;

  index= params->index;
  E= params->E; J= params->J;
  pot= calc_pot_within_interval(r, p);
  
  return(2.0 * SQR(r) * (E - pot - PHI_S(r, index))- SQR(J));
};

/* Calculates the square of vr! */
double calc_vr(double r, long index, double E, double J) {
  long k;
  struct Interval inter;
  struct calc_vr_params params;
  double vr;

  if (SEARCH_GRID) {
    inter= search_grid_get_interval(r_grid, r);
  } else {
    inter.min= 0;
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
//#define FUNC(j, k, E, J) (2.0 * SQR(star[(k)].r) * ((E) - (star[(k)].phi + PHI_S(star[k].r, j))) - SQR(J))
double find_root_vr(long index, long k, double E, double J) {
  struct calc_vr_params p;
  int status, not_converged;
  double r_low=-1.0, r_high=-1.0, apsis, prev_apsis= -1.0;
  long iter;
  gsl_function F;

  not_converged= 1;
  iter= APSIDES_MAX_ITER;
  p.index= index; p.k= k;
  p.E= E; p.J= J;
  F.function= &calc_vr_within_interval;
  F.params= &p;

  /* check if the values of F at the interval boundaries have different signs */
  if (GSL_SIGN(GSL_FN_EVAL(&F, star[k].r))==GSL_SIGN(GSL_FN_EVAL(&F, star[k+1].r))) {
    eprintf("The signs of F[k] and F[k+1] are the same!\n");
    eprintf("k= %li, F[k]= %g, F[k+1]= %g, r[k]= %g, r[k+1]= %g\n", k, 
        GSL_FN_EVAL(&F, star[k].r), GSL_FN_EVAL(&F, star[k+1].r), star[k].r, star[k+1].r);
    eprintf("FUNC[k]= %g, FUNC[k+1]= %g\n", FUNC(index, k, E, J), FUNC(index, k+1, E, J));
    eprintf("DBL_EPSILON= %g, r[k]= %g\n", DBL_EPSILON, star[k].r);
    exit(1);
  }

  status= gsl_root_fsolver_set(q_root, &F, star[k].r, star[k+1].r);
  if (status) {
    eprintf("Initialization of root solver failed! Error Code: %i\n", status);
  exit(1);
  } else {
    //dprintf("Initialization successful!\n");
  };

  while(not_converged && iter) {
    status= gsl_root_fsolver_iterate(q_root);
    if (!status) {
      r_low= gsl_root_fsolver_x_lower(q_root);
      r_high= gsl_root_fsolver_x_upper(q_root);
      not_converged= (gsl_root_test_interval(r_low, r_high, APSIDES_PRECISION, APSIDES_PRECISION)==GSL_CONTINUE);
    } else {
      if (status== GSL_EBADFUNC) {
        eprintf("Iteration encountered a singular point, i.e. vr= Inf or NaN!\n");
        exit(1);
      } else {
        eprintf("Iteration stopped due to an unknown error! Error Code: %i\n", status);
        exit(1);
      };
    };
    iter--;
    if (iter<10) {
      dprintf("Iterated %li times and did not get apside error down to %g!\n", APSIDES_MAX_ITER-iter, APSIDES_PRECISION);
      dprintf("The current interval is [%.16g,%.16g]\n", r_low, r_high);
      dprintf("Interval width is %.14g\n", r_high-r_low);
      dprintf("Distance between the stars is %g\n", star[k].r- star[k+1].r);
      dprintf("Values of vr range from %g to %g\n", GSL_FN_EVAL(&F, r_low), GSL_FN_EVAL(&F, r_high));
      if (prev_apsis< 0.) {
	dprintf("Consider now APSIDES_CONVERGENCE= %g.\n", APSIDES_CONVERGENCE);
	prev_apsis= gsl_root_fsolver_root(q_root);
      } else {
	apsis= gsl_root_fsolver_root(q_root);
	not_converged=  not_converged && 
          (gsl_root_test_delta(apsis, prev_apsis, APSIDES_CONVERGENCE, APSIDES_CONVERGENCE)==GSL_CONTINUE);
	prev_apsis= apsis;
      };
    };
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
    dprintf("Distance between the stars is %g\n", star[k].r- star[k+1].r);
    dprintf("Values of vr range from %g to %g\n", GSL_FN_EVAL(&F, r_low), GSL_FN_EVAL(&F, r_high));
    exit(1);
  };


  if (r_high-r_low>APSIDES_PRECISION+APSIDES_PRECISION*MIN(r_high, r_low)) 
    dprintf("Wrong assumption!!! delta_r=%g, prec=%g\n", r_high-r_low, 
APSIDES_PRECISION+APSIDES_PRECISION*MIN(r_high,r_low));

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
    exit_cleanly(-1);
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
