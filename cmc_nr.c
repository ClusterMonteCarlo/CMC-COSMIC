/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

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
//	rtb = f < 0.0 ? (rdx = -1, fdx = 1, dx = x2 - x1, x1) : (rdx = 0, fdx = -1, dx = x1 - x2, x2);
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

/* Find Zero OF Q */
/* Finding the root of the function defined in 'func' by bisection.  */
/* Requires the root to be bracketed between x1 and x2. */
#define JMAX 500
#define FUNC(k, E, J) \
	(2.0 * (E - star[k].phi) - SQR(J / star[k].r))

long FindZero_Q(long x1, long x2, double E, double J)
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
