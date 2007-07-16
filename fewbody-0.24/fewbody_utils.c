/* -*- linux-c -*- */
/* fewbody_utils.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "fewbody.h"

/* allocate a vector */
double *fb_malloc_vector(int n)
{
	return((double *) malloc(n * sizeof(double)));
}

/* allocate a matrix */
double **fb_malloc_matrix(int nr, int nc)
{
	int i;
	double **m;

	m = (double **) malloc(nr * sizeof(double *));

	m[0] = (double *) malloc(nr * nc * sizeof(double));

	for (i=1; i<nr; i++) {
		m[i] = m[i-1] + nc;
	}

	return(m);
}

/* free a vector */
void fb_free_vector(double *v)
{
	free(v);
}

/* free a matrix */
void fb_free_matrix(double **m)
{
	free(m[0]);
	free(m);
}

/* a fast square function */
double fb_sqr(double x)
{
	return(x*x);
}

/* a fast cube function */
double fb_cub(double x)
{
	return(x*x*x);
}

/* the dot product of two vectors */
double fb_dot(double x[3], double y[3])
{
	return(x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

/* the modulus of a vector */
double fb_mod(double x[3])
{
	return(sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
}

/* the cross product of two vectors */
int fb_cross(double x[3], double y[3], double z[3])
{
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];

	return(0);
}

/* something to calculate the angular momentum */
int fb_angmom(fb_obj_t *star, int nstar, double L[3])
{
	int i, k;
	double **l;

	l = fb_malloc_matrix(nstar, 3);

	for (k=0; k<3; k++) {
		L[k] = 0.0;
	}
	
	for (i=0; i<nstar; i++) {
		fb_cross(star[i].x, star[i].v, l[i]);
		for (k=0; k<3; k++) {
			L[k] += star[i].m * l[i][k];
		}
	}

	fb_free_matrix(l);
	
	return(0);
}

/* something to calculate the internal angular momentum */
void fb_angmomint(fb_obj_t *star, int nstar, double L[3])
{
	int i, k;

	for (k=0; k<3; k++) {
		L[k] = 0.0;
	}
	
	for (i=0; i<nstar; i++) {
		for (k=0; k<3; k++) {
			L[k] += star[i].Lint[k];
		}
	}
}

/* something to calculate the internal energy */
double fb_einttot(fb_obj_t *star, int nstar)
{
	int i;
	double eint=0.0;

	for (i=0; i<nstar; i++) {
		eint += star[i].Eint;
	}
     
	return(eint);
}

/* calculates the total potential energy of the system */
double fb_petot(fb_obj_t *star, int nstar)
{
	int i, j, k;
	double pe=0.0, r[3];

	for (i=0; i<nstar; i++) {
		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				r[k] = star[j].x[k] - star[i].x[k];
			}
			pe += -star[i].m * star[j].m / fb_mod(r);
		}
	}
	
	return(pe);
}

/* calculates the total kinetic energy of the system */
double fb_ketot(fb_obj_t *star, int nstar)
{
	int i;
	double ke=0.0;

	for (i=0; i<nstar; i++) {
		ke += 0.5 * star[i].m * fb_dot(star[i].v, star[i].v);
	}
     
	return(ke);
}

/* calculates the potential energy of the bound members of the system */
double fb_outerpetot(fb_obj_t **obj, int nobj)
{
	int i, j, k;
	double pe=0.0, r[3];

	for (i=0; i<nobj; i++) {
		for (j=i+1; j<nobj; j++) {
			for (k=0; k<3; k++) {
				r[k] = obj[j]->x[k] - obj[i]->x[k];
			}
			pe += -obj[i]->m * obj[j]->m / fb_mod(r);
		}
	}
	
	return(pe);
}

/* calculates the kinetic energy of the bound members of the system */
double fb_outerketot(fb_obj_t **obj, int nobj)
{
	int i;
	double ke=0.0;

	for (i=0; i<nobj; i++) {
		ke += 0.5 * obj[i]->m * fb_dot(obj[i]->v, obj[i]->v);
	}
     
	return(ke);
}

/* solve the Kepler equation for the eccentric anomaly, given the mean anomaly and eccentricity */
double fb_kepler(double e, double mean_anom)
{
	int status, iter;
	double ecc_anom, params[2];
	gsl_function F;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	
	/* set up the root solver */
	F.function = &fb_keplerfunc;
	F.params = &params;

	/* set the parameters */
	params[0] = e;
	params[1] = mean_anom;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, 0.0, 2.0*FB_CONST_PI);
	
	/* get eccentric anomaly by root-finding */
	iter = 0;
	do {
		iter++;
		gsl_root_fsolver_iterate(s);
		status = gsl_root_test_interval(gsl_root_fsolver_x_lower(s), gsl_root_fsolver_x_upper(s), \
						FB_ROOTSOLVER_ABS_ACC, FB_ROOTSOLVER_REL_ACC);
	} while (status == GSL_CONTINUE && iter < FB_ROOTSOLVER_MAX_ITER);

	if (iter >= FB_ROOTSOLVER_MAX_ITER) {
		fprintf(stderr, "Root finder failed to converge.\n");
		exit(1);
	}

	/* we've got the root */
	ecc_anom = gsl_root_fsolver_root(s);
	
	/* free memory associated with root solver */
	gsl_root_fsolver_free(s);

	return(ecc_anom);
}

/* the Kepler function for the root finder */
double fb_keplerfunc(double ecc_anom, void *params)
{
	double e, mean_anom;
	
	e = ((double *)params)[0];
	mean_anom = ((double *)params)[1];
	
	return(ecc_anom - e * sin(ecc_anom) - mean_anom);
}

/* calculate the relative tidal acceleration */
double fb_reltide(fb_obj_t *bin, fb_obj_t *single, double r)
{
	double arel, atid;

	arel = bin->m / fb_sqr(bin->a * (1.0 + bin->e));
	/* Factor in numerator is (single->m + bin->m) instead of single->m for the case of 
	   small single->m, for which we want to at some point resolve "bin" to get the 
	   motion of "single".  This will eventually be fixed when we we use full collapsed 
	   binary member positions in the integration, and not just the CM of the binary. */
	atid = 2.0 * (single->m + bin->m) / fb_cub(r) * bin->a * (1.0 + bin->e);

	return(atid/arel);
}
