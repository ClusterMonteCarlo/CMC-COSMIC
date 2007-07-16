/* -*- linux-c -*- */
/* fewbody_nonks.c

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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "fewbody.h"

/* the derivatives function for the GSL ODE integrator */
#define FB_FM(i, j, k) fm[nstar*3*i + 3*j + k]
int fb_nonks_func(double t, const double *y, double *f, void *params)
{
	int i, j, k, nstar;
	double r[3], *m, *fm, val;

	nstar = (*(fb_nonks_params_t *) params).nstar;
	m = (*(fb_nonks_params_t *) params).m;

	fm = fb_malloc_vector(nstar * nstar * 3);

	/* calculate the matrix */
	for (i=0; i<nstar; i++) {
		for (j=0; j<i; j++) {
			for (k=0; k<3; k++) {
				FB_FM(i, j, k) = -FB_FM(j, i, k);
			}
		}

		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				r[k] = y[j*6+k] - y[i*6+k];	
			}
			
			val = 1.0 / fb_cub(fb_mod(r));

			for (k=0; k<3; k++) {
				FB_FM(i, j, k) = val * r[k];
			}
		}
	}

	/* calculate derivatives */
	for (i=0; i<nstar; i++) {
		for (k=0; k<3; k++) {
			f[i*6+k] = y[i*6+k+3];
			f[i*6+k+3] = 0.0;
		}
		
		for (j=0; j<i; j++) {
			for (k=0; k<3; k++) {
				f[i*6+k+3] += m[j] * FB_FM(i, j, k);
			}
		}

		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				f[i*6+k+3] += m[j] * FB_FM(i, j, k);
			}
		}
	}
	
	fb_free_vector(fm);

	return(GSL_SUCCESS);
}
#undef FB_FM

/* the Jacobian for the GSL ODE integrator */
int fb_nonks_jac(double t, const double *y, double *dfdy, double *dfdt, void *params)
{
	unsigned int i, j, k, a, b, kk;
	int nstar;
	double r[3], *m, val;
	gsl_matrix_view dfdy_mat;
	gsl_matrix *matrix;

	nstar = (*(fb_nonks_params_t *) params).nstar;
	m = (*(fb_nonks_params_t *) params).m;

	/* allocate matrices */
	dfdy_mat = gsl_matrix_view_array(dfdy, nstar*6, nstar*6);
	matrix = &dfdy_mat.matrix;
	
	/* set dfdt to zero */
	for (j=0; j<nstar*6; j++) {
		dfdt[j] = 0.0;
	}

	/* set the matrix to zero to begin with */
	gsl_matrix_set_zero(matrix);

	/* then set the actual values */
	for (i=0; i<nstar; i++) {
		for (a=0; a<3; a++) {
			gsl_matrix_set(matrix, i+a, i+a+3, 1.0);

			for (k=0; k<nstar; k++) {
				for (b=0; b<3; b++) {
					val = 0.0;
					for (j=0; j<nstar; j++) {
						if (j != i) {
							for (kk=0; kk<3; kk++) {
								r[kk] = y[i*6+kk] - y[j*6+kk];
							}
							val -= m[j] * ((double) (FB_DELTA(i, k) - FB_DELTA(j, k))) / fb_cub(fb_mod(r)) * \
								(((double) FB_DELTA(a, b)) - 3.0 * (y[i*6+a] - y[j*6+a])*(y[i*6+b] - y[j*6+b])/fb_dot(r, r));
						}
					}
					gsl_matrix_set(matrix, i+a+3, k+b, val);
				}
			}
		}
	}
	
	return(GSL_SUCCESS);
}

void fb_euclidean_to_nonks(fb_obj_t **star, double *y, int nstar)
{
	int i, j;
	
	for (i=0; i<nstar; i++) {
		for (j=0; j<3; j++) {
			y[i*6+j] = star[i]->x[j];
			y[i*6+j+3] = star[i]->v[j];
		}
	}
}

void fb_nonks_to_euclidean(double *y, fb_obj_t **star, int nstar)
{
	int i, j;

	for (i=0; i<nstar; i++) {
		for (j=0; j<3; j++) {
			star[i]->x[j] = y[i*6+j];
			star[i]->v[j] = y[i*6+j+3];
		}
	}	
}
