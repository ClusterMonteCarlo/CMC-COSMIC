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
#define FB_REL(i, j, k) fmr[nstar*3*i + 3*j + k]
#define pi2 9.869604401089359

int fb_nonks_func(double t, const double *y, double *f, void *params)
{
	int i, j, k, nstar, PN1, PN2, PN25, PN3, PN35;
	double n[3],r[3],v[3], *m, *fm, *fmr, val, rdot, SM, nu;
	double  A=0.,B=0.,A2=0.,B2=0.,A4=0.,B4=0.,A5=0.,B5=0.,A6=0.,B6=0.,A7=0.,B7=0.;
	double clight, clight2, clight4, clight5, R, Mclus, aclus;
	double SM2, SM3;
	double rdot2,rdot3,rdot4,rdot6,v2,v4,v6;
	double r_mod, r_mod2, r_mod3;
	fb_units_t units;
	
	nstar = (*(fb_nonks_params_t *) params).nstar;
	m = (*(fb_nonks_params_t *) params).m;
	PN1 = (*(fb_nonks_params_t *) params).PN1;
	PN2 = (*(fb_nonks_params_t *) params).PN2;
	PN25 = (*(fb_nonks_params_t *) params).PN25;
	PN3 = (*(fb_nonks_params_t *) params).PN3;
	PN35 = (*(fb_nonks_params_t *) params).PN35;
	units = (*(fb_nonks_params_t *) params).units;

	clight = FB_CONST_C / units.v;
	clight2 = fb_sqr(clight);
	clight4 = fb_sqr(clight2);
	clight5 = clight4 * clight;

	fm = fb_malloc_vector(nstar * nstar * 3);
	fmr = fb_malloc_vector(nstar * nstar * 3);

	/* calculate the matrix */
	for (i=0; i<nstar; i++) {
		for (j=0; j<i; j++) {
			for (k=0; k<3; k++) {
				FB_FM(i, j, k) = -FB_FM(j, i, k);
				FB_REL(i, j, k) = -FB_REL(j, i, k);
			}
		}
		for (j=i+1; j<nstar; j++) {
			
			/* First, precalculate a lot of quantities */
			SM = m[i] + m[j];
			SM2 = SM*SM;
			nu = m[i]*m[j]/(SM2);
			SM3 = SM2*SM;
	
		
			for (k=0; k<3; k++) {
				r[k] = y[j*6+k] - y[i*6+k];	
				v[k] = y[j*6+k+3] - y[i*6+k+3]; 
			}

			rdot = 0;
			r_mod = fb_mod(r);
			r_mod2 = r_mod*r_mod;
			r_mod3 = r_mod2*r_mod;

			for (k=0; k<3; k++) {
				n[k] = r[k]/r_mod; 
				rdot += n[k]*v[k]; 
			}

			v2 = fb_mod(v)*fb_mod(v);
			v4 = v2*v2;
			v6 = v4*v2;
			rdot2 = rdot*rdot;
			rdot3 = rdot2*rdot;
			rdot4 = rdot2*rdot2;
			rdot6 = rdot4*rdot2;
			
			if (PN1) {
				A2 = (-3*rdot2*nu/2 + v2 + 3*nu*v2 - \
				      SM*(4+2*nu)/r_mod)/clight2;
				B2 = (-4 + 2*nu)*rdot/clight2;
			} 

			if (PN2) {
				A4 = (15*rdot4*nu/8-45*rdot4*nu*nu/8-\
				      9*rdot2*nu*v2/2+6*rdot2*nu*nu*v2+\
				      3*nu*v4-4*nu*nu*v2*\
				      v2+SM*(-2*rdot2-25*rdot2*nu-2*rdot2*nu*nu-\
							    13*nu*v2/2+2*nu*nu*v2)\
				      /r_mod+SM2*(9+87*nu/4)/r_mod2)/clight4;
				B4 = (9*rdot3*nu/2+3*rdot3*nu*nu-15*rdot*nu*v2/\
				      2-2*rdot*nu*nu*v2+SM*(2*rdot+41*rdot*nu/2+4*rdot*nu*nu)/\
				      r_mod)/clight4;
			}

			if (PN25) {
				A5 = (-24*rdot*nu*v2*SM/(5*r_mod)-\
				      136*rdot*nu*SM2/(15*r_mod2))/clight5;
				B5 = (8*nu*v2*SM/(5*r_mod) + \
				      24*nu*SM2/(5*r_mod2))/clight5;
			}
			
			
			if (PN3) {
				A6 = -((16+(1399/12-41*pi2/16)*nu+35.5*nu*nu)*SM3/r_mod3 +\
				       nu*v2*SM2*(20827/840+123*pi2/64-nu*nu)/r_mod2 - \
				       rdot2*SM2*(1+(22717/168+615/64*pi2)*nu+11*nu*nu/8-7*nu*nu*nu)/r_mod2 - \
				       0.25*nu*v6*(11-49*nu+52*nu*nu) + \
				       35*rdot6*nu*(1-5*nu+5*nu*nu)/16 -\
				       0.25*nu*SM*v4*(75+32*nu-40*nu*nu)/r_mod - \
				       0.5*nu*rdot4*SM*(158-69*nu-60*nu*nu)/r_mod + \
				       nu*SM*rdot2*v2*(121-16*nu-20*nu*nu)/r_mod  + \
				       3*nu*v4*rdot2*(20-79*nu+60*nu*nu)/8 -\
				       15*nu*rdot4*v2*(4-18*nu+17*nu*nu)/8 )/clight4/clight2;
				B6 = -rdot*((4+((5849/840)+(123/32)*pi2)*nu-25*nu*nu-8*nu*nu*nu)*SM2/\
					    r_mod2+nu*v4*(65-152*nu-48*nu*nu)/8+\
					    15*nu*rdot4*(3-8*nu-2*nu*nu)/8+nu*(15+27*nu+10*nu*nu)*v2*\
					    SM/r_mod-nu*SM*rdot2*(329+177*nu+108*nu*nu)/r_mod/6-\
					    0.75*nu*rdot2*v2*(16-37*nu-16*nu*nu))/clight4/clight2;
			}

			if (PN35) {
				A7 =  1.6*nu*SM*rdot*(23*SM2*(43+14*nu)/r_mod2/14+\
						      3*v4*(61+70*nu)/\
						      28+70*rdot4+SM*v2*(519-1267*nu)/42/r_mod +\
						      SM*rdot2*(147+188*nu)/4/r_mod-\
						      15*rdot2*v2*(19+2*nu)/4)/r_mod/clight5/clight2;
				B7 = -1.6*nu*SM*(SM2*(1325+546*nu)/r_mod2/42+\
						 v4*(313+42*nu)/28+75*rdot4-\
						 SM*v2*(205+777*nu)/r_mod/42   +\
						 SM*rdot2*(205+424*nu)/r_mod/12-\
						 0.75*rdot2*v2*(113+2*nu))/r_mod/clight5/clight2;
			}

			A = A2 + A4 + A5 + A6 + A7;
			B = B2 + B4 + B5 + B6 + B7;


			val = 1.0 / r_mod3;
			for (k=0; k<3; k++) {
				FB_FM(i, j, k) = val * r[k];
				FB_REL(i, j, k) = (A*n[k]+B*v[k]) / r_mod2;
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
				f[i*6+k+3] += m[j] * FB_FM(i, j, k) + m[j] * FB_REL(i,j,k);
			}
		}

		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				f[i*6+k+3] += m[j] * FB_FM(i, j, k) + m[j] * FB_REL(i,j,k);
			}
		}
	}

	fb_free_vector(fm);
	fb_free_vector(fmr);

	return(GSL_SUCCESS);
}
#undef FB_FM
#undef FB_REL

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
