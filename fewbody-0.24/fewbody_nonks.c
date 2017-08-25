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
	//PAU int i, j, k, nstar;
	//PAU double r[3], *m, *fm, val;
	int i, j, k, nstar, PN1, PN2, PN25, PN3, PN35;
	double n[3],r[3],v[3], *m, *fm, *fmr, val, rdot, SM, nu, A, B, A2, B2, A4, B4, A5, B5,  A6, B6, A7, B7;
	double clight, clight2, clight4, clight5, R, Mclus, aclus;
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

			SM = m[i] + m[j];
			nu = m[i]*m[j]/(SM*SM);
			for (k=0; k<3; k++) {
				r[k] = y[j*6+k] - y[i*6+k];	
				v[k] = y[j*6+k+3] - y[i*6+k+3]; 
			}
			
			/* this is not necessary, since collisions are handled properly elsewhere
			   in the code */
			/* if(fb_mod(r)<2*SM/clight2){ */
/* 				collisions++; */
/* 				if(collisions==1)printf("%dst collision at %6.16g\n",collisions,t); */
/* 				if(collisions==2)printf("%dnd collision at %6.16g\n",collisions,t); */
/* 				if(collisions==3)printf("%drd collision at %6.16g\n",collisions,t); */
/* 				if(collisions>3)printf("%dth collision at %6.16g\n",collisions,t); */
/* 			  if(nstar-collisions<=1) exit(-1);	 */
/* 			} */
			
			for (k=0; k<3; k++) {
				n[k] = r[k]/fb_mod(r); 
			}
			
			rdot = 0;
			for (k=0; k<3; k++) {
				rdot += n[k]*v[k]; 
			}
			
			if (PN1) {
				A2 = (-3*rdot*rdot*nu/2 + fb_sqr(fb_mod(v)) + 3*nu*fb_sqr(fb_mod(v)) - \
				      SM*(4+2*nu)/fb_mod(r))/clight2;
				B2 = (-4 + 2*nu)*rdot/clight2;
			} else {
				A2 = 0.0;
				B2 = 0.0;
			}

			if (PN2) {
				A4 = (15*rdot*rdot*rdot*rdot*nu/8-45*rdot*rdot*rdot*rdot*nu*nu/8-\
				      9*rdot*rdot*nu*fb_sqr(fb_mod(v))/2+6*rdot*rdot*nu*nu*fb_sqr(fb_mod(v))+\
				      3*nu*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))-4*nu*nu*fb_sqr(fb_mod(v))*\
				      fb_sqr(fb_mod(v))+SM*(-2*rdot*rdot-25*rdot*rdot*nu-2*rdot*rdot*nu*nu-\
							    13*nu*fb_sqr(fb_mod(v))/2+2*nu*nu*fb_sqr(fb_mod(v)))\
				      /fb_mod(r)+SM*SM*(9+87*nu/4)/fb_sqr(fb_mod(r)))/clight4;
				B4 = (9*rdot*rdot*rdot*nu/2+3*rdot*rdot*rdot*nu*nu-15*rdot*nu*fb_sqr(fb_mod(v))/\
				      2-2*rdot*nu*nu*fb_sqr(fb_mod(v))+SM*(2*rdot+41*rdot*nu/2+4*rdot*nu*nu)/\
				      fb_mod(r))/clight4;
			} else {
				A4 = 0.0;
				B4 = 0.0;
			}

			if (PN25) {
				A5 = (-24*rdot*nu*fb_sqr(fb_mod(v))*SM/(5*fb_mod(r))-\
				      136*rdot*nu*SM*SM/(15*fb_sqr(fb_mod(r))))/clight5;
				B5 = (8*nu*fb_sqr(fb_mod(v))*SM/(5*fb_mod(r)) + \
				      24*nu*SM*SM/(5*fb_sqr(fb_mod(r))))/clight5;
			} else {
				A5 = 0.0;
				B5 = 0.0;
			}
			
			if (PN3) {
				A6 = -((16+(1399/12-41*pi2/16)*nu+35.5*nu*nu)*SM*SM*SM/fb_mod(r)/fb_mod(r)/fb_mod(r)+\
				       nu*fb_sqr(fb_mod(v))*SM*SM*(20827/840+123*pi2/64-nu*nu)/fb_sqr(fb_mod(r)) - \
				       rdot*rdot*SM*SM*(1+(22717/168+615/64*pi2)*nu+11*nu*nu/8-7*nu*nu*nu)/fb_sqr(fb_mod(r)) - \
				       0.25*nu*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*(11-49*nu+52*nu*nu) + \
				       35*rdot*rdot*rdot*rdot*rdot*rdot*nu*(1-5*nu+5*nu*nu)/16 -\
				       0.25*nu*SM*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*(75+32*nu-40*nu*nu)/fb_mod(r) - \
				       0.5*nu*rdot*rdot*rdot*rdot*SM*(158-69*nu-60*nu*nu)/fb_mod(r) + \
				       nu*SM*rdot*rdot*fb_sqr(fb_mod(v))*(121-16*nu-20*nu*nu)/fb_mod(r)  + \
				       3*nu*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*rdot*rdot*(20-79*nu+60*nu*nu)/8 -\
				       15*nu*rdot*rdot*rdot*rdot*fb_sqr(fb_mod(v))*(4-18*nu+17*nu*nu)/8 )/clight4/clight2;
				B6 = -rdot*((4+((5849/840)+(123/32)*pi2)*nu-25*nu*nu-8*nu*nu*nu)*SM*SM/\
					    fb_sqr(fb_mod(r))+nu*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*(65-152*nu-48*nu*nu)/8+\
					    15*nu*rdot*rdot*rdot*rdot*(3-8*nu-2*nu*nu)/8+nu*(15+27*nu+10*nu*nu)*fb_sqr(fb_mod(v))*\
					    SM/fb_mod(r)-nu*SM*rdot*rdot*(329+177*nu+108*nu*nu)/fb_mod(r)/6-\
					    0.75*nu*rdot*rdot*fb_sqr(fb_mod(v))*(16-37*nu-16*nu*nu))/clight4/clight2;
			} else {
				A6 = 0.0;
				B6 = 0.0;
			}

			if (PN35) {
				A7 =  1.6*nu*SM*rdot*(23*SM*SM*(43+14*nu)/fb_sqr(fb_mod(r))/14+\
						      3*fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*(61+70*nu)/\
						      28+70*rdot*rdot*rdot*rdot+SM*fb_sqr(fb_mod(v))*(519-1267*nu)/42/fb_mod(r) +\
						      SM*rdot*rdot*(147+188*nu)/4/fb_mod(r)-\
						      15*rdot*rdot*fb_sqr(fb_mod(v))*(19+2*nu)/4)/fb_mod(r)/clight5/clight2;
				B7 = -1.6*nu*SM*(SM*SM*(1325+546*nu)/fb_sqr(fb_mod(r))/42+\
						 fb_sqr(fb_mod(v))*fb_sqr(fb_mod(v))*(313+42*nu)/28+75*rdot*rdot*rdot*rdot-\
						 SM*fb_sqr(fb_mod(v))*(205+777*nu)/fb_mod(r)/42   +\
						 SM*rdot*rdot*(205+424*nu)/fb_mod(r)/12-\
						 0.75*rdot*rdot*fb_sqr(fb_mod(v))*(113+2*nu))/fb_mod(r)/clight5/clight2;
			} else {
				A7 = 0.0;
				B7 = 0.0;
			}

			A = A2 + A4 + A5 + A6 + A7;
			B = B2 + B4 + B5 + B6 + B7;

			/* fprintf(stdout, "%g\n", fb_mod(v)/clight); */

			val = 1.0 / fb_cub(fb_mod(r));

			for (k=0; k<3; k++) {
				FB_FM(i, j, k) = val * r[k];
				FB_REL(i, j, k) = (A*n[k]+B*v[k]) / fb_sqr(fb_mod(r)) / SM;
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
				// PAU f[i*6+k+3] += m[j] * FB_FM(i, j, k);
				f[i*6+k+3] += m[j] * FB_FM(i, j, k) + m[j] * FB_REL(i,j,k);
			}
		}

		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				// PAU f[i*6+k+3] += m[j] * FB_FM(i, j, k);
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
	for (j=0; j< (unsigned int) nstar*6; j++) {
		dfdt[j] = 0.0;
	}

	/* set the matrix to zero to begin with */
	gsl_matrix_set_zero(matrix);

	/* then set the actual values */
	for (i=0; i< (unsigned int) nstar; i++) {
		for (a=0; a<3; a++) {
			gsl_matrix_set(matrix, i+a, i+a+3, 1.0);

			for (k=0; k< (unsigned int) nstar; k++) {
				for (b=0; b<3; b++) {
					val = 0.0;
					for (j=0; j< (unsigned int) nstar; j++) {
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
