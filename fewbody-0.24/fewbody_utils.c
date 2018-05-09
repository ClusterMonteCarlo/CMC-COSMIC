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

/* same, but for the relativistic internal nergy */
double fb_einttot_rel(fb_obj_t *star, int nstar)
{
	int i;
	double eint=0.0;

	for (i=0; i<nstar; i++) {
		eint += star[i].Eint_rel;
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
/* calculates the energy lost to gravitational waves this timestep */ 
double fb_dedt_gw(fb_obj_t *star, int nstar, fb_units_t units, fb_nonks_params_t nonks_params)
{
	/* If no 2.5 pN terms, no need */
	if(nonks_params.PN25 == 0)
		return(0.0);

	int i,j,k;
	double dedt=0.0,temp=0.0;
	double clight = FB_CONST_C / units.v;
	double clight5 = fb_cub(clight)*fb_sqr(clight);
	double m1,m2,r,v,m,mu,nu,rdot,v12v12,n12v12;
	double r_vec[3],n[3],v_vec[3];


	for (i=0; i<nstar; i++) {
		for (j=i+1; j<nstar; j++) {
				for (k=0; k<3; k++) {
					r_vec[k] = (star[i].x[k] - star[j].x[k]) * units.l * (FB_CONST_C*FB_CONST_C/FB_CONST_G) / FB_CONST_MSUN;
					v_vec[k] = (star[i].v[k] - star[j].v[k]) * units.v / FB_CONST_C;
					n[k] = r_vec[k];
				}
				m = (star[i].m + star[j].m) * units.m / FB_CONST_MSUN;
				mu = (star[i].m * star[j].m)/(star[i].m + star[j].m) * units.m / FB_CONST_MSUN;
				nu = mu / m;
				r = fb_mod(r_vec);
				v = fb_mod(v_vec);
				n[0] /= r; n[1] /= r; n[2] /= r;
				rdot = fb_dot(n,v_vec);

				temp = -8. *mu* nu * m*m*m * (4*v*v - 11.*rdot*rdot/3) / (r*r*r*r) / 5.; 
				dedt += temp;

		}
	}

	return dedt*fb_cub(FB_CONST_C)*fb_sqr(FB_CONST_C/units.v)*units.t/(FB_CONST_G*units.m);
}


void fb_n_ecc(fb_obj_t *obj1, fb_obj_t *obj2, double *a, double *e, fb_units_t units)
{
	double x1[3],x2[3],v1[3],v2[3],xrel[3],vrel[3],l1[3],l2[3],l[3],L[3],A[3];
    double x_com[3], v_com[3];
	double E; 
	int i;

	for (i=0; i<3; i++) {
        x_com[i] = (obj1->m * obj1->x[i] + obj2->m * obj2->x[i]) / (obj1->m + obj2->m);
        v_com[i] = (obj1->m * obj1->v[i] + obj2->m * obj2->v[i]) / (obj1->m + obj2->m);
		x1[i] = obj1->x[i] - x_com[i];
		x2[i] = obj2->x[i] - x_com[i];
		v1[i] = obj1->v[i] - v_com[i];
		v2[i] = obj2->v[i] - v_com[i];
		xrel[i] = obj1->x[i] - obj2->x[i];
		vrel[i] = obj1->v[i] - obj2->v[i];
	}

	/* compute the orbital energy and semi-major axis at merger*/
	E = 0.5 * (obj1->m * fb_dot(v1,v1) + obj2->m * fb_dot(v2,v2)) - 
		obj1->m * obj2->m/fb_mod(xrel);
	*a = -obj1->m * obj2->m / (2.0 * E);

	/* set internal angular momentum */
	fb_cross(x1, v1, l1);
	fb_cross(x2, v2, l2);

	/* compute angular momenta for LRL vector*/
	for (i=0; i<3; i++) {
		L[i] = obj1->m * l1[i] + obj2->m * l2[i];
		l[i] = L[i] * (obj1->m + obj2->m)/(obj1->m * obj2->m);
	}
	
	/* -A = l x v + G M \hat r */
	fb_cross(vrel, l, A);
	for (i=0; i<3; i++) {
		A[i] -= (obj1->m + obj2->m) * xrel[i]/fb_mod(xrel);
	}
	
	/* magnitude of A gives the eccentricity at merger*/
	*e = fb_mod(A)/(obj1->m + obj2->m);
}

void fb_check_ecc_for_inspiral(fb_obj_t *obj1, fb_obj_t *obj2, double sep_M, fb_units_t units)
{
    double a, e, last_sep;

    /* if the last seperation from these two is the same, then its the same binary */
	if(obj1->last_sepM == obj2->last_sepM) {
		last_sep = obj1->last_sepM;
    } else { /* otherwise update it and stop */
        obj1->last_sepM = sep_M;
        obj2->last_sepM = sep_M;

        return;
    }

    /* update the last seperation */
	obj1->last_sepM = sep_M;
	obj2->last_sepM = sep_M;
		
    /* update the last a,e at 100M */
	if (last_sep > 100 && sep_M < 100){
        fb_n_ecc(obj1,obj2,&a,&e,units);

		obj1->e_100M[obj1->ncoll - 1] = e; 
		obj2->e_100M[obj1->ncoll - 1] = e;
		obj1->a_100M[obj1->ncoll - 1] = a; 
		obj2->a_100M[obj1->ncoll - 1] = a;
	}

    /* update the last a,e at 500M */
	if (last_sep > 50 && sep_M < 50){
        fb_n_ecc(obj1,obj2,&a,&e,units);

		obj1->e_50M[obj1->ncoll - 1] = e; 
		obj2->e_50M[obj1->ncoll - 1] = e;
		obj1->a_50M[obj1->ncoll - 1] = a; 
		obj2->a_50M[obj1->ncoll - 1] = a;
	}

    /* update the last a,e at 500M */
	if (last_sep > 500 && sep_M < 500){
        fb_n_ecc(obj1,obj2,&a,&e,units);

		obj1->e_500M[obj1->ncoll - 1] = e; 
		obj2->e_500M[obj1->ncoll - 1] = e;
		obj1->a_500M[obj1->ncoll - 1] = a; 
		obj2->a_500M[obj1->ncoll - 1] = a;
	}

}

/* computes the pair-wise energy between two stars in the center-of-mass frame*/ 
double fb_pn_ecc_t(fb_obj_t *star1, fb_obj_t *star2, fb_units_t units, fb_nonks_params_t nonks_params)
{
	double m, mu, nu, e, h, e2;
	double pi2 = FB_CONST_PI*FB_CONST_PI;
	double nu2,nu3;

	m = (star1->m + star2->m)*units.m / FB_CONST_MSUN;
	mu = (star1->m * star2->m)/(star1->m + star2->m) * units.m / FB_CONST_MSUN;
	nu = mu / m;
	nu2 = nu*nu;
	nu3 = nu2*nu;

	e = E_rel(star1, star2, units, nonks_params)*m/mu;
	h = J_rel(star1, star2, units, nonks_params)*m*m/mu;

	if(e > 0.)
		return -1;

	e = fabs(e);
	e2 = 1-2*e*h*h;
	if (nonks_params.PN1 == 1)
		e2 += - 4*(1-nu)*e + (17-7*nu)*e*e*h*h; 
	if (nonks_params.PN2 == 1)
		e2 += 2*(2+nu+5*nu*nu)*e*e + (17-11*nu)*e/h/h - (112  -47*nu + 16*nu*nu)*e*e*e*h*h 
			+ 3*(2*nu-5)*(1-2*e*h*h)*pow(2*e,1.5)/h; 
	if (nonks_params.PN3 == 1)
		e2 += (8*e*e*e/6720.) * (23520. - 464800*nu + 179760*nu2 + 16800*nu3 
			- 2520.*sqrt(2*e*h*h)*(265-193*nu+46*nu2) 
			- 525*(2*e*h*h)*(-528. + 200*nu - 77*nu2 + 24*nu3) 
			- 6*(73920 - 260272*nu + 4305*pi2*nu + 61040*nu2)/(2*e*h*h)
			+ 70*(16380 - 19964*nu + 123*pi2*nu + 3240*nu2) / sqrt(2*e*h*h)
			+ 8*(53760 - 176024*nu + 4305*pi2*nu + 15120*nu2) / (4*e*e*h*h*h*h)
			- 70*(10080 - 13952*nu + 123*pi2*nu + 1440*nu2) / pow(2*e*h*h,1.5));

	if(e2 >= 0)
		return sqrt(e2);
	else	
		return 0;
}


/* Computes the relativistic energy in fewbody code units*/
double fb_E_rel(fb_obj_t *star1, fb_obj_t *star2, fb_units_t units, fb_nonks_params_t nonks_params)
{
	double mass,energy,mu;
	mass = (star1->m + star2->m)*units.m / FB_CONST_MSUN;
	mu = (star1->m * star2->m)/(star1->m + star2->m) * units.m / FB_CONST_MSUN;
	energy = E_rel(star1, star2, units, nonks_params) *mass;

	return energy*fb_sqr(FB_CONST_C/units.v)*(FB_CONST_MSUN/units.m);
}

/* computes the pair-wise energy between two stars in the center-of-mass frame 
 *
 * Returns in dimmensionless units (E/M), where M is the total mass of the 
 * pair*/ 
double E_rel(fb_obj_t *star1, fb_obj_t *star2, fb_units_t units, fb_nonks_params_t nonks_params)
{
	double r_vec[3],v_vec[3],n[3],r,v,m,mu,nu,rdot;
	double rdot2,rdot4,rdot6,v2,v4,v6,v8,nu2,nu3;
	double energy = 0.;
	int k;

	double clight = FB_CONST_C / units.v;
	double clight2 = fb_sqr(clight);
	double clight4 = fb_sqr(clight2);
	double pi2 = FB_CONST_PI*FB_CONST_PI;

	for (k=0; k<3; k++) {
		r_vec[k] = (star1->x[k] - star2->x[k]) * units.l * (FB_CONST_C*FB_CONST_C/FB_CONST_G) / FB_CONST_MSUN;
		v_vec[k] = (star1->v[k] - star2->v[k]) * units.v / FB_CONST_C;
		n[k] = r_vec[k];
	}

	m = (star1->m + star2->m)*units.m / FB_CONST_MSUN;
	mu = (star1->m * star2->m)/(star1->m + star2->m) * units.m / FB_CONST_MSUN;
	nu = mu / m;
	nu2 = nu*nu;
	nu3 = nu2*nu; 
	r = fb_mod(r_vec);
	v = fb_mod(v_vec);
	v2 = v*v;
	v4 = v2*v2;
	v6 = v4*v2;
	v8 = v4*v4;
	n[0] /= r; n[1] /= r; n[2] /= r;
	rdot = fb_dot(v_vec,n);
	rdot2 = rdot*rdot;
	rdot4 = rdot2*rdot2;
	rdot6 = rdot4*rdot2;

	energy += (v*v / 2) - (m / r); 
	if (nonks_params.PN1 == 1)
		energy +=  (3*v4/8 - 9*nu*v4/8 + 
				m*(rdot2*nu/2 + 3*v2/2 + nu*v2/2)/r + m*m/(2*r*r));
	
	if (nonks_params.PN2 == 1)
		energy +=  (5*v6/16 - 35*nu*v6/16 + 65*nu2*v6/16 + m*(-3*rdot4*nu/8
					+ 9*rdot4*nu2/8 + rdot2*nu*v2/4 - 15*rdot2*nu2*v2/4 
					+ 21*v4/8 - 23*nu*v4/8 - 27*nu2*v4/8)/r + m*m*(rdot2/2 
					+ 69*rdot2*nu/8 + 3*rdot2*nu2/2 + 7*v2/4 - 55*nu*v2/8 
					+ nu2*v2/2)/(r*r) + m*m*m*(-15*nu/4 - 0.5)/(r*r*r));

	if (nonks_params.PN3 == 1)
		energy += (35 - 413*nu + 1666*nu2 - 2261*nu3)*v8/128. + m*(55 - 215*nu + 116*nu2 
					+ 325*nu3)*v6/r/16. + m*(5*nu - 25*nu2 + 25*nu3)*rdot6/r/16.
					- m*(21*nu + 75*nu2 - 375*nu3)*v4*rdot2/r/16. - m*(9*nu - 84*nu2
					+ 165*nu3)*v2*rdot4/r/16. + m*m*(135 - 194*nu + 406*nu2 - 108*nu3)
					*v4/r/r/16. + m*m*(12 + 248*nu - 815*nu2 - 324*nu3)*v2*rdot2/r/r/16.
					- m*m*(731*nu - 492*nu2 - 288*nu3)*rdot4/r/r/48. + m*m*m*(2800 - (53976 - 
					1435*pi2)*nu - 11760*nu2 + 1120*nu3)*v2/r/r/r/2240. + m*m*m*(3360 + (
					18568 - 4305*pi2)*nu + 28560*nu2 + 7840*nu3)*rdot2/r/r/r/2240.
					+ m*m*m*m*(315 + 18469*nu)/r/r/r/r/840.;
					
	energy *= mu;
	return energy / m ;// * FB_CONST_MSUN / units.m ;

}

/* computes the pair-wise angular momentum between two stars in the center-of-mass frame
 *
 * Returns in dimmensionless units (J/M^2), where M is the total mass of the 
 * pair*/ 
double J_rel(fb_obj_t *star1, fb_obj_t *star2, fb_units_t units, fb_nonks_params_t nonks_params)
{
	double r_vec[3],v_vec[3],n[3],z[3],r,v,m,nu,mu,rdot;
	double rdot2,rdot4,v2,v4,v6;
	double J = 0.;
	int k;

	double clight = FB_CONST_C / units.v;
	double clight2 = fb_sqr(clight);
	double clight4 = fb_sqr(clight2);
	double pi2 = FB_CONST_PI*FB_CONST_PI;
	double nu2,nu3;

	for (k=0; k<3; k++) {
		r_vec[k] = (star1->x[k] - star2->x[k]) * units.l * (FB_CONST_C*FB_CONST_C/FB_CONST_G) / FB_CONST_MSUN;
		v_vec[k] = (star1->v[k] - star2->v[k]) * units.v / FB_CONST_C;
		n[k] = r_vec[k];
	}
	
	m = (star1->m + star2->m)*units.m / FB_CONST_MSUN;
	mu = (star1->m * star2->m)/(star1->m + star2->m) * units.m / FB_CONST_MSUN;
	nu = mu / m;
	nu2 = nu*nu;
	nu3 = nu2*nu;
	r = fb_mod(r_vec);
	v = fb_mod(v_vec);
	v2 = v*v;
	v4 = v2*v2;
	v6 = v4*v2;
	n[0] /= r; n[1] /= r; n[2] /= r;
	rdot = fb_dot(v_vec,n);
	rdot2 = rdot*rdot;
	rdot4 = rdot2*rdot2;

	J += 1.; 
	if (nonks_params.PN1 == 1)
		J +=  ((1-3*nu)*v2/2 + m*(3+nu)/r);
	
	if (nonks_params.PN2 == 1)
		J += (3*v4/8 - 21*nu*v4/8 + 39*nu*nu*v4/8 + m*(-rdot2*nu - 5*rdot2*nu*nu/2
			   + 7*v2/2 - 5*nu*v2 - 9*nu*nu*v2/2)/r + m*m*(3.5 - 41*nu/4 + nu*nu)/(r*r));

	if (nonks_params.PN3 == 1)
		J += (5./2 - (5199./280. - 41.*pi2/32.)*nu - 7*nu2 + nu3) * (m*m*m/r/r/r)
			+ (5 - 59*nu + 238*nu2 - 323*nu3)*v6/16.
			+ (135 - 322*nu + 315*nu2 - 108*nu3)*(m*m/r/r)*v2/12.
			+ (12 - 287*nu - 951*nu2 - 324*nu3)*(m*m/r/r)*rdot2/24.
			+ (33 - 142*nu + 106*nu2 + 195*nu3)*m*v4/r/8.
			- nu*(12 - 7*nu - 75*nu2)*m*v2*rdot2/r/4
			+ 3*nu*(2 - 2*nu - 11*nu2)*m*rdot4/r/8.;

	fb_cross(r_vec,v_vec,z);
	J *= fb_mod(z)*mu;

	return J / m / m;
}

/*Computes the pairwise seperation of two particles in units of the total mass*/
double fb_compute_distance_in_M(fb_obj_t *star1, fb_obj_t *star2, fb_units_t units)
{
	int k;
	double r[3], m, rm;

	for (k=0; k<3; k++){
		r[k] = (star1->x[k] - star2->x[k]) * units.l * (FB_CONST_C*FB_CONST_C/FB_CONST_G);
	}

	m = (star1->m + star2->m)*units.m;
	rm = fb_mod(r) / m;

	return(rm);
}

/* calculates the total energy of the system, including relativistic terms */ 

double fb_Etot_rel(fb_obj_t *star, int nstar, fb_units_t units, fb_nonks_params_t nonks_params)
{
	double energy=0., mass, temp;
	int i,j;
	for (i=0; i<nstar; i++) {
		for (j=i+1; j<nstar; j++) {
			temp = fb_E_rel(&star[i], &star[j], units, nonks_params);
			energy += temp; 
		}
	}

	return energy;
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
