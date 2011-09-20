/* -*- linux-c -*- */
/* fewbody_hier.c

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
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "fewbody.h"

/* allocate memory for a hier_t */
void fb_malloc_hier(fb_hier_t *hier)
{
	int i;

	hier->hi = (int *) malloc(hier->nstarinit * sizeof(int)) - 1;
	fb_create_indices(hier->hi, hier->nstarinit);
	hier->narr = (int *) malloc((hier->nstarinit-1)*sizeof(int)) - 2;
	/* change from malloc to calloc to get rid of harmless valgrind errors...
	   the errors occur in fb_normalize(), where uninitialized memory is normalized - 
	   these are clearly harmless */
	hier->hier = (fb_obj_t *) calloc(hier->hi[hier->nstarinit] + 1, sizeof(fb_obj_t));
	for (i=0; i<hier->hi[hier->nstarinit]+1; i++) {
		hier->hier[i].ncoll = 0;
		hier->hier[i].id = (long *) malloc(hier->nstarinit * sizeof(long));
	}
	hier->obj = (fb_obj_t **) malloc(hier->nstarinit * sizeof(fb_obj_t *));
}

/* initialize to a flat hier */
void fb_init_hier(fb_hier_t *hier)
{
	int i;
	
	hier->nobj = hier->nstar;
	
	for (i=2; i<=hier->nstar; i++) {
		hier->narr[i] = 0;
	}
	
	for (i=0; i<hier->nstar; i++) {
		hier->obj[i] = &(hier->hier[hier->hi[1]+i]);
	}
}

/* free memory */
void fb_free_hier(fb_hier_t hier)
{
	int i;

	for (i=0; i<hier.hi[hier.nstarinit]+1; i++) {
		free(hier.hier[i].id);
	}
	free(hier.hi+1);
	free(hier.narr+2);
	free(hier.hier);
	free(hier.obj);
}

/* trickle down hier */
void fb_trickle(fb_hier_t *hier, double t)
{
	int i, j;
	
	for (i=hier->nstar; i>=2; i--) {
		for (j=0; j<hier->narr[i]; j++) {
			fb_downsync(&(hier->hier[hier->hi[i] + j]), t);
		}
	}
}

/* trickle up hier */
void fb_elkcirt(fb_hier_t *hier, double t)
{
	int i, j;
	
	for (i=2; i<=hier->nstar; i++) {
		for (j=0; j<hier->narr[i]; j++) {
			fb_upsync(&(hier->hier[hier->hi[i] + j]), t);
		}
	}
}

/* created the index array */
int fb_create_indices(int *hi, int nstar)
{
	int i, j=0;
	
	for (i=1; i<=nstar; i++) {
		hi[i] = j;
		j += nstar/i;
	}
	
	return(0);
}

/* how many members are in the immediate hierarchy */
int fb_n_hier(fb_obj_t *obj)
{
	if (obj == NULL) {
		return(0);
	} else if ((obj->obj[0] == NULL) && (obj->obj[1] == NULL)) {
		return(1);
	} else if ((obj->obj[0]->obj[0] == NULL) && (obj->obj[0]->obj[1] == NULL) && (obj->obj[1]->obj[0] == NULL) && (obj->obj[1]->obj[1] == NULL)) {
		return(2);
	} else if (((obj->obj[0]->obj[0] == NULL) && (obj->obj[0]->obj[1] == NULL)) || ((obj->obj[1]->obj[0] == NULL) && (obj->obj[1]->obj[1] == NULL))) {
		return(3);
	} else {
		return(4);
	}
}

/* print the hierarchy information */
char *fb_sprint_hier(fb_hier_t hier, char string[FB_MAX_STRING_LENGTH])
{
	int i;
	
	snprintf(string, FB_MAX_STRING_LENGTH, "nstar=%d nobj=%d: ", hier.nstar, hier.nobj);
	for (i=0; i<hier.nobj; i++) {
		snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), " %s", hier.obj[i]->idstring);
	}
	
	return(string);
}

/* print the hierarchy information in a nice human-readable format */
char *fb_sprint_hier_hr(fb_hier_t hier, char string[FB_MAX_STRING_LENGTH])
{
	int i;

	/* zero out string so strlen(string) returns 0 */
	string[0] = '\0';
	
	/* loop through objects */
	for (i=0; i<hier.nobj; i++) {
		/* put a dash inbetween labels */
		if (i > 0) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "-");
		}
		/* a name for each type of hierarchy */
		if (hier.obj[i]->n == 1) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "single");
		} else if (hier.obj[i]->n == 2) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "binary");
		} else if (hier.obj[i]->n == 3) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "triple");
		} else if (hier.obj[i]->n == 4) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "quadruple");
		} else if (hier.obj[i]->n == 5) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "quintuple");
		} else if (hier.obj[i]->n == 6) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "sextuple");
		} else if (hier.obj[i]->n == 7) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "septuple");
		} else if (hier.obj[i]->n == 8) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "octuple");
		} else if (hier.obj[i]->n == 9) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "nontuple");
		} else if (hier.obj[i]->n == 10) {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "decituple");
		} else {
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "%d", hier.obj[i]->n);	
			snprintf(&(string[strlen(string)]), FB_MAX_STRING_LENGTH-strlen(string), "tuple");
		}
	}

	return(string);
}

/* merge the object's properties up---calculate the binary's properties from the
   stars' properties */
void fb_upsync(fb_obj_t *obj, double t)
{
	int i;
	double m0, m1, x0[3], x1[3], xrel[3], v0[3], v1[3], vrel[3], E, l0[3], l1[3], l[3];
	double ypp[3], psi, ecc_anom, L[3], A[3];

	/* a little bit of paranoia here */
	obj->ncoll = 0;
	obj->id[0] = -1;

	/* set time of upsync */
	obj->t = t;

	/* create idstring */
	snprintf(obj->idstring, FB_MAX_STRING_LENGTH, "[%s %s]", obj->obj[0]->idstring, obj->obj[1]->idstring);

	/* update number of stars in hierarchy */
	obj->n = obj->obj[0]->n + obj->obj[1]->n;
	
	/* update mass */
	m0 = obj->obj[0]->m;
	m1 = obj->obj[1]->m;
	obj->m = m0 + m1;

	/* update position and velocity */
	for (i=0; i<3; i++) {
		obj->x[i] = (m0*(obj->obj[0]->x[i]) + m1*(obj->obj[1]->x[i])) / (obj->m);
		obj->v[i] = (m0*(obj->obj[0]->v[i]) + m1*(obj->obj[1]->v[i])) / (obj->m);
	}
	
	/* update semimajor axis */
	for (i=0; i<3; i++) {
		x0[i] = obj->obj[0]->x[i] - obj->x[i];
		x1[i] = obj->obj[1]->x[i] - obj->x[i];
		xrel[i] = obj->obj[0]->x[i] - obj->obj[1]->x[i];
		v0[i] = obj->obj[0]->v[i] - obj->v[i];
		v1[i] = obj->obj[1]->v[i] - obj->v[i];
		vrel[i] = obj->obj[0]->v[i] - obj->obj[1]->v[i];
	}
	
	E = 0.5*(m0*fb_dot(v0, v0) + m1*fb_dot(v1, v1)) - m0*m1/fb_mod(xrel);
	
	if (E >= 0.0) {
		fprintf(stderr, "fb_upsync: E = %g >= 0!\n", E);
		exit(1);
	}
	
	obj->a = -m0 * m1 / (2.0 * E);

	/* update angular momentum, Runge-Lenz vector, and eccentricity (using the Runge-Lenz vector) */
	fb_cross(x0, v0, l0);
	fb_cross(x1, v1, l1);
	
	for (i=0; i<3; i++) {
		L[i] = m0 * l0[i] + m1 * l1[i];
		l[i] = L[i] * (m0+m1)/(m0*m1);
	}
	
	/* -A = l x v + G M \hat r */
	fb_cross(vrel, l, A);
	for (i=0; i<3; i++) {
		A[i] -= (m0+m1) * xrel[i]/fb_mod(xrel);
	}
	
	obj->e = fb_mod(A)/(m0+m1);
	
	/* must be careful with circular orbits, which have a null A vector */
	if (obj->e == 0.0) {
		for (i=0; i<3; i++) {
			A[i] = xrel[i];
		}
	}
	
	/* define coord system */
	for (i=0; i<3; i++) {
		obj->Lhat[i] = L[i]/fb_mod(L); /* z-direction */
		obj->Ahat[i] = A[i]/fb_mod(A); /* x-direction */
	}
	fb_cross(obj->Lhat, obj->Ahat, ypp);
	
	/* and set mean anomaly */
	psi = atan2(fb_dot(x0, ypp), fb_dot(x0, obj->Ahat));
	ecc_anom = acos((obj->e + cos(psi)) / (1.0 + obj->e * cos(psi)));
	if (psi < 0.0) {
		ecc_anom = 2.0 * FB_CONST_PI - ecc_anom;
	}
	obj->mean_anom = ecc_anom - obj->e * sin(ecc_anom);
}

/* randomly orient binary */
void fb_randorient(fb_obj_t *obj, gsl_rng *rng, struct rng_t113_state *curr_st)
{
	int i;
	double theta, phi, omega, xp[3], yp[3], zp[3];

	/* random angles */
/*
	theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0);
	phi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
	omega = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
	obj->mean_anom = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
*/

	//MPI2: Changing this temporarily? to use global random state.
	theta = acos(2.0 * rng_t113_dbl_new(curr_st) - 1.0);
	phi = 2.0 * FB_CONST_PI * rng_t113_dbl_new(curr_st);
	omega = 2.0 * FB_CONST_PI * rng_t113_dbl_new(curr_st);
	obj->mean_anom = 2.0 * FB_CONST_PI * rng_t113_dbl_new(curr_st);

	/* set up coordinate transformations */
	zp[0] = -sin(theta) * sin(phi);
	zp[1] = sin(theta) * cos(phi);
	zp[2] = cos(theta);
	
	xp[0] = cos(phi);
	xp[1] = sin(phi);
	xp[2] = 0.0;
	
	fb_cross(zp, xp, yp);
	
	for (i=0; i<3; i++) {
		obj->Lhat[i] = zp[i];
		obj->Ahat[i] = xp[i] * cos(omega) + yp[i] * sin(omega);
	}
}

/* merge the object's properties down---calculate the objs' properties from the
   binary's properties */
void fb_downsync(fb_obj_t *obj, double t)
{
	int i;
	double xpp[3], ypp[3], zpp[3], er[3], epsi[3];
	double a, e, m0, m1, omega, mean_anom, ecc_anom, psi, r0, r1, L, psidot, E, r0dot, r0dot2, r1dot;

	/* we're assuming here that the objs' masses are already set */
	a = obj->a;
	e = obj->e;
	m0 = obj->obj[0]->m;
	m1 = obj->obj[1]->m;
	
	/* angular frequency */
	omega = sqrt(obj->m/fb_cub(a));
	/* mean anomaly, between 0 and 2PI */
	mean_anom = (obj->mean_anom + omega * (t - obj->t)) / (2.0 * FB_CONST_PI);
	mean_anom = (mean_anom - floor(mean_anom)) * 2.0 * FB_CONST_PI;
	/* eccentric anomaly, from solving the Kepler equation */
	ecc_anom = fb_kepler(e, mean_anom);
	/* true anomaly, between 0 and 2PI */
	psi = acos((cos(ecc_anom) - e) / (1.0 - e * cos(ecc_anom)));
	/* this step is necessary because acos() returns a value between 0 and PI */
	if (ecc_anom > FB_CONST_PI) {
		psi = 2.0 * FB_CONST_PI - psi;
	}
	
	/* define vectors, based on angular momentum and Runge-Lenz vectors */
	for (i=0; i<3; i++) {
		zpp[i] = obj->Lhat[i];
		xpp[i] = obj->Ahat[i];
	}
	fb_cross(zpp, xpp, ypp);
	
	for (i=0; i<3; i++) {
		er[i] = xpp[i] * cos(psi) + ypp[i] * sin(psi);
		epsi[i] = ypp[i] * cos(psi) - xpp[i] * sin(psi);
	}
	
	/* determine binary params */
	r0 = a * (1.0 - fb_sqr(e)) / ((1.0 + e*cos(psi)) * (1.0 + m0/m1));
	r1 = a * (1.0 - fb_sqr(e)) / ((1.0 + e*cos(psi)) * (1.0 + m1/m0));

	L = m0 * m1 * sqrt((m0+m1)*a*(1.0-fb_sqr(e))) / (m0 + m1);
	psidot = L / (m0 * fb_sqr(r0) + m1 * fb_sqr(r1));
	E = -m0 * m1 / (2.0 * a);

	/* must be careful for circular orbits */
	r0dot2 = (E + m0*m1/(r0+r1)) * 2.0 / (m0 * (1.0+m0/m1)) - fb_sqr(r0*psidot);
	if (ecc_anom > FB_CONST_PI) {
		r0dot = -(r0dot2<=0.0?0.0:sqrt(r0dot2));
	} else {
		r0dot = (r0dot2<=0.0?0.0:sqrt(r0dot2));
	}
	r1dot = m0 * r0dot / m1;
	
	for (i=0; i<3; i++) {
		obj->obj[0]->x[i] = r0 * er[i] + obj->x[i];
		obj->obj[0]->v[i] = r0dot * er[i] + r0 * psidot * epsi[i] + obj->v[i];
		obj->obj[1]->x[i] = -r1 * er[i] + obj->x[i];
		obj->obj[1]->v[i] = -r1dot * er[i] - r1 * psidot * epsi[i] + obj->v[i];
	}
}

/* copy one object to another, being careful about any pointers */
void fb_objcpy(fb_obj_t *obj1, fb_obj_t *obj2)
{
	int i;
	long *longptr;
	
	for (i=0; i<obj2->ncoll; i++) {
		obj1->id[i] = obj2->id[i];
	}

	/* prevent any dangling pointers */
	longptr = obj1->id;
	*obj1 = *obj2;
	obj1->id = longptr;
}
