/* -*- linux-c -*- */
/* fewbody_int.c

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
#include <gsl/gsl_rng.h>
#include "fewbody.h"

/* allocate memory for ks_params */
void fb_malloc_ks_params(fb_ks_params_t *ks_params)
{
	ks_params->m = fb_malloc_vector(ks_params->nstar);
	ks_params->M = fb_malloc_vector(ks_params->kstar);
	ks_params->amat = fb_malloc_matrix(ks_params->nstar, ks_params->kstar);
	ks_params->Tmat = fb_malloc_matrix(ks_params->kstar, ks_params->kstar);
}

/* initialize ks_params; assumes ks_params is already malloc()ed */
void fb_init_ks_params(fb_ks_params_t *ks_params, fb_hier_t hier)
{
	int i, j, k;
	double *y;

	/* exit if hier is not consistent with ks_params */
	if (ks_params->nstar != hier.nobj) {
		fprintf(stderr, "fb_init_ks_params(): ks_params->nstar != hier.nobj: ks_params->nstar=%d  hier.nobj=%d\n", \
			ks_params->nstar, hier.nobj);
		exit(1);
	}

	/* first set the mass matrix */
	for (i=0; i<hier.nobj; i++) {
		ks_params->m[i] = hier.obj[i]->m;
	}
	
	/* calculate the M_k */
	k = -1;
	for (i=0; i<hier.nobj-1; i++) {
		for (j=i+1; j<hier.nobj; j++) {
			k++;
			ks_params->M[k] = ks_params->m[i] * ks_params->m[j];
		}
	}

	/* calculate and set the a and T matrices */
	fb_calc_amat(ks_params->amat, ks_params->nstar, ks_params->kstar);
	fb_calc_Tmat(ks_params->amat, ks_params->m, ks_params->Tmat, ks_params->nstar, ks_params->kstar);

	/* set Einit */
	y = fb_malloc_vector(8*ks_params->kstar+1);
	fb_euclidean_to_ks(hier.obj, y, ks_params->nstar, ks_params->kstar);
	ks_params->Einit = fb_ks_Einit(y, *ks_params);

	fb_free_vector(y);
}

/* free memory for ks_params */
void fb_free_ks_params(fb_ks_params_t ks_params)
{
	fb_free_vector(ks_params.m);
	fb_free_vector(ks_params.M);
	fb_free_matrix(ks_params.amat);
	fb_free_matrix(ks_params.Tmat);
}

/* allocate memory for nonks_params */
void fb_malloc_nonks_params(fb_nonks_params_t *nonks_params)
{
	nonks_params->m = fb_malloc_vector(nonks_params->nstar);
}

/* initialize nonks_params; assumes nonks_params is already malloc()ed */
void fb_init_nonks_params(fb_nonks_params_t *nonks_params, fb_hier_t hier)
{
	int i;

	/* exit if hier is not consistent with nonks_params */
	if (nonks_params->nstar != hier.nobj) {
		fprintf(stderr, "fb_init_nonks_params(): nonks_params->nstar != hier.nobj: nonks_params->nstar=%d  hier.nobj=%d\n", \
			nonks_params->nstar, hier.nobj);
		exit(1);
	}

	/* set the mass vector */
	for (i=0; i<hier.nobj; i++) {
		nonks_params->m[i] = hier.obj[i]->m;
	}
}

/* free memory for ks_params */
void fb_free_nonks_params(fb_nonks_params_t nonks_params)
{
	fb_free_vector(nonks_params.m);
}

