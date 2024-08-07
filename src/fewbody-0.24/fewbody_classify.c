/* -*- linux-c -*- */
/* fewbody_classify.c

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
#include "fewbody.h"

/* classify the stars into hierarchies; i.e., build the binary tree */
// PAU int fb_classify(fb_hier_t *hier, double t, double tidaltol)
int fb_classify(fb_hier_t *hier, double t, double tidaltol, double speedtol, fb_units_t units, fb_input_t input)
{
	int i, j, k, n, isave[2], cont=1;
	double a, amin, E, xrel[3], v0[3], v1[3], vcm[3], vrel[3], ftid;

	/* initialize to flat hier */
	fb_init_hier(hier);

	/* first build the hierarchy */
	while (cont) {
		amin = FB_AMIN;
		cont = 0;
		for (j=0; j<hier->nobj; j++) {
			for (k=j+1; k<hier->nobj; k++) {
				for (i=0; i<3; i++) {
					xrel[i] = hier->obj[j]->x[i] - hier->obj[k]->x[i];
					vcm[i] = (hier->obj[j]->m * hier->obj[j]->v[i] + hier->obj[k]->m * hier->obj[k]->v[i]) /\
						(hier->obj[j]->m + hier->obj[k]->m);
					v0[i] = hier->obj[j]->v[i] - vcm[i];
					v1[i] = hier->obj[k]->v[i] - vcm[i];
				}
				
				E = 0.5 * (hier->obj[j]->m * fb_dot(v0, v0) + hier->obj[k]->m * fb_dot(v1, v1)) - \
					hier->obj[j]->m * hier->obj[k]->m / fb_mod(xrel);
				
				if (E < 0.0) {
					a = -hier->obj[j]->m * hier->obj[k]->m / (2.0 * E);
					if (a < amin) {
						amin = a;
						cont = 1;
						isave[0] = j;
						isave[1] = k;
					}
				}
			}
		}
		
		/* did we find a binary? */
		if (cont) {
			/* swap indices so object on the left contains more stars, or, if
			   they are each single, so object on left has smaller id */
			if (hier->obj[isave[0]]->n == 1 && hier->obj[isave[1]]->n == 1) {
				if (hier->obj[isave[0]]->id[0] > hier->obj[isave[1]]->id[0]) {
					j = isave[0];
					isave[0] = isave[1];
					isave[1] = j;
				}
			} else {
				if (hier->obj[isave[0]]->n < hier->obj[isave[1]]->n) {
					j = isave[0];
					isave[0] = isave[1];
					isave[1] = j;
				}
			}

			/* create new hierarchical object */
			n = hier->obj[isave[0]]->n + hier->obj[isave[1]]->n;
			hier->hier[hier->hi[n]+hier->narr[n]].obj[0] = hier->obj[isave[0]];
			hier->hier[hier->hi[n]+hier->narr[n]].obj[1] = hier->obj[isave[1]];
			fb_upsync(&(hier->hier[hier->hi[n]+hier->narr[n]]), t);
			hier->obj[isave[0]] = &(hier->hier[hier->hi[n]+hier->narr[n]]);
			hier->narr[n]++;

			/* delete extraneous object */
			hier->obj[isave[1]] = hier->obj[hier->nobj-1];
			hier->nobj--;
		}
	}
	
	/******************************/
	/* now start performing tests */
	/******************************/

	/* make sure all pairs have positive relative velocities */
	for (j=0; j<hier->nobj; j++) {
		for (k=j+1; k<hier->nobj; k++) {
			for (i=0; i<3; i++) {
				xrel[i] = hier->obj[j]->x[i] - hier->obj[k]->x[i];
				vrel[i] = hier->obj[j]->v[i] - hier->obj[k]->v[i];
			}
			
			if (fb_dot(vrel, xrel)/fb_mod(xrel) < 0.0) {
				fb_dprintf("fewbody: classify(): not done: vrel < 0\n");
				return(0);
			}
		}
	}

	/* check tidal perturbations for the unbound objects */
	for (j=0; j<hier->nobj; j++) {
		if (hier->obj[j]->n > 1) {
			/* sum over other objects to get relative tide */
			ftid = 0.0;
			for (k=0; k<hier->nobj; k++) {
				if (k != j) {
					for (i=0; i<3; i++) {
						xrel[i] = hier->obj[j]->x[i] - hier->obj[k]->x[i];
					}
					ftid += fb_reltide(hier->obj[j], hier->obj[k], fb_mod(xrel));
				}
			}
			if (ftid >= tidaltol) {
				fb_dprintf("fewbody: classify(): unbound: tidal tolerance exceeded; encounter not done.\n");
				return(0);
			}			
		}
	}
	
	/* check stability of bound hierarchies */
	for (i=2; i<=hier->nstar; i++) {
		for (j=0; j<hier->narr[i]; j++) {
			// if (!fb_is_stable(&(hier->hier[hier->hi[i]+j]))) {
			if (!fb_is_stable(&(hier->hier[hier->hi[i]+j]), speedtol, units, input)) {
				fb_dprintf("fewbody: classify(): unstable hierarchy: i=%d hier->narr[i]=%d j=%d\n", 
					   i, hier->narr[i], j);
				return(0);
			}
		}
	}
	
	/* it is possible that the system has passed all the above tests, but its
	   outer objects are bound */
	if (fb_outerketot(hier->obj, hier->nobj) + fb_outerpetot(hier->obj, hier->nobj) < 0.0) {
		fb_dprintf("fewbody: classify(): outer objects are bound; encounter not done.\n");
		return(0);
	}

	/* all done! */
	fb_dprintf("fewbody: classify(): encounter done!\n");
	return(1);
}

/* check the stability of an arbitrary hierarchical object */
// PAU int fb_is_stable(fb_obj_t *obj)
int fb_is_stable(fb_obj_t *obj, double speedtol, fb_units_t units, fb_input_t input)
{
	if (fb_n_hier(obj) == 2) {
		// PAU return(fb_is_stable_binary(obj));
		return(fb_is_stable_binary(obj, speedtol, units, input));
	} else if (fb_n_hier(obj) == 3) {
		return(fb_is_stable_triple(obj));
	} else if (fb_n_hier(obj) == 4) {
		return(fb_is_stable_quad(obj));
	} else {
		return(1);
	}
}

int fb_is_stable_binary(fb_obj_t *obj, double speedtol, fb_units_t units, fb_input_t input)
{
	double vrelperi;
	double clight;

	clight = FB_CONST_C / units.v;
	vrelperi = sqrt(2.0 * obj->m / (obj->a * (1.0 - obj->e)) - obj->m / obj->a);

	/* test for collision at pericenter */
	// PAU if (fb_is_collision(obj->a * (1.0 - obj->e), obj->obj[0]->R, obj->obj[1]->R)) {
	if (fb_is_collision(obj->a * (1.0 - obj->e), obj->obj[0]->R, obj->obj[1]->R, obj->obj[0]->m, obj->obj[1]->m, obj->obj[0]->k_type, obj->obj[1]->k_type, units.m, units.l, input.BHNS_TDE_FLAG) || 
		vrelperi / clight >= speedtol) {
		return(0);
	} else {
		return(1);
	}
}

int fb_is_stable_triple(fb_obj_t *obj)
{
	if (fb_n_hier(obj->obj[0]) == 1) {
		return(fb_mardling(obj, 1, 0));
	} else {
		return(fb_mardling(obj, 0, 1));
	}
}

int fb_is_stable_quad(fb_obj_t *obj)
{
	int ib, is;
	
	/* use the wider binary as the inner binary */
	if (obj->obj[0]->a >= obj->obj[1]->a) {
		ib = 0;
		is = 1;
	} else {
		ib = 1;
		is = 0;
	}

	if (fb_mardling(obj, ib, is)) {
		return(1);
	} else {
		return(0);
	}
}

/* the mardling criterion for the stability of triples or quadruples */
int fb_mardling(fb_obj_t *obj, int ib, int is)
{
	int i;
	double C=2.8, ain, aout, qout, eout, Rpout, inc, l0[3], l1[3], Lout[3], Lin[3], f;
	double a2, fquad;

	/* set useful variables */
	ain = obj->obj[ib]->a;
	qout = obj->obj[is]->m / obj->obj[ib]->m;
	aout = obj->a;
	eout = obj->e;
	Rpout = aout * (1.0 - eout);

	/* calculate the extra ad hoc factor for quadruples if this is a quad */
	if ((obj->obj[is]->obj[0] != NULL) && (obj->obj[is]->obj[1] != NULL)) {
		a2 = obj->obj[is]->a;
		fquad = 0.1 * a2/ain;
	} else {
		fquad = 0.0;
	}

	/* calculate inclination */
	fb_cross(obj->obj[0]->x, obj->obj[0]->v, l0);
	fb_cross(obj->obj[1]->x, obj->obj[1]->v, l1);
	
	for (i=0; i<3; i++) {
		Lout[i] = obj->obj[0]->m * l0[i] + obj->obj[1]->m * l1[i];
	}

	fb_cross(obj->obj[ib]->obj[0]->x, obj->obj[ib]->obj[0]->v, l0);
	fb_cross(obj->obj[ib]->obj[1]->x, obj->obj[ib]->obj[1]->v, l1);
	
	for (i=0; i<3; i++) {
		Lin[i] = obj->obj[ib]->obj[0]->m * l0[i] + obj->obj[ib]->obj[1]->m * l1[i];
	}

	/* acos returns a value between 0 and PI */
	inc = acos(fb_dot(Lin, Lout)/(fb_mod(Lin)*fb_mod(Lout)));
	
	/* may be unconditionally unstable due to mass ratio */
	if (qout > 5.0) {
		return(0);
	}
	
	/* ad hoc inclination factor */
	f = 1.0 - 0.3 * inc / FB_CONST_PI + fquad;
	
	/* otherwise use usual Mardling stability criterion */
	if (Rpout >= C*f*pow((1.0+qout)*(1.0+eout)/sqrt(1.0-eout), 0.4)*ain) {
		fb_dprintf("fewbody: mardling(): stable triple or quadruple: Rpout/Rpout,crit=%.6g eout=%.6g inc=%.6g degrees fquad=%.6g\n", 
			   Rpout/(C*f*pow((1.0+qout)*(1.0+eout)/sqrt(1.0-eout), 0.4)*ain),
			   eout, inc * 180.0/FB_CONST_PI, fquad);
		return(1);
	} else {
		fb_dprintf("fewbody: mardling(): unstable triple or quadruple: Rpout/Rpout,crit=%.6g eout=%.6g inc=%.6g degrees fquad=%.6g\n", 
			   Rpout/(C*f*pow((1.0+qout)*(1.0+eout)/sqrt(1.0-eout), 0.4)*ain),
			   eout, inc * 180.0/FB_CONST_PI, fquad);
		return(0);
	}
}
