/* -*- linux-c -*- */
/* fewbody_collapse.c

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

/* build the binary tree, subject to tidal criterion */
int fb_collapse(fb_hier_t *hier, double t, double tidaltol, double speedtol, fb_units_t units, fb_nonks_params_t nonks_params, fb_input_t input)
{
	int i, j, k, n, isave[2], cont=1, retval=0;
	double a, amin, E, xrel[3], v0[3], v1[3], vcm[3], ftid, E_pN;

	/* we only need the even terms here, since those are the conservative energy terms */
	int PN_on = nonks_params.PN1 + nonks_params.PN2 + nonks_params.PN3;

	/* first find the tightest binary and test to see whether it is unperturbed */
	while (cont) {
		amin = FB_AMIN;
		cont = 0;
		for (j=0; j<hier->nobj; j++) {
			for (k=j+1; k<hier->nobj; k++) {
				for (i=0; i<3; i++) {
					xrel[i] = hier->obj[j]->x[i] - hier->obj[k]->x[i];
					vcm[i] = (hier->obj[j]->m * hier->obj[j]->v[i] + \
						  hier->obj[k]->m * hier->obj[k]->v[i]) /\
						(hier->obj[j]->m + hier->obj[k]->m);
					v0[i] = hier->obj[j]->v[i] - vcm[i];
					v1[i] = hier->obj[k]->v[i] - vcm[i];
				}
				
				E = 0.5 * (hier->obj[j]->m * fb_dot(v0, v0) + hier->obj[k]->m * fb_dot(v1, v1)) - \
					hier->obj[j]->m * hier->obj[k]->m / fb_mod(xrel);

				E_pN = fb_E_rel(hier->obj[j], hier->obj[k], units, nonks_params);

				if (E < 0.0 && (!PN_on || E_pN < 0.0)) {
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
			/* swap indices so object on the left contains more stars */
			if (hier->obj[isave[0]]->n < hier->obj[isave[1]]->n) {
				j = isave[0];
				isave[0] = isave[1];
				isave[1] = j;
			}
			
			/* temporarily store new hierarchical object */
			n = hier->obj[isave[0]]->n + hier->obj[isave[1]]->n;
			hier->hier[hier->hi[n]+hier->narr[n]].obj[0] = hier->obj[isave[0]];
			hier->hier[hier->hi[n]+hier->narr[n]].obj[1] = hier->obj[isave[1]];
			fb_upsync(&(hier->hier[hier->hi[n]+hier->narr[n]]), t);
			
			/* test for tidal perturbation */
			ftid = 0.0;
			for (j=0; j<hier->nobj; j++) {
				if (j != isave[0] && j != isave[1]) {
					for (k=0; k<3; k++) {
						xrel[k] = hier->obj[j]->x[k] - hier->hier[hier->hi[n]+hier->narr[n]].x[k];
					}
					ftid += fb_reltide(&(hier->hier[hier->hi[n]+hier->narr[n]]), hier->obj[j], fb_mod(xrel));
				}
			}
			
			/* can't collapse if tidal force is too large */
			if (ftid >= tidaltol) {
				cont = 0;
			}
				
			/* can't collapse if the object is not stable */
			if (!fb_is_stable(&(hier->hier[hier->hi[n]+hier->narr[n]]), speedtol, units, input)) {
				cont = 0;
			}

			/* also can't collapse if internal tide is ever too large, since then we can't do
			   the orbits analytically */
			if (n > 2) {
				if (hier->hier[hier->hi[n]+hier->narr[n]].obj[0]->n >= 2) {
					if (fb_reltide(hier->hier[hier->hi[n]+hier->narr[n]].obj[0], hier->hier[hier->hi[n]+hier->narr[n]].obj[1], hier->hier[hier->hi[n]+hier->narr[n]].a*(1.0-hier->hier[hier->hi[n]+hier->narr[n]].e)) >= tidaltol) {
						fb_dprintf("fewbody: collapse(): not collapsing n=%d hierarchy due to tide\n", n);
						cont = 0;
					}
				}
				if (hier->hier[hier->hi[n]+hier->narr[n]].obj[1]->n >= 2) {
					if (fb_reltide(hier->hier[hier->hi[n]+hier->narr[n]].obj[1], hier->hier[hier->hi[n]+hier->narr[n]].obj[0], hier->hier[hier->hi[n]+hier->narr[n]].a*(1.0-hier->hier[hier->hi[n]+hier->narr[n]].e)) >= tidaltol) {
						fb_dprintf("fewbody: collapse(): not collapsing n=%d hierarchy due to tide\n", n);
						cont = 0;
					}
				}
			}

			/* can we actually collapse? */
			if (cont) {
				retval = 1;
				
				/* create new hierarchical object */
				hier->obj[isave[0]] = &(hier->hier[hier->hi[n]+hier->narr[n]]);
				hier->narr[n]++;
				
				/* delete extraneous object */
				hier->obj[isave[1]] = hier->obj[hier->nobj-1];
				hier->nobj--;
			}
		}
	}

	return(retval);
}

/* expand the tree if a tide is exceeded */
int fb_expand(fb_hier_t *hier, double t, double tidaltol)
{
	int i, j, k, cont=1, retval=0, n;
	double xrel[3], ftid;
	fb_obj_t *obj1ptr, *obj2ptr;
	
	while (cont) {
		cont = 0;
		for (i=0; i<hier->nobj; i++) {
			if (hier->obj[i]->n > 1) {
				ftid = 0.0;
				for (j=0; j<hier->nobj; j++) {
					if (j != i) {
						for (k=0; k<3; k++) {
							xrel[k] = hier->obj[i]->x[k] - hier->obj[j]->x[k];
						}
						ftid += fb_reltide(hier->obj[i], hier->obj[j], fb_mod(xrel));
					}
				}
				if (ftid >= tidaltol) {
					cont = 1;
					break;
				}
			}
		}

		/* expand stuff */
		if (cont) {
			retval = 1;
			n = hier->obj[i]->n;

			/* store the pointers of the two objects to swap */
			obj1ptr = hier->obj[i];
			obj2ptr = &(hier->hier[hier->hi[n]+hier->narr[n]-1]);

			/* make new obj */
			fb_downsync(hier->obj[i], t);			
			hier->obj[hier->nobj] = hier->obj[i]->obj[1];
			hier->obj[i] = hier->obj[i]->obj[0];
			hier->nobj++;
			
			/* musical pointers */
			for (i=0; i<hier->nobj; i++) {
				if (hier->obj[i] == obj2ptr) {
					hier->obj[i] = obj1ptr;
				}
			}
			if (n < hier->nstar) {
				for (j=n+1; j<=hier->nstar; j++) {
					for (i=0; i<hier->narr[j]; i++) {
						if (hier->hier[hier->hi[j]+i].obj[0] == obj2ptr) {
							hier->hier[hier->hi[j]+i].obj[0] = obj1ptr;
						}
						if (hier->hier[hier->hi[j]+i].obj[1] == obj2ptr) {
							hier->hier[hier->hi[j]+i].obj[1] = obj1ptr;
						}
					}
				}
			}
			fb_objcpy(obj1ptr, obj2ptr);
			hier->narr[n]--;
		}
	}
	
	return(retval);
}
