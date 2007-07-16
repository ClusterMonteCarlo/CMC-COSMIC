/* -*- linux-c -*- */
/* fewbody_coll.c

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

/* the main collision criterion */
int fb_is_collision(double r, double R1, double R2)
{
	if (r < R1 + R2) {
		return(1);
	} else {
		return(0);
	}
}

int fb_collide(fb_hier_t *hier, double f_exp)
{
	int i, j=-1, k, retval=0, cont=1;
	double R[3], peinit;

	/* this is a non-recursive way to perform a recursive operation: keep going until there are no more
	   mergers */
	while (cont) {
		cont = 0;
		for (i=0; i<hier->nstar-1; i++) {
			for (j=i+1; j<hier->nstar; j++) {
				/* calculate relative separation */
				for (k=0; k<3; k++) {
					R[k] = hier->hier[hier->hi[1]+i].x[k] - hier->hier[hier->hi[1]+j].x[k];
				}

				/* test collision criterion */
				if (fb_is_collision(fb_mod(R), hier->hier[hier->hi[1]+i].R, hier->hier[hier->hi[1]+j].R)) {
					cont = 1;
					/* break out of the double loop if there is a collision, so we can merge
					   the stars immediately */
					break;
				}
			}
			/* break out of the double loop if there is a collision, so we can merge
			   the stars immediately */
			if (cont) {
				break;
			}
		}
		
		/* merge two stars if necessary */
		if (cont) {
			/* return 1 if there is a collision */
			retval = 1;

			/* calculate the potential energy before the collision, since we're going to need
			   to account for the change in potential energy, and put it in Eint */
			peinit = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);			
			
			/* do the actual merger */
			fb_dprintf("fewbody: collide(): merging stars: i=%d j=%d\n", i, j);
			fb_merge(&(hier->hier[hier->hi[1]+i]), &(hier->hier[hier->hi[1]+j]), hier->nstarinit, f_exp);
			fb_objcpy(&(hier->hier[hier->hi[1]+j]), &(hier->hier[hier->hi[1]+hier->nstar-1]));
			hier->nstar--;

			/* calculate the difference in potential energy before and after the collision, and put it
			   in Eint for accounting */
			hier->hier[hier->hi[1]+i].Eint += peinit - fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);
		}
	}

	return(retval);
}

void fb_merge(fb_obj_t *obj1, fb_obj_t *obj2, int nstarinit, double f_exp)
{
	int i;
	double x1[3], x2[3], v1[3], v2[3], l1[3], l2[3];
	fb_obj_t tmpobj;

	/* sanity check */
	if (obj1->n != 1 || obj2->n != 1) {
		fprintf(stderr, "fb_merge: trying to merge an object that isn't single?!\n");
		exit(1);
	}

	/* need some temporary storage here */
	tmpobj.id = (long *) malloc(nstarinit * sizeof(long));

	/* merge id's */
	tmpobj.ncoll = obj1->ncoll + obj2->ncoll;
	for (i=0; i<obj1->ncoll; i++) {
		tmpobj.id[i] = obj1->id[i];
	}
	for (i=0; i<obj2->ncoll; i++) {
		tmpobj.id[obj1->ncoll + i] = obj2->id[i];
	}

	/* create idstring */
	snprintf(tmpobj.idstring, FB_MAX_STRING_LENGTH, "%s:%s", obj1->idstring, obj2->idstring);

	/* assume no mass loss */
	tmpobj.m = obj1->m + obj2->m;

	/* this is just a simple prescription */
	tmpobj.R = f_exp * (obj1->R + obj2->R);

	/* set new position and velocity, calculate relative positions and velocities */
	for (i=0; i<3; i++) {
		tmpobj.x[i] = (obj1->m * obj1->x[i] + obj2->m * obj2->x[i]) / tmpobj.m;
		tmpobj.v[i] = (obj1->m * obj1->v[i] + obj2->m * obj2->v[i]) / tmpobj.m;
		x1[i] = obj1->x[i] - tmpobj.x[i];
		x2[i] = obj2->x[i] - tmpobj.x[i];
		v1[i] = obj1->v[i] - tmpobj.v[i];
		v2[i] = obj2->v[i] - tmpobj.v[i];
	}

	/* set internal energy, using the difference in kinetic energy; the difference in potential energy
	   depends on the positions of the other stars, and will be calculated later and added to Eint */
	tmpobj.Eint = obj1->Eint + obj2->Eint +
		0.5 * (obj1->m * fb_dot(obj1->v, obj1->v) + obj2->m * fb_dot(obj2->v, obj2->v)) -
		0.5 * tmpobj.m * fb_dot(tmpobj.v, tmpobj.v);

	/* set internal angular momentum */
	fb_cross(x1, v1, l1);
	fb_cross(x2, v2, l2);
	for (i=0; i<3; i++) {
		tmpobj.Lint[i] = obj1->Lint[i] + obj2->Lint[i] + obj1->m * l1[i] + obj2->m * l2[i];
	}
	
	/* and better set these, too... */
	tmpobj.n = 1;
	
	tmpobj.obj[0] = NULL;
	tmpobj.obj[1] = NULL;

	/* finally, copy over the merger from temporary storage */
	fb_objcpy(obj1, &tmpobj);

	free(tmpobj.id);
}
