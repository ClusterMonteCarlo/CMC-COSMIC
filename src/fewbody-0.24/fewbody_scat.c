/* -*- linux-c -*- */
/* fewbody_scat.c

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

/* move objects in from infinity analytically along a hyperbolic orbit */
void fb_init_scattering(fb_obj_t *obj[2], double vinf, double b, double rtid)
{
	double m0, m1, M, r, rperi, x, y, vx, vy, qa, qb, qc, rad;
	
	/* some useful variables */
	m0 = obj[0]->m;
	m1 = obj[1]->m;
	M = m0 + m1;

	/* pericenter */
	rperi = (vinf==0.0?0.0:(M/fb_sqr(vinf)*(sqrt(1.0+fb_sqr(b*fb_sqr(vinf)/M))-1.0)));

	/* make sure r>=rperi, otherwise analytically moving the obj's below will give NANs */
	r = FB_MAX(rtid, rperi);

	/* solve for y */
	qa = M*M + fb_sqr(b*vinf*vinf);
	qb = 2.0*(M*r-fb_sqr(b*vinf))*b*vinf*vinf;
	qc = -2.0*fb_sqr(b*vinf)*M*r + fb_sqr(fb_sqr(b*vinf));
	rad = qb*qb-4.0*qa*qc;
	y = (-qb+(rad<=0.0?0.0:sqrt(rad)))/(2.0*qa);

	/* x can sometimes be zero, so need to worry about the square root here */
	x = (r*r-y*y<=0.0?0.0:sqrt(r*r-y*y));

	/* determine v_x from quadratic */
	if (x > 0.0) {
		qa = 1.0 + fb_sqr(y/x);
		qb = 2.0*y*b*vinf/(x*x);
		qc = fb_sqr(b*vinf/x) - 2.0*M/r - vinf*vinf;
		rad = qb*qb-4.0*qa*qc;
		vx = (-qb-(rad<=0.0?0.0:sqrt(rad)))/(2.0*qa);
		vy = (b*vinf+y*vx)/x;
	} else { /* x=0 */
		vx = -vinf;
		vy = 0.0;
	}

	/* set the positions and velocities */
	obj[0]->x[0] = - x / (1.0 + m0/m1);
	obj[0]->x[1] = - y / (1.0 + m0/m1);
	obj[0]->x[2] = 0.0;
	
	obj[0]->v[0] = - vx / (1.0 + m0/m1);
	obj[0]->v[1] = - vy / (1.0 + m0/m1);
	obj[0]->v[2] = 0.0;

	obj[1]->x[0] = x / (1.0 + m1/m0);
	obj[1]->x[1] = y / (1.0 + m1/m0);
	obj[1]->x[2] = 0.0;
	
	obj[1]->v[0] = vx / (1.0 + m1/m0);
	obj[1]->v[1] = vy / (1.0 + m1/m0);
	obj[1]->v[2] = 0.0;
}

/* normalize all variables to our system of units */
void fb_normalize(fb_hier_t *hier, fb_units_t units)
{
	int i, k;

	/* just normalize everything, since it can't hurt */
	for (i=hier->hi[1]; i<=hier->hi[hier->nstarinit]; i++) {
		hier->hier[i].m /= units.m;
		hier->hier[i].R /= units.l;
		hier->hier[i].Eint /= units.E;
		hier->hier[i].a /= units.l;
		hier->hier[i].t /= units.t;
		for (k=0; k<3; k++) {
			hier->hier[i].x[k] /= units.l;
			hier->hier[i].v[k] /= units.v;
			hier->hier[i].Lint[k] /= units.m * units.l * units.v;
		}
	}
}
