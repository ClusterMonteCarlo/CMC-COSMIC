/* -*- linux-c -*- */
/* fewbody_io.c

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

/* print the version */
void fb_print_version(FILE *stream)
{
	fprintf(stream, "** Fewbody %s (%s) [%s] **\n", FB_VERSION, FB_NICK, FB_DATE);
}

/* print the output in Starlab story format */
void fb_print_story(fb_obj_t *star, int nstar, double t, char *logentry)
{
	int i, j;
	double mtot, r[3], v[3], E, L[3], Lint[3];
	
	/* calculate total mass and the motion of the center of mass (should be zero) */
	mtot = 0.0;
	for (j=0; j<3; j++) {
		r[j] = 0.0;
		v[j] = 0.0;
	}
	for (i=0; i<nstar; i++) {
		mtot += star[i].m;
		for (j=0; j<3; j++) {
			r[j] += star[i].m * star[i].x[j];
			v[j] += star[i].m * star[i].v[j];
		}
	}
	for (j=0; j<3; j++) {
		r[j] /= mtot;
		v[j] /= mtot;
	}

	E = fb_petot(star, nstar) + fb_ketot(star, nstar) + fb_einttot(star, nstar);
	fb_angmom(star, nstar, L);
	fb_angmomint(star, nstar, Lint);
	for (j=0; j<3; j++) {
		L[j] += Lint[j];
	}
	
	fprintf(stdout, "(Particle\n");
	fprintf(stdout, "  N  =  %d\n", nstar);
	
	fprintf(stdout, "(Log\n");
	fprintf(stdout, "%s", logentry);
	logentry[0] = '\0';
	fprintf(stdout, ")Log\n");
	
	fprintf(stdout, "(Dynamics\n");
	fprintf(stdout, "  system_time  =  %.9g\n", t);
	fprintf(stdout, "  t  =  %.9g\n", t);
	fprintf(stdout, "  m  =  %.9g\n", mtot);
	fprintf(stdout, "  r  =  %.9g  %.9g  %.9g\n", r[0], r[1], r[2]);
	fprintf(stdout, "  v  =  %.9g  %.9g  %.9g\n", v[0], v[1], v[2]);
	/* what to do here? */
	fprintf(stdout, "  R_eff  =  %.9g\n", 0.0);
	fprintf(stdout, "  E  =  %.9g\n", E);
	fprintf(stdout, "  L  =  %.9g  %.9g  %.9g\n", L[0], L[1], L[2]);
	fprintf(stdout, ")Dynamics\n");
	
	fprintf(stdout, "(Hydro\n");
	fprintf(stdout, ")Hydro\n");
	
	fprintf(stdout, "(Star\n");
	fprintf(stdout, ")Star\n");
	
	for (i=0; i<nstar; i++) {
		fprintf(stdout, "(Particle\n");
		fprintf(stdout, "  i  =  %d\n", i+1);
		fprintf(stdout, "  N  =  %d\n", 1);
		fprintf(stdout, "(Dynamics\n");
		fprintf(stdout, "  t  =  %.9g\n", t);
		fprintf(stdout, "  m  =  %.9g\n", star[i].m);
		fprintf(stdout, "  r  =  %.9g  %.9g  %.9g\n", star[i].x[0], star[i].x[1], star[i].x[2]);
		fprintf(stdout, "  v  =  %.9g  %.9g  %.9g\n", star[i].v[0], star[i].v[1], star[i].v[2]);
		fprintf(stdout, "  R_eff  =  %.9g\n", star[i].R);
		fprintf(stdout, ")Dynamics\n");
		fprintf(stdout, ")Particle\n");
	}
	
	fprintf(stdout, ")Particle\n");
}
