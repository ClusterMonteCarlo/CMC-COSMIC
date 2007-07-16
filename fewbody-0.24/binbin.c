/* -*- linux-c -*- */
/* binbin.c

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
#include <time.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "fewbody.h"
#include "binbin.h"

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  binbin [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -m --m00 <m00/MSUN>          : set mass of star 0 of binary 0 [%.6g]\n", FB_M00/FB_CONST_MSUN);
	fprintf(stream, "  -n --m01 <m01/MSUN>          : set mass of star 1 of binary 0 [%.6g]\n", FB_M01/FB_CONST_MSUN);
	fprintf(stream, "  -o --m10 <m10/MSUN>          : set mass of star 0 of binary 1 [%.6g]\n", FB_M10/FB_CONST_MSUN);
	fprintf(stream, "  -p --m11 <m11/MSUN>          : set mass of star 1 of binary 1 [%.6g]\n", FB_M11/FB_CONST_MSUN);
	fprintf(stream, "  -r --r00 <r00/RSUN>          : set radius of star 0 of binary 0 [%.6g]\n", FB_R00/FB_CONST_RSUN);
	fprintf(stream, "  -g --r01 <r01/RSUN>          : set radius of star 1 of binary 0 [%.6g]\n", FB_R01/FB_CONST_RSUN);
	fprintf(stream, "  -i --r10 <r10/RSUN>          : set radius of star 0 of binary 1 [%.6g]\n", FB_R10/FB_CONST_RSUN);
	fprintf(stream, "  -j --r11 <r11/RSUN>          : set radius of star 1 of binary 1 [%.6g]\n", FB_R11/FB_CONST_RSUN);
	fprintf(stream, "  -a --a0 <a0/AU>              : set semimajor axis of binary 0 [%.6g]\n", FB_A0/FB_CONST_AU);
	fprintf(stream, "  -q --a1 <a1/AU>              : set semimajor axis of binary 1 [%.6g]\n", FB_A1/FB_CONST_AU);
	fprintf(stream, "  -e --e0 <e0>                 : set eccentricity of binary 0 [%.6g]\n", FB_E0);
	fprintf(stream, "  -f --e1 <e1>                 : set eccentricity of binary 1 [%.6g]\n", FB_E1);
	fprintf(stream, "  -v --vinf <vinf/v_crit>      : set velocity at infinity [%.6g]\n", FB_VINF);
	fprintf(stream, "  -b --b <b/(a0+a1)>           : set impact parameter [%.6g]\n", FB_B);
	fprintf(stream, "  -t --tstop <tstop/t_dyn>     : set stopping time [%.6g]\n", FB_TSTOP);
	fprintf(stream, "  -D --dt <dt/t_dyn>           : set approximate output dt [%.6g]\n", FB_DT);
	fprintf(stream, "  -c --tcpustop <tcpustop/sec> : set cpu stopping time [%.6g]\n", FB_TCPUSTOP);
	fprintf(stream, "  -A --absacc <absacc>         : set integrator's absolute accuracy [%.6g]\n", FB_ABSACC);
	fprintf(stream, "  -R --relacc <relacc>         : set integrator's relative accuracy [%.6g]\n", FB_RELACC);
	fprintf(stream, "  -N --ncount <ncount>         : set number of integration steps between calls\n");
	fprintf(stream, "                                 to fb_classify() [%d]\n", FB_NCOUNT);
	fprintf(stream, "  -z --tidaltol <tidaltol>     : set tidal tolerance [%.6g]\n", FB_TIDALTOL);
	fprintf(stream, "  -x --fexp <f_exp>            : set expansion factor of merger product [%.6g]\n", FB_FEXP);
	fprintf(stream, "  -k --ks                      : turn K-S regularization on or off [%d]\n", FB_KS);
	fprintf(stream, "  -s --seed                    : set random seed [%ld]\n", FB_SEED);
	fprintf(stream, "  -d --debug                   : turn on debugging\n");
	fprintf(stream, "  -V --version                 : print version info\n");
	fprintf(stream, "  -h --help                    : display this help text\n");
}

/* calculate the units used */
int calc_units(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[0]->obj[0]->m * obj[0]->obj[1]->m / obj[0]->a + \
			 obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	units->l = obj[0]->a + obj[1]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
	
	return(0);
}

/* the main attraction */
int main(int argc, char *argv[])
{
	int i, j;
	unsigned long int seed;
	double m00, m01, m10, m11, r00, r01, r10, r11, a0, a1, e0, e1;
	double rtid, vinf, b, m0, m1, M, mu, Ei, E, Lint[3], Li[3], l0[3], l1[3], L[3], r[3], t;
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;
	const char *short_opts = "m:n:o:p:r:g:i:j:a:q:e:f:v:b:t:D:c:A:R:N:z:x:k:s:dVh";
	const struct option long_opts[] = {
		{"m00", required_argument, NULL, 'm'},
		{"m01", required_argument, NULL, 'n'},
		{"m10", required_argument, NULL, 'o'},
		{"m11", required_argument, NULL, 'p'},
		{"r00", required_argument, NULL, 'r'},
		{"r01", required_argument, NULL, 'g'},
		{"r10", required_argument, NULL, 'i'},
		{"r11", required_argument, NULL, 'j'},
		{"a0", required_argument, NULL, 'a'},
		{"a1", required_argument, NULL, 'q'},
		{"e0", required_argument, NULL, 'e'},
		{"e1", required_argument, NULL, 'f'},
		{"vinf", required_argument, NULL, 'v'},
		{"b", required_argument, NULL, 'b'},
		{"tstop", required_argument, NULL, 't'},
		{"dt", required_argument, NULL, 'D'},
		{"tcpustop", required_argument, NULL, 'c'},
		{"absacc", required_argument, NULL, 'A'},
		{"relacc", required_argument, NULL, 'R'},
		{"ncount", required_argument, NULL, 'N'},
		{"tidaltol", required_argument, NULL, 'z'},
		{"fexp", required_argument, NULL, 'x'},
		{"ks", required_argument, NULL, 'k'},
		{"seed", required_argument, NULL, 's'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};

	/* set parameters to default values */
	m00 = FB_M00;
	m01 = FB_M01;
	m10 = FB_M10;
	m11 = FB_M11;
	r00 = FB_R00;
	r01 = FB_R01;
	r10 = FB_R10;
	r11 = FB_R11;
	a0 = FB_A0;
	a1 = FB_A1;
	e0 = FB_E0;
	e1 = FB_E1;
	vinf = FB_VINF;
	b = FB_B;
	input.ks = FB_KS;
	input.tstop = FB_TSTOP;
	input.Dflag = 0;
	input.dt = FB_DT;
	input.tcpustop = FB_TCPUSTOP;
	input.absacc = FB_ABSACC;
	input.relacc = FB_RELACC;
	input.ncount = FB_NCOUNT;
	input.tidaltol = FB_TIDALTOL;
	input.fexp = FB_FEXP;
	seed = FB_SEED;
	fb_debug = FB_DEBUG;
	
	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'm':
			m00 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'n':
			m01 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'o':
			m10 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'p':
			m11 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'r':
			r00 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'g':
			r01 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'i':
			r10 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'j':
			r11 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'a':
			a0 = atof(optarg) * FB_CONST_AU;
			break;
		case 'q':
			a1 = atof(optarg) * FB_CONST_AU;
			break;
		case 'e':
			e0 = atof(optarg);
			if (e0 >= 1.0) {
				fprintf(stderr, "e0 must be less than 1\n");
				return(1);
			}
			break;
		case 'f':
			e1 = atof(optarg);
			if (e1 >= 1.0) {
				fprintf(stderr, "e1 must be less than 1\n");
				return(1);
			}
			break;
		case 'v':
			vinf = atof(optarg);
			if (vinf < 0.0) {
				fprintf(stderr, "vinf must be non-negative\n");
				return(1);
			}
			break;
		case 'b':
			b = atof(optarg);
			if (b < 0.0) {
				fprintf(stderr, "b must be non-negative\n");
				return(1);
			}
			break;
		case 't':
			input.tstop = atof(optarg);
			break;
		case 'D':
			input.Dflag = 1;
			input.dt = atof(optarg);
			break;
		case 'c':
			input.tcpustop = atof(optarg);
			break;
		case 'A':
			input.absacc = atof(optarg);
			break;
		case 'R':
			input.relacc = atof(optarg);
			break;
		case 'N':
			input.ncount = atoi(optarg);
			break;
		case 'z':
			input.tidaltol = atof(optarg);
			break;
		case 'x':
			input.fexp = atof(optarg);
			break;
		case 'k':
			input.ks = atoi(optarg);
			break;
		case 's':
			seed = atol(optarg);
			break;
		case 'd':
			fb_debug = 1;
			break;
		case 'V':
			fb_print_version(stdout);
			return(0);
		case 'h':
			fb_print_version(stdout);
			fprintf(stdout, "\n");
			print_usage(stdout);
			return(0);
		default:
			break;
		}
	}
	
	/* check to make sure there was nothing crazy on the command line */
	if (optind < argc) {
		print_usage(stdout);
		return(1);
	}

	/* initialize a few things for integrator */
	t = 0.0;
	hier.nstarinit = 4;
	hier.nstar = 4;
	fb_malloc_hier(&hier);
	fb_init_hier(&hier);

	/* put stuff in log entry */
	snprintf(input.firstlogentry, FB_MAX_LOGENTRY_LENGTH, "  command line:");
	for (i=0; i<argc; i++) {
		snprintf(&(input.firstlogentry[strlen(input.firstlogentry)]), 
			 FB_MAX_LOGENTRY_LENGTH-strlen(input.firstlogentry), " %s", argv[i]);
	}
	snprintf(&(input.firstlogentry[strlen(input.firstlogentry)]), 
		 FB_MAX_LOGENTRY_LENGTH-strlen(input.firstlogentry), "\n");
	
	/* print out values of paramaters */
	fprintf(stderr, "PARAMETERS:\n");
	fprintf(stderr, "  ks=%d  seed=%ld\n", input.ks, seed);
	fprintf(stderr, "  a0=%.6g AU  e0=%.6g  m00=%.6g MSUN  m01=%.6g MSUN  r00=%.6g RSUN  r01=%.6g RSUN\n", \
		a0/FB_CONST_AU, e0, m00/FB_CONST_MSUN, m01/FB_CONST_MSUN, r00/FB_CONST_RSUN, r01/FB_CONST_RSUN);
	fprintf(stderr, "  a1=%.6g AU  e1=%.6g  m10=%.6g MSUN  m11=%.6g MSUN  r10=%.6g RSUN  r11=%.6g RSUN\n", \
		a1/FB_CONST_AU, e1, m10/FB_CONST_MSUN, m11/FB_CONST_MSUN, r10/FB_CONST_RSUN, r11/FB_CONST_RSUN);
	fprintf(stderr, "  vinf=%.6g  b=%.6g  tstop=%.6g  tcpustop=%.6g\n", \
		vinf, b, input.tstop, input.tcpustop);
	fprintf(stderr, "  tidaltol=%.6g  abs_acc=%.6g  rel_acc=%.6g  ncount=%d  fexp=%.6g\n\n", \
		input.tidaltol, input.absacc, input.relacc, input.ncount, input.fexp);

	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, seed);

	/* create binaries */
	hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
	hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
	hier.hier[hier.hi[2]+0].t = t;
	hier.hier[hier.hi[2]+1].obj[0] = &(hier.hier[hier.hi[1]+2]);
	hier.hier[hier.hi[2]+1].obj[1] = &(hier.hier[hier.hi[1]+3]);
	hier.hier[hier.hi[2]+1].t = t;

	/* give the objects some properties */
	for (j=0; j<hier.nstar; j++) {
		hier.hier[hier.hi[1]+j].ncoll = 1;
		hier.hier[hier.hi[1]+j].id[0] = j;
		snprintf(hier.hier[hier.hi[1]+j].idstring, FB_MAX_STRING_LENGTH, "%d", j);
		hier.hier[hier.hi[1]+j].n = 1;
		hier.hier[hier.hi[1]+j].obj[0] = NULL;
		hier.hier[hier.hi[1]+j].obj[1] = NULL;
		hier.hier[hier.hi[1]+j].Eint = 0.0;
		hier.hier[hier.hi[1]+j].Lint[0] = 0.0;
		hier.hier[hier.hi[1]+j].Lint[1] = 0.0;
		hier.hier[hier.hi[1]+j].Lint[2] = 0.0;
	}

	hier.hier[hier.hi[1]+0].R = r00;
	hier.hier[hier.hi[1]+1].R = r01;
	hier.hier[hier.hi[1]+2].R = r10;
	hier.hier[hier.hi[1]+3].R = r11;

	hier.hier[hier.hi[1]+0].m = m00;
	hier.hier[hier.hi[1]+1].m = m01;
	hier.hier[hier.hi[1]+2].m = m10;
	hier.hier[hier.hi[1]+3].m = m11;

	hier.hier[hier.hi[2]+0].m = m00 + m01;
	hier.hier[hier.hi[2]+1].m = m10 + m11;

	hier.hier[hier.hi[2]+0].a = a0;
	hier.hier[hier.hi[2]+1].a = a1;
	
	hier.hier[hier.hi[2]+0].e = e0;
	hier.hier[hier.hi[2]+1].e = e1;

	hier.obj[0] = &(hier.hier[hier.hi[2]+0]);
	hier.obj[1] = &(hier.hier[hier.hi[2]+1]);
	hier.obj[2] = NULL;
	hier.obj[3] = NULL;

	/* get the units and normalize */
	calc_units(hier.obj, &units);
	fb_normalize(&hier, units);
	
	fprintf(stderr, "UNITS:\n");
	fprintf(stderr, "  v=v_crit=%.6g km/s  l=%.6g AU  t=t_dyn=%.6g yr\n", \
		units.v/1.0e5, units.l/FB_CONST_AU, units.t/FB_CONST_YR);
	fprintf(stderr, "  M=%.6g M_sun  E=%.6g erg\n\n", units.m/FB_CONST_MSUN, units.E);
	
	/* move hierarchies analytically in from infinity along hyperbolic orbit */
	fb_dprintf("moving hierarchies analytically in from infinity...\n");

	m0 = hier.obj[0]->m;
	m1 = hier.obj[1]->m;
	M = m0 + m1;
	mu = m0 * m1 / M;

	Ei = 0.5 * mu * fb_sqr(vinf);

	a0 = hier.obj[0]->a;
	a1 = hier.obj[1]->a;

	e0 = hier.obj[0]->e;
	e1 = hier.obj[1]->e;
	
	m00 = hier.obj[0]->obj[0]->m;
	m01 = hier.obj[0]->obj[1]->m;
	m10 = hier.obj[1]->obj[0]->m;
	m11 = hier.obj[1]->obj[1]->m;

	rtid = pow(2.0*(m0+m1)/input.tidaltol, 1.0/3.0) * \
		FB_MAX(pow(m0, -1.0/3.0)*a0*(1.0+e0), pow(m1, -1.0/3.0)*a1*(1.0+e1));

	fb_init_scattering(hier.obj, vinf, b, rtid);
	
	/* and check to see that we conserved energy and angular momentum */
	fb_cross(hier.obj[0]->x, hier.obj[0]->v, l0);
	fb_cross(hier.obj[1]->x, hier.obj[1]->v, l1);
	
	for (i=0; i<3; i++) {
		L[i] = (m0 * l0[i] + m1 * l1[i]);
		r[i] = hier.obj[1]->x[i] - hier.obj[0]->x[i];
	}

	E = - m0 * m1 / fb_mod(r) + 0.5 * (m0 * fb_dot(hier.obj[0]->v, hier.obj[0]->v) + \
					m1 * fb_dot(hier.obj[1]->v, hier.obj[1]->v));

	fb_dprintf("L0=%.6g DeltaL/L0=%.6g DeltaL=%.6g\n", mu*b*vinf, fb_mod(L)/(mu*b*vinf)-1.0, fb_mod(L)-mu*b*vinf);
	fb_dprintf("E0=%.6g DeltaE/E0=%.6g DeltaE=%.6g\n\n", Ei, E/Ei-1.0, E-Ei);

	/* trickle down the binary properties, then back up */
	for (j=0; j<2; j++) {
		fb_dprintf("obj[%d]->a=%e\n", j, hier.obj[j]->a);
		fb_dprintf("obj[%d]->e=%e\n", j, hier.obj[j]->e);
		fb_dprintf("obj[%d]->m=%e\n", j, hier.obj[j]->m);
		fb_randorient(&(hier.hier[hier.hi[2]+j]), rng);
		fb_downsync(&(hier.hier[hier.hi[2]+j]), t);
		fb_upsync(&(hier.hier[hier.hi[2]+j]), t);
		fb_dprintf("obj[%d]->a=%e\n", j, hier.obj[j]->a);
		fb_dprintf("obj[%d]->e=%e\n", j, hier.obj[j]->e);
		fb_dprintf("obj[%d]->m=%e\n", j, hier.obj[j]->m);
	}
	
	/* store the initial energy and angular momentum*/
	Ei = fb_petot(&(hier.hier[hier.hi[1]]), hier.nstar) + fb_ketot(&(hier.hier[hier.hi[1]]), hier.nstar) +
		fb_einttot(&(hier.hier[hier.hi[1]]), hier.nstar);
	fb_angmom(&(hier.hier[hier.hi[1]]), hier.nstar, Li);
	fb_angmomint(&(hier.hier[hier.hi[1]]), hier.nstar, Lint);
	for (j=0; j<3; j++) {
		Li[j] += Lint[j];
	}

	/* integrate along */
	fb_dprintf("calling fewbody()...\n");
	
	/* call fewbody! */
	retval = fewbody(input, &hier, &t);

	/* print information to screen */
	fprintf(stderr, "OUTCOME:\n");
	if (retval.retval == 1) {
		fprintf(stderr, "  encounter complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
			t, t * units.t/FB_CONST_YR,
			fb_sprint_hier(hier, string1),
			fb_sprint_hier_hr(hier, string2));
	} else {
		fprintf(stderr, "  encounter NOT complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
			t, t * units.t/FB_CONST_YR,
			fb_sprint_hier(hier, string1),
			fb_sprint_hier_hr(hier, string2));
	}

	fb_dprintf("there were %ld integration steps\n", retval.count);
	fb_dprintf("fb_classify() was called %ld times\n", retval.iclassify);
	
	fprintf(stderr, "FINAL:\n");
	fprintf(stderr, "  t_final=%.6g (%.6g yr)  t_cpu=%.6g s\n", \
		t, t*units.t/FB_CONST_YR, retval.tcpu);

	fprintf(stderr, "  L0=%.6g  DeltaL/L0=%.6g  DeltaL=%.6g\n", fb_mod(Li), retval.DeltaLfrac, retval.DeltaL);
	fprintf(stderr, "  E0=%.6g  DeltaE/E0=%.6g  DeltaE=%.6g\n", Ei, retval.DeltaEfrac, retval.DeltaE);
	fprintf(stderr, "  Rmin=%.6g (%.6g RSUN)  Rmin_i=%d  Rmin_j=%d\n", \
		retval.Rmin, retval.Rmin*units.l/FB_CONST_RSUN, retval.Rmin_i, retval.Rmin_j);
	fprintf(stderr, "  Nosc=%d (%s)\n", retval.Nosc, (retval.Nosc>=1?"resonance":"non-resonance"));
	
	/* free GSL stuff */
	gsl_rng_free(rng);

	/* free our own stuff */
	fb_free_hier(hier);

	/* done! */
	return(0);
}
