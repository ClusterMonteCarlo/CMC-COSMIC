/* -*- linux-c -*- */
/* triplebin.c

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
#include "triplebin.h"

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  triplebin [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -m --m000 <m000/MSUN>        : set mass of star 0 of inner binary of triple [%.6g]\n", FB_M000/FB_CONST_MSUN);
	fprintf(stream, "  -n --m001 <m001/MSUN>        : set mass of star 1 of inner binary of triple [%.6g]\n", FB_M001/FB_CONST_MSUN);
	fprintf(stream, "  -o --m01 <m01/MSUN>          : set mass of outer star of triple [%.6g]\n", FB_M01/FB_CONST_MSUN);
	fprintf(stream, "  -p --m10 <m10/MSUN>          : set mass of star 0 of binary [%.6g]\n", FB_M10/FB_CONST_MSUN);
	fprintf(stream, "  -l --m11 <m11/MSUN>          : set mass of star 1 of binary [%.6g]\n", FB_M11/FB_CONST_MSUN);
	fprintf(stream, "  -r --r000 <r000/RSUN>        : set radius of star 0 of inner binary of triple [%.6g]\n", FB_R000/FB_CONST_RSUN);
	fprintf(stream, "  -g --r001 <r001/RSUN>        : set radius of star 1 of inner binary of triple [%.6g]\n", FB_R001/FB_CONST_RSUN);
	fprintf(stream, "  -i --r01 <r01/RSUN>          : set radius of outer star of triple [%.6g]\n", FB_R01/FB_CONST_RSUN);
	fprintf(stream, "  -j --r10 <r10/RSUN>          : set radius of star 0 of binary [%.6g]\n", FB_R10/FB_CONST_RSUN);
	fprintf(stream, "  -u --r11 <r11/RSUN>          : set radius of star 1 of binary [%.6g]\n", FB_R11/FB_CONST_RSUN);
	fprintf(stream, "  -a --a00 <a00/AU>            : set inner semimajor axis of triple [%.6g]\n", FB_A00/FB_CONST_AU);
	fprintf(stream, "  -Q --a0 <a0/AU>              : set outer semimajor axis of triple [%.6g]\n", FB_A0/FB_CONST_AU);
	fprintf(stream, "  -q --a1 <a1/AU>              : set semimajor axis of binary [%.6g]\n", FB_A1/FB_CONST_AU);
	fprintf(stream, "  -e --e00 <e00>               : set inner eccentricity of triple [%.6g]\n", FB_E00);
	fprintf(stream, "  -F --e0 <e0>                 : set outer eccentricity of triple [%.6g]\n", FB_E0);
	fprintf(stream, "  -f --e1 <e1>                 : set eccentricity of binary [%.6g]\n", FB_E1);
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
	double mu, m0, m1, m00, m01, m10, m11, m000, m001, a00, a0, a1;

	m0 = obj[0]->m;
	m1 = obj[1]->m;
	m00 = obj[0]->obj[0]->m;
	m01 = obj[0]->obj[1]->m;
	m10 = obj[1]->obj[0]->m;
	m11 = obj[1]->obj[1]->m;
	m000 = obj[0]->obj[0]->obj[0]->m;
	m001 = obj[0]->obj[0]->obj[1]->m;
	
	mu = m0*m1/(m0+m1);

	a0 = obj[0]->a;
	a1 = obj[1]->a;
	a00 = obj[0]->obj[0]->a;
	
	/* note that the unit of velocity here is set to an *approximate* value of v_c */
	units->v = sqrt(FB_CONST_G/mu*(m10*m11/a1 + m00*m01/a0 + m000*m001/a00));
	units->l = a0 + a1;
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
	double rtid, m000, m001, m00, m01, m10, m11, r000, r001, r01, r10, r11, a00, a0, a1, e00, e0, e1;
	double vcrit, vinf, b, m0, m1, M, mu, Ei, E, Lint[3], Li[3], l0[3], l1[3], L[3], r[3], t;
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;
	const char *short_opts = "m:n:o:p:l:r:g:i:j:u:a:Q:q:e:F:f:v:b:t:D:c:A:R:N:z:x:k:s:dVh";
	const struct option long_opts[] = {
		{"m000", required_argument, NULL, 'm'},
		{"m001", required_argument, NULL, 'n'},
		{"m01", required_argument, NULL, 'o'},
		{"m10", required_argument, NULL, 'p'},
		{"m11", required_argument, NULL, 'l'},
		{"r000", required_argument, NULL, 'r'},
		{"r001", required_argument, NULL, 'g'},
		{"r01", required_argument, NULL, 'i'},
		{"r10", required_argument, NULL, 'j'},
		{"r11", required_argument, NULL, 'u'},
		{"a00", required_argument, NULL, 'a'},
		{"a0", required_argument, NULL, 'Q'},
		{"a1", required_argument, NULL, 'q'},
		{"e00", required_argument, NULL, 'e'},
		{"e0", required_argument, NULL, 'F'},
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
	m000 = FB_M000;
	m001 = FB_M001;
	m01 = FB_M01;
	m10 = FB_M10;
	m11 = FB_M11;
	r000 = FB_R000;
	r001 = FB_R001;
	r01 = FB_R01;
	r10 = FB_R10;
	r11 = FB_R11;
	a00 = FB_A00;
	a0 = FB_A0;
	a1 = FB_A1;
	e00 = FB_E00;
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
			m000 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'n':
			m001 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'o':
			m01 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'p':
			m10 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'l':
			m11 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'r':
			r000 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'g':
			r001 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'i':
			r01 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'j':
			r10 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'u':
			r11 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'a':
			a00 = atof(optarg) * FB_CONST_AU;
			break;
		case 'Q':
			a0 = atof(optarg) * FB_CONST_AU;
			break;
		case 'q':
			a1 = atof(optarg) * FB_CONST_AU;
			break;
		case 'e':
			e00 = atof(optarg);
			if (e00 >= 1.0) {
				fprintf(stderr, "e00 must be less than 1\n");
				return(1);
			}
			break;
		case 'F':
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
	hier.nstarinit = 5;
	hier.nstar = 5;
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
	fprintf(stderr, "  a00=%.6g AU  e00=%.6g  m000=%.6g MSUN  m001=%.6g MSUN  r000=%.6g RSUN  r001=%.6g RSUN\n", \
		a00/FB_CONST_AU, e00, m000/FB_CONST_MSUN, m001/FB_CONST_MSUN, r000/FB_CONST_RSUN, r001/FB_CONST_RSUN);
	fprintf(stderr, "  a0=%.6g AU  e0=%.6g  m01=%.6g MSUN  r01=%.6g RSUN\n", \
		a0/FB_CONST_AU, e0, m01/FB_CONST_MSUN, r01/FB_CONST_RSUN);
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

	/* create hierarchies */
	hier.narr[2] = 2;
	hier.narr[3] = 1;
	/* inner binary of triple */
	hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
	hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
	hier.hier[hier.hi[2]+0].t = t;
	/* outer binary of triple */
	hier.hier[hier.hi[3]+0].obj[0] = &(hier.hier[hier.hi[2]+0]);
	hier.hier[hier.hi[3]+0].obj[1] = &(hier.hier[hier.hi[1]+2]);
	hier.hier[hier.hi[3]+0].t = t;
	/* the other binary */
	hier.hier[hier.hi[2]+1].obj[0] = &(hier.hier[hier.hi[1]+3]);
	hier.hier[hier.hi[2]+1].obj[1] = &(hier.hier[hier.hi[1]+4]);
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

	hier.hier[hier.hi[1]+0].R = r000;
	hier.hier[hier.hi[1]+1].R = r001;
	hier.hier[hier.hi[1]+2].R = r01;
	hier.hier[hier.hi[1]+3].R = r10;
	hier.hier[hier.hi[1]+4].R = r11;

	hier.hier[hier.hi[1]+0].m = m000;
	hier.hier[hier.hi[1]+1].m = m001;
	hier.hier[hier.hi[1]+2].m = m01;
	hier.hier[hier.hi[1]+3].m = m10;
	hier.hier[hier.hi[1]+4].m = m11;

	hier.hier[hier.hi[2]+0].m = m000 + m001;
	hier.hier[hier.hi[3]+0].m = m000 + m001 + m01;
	hier.hier[hier.hi[2]+1].m = m10 + m11;

	hier.hier[hier.hi[2]+0].a = a00;
	hier.hier[hier.hi[3]+0].a = a0;
	hier.hier[hier.hi[2]+1].a = a1;
	
	hier.hier[hier.hi[2]+0].e = e00;
	hier.hier[hier.hi[3]+0].e = e0;
	hier.hier[hier.hi[2]+1].e = e1;

	hier.nobj = 2;
	hier.obj[0] = &(hier.hier[hier.hi[3]+0]);
	hier.obj[1] = &(hier.hier[hier.hi[2]+1]);
	hier.obj[2] = NULL;
	hier.obj[3] = NULL;
	hier.obj[4] = NULL;

	/* get the units and normalize */
	calc_units(hier.obj, &units);
	fb_normalize(&hier, units);
	
	/* temporarily place hierarchies at origin so we can calculate v_c */
	for (i=0; i<2; i++) {
		for (j=0; j<3; j++) {
			hier.obj[i]->x[j] = 0.0;
			hier.obj[i]->v[j] = 0.0;
		}
	}

	/* randomize binary orientations and downsync */
	fb_randorient(&(hier.hier[hier.hi[3]+0]), rng);
	fb_downsync(&(hier.hier[hier.hi[3]+0]), t);
	fb_randorient(&(hier.hier[hier.hi[2]+0]), rng);
	fb_downsync(&(hier.hier[hier.hi[2]+0]), t);
	fb_randorient(&(hier.hier[hier.hi[2]+1]), rng);
	fb_downsync(&(hier.hier[hier.hi[2]+1]), t);

	/* calculate v_c numerically, since it can't be calculated analytically since we have a triple */
	vcrit = sqrt(-(hier.obj[0]->m + hier.obj[1]->m)/(hier.obj[0]->m * hier.obj[1]->m) * 2.0 * \
		(fb_petot(&(hier.hier[hier.hi[1]+0]), 3) + fb_ketot(&(hier.hier[hier.hi[1]+0]), 3) + \
		 fb_petot(&(hier.hier[hier.hi[1]+3]), 2) + fb_ketot(&(hier.hier[hier.hi[1]+3]), 2)));
	
	fprintf(stderr, "UNITS:\n");
	fprintf(stderr, "  v_crit=%.6g km/s  v=%.6g km/s  l=%.6g AU  t=t_dyn=%.6g yr\n", \
		vcrit*units.v/1.0e5, units.v/1.0e5, units.l/FB_CONST_AU, units.t/FB_CONST_YR);
	fprintf(stderr, "  M=%.6g M_sun  E=%.6g erg\n\n", units.m/FB_CONST_MSUN, units.E);
	
	/* move hierarchies analytically in from infinity along hyperbolic orbit */
	fb_dprintf("moving hierarchies analytically in from infinity...\n");

	m0 = hier.obj[0]->m;
	m1 = hier.obj[1]->m;
	M = m0 + m1;
	mu = m0 * m1 / M;

	Ei = 0.5 * mu * fb_sqr(vinf*vcrit);

	a0 = hier.obj[0]->a;
	a1 = hier.obj[1]->a;

	e0 = hier.obj[0]->e;
	e1 = hier.obj[1]->e;
	
	m00 = hier.obj[0]->obj[0]->m;
	m01 = hier.obj[0]->obj[1]->m;
	m10 = hier.obj[1]->obj[0]->m;
	m11 = hier.obj[1]->obj[1]->m;

	rtid = pow(2.0*m0*m1/input.tidaltol, 1.0/3.0) * \
		FB_MAX(pow(m00*m01, -1.0/3.0)*a0*(1.0+e0), pow(m10*m11, -1.0/3.0)*a1*(1.0+e1));
	fb_init_scattering(hier.obj, vinf*vcrit, b, rtid);
	
	/* and check to see that we conserved energy and angular momentum */
	fb_cross(hier.obj[0]->x, hier.obj[0]->v, l0);
	fb_cross(hier.obj[1]->x, hier.obj[1]->v, l1);
	
	for (i=0; i<3; i++) {
		L[i] = (m0 * l0[i] + m1 * l1[i]);
		r[i] = hier.obj[1]->x[i] - hier.obj[0]->x[i];
	}

	E = - m0 * m1 / fb_mod(r) + 0.5 * (m0 * fb_dot(hier.obj[0]->v, hier.obj[0]->v) + \
					m1 * fb_dot(hier.obj[1]->v, hier.obj[1]->v));

	fb_dprintf("L0=%.6g DeltaL/L0=%.6g DeltaL=%.6g\n", mu*b*vinf*vcrit, fb_mod(L)/(mu*b*vinf*vcrit)-1.0, fb_mod(L)-mu*b*vinf*vcrit);
	fb_dprintf("E0=%.6g DeltaE/E0=%.6g DeltaE=%.6g\n\n", Ei, E/Ei-1.0, E-Ei);

	/* trickle down properties */
	fb_trickle(&hier, t);

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
	
	/* free GSL stuff */
	gsl_rng_free(rng);

	/* free our own stuff */
	fb_free_hier(hier);

	/* done! */
	return(0);
}
