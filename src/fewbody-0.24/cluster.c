/* -*- linux-c -*- */
/* cluster.c

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
#include <gsl/gsl_roots.h>
#include "fewbody.h"
#include "cluster.h"

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  cluster [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -n --n <n>                    : set number of stars in cluster [%d]\n", FB_N);
	fprintf(stream, "  -m --m <m/MSUN>               : set mass of each star [%.6g]\n", FB_M/FB_CONST_MSUN);
	fprintf(stream, "  -r --r <r/RSUN>               : set radius of each star [%.6g]\n", FB_R/FB_CONST_RSUN);
	fprintf(stream, "  -S --sigma <sigma/(km/s)>     : set central velocity dispersion [%.6g]\n", FB_SIGMA/1.0e5);
	fprintf(stream, "  -T --rmax <rmax/PARSEC>       : set truncation radius [%.6g]\n", FB_RMAX/FB_CONST_PARSEC);
	fprintf(stream, "  -t --tstop <tstop/t_dyn>      : set stopping time [%.6g]\n", FB_TSTOP);
	fprintf(stream, "  -P --tphysstop <tphysstop/yr> : set physical stopping time [%.6g]\n", FB_TPHYSSTOP/FB_CONST_YR);
	fprintf(stream, "  -D --dt <dt/t_dyn>            : set approximate output dt [%.6g]\n", FB_DT);
	fprintf(stream, "  -c --tcpustop <tcpustop/sec>  : set cpu stopping time [%.6g]\n", FB_TCPUSTOP);
	fprintf(stream, "  -A --absacc <absacc>          : set integrator's absolute accuracy [%.6g]\n", FB_ABSACC);
	fprintf(stream, "  -R --relacc <relacc>          : set integrator's relative accuracy [%.6g]\n", FB_RELACC);
	fprintf(stream, "  -N --ncount <ncount>          : set number of integration steps between calls\n");
	fprintf(stream, "                                  to fb_classify() [%d]\n", FB_NCOUNT);
	fprintf(stream, "  -z --tidaltol <tidaltol>      : set tidal tolerance [%.6g]\n", FB_TIDALTOL);
	fprintf(stream, "  -x --fexp <f_exp>             : set expansion factor of merger product [%.6g]\n", FB_FEXP);
	fprintf(stream, "  -k --ks                       : turn K-S regularization on or off [%d]\n", FB_KS);
	fprintf(stream, "  -s --seed                     : set random seed [%ld]\n", FB_SEED);
	fprintf(stream, "  -d --debug                    : turn on debugging\n");
	fprintf(stream, "  -V --version                  : print version info\n");
	fprintf(stream, "  -h --help                     : display this help text\n");
}

/* calculate the units used (N-body units) */
void calc_units(fb_hier_t hier, fb_units_t *units)
{
	int i, j, k;
	double petot, ketot, r[3];

	/* the unit of mass is defined so that M_tot=1 */
	units->m = 0.0;
	for (i=0; i<hier.nstar; i++) {
		units->m += hier.hier[hier.hi[1]+i].m;
	}

	/* calculate total potential and kinetic energy */
	petot = 0.0;
	for (i=0; i<hier.nstar-1; i++) {
		for (j=i+1; j<hier.nstar; j++) {
			for (k=0; k<3; k++) {
				r[k] = hier.hier[hier.hi[1]+j].x[k] - hier.hier[hier.hi[1]+i].x[k];
			}
			petot += - FB_CONST_G * hier.hier[hier.hi[1]+i].m * hier.hier[hier.hi[1]+j].m / fb_mod(r);
		}
	}

	ketot = 0.0;
	for (i=0; i<hier.nstar; i++) {
		ketot += 0.5 * hier.hier[hier.hi[1]+i].m * fb_dot(hier.hier[hier.hi[1]+i].v, hier.hier[hier.hi[1]+i].v);
	}
	
	/* the unit of energy is defined so that E=-1/4 */
	units->E = -4.0 * (petot + ketot);
	
	/* with G=1, all other units are derived */
	units->t = FB_CONST_G * pow(units->m, 2.5) * pow(units->E, -1.5);
	units->l = FB_CONST_G * fb_sqr(units->m) / units->E;
	units->v = units->l / units->t;
}

/* f_0-f(v/v_esc) for the Plummer model */
double fv(double v, void *params)
{
	double f;
	
	f = ((double *)params)[0];
	return(f - 512.0/(7.0*FB_CONST_PI) * (sqrt(1.0-fb_sqr(v))*(-7.0*v/256.0 + 121.0*fb_cub(v)/384.0 - 263.0*fb_sqr(v)*fb_cub(v)/480.0 + 31.0*v*fb_sqr(fb_cub(v))/80.0 - fb_cub(fb_cub(v))/10.0) + 7.0*asin(v)/256.0));
}

/* get speed from Plummer distribution */
double vf(double f)
{
	int status, iter;
	double v, params[1];
	gsl_function F;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	
	/* set up the root solver */
	F.function = &fv;
	F.params = &params;

	/* set the parameters */
	params[0] = f;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, 0.0, 1.0);
	
	/* get v/v_esc by root-finding */
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
	v = gsl_root_fsolver_root(s);
	
	/* free memory associated with root solver */
	gsl_root_fsolver_free(s);

	return(v);
}

/* the main attraction */
int main(int argc, char *argv[])
{
	int i, j, n;
	unsigned long int seed;
	double tphysstop, m, r, Ei, Li[3], Lint[3], t, M, a, sigma, rmax, Xmax, radius, v, theta, phi, xcm[3], vcm[3];
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;
	const char *short_opts = "n:m:r:S:T:t:P:D:c:A:R:N:z:x:k:s:dVh";
	const struct option long_opts[] = {
		{"n", required_argument, NULL, 'n'},
		{"m", required_argument, NULL, 'm'},
		{"r", required_argument, NULL, 'r'},
		{"sigma", required_argument, NULL, 'S'},
		{"rmax", required_argument, NULL, 'T'},
		{"tstop", required_argument, NULL, 't'},
		{"tphysstop", required_argument, NULL, 'P'},
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
	n = FB_N;
	m = FB_M;
	r = FB_R;
	sigma = FB_SIGMA;
	rmax = FB_RMAX;
	input.ks = FB_KS;
	input.tstop = FB_TSTOP;
	tphysstop = FB_TPHYSSTOP;
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
		case 'n':
			n = atoi(optarg);
			break;
		case 'm':
			m = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'r':
			r = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'S':
			sigma = atof(optarg) * 1.0e5;
			break;
		case 'T':
			rmax = atof(optarg) * FB_CONST_PARSEC;
			break;
		case 't':
			input.tstop = atof(optarg);
			break;
		case 'P':
			tphysstop = atof(optarg) * FB_CONST_YR;
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

	/* set up parameters for Plummer model */
	M = ((double) n) * m;
	a = FB_CONST_G * M / (6.0*fb_sqr(sigma));
	Xmax = pow(1.0+fb_sqr(a/rmax), -1.5);
	
	/* initialize a few things for integrator */
	t = 0.0;
	hier.nstarinit = n;
	hier.nstar = n;
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
	fprintf(stderr, "  M=%.6g MSUN  sigma=%.6g km/s  a=%.6g PARSEC  rmax=%.6g PARSEC  rmax/a=%.6g\n",
		M/FB_CONST_MSUN, sigma/1.0e5, a/FB_CONST_PARSEC, rmax/FB_CONST_PARSEC, rmax/a);
	fprintf(stderr, "  n=%d  m=%.6g MSUN  r=%.6g RSUN  tstop=%.6g  tphysstop=%.6g yr  tcpustop=%.6g\n", \
		n, m/FB_CONST_MSUN, r/FB_CONST_RSUN, input.tstop, tphysstop/FB_CONST_YR, input.tcpustop);
	fprintf(stderr, "  tidaltol=%.6g  abs_acc=%.6g  rel_acc=%.6g  ncount=%d  fexp=%.6g\n\n", \
		input.tidaltol, input.absacc, input.relacc, input.ncount, input.fexp);

	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, seed);

	/* prepare to calculate the center of mass position and velocity */
	for (i=0; i<3; i++) {
		xcm[i] = 0.0;
		vcm[i] = 0.0;
	}

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
		hier.hier[hier.hi[1]+j].m = m;
		hier.hier[hier.hi[1]+j].R = r;
		
		/* draw radius from Plummer distribution */
		radius = a / sqrt(pow(Xmax*gsl_rng_uniform(rng), -2.0/3.0)-1.0);
		theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0);
		phi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
		hier.hier[hier.hi[1]+j].x[0] = radius * sin(theta) * cos(phi);
		hier.hier[hier.hi[1]+j].x[1] = radius * sin(theta) * sin(phi);
		hier.hier[hier.hi[1]+j].x[2] = radius * cos(theta);

		/* draw velocity from Plummer distribution */
		v = vf(gsl_rng_uniform(rng)) * sqrt(2.0*FB_CONST_G*M/(a*sqrt(1.0+fb_sqr(radius/a))));
		theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0);
		phi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
		hier.hier[hier.hi[1]+j].v[0] = v * sin(theta) * cos(phi);
		hier.hier[hier.hi[1]+j].v[1] = v * sin(theta) * sin(phi);
		hier.hier[hier.hi[1]+j].v[2] = v * cos(theta);
		
		/* calculate center of mass position and velocity */
		for (i=0; i<3; i++) {
			xcm[i] += hier.hier[hier.hi[1]+j].x[i] / ((double) n);
			vcm[i] += hier.hier[hier.hi[1]+j].v[i] / ((double) n);
		}
	}
	
	/* transform to center of mass frame */
	for (j=0; j<hier.nstar; j++) {
		for (i=0; i<3; i++) {
			hier.hier[hier.hi[1]+j].x[i] -= xcm[i];
			hier.hier[hier.hi[1]+j].v[i] -= vcm[i];
		}
	}

	/* get the units and normalize */
	calc_units(hier, &units);
	fb_normalize(&hier, units);
	
	/* reset input.tstop if necessary, based on the value of tphysstop */
	if (tphysstop/units.t < input.tstop) {
		input.tstop = tphysstop/units.t;
		fb_dprintf("decreasing tstop in accord with tphysstop: tstop=%.6g\n", input.tstop);
	}

	fb_dprintf("virial ratio: q=%.6g\n", -fb_ketot(&(hier.hier[hier.hi[1]]), hier.nstar)/fb_petot(&(hier.hier[hier.hi[1]]), hier.nstar));

	fprintf(stderr, "UNITS:\n");
	fprintf(stderr, "  v=%.6g km/s  l=%.6g AU  t=t_dyn=%.6g yr\n", \
		units.v/1.0e5, units.l/FB_CONST_AU, units.t/FB_CONST_YR);
	fprintf(stderr, "  M=%.6g M_sun  E=%.6g erg\n\n", units.m/FB_CONST_MSUN, units.E);

	/* calculate the initial energy and angular momentum, for bookkeeping */
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
	// PAU retval = fewbody(input, &hier, &t);
	retval = fewbody(input, units, &hier, &t, rng);

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
