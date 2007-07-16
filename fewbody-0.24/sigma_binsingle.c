/* -*- linux-c -*- */
/* sigma_binsingle.c

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
#include "sigma_binsingle.h"

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  sigma_binsingle [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -m --m0 <m0/MSUN>             : set mass of single star [%.6g]\n", FB_M0/FB_CONST_MSUN);
	fprintf(stream, "  -n --m10 <m10/MSUN>           : set mass of star 0 of binary [%.6g]\n", FB_M10/FB_CONST_MSUN);
	fprintf(stream, "  -o --m11 <m11/MSUN>           : set mass of star 1 of binary [%.6g]\n", FB_M11/FB_CONST_MSUN);
	fprintf(stream, "  -r --r0 <r0/RSUN>             : set radius of single star [%.6g]\n", FB_R0/FB_CONST_RSUN);
	fprintf(stream, "  -g --r10 <r10/RSUN>           : set radius of star 0 of binary [%.6g]\n", FB_R10/FB_CONST_RSUN);
	fprintf(stream, "  -i --r11 <r11/RSUN>           : set radius of star 1 of binary [%.6g]\n", FB_R11/FB_CONST_RSUN);
	fprintf(stream, "  -a --a1 <a1/AU>               : set semimajor axis of binary [%.6g]\n", FB_A1/FB_CONST_AU);
	fprintf(stream, "  -e --e1 <e1>                  : set eccentricity of binary 0 [%.6g]\n", FB_E1);
	fprintf(stream, "  -v --vinf <vinf/v_crit>       : set velocity at infinity [%.6g]\n", FB_VINF);
	fprintf(stream, "  -P --precision <dsigma/sigma> : set fractional precision on cross\n");
	fprintf(stream, "                                  section [%.6g]\n", FB_PRECISION);
	fprintf(stream, "  -t --tstop <tstop/t_dyn>      : set stopping time [%.6g]\n", FB_TSTOP);
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

/* calculate the units used */
int calc_units(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	units->l = obj[1]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
	
	return(0);
}

/* the main attraction */
int main(int argc, char *argv[])
{
	int i, j, k, res, err, sid, bid, precisionmet=0;
	long l, nb;
	unsigned long int seed;
	double m0, m10, m11, r0, r10, r11, a1, e1, precision;
	double rtid, vinf, rperi, b, b0, bmin, bmax, db, bxlast, t, vc;
	fb_sigma_t sigmacurr[11], sigmaprev[11];
	char outcome[11][FB_MAX_STRING_LENGTH];
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;
	const char *short_opts = "m:n:o:r:g:i:a:e:v:P:t:D:c:A:R:N:z:x:k:s:dVh";
	const struct option long_opts[] = {
		{"m0", required_argument, NULL, 'm'},
		{"m10", required_argument, NULL, 'n'},
		{"m11", required_argument, NULL, 'o'},
		{"r0", required_argument, NULL, 'r'},
		{"r10", required_argument, NULL, 'g'},
		{"r11", required_argument, NULL, 'i'},
		{"a1", required_argument, NULL, 'a'},
		{"e1", required_argument, NULL, 'e'},
		{"vinf", required_argument, NULL, 'v'},
		{"precision", required_argument, NULL, 'P'},
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
	m0 = FB_M0;
	m10 = FB_M10;
	m11 = FB_M11;
	r0 = FB_R0;
	r10 = FB_R10;
	r11 = FB_R11;
	a1 = FB_A1;
	e1 = FB_E1;
	vinf = FB_VINF;
	precision = 1.0e-2;
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
			m0 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'n':
			m10 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'o':
			m11 = atof(optarg) * FB_CONST_MSUN;
			break;
		case 'r':
			r0 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'g':
			r10 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'i':
			r11 = atof(optarg) * FB_CONST_RSUN;
			break;
		case 'a':
			a1 = atof(optarg) * FB_CONST_AU;
			break;
		case 'e':
			e1 = atof(optarg);
			if (e1 >= 1.0) {
				fprintf(stderr, "e0 must be less than 1\n");
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
		case 'P':
			precision = atof(optarg);
			if (precision <= 0.0) {
				fprintf(stderr, "precision must be > 0\n");
				return(1);
			} else if (precision > 1.0) {
				fprintf(stderr, "precision must be < 1\n");
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
	hier.nstarinit = 3;
	fb_malloc_hier(&hier);

	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, seed);

	/* zero everything out before summing */
	for (i=0; i<11; i++) {
		sigmacurr[i].sigma_r = 0.0;
		sigmacurr[i].dsigma2_r_minus = 0.0;
		sigmacurr[i].dsigma2_r_plus = 0.0;
		sigmacurr[i].sigma_nr = 0.0;
		sigmacurr[i].dsigma2_nr_minus = 0.0;
		sigmacurr[i].dsigma2_nr_plus = 0.0;
		sigmaprev[i].sigma_r = 0.0;
		sigmaprev[i].dsigma2_r_minus = 0.0;
		sigmaprev[i].dsigma2_r_plus = 0.0;
		sigmaprev[i].sigma_nr = 0.0;
		sigmaprev[i].dsigma2_nr_minus = 0.0;
		sigmaprev[i].dsigma2_nr_plus = 0.0;
	}
	/* loop through and calcluate cross sections */
	nb = (long) 1.0/precision;
	while (!precisionmet) {
		/* select b0 that gives r_p=2a */
		vc = sqrt(FB_CONST_G * (m0+m10+m11) / (m0*(m10+m11)) * (m10*m11/a1));
		rperi = 2.0 * a1;
		/* b0 is in code units here (vinf is in units of vc) */
		b0 = sqrt(fb_sqr(rperi) + 2.0 * FB_CONST_G * (m0+m10+m11) * rperi / fb_sqr(vinf*vc)) / a1;
		db = b0/((double) nb);
		bxlast = GSL_POSINF;

		l = 0;
		bmin = ((double) l) * db;
		bmax = bmin + db;
		while (bmax <= b0 || bmax <= 2.0*bxlast) {
			bmin = ((double) l) * db;
			bmax = bmin + db;
			/* choose b uniformly in area from bmin to bmax */
			b = sqrt(gsl_rng_uniform(rng) * (fb_sqr(bmax)-fb_sqr(bmin)) + fb_sqr(bmin));
			
			/* do scattering experiment */
			/* initialize a few things for integrator */
			t = 0.0;
			hier.nstar = 3;
			fb_init_hier(&hier);
			
			/* create binary */
			hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+1]);
			hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+2]);
			hier.hier[hier.hi[2]+0].t = t;
			
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
			
			hier.hier[hier.hi[1]+0].R = r0;
			hier.hier[hier.hi[1]+1].R = r10;
			hier.hier[hier.hi[1]+2].R = r11;
			
			hier.hier[hier.hi[1]+0].m = m0;
			hier.hier[hier.hi[1]+1].m = m10;
			hier.hier[hier.hi[1]+2].m = m11;
			
			hier.hier[hier.hi[2]+0].m = m10 + m11;
			
			hier.hier[hier.hi[2]+0].a = a1;
			hier.hier[hier.hi[2]+0].e = e1;
			
			hier.obj[0] = &(hier.hier[hier.hi[1]+0]);
			hier.obj[1] = &(hier.hier[hier.hi[2]+0]);
			hier.obj[2] = NULL;
			
			/* get the units and normalize */
			calc_units(hier.obj, &units);
			fb_normalize(&hier, units);
			
			/* move hierarchies analytically in from infinity along hyperbolic orbit */
			rtid = pow(2.0*(hier.obj[0]->m + hier.obj[1]->m) / (hier.obj[1]->m * input.tidaltol), 1.0/3.0) * 
				hier.obj[1]->a * (1.0 + hier.obj[1]->e);
			fb_init_scattering(hier.obj, vinf, b, rtid);
			
			/* trickle down the binary properties, then back up */
			fb_randorient(&(hier.hier[hier.hi[2]+0]), rng);
			fb_downsync(&(hier.hier[hier.hi[2]+0]), t);
			/* fb_upsync(&(hier.hier[hier.hi[2]+0]), t); */
			
			/* call fewbody! */
			retval = fewbody(input, &hier, &t);

			k = -1;
			res = 0;
			err = 0;

			/* analyze outcome */
			if (retval.retval == 1 && 
			    FB_MIN(fabs(retval.DeltaEfrac), fabs(retval.DeltaE)) <= 1.0e-4 &&
			    FB_MIN(fabs(retval.DeltaLfrac), fabs(retval.DeltaL)) <= 1.0e-4) {
				/* check for resonance */
				if (retval.Nosc >= 1) {
					res = 1;
				} else {
					res = 0;
				}
				/* classify outcome */
				if (hier.nstar == 3) {
					if (hier.nobj == 1) {
						/* triples physically forbidden */
						err = 1;
					} else if (hier.nobj == 2) {
						if (fb_n_hier(hier.obj[0])==2) {
							bid = 0;
							sid = 1;
						} else {
							bid = 1;
							sid = 0;
						}
						if (hier.obj[sid]->id[0] == 0) {
							k = 0;
						} else if (hier.obj[sid]->id[0] == 1) {
							k = 1;
						} else {
							k = 2;
						}
					} else { /* hier.nobj == 3 */
						k = 3;
					}
                                } else if (hier.nstar == 2) {
					if (hier.nobj == 1) {
						if ((hier.obj[0]->obj[0]->ncoll==1 && hier.obj[0]->obj[0]->id[0]==0) ||
						    (hier.obj[0]->obj[1]->ncoll==1 && hier.obj[0]->obj[1]->id[0]==0)) {
							k = 4;
						} else if ((hier.obj[0]->obj[0]->ncoll==1 && hier.obj[0]->obj[0]->id[0]==1) ||
							   (hier.obj[0]->obj[1]->ncoll==1 && hier.obj[0]->obj[1]->id[0]==1)) {
							k = 5;
						} else {
							k = 6;
						}
					} else { /* hier.nobj == 2 */
						if ((hier.obj[0]->ncoll==1 && hier.obj[0]->id[0]==0) ||
						    (hier.obj[1]->ncoll==1 && hier.obj[1]->id[0]==0)) {
							k = 7;
						} else if ((hier.obj[0]->ncoll==1 && hier.obj[0]->id[0]==1) ||
							   (hier.obj[1]->ncoll==1 && hier.obj[1]->id[0]==1)) {
							k = 8;
						} else {
							k = 9;
						}
					}
				} else if (hier.nstar == 1) {
					k = 10;
				} else {
					fprintf(stderr, "ERROR: hier.nstar!=1,2,3!\n");
					exit(1);
				}
                        } else { /* bad outcome */
				err = 1;
			}

			/* tally up cross sections */
			if (!err) {
				for (i=0; i<11; i++) {
					if (i==k) {
						if (res) {
							sigmacurr[i].sigma_r += FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin));
							sigmacurr[i].dsigma2_r_minus += fb_sqr(0.1587 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
							sigmacurr[i].dsigma2_nr_plus += fb_sqr(0.8413/21.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
						} else {
							sigmacurr[i].sigma_nr += FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin));
							sigmacurr[i].dsigma2_nr_minus += fb_sqr(0.1587 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));	
							sigmacurr[i].dsigma2_r_plus += fb_sqr(0.8413/21.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
						}
					} else {
						sigmacurr[i].dsigma2_r_plus += fb_sqr(0.8413/21.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
						sigmacurr[i].dsigma2_nr_plus += fb_sqr(0.8413/21.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
					}
				}
			} else {
				for (i=0; i<11; i++) {
					sigmacurr[i].dsigma2_r_plus += fb_sqr(0.8413/22.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
					sigmacurr[i].dsigma2_nr_plus += fb_sqr(0.8413/22.0 * FB_CONST_PI * (fb_sqr(bmax)-fb_sqr(bmin)));
				}	
			}
			
			/* if the outcome was not a weak flyby, then record last impact parameter */
			if (!err && (k!=0 || (k==0 && res))) {
				bxlast = bmax;
			}
			
			/*
			fprintf(stderr, "outcome:  retval=%d  %s  (%s)\n", retval.retval, 
				fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2));
			fprintf(stderr, "  DeltaL/L0=%.6g  DeltaL=%.6g\n", retval.DeltaLfrac, retval.DeltaL);
			fprintf(stderr, "  DeltaE/E0=%.6g  DeltaE=%.6g\n", retval.DeltaEfrac, retval.DeltaE);
			fprintf(stderr, "  l=%ld b0=%g bmin=%g bmax=%g bxlast=%g k=%d res=%d err=%d\n", l, b0, bmin, bmax, bxlast, k, res, err);
			*/

			l++;
		}
		
		/* add current and previous sigmas, and transfer current to previous */
		
		/* increase resolution geometrically */
		nb *= 2;
		
		precisionmet = 1;
	}

	snprintf(outcome[0], FB_MAX_STRING_LENGTH, "[1 2] 0 (preservation)");
	snprintf(outcome[1], FB_MAX_STRING_LENGTH, "[0 2] 1 (exchange_1)");
	snprintf(outcome[2], FB_MAX_STRING_LENGTH, "[0 1] 2 (exchange_2)");
	snprintf(outcome[3], FB_MAX_STRING_LENGTH, "0 1 2 (ionization)");
	snprintf(outcome[4], FB_MAX_STRING_LENGTH, "[1:2 0] (merger_binary_12)");
	snprintf(outcome[5], FB_MAX_STRING_LENGTH, "[0:2 1] (merger_binary_02)");
	snprintf(outcome[6], FB_MAX_STRING_LENGTH, "[0:1 2] (merger_binary_01)");
	snprintf(outcome[7], FB_MAX_STRING_LENGTH, "1:2 0 (merger_12)");
	snprintf(outcome[8], FB_MAX_STRING_LENGTH, "0:2 1 (merger_02)");
	snprintf(outcome[9], FB_MAX_STRING_LENGTH, "0:1 2 (merger_01)");
	snprintf(outcome[10], FB_MAX_STRING_LENGTH, "0:1:2 (triple_merger)");

	/* precision has now been met, so print out results */
	fprintf(stdout, "outcome  sigma_r  dsigma_r_minus  dsigma_r_plus  sigma_nr  dsigma_nr_minus  dsigma_nr_plus\n");
	for (i=0; i<11; i++) {
		fprintf(stdout, "%27s  %6.6g  %6.6g  %6.6g  %6.6g  %6.6g  %6.6g\n", outcome[i], 
			sigmacurr[i].sigma_r/FB_CONST_PI, 
			sqrt(sigmacurr[i].dsigma2_r_minus)/FB_CONST_PI, 
			sqrt(sigmacurr[i].dsigma2_r_plus)/FB_CONST_PI,
			sigmacurr[i].sigma_nr/FB_CONST_PI, 
			sqrt(sigmacurr[i].dsigma2_nr_minus)/FB_CONST_PI, 
			sqrt(sigmacurr[i].dsigma2_nr_plus)/FB_CONST_PI);
	}

	/* free GSL stuff */
	gsl_rng_free(rng);

	/* free our own stuff */
	fb_free_hier(hier);

	/* done! */
	return(0);
}
