/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief calculate the units used
*
* @param obj[2] ?
* @param bs_units ?
*/
void bs_calcunits(fb_obj_t *obj[2], fb_units_t *bs_units)
{
	bs_units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	bs_units->l = obj[1]->a;
	bs_units->t = bs_units->l / bs_units->v;
	bs_units->m = bs_units->l * fb_sqr(bs_units->v) / FB_CONST_G;
	bs_units->E = bs_units->m * fb_sqr(bs_units->v);
}

/**
* @brief ?
*
* @param t ?
* @param ksin ?
* @param kbin ?
* @param W ?
* @param bmax ?
* @param hier ?
* @param rng ?
*
* @return ?
*/
fb_ret_t binsingle(double *t, long ksin, long kbin, double W, double bmax, fb_hier_t *hier, gsl_rng *rng)
{
	int j;
	long jbin;
	double vc, b, rtid, m0, m1, a1, e1, m10, m11;
	fb_units_t fb_units;
	fb_input_t input;
	fb_ret_t retval;
	
	/* a useful definition */
	jbin = star[kbin].binind;

	/* v_inf should be in units of v_crit */
#ifdef USE_MPI
	vc = sqrt(binary[jbin].m1 * binary[jbin].m2 * (star_m[get_global_idx(kbin)] + star_m[get_global_idx(ksin)]) / \
		  (binary[jbin].a * star_m[get_global_idx(kbin)] * star_m[get_global_idx(ksin)] * ((double) clus.N_STAR)));
#else
	vc = sqrt(binary[jbin].m1 * binary[jbin].m2 * (star[kbin].m + star[ksin].m) / \
		  (binary[jbin].a * star[kbin].m * star[ksin].m * ((double) clus.N_STAR)));
#endif

#ifndef USE_MPI
	curr_st = &st[findProcForIndex(ksin)];
#endif
	b = sqrt(rng_t113_dbl_new(curr_st)) * bmax / binary[jbin].a;
	/* b should be in units of a */
	//b = sqrt(rng_t113_dbl()) * bmax / binary[jbin].a;
				
	/* set parameters */
	input.ks = 0;
	input.tstop = 1.0e7;
	input.Dflag = 0;
	input.dt = 0.0;
	input.tcpustop = 60.0;
	input.absacc = 1.0e-9;
	input.relacc = 1.0e-9;
	input.ncount = 500;
	input.tidaltol = 1.0e-5;
	input.firstlogentry[0] = '\0';
	input.fexp = 1.0;
	fb_debug = 0;
	
	/* initialize a few things for integrator */
	*t = 0.0;
	hier->nstar = 3;
	fb_init_hier(hier);

	/* create binary */
	hier->hier[hier->hi[2]+0].obj[0] = &(hier->hier[hier->hi[1]+1]);
	hier->hier[hier->hi[2]+0].obj[1] = &(hier->hier[hier->hi[1]+2]);
	hier->hier[hier->hi[2]+0].t = *t;

	/* give the objects some properties */
	for (j=0; j<hier->nstar; j++) {
		hier->hier[hier->hi[1]+j].ncoll = 1;
		sprintf(hier->hier[hier->hi[1]+j].idstring, "%d", j);
		hier->hier[hier->hi[1]+j].n = 1;
		hier->hier[hier->hi[1]+j].obj[0] = NULL;
		hier->hier[hier->hi[1]+j].obj[1] = NULL;
		hier->hier[hier->hi[1]+j].Lint[0] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[1] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[2] = 0.0;
	}
	
	hier->hier[hier->hi[1]+0].id[0] = star[ksin].id;
	hier->hier[hier->hi[1]+1].id[0] = binary[jbin].id1;
	hier->hier[hier->hi[1]+2].id[0] = binary[jbin].id2;

	if (SS_COLLISION) {
		hier->hier[hier->hi[1]+0].R = star[ksin].rad * units.l;
		hier->hier[hier->hi[1]+1].R = binary[jbin].rad1 * units.l;
		hier->hier[hier->hi[1]+2].R = binary[jbin].rad2 * units.l;
	} else {
		hier->hier[hier->hi[1]+0].R = 0.0;
		hier->hier[hier->hi[1]+1].R = 0.0;
		hier->hier[hier->hi[1]+2].R = 0.0;
	}

#ifdef USE_MPI
	hier->hier[hier->hi[1]+0].m = star_m[get_global_idx(ksin)] * units.mstar;
#else
	hier->hier[hier->hi[1]+0].m = star[ksin].m * units.mstar;
#endif
	hier->hier[hier->hi[1]+1].m = binary[jbin].m1 * units.mstar;
	hier->hier[hier->hi[1]+2].m = binary[jbin].m2 * units.mstar;

	hier->hier[hier->hi[2]+0].m = hier->hier[hier->hi[1]+1].m + hier->hier[hier->hi[1]+2].m;

	hier->hier[hier->hi[1]+0].Eint = star[ksin].Eint * units.E;
	hier->hier[hier->hi[1]+1].Eint = binary[jbin].Eint1 * units.E;
	hier->hier[hier->hi[1]+2].Eint = binary[jbin].Eint2 * units.E;

	hier->hier[hier->hi[2]+0].a = binary[jbin].a * units.l;
	hier->hier[hier->hi[2]+0].e = binary[jbin].e;

	hier->obj[0] = &(hier->hier[hier->hi[1]+0]);
	hier->obj[1] = &(hier->hier[hier->hi[2]+0]);
	hier->obj[2] = NULL;

	/* logging */
	parafprintf(binintfile, "********************************************************************************\n");
	parafprintf(binintfile, "type=BS t=%.9g\n", TotalTime);
	parafprintf(binintfile, "params: b=%g v=%g\n", b, W/vc);
	/* set units to 1 since we're already in CGS */
	fb_units.v = fb_units.l = fb_units.t = fb_units.m = fb_units.E = 1.0;
	parafprintf(binintfile, "input: ");
	binint_log_obj(hier->obj[0], fb_units);
	parafprintf(binintfile, "input: ");
	binint_log_obj(hier->obj[1], fb_units);

	/* get the units and normalize */
	bs_calcunits(hier->obj, &fb_units);
	fb_normalize(hier, fb_units);
	
	/* move objects in from infinity along hyperbolic orbit */
	m0 = hier->obj[0]->m;
	m1 = hier->obj[1]->m;
	a1 = hier->obj[1]->a;
	e1 = hier->obj[1]->e;
	m10 = hier->obj[1]->obj[0]->m;
	m11 = hier->obj[1]->obj[1]->m;

	rtid = pow(2.0*m0*m1/input.tidaltol, 1.0/3.0) * pow(m10*m11, -1.0/3.0)*a1*(1.0+e1);

	fb_init_scattering(hier->obj, W/vc, b, rtid);
	
#ifndef USE_MPI
	curr_st = &st[findProcForIndex(ksin)];
#endif

	/* trickle down the binary properties, then back up */
	fb_randorient(&(hier->hier[hier->hi[2]+0]), rng, curr_st);
	fb_downsync(&(hier->hier[hier->hi[2]+0]), *t);
	fb_upsync(&(hier->hier[hier->hi[2]+0]), *t);
	
	/* call fewbody! */
	retval = fewbody(input, hier, t);

	/* and return */
	return(retval);
}
