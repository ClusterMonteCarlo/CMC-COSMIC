/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

/* calculate the units used */
void bs_calcunits(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	units->l = obj[1]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
}

/* the input units are CGS */
fb_ret_t binsingle(double *t, double m0, double m10, double m11, double r0, double r10, double r11, double a1, double e1, double vinf, double b, fb_units_t *units, fb_hier_t *hier, gsl_rng *rng)
{
	int j;
	double rtid, m1;
	fb_input_t input;
	fb_ret_t retval;

	/* sanity check */
	if (m0 <= 0.0 || m10 <= 0.0 || m11 <= 0.0) {
		eprintf("unphysical mass: m0=%g m10=%g m11=%g\n", m0, m10, m11);
		exit_cleanly(1);
	} else if (r0 < 0.0 || r10 < 0.0 || r11 < 0.0) {
		eprintf("unphysical radius: r0=%g r10=%g r11=%g\n", r0, r10, r11);
		exit_cleanly(1);
	} else if (a1 <= 0.0) {
		eprintf("unphysical semimajor axis: a1=%g\n", a1);
		exit_cleanly(1);
	} else if (e1 < 0.0 || e1 >= 1.0) {
		eprintf("unphysical eccentricity: e1=%g\n", e1);
		exit_cleanly(1);
	} else if (vinf < 0.0 || b < 0.0) {
		eprintf("unphysical scattering parameters: vinf=%g b=%g\n", vinf, b);
		exit_cleanly(1);
	}

	/* set parameters */
	input.ks = 0;
	input.tstop = 1.0e9;
	input.Dflag = 0;
	input.dt = 0.0;
	input.tcpustop = 60.0;
	input.absacc = 1.0e-9;
	input.relacc = 1.0e-9;
	input.ncount = 500;
	input.tidaltol = 1.0e-5;
	input.firstlogentry[0] = '\0';
	input.fexp = 3.0;
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
		hier->hier[hier->hi[1]+j].id[0] = j;
		sprintf(hier->hier[hier->hi[1]+j].idstring, "%d", j);
		hier->hier[hier->hi[1]+j].n = 1;
		hier->hier[hier->hi[1]+j].obj[0] = NULL;
		hier->hier[hier->hi[1]+j].obj[1] = NULL;
		hier->hier[hier->hi[1]+j].Eint = 0.0;
		hier->hier[hier->hi[1]+j].Lint[0] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[1] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[2] = 0.0;
	}

	hier->hier[hier->hi[1]+0].R = r0;
	hier->hier[hier->hi[1]+1].R = r10;
	hier->hier[hier->hi[1]+2].R = r11;

	hier->hier[hier->hi[1]+0].m = m0;
	hier->hier[hier->hi[1]+1].m = m10;
	hier->hier[hier->hi[1]+2].m = m11;

	hier->hier[hier->hi[2]+0].m = m10 + m11;

	hier->hier[hier->hi[2]+0].a = a1;
	hier->hier[hier->hi[2]+0].e = e1;

	hier->obj[0] = &(hier->hier[hier->hi[1]+0]);
	hier->obj[1] = &(hier->hier[hier->hi[2]+0]);
	hier->obj[2] = NULL;

	/* get the units and normalize */
	bs_calcunits(hier->obj, units);
	fb_normalize(hier, *units);
	
	/* move objects in from infinity along hyperbolic orbit */
	m0 = hier->obj[0]->m;
	m1 = hier->obj[1]->m;
	a1 = hier->obj[1]->a;
	e1 = hier->obj[1]->e;
	m10 = hier->obj[1]->obj[0]->m;
	m11 = hier->obj[1]->obj[1]->m;

	rtid = pow(2.0*m0*m1/input.tidaltol, 1.0/3.0) * pow(m10*m11, -1.0/3.0)*a1*(1.0+e1);

	fb_init_scattering(hier->obj, vinf, b, rtid);
	
	/* trickle down the binary properties, then back up */
	fb_randorient(&(hier->hier[hier->hi[2]+0]), rng);
	fb_downsync(&(hier->hier[hier->hi[2]+0]), *t);
	fb_upsync(&(hier->hier[hier->hi[2]+0]), *t);
	
	/* call fewbody! */
	retval = fewbody(input, hier, t);

	/* and return */
	return(retval);
}
