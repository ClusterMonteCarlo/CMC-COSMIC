/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

/* calculate the units used */
void bb_calcunits(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[0]->obj[0]->m * obj[0]->obj[1]->m / obj[0]->a + \
			 obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	units->l = obj[0]->a + obj[1]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
}

/* the main attraction */
fb_ret_t binbin(double *t, double m00, double m01, double m10, double m11, double r00, double r01, double r10, double r11, double a0, double a1, double e0, double e1, double vinf, double b, fb_units_t *units, fb_hier_t *hier, gsl_rng *rng)
{
	int j;
	double rtid, m0, m1;
	fb_input_t input;
	fb_ret_t retval;

	/* sanity check */
	if (m00 <= 0.0 || m01 <= 0.0 || m10 <= 0.0 || m11 <= 0.0) {
		eprintf("unphysical mass: m00=%g m01=%g m10=%g m11=%g\n", m00, m01, m10, m11);
		exit_cleanly(1);
	} else if (r00 < 0.0 || r01 < 0.0 || r10 < 0.0 || r11 < 0.0) {
		eprintf("unphysical radius: r00=%g r01=%g r10=%g r11=%g\n", r00, r01, r10, r11);
		exit_cleanly(1);
	} else if (a0 <= 0.0 || a1 <= 0.0) {
		eprintf("unphysical semimajor axis: a0=%g a1=%g\n", a0, a1);
		exit_cleanly(1);
	} else if (e0 < 0.0 || e0 >= 1.0 || e1 < 0.0 || e1 >= 1.0) {
		eprintf("unphysical eccentricity: e0=%g e1=%g\n", e0, e1);
		exit_cleanly(1);
	} else if (vinf < 0.0 || b < 0.0) {
		eprintf("unphysical scattering parameters: vinf=%g b=%g\n", vinf, b);
		exit_cleanly(1);
	}

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
	input.fexp = 3.0;
	fb_debug = 0;
	
	/* initialize a few things for integrator */
	*t = 0.0;
	hier->nstar = 4;
	fb_init_hier(hier);

	/* create binaries */
	hier->hier[hier->hi[2]+0].obj[0] = &(hier->hier[hier->hi[1]+0]);
	hier->hier[hier->hi[2]+0].obj[1] = &(hier->hier[hier->hi[1]+1]);
	hier->hier[hier->hi[2]+0].t = *t;
	hier->hier[hier->hi[2]+1].obj[0] = &(hier->hier[hier->hi[1]+2]);
	hier->hier[hier->hi[2]+1].obj[1] = &(hier->hier[hier->hi[1]+3]);
	hier->hier[hier->hi[2]+1].t = *t;

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

	hier->hier[hier->hi[1]+0].R = r00;
	hier->hier[hier->hi[1]+1].R = r01;
	hier->hier[hier->hi[1]+2].R = r10;
	hier->hier[hier->hi[1]+3].R = r11;

	hier->hier[hier->hi[1]+0].m = m00;
	hier->hier[hier->hi[1]+1].m = m01;
	hier->hier[hier->hi[1]+2].m = m10;
	hier->hier[hier->hi[1]+3].m = m11;

	hier->hier[hier->hi[2]+0].m = m00 + m01;
	hier->hier[hier->hi[2]+1].m = m10 + m11;

	hier->hier[hier->hi[2]+0].a = a0;
	hier->hier[hier->hi[2]+1].a = a1;
	
	hier->hier[hier->hi[2]+0].e = e0;
	hier->hier[hier->hi[2]+1].e = e1;

	hier->obj[0] = &(hier->hier[hier->hi[2]+0]);
	hier->obj[1] = &(hier->hier[hier->hi[2]+1]);
	hier->obj[2] = NULL;
	hier->obj[3] = NULL;

	/* get the units and normalize */
	bb_calcunits(hier->obj, units);
	fb_normalize(hier, *units);
	
	m0 = hier->obj[0]->m;
	m1 = hier->obj[1]->m;
	a0 = hier->obj[0]->a;
	a1 = hier->obj[1]->a;
	e0 = hier->obj[0]->e;
	e1 = hier->obj[1]->e;
	m00 = hier->obj[0]->obj[0]->m;
	m01 = hier->obj[0]->obj[1]->m;
	m10 = hier->obj[1]->obj[0]->m;
	m11 = hier->obj[1]->obj[1]->m;

	rtid = pow(2.0*m0*m1/input.tidaltol, 1.0/3.0) * \
		FB_MAX(pow(m00*m01, -1.0/3.0)*a0*(1.0+e0), pow(m10*m11, -1.0/3.0)*a1*(1.0+e1));

	fb_init_scattering(hier->obj, vinf, b, rtid);
	
	/* trickle down the binary properties, then back up */
	for (j=0; j<2; j++) {
		fb_randorient(&(hier->hier[hier->hi[2]+j]), rng);
		fb_downsync(&(hier->hier[hier->hi[2]+j]), *t);
		fb_upsync(&(hier->hier[hier->hi[2]+j]), *t);
	}
	
	/* call fewbody! */
	retval = fewbody(input, hier, t);

	/* and return */
	return(retval);
}
