/* -*- linux-c -*- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <fitsio.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#include "cmc_vars.h"

void load_fits_file_data(void)
{
	long i, j;

	newstarid = 0;

	/* set units */
	units_set();

	/* copy everything over from cfd */
	for (i=0; i<=cfd->NOBJ+1; i++) {
		star[i].id = cfd->obj_id[i];
		if (star[i].id > newstarid) {
			newstarid = star[i].id;
		}
		star[i].k = cfd->obj_k[i];
		star[i].m = cfd->obj_m[i] * ((double) clus.N_STAR);
		star[i].rad = cfd->obj_Reff[i];
		star[i].r = cfd->obj_r[i];
		star[i].vr = cfd->obj_vr[i];
		star[i].vt = cfd->obj_vr[i];
		star[i].binind = cfd->obj_binind[i];
	}

	for (i=0; i<=cfd->NBINARY; i++) {
		j = cfd->bs_index[i];
		binary[j].id1 = cfd->bs_id1[i];
		binary[j].id2 = cfd->bs_id2[i];
		if (binary[j].id1 > newstarid) {
			newstarid = binary[j].id1;
		}
		if (binary[j].id2 > newstarid) {
			newstarid = binary[j].id2;
		}
		binary[j].rad1 = cfd->bs_Reff1[i];
		binary[j].rad2 = cfd->bs_Reff2[i];
		binary[j].m1 = cfd->bs_m1[i] * ((double) clus.N_STAR);
		binary[j].m2 = cfd->bs_m2[i] * ((double) clus.N_STAR);
		binary[j].a = cfd->bs_a[i];
		binary[j].e = cfd->bs_e[i];
		binary[j].inuse = 1;
	}

	/* some assignments so the code won't break */
	star[0].r = ZERO; 
	star[clus.N_STAR + 1].r = SF_INFINITY;
	Mtotal = 1.0;

	/* some normalization that set_imf() used to take care of */
	totmass = clus.N_STAR;
	avemass = 1.0;
	initial_total_mass = totmass; //+ star[0].m;

        /* check if there is a black hole in the file */
        if (OVERRIDE_CENTRAL_MASS) {
		cenma.m = CENTRAL_MASS * MSUN / units.mstar / avemass;
        } else {
		cenma.m = star[0].m / avemass;
        }
}
