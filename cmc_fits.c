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
	for (i=0; i<=cfd.NOBJ+1; i++) {
		star[i].id = cfd.obj_id[i];
		if (star[i].id > newstarid) {
			newstarid = star[i].id;
		}
		star[i].se_k = cfd.obj_k[i];
		star[i].m = cfd.obj_m[i] * ((double) clus.N_STAR);
		star[i].rad = cfd.obj_Reff[i];
		star[i].r = cfd.obj_r[i];
		star[i].vr = cfd.obj_vr[i];
		star[i].vt = cfd.obj_vt[i];
		star[i].binind = cfd.obj_binind[i];
		/*Sourav: putting creation time and lifetime as a variable*/
		/*Sourav: this is ongoing changes*/
		if (!star[i].binind){
			if (STAR_AGING_SCHEME==1 ||STAR_AGING_SCHEME==3){
				if(PREAGING){
					star[i].createtime = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
				} else {
					star[i].createtime = 0.0;
				}	
				star[i].lifetime = pow(10.0,9.921)*pow((star[i].m * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			}
			else {
				star[i].createtime = 0.0;
				star[i].lifetime = GSL_POSINF;
			}
		}

	}

	for (i=0; i<=cfd.NBINARY; i++) {
		j = cfd.bs_index[i];
		binary[j].id1 = cfd.bs_id1[i];
		binary[j].id2 = cfd.bs_id2[i];
		if (binary[j].id1 > newstarid) {
			newstarid = binary[j].id1;
		}
		if (binary[j].id2 > newstarid) {
			newstarid = binary[j].id2;
		}
		binary[j].rad1 = cfd.bs_Reff1[i];
		binary[j].rad2 = cfd.bs_Reff2[i];
		binary[j].m1 = cfd.bs_m1[i] * ((double) clus.N_STAR);
		binary[j].m2 = cfd.bs_m2[i] * ((double) clus.N_STAR);
		binary[j].a = cfd.bs_a[i];
		binary[j].e = cfd.bs_e[i];
		binary[j].inuse = 1;
		/*Sourav: assign lifetimes to the binary components*/
		if (STAR_AGING_SCHEME==1 ||STAR_AGING_SCHEME==3){
			if (PREAGING){
				binary[j].createtime_m1 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
				binary[j].createtime_m2 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			}else {
				binary[j].createtime_m1 = 0.0;
				binary[j].createtime_m2 = 0.0;
			}
			binary[j].lifetime_m1 = pow(10.0,9.921)*pow((binary[j].m1 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			binary[j].lifetime_m2 = pow(10.0,9.921)*pow((binary[j].m2 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
		}
		else {
			binary[j].createtime_m1 = 0.0;
			binary[j].createtime_m2 = 0.0;
			binary[j].lifetime_m1 = GSL_POSINF;
			binary[j].lifetime_m2 = GSL_POSINF;
		}
	}

	/* some assignments so the code won't break */
	star[0].r = ZERO; 
	star[clus.N_STAR + 1].r = SF_INFINITY;
	Mtotal = 1.0;

	// central mass business, read in from file
	// I believe the normalization should be correct, from above
	cenma.m = star[0].m;
        cenma.m_new= star[0].m;
	star[0].m = 0.0;

#ifdef USE_MPI
	if(myid==0)
#endif
	{
		if (BH_R_DISRUPT_NB> 0) {
			diaprintf("R_disrupt in NB-units for all stars, Rdisr=%lg\n", BH_R_DISRUPT_NB);
		} else {
			diaprintf("R_disrupt for a solar mass star in NB-units, Rdisr=%lg\n", 
					pow(2.*cenma.m/MSUN*units.mstar, 1./3.)*RSUN/units.l);
		};
	}
}
