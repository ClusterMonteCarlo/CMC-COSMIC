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

void load_binary_data(int d_k, int s_k) {
		binary[d_k].id1 = cfd.bs_id1[s_k];
		binary[d_k].id2 = cfd.bs_id2[s_k];
		if (binary[d_k].id1 > newstarid) {
			newstarid = binary[d_k].id1;
		}
		if (binary[d_k].id2 > newstarid) {
			newstarid = binary[d_k].id2;
		}
		binary[d_k].rad1 = cfd.bs_Reff1[s_k];
		binary[d_k].rad2 = cfd.bs_Reff2[s_k];
		binary[d_k].m1 = cfd.bs_m1[s_k] * ((double) clus.N_STAR);
		binary[d_k].m2 = cfd.bs_m2[s_k] * ((double) clus.N_STAR);
		binary[d_k].a = cfd.bs_a[s_k];
		binary[d_k].e = cfd.bs_e[s_k];
		binary[d_k].inuse = 1;
		/*Sourav: assign lifetimes to the binary components*/
		if (STAR_AGING_SCHEME==1 ||STAR_AGING_SCHEME==3){
			if (PREAGING){
				binary[d_k].createtime_m1 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
				binary[d_k].createtime_m2 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			}else {
				binary[d_k].createtime_m1 = 0.0;
				binary[d_k].createtime_m2 = 0.0;
			}
			binary[d_k].lifetime_m1 = pow(10.0,9.921)*pow((binary[d_k].m1 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			binary[d_k].lifetime_m2 = pow(10.0,9.921)*pow((binary[d_k].m2 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
		}
		else {
			binary[d_k].createtime_m1 = 0.0;
			binary[d_k].createtime_m2 = 0.0;
			binary[d_k].lifetime_m1 = GSL_POSINF;
			binary[d_k].lifetime_m2 = GSL_POSINF;
		}
}

void load_fits_file_data(void)
{
	long i, g_i=0;
	int b_i=1;

	newstarid = 0;

	/* set units */
	units_set();

	//MPI: Just to load the 0th element (never used) with whatever is in the cfd struct.
	load_binary_data(0, 0);

#ifdef USE_MPI
	//MPI3: Copying only the data each process needs.
	for (i=0; i<=End[myid] - Start[myid]+1; i++) {
		//MPI3: Getting global index to read only stars that belong to this processor.
		g_i = get_global_idx(i);
#else
	/* copy everything over from cfd */
	for (i=0; i<=cfd.NOBJ+1; i++) {
		g_i = i;
#endif
		star[i].id = cfd.obj_id[g_i];
		if (star[i].id > newstarid) {
			newstarid = star[i].id;
		}
		star[i].se_k = cfd.obj_k[g_i];
		star[i].rad = cfd.obj_Reff[g_i];
		star[i].vr = cfd.obj_vr[g_i];
		star[i].vt = cfd.obj_vt[g_i];

      if(cfd.obj_binind[g_i])
      {
#ifdef USE_MPI
         load_binary_data(b_i, cfd.obj_binind[g_i]);
         star[i].binind = b_i;
#else
			load_binary_data(cfd.obj_binind[g_i], cfd.obj_binind[g_i]);
			star[i].binind = cfd.obj_binind[g_i];
#endif
         b_i++;
      }

		star[i].m = cfd.obj_m[g_i] * ((double) clus.N_STAR);
		star[i].r = cfd.obj_r[g_i];
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

#ifdef USE_MPI
	num_bin = b_i-1;
#else
	if((b_i-1)!=cfd.NOBJ) eprintf("Binary number inconsistent: in binary array = %d, in cfd = %d\n", b_i, cfd.NOBJ);
#endif

#ifdef USE_MPI
	//MPI3: All procs read data for global arrays.
	for (i=0; i<=cfd.NOBJ+1; i++) {
		star_m[i] = cfd.obj_m[i] * ((double) clus.N_STAR);
		star_r[i] = cfd.obj_r[i];
	}
#endif

	/* some assignments so the code won't break */
#ifdef USE_MPI
	star_r[0] = ZERO; 
	star_r[clus.N_STAR + 1] = SF_INFINITY;
#else
	star[0].r = ZERO; 
	star[clus.N_STAR + 1].r = SF_INFINITY;
#endif
	Mtotal = 1.0;

	// central mass business, read in from file
	// I believe the normalization should be correct, from above
#ifdef USE_MPI
	cenma.m = star_m[0];
	cenma.m_new= star_m[0];
	star_m[0] = 0.0;
#else
	cenma.m = star[0].m;
	cenma.m_new= star[0].m;
	star[0].m = 0.0;
#endif

	if (BH_R_DISRUPT_NB> 0) {
		diaprintf("R_disrupt in NB-units for all stars, Rdisr=%lg\n", BH_R_DISRUPT_NB);
	} else {
		diaprintf("R_disrupt for a solar mass star in NB-units, Rdisr=%lg\n", 
				pow(2.*cenma.m/MSUN*units.mstar, 1./3.)*RSUN/units.l);
	};
}
