#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap/bse_wrap.h"

void zero_out_array(double* ptr, int size)
{
	int i;
	for( i=0; i<size; i++)
		ptr[i] = 0.0;
}

void stellar_evolution_init(void){  
	double tphysf, dtp, vs[12];
	int i, j=0;
	long k, kb;
       struct rng_t113_state temp_state;
	binary_t tempbinary;

	/* SSE */
	/* bse_set_hewind(0.5); */

	/* BSE */
	bse_set_neta(BSE_NETA);
	bse_set_bwind(BSE_BWIND);
	bse_set_hewind(BSE_HEWIND);
	bse_set_windflag(BSE_WINDFLAG);
	bse_set_alpha1(BSE_ALPHA1); /* FIXME: is 3 too high? (normally 1.0) */
	bse_set_lambda(BSE_LAMBDA);
	bse_set_ceflag(BSE_CEFLAG);
	bse_set_tflag(BSE_TFLAG);
	bse_set_ifflag(BSE_IFFLAG);
	bse_set_wdflag(BSE_WDFLAG);
	bse_set_bhflag(BSE_BHFLAG);
	bse_set_nsflag(BSE_NSFLAG);
	bse_set_mxns(BSE_MXNS); //3 if nsflag=1 or 2, 1.8 if nsflag=0 (see evolv2.f)
	bse_set_idum(BSE_IDUM);
	bse_set_pts1(0.05);
	bse_set_pts2(0.01);
	bse_set_pts3(0.02);
	bse_set_sigma(BSE_SIGMA);
	bse_set_beta(BSE_BETA);
	bse_set_xi(1.0);
	bse_set_acc2(1.5);
	bse_set_epsnov(0.001);
	bse_set_eddfac(BSE_EDDFAC); /* (normally 1.0) */
	bse_set_gamma(BSE_GAMMA);
  bse_set_merger(-1.0);
  
	/* set parameters relating to metallicity */
	zpars = (double *) malloc(20 * sizeof(double));
	bse_zcnsts(&METALLICITY, zpars);

	/* set collisions matrix */
	bse_instar();
	dprintf("se_init: %g %g %g %d %g %g %d %d %d %d %d %d %g %d %g %g %g %g\n", BSE_NETA, BSE_BWIND, BSE_HEWIND, BSE_WINDFLAG, BSE_ALPHA1, BSE_LAMBDA, BSE_CEFLAG, BSE_TFLAG, BSE_IFFLAG, BSE_WDFLAG, BSE_BHFLAG, BSE_NSFLAG, BSE_MXNS, BSE_IDUM, BSE_SIGMA, BSE_BETA, BSE_EDDFAC, BSE_GAMMA);

#ifdef USE_MPI 
	for (k=1; k<=mpiEnd-mpiBegin+1; k++) {
		j = get_global_idx(k);
#else
	/* set initial properties of stars */
	for (k=1; k<=clus.N_MAX; k++) {
#endif

		if (star[k].binind == 0) { /* single star */
#ifdef USE_MPI
			star[k].se_mass = star_m[j] * units.mstar / MSUN;
#else
			star[k].se_mass = star[k].m * units.mstar / MSUN;
#endif
			/* setting the type */
			if(star[k].se_mass <= 0.7){
				star[k].se_k = 0;
			} else {
				star[k].se_k = 1;
			}
			star[k].se_mt = star[k].se_mass;
			star[k].se_ospin = 0.0;
			star[k].se_epoch = 0.0;
			star[k].se_tphys = 0.0;

			/* evolve slightly (1 year) for initial radii */
			tphysf = 1.0e-6;
			dtp = tphysf - star[k].se_tphys;
			dtp = 0.0;
#ifdef USE_MPI
			DMse += star_m[j] * madhoc;
#else
			DMse_mimic[findProcForIndex(k)] += star[k].m * madhoc;
#endif
			/* Update star id for pass through. */
			bse_set_id1_pass(star[k].id);
			bse_set_id2_pass(0);
			tempbinary.bse_mass0[0] = star[k].se_mass;
			tempbinary.bse_mass0[1] = 0.0;
			tempbinary.bse_kw[0] = star[k].se_k;
			tempbinary.bse_kw[1] = 15;
			tempbinary.bse_mass[0] = star[k].se_mt;
			tempbinary.bse_mass[1] = 0.0;
			tempbinary.bse_radius[0] = star[k].se_radius;
			tempbinary.bse_radius[1] = 0.0;
			tempbinary.bse_lum[0] = star[k].se_lum;
			tempbinary.bse_lum[1] = 0.0;
			tempbinary.bse_massc[0] = star[k].se_mc;
			tempbinary.bse_massc[1] = 0.0;
			tempbinary.bse_radc[0] = star[k].se_rc;
			tempbinary.bse_radc[1] = 0.0;
			tempbinary.bse_menv[0] = star[k].se_menv;
			tempbinary.bse_menv[1] = 0.0;
			tempbinary.bse_renv[0] = star[k].se_renv;
			tempbinary.bse_renv[1] = 0.0;
			tempbinary.bse_ospin[0] = star[k].se_ospin;
			tempbinary.bse_ospin[1] = 0.0;
			/*
				tempbinary.bse_B_0[0] = star[k].se_B_0;
				tempbinary.bse_B_0[1] = 0.0;
				tempbinary.bse_bacc[0] = star[k].se_bacc;
				tempbinary.bse_bacc[1] = 0.0;
				tempbinary.bse_tacc[0] = star[k].se_tacc;
				tempbinary.bse_tacc[1] = 0.0;
			 */
			tempbinary.bse_epoch[0] = star[k].se_epoch;
			tempbinary.bse_epoch[1] = 0.0;
			tempbinary.bse_tms[0] = star[k].se_tms;
			tempbinary.bse_tms[1] = 0.0;
			tempbinary.bse_tb = 0.0;
			tempbinary.e = 0.0;
			/*
				bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
				&(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
				&(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
				&(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
			 */
                        get_rng_t113(&temp_state);
                        bse_set_taus113state(temp_state, 0);
			bse_evolv2(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
					&(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
					&(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
					&(tempbinary.bse_ospin[0]),
					&(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
					&(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
					&(tempbinary.bse_tb), &(tempbinary.e), vs);

                        temp_state=bse_get_taus113state();
                        set_rng_t113(temp_state);
			//MPI2: Since the BSE rng is not yet parallelized, disabling kicks by setting vs[] to zero.
			//zero_out_array(vs, 12);

			star[k].se_mass = tempbinary.bse_mass0[0];
			star[k].se_k = tempbinary.bse_kw[0];
			star[k].se_mt = tempbinary.bse_mass[0];
			star[k].se_radius = tempbinary.bse_radius[0];
			star[k].se_lum = tempbinary.bse_lum[0];
			star[k].se_mc = tempbinary.bse_massc[0];
			star[k].se_rc = tempbinary.bse_radc[0];
			star[k].se_menv = tempbinary.bse_menv[0];
			star[k].se_renv = tempbinary.bse_renv[0];
			star[k].se_ospin = tempbinary.bse_ospin[0];
			/*
				star[k].se_B_0 = tempbinary.bse_B_0[0];
				star[k].se_bacc = tempbinary.bse_bacc[0];
				star[k].se_tacc = tempbinary.bse_tacc[0];
			 */
			star[k].se_epoch = tempbinary.bse_epoch[0];
			star[k].se_tms = tempbinary.bse_tms[0];

			star[k].rad = star[k].se_radius * RSUN / units.l;
#ifdef USE_MPI
			star_m[j] = star[k].se_mt * MSUN / units.mstar;
			DMse -= star_m[j] * madhoc;
#else
			star[k].m = star[k].se_mt * MSUN / units.mstar;
			DMse_mimic[findProcForIndex(k)] -= star[k].m * madhoc;
#endif
			/* birth kicks */
			if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[2]*vs[2]) != 0.0) {
				//dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
			}
			star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);

#ifndef USE_MPI
			curr_st = &st[findProcForIndex(k)];
			//printf("starid = %d proc = %d start = %d end = %d\n", k, findProcForIndex(k), Start[findProcForIndex(k)], End[findProcForIndex(k)]);
#endif

			vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
			//star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
			set_star_EJ(k);
		} else if (star[k].binind > 0) { /* binary */
			star[k].se_k = NOT_A_STAR; /* just for safety */
			kb = star[k].binind;
			binary[kb].bse_mass0[0] = binary[kb].m1 * units.mstar / MSUN;
			binary[kb].bse_mass0[1] = binary[kb].m2 * units.mstar / MSUN;
			for (i=0; i<=1; i++) {
				if(binary[kb].bse_mass0[i] <= 0.7){
					binary[kb].bse_kw[i] = 0;
				} else {
					binary[kb].bse_kw[i] = 1;
				}
				binary[kb].bse_mass[i] = binary[kb].bse_mass0[i];
				binary[kb].bse_ospin[i] = 0.0;
				binary[kb].bse_epoch[i] = 0.0;
			}
			binary[kb].bse_tphys = 0.0;

			/* set binary orbital period (in days) from a */
			binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;

			/* evolve slightly (1 year) for initial radii */
			tphysf = 1.0e-6;
			dtp = tphysf - binary[kb].bse_tphys;
			dtp = 0.0;
#ifdef USE_MPI
			DMse += (binary[kb].m1 + binary[kb].m2) * madhoc;
#else
			DMse_mimic[findProcForIndex(k)] += (binary[kb].m1 + binary[kb].m2) * madhoc;
#endif
			/* Update star id for pass through. */
			bse_set_id1_pass(binary[kb].id1);
			bse_set_id2_pass(binary[kb].id2);
                        get_rng_t113(&temp_state);
                        bse_set_taus113state(temp_state, 0);
			bse_evolv2(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
					&(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
					&(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]), &(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
					&(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
					&(binary[kb].bse_tb), &(binary[kb].e), vs);
                        temp_state=bse_get_taus113state();
                        set_rng_t113(temp_state);
			//MPI2: Since the BSE rng is not yet parallelized, disabling kicks by setting vs[] to zero.
			//zero_out_array(vs, 12);

			handle_bse_outcome(k, kb, vs, tphysf);
		} else {
			eprintf("totally confused!\n");
			exit_cleanly(-1, __FUNCTION__);
		}
	}

#ifdef USE_MPI
	//if(myid==0)
		MPI_Allreduce(MPI_IN_PLACE, &DMse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//else
	//	MPI_Allreduce(&DMse, &DMse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	for(i=0; i<procs; i++)
		DMse += DMse_mimic[i];
#endif
}

/* note that this routine is called after perturb_stars() and get_positions() */
void do_stellar_evolution(gsl_rng *rng)
{
	long k, kb;
	int g_k;
	int kprev;
	double dtp, tphysf, vs[12], VKO;
        struct rng_t113_state temp_state;
	binary_t tempbinary;
  bse_set_merger(-1.0);
	/* double vk, theta; */

	//MPI2: Running till N_MAX_NEW+1 following the bugfix of incrementing N_MAX_NEW after SE.
	//MPI3: Why to run till N_MAX_NEW+1?! I guess it is to account for the sentinel. But now, we have no sentinel.
#ifdef USE_MPI
	for(k=1; k<=clus.N_MAX_NEW; k++){ 
		g_k = get_global_idx(k);
#else
	for(k=1; k<=clus.N_MAX_NEW+1; k++){ 
#endif
		if (star[k].binind == 0) { /* single star */
			tphysf = TotalTime / MEGA_YEAR;
			dtp = tphysf;
			dtp = 0.0;
			kprev = star[k].se_k;	  

#ifdef USE_MPI
			if (star_m[get_global_idx(k)]<=DBL_MIN && star[k].vr==0. && star[k].vt==0. && star[k].E==0. && star[k].J==0.){ //ignoring zeroed out stars
				dprintf ("zeroed out star: skipping SE:\n"); 
				dprintf ("k=%ld m=%g r=%g phi=%g vr=%g vt=%g E=%g J=%g\n", k, star_m[k], star_r[k], star_phi[k], star[k].vr, star[k].vt, star[k].E, star[k].J);	
#else
			if (star[k].m<=DBL_MIN && star[k].vr==0. && star[k].vt==0. && star[k].E==0. && star[k].J==0.){ //ignoring zeroed out stars
					dprintf ("zeroed out star: skipping SE:\n"); 
					dprintf ("k=%ld m=%g r=%g phi=%g vr=%g vt=%g E=%g J=%g\n", k, star[k].m, star[k].r, star[k].phi, star[k].vr, star[k].vt, star[k].E, star[k].J);	
#endif
			} else {
#ifdef USE_MPI
				DMse += star_m[g_k] * madhoc;
#else
				DMse_mimic[findProcForIndex(k)] += star[k].m * madhoc;
#endif
				/* Update star id for pass through. */
				bse_set_id1_pass(star[k].id);
				bse_set_id2_pass(0);
				tempbinary.bse_mass0[0] = star[k].se_mass;
				tempbinary.bse_mass0[1] = 0.0;
				tempbinary.bse_kw[0] = star[k].se_k;
				tempbinary.bse_kw[1] = 15;
				tempbinary.bse_mass[0] = star[k].se_mt;
				tempbinary.bse_mass[1] = 0.0;
				tempbinary.bse_radius[0] = star[k].se_radius;
				tempbinary.bse_radius[1] = 0.0;
				tempbinary.bse_lum[0] = star[k].se_lum;
				tempbinary.bse_lum[1] = 0.0;
				tempbinary.bse_massc[0] = star[k].se_mc;
				tempbinary.bse_massc[1] = 0.0;
				tempbinary.bse_radc[0] = star[k].se_rc;
				tempbinary.bse_radc[1] = 0.0;
				tempbinary.bse_menv[0] = star[k].se_menv;
				tempbinary.bse_menv[1] = 0.0;
				tempbinary.bse_renv[0] = star[k].se_renv;
				tempbinary.bse_renv[1] = 0.0;
				tempbinary.bse_ospin[0] = star[k].se_ospin;
				tempbinary.bse_ospin[1] = 0.0;
				/*
					tempbinary.bse_B_0[0] = star[k].se_B_0;
					tempbinary.bse_B_0[1] = 0.0;
					tempbinary.bse_bacc[0] = star[k].se_bacc;
					tempbinary.bse_bacc[1] = 0.0;
					tempbinary.bse_tacc[0] = star[k].se_tacc;
					tempbinary.bse_tacc[1] = 0.0;
				 */
				tempbinary.bse_epoch[0] = star[k].se_epoch;
				tempbinary.bse_epoch[1] = 0.0;
				tempbinary.bse_tms[0] = star[k].se_tms;
				tempbinary.bse_tms[1] = 0.0;
				tempbinary.bse_tb = 0.0;
				tempbinary.e = 0.0;
				/*
					bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
					&(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
					&(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
					&(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
				 */
                                get_rng_t113(&temp_state);
                                bse_set_taus113state(temp_state, 0);
				bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
						&(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
						&(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
						&(tempbinary.bse_ospin[0]), 
						&(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
						&(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
						&(tempbinary.bse_tb), &(tempbinary.e), vs);
                                temp_state=bse_get_taus113state();
                                set_rng_t113(temp_state);

				//MPI2: Since the BSE rng is not yet parallelized, disabling kicks by setting vs[] to zero.
				//zero_out_array(vs, 12);

				star[k].se_mass = tempbinary.bse_mass0[0];
				star[k].se_k = tempbinary.bse_kw[0];
				star[k].se_mt = tempbinary.bse_mass[0];
				star[k].se_radius = tempbinary.bse_radius[0];
				star[k].se_lum = tempbinary.bse_lum[0];
				star[k].se_mc = tempbinary.bse_massc[0];
				star[k].se_rc = tempbinary.bse_radc[0];
				star[k].se_menv = tempbinary.bse_menv[0];
				star[k].se_renv = tempbinary.bse_renv[0];
				star[k].se_ospin = tempbinary.bse_ospin[0];
				/*
					star[k].se_B_0 = tempbinary.bse_B_0[0];
					star[k].se_bacc = tempbinary.bse_bacc[0];
					star[k].se_tacc = tempbinary.bse_tacc[0];
				 */
				star[k].se_epoch = tempbinary.bse_epoch[0];
				star[k].se_tms = tempbinary.bse_tms[0];

				star[k].rad = star[k].se_radius * RSUN / units.l;
#ifdef USE_MPI
				star_m[g_k] = star[k].se_mt * MSUN / units.mstar;
				DMse -= star_m[g_k] * madhoc;
#else
				star[k].m = star[k].se_mt * MSUN / units.mstar;
				DMse_mimic[findProcForIndex(k)] -= star[k].m * madhoc;
#endif     

				/* birth kicks */
				if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
					//dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
				}
				star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);


#ifndef USE_MPI
				curr_st = &st[findProcForIndex(k)];
#endif

				vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
				//star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
				set_star_EJ(k);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
				/* birth kicks */
				if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
					//   dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
					dprintf("birth kick(iso): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g VKO=%.18g type=%d star_id=%ld Pi=%g\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[k].se_k, star[k].id,star[k].se_ospin);
				}

				/* WD birth kicks, just in case they exist */
				/* if ((star[k].se_k >= 10 && star[k].se_k <= 12) && star[k].se_k != kprev) { */
				/* vk = 2.0e5 / (units.l/units.t); */
				/* 	theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0); */
				/* 	star[k].vr += cos(theta) * vk; */
				/* 	star[k].vt += sin(theta) * vk; */
				/* 	set_star_EJ(k); */
				/* } */
			}
		} else { /* binary */
			tphysf = TotalTime / MEGA_YEAR;
			kb = star[k].binind;
			dtp = tphysf - binary[kb].bse_tphys;
			dtp = 0.0;
#ifdef USE_MPI
			if (star_m[g_k]<=DBL_MIN && binary[kb].a==0. && binary[kb].e==0. && binary[kb].m1==0. && binary[kb].m2==0.){ //ignoring zeroed out binaries
				if(myid==0) {
					dprintf ("zeroed out star: skipping SE:\n");	
					dprintf ("k=%ld kb=%ld m=%g m1=%g m2=%g a=%g e=%g r=%g\n", k, kb, star_m[k], binary[kb].m1, binary[kb].m2, binary[kb].a, binary[kb].e, star_r[k]);
				}
#else
			if (star[k].m<=DBL_MIN && binary[kb].a==0. && binary[kb].e==0. && binary[kb].m1==0. && binary[kb].m2==0.){ //ignoring zeroed out binaries
						dprintf ("zeroed out star: skipping SE:\n");	
						dprintf ("k=%ld kb=%ld m=%g m1=%g m2=%g a=%g e=%g r=%g\n", k, kb, star[k].m, binary[kb].m1, binary[kb].m2, binary[kb].a, binary[kb].e, star[k].r);
#endif
			} else {
				/* set binary orbital period (in days) from a */
				binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;
#ifdef USE_MPI
				DMse += (binary[kb].m1 + binary[kb].m2) * madhoc;
#else
				DMse_mimic[findProcForIndex(k)] += (binary[kb].m1 + binary[kb].m2) * madhoc;
#endif
				/* Update star id for pass through. */
				bse_set_id1_pass(binary[kb].id1);
				bse_set_id2_pass(binary[kb].id2);
                                get_rng_t113(&temp_state);
                                bse_set_taus113state(temp_state, 0);
				bse_evolv2_safely(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
						&(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
						&(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]), &(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
						&(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
						&(binary[kb].bse_tb), &(binary[kb].e), vs);
                                temp_state=bse_get_taus113state();
                                set_rng_t113(temp_state);

				//MPI2: Since the BSE rng is not yet parallelized, disabling kicks by setting vs[] to zero.
				//zero_out_array(vs, 12);

				if(isnan(binary[kb].bse_radius[0])){
					fprintf(stderr, "An isnan occured for r1 cmc_stellar_evolution.c\n");
					fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g\n", binary[kb].bse_tphys, tphysf, binary[kb].bse_kw[0], binary[kb].bse_kw[1], binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_radius[0], binary[kb].bse_radius[1], binary[kb].bse_lum[0], binary[kb].bse_lum[1], binary[kb].bse_tb);
					fprintf(stderr, "k= %ld kb=%ld star_id=%ld bin_id1=%ld bin_id2=%ld \n", k, kb, star[k].id, binary[kb].id1, binary[kb].id2);
					exit(1);
				} 
				if(isnan(binary[kb].bse_radius[1])){
					fprintf(stderr, "An isnan occured for r2 cmc_stellar_evolution.c\n");
					fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", binary[kb].bse_tphys, tphysf, binary[kb].bse_kw[0], binary[kb].bse_kw[1], binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_radius[0], binary[kb].bse_radius[1], binary[kb].bse_lum[0], binary[kb].bse_lum[1]);
					fprintf(stderr, "k= %ld kb=%ld star_id=%ld bin_id1=%ld bin_id2=%ld \n", k, kb, star[k].id, binary[kb].id1, binary[kb].id2);
					exit(1);
				}
				handle_bse_outcome(k, kb, vs, tphysf);
			}
		}
	}

#ifdef USE_MPI
	double temp = 0.0;

	MPI_Status stat;
	int i;
	
	//MPI2: Avoiding reduce to improve accuracy.
	if(myid!=0)
		MPI_Send(&DMse, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	else
		for(i=1;i<procs;i++)
		{
			MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
			DMse += temp;
		}
#else
	int i;
	for(i=0; i<procs; i++)
		DMse += DMse_mimic[i];
#endif
}

void write_stellar_data(void){
  long k, kb;
  FILE *stel_file;
  char filename[1024];
  
  se_file_counter++;
  
  /* single star info */
  sprintf(filename, "%s_stellar_info.%05d.dat", outprefix, se_file_counter);
  stel_file = fopen(filename, "w");
  if (stel_file==NULL){
    fprintf(stderr,
	    "file cannot be opened to write stellar info\n");
    return;
  }
  fprintf(stel_file, "# time (Myr): %e\n",
	  TotalTime/MEGA_YEAR);
  fprintf(stel_file, "# time (FP):  %e\n", TotalTime);
  fprintf(stel_file,
	  "#  id        mass        radius     luminosity  type\n");
  fprintf(stel_file,
	  "#======= ============ ============ ============ ====\n");
  for(k=1; k<=clus.N_MAX; k++){
    fprintf(stel_file, "%08ld ", k);
    fprintf(stel_file, "%e ", star[k].se_mt);
    fprintf(stel_file, "%e ", star[k].se_radius);
    fprintf(stel_file, "%e ", star[k].se_lum);
    fprintf(stel_file, "%d ", star[k].se_k);
    fprintf(stel_file, "\n");
  }
  fclose(stel_file);
  
  /* binary star info */
  sprintf(filename, "%s_binary_stellar_info.%05d.dat", outprefix, se_file_counter);
  stel_file = fopen(filename, "w");
  if (stel_file==NULL){
    fprintf(stderr,
	    "file cannot be opened to write binary stellar info\n");
    return;
  }
  fprintf(stel_file, "# time (Myr): %e\n",
	  TotalTime/MEGA_YEAR);
  fprintf(stel_file, "# time (FP):  %e\n", TotalTime);
  fprintf(stel_file, "#1:id1 #2:id2 #3:M1[MSUN] #4:M2 #5:R1[RSUN] #6:R2 #7:k1 #8:k2 #9:Porb[day] #10:e #11:L1[LSUN] #12:L2 #13:Mcore1[MSUN] #14:Mcore2 #15:Rcore1[RSUN] #16:Rcore2 #17:Menv1[MSUN] #18:Menv2 #19:Renv1[RSUN] #20:Renv2 #21:Tms1[MYR] #22:Tms2 #23:Mdot1[MSUN/YR] #24:Mdot2 #25:R1/ROL1 #26:R2/ROL2\n");
  for(k=1; k<=clus.N_MAX; k++){
    if (star[k].binind) {
      kb = star[k].binind;
      fprintf(stel_file, "%08ld %08ld %g %g %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
	      binary[kb].id1, binary[kb].id2,
	      binary[kb].bse_mass[0], binary[kb].bse_mass[1],
	      binary[kb].bse_radius[0], binary[kb].bse_radius[1],
	      binary[kb].bse_kw[0], binary[kb].bse_kw[1],
	      binary[kb].bse_tb, binary[kb].e,
	      binary[kb].bse_lum[0], binary[kb].bse_lum[1],
	      binary[kb].bse_massc[0], binary[kb].bse_massc[1],
	      binary[kb].bse_radc[0], binary[kb].bse_radc[1],
	      binary[kb].bse_menv[0], binary[kb].bse_menv[1],
	      binary[kb].bse_renv[0], binary[kb].bse_renv[1],
	      binary[kb].bse_tms[0], binary[kb].bse_tms[1],
	      binary[kb].bse_bcm_dmdt[0], binary[kb].bse_bcm_dmdt[1],
	      binary[kb].bse_bcm_radrol[0], binary[kb].bse_bcm_radrol[1]);
    }
  }
  fclose(stel_file);
}

void handle_bse_outcome(long k, long kb, double *vs, double tphysf)
{
  int j, g_k;
  long knew, knewp;
  double dtp, VKO;
  
#ifndef USE_MPI
  curr_st = &st[findProcForIndex(k)];
#endif

  if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0 && binary[kb].bse_tb > 0.0) {
    /* normal evolution */
    binary[kb].rad1 = binary[kb].bse_radius[0] * RSUN / units.l;
    binary[kb].rad2 = binary[kb].bse_radius[1] * RSUN / units.l;
    binary[kb].m1 = binary[kb].bse_mass[0] * MSUN / units.mstar;
    binary[kb].m2 = binary[kb].bse_mass[1] * MSUN / units.mstar;
#ifdef USE_MPI
	 g_k = get_global_idx(k);
    star_m[g_k] = binary[kb].m1 + binary[kb].m2;
    DMse -= star_m[g_k] * madhoc;
#else
    star[k].m = binary[kb].m1 + binary[kb].m2;
    DMse_mimic[findProcForIndex(k)] -= star[k].m * madhoc;
#endif
    binary[kb].a = pow((binary[kb].bse_mass[0]+binary[kb].bse_mass[1])*sqr(binary[kb].bse_tb/365.25), 1.0/3.0)
      * AU / units.l;
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    if (vs[0] <= 0.0) {
        star[k].vr += 0.0;
        star[k].vt += 0.0;
        set_star_EJ(k);
	VKO = 0.0;
    } else if (vs[4] <= 0.0) {
    /* If one kick */
       star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
       //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(k);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else {
    /* Two kicks */
       star[k].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[k].vt),vs[5],vs[6], curr_st);//will mean a different angle for each kick, should be ok as this is a randomised process...
       //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(k);
       VKO = sqrt(vs[1]*vs[1] + vs[2]*vs[2] + vs[3]*vs[3]) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    }
if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
	dprintf("birth kick(bin): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g VK=%.18g id1=%ld id2=%ld pid1=%ld pid2=%ld type1=%d type2=%d \n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,binary[kb].id1,binary[kb].id2,binary[kb].id1,binary[kb].id2, binary[kb].bse_kw[0], binary[kb].bse_kw[1]);
    }
    /* extract some binary info from BSE's bcm array */
    j = 1;
    while (bse_get_bcm(j, 1) >= 0.0) {
      j++;
    }
    j--;
    if (j >= 1) {
      if (fabs((binary[kb].bse_tphys - bse_get_bcm(j,1))/binary[kb].bse_tphys) >= 1.0e-6) {
	wprintf("binary[kb].bse_tphys=%g bcmtime=%g\n", binary[kb].bse_tphys, bse_get_bcm(j,1));
	/* exit_cleanly(-1); */
      }
      binary[kb].bse_bcm_dmdt[0] = bse_get_bcm(j, 14);
      binary[kb].bse_bcm_dmdt[1] = bse_get_bcm(j, 28);
      binary[kb].bse_bcm_radrol[0] = bse_get_bcm(j, 15);
      binary[kb].bse_bcm_radrol[1] = bse_get_bcm(j, 29);
    } else {
      eprintf("Could not extract BSE bcm info!  Input dtp not exactly equal to tphysf-tphys?");
      exit_cleanly(-1, __FUNCTION__);
    }
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0) {
    /* disruption with both stars "intact" */
    //dprintf("binary disrupted via BSE with both stars intact\n");
    knew = create_star(k, 1);
    knewp = create_star(k, 1);
    cp_binmemb_to_star(k, 0, knew);
    cp_binmemb_to_star(k, 1, knewp);
#ifdef USE_MPI
	 //MPI3: There are going to be problems when new stars are created! _m refers to global array, what index should be used?
    DMse -= (star_m[knew] + star_m[knewp]) * madhoc;
#else
    DMse_mimic[findProcForIndex(k)] -= (star[knew].m + star[knewp].m) * madhoc;
#endif

    fprintf(semergedisruptfile, "t=%g disruptboth id1=%ld(m1=%g) id2=%ld(m2=%g) (r=%g)\n", 
	    TotalTime, 
	    star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN,
	    star[knewp].id, star[knewp].m * units.mstar / FB_CONST_MSUN,
	    star[k].r);

    destroy_obj(k);
    /* in this case vs is relative speed between stars at infinity */
/*    star[knew].vr += star[knewp].m/(star[knew].m+star[knewp].m) * vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += star[knewp].m/(star[knew].m+star[knewp].m) * sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    star[knewp].vr += -star[knew].m/(star[knew].m+star[knewp].m) * vs[2] * 1.0e5 / (units.l/units.t);
    star[knewp].vt += -star[knew].m/(star[knew].m+star[knewp].m) * sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knewp);
*/
    /* vs is now the recoil speed of both runaway stars owing to SN disruption */
    if (vs[4]>0.0 && vs[8]<=0.0 ) {
      /* 1 kick occured and it disrupted the system */
        if (vs[0]==1) {
      /* Star knew went SN */
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);

            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            set_star_EJ(knew);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knewp].vr += vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);
            //star[knewp].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t); 
            set_star_EJ(knewp);
        } else {
       /* Star knewp went SN */
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knew].vr += vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            //star[knew].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t); //minus at front?
            set_star_EJ(knew);
        }
    } else if ((vs[4]>0.0 && vs[8]>0.0) && (vs[4] == vs[8])) {
      /* Two SNe and the 2nd one disrupts the system. 
         Thus the primary receives vs[1-3] and vs[9-11] secondary receives vs[1-3] and vs[5-7] */
        if (vs[0]==1) {
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            star[knewp].vt = star[knew].vt; //update, now seperate companion tangent velocity, after first SN when still part of binary here.
            vt_add_kick(&(star[knew].vt),vs[9],vs[10], curr_st);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[9]*vs[9]+vs[10]*vs[10]+vs[11]*vs[11]);
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);//finalise companion vt...

            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        } else {
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            star[knew].vt = star[knewp].vt;
            vt_add_kick(&(star[knewp].vt),vs[9],vs[10], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[9]*vs[9]+vs[10]*vs[10]+vs[11]*vs[11]);
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
        }
    } else if ((vs[4]>0.0 && vs[8]>0.0) && (vs[4] != vs[8])) {
      /* Two SNe and the 1st one disrupts the system.
         Primary feels vs[1-3] and secondary feels vs[5-7] and vs[9-11]. */
        if (vs[0]==1) {
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            set_star_EJ(knew);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knewp].vr += vs[7] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);
            vt_add_kick(&(star[knewp].vt),vs[9],vs[10], curr_st);
            //star[knewp].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        } else {
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
				VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knew].vr += vs[7] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            vt_add_kick(&(star[knew].vt),vs[9],vs[10], curr_st);
            //star[knew].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
        }
    } else {
      /* No kick */
        star[knew].vr += 0.0;
        star[knew].vt += 0.0;
        set_star_EJ(knew);
        VKO = 0.0;
        star[knewp].vr += 0.0;
        star[knewp].vt += 0.0;
        set_star_EJ(knewp);
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
      dprintf("birth kick(iso1.2): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g id1=%ld id2=%ld pid1=%ld pid2=%ld type1=%d type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k,star[knewp].se_k);
    }
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] == 0.0) {
    /* secondary star gone */
    //dprintf("binary disrupted via BSE with first star intact\n");
    knew = create_star(k, 1);
    cp_binmemb_to_star(k, 0, knew);

    fprintf(semergedisruptfile, "t=%g disrupt1 idr=%ld(mr=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g)\n", 
	    TotalTime, 
	    star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN,
	    binary[kb].id1, binary[kb].m1 * units.mstar / FB_CONST_MSUN,
	    binary[kb].id2, binary[kb].m2 * units.mstar / FB_CONST_MSUN,
	    star[k].r);

    destroy_obj(k);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }

/*    star[knew].vr += vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew); */
    /* Update velocity owing to kick(s?) - could be overkill here... */
    if (vs[0]>0.0 && vs[4]<=0.0 ) {
       /* One SN occured */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]==vs[4])) {
       /* 1 SN occured disrupted system then one star killed itself */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]!=vs[4])) {
       /* 2 SNe then system mergers, star feels both kicks (i.e. it is COM) */
       /* NOTE: merger remnants of double compact (NS or BH) binaries do not receive kicks in BSE as of yet!
                When it is introduced may have to put new else if statement in here collecting all (different) vs[]'s. */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    } else {
      /* No kick - going to set_star_EJ seems overkill (here and elsewhere). */
       star[knew].vr += 0.0;
       star[knew].vt += 0.0;
       set_star_EJ(knew);
       VKO = 0.0;
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
	dprintf("birth kick(iso2): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g id1=%ld id2=na pid1=%ld pid2=%ld type1=%d type2=15\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k);
    }
    /* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
    dtp = tphysf - star[knew].se_tphys;
    dtp = 0.0;
  /* Update star id for pass through. */
    bse_set_id1_pass(star[knew].id);
    bse_set_id2_pass(0);
    /* HAVE REMOVED SAFE SINGLE EVOLVE BELIEVING THAT NOW SELF-CONSISTENT MERGERS SHOULD ALWAYS RESULT
    bse_evolv1_safely(&(star[knew].se_k), &(star[knew].se_mass), &(star[knew].se_mt), &(star[knew].se_radius), 
		      &(star[knew].se_lum), &(star[knew].se_mc), &(star[knew].se_rc), &(star[knew].se_menv), 
		      &(star[knew].se_renv), &(star[knew].se_ospin), &(star[knew].se_epoch), &(star[knew].se_tms), 
		      &(star[knew].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
    */
    
    star[knew].rad = star[knew].se_radius * RSUN / units.l;
#ifdef USE_MPI
    star_m[knew] = star[knew].se_mt * MSUN / units.mstar;
    DMse -= star_m[knew] * madhoc;
#else
    star[knew].m = star[knew].se_mt * MSUN / units.mstar;
    DMse_mimic[findProcForIndex(k)] -= star[knew].m * madhoc;
#endif

    /* birth kicks */
    /* ALSO REMOVE BELOW AS WE NOW WONT PRODUCE ANOTHER KICK
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
    vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
    //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
	dprintf("birth kick(iso3): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VKO=%.18g, id=%ld id1=%ld id2=%ld type1=%d type2=15\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k);
    }
    */
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] != 0.0) {
    /* primary star gone */
    //dprintf("binary disrupted via BSE with second star intact\n");
    knew = create_star(k, 1);
    cp_binmemb_to_star(k, 1, knew);

    fprintf(semergedisruptfile, "t=%g disrupt2 idr=%ld(mr=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g)\n", 
	    TotalTime, 
	    star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN,
	    binary[kb].id1, binary[kb].m1 * units.mstar / FB_CONST_MSUN,
	    binary[kb].id2, binary[kb].m2 * units.mstar / FB_CONST_MSUN,
	    star[k].r);
    
    destroy_obj(k);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
/*    star[knew].vr += vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew); */
    /* Update velocity owing to kick(s?) - could be overkill here... */
    if (vs[0]>0.0 && vs[4]<=0.0 ) {
       /* One SN occured */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]==vs[4])) {
       /* 1 SN occured disrupted system then one star killed itself */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]!=vs[4])) {
       /* 2 SNe then system mergers, star feels both kicks (i.e. it is COM) */
       /* NOTE: merger remnants of double compact (NS or BH) binaries do not receive kicks in BSE as of yet!
                When it is introduced may have to put new else if statement in here collecting all (different) vs[]'s. */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    } else {
      /* No kick - going to set_star_EJ seems overkill (here and elsewhere). */
       star[knew].vr += 0.0;
       star[knew].vt += 0.0;
       set_star_EJ(knew);
       VKO = 0.0;
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
	dprintf("birth kick(iso4): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VKO=%.18g id=%ld id1=%ld id2=%ld type1=15 type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO, star[knew].id, binary[kb].id1, binary[kb].id2, star[knew].se_k);
    }
    
    /* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
    dtp = tphysf - star[knew].se_tphys;
    dtp = 0.0;
  /* Update star id for pass through. */
    /* REMOVED AS PART OF EVOLV1 PHASE OUT
    bse_set_id1_pass(star[knew].id);
    bse_set_id2_pass(0);
    bse_evolv1_safely(&(star[knew].se_k), &(star[knew].se_mass), &(star[knew].se_mt), &(star[knew].se_radius), 
		      &(star[knew].se_lum), &(star[knew].se_mc), &(star[knew].se_rc), &(star[knew].se_menv), 
		      &(star[knew].se_renv), &(star[knew].se_ospin), &(star[knew].se_epoch), &(star[knew].se_tms), 
		      &(star[knew].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
    */
    
    star[knew].rad = star[knew].se_radius * RSUN / units.l;
#ifdef USE_MPI
    star_m[knew] = star[knew].se_mt * MSUN / units.mstar;
    DMse -= star_m[knew] * madhoc; 
#else
    star[knew].m = star[knew].se_mt * MSUN / units.mstar;
    DMse_mimic[findProcForIndex(k)] -= star[knew].m * madhoc; 
#endif

    /* birth kicks */
    /* REMOVED AGAIN NOT TO ADD KICK AGAIN
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
    vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
    //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
	dprintf("birth kick(iso6): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g type1=15 type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].se_k);
    }
    */
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] == 0.0) {
    /* both stars gone */
    //dprintf("binary disrupted via BSE with no stars intact\n");
    destroy_obj(k);
  } else {
    dprintf("unhandled binary outcome!\n");
    dprintf("bse_mass0=%g bse_mass1=%g tb=%g\n", 
	    binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_tb);
    exit_cleanly(-1, __FUNCTION__);
  }
}

/* k=star index, kbi=0 or 1, knew=index of new star */
void cp_binmemb_to_star(long k, int kbi, long knew)
{
  long kb;
  
  kb = star[k].binind;
#ifdef USE_MPI
  int g_k = get_global_idx(k);
  star_r[knew] = star_r[g_k];
  star_m[knew] = binary[kb].bse_mass[kbi] * MSUN / units.mstar;
  star_phi[knew] = star_phi[g_k];
#else
  int g_k = k;
  /* and set the stars' dynamical properties */
  star[knew].r = star[k].r;
  star[knew].m = binary[kb].bse_mass[kbi] * MSUN / units.mstar;
  star[knew].phi = star[k].phi;
#endif
  star[knew].vr = star[k].vr;
  star[knew].vt = star[k].vt;
  set_star_EJ(knew);
  set_star_news(knew);
  set_star_olds(knew);
  /* mark stars as interacted so they don't undergo E_CONS mode stuff */
  star[knew].interacted = 1;
  if (kbi == 0) {
    star[knew].Eint = binary[kb].Eint1;
    star[knew].id = binary[kb].id1;
  } else {
    star[knew].Eint = binary[kb].Eint2;
    star[knew].id = binary[kb].id2;
  }
  star[knew].rad = binary[kb].bse_radius[kbi] * RSUN / units.l;
  star[knew].se_mass = binary[kb].bse_mass0[kbi]; /* initial mass (at curent epoch?) */
  star[knew].se_k = binary[kb].bse_kw[kbi];
  star[knew].se_mt = binary[kb].bse_mass[kbi]; /* current mass */
  star[knew].se_ospin = binary[kb].bse_ospin[kbi];
  star[knew].se_epoch = binary[kb].bse_epoch[kbi];
  star[knew].se_tphys = binary[kb].bse_tphys;
  star[knew].se_radius = binary[kb].bse_radius[kbi];
  star[knew].se_lum = binary[kb].bse_lum[kbi];
  star[knew].se_mc = binary[kb].bse_massc[kbi];
  star[knew].se_rc = binary[kb].bse_radc[kbi];
  star[knew].se_menv = binary[kb].bse_menv[kbi];
  star[knew].se_renv = binary[kb].bse_renv[kbi];
  star[knew].se_tms = binary[kb].bse_tms[kbi];
  //Sourav: toy rejuvenation- variables updating
  if (kbi==0){
    star[knew].createtime = binary[kb].createtime_m1;
    star[knew].lifetime = binary[kb].lifetime_m1;
  } else {
    star[knew].createtime = binary[kb].createtime_m2;
    star[knew].lifetime = binary[kb].lifetime_m2;
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
void cp_SEvars_to_newstar(long oldk, int kbi, long knew)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    star[knew].se_mass = star[oldk].se_mass;
    star[knew].se_k = star[oldk].se_k;
    star[knew].se_mt = star[oldk].se_mt;
    star[knew].se_ospin = star[oldk].se_ospin;
    star[knew].se_epoch = star[oldk].se_epoch;
    star[knew].se_tphys = star[oldk].se_tphys;
    star[knew].se_radius = star[oldk].se_radius;
    star[knew].se_lum = star[oldk].se_lum;
    star[knew].se_mc = star[oldk].se_mc;
    star[knew].se_rc = star[oldk].se_rc;
    star[knew].se_menv = star[oldk].se_menv;
    star[knew].se_renv = star[oldk].se_renv;
    star[knew].se_tms = star[oldk].se_tms;
    //Sourav: toy rejuvenation- updating the createtime and lifetime
    star[knew].createtime = star[oldk].createtime;
    star[knew].lifetime = star[oldk].lifetime;
    star[knew].rad = star[oldk].rad;
  } else { /* star comes from input binary */
    star[knew].se_mass = binary[kb].bse_mass0[kbi];
    star[knew].se_k = binary[kb].bse_kw[kbi];
    star[knew].se_mt = binary[kb].bse_mass[kbi];
    star[knew].se_ospin = binary[kb].bse_ospin[kbi];
    star[knew].se_epoch = binary[kb].bse_epoch[kbi];
    star[knew].se_tphys = binary[kb].bse_tphys;
    star[knew].se_radius = binary[kb].bse_radius[kbi];
    star[knew].se_lum = binary[kb].bse_lum[kbi];
    star[knew].se_mc = binary[kb].bse_massc[kbi];
    star[knew].se_rc = binary[kb].bse_radc[kbi];
    star[knew].se_menv = binary[kb].bse_menv[kbi];
    star[knew].se_renv = binary[kb].bse_renv[kbi];
    star[knew].se_tms = binary[kb].bse_tms[kbi];
    //Sourav: toy rejuvenation- updating the rejuv variables for two cases- mass1 and mass2 
    if (kbi==0){
    	star[knew].createtime = binary[kb].createtime_m1;
    	star[knew].lifetime = binary[kb].lifetime_m1;
        star[knew].rad = binary[kb].rad1;
    } else {
	star[knew].createtime = binary[kb].createtime_m2;
    	star[knew].lifetime = binary[kb].lifetime_m2;	
        star[knew].rad = binary[kb].rad2;
    } 
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
void cp_m_to_newstar(long oldk, int kbi, long knew)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
#ifdef USE_MPI
    star_m[get_global_idx(knew)] = star_m[get_global_idx(oldk)];
#else
    star[knew].m = star[oldk].m;
#endif
  } else { /* star comes from input binary */
    if (kbi == 0) {
#ifdef USE_MPI
      star_m[get_global_idx(knew)] = binary[kb].m1; //should this be multiplied by MSUN/units.mstar ?
    } else {
      star_m[get_global_idx(knew)] = binary[kb].m2;//should this be multiplied by MSUN/units.mstar ?
#else
      star[knew].m = binary[kb].m1; //should this be multiplied by MSUN/units.mstar ?
    } else {
      star[knew].m = binary[kb].m2;//should this be multiplied by MSUN/units.mstar ?
#endif
    }
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
void cp_SEvars_to_star(long oldk, int kbi, star_t *target_star)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    target_star->rad = star[oldk].rad; //PDK addition
    target_star->se_mass = star[oldk].se_mass;
    target_star->se_k = star[oldk].se_k;
    target_star->se_mt = star[oldk].se_mt;
    target_star->se_ospin = star[oldk].se_ospin;
    target_star->se_epoch = star[oldk].se_epoch;
    target_star->se_tphys = star[oldk].se_tphys;
    target_star->se_radius = star[oldk].se_radius;
    target_star->se_lum = star[oldk].se_lum;
    target_star->se_mc = star[oldk].se_mc;
    target_star->se_rc = star[oldk].se_rc;
    target_star->se_menv = star[oldk].se_menv;
    target_star->se_renv = star[oldk].se_renv;
    target_star->se_tms = star[oldk].se_tms;
    //Sourav: toy rejuvenation- updating rejuv variables
    target_star->createtime = star[oldk].createtime;
    target_star->lifetime = star[oldk].lifetime;
  } else { /* star comes from input binary */
    target_star->se_mass = binary[kb].bse_mass0[kbi];
    target_star->se_k = binary[kb].bse_kw[kbi];
    target_star->se_mt = binary[kb].bse_mass[kbi];
    target_star->se_ospin = binary[kb].bse_ospin[kbi];
    target_star->se_epoch = binary[kb].bse_epoch[kbi];
    target_star->se_tphys = binary[kb].bse_tphys;
    target_star->se_radius = binary[kb].bse_radius[kbi];
    target_star->se_lum = binary[kb].bse_lum[kbi];
    target_star->se_mc = binary[kb].bse_massc[kbi];
    target_star->se_rc = binary[kb].bse_radc[kbi];
    target_star->se_menv = binary[kb].bse_menv[kbi];
    target_star->se_renv = binary[kb].bse_renv[kbi];
    target_star->se_tms = binary[kb].bse_tms[kbi];
    //Sourav: toy rejuvenation- updating rejuv variables for two cases mass1 and mass2
    if (kbi==1){
   	target_star->createtime = binary[kb].createtime_m1;
	target_star->lifetime = binary[kb].lifetime_m1;
        target_star->rad = binary[kb].rad1; // PDK addition
    } 
    else {
	target_star->createtime = binary[kb].createtime_m2;
	target_star->lifetime = binary[kb].lifetime_m2;
        target_star->rad = binary[kb].rad2; // PDK addition
    }
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
void cp_m_to_star(long oldk, int kbi, star_t *target_star)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
#ifdef USE_MPI
    target_star->m = star_m[get_global_idx(oldk)];
#else
    target_star->m = star[oldk].m;
#endif
  } else { /* star comes from input binary */
    if (kbi == 0) {
      target_star->m = binary[kb].m1;//should this be multiplied by MSUN/units.mstar ?
    } else {
      target_star->m = binary[kb].m2;//should this be multiplied by MSUN/units.mstar ?
    }
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
/* set everything except tb */
void cp_SEvars_to_newbinary(long oldk, int oldkbi, long knew, int kbinew)
{
	long kbold, kbnew;

	kbold = star[oldk].binind;
	kbnew = star[knew].binind;

	if (oldkbi == -1) { /* star comes from input single star */
		binary[kbnew].bse_mass0[kbinew] = star[oldk].se_mass;
		binary[kbnew].bse_kw[kbinew] = star[oldk].se_k;
		binary[kbnew].bse_mass[kbinew] = star[oldk].se_mt;
		binary[kbnew].bse_ospin[kbinew] = star[oldk].se_ospin;
		binary[kbnew].bse_epoch[kbinew] = star[oldk].se_epoch;
		binary[kbnew].bse_tphys = star[oldk].se_tphys; /* tphys should be the same for both input stars so this should be OK */
		binary[kbnew].bse_radius[kbinew] = star[oldk].se_radius;
		binary[kbnew].bse_lum[kbinew] = star[oldk].se_lum;
		binary[kbnew].bse_massc[kbinew] = star[oldk].se_mc;
		binary[kbnew].bse_radc[kbinew] = star[oldk].se_rc;
		binary[kbnew].bse_menv[kbinew] = star[oldk].se_menv;
		binary[kbnew].bse_renv[kbinew] = star[oldk].se_renv;
		binary[kbnew].bse_tms[kbinew] = star[oldk].se_tms;
		//Sourav: toy rejuv- updating rejuv variables to the binary member from the single star
		if (kbinew==0){
			binary[kbnew].createtime_m1 = star[oldk].createtime;
			binary[kbnew].lifetime_m1 = star[oldk].lifetime;
			binary[kbnew].rad1 = star[oldk].rad; // PDK addition
		} else {
			binary[kbnew].createtime_m2 = star[oldk].createtime;
			binary[kbnew].lifetime_m2 = star[oldk].lifetime;
			binary[kbnew].rad2 = star[oldk].rad; // PDK addition
		}
	} else { /* star comes from input binary */
		binary[kbnew].bse_mass0[kbinew] = binary[kbold].bse_mass0[oldkbi];
		binary[kbnew].bse_kw[kbinew] = binary[kbold].bse_kw[oldkbi];
		binary[kbnew].bse_mass[kbinew] = binary[kbold].bse_mass[oldkbi];
		binary[kbnew].bse_ospin[kbinew] = binary[kbold].bse_ospin[oldkbi];
		binary[kbnew].bse_epoch[kbinew] = binary[kbold].bse_epoch[oldkbi];
		binary[kbnew].bse_tphys = binary[kbold].bse_tphys;
		binary[kbnew].bse_radius[kbinew] = binary[kbold].bse_radius[oldkbi];
		binary[kbnew].bse_lum[kbinew] = binary[kbold].bse_lum[oldkbi];
		binary[kbnew].bse_massc[kbinew] = binary[kbold].bse_massc[oldkbi];
		binary[kbnew].bse_radc[kbinew] = binary[kbold].bse_radc[oldkbi];
		binary[kbnew].bse_menv[kbinew] = binary[kbold].bse_menv[oldkbi];
		binary[kbnew].bse_renv[kbinew] = binary[kbold].bse_renv[oldkbi];
		binary[kbnew].bse_tms[kbinew] = binary[kbold].bse_tms[oldkbi];
		//Sourav: toy rejuv- updating rejuv variables to binary members from binary members
		//There can be four cases. m1,2(old)->m1,2(new) 
		if(kbinew==0){
			if (oldkbi==0){
				binary[kbnew].createtime_m1 = binary[kbold].createtime_m1;
				binary[kbnew].lifetime_m1 = binary[kbold].lifetime_m1;
				binary[kbnew].rad1 = binary[kbold].rad1; // PDK addition
			} else {
				binary[kbnew].createtime_m1 = binary[kbold].createtime_m2;
				binary[kbnew].lifetime_m1 = binary[kbold].lifetime_m2;
				binary[kbnew].rad1 = binary[kbold].rad2; // PDK addition
			}
		} else {
			if (oldkbi==0){
				binary[kbnew].createtime_m2 = binary[kbold].createtime_m1;
				binary[kbnew].lifetime_m2 = binary[kbold].lifetime_m1;
				binary[kbnew].rad2 = binary[kbold].rad1; // PDK addition
			} else {
				binary[kbnew].createtime_m2 = binary[kbold].createtime_m2;
				binary[kbnew].lifetime_m2 = binary[kbold].lifetime_m2;
				binary[kbnew].rad2 = binary[kbold].rad2; // PDK addition
			}
		}
	}
}


/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
/* set everything except tb */
void cp_starSEvars_to_binmember(star_t instar, long binindex, int bid)
{
  binary[binindex].bse_mass0[bid] = instar.se_mass;
  binary[binindex].bse_kw[bid] = instar.se_k;
  binary[binindex].bse_mass[bid] = instar.se_mt;
  binary[binindex].bse_ospin[bid] = instar.se_ospin;
  binary[binindex].bse_epoch[bid] = instar.se_epoch;
  binary[binindex].bse_tphys = instar.se_tphys; /* tphys should be the same for both input stars so this should be OK */
  binary[binindex].bse_radius[bid] = instar.se_radius;
  binary[binindex].bse_lum[bid] = instar.se_lum;
  binary[binindex].bse_massc[bid] = instar.se_mc;
  binary[binindex].bse_radc[bid] = instar.se_rc;
  binary[binindex].bse_menv[bid] = instar.se_menv;
  binary[binindex].bse_renv[bid] = instar.se_renv;
  binary[binindex].bse_tms[bid] = instar.se_tms;
  //Sourav: toy rejuv- updating rejuv variables from a single to a binary member
  if (bid==0){
    binary[binindex].createtime_m1 = instar.createtime;
    binary[binindex].lifetime_m1 = instar.lifetime;
    binary[binindex].rad1 = instar.rad; // PDK addition
  } else {
    binary[binindex].createtime_m2 = instar.createtime;
    binary[binindex].lifetime_m2 = instar.lifetime;
    binary[binindex].rad1 = instar.rad; // PDK addition
  }
}

void cp_starmass_to_binmember(star_t instar, long binindex, int bid)
{
  if (bid == 0) {
    binary[binindex].m1 = instar.m;
  } else {
    binary[binindex].m2 = instar.m;
  }
}
