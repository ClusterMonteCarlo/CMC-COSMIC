/* -*- linux-c -*- */
/* trivial */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <sys/times.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#define _MAIN_
#include "cmc_vars.h"

#ifdef USE_CUDA
#include "cuda/cmc_cuda.h"
#endif

int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i, j;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	/* set some important global variables */
	quiet = 0;
	debug = 0;
	initial_total_mass = 1.0;
	newstarid = 0;
	cenma.E = 0.0;

	/* trap signals */
	trap_sigs();

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	
	/* parse input */
	parser(argc, argv, rng);

	/* print version information to log file */
	print_version(logfile);

        /* initialize the Search_Grid r_grid */
        if (SEARCH_GRID) {
          r_grid= search_grid_initialize(SG_POWER_LAW_EXPONENT, \
              SG_MATCH_AT_FRACTION, SG_STARSPERBIN, SG_PARTICLE_FRACTION);
        };

        /* initialize the root finder algorithm */
        q_root = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

#ifdef DEBUGGING
        /* create a new hash table for the star ids */
        star_ids= g_hash_table_new(g_int_hash, g_int_equal);
        /* load a trace list */
        //load_id_table(star_ids, "trace_list");
#endif
	/* Set up initial conditions */
	load_fits_file_data();

	TidalMassLoss = 0.0;
	Etidal = 0.0;
	N_bb = 0;			/* number of bin-bin interactions */
	N_bs = 0;
	E_bb = 0.0;
	E_bs = 0.0;
	Echeck = 0; 		
	se_file_counter = 0; 	
	snap_num = 0; 		
	bh_snap_num = 0;
	StepCount = 0; 		
	tcount = 1;
	TotalTime = 0.0;
	reset_rng_t113(IDUM);

	/* Meagan - 3bb */
	N3bbformed = 0;
	delta_E_3bb = 0.0;

	/* Meagan - bh counters */
	bhsingle = 0;
	bhbinary = 0;
	bhbh = 0;
	bhnonbh = 0;
	bh13 = 0;
	bh10 = 0;
	bh11 = 0;
	bh12 = 0;
	bhwd = 0;
	bhstar = 0;
	bh01=0;
	bh26=0;
	bh7=0;
	bh89=0;
	esc_bhsingle = 0;
	esc_bhbinary = 0;
	esc_bhbh = 0;
	esc_bhnonbh = 0;
	esc_bh13 = 0;
	esc_bh10 = 0;
	esc_bh11 = 0;
	esc_bh12 = 0;
	esc_bhwd = 0;
	esc_bhstar = 0;
	esc_bh01 = 0;
	esc_bh26 = 0;
	esc_bh7 = 0;
	esc_bh89 = 0;
	esc_bhsingle_tot = 0;
	esc_bhbinary_tot = 0;
	esc_bhbh_tot = 0;
	esc_bhnonbh_tot = 0;
	esc_bh13_tot = 0;
	esc_bh10_tot = 0;
	esc_bh11_tot = 0;
	esc_bh12_tot = 0;
	esc_bhwd_tot = 0;
	esc_bhstar_tot = 0;
	esc_bh01_tot = 0;
	esc_bh26_tot = 0;
	esc_bh7_tot = 0;
	esc_bh89_tot = 0;

	/* binary remainders */
	clus.N_MAX = clus.N_STAR;
	N_b = clus.N_BINARY;
	calc_sigma_r();
	central_calculate();
        if (WRITE_EXTRA_CORE_INFO) {
          no_remnants= no_remnants_core(6);
        }
	M_b = 0.0;
	E_b = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		j = star[i].binind;
		if (j && binary[j].inuse) {
			M_b += star[i].m;
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}
	
	/* print out binary properties to a file */
	print_initial_binaries();
	
	orbit_r = R_MAX;
	potential_calculate();
	total_bisections= 0;
	if (SEARCH_GRID) 
		search_grid_update(r_grid);
	
	/* compute energy initially */
	star[0].E = star[0].J = 0.0;
	ComputeEnergy();
	star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;

	/* Noting the total initial energy, in order to set termination energy. */
	Etotal.ini = Etotal.tot;
	
	Etotal.New = 0.0;
	Eescaped = 0.0;
	Jescaped = 0.0;
	Eintescaped = 0.0;
	Ebescaped = 0.0;
	Eoops = 0.0;
	
	comp_mass_percent();
	comp_multi_mass_percent();

	/* If we don't set it here, new stars created by breaking binaries (BSE) will
	 * end up in the wrong place */
	clus.N_MAX_NEW = clus.N_MAX;
	/* initialize stellar evolution things */
	DMse = 0.0;
	if (STELLAR_EVOLUTION > 0) {
		stellar_evolution_init();
	}
	
	update_vars();
	
	times(&tmsbufref);

	/* calculate central quantities */
	central_calculate();
        if (WRITE_EXTRA_CORE_INFO) {
          no_remnants= no_remnants_core(6);
        }

	/* calculate dynamical quantities */
	clusdyn_calculate();

	/* Printing Results for initial model */
	print_results();

	/* print handy script for converting output files to physical units */
	print_conversion_script();
	
#ifdef USE_CUDA
	cuInitialize();
#endif

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) {
		/* calculate central quantities */
		central_calculate();
		
		/* Get new time step */
		Dt = GetTimeStep(rng);
		
		/*dprintf ("before tide: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		     star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
	     		star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		/* if tidal mass loss in previous time step is > 5% reduce PREVIOUS timestep by 20% */
		if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
			diaprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
			Dt = Prev_Dt * 0.8;
		} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
			diaprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
			Dt = Prev_Dt * 1.1;
		}
		
		/*dprintf ("after tide before dynamics: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		     star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
	     		star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/
			
		TotalTime += Dt;

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		clus.N_MAX_NEW = clus.N_MAX;

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		if (PERTURB > 0) {
			dynamics_apply(Dt, rng);
		}

		/*dprintf ("after dynamics before SE: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		     star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
	     		star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
		   will disappear! */
		clus.N_MAX_NEW++;

		/* evolve stars up to new time */
		DMse = 0.0;
		if (STELLAR_EVOLUTION > 0) {
			do_stellar_evolution(rng);
		}

		/*dprintf ("after SE: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		     star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
	     		star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/
		
		Prev_Dt = Dt;

		/* some numbers necessary to implement Stodolkiewicz's
	 	 * energy conservation scheme */
		for (i = 1; i <= clus.N_MAX_NEW; i++) {
			/* saving velocities */
			star[i].vtold = star[i].vt;
			star[i].vrold = star[i].vr;
			
			/* the following will get updated after sorting and
			 * calling potential_calculate(), needs to be saved 
			 * now */  
			star[i].Uoldrold = star[i].phi + PHI_S(star[i].r, i);
			
			/* Unewrold will be calculated after 
			 * potential_calculate() using [].rOld
			 * Unewrnew is [].phi after potential_calculate() */
		}

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		DMrejuv = 0.0;
		if (STAR_AGING_SCHEME > 0) {
			for (i=1; i<=clus.N_MAX; i++){
				remove_old_star(TotalTime, i);
			}	
		}
		
		/* this calls get_positions() */
		tidally_strip_stars();

		/* more numbers necessary to implement Stodolkiewicz's
	 	 * energy conservation scheme */
		for (i = 1; i <= clus.N_MAX_NEW; i++) {
			/* the following cannot be calculated after sorting 
			 * and calling potential_calculate() */
			star[i].Uoldrnew = potential(star[i].rnew) + PHI_S(star[i].rnew, i);
		}

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		ComputeIntermediateEnergy();

		/* Sorting stars by radius. The 0th star at radius 0 
		   and (N_STAR+1)th star at SF_INFINITY are already set earlier.
		 */
		qsorts(star+1,clus.N_MAX_NEW);

		potential_calculate();
                if (SEARCH_GRID)
                  search_grid_update(r_grid);

		set_velocities3();

		comp_mass_percent();
		comp_multi_mass_percent();

		ComputeEnergy();

		/* reset interacted flag */
		for (i = 1; i <= clus.N_MAX; i++) {
			star[i].interacted = 0;
			/* Meagan - 3bb */
			star[i].threebb_interacted = 0;

		}
		
		/* update variables, then print */
		update_vars();
		tcount++;
		
		/* calculate dynamical quantities */
		clusdyn_calculate();
                if (WRITE_EXTRA_CORE_INFO) {
                  no_remnants= no_remnants_core(6);
                }

		print_results();
		/* take a snapshot, we need more accurate 
		 * and meaningful criterion 
		 */
		if(tcount%SNAPSHOT_DELTACOUNT==0) {
			print_2Dsnapshot();
			if (WRITE_STELLAR_INFO) {
				write_stellar_data();
			}
		}		
		// Meagan - bh snapshot
		if(tcount%BH_SNAPSHOT_DELTACOUNT==0) {
			print_bh_snapshot();
		}		

	} /* End FOR (time step iteration loop) */

	times(&tmsbuf);
	dprintf("Usr time = %.6e ", (double)
		(tmsbuf.tms_utime-tmsbufref.tms_utime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time = %.6e\n", (double)
		(tmsbuf.tms_stime-tmsbufref.tms_stime)/sysconf(_SC_CLK_TCK));
	dprintf("Usr time (ch) = %.6e ", (double)
		(tmsbuf.tms_cutime-tmsbufref.tms_cutime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time (ch)= %.6e seconds\n", (double)
		(tmsbuf.tms_cstime-tmsbufref.tms_cstime)/sysconf(_SC_CLK_TCK));
	
        printf("The total number of bisections is %li\n", total_bisections);
	/* free RNG */
	gsl_rng_free(rng);
	
#ifdef USE_CUDA
	cuCleanUp();
#endif

	/* flush buffers before returning */
	close_buffers();
	free_arrays();
        if (SEARCH_GRID)
          search_grid_free(r_grid);
#ifdef DEBUGGING
        g_hash_table_destroy(star_ids);
        g_array_free(id_array, TRUE);
#endif
	return(0);
}
