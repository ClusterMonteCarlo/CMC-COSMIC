/* -*- linux-c -*- */

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
#include "cmc_vars.h"

int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	/* set debugging before toggle_debugging can ever be called */
	debug = 0;

	/* Catch some signals */
	signal(SIGINT, exit_cleanly);
	signal(SIGTERM, exit_cleanly);
	signal(SIGQUIT, exit_cleanly);
	signal(SIGUSR1, toggle_debugging);
	
	/* override GSL error handler */
	gsl_set_error_handler(&sf_gsl_errhandler);
	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	
	ReadSnapshot = 0;	/* this becomes 1 if a snapshot is used for restart */

	/* set other variables */
	newstarid = 0;

	/* parses input file, allocates memory, initializes some variables, 
	 * and readies file I/O */
	if (parser(argc, argv, rng) == 0) {
		return 0;
	}

	/* Set up initial conditions, possibly from a previous snapshot */
	/* ReadSnapshot is supposed to be set here if necessary         */
	read_fits_file_data(INPUT_FILE);

	if(!ReadSnapshot){
		/* if the file read is snapshot the following is skipped 
		 * since reading a snapshots sets appropiate values already */
		TidalMassLoss = 0.0;
		Etidal = 0.0;
		N_bb = 0;			/* number of bin-bin interactions */
		N_bs = 0;
		E_bb = 0.0;
		E_bs = 0.0;
		Echeck = 0; 		
		se_file_counter = 0; 	
		snap_num = 0; 		
		StepCount = 0; 		
		tcount = 1;
		TotalTime = 0.0;
		reset_rng_t113(IDUM);
		/* for the sub zone business, pericenter of stars are used to 
		 * determine if they will be relaxed, initially all pericenters 
		 * are set to 0, so all stars will be relaxed */
		for(i=0; i<=clus.N_STAR+1; i++){
			star[i].r_peri = 0.0;
		}
		
		/* assign star id's */
		for(i=1; i<=clus.N_STAR_NEW; i++) {
			star[i].id = star_get_id_new();
		}
		
		/* assign binaries */
		assign_binaries();
		/* assign_binaries_test(); */

		orbit_r = R_MAX;
		potential_calculate();

		ComputeEnergy();
		Etotal.ini = Etotal.tot;   /* Noting the total initial energy, in order to set termination energy. */

		Etotal.New = 0.0;
		Eescaped = 0.0;
		Jescaped = 0.0;
		Eintescaped = 0.0;
		Ebescaped = 0.0;
		Eoops = 0.0;
		
		comp_mass_percent();
		comp_multi_mass_percent();
		
		/* initialize stellar evolution things */
#ifdef SE
		if (STELLAR_EVOLUTION > 0) {
			dprintf("Initialization for stellar evolution...\n");
			stellar_evolution_init();
			dprintf("\tDone.\n");
		}
#endif
		
		sub.count = 0;
		update_vars();
	} else {
		ComputeEnergy2();
	}
	
	times(&tmsbufref);

	/* calculate central quantities */
	central_calculate();

	/* Printing Results for initial model */
	print_results();
	
	/* take an initial snapshot, just for fun */
//	if (CHECKPOINT_PERIOD) {
//		chkpnt_fits(rng);
//	}
	
	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (tcount<=500000) {
#ifdef TRY1
		static ch1=0, ch2=0;

		if(ch2==0 && cenma.m>0){
			ch2 = 1;
		}
#endif
		/* DEBUG
		for (i=0; i<N_STAR_DIM; i++) {
			fprintf(stdout, "i=%ld m=%g rad=%g binind=%ld",
				i, star[i].m * units.mstar / MSUN, star[i].rad * units.l / RSUN, star[i].binind);
			if (star[i].binind) {
				fprintf(stdout, " m1=%g m2=%g rad1=%g rad2=%g",
					binary[star[i].binind].m1 * units.mstar / MSUN, binary[star[i].binind].m2 * units.mstar / MSUN,
					binary[star[i].binind].rad1 * units.l / RSUN, binary[star[i].binind].rad2 * units.l / RSUN);
			}
			fprintf(stdout, "\n");
		}
		exit_cleanly(1);
		DEBUG */

		/*write stellar data from time to time */
		/*if (((tcount-1) % 20 == 0) && (STELLAR_EVOLUTION > 0)){
			write_stellar_data();
		}*/

		/* Check for END of simulation */
		if (CheckStop(tmsbufref) != 0)
			break;

		/* calculate central quantities */
		central_calculate();
		
		/* Get new time step */
		Dt = GetTimeStep(rng);

		/* Make sure stellar mass loss is not greater than 
		 * 1% in the timestep */
		/* FIXME add stuff here FIXME */

		/* if tidal mass loss in previous time step is > 5% reduce 
		   PREVIOUS timestep by 20% */
		if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
			diaprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
			Dt = Prev_Dt * 0.8;
		} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
			diaprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
			Dt = Prev_Dt * 1.1;
		}

		TotalTime = TotalTime + Dt;

		setup_sub_time_step();

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		clus.N_MAX_NEW = clus.N_MAX;

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		if (PERTURB > 0) {
			dynamics_apply(Dt, rng);
		}

		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
		   will disappear! */
		clus.N_MAX_NEW++;

#ifdef SE
		if (STELLAR_EVOLUTION > 0) {
			do_stellar_evolution();
		}
#endif

		Prev_Dt = Dt;

		/* some numbers necessary to implement Stodolkiewicz's
	 	 * energy conservation scheme */
		if(E_CONS==2 || E_CONS==3){
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
		}

		sniff_stars();

		/* more numbers necessary to implement Stodolkiewicz's
	 	 * energy conservation scheme */
		if(E_CONS==2 || E_CONS==3){
			for (i = 1; i <= clus.N_MAX_NEW; i++) {
				/* the following cannot be calculated after sorting 
				 * and calling potential_calculate() */
				star[i].Uoldrnew = potential(star[i].rnew) + PHI_S(star[i].rnew, i);
			}
		}

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		ComputeIntermediateEnergy();

		/* Sorting stars by radius. The 0th star at radius 0 
		   and (N_STAR+1)th star at SF_INFINITY are already set earlier.
		 */
		//mqsort(star, 1, clus.N_MAX_NEW);
		qsorts(star+1,clus.N_MAX_NEW);

		potential_calculate();

		if(E_CONS==2){
			set_velocities();
		} else if (E_CONS == 3) {
			set_velocities3();
		}

		comp_mass_percent();
		comp_multi_mass_percent();

		/* Recompute Energy. Uses the specified ECONS_MODE */
		RecomputeEnergy();

		/* reset interacted flag */
		for (i = 1; i <= clus.N_MAX; i++) {
			star[i].interacted = 0;
		}
		
		/* update variables, then print */
		update_vars();
		tcount++;
		
		print_results();
		/* take a snapshot, we need more accurate 
		 * and meaningful criterion 
		 */
		if(CHECKPOINT_PERIOD && (tcount%CHECKPOINT_PERIOD==0)) {
			chkpnt_fits(rng);
		}
		if(SNAPSHOT_PERIOD && (tcount%SNAPSHOT_PERIOD==0)) {
			print_2Dsnapshot();
		}
	} /* End FOR (time step iteration loop) */

	times(&tmsbuf);
	fprintf(stderr, "Usr time = %.6e ", (double)
	(tmsbuf.tms_utime-tmsbufref.tms_utime)/sysconf(_SC_CLK_TCK));
	fprintf(stderr, "Sys time = %.6e\n", (double)
	(tmsbuf.tms_stime-tmsbufref.tms_stime)/sysconf(_SC_CLK_TCK));
	fprintf(stderr, "Usr time (ch) = %.6e ", (double)
	(tmsbuf.tms_cutime-tmsbufref.tms_cutime)/sysconf(_SC_CLK_TCK));
	fprintf(stderr, "Sys time (ch)= %.6e seconds\n", (double)
	(tmsbuf.tms_cstime-tmsbufref.tms_cstime)/sysconf(_SC_CLK_TCK));

	/* free RNG */
	gsl_rng_free(rng);
	
	/* flush buffers before returning */
	close_buffers();
	free_arrays();
	return(0);
}
