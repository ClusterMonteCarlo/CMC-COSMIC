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

void take_snapshot(void);

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
		for(i=0; i<=clus.N_STAR_NEW; i++) {
			star[i].id = i;
		}
		for (i=clus.N_STAR_NEW+1; i<N_STAR_DIM; i++) {
			star[i].id = -1;
		}
		
		/* assign binaries */
		assign_binaries();

		/* Computing the gravity at the star locations provided in sr[].
		   Results returned in star[].gravity. Returns N_MAX. Also computes
		   Rtidal using Mtotal and orbit_r. */
		orbit_r = R_MAX;
		errstat = gravity();
		Etotal.ini = Etotal.tot;   /* Noting the total initial energy, 
							in order to set termination energy. */
		ComputeEnergy();
		Etotal.New = 0.0; Eescaped = 0.0; Ebescaped = 0.0; Jescaped = 0.0;
		
		comp_mass_percent();
		
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

	/* Printing Results for initial model */
	print_results();
	
	/* FIXME: FITS snapshotting not working yet? */
	/* 
	if (DUMPS) {
		sshot_fits(rng);
	}
	*/
	
	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (tcount<=500000) {
#ifdef TRY1
		static ch1=0, ch2=0;

		if(ch2==0 && cenma.m>0){
			DT_FACTOR /=10.0;
			ch2 = 1;
		}
#endif
		
		
		/*write stellar data from time to time */
		/*if (((tcount-1) % 20 == 0) && (STELLAR_EVOLUTION > 0)){
			write_stellar_data();
		}*/

		/* Check for END of simulation */
		if (CheckStop() != 0)
			break;

		/* Get new time step */
		Dt = GetTimeStep();

		/* Make sure stellar mass loss is not greater than 
		 * 1% in the timestep */
		/* FIXME add stuff here FIXME */

		/* if tidal mass loss in previous time step is > 5% reduce 
		   PREVIOUS timestep by 20% */
		if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
			dprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
			Dt = Prev_Dt * 0.8;
		} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
			dprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
			Dt = Prev_Dt * 1.1;
		}

		Dt /= DT_FACTOR;
		TotalTime = TotalTime + Dt;

		setup_sub_time_step();
		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		if (PERTURB > 0) {
			if (ORIGINAL_PERTURB_STARS) {
				/* perturb_stars(Dt); */
				perturb_stars_new(Dt, rng);
			} else {
				/* perturb_stars_fewbody(Dt, rng); */
				perturb_stars_new(Dt, rng);
			}
		}

#ifdef SE
		if (STELLAR_EVOLUTION > 0) {
			do_stellar_evolution();
		}
#endif

		Prev_Dt = Dt;

		sniff_stars();

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

		/* Computing the gravity at the star locations provided in sr[].
		   Results returned in star[].gravity. Computes Mtotal 
		   (total mass).  Returns N_MAX. 
		   This is the  only place where N_MAX can change.
		   Also computes new Rtidal using Mtotal and orbit_r */
		errstat = gravity();

		comp_mass_percent();

		/* Recompute Energy. Uses the specified ECONS_MODE */
		RecomputeEnergy();
		
		/* update variables, then print */
		update_vars();
		tcount++;
		print_results();
		/* FIXME: FITS snapshotting not working yet? */
		/* if(tcount%1==0 && DUMPS) sshot_fits(rng); */
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
