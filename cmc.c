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
#include <string.h>

#ifdef USE_CUDA
#include "cuda/cmc_cuda.h"
#endif

#define PROC 4 //to mimic rng of the serial version to the parallel version. This has to be changed every time the no.of processors you want to simulate is changed.

int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i, j;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	//MPI2: Time variables to compute total time.
	double timeS, timeE;
	//Temp file handle for debugging
	FILE *ftest;
	char num[5],filename[20], tempstr[20];

#ifdef USE_MPI
	//MPI2: Some code from the main branch might have been removed in the MPI version. Please check.
	//MPI2: At some point, change all loops to run between 2 variables: Begin to End, which will be decided based on if it is MPI or not.
	//int myid, procs; //Declared in cmc.h to make it available across all files
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	//int mpiBegin, mpiEnd; //declared in cmc_vars.h
	//int *mpiDisp, *mpiLen; //declared in cmc_vars.h
	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));
#else
	//the variable(Macro) PROC has to be changed everytime the no.of processors are changed.
	procs = PROC;
	//MPI2: Compared random numbers of serial and parallel. Seem to match perfectly.
#endif

	Start = (int *) malloc(procs * sizeof(int));
	End = (int *) malloc(procs * sizeof(int));

	create_timing_files();

	/* set some important global variables */
	set_global_vars1();

	/* trap signals */
	trap_sigs(); // captures i/o signals to close

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);

	//MPI2: Storing time to compute total time.
	//timeStart2(&timeS);

	get_star_data(argc, argv, rng);

	set_rng_states();

	mpiInitBcastGlobArrays();

	calc_sigma_new();

	/* calculate central quantities */
	calc_central_new();

	bin_vars_calculate();

	/* print out binary properties to a file */
	//MPI2: skipping outputs for initial steps of MPI
	//print_initial_binaries();

	orbit_r = R_MAX;

#ifdef USE_MPI
	//MPI2: Calculating indices which will be used in all loops till beginning of the main loop. The value 20 depends on the p value used in calc_sigma_new()
   mpiFindIndicesCustom( clus.N_MAX, 20, myid, &mpiBegin, &mpiEnd );
#endif

	calc_potential_new();

	//Calculating disp and len for mimcking parallel rng.
	findLimits( clus.N_MAX, 20 );

	total_bisections= 0;

	/*
		Skipping search grid for MPI
		if (SEARCH_GRID) 
		search_grid_update(r_grid);
	 */	

	/* compute energy initially */
	//MPI2: do on root.
	//Bharath: doesnt have to be done on root since all procs have the star array, and they probably need the 0th star's values.
	star[0].E = star[0].J = 0.0;

	compute_energy_new();

	//MPI2: The following line has been moved after load_fits_file_data() - above
	//	star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;

/*
	//MPI2: ignore for now
	comp_mass_percent(); //used for diagnostics. needs neighouring particles to compute cum. sum.
	//MPI2: ignore for now
	comp_multi_mass_percent();
*/

	/* If we don't set it here, new stars created by breaking binaries (BSE) will
	 * end up in the wrong place */
	clus.N_MAX_NEW = clus.N_MAX;
	/* initialize stellar evolution things */
	DMse = 0.0;

	if (STELLAR_EVOLUTION > 0) {
		//MPI2: Ignore binaries. For binaries, there is going to be a problem shuffling binary array. To be thought about later.
		stellar_evolution_init(); //whole stellar evol. part does not need any data from other particles
	}

	//MPI2: Binaries. Ignoring for now.
	update_vars(); //might need communication for bin. index array. needs cum.sum.

	//root node & Bcast? Does not have to be done on root as per Stefan.	
	times(&tmsbufref);

	/* calculate central quantities */
	calc_central_new();

	/* can skip for MPI
		if (WRITE_EXTRA_CORE_INFO) {
		no_remnants= no_remnants_core(6);
		}
	 */

	calc_clusdyn_new();

	/* Printing Results for initial model */
	//skip outputs for MPI
	//print_results();

	/* print handy script for converting output files to physical units */
	//skip outputs for MPI
	//print_conversion_script();

#ifdef USE_CUDA
	cuInitialize();
#endif

	timeStart2(&timeS);

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) 
	{
		/* calculate central quantities */
		calc_central_new();

		calc_timestep(rng);
		//MPI2: Setting the timestep (hardcoding) for testing.
		if(myid==0)
			printf("Dt = %.14g\n", Dt);

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		clus.N_MAX_NEW = clus.N_MAX;

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		//MPI2: Tested for outputs: vr, vt, E and J. Tests performed with same seed for rng of all procs. Check done only for proc 0's values as others cant be tested due to rng. Must test after rng is replaced.
		if (PERTURB > 0)
			dynamics_apply(Dt, rng);


		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
			will disappear! */
		clus.N_MAX_NEW++;

		/* evolve stars up to new time */
		DMse = 0.0;

		//MPI2: Tested for outputs: rad, m E. Check if rng is used at all. Testing done only for proc 0.
		if (STELLAR_EVOLUTION > 0)
			do_stellar_evolution(rng);

		Prev_Dt = Dt;

		/* some numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation1();

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		toy_rejuvenation();

		/* this calls get_positions() */
		tidally_strip_stars1();

#ifdef USE_MPI
		strcpy(filename, "test_rng_par");
		strcpy(tempstr, filename);
		sprintf(num, "%d", myid);
		strcat(tempstr, num);
		strcat(tempstr, ".dat");
		for( i = 0; i < procs; i++ )
		{
			if(myid == i)
			{
				//printf("Start[i]=%d\tend=\%d\n", Start[i], End[i]);
				ftest = fopen( tempstr, "w" );
				for( j = Start[i]; j <= End[i]; j++ )
					fprintf(ftest, "%ld\t%.18g\n", j, star[j].vrnew );
				fclose(ftest);
			}
		}
		if(myid==0)
			system("./process.sh");
#else
		strcpy(tempstr, "test_rng_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].vrnew );
		fclose(ftest);
#endif

		/* more numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation2();

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		ComputeIntermediateEnergy();

		pre_sort_comm();

		tidally_strip_stars2();

/*
#ifdef USE_MPI
	printf("myid = %d\trn = %ld\n",myid, rng_t113_int_new(curr_st)) ;
#else
	for(i=0; i<procs; i++)
		printf("i = %d\trn = %ld\n",i, rng_t113_int_new(&st[i])) ;
#endif
*/
		strcpy(funcName, "qsorts");
		static double timeTotLoc;
		timeStart();

#ifdef USE_MPI
		if(myid==0)
#endif
			/* Sorting stars by radius. The 0th star at radius 0 
				and (N_STAR+1)th star at SF_INFINITY are already set earlier.
			 */
			//MPI2: Only running till N_MAX for now since no new stars are created, later has to be changed to N_MAX_NEW.
			//MPI2: Changin to N_MAX+1, later change to N_MAX_NEW+1
			qsorts(star+1,clus.N_MAX); //parallel sort

		timeEnd(fileTime, funcName, &timeTotLoc);

		calc_potential_new2();

#ifdef USE_MPI
		//MPI2: Calculating new indices which will be used in all loops till end of next timestep (qsorts).
		mpiFindIndicesCustom( clus.N_MAX, 20, myid, &mpiBegin, &mpiEnd );
#endif

		//Calculating Start and End values for each processor for mimcking parallel rng.
		findLimits( clus.N_MAX, 20 );

		set_velocities3();

		post_sort_comm();

		//commenting out for MPI
		/*
			if (SEARCH_GRID)
			search_grid_update(r_grid);
		 */


		/*
		//MPI2: Ignoring for MPI
#ifndef USE_MPI
comp_mass_percent();
comp_multi_mass_percent();
#endif
		 */

		compute_energy_new();

#ifdef USE_MPI
		for (i=mpiBegin; i<=mpiEnd; i++)
#else
		for (i = 1; i <= clus.N_MAX; i++) 
#endif
			/* reset interacted flag */
			star[i].interacted = 0;

		//MPI2: Binaries. Ignore for now.
		/* update variables, then print */
		update_vars();

		tcount++;

		calc_clusdyn_new();

		//commenting out for MPI
		/*
			if (WRITE_EXTRA_CORE_INFO) {
			no_remnants= no_remnants_core(6);
			}
		 */

		//commenting out for MPI
		//print_results();
		/* take a snapshot, we need more accurate 
		 * and meaningful criterion 
		 */
		if(tcount%SNAPSHOT_DELTACOUNT==0) {
			print_2Dsnapshot();
			if (WRITE_STELLAR_INFO) {
				write_stellar_data();
			}
		}		

	} /* End WHILE (time step iteration loop) */

	timeEnd2(fileTime, "TotalTime", &timeS, &timeE, &timeT);


#ifdef USE_MPI
	if(myid==0)
#endif
	{
		ftest = fopen("mpi_globvar.dat","w");
		fprintf(ftest, "clus.N_MAX_NEW=%ld\nclus.N_MAX=%ld\nTidalMassLoss=%g\nOldTidalMassLoss=%g\nPrev_Dt=%g\nEescaped=%g\nJescaped=%g\nEintescaped=%g\nEbescaped=%g\nMtotal=%g\ninitial_total_mass=%g\nDMse=%g\nDMrejuv=%g\ncenma.m=%g\ncenma.m_new=%g\ncenma.E=%g\ntcount=%ld\nStepCount=%ld\nsnap_num=%ld\nEcheck=%ld\nnewstarid=%ld\n",
				clus.N_MAX_NEW,
				clus.N_MAX,
				TidalMassLoss,
				OldTidalMassLoss,
				Prev_Dt,
				Eescaped,
				Jescaped,
				Eintescaped,
				Ebescaped, 
				Mtotal, 
				initial_total_mass,
				DMse,
				DMrejuv,
				cenma.m,
				cenma.m_new,
				cenma.E,
				tcount,
				StepCount,
				snap_num,
				Echeck, 
				newstarid);
		fclose(ftest);
	}

	//root node?
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

#ifdef USE_MPI
	free(mpiDisp);
	free(mpiLen);
	free(curr_st);
#else
	free(st);
#endif
	free(Start);
	free(End);

	if (SEARCH_GRID)
		search_grid_free(r_grid);
#ifdef DEBUGGING
	g_hash_table_destroy(star_ids);
	g_array_free(id_array, TRUE);
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return(0);
}
