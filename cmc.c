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

/*
#ifndef SIMUL
#define SIMUL 1 //to mimic rng of the serial version to the parallel version. This has to be given by "make simul=<PROCS>" while compiing.
#endif
*/

int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i, j;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	//Temp file handle for debugging
	FILE *ftest;
	//MPI2: Some variables to assist debugging
	char num[5],filename[20], tempstr[20];

#ifdef USE_MPI
	//MPI2: Some code from the main branch might have been removed in the MPI version. Please check.
	//MPI2: At some point, change all loops to run between 2 variables: Begin to End, which will be decided based on if it is MPI or not.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));

	if(myid==0)
		new_size = (int *)calloc(procs, sizeof(int)); 
#else
	procs = 1;
	created_star_dyn_node = (int *) calloc(procs, sizeof(int));
	created_star_se_node = (int *) calloc(procs, sizeof(int));
#endif

	Start = (int *) calloc(procs, sizeof(int));
	End = (int *) calloc(procs, sizeof(int));

	create_timing_files();

	/* set some important global variables */
	set_global_vars1();

	/* trap signals */
	trap_sigs(); // captures i/o signals to close

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);

	get_star_data(argc, argv, rng);

	N_b_OLD = N_b;
	N_b_NEW = N_b;

	set_rng_states();

	mpiInitBcastGlobArrays();

	calc_sigma_new();

	/* calculate central quantities */
	calc_central_new();

	/* print out binary properties to a file */
	//MPI2: skipping outputs for initial steps of MPI
	//print_initial_binaries();

	orbit_r = R_MAX;

	calc_potential_new();

	//MPI2: Setting this because in the MPI version only _r is set, and will create problem while sorting.
	star[clus.N_MAX + 1].r = SF_INFINITY;
	star[clus.N_MAX + 1].phi = 0.0;

	/* MPI2: Calculating disp and len for mimcking parallel rng */
	findLimits( clus.N_MAX, 20 );

	alloc_bin_buf();

	distr_bin_data();

	bin_vars_calculate();

	total_bisections= 0;

	/*
		Skipping search grid for MPI
		if (SEARCH_GRID) 
		search_grid_update(r_grid);
	 */	

	/* compute energy initially */
	star[0].E = star[0].J = 0.0;

	compute_energy_new();

	set_energy_vars();

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

#ifndef USE_MPI
	for(i=0; i<procs; i++)
		created_star_se_node[i] = 0;
#endif

	//MPI2: Binaries. Ignoring for now.
	update_vars(); //might need communication for bin. index array. needs cum.sum.

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
	// MPI2: Not file outputs for MPI
	print_results();

	/* print handy script for converting output files to physical units */
	//skip outputs for MPI
	//print_conversion_script();

/*
#ifdef USE_CUDA
	cuInitialize();
#endif

#ifdef USE_MPI
		strcpy(filename2, "test_E_par");
		strcpy(tempstr2, filename2);
		sprintf(num2, "%d", myid);
		strcat(tempstr2, num2);
		strcat(tempstr2, ".dat");
#else
		strcpy(tempstr2, "test_E_ser.dat");
#endif
		ftest2 = fopen( tempstr2, "w" );
*/
	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) 
	{
#ifndef USE_MPI
		for(i=0; i<procs; i++)
		{
			created_star_dyn_node[i] = 0;
			created_star_se_node[i] = 0;
		}
#endif

#ifdef USE_MPI
		Eescaped_old = Eescaped;
		Jescaped_old = Jescaped;
		Eintescaped_old = Eintescaped;
		Ebescaped_old = Ebescaped;
		TidalMassLoss_old = TidalMassLoss;
		Etidal_old = Etidal;

		Eescaped = 0.0;
		Jescaped = 0.0;
		Eintescaped = 0.0;
		Ebescaped = 0.0;
		TidalMassLoss = 0.0;
		Etidal = 0.0;
#endif

		/* calculate central quantities */
		calc_central_new();

		calc_timestep(rng);
#ifdef USE_MPI
		if(myid==0){
			strcpy(tempstr, "test_dt_par.dat");
			ftest = fopen( tempstr, "a" );
			fprintf(ftest, "%ld\t%.18g\n", tcount, Dt );
			fclose(ftest);
		}
#else
		strcpy(tempstr, "test_dt_ser.dat");
		ftest = fopen( tempstr, "a" );
		fprintf(ftest, "%ld\t%.18g\n", tcount, Dt );
		fclose(ftest);
#endif
		//Dt = 0.0001065233252096;
#ifdef USE_MPI
if(myid==0)
#endif
		printf("FINAL TIMESTEP = %.14g\n", Dt);

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		clus.N_MAX_NEW = clus.N_MAX;

/*
#ifdef USE_MPI
		printf("id = %d\tBefore dyn_apply N_MAX = %ld\tN_MAX_NEW = %ld\n",myid, clus.N_MAX, clus.N_MAX_NEW);
#else
		printf("Before dyn_apply N_MAX = %ld\tN_MAX_NEW = %ld\n", clus.N_MAX, clus.N_MAX_NEW);
#endif
*/

//printf("0\ncenma.E=%.18g\nEescaped=%.18g\nEbescaped=%.18g\nEintescaped=%.18g\n", cenma.E, Eescaped, Ebescaped, Eintescaped );
		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		//MPI2: Tested for outputs: vr, vt, E and J. Tests performed with same seed for rng of all procs. Check done only for proc 0's values as others cant be tested due to rng. Must test after rng is replaced.
		if (PERTURB > 0)
			dynamics_apply(Dt, rng);

/*
#ifndef USE_MPI
		for(i=0; i<procs; i++)
			printf("node %ld=%d, %d\t", i, created_star_dyn_node[i], created_star_se_node[i]);
		printf("\n");
#endif

#ifdef USE_MPI
		printf("id = %d\tAfter dyn_apply N_MAX = %ld\tN_MAX_NEW = %ld\n",myid, clus.N_MAX, clus.N_MAX_NEW);
#else
		printf("After dyn_apply N_MAX = %ld\tN_MAX_NEW = %ld\n", clus.N_MAX, clus.N_MAX_NEW);
#endif
*/

		/* evolve stars up to new time */
		DMse = 0.0;

		//MPI2: Tested for outputs: rad, m E. Check if rng is used at all. Testing done only for proc 0.
		if (STELLAR_EVOLUTION > 0)
			do_stellar_evolution(rng);

		Prev_Dt = Dt;

		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
			will disappear! */
		clus.N_MAX_NEW++;

		/* some numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation1();

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		toy_rejuvenation();

/*
#ifndef USE_MPI
		for(i=0; i<procs; i++)
			printf("node %ld=%d, %d\t", i, created_star_dyn_node[i], created_star_se_node[i]);
		printf("\n");
#endif

#ifdef USE_MPI
		printf("id = %d\tAfter SE N_MAX = %ld\tN_MAX_NEW = %ld\n",myid, clus.N_MAX, clus.N_MAX_NEW);
#else
		printf("After SE N_MAX = %ld\tN_MAX_NEW = %ld\n", clus.N_MAX, clus.N_MAX_NEW);
#endif
*/

//printf("1\ncenma.E=%.18g\nEescaped=%.18g\nEbescaped=%.18g\nEintescaped=%.18g\n", cenma.E, Eescaped, Ebescaped, Eintescaped );
		/* this calls get_positions() */
		tidally_strip_stars1();

//printf("rmax = %.18g\n",rmax);

		/* more numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation2();

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		ComputeIntermediateEnergy();
//printf("2\ncenma.E=%.18g\nEescaped=%.18g\nEbescaped=%.18g\nEintescaped=%.18g\n", cenma.E, Eescaped, Ebescaped, Eintescaped );

/*
		for(i=clus.N_MAX+2; i<=clus.N_MAX_NEW; i++) 
			if(star[i].binind>0)
				printf("myid = %d\t%ld\t%.18g\n", myid, i, binary[star[i].binind].a );
*/

		pre_sort_comm();

		collect_bin_data();		

#ifdef USE_MPI
	if(myid==0){
		strcpy(tempstr, "test_var_par.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= N_b_NEW+5; i++ )
			//if(star[i].binind>0)
				fprintf(ftest, "%ld\t%.18g\n", i, binary[i].a );
		fclose(ftest);
	}
#else
	//printf("Nmaxnew = %ld\n", clus.N_MAX_NEW);
	strcpy(tempstr, "test_var_ser.dat");
	ftest = fopen( tempstr, "w" );
	for( i = 1; i <= N_b_NEW+5; i++ )
		//if(star[i].binind>0)
			fprintf(ftest, "%ld\t%.18g\n", i, binary[i].a );
	fclose(ftest);
#endif

		tidally_strip_stars2();


	if(myid==0){
		strcpy(tempstr, "test_rng_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX_NEW; i++ )
			//if(star[i].binind>0)
				fprintf(ftest, "%ld\t%.18g\n", i, star[i].r);
		fclose(ftest);
	}

		for( i = 1; i <= clus.N_MAX_NEW; i++ )
			if(star[i].r==0)
				printf("######## FOUND BUG!!!!!!! %d\n\n", i);

//printf("3\ncenma.E=%.18g\nEescaped=%.18g\nEbescaped=%.18g\nEintescaped=%.18g\n", cenma.E, Eescaped, Ebescaped, Eintescaped );
		qsorts_new();

#ifdef USE_MPI
	if(myid==0){
		strcpy(tempstr, "test_rng_par.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			//if(star[i].binind>0)
				fprintf(ftest, "%ld\t%.18g\n", i, star[i].r);
		fclose(ftest);
	}
#else
	//printf("Nmaxnew = %ld\n", clus.N_MAX_NEW);
	strcpy(tempstr, "test_rng_ser.dat");
	ftest = fopen( tempstr, "w" );
	for( i = 1; i <= clus.N_MAX; i++ )
		//if(star[i].binind>0)
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].r );
	fclose(ftest);
#endif

		calc_potential_new2();

		//Calculating Start and End values for each processor for mimcking parallel rng.
		findLimits( clus.N_MAX, 20 );

		set_velocities3();

		post_sort_comm();

		distr_bin_data();
//printf("4\ncenma.E=%.18g\nEescaped=%.18g\nEbescaped=%.18g\nEintescaped=%.18g\n", cenma.E, Eescaped, Ebescaped, Eintescaped );

		//commenting out for MPI
		/*
			if (SEARCH_GRID)
			search_grid_update(r_grid);
	
		//MPI2: Ignoring for MPI
		comp_mass_percent();
		comp_multi_mass_percent();
		 */

		compute_energy_new();

//printf("ETOTAL = %.14g\n\n", Etotal.tot);

		reset_interaction_flags();

		//MPI2: Binaries. Ignore for now.
		/* update variables, then print */
		update_vars();

		tcount++;

		calc_clusdyn_new();

#ifdef USE_MPI
		strcpy(filename, "test_out_par");
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
					//if(star[j].binind>0)
						//fprintf(ftest, "%ld\t%.18g\n", j, binary[star[j].binind].a);
					fprintf(ftest, "%ld\t%.18g\n", j, star_r[j]);
				fclose(ftest);
			}
		}
		if(myid==0)
			system("./process2.sh");
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			//if(star[i].binind>0)
				//fprintf(ftest, "%ld\t%.18g\n", i, binary[star[i].binind].a);
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].r);
		fclose(ftest);
#endif

		//MPI2: Commenting out for MPI
		/*
			if (WRITE_EXTRA_CORE_INFO) {
			no_remnants= no_remnants_core(6);
			}
		 */

#ifdef USE_MPI
if(myid==0)
#endif
	printf("========> N_MAX =  %ld\n", clus.N_MAX);

		//MPI2: Commenting out for MPI
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

	} /* End WHILE (time step iteration loop) */

	times(&tmsbuf);

	dprintf("Usr time = %.6e ", (double)
			(tmsbuf.tms_utime-tmsbufref.tms_utime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time = %.6e\n", (double)
			(tmsbuf.tms_stime-tmsbufref.tms_stime)/sysconf(_SC_CLK_TCK));
	dprintf("Usr time (ch) = %.6e ", (double)
			(tmsbuf.tms_cutime-tmsbufref.tms_cutime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time (ch)= %.6e seconds\n", (double)
			(tmsbuf.tms_cstime-tmsbufref.tms_cstime)/sysconf(_SC_CLK_TCK));

#ifdef USE_MPI
	if(myid==0)
#endif
	printf("The total number of bisections is %li\n", total_bisections);

	/* free RNG */
	gsl_rng_free(rng);

#ifdef USE_CUDA
	cuCleanUp();
#endif

	/* flush buffers before returning */
	close_buffers();
	free_arrays();

/*
#ifdef USE_MPI
	if(myid==0)
		system("./process2.sh");
#endif
	//fclose(ftest2);
*/

#ifdef USE_MPI
	free(mpiDisp);
	free(mpiLen);
	free(curr_st);
	free(binary_buf);
	free(num_bin_buf);
#else
	free(st); //commenting because it throws some error
#endif

	//free(Start);
	//free(End);


	//free(multi_mass_r[i]);
/*
	free(multi_mass_r);
	free(sigma_array.r);
	free(sigma_array.sigma);
	free(mass_bins);
*/
	if(SNAPSHOT_WINDOWS)
		free(SNAPSHOT_WINDOWS);
/*
	if(SNAPSHOT_WINDOW_UNITS)
		free(SNAPSHOT_WINDOW_UNITS);
*/
	if (SEARCH_GRID)
		search_grid_free(r_grid);

	if(zpars)
		free(zpars);

#ifdef DEBUGGING
	g_hash_table_destroy(star_ids);
	g_array_free(id_array, TRUE);
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return(0);
}
