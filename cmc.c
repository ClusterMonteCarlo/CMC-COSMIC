/* -*- linux-c -*- */
/* trivial */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#define _MAIN_
#include "cmc_vars.h"
#include <string.h>

#ifdef USE_CUDA
#include "cuda/cmc_cuda.h"
#endif

int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	//MPI3: There might be overhead if timing is done due to MPI_Barriers.
	double tmpTimeStart, tmpTimeStart_full;
	double t_full=0.0, t_init=0.0, t_sort=0.0, t_ener=0.0, t_se=0.0, t_dyn=0.0, t_orb=0.0, t_oth=0.0, t_filemer=0.0;

#ifdef USE_MPI
	//MPI2: Some code from the main branch might have been removed in the MPI version. Please check.
	//MPI2: At some point, change all loops to run between 2 variables: Begin to End, which will be decided based on if it is MPI or not.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	//MPI2: Used for Gathering particles on root node before sorting
	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));

	if(myid==0)
	{
		disp = (int*) malloc(procs * sizeof(int));
		len = (int*) malloc(procs * sizeof(int));
		new_size = (long*) malloc(procs * sizeof(long));
	}

/*
//code for using gdb with mpi.
//http://www.open-mpi.org/faq/?category=debugging#serial-debuggers
{
    int ii = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == ii)
        sleep(5);
}
*/

#else
	procs = 1;
	created_star_dyn_node = (int *) calloc(procs, sizeof(int));
	created_star_se_node = (int *) calloc(procs, sizeof(int));
	DMse_mimic = (double *) calloc(procs, sizeof(double));
#endif

	tmpTimeStart_full = timeStartSimple();
	tmpTimeStart = timeStartSimple();
	create_timing_files();

	/* set some important global variables */
	set_global_vars1();

	/* trap signals */
	trap_sigs(); //captures i/o signals to close

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);

	//Currently doing on all nodes to avoid broadcasting of global variables. Later, might have to be split into 2 or more functions, and store all global variables into a structure for easy broadcast.
	/* parse input */
	parser(argc, argv, rng); //to do parallel i/o

	Start = (int *) calloc(procs, sizeof(int));
	End = (int *) calloc(procs, sizeof(int));

	/* MPI2: Calculating Start and End for mimcking parallel rng */
	findLimits( cfd.NOBJ, 20 );

	//MPI3: Allocate N/procs + 10% for each node. Also allocate a separate buffer array for receiving ghost particles. File I/O, each process takes its slice of data. Also, assemble the global arrays - _m, and _r.
	get_star_data(argc, argv, rng);

	Start = (int *) calloc(procs, sizeof(int));
	End = (int *) calloc(procs, sizeof(int));

	findLimits( clus.N_MAX, 20 );

	N_b_OLD = N_b;
	N_b_NEW = N_b;

	set_rng_states();

	//MPI3: Need to check if this function and the sigma variable is reqd at all! Use the buffer for collecting ghost particles.
	//MPI3: Ignoring for now. Parallelize after verifying if the variable is reqd at all.
	//calc_sigma_new();

	/* calculate central quantities */
	calc_central_new();

	/* print out binary properties to a file */
	//MPI2: skipping  file outputs for now.
	//print_initial_binaries();

	orbit_r = R_MAX;

	calc_potential_new();

	//MPI2: Setting this because in the MPI version only _r is set, and will create problem while sorting.
#ifndef USE_MPI
	star[clus.N_MAX + 1].r = SF_INFINITY;
	star[clus.N_MAX + 1].phi = 0.0;
#endif

	//MPI3: Temporarily avoiding binary stuff.
	if(N_b !=0)
	{
		alloc_bin_buf();

		distr_bin_data();

		bin_vars_calculate();
	}

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


#ifndef USE_MPI
	for(i=0; i<procs; i++)
		DMse_mimic[i] = 0.0;
#endif

	if (STELLAR_EVOLUTION > 0) {
		stellar_evolution_init();
	}

#ifndef USE_MPI
	for(i=0; i<procs; i++)
		created_star_se_node[i] = 0;
#endif

	//OPT: M_b E_b calculated twice? Check for redundancy.
	if(N_b!=0)
		update_vars();

	times(&tmsbufref);

	//OPT: Check for redundancy. Ask Stefan
	/* calculate central quantities */
	calc_central_new();

	/* can skip for MPI
		if (WRITE_EXTRA_CORE_INFO) {
		no_remnants= no_remnants_core(6);
		}
	 */

	calc_clusdyn_new();

	//MPI3: Commenting out outputs for now.
	/* Printing Results for initial model */
	//print_results();

	/* print handy script for converting output files to physical units */
	//skip outputs for MPI
	//print_conversion_script();

#ifdef USE_CUDA
	cuInitialize();
#endif

#ifdef USE_MPI
	if (STELLAR_EVOLUTION > 0) {
		mpiFindDispAndLenCustom( clus.N_MAX, 20, mpiDisp, mpiLen );

		for(i=0;i<procs;i++)
			mpiLen[i] *= sizeof(double); 

		//MPI3: THis needs to be done only if SE is on? If yes, put a condition.
		//MPI3: Can be replaced by AllGatherv. For now, keep it the way it is. The same problem will reappear after sorting, where we need to do the same for the 3 global arrays, r, phi, and m. Have to think abt a novel solution.
		MPI_Status stat;
		if(myid!=0)
			MPI_Send(&star_m[mpiDisp[myid]], mpiLen[myid], MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		else
			for(i=1;i<procs;i++)
				MPI_Recv(&star_m[mpiDisp[i]], mpiLen[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &stat);
		MPI_Bcast(star_m, clus.N_MAX+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
#endif
	timeEndSimple(tmpTimeStart, &t_init);

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) 
	{

		tmpTimeStart = timeStartSimple();
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

		//=========================================
		//MPI3: COMPLETED PARALLELIZATION TILL HERE>
		//=========================================
		/*if(myid==1)
			for(i=0; i<procs; i++)
				printf("%ld %d %d\n", i, Start[i], End[i]);
		*/
		//=========================================

		calc_timestep(rng);

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
#ifdef USE_MPI
		clus.N_MAX_NEW = mpiEnd-mpiBegin+1;
#else
		clus.N_MAX_NEW = clus.N_MAX;
#endif
		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		//MPI2: Tested for outputs: vr, vt, E and J. Tests performed with same seed for rng of all procs. Check done only for proc 0's values as others cant be tested due to rng. Must test after rng is replaced.
		if (PERTURB > 0)
			dynamics_apply(Dt, rng);
		timeEndSimple(tmpTimeStart, &t_dyn);

		tmpTimeStart = timeStartSimple();
		//MPI2: Tested for outputs: rad, m E. Check if rng is used at all. Testing done only for proc 0.
		if (STELLAR_EVOLUTION > 0)
			do_stellar_evolution(rng);
		timeEndSimple(tmpTimeStart, &t_se);

		tmpTimeStart = timeStartSimple();
		Prev_Dt = Dt;

		/* evolve stars up to new time */
		DMse = 0.0;
#ifndef USE_MPI
		for(i=0; i<procs; i++)
			DMse_mimic[i] = 0.0;
#endif

		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
			will disappear! */
#ifndef USE_MPI
		clus.N_MAX_NEW++;
#endif

		/* some numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation1();

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		toy_rejuvenation();
		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		/* this calls get_positions() */
		new_orbits_calculate();
		timeEndSimple(tmpTimeStart, &t_orb);

		tmpTimeStart = timeStartSimple();
		/* more numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation2();

		tidally_strip_stars();

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		ComputeIntermediateEnergy();

		//pre_sort_comm();

		//collect_bin_data();		

		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		qsorts_new();

		post_sort_comm();
		timeEndSimple(tmpTimeStart, &t_sort);

		tmpTimeStart = timeStartSimple();
		calc_potential_new();

		//Calculating Start and End values for each processor for mimcking parallel rng.
		findLimits( clus.N_MAX, 20 );
		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		energy_conservation3();
		timeEndSimple(tmpTimeStart, &t_ener);

		//distr_bin_data();

		//commenting out for MPI
		/*
			if (SEARCH_GRID)
			search_grid_update(r_grid);
	
		//MPI2: Ignoring for MPI
		comp_mass_percent();
		comp_multi_mass_percent();
		 */

		tmpTimeStart = timeStartSimple();
		compute_energy_new();

		reset_interaction_flags();

		//MPI2: Commenting out for MPI
		/* update variables, then print */
		//update_vars();

		tcount++;

		calc_clusdyn_new();

		//MPI2: Commenting out for MPI
		/*
			if (WRITE_EXTRA_CORE_INFO) {
			no_remnants= no_remnants_core(6);
			}
		 */

/* TESTING FOR KEVIN */
/*
#ifdef USE_MPI
		// Only proc with id 0 prints out.
		if(myid==0)
		{
			strcpy(tempstr, "test_out_par.dat");
			ftest = fopen( tempstr, "w" );
			for( i = 1; i <= clus.N_MAX; i++ )
				fprintf(ftest, "%ld\t%.18g\n",i, star_r[i]);
			fclose(ftest);
		}
		MPI_Barrier(MPI_COMM_WORLD);
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].r);
		fclose(ftest);
#endif
*/
		print_results();
		timeEndSimple(tmpTimeStart, &t_oth);

		/* take a snapshot, we need more accurate 
		 * and meaningful criterion 
		 */
/*
#ifdef USE_MPI
	if(myid==0)
#endif
		if(tcount%SNAPSHOT_DELTACOUNT==0) {
			print_2Dsnapshot();
			if (WRITE_STELLAR_INFO) {
				write_stellar_data();
			}
		}		
*/
	} /* End WHILE (time step iteration loop) */

	tmpTimeStart = timeStartSimple();
	mpi_files_merge();
	timeEndSimple(tmpTimeStart, &t_filemer);

	times(&tmsbuf);
	timeEndSimple(tmpTimeStart_full, &t_full);

#ifdef USE_MPI
	if(myid==0)
#endif
	{
		printf("******************************************************************************\n");
		printf("Time Taken:\n------------------------\n%.4lf:\tInitialization\n%.4lf:\tDynamics\n%.4lf:\tStellar Evolution\n%.4lf:\tOrbit Calculation\n%.4lf:\tSorting\n%.4lf:\tEnergy Conservation\n%.4lf:\tOthers\n%.4lf:\tFiles Merge\n------------------------\n%.4lf:\tTotal\n------------------------\n", t_init, t_dyn, t_se, t_orb, t_sort, t_ener, t_oth, t_filemer, t_full);

/*
		printf("%.4lf seconds of processing\n", t_full);
		printf("%.4lf seconds took by initialization\n", t_init);
		printf("%.4lf seconds took by sorting\n", t_sort);
		printf("%.4lf seconds took by energy conservation\n", t_ener);
		printf("%.4lf seconds took by dynamics\n", t_dyn);
		printf("%.4lf seconds took by SE\n", t_se);
		printf("%.4lf seconds took by orb. calc.\n", t_orb);
		printf("%.4lf seconds took by others\n", t_oth);
		printf("%.4lf seconds took by file merge\n", t_filemer);
*/
	}

	dprintf("Usr time = %.6e ", (double)
			(tmsbuf.tms_utime-tmsbufref.tms_utime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time = %.6e\n", (double)
			(tmsbuf.tms_stime-tmsbufref.tms_stime)/sysconf(_SC_CLK_TCK));
	dprintf("Usr time (ch) = %.6e ", (double)
			(tmsbuf.tms_cutime-tmsbufref.tms_cutime)/sysconf(_SC_CLK_TCK));
	dprintf("Sys time (ch)= %.6e seconds\n", (double)
			(tmsbuf.tms_cstime-tmsbufref.tms_cstime)/sysconf(_SC_CLK_TCK));

	/* free RNG */
	//gsl_rng_free(rng);

#ifdef USE_CUDA
	cuCleanUp();
#endif

	/* flush buffers before returning */
	close_buffers();
	//free_arrays();

/*
#ifdef USE_MPI
	free(mpiDisp);
	free(mpiLen);
	free(curr_st);
	free(binary_buf);
	free(num_bin_buf);
#else
	free(st); //commenting because it throws some error
#endif
*/
	free(Start);
	free(End);

	//if (SEARCH_GRID)
	//	search_grid_free(r_grid);

	//if(zpars)
	//	free(zpars);

#ifdef DEBUGGING
	g_hash_table_destroy(star_ids);
	g_array_free(id_array, TRUE);
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return(0);
}
/*
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
				for( j = 1; j <= End[i]-Start[i]+1; j++ )
				//for( j = mpiBegin; j <= mpiEnd; j++ )
				//for( j = 1; j <= clus.N_MAX; j++ )
					fprintf(ftest, "%ld\t%.18g\n", mpiBegin+j-1, star[j].vtnew);
					//fprintf(ftest, "%ld\t%ld\n", mpiBegin+j-1, star[j].id);
					//fprintf(ftest, "%ld\t%.18g\n", j, star_m[j]);
				fclose(ftest);
			}
		}
		if(myid==0)
			system("./process.sh");
MPI_Barrier(MPI_COMM_WORLD);
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].vtnew);
			//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
			//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
		fclose(ftest);
#endif
//return;

*/

/*
#ifdef USE_MPI 
   printf("id = %d\trng = %li\n", myid, rng_t113_int_new(curr_st));
#else
   for(i=0; i<procs; i++)
      printf("i = %d\trng = %li\n", i, rng_t113_int_new(&st[i]));
#endif
*/

/*
#ifdef USE_MPI
		if(myid==2)
		{
			strcpy(tempstr, "test_out_par.dat");
			ftest = fopen( tempstr, "w" );
			for( i = 1; i <= clus.N_MAX; i++ )
				fprintf(ftest, "%.18g\n", star_r[i]);
			//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
			//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
			fclose(ftest);
		}
#endif
MPI_Barrier(MPI_COMM_WORLD);
*/

/*
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
				{
					//if(star[j].binind>0)
						//fprintf(ftest, "%ld\t%.18g\n", j, binary[star[j].binind].a);
				if(star[j].id <= 0)
					fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\n", j, star_r[j], star[j].id, binary[star[j].binind].id1, binary[star[j].binind].id2);
				else
					fprintf(ftest, "%ld\t%.18g\t%ld\n", j, star_r[j], star[j].id);
				}
				fclose(ftest);
			}
		}
		if(myid==0)
			system("./process.sh");
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
		{
			//if(star[i].binind>0)
				//fprintf(ftest, "%ld\t%.18g\n", i, binary[star[i].binind].a);
		if(star[i].id <= 0)
			fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\n", i, star[i].r, star[i].id, binary[star[i].binind].id1, binary[star[i].binind].id2);
		else
			fprintf(ftest, "%ld\t%.18g\t%ld\n", i, star[i].r, star[i].id);
		}
		fclose(ftest);
#endif
*/

/*

#ifdef USE_MPI
	//printf("%ld----------------------\n",clus.N_MAX);
	if(myid==1)
	{
		strcpy(tempstr, "test_out_par.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			fprintf(ftest, "%ld\t%.18g\n",i, star_m[i]);
		//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
		//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
		fclose(ftest);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#else
	strcpy(tempstr, "test_out_ser.dat");
	ftest = fopen( tempstr, "w" );
	for( i = 1; i <= clus.N_MAX; i++ )
		fprintf(ftest, "%ld\t%.18g\n", i, star[i].m);
	//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
	//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
	fclose(ftest);
#endif
return;

 
#ifdef USE_MPI
		if(myid==0)
		{
			ftest = fopen("mpi_globvar.dat","w");
#else
		{
			ftest = fopen("ser_globvar.dat","w");
#endif
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
#ifdef USE_MPI
		if(myid==2)
		{
			strcpy(tempstr, "test_out_par.dat");
			ftest = fopen( tempstr, "w" );
			for( i = 1; i <= clus.N_MAX; i++ )
				fprintf(ftest, "%ld\t%.18g\n", i, star_r[i]);
			//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
			//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
			fclose(ftest);
		}

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
				for( j = 1; j <= End[i]-Start[i]+1; j++ )
					//for( j = mpiBegin; j <= mpiEnd; j++ )
					//for( j = 1; j <= clus.N_MAX; j++ )
					fprintf(ftest, "%ld\t%.18g\n", mpiBegin+j-1, star[j].E);
					//fprintf(ftest, "%ld\t%ld\n", mpiBegin+j-1, star[j].id);
					//fprintf(ftest, "%ld\t%.18g\n", j, star_m[j]);
				fclose(ftest);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==0)
		{
			char process_str[30];
			sprintf(process_str, "./process.sh %d", procs);
			system(process_str);
		}
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
			fprintf(ftest, "%ld\t%.18g\n", i, star[i].E);
			//fprintf(ftest, "%ld\t%.18g\t\n", i, star[i].r);
			//fprintf(ftest, "%ld\t%ld\t\n", i, star[i].id);
		fclose(ftest);
#endif

*/
