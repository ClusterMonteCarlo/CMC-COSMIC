/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */
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


/**
* @brief The main function.
* @details
Monte Carlo (MC) methods calculate the dynamical evolution of a collisional system of N stars in the Fokker-Planck approximation, which applies when the evolution of the cluster is dominated by two-body relaxation, and the relaxation time is much larger than the dynamical time. In practice, further assumptions of spherical symmetry and dynamical equilibrium have to be made. The Henon MC technique (Henon 1971), which is based on orbit averaging, represents a balanced compromise between realism and speed. The MC method allows for a star-by-star realization of the cluster, with its N particles representing the N stars in the cluster. Integration is done on the relaxation timescale, and the total computational cost scales as O(NlogN) (Henon 1971).

Our code here is based on the Henon-type MC cluster evolution code CMC (Cluster Monte Carlo), developed over many years by Joshi et al. (2000, 2001), Fregeau et al. (2003), Fregeau & Rasio (2007), Chatterjee et al. (2010), and Umbreit et al. (2012). CMC includes a detailed treatment of strong binary star interactions and physical stellar collisions (Fregeau & Rasio 2007), as well as an implementation of single and binary star evolution (Chatterjee et al. 2010) and the capability of handling the dynamics around a central massive black hole (Umbreit et al. 2012).

Code Overview:\n
The principal data structure is an array of structures of type star_t. One or more of these stars can be binaries. Binaries are stored in a separate array of structures of type binary_t. If an element in the star array is a binary, it's binind property/variable is assigned a non-zero value. Moreover, the value of this variable is assigned in such a way that it indicates the index of the binary array which holds the properties of this binary.

CMC can be run both serially as well as parallely.\n

Parallelization Strategy:\n
There are various routines which have varying dependencies and accesses between elements of the data structures. To minmize inter-process communication required by these routines as much as possible, we partition the data such that the number of stars held by each processor is a multiple of MIN_CHUNK_SIZE (input parameter, defaults to 20, refer http://adsabs.harvard.edu/abs/2013ApJS..204...15P for a justification for the choice of this value, and detailed explanation). The remainder of stars which is < MIN_CHUNK_SIZE after division goes to the last processor.

Most routines do computations on a star by star basis, and access only the properties of the star being processed. However, the values of gravitational potential, masses of the stars and radial positions of the stars are accessed in a complex, data-dependent manner throughout the code. Hence, we store these values in separate arrays and duplicate/copy these across all processors. Since these properties change during a timestep, we synchronize these arrays at the end of each time step.

Like mentioned above, most routines perform computations on a star by star basis. The only routine that differs from this pattern is the sorting routine that sorts all the stars by their radial positions. We use a parallel sorting algorithm called Sample Sort (see http://adsabs.harvard.edu/abs/2013ApJS..204...15P for a detailed discussion).

The parallel sorting routine does the sorting and automatically redistributes the data to the processors, however, there is no control over the number of stars that will end up on each processor. So, after the sorting takes place, we exchange data between processors to adhere to the data partitioning scheme mentioned above.
*
* @param argc Input argument count
* @param argv[] Input argument list:
USAGE:
  <path>/cmc [options...] <input_file> <output_file_prefix>
  mpirun -np <procs> <path>/cmc  <input_file> <output_file_prefix>

OPTIONS:
  -q --quiet   : do not print diagnostic info to stdout
  -d --debug   : turn on debugging
  -V --version : print version info
  -h --help    : display this help text
  -s --streams	:Run with multiple random streams. To mimic the parallel version with the given number of processors
*
* @return indicates how the program exited. 0 implies normal exit, and nonzero value implies abnormal termination.
*/
int main(int argc, char *argv[])
{
	struct tms tmsbuf, tmsbufref;
	long i;
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	//MPI: Variables used for measuring timing if the input parameter TIMER is set to 1. Caveat: There might be some small overhead if the parallel code is timed due to the barriers (synchronization points) used.
	double tmpTimeStart, tmpTimeStart_init, tmpTimeStart_full;
	double t_full=0.0, t_init=0.0, t_cen_calc=0.0, t_timestep=0.0, t_dyn=0.0, t_se=0.0, t_orb=0.0, t_tid_str=0.0, t_sort=0.0, t_postsort_comm=0.0, t_pot_cal=0.0, t_ener_con3=0.0, t_calc_io_vars1=0.0, t_calc_io_vars2=0.0, t_io=0.0, t_io_ignore=0.0, t_comp_ener=0.0, t_upd_vars=0.0, t_oth=0.0;
	t_sort_only=0.0;
	t_sort_lb=0.0;
	t_sort_lsort1=0.0;
	t_sort_splitters=0.0;
	t_sort_a2a=0.0;
	t_sort_lsort2=0.0;
	t_sort_oth=0.0;
	t_comm=0.0;

#ifdef USE_MPI
	//MPI: Some code from the main branch might have been removed in the MPI version. Please check.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* MPI: These variables are used for storing data partitioning related information in the parallel version. These are in particular useful for some MPI communication calls. */
	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));
#else
	//If the serial version is being compiled, setting procs to 1.
	procs = 1;
#endif

	/* Starting timer to measure the overall time taken */
#ifdef USE_MPI
	tmpTimeStart_init = MPI_Wtime();
#else
	struct timeval tim;
	gettimeofday(&tim, NULL);
	tmpTimeStart_init=tim.tv_sec+(tim.tv_usec/1000000.0);
#endif


	/* sets some important global variables */
	set_global_vars1();

	/* captures i/o signals to close */
	trap_sigs();

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	

	/* parse input parameter file, and read input data */
	parser(argc, argv, rng);

	/* MPI: These variables are used for storing data partitioning related information in the parallel version. These arrays store the start and end indices in the global array that each processor is responsible for processing. In the serial version, these are used to mimic the parallel version to obtain comparable results. */
	Start = (int *) calloc(procs, sizeof(int));
	End = (int *) calloc(procs, sizeof(int));

	if(USE_TT_FILE)
		load_tidal_tensor();

#ifndef USE_MPI
/*
MPI: These are arrays used for mimicking the parallel version. In the parallel version new stars created are stored at the end of the local star array of each processor. In the serial version, the same is done i.e. they are placed at the end of the star array beyond the sentinel (which is a nullified star indicating the end of old stars and beginning of newly created ones). However, in order for the serial version to mimic the parallel, it is essential to know which node would have created these stars in a corresponding parallel run so as to draw a random number from the appropriate stream. Although this is trivial for the old stars since it is a just function of the index of the star, for newly created stars it is not, since they are stored at the end of the array and are mixed up.

We use these two arrays to store the number of stars created by each node during a timestep. Stars can be created in two routines - dynamics, or stellar evolution, so we have an array for each of the two. These arrays are as large as the number of processors and their elements are used to store the number or stars that would have been created by each processor in a corresponding parallel run. With this supplementary information, we can deduce which proc would be having a given star, in the parallel version, and draw a random number from an appropriate stream.
*/
	created_star_dyn_node = (int *) calloc(procs, sizeof(int));
	created_star_se_node = (int *) calloc(procs, sizeof(int));
	/* MPI: Used to mimic the DMse quantity of the parallel version. */
	DMse_mimic = (double *) calloc(procs, sizeof(double));
#endif

	if(RESTART_TCOUNT == 0){
		/* MPI: Populating Start and End arrays based on data partitioning scheme */
		findLimits( cfd.NOBJ, MIN_CHUNK_SIZE );

		/* copy the input data which has been read from the input file into the cfd struct in parser() into the star and binary data structures */
		get_star_data(argc, argv, rng);

		/* Initializes rng states for parallel version and serial mimic */
		set_rng_states();
	} else {
		/*Load the restart file directly into star and binary array*/
		load_restart_file();

		/*Find the limits on each MPI process, but using the actual star numbers*/
		findLimits( clus.N_MAX, MIN_CHUNK_SIZE );

		/*Should already be sorted, but new stars need to be accounted for*/
		qsorts_new();

		/*Communicate that information to the global arrays*/
		post_sort_comm();

	}

//Probably not needed anymore
//	N_b_OLD = N_b;
//	N_b_NEW = N_b;

	if(USE_TT_FILE){
		load_tidal_tensor();
		orbit_r = compute_tidal_boundary();
	} else {
		orbit_r = R_MAX;
	}

	/* compute the potential */
	calc_potential_new();

	calc_sigma_new();

	/* calculate central quantities */
	central_calculate();

	//MPI: Setting this because in the MPI version only _r is set, and will create problem while sorting.
#ifndef USE_MPI
	star[clus.N_MAX + 1].r = SF_INFINITY;
	star[clus.N_MAX + 1].phi = 0.0;
#endif

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

	/* Calculates some global binary variables - total binary mass,and E. */
	bin_vars_calculate();

	/*
		Skipping search grid for MPI
		if (SEARCH_GRID) 
		search_grid_update(r_grid);
	 */	

	/* compute energy initially */
	star[0].E = star[0].J = 0.0;

    set_energy_vars();

	compute_energy_new();


	/* If we don't set it here, new stars created by breaking binaries (BSE) will
	 * end up in the wrong place */
#ifdef USE_MPI
    clus.N_MAX_NEW = mpiEnd-mpiBegin+1;
#else
    clus.N_MAX_NEW = clus.N_MAX;
#endif

	comp_mass_percent();

	comp_multi_mass_percent();

	/* initialize stellar evolution things */
	DMse = 0.0;


#ifndef USE_MPI
	for(i=0; i<procs; i++)
		DMse_mimic[i] = 0.0;
#endif


	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);
	if (STELLAR_EVOLUTION > 0) {
	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);
		if(RESTART_TCOUNT == 0 )
			stellar_evolution_init();
		else
			restart_stellar_evolution();
	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);
	}
	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);


#ifndef USE_MPI
	for(i=0; i<procs; i++)
		created_star_se_node[i] = 0;
#endif

	//OPT: M_b E_b calculated twice, also in bin_vars_calculate? Check for redundancy.
	update_vars();

	times(&tmsbufref);

	//OPT: Check for redundancy. Ask Stefan
	/* calculate central quantities */
	central_calculate();

	if (WRITE_EXTRA_CORE_INFO) {
		no_remnants= no_remnants_core(6);
	}

	calc_clusdyn_new();

	if(RESTART_TCOUNT == 0){
	/* print out binary properties to a file */
	//print_initial_binaries();
	
		/* Printing Results for initial model */
		print_results();

		/* print handy script for converting output files to physical units */
		print_conversion_script();

	#ifdef USE_CUDA
		cuInitialize();
	#endif

	#ifdef USE_MPI
		tmpTimeStart = timeStartSimple();
		//MPI: Apparently there are changes to the masses of the stars in stellar_evolution_init, so here before we start the timestep loop, we synchronize the mass array across all processors.
		if (STELLAR_EVOLUTION > 0) {
			mpiFindDispAndLenCustom( clus.N_MAX, MIN_CHUNK_SIZE, mpiDisp, mpiLen );

			for(i=0;i<procs;i++)
				mpiLen[i] *= sizeof(double); 

			//MPI: Can be replaced by AllGatherv similar to the way it's done after sorting in post_sort_comm.
			MPI_Status stat;
			if(myid!=0)
				MPI_Send(&star_m[mpiDisp[myid]], mpiLen[myid], MPI_BYTE, 0, 0, MPI_COMM_WORLD);
			else
				for(i=1;i<procs;i++)
					MPI_Recv(&star_m[mpiDisp[i]], mpiLen[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &stat);
			MPI_Bcast(star_m, clus.N_MAX+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		timeEndSimple(tmpTimeStart, &t_comm);
	}
	#endif

#ifdef USE_MPI
	t_init = MPI_Wtime() - tmpTimeStart_init;
	tmpTimeStart_full = MPI_Wtime();
#else
	gettimeofday(&tim, NULL);
	t_init += tim.tv_sec+(tim.tv_usec/1000000.0) - tmpTimeStart_init;
	tmpTimeStart_full=tim.tv_sec+(tim.tv_usec/1000000.0);
#endif

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop() == 0) 
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
//MPI: These are some global variables that are update at various places during the timestep, and towards the end need to be summed up across all processors. So, we store the values of these variables from the previous timestep into corresponding _old variables, and reset the actual variables to zero. At the end of the timestep, we cumulate/reduce the actual variables across processors and finally add them to the _old value i.e. total value of the variable from the previous timestep to obtain the updated values for these variables.
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
		timeEndSimple(tmpTimeStart, &t_oth);

		/* calculate central quantities */
		tmpTimeStart = timeStartSimple();
		central_calculate();
		timeEndSimple(tmpTimeStart, &t_cen_calc);

		/* calculate timestep */
		tmpTimeStart = timeStartSimple();
		calc_timestep(rng);
		timeEndSimple(tmpTimeStart, &t_timestep);

		/* If we're using a tidal tensor, advance to the new orbit of 
		 * the cluster */
		if(USE_TT_FILE)
			orbit_r = compute_tidal_boundary();

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		tmpTimeStart = timeStartSimple();
		//MPI: Note that N_MAX_NEW is the number of local stars (including newly created ones) in the parallel version, whereas in the serial it refers to the total number of stars including newly created ones. N_MAX is the total number of stars excluding newly created ones in both serial and parallel versions.
#ifdef USE_MPI
		clus.N_MAX_NEW = mpiEnd-mpiBegin+1;
#else
		clus.N_MAX_NEW = clus.N_MAX;
#endif
		timeEndSimple(tmpTimeStart, &t_oth);

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		tmpTimeStart = timeStartSimple();
		if (PERTURB > 0)
			dynamics_apply(Dt, rng);
		timeEndSimple(tmpTimeStart, &t_dyn);

		tmpTimeStart = timeStartSimple();
	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);
		if (STELLAR_EVOLUTION > 0)
			do_stellar_evolution(rng);
		timeEndSimple(tmpTimeStart, &t_se);
	printf("GOT HERE=%d in %s on proc=%d\n",__LINE__,__FILE__,myid);

		tmpTimeStart = timeStartSimple();
		Prev_Dt = Dt;

		/* evolve stars up to new time */
		DMse = 0.0;
#ifndef USE_MPI
		for(i=0; i<procs; i++)
			DMse_mimic[i] = 0.0;
#endif

#ifndef USE_MPI
		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
			will disappear! */
		/* This really has to come after SE otherwise merger products will disappear. */
		clus.N_MAX_NEW++;
#endif

		/* some numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		energy_conservation1();

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		toy_rejuvenation();
		timeEndSimple(tmpTimeStart, &t_oth);

		/* this calls get_positions() */
		tmpTimeStart = timeStartSimple();
		new_orbits_calculate();
		timeEndSimple(tmpTimeStart, &t_orb);

		/* more numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		tmpTimeStart = timeStartSimple();
		energy_conservation2();
		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		if(TIDALLY_STRIP_STARS > 0)
			tidally_strip_stars();
		timeEndSimple(tmpTimeStart, &t_tid_str);

		/* Compute Intermediate Energies of stars. 
		 * Also transfers new positions and velocities from srnew[], 
		 * svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
		 */
		tmpTimeStart = timeStartSimple();
		ComputeIntermediateEnergy();
		timeEndSimple(tmpTimeStart, &t_oth);

		/* sort stars by radial positions */
		tmpTimeStart = timeStartSimple();
		qsorts_new();
		timeEndSimple(tmpTimeStart, &t_sort);

		tmpTimeStart = timeStartSimple();
		post_sort_comm();
		timeEndSimple(tmpTimeStart, &t_postsort_comm);

		/* compute the potential */
		tmpTimeStart = timeStartSimple();
		calc_potential_new();
		timeEndSimple(tmpTimeStart, &t_pot_cal);

		//MPI: Since N_MAX is updated during potential calculation, we update Start and End values for each processor based on new number of stars.
		tmpTimeStart = timeStartSimple();
		findLimits( clus.N_MAX, MIN_CHUNK_SIZE );
		timeEndSimple(tmpTimeStart, &t_oth);

		tmpTimeStart = timeStartSimple();
		energy_conservation3();
		timeEndSimple(tmpTimeStart, &t_ener_con3);

		//commenting out for MPI
		/*
			if (SEARCH_GRID)
			search_grid_update(r_grid);
		 */
	
		tmpTimeStart = timeStartSimple();
		comp_mass_percent();
		timeEndSimple(tmpTimeStart, &t_calc_io_vars1);

		tmpTimeStart = timeStartSimple();
		comp_multi_mass_percent();
		timeEndSimple(tmpTimeStart, &t_calc_io_vars2);

		tmpTimeStart = timeStartSimple();
		compute_energy_new();
		timeEndSimple(tmpTimeStart, &t_comp_ener);

		tmpTimeStart = timeStartSimple();
		reset_interaction_flags();
		timeEndSimple(tmpTimeStart, &t_oth);

		/* update variables, then print */
		tmpTimeStart = timeStartSimple();
		update_vars();
		timeEndSimple(tmpTimeStart, &t_upd_vars);

		tmpTimeStart = timeStartSimple();
		calc_clusdyn_new();
		timeEndSimple(tmpTimeStart, &t_oth);

        if (WRITE_EXTRA_CORE_INFO) {
            no_remnants= no_remnants_core(6);
        }


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


		tmpTimeStart = timeStartSimple();
		print_results();

		print_snapshot_windows();

		tcount++;


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
		timeEndSimple(tmpTimeStart, &t_io);

		tmpTimeStart = timeStartSimple();
		//Print out timer file
		if(TIMER)
		{
			rootfprintf(timerfile, "%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", tcount, t_cen_calc, t_timestep, t_dyn, t_se, t_orb, t_tid_str, t_sort, t_postsort_comm, t_pot_cal, t_ener_con3, t_calc_io_vars1, t_calc_io_vars2, t_comp_ener, t_upd_vars, t_io, t_io_ignore, t_oth, t_sort_lsort1, t_sort_splitters, t_sort_a2a, t_sort_lsort2, t_sort_oth, t_sort_lb, t_sort_only);
		}
		timeEndSimple(tmpTimeStart, &t_io_ignore);

	    update_tspent(tmsbufref);

		if(CheckCheckpoint())
			save_restart_file();

	} /* End time step iteration loop */

	save_restart_file();

	times(&tmsbuf);
#ifdef USE_MPI
	t_full = MPI_Wtime() - tmpTimeStart_full;
#else
	gettimeofday(&tim, NULL);
	t_full += tim.tv_sec+(tim.tv_usec/1000000.0) - tmpTimeStart_full;
#endif

	//Print overall timings to stdout
	rootprintf("******************************************************************************\n");
	rootprintf("Time for Initialization: %.4lf\n", t_init);
	rootprintf("Total Time Taken: %.4lf\n", t_full);
	rootprintf("Time for Communication: %.4lf\n", t_comm);
	rootprintf("******************************************************************************\n");

	//Print overall timings to file
	if(TIMER)
	{
		rootfprintf(timerfile, "******************************************************************************\n");
		rootfprintf(timerfile, "Time for Initialization: %.4lf\n", t_init);
		rootfprintf(timerfile, "Total Time Taken: %.4lf\n", t_full);
		rootfprintf(timerfile, "Time for Communication: %.4lf\n", t_comm);
		rootfprintf(timerfile, "******************************************************************************\n");
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
	free_arrays();

#ifdef USE_MPI
	free(mpiDisp);
	free(mpiLen);
	free(curr_st);
//Probably not needed anymore
//	free(binary_buf);
//	free(num_bin_buf);
#else
	free(st);
#endif

	free(Start);
	free(End);

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

#ifdef USE_MPI
	int j;
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
			for( j = 1; j <= mpiEnd-mpiBegin+1; j++ )
				fprintf(ftest, "%d\t%.18g\n", get_global_idx(j), sigma_array.sigma[j]);
			fclose(ftest);
		}
	}
	if(myid==0)
	{
		char process_str[30];
		sprintf(process_str, "./process.sh %d", procs);
		system(process_str);
	}
#else
	int j;
	strcpy(tempstr, "test_out_ser.dat");
	ftest = fopen( tempstr, "w" );
	for( j = 1; j <= clus.N_MAX; j++ )
	fprintf(ftest, "%d\t%.18g\n", j, sigma_array.sigma[j]);
	fclose(ftest);
#endif

*/
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

//Testing MPI-IO
/*
#ifdef USE_MPI
    char mpi_iotest_wrbuf[1000000];
    char mpi_iotest_buf[100];
    int mpi_iotest_len=0, foutoffset=0, cum_offset=0, tot_offset=0;
    MPI_File fh;
    MPI_Status mpistat;

    //MPI3: Opening and closing file once in case it exists to delete it. Then opening again. Should find a better solution later.
    MPI_File_delete ("mpiiotest.log", MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, "mpiiotest.log", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//---------------------
    mpi_iotest_wrbuf[0] = '\0';
    mpi_iotest_len=0;
    foutoffset=0;

    parafprintf(iotest, "My name is %d\n", myid);
    //   sprintf(mpi_iotest_buf, "My name is %d\n", myid);
    //   strcat(mpi_iotest_wrbuf, mpi_iotest_buf);
    //   mpi_iotest_len += strlen(mpi_iotest_buf);
    if(myid==2)
        parafprintf(iotest, "My name is not at all fugu\n");
    else
        parafprintf(iotest, "My name is also fugu\n");
    //   strcat(mpi_iotest_wrbuf, mpi_iotest_buf);
    //   mpi_iotest_len += strlen(mpi_iotest_buf);
    MPI_Exscan(&mpi_iotest_len, &foutoffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    foutoffset += cum_offset;

    MPI_File_write_at_all(fh, foutoffset, mpi_iotest_wrbuf, mpi_iotest_len, MPI_CHAR, &mpistat);
    MPI_Allreduce (&mpi_iotest_len, &tot_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    cum_offset += tot_offset;

    //printf("----->%d strln=%d foutoffset=%d cum_off=%d\n", myid, lenstr, foutoffset, cum_offset);
    MPI_File_close(&fh);
#endif
*/
