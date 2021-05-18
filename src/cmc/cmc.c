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
#include "hdf5.h"
#include "hdf5_hl.h"

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

	//MPI: Some code from the main branch might have been removed in the MPI version. Please check.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* MPI: These variables are used for storing data partitioning related information in the parallel version. These are in particular useful for some MPI communication calls. */
	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));

	/* Starting timer to measure the overall time taken */
	tmpTimeStart_init = MPI_Wtime();

	/* sets some important global variables */
	set_global_vars1();

        /* sets some important global logging variables */
        //set_global_logfile_vars();

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

	/* Load the optional tidal tensor or dynamical friction files */
	if(USE_TT_FILE){
		load_tidal_tensor();
		orbit_r = compute_tidal_boundary();
	} else {
		orbit_r = R_MAX;
	}

	if(USE_DF_CUTOFF)
		load_dynamical_friction_data();

	/* compute the potential */
	calc_potential_new();

	calc_sigma_new();

	/* calculate central quantities */
	central_calculate();

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


	/* compute energy initially */
	star[0].E = star[0].J = 0.0;

	set_energy_vars();

	compute_energy_new();


	/* If we don't set it here, new stars created by breaking binaries (BSE) will
	 * end up in the wrong place */
	clus.N_MAX_NEW = mpiEnd-mpiBegin+1;

	comp_mass_percent();

	comp_multi_mass_percent();

	/* initialize stellar evolution things */
	DMse = 0.0;

	if (STELLAR_EVOLUTION > 0) {
		if(RESTART_TCOUNT == 0 )
			stellar_evolution_init();
		else
			restart_stellar_evolution();
	}

	update_vars();

	times(&tmsbufref);

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

	t_init = MPI_Wtime() - tmpTimeStart_init;
	tmpTimeStart_full = MPI_Wtime();

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop() == 0) 
	{

		tmpTimeStart = timeStartSimple();

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
		clus.N_MAX_NEW = mpiEnd-mpiBegin+1;
		timeEndSimple(tmpTimeStart, &t_oth);

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */
		tmpTimeStart = timeStartSimple();
		if (PERTURB > 0)
			dynamics_apply(Dt, rng);
		timeEndSimple(tmpTimeStart, &t_dyn);

		tmpTimeStart = timeStartSimple();
		if (STELLAR_EVOLUTION > 0)
			do_stellar_evolution(rng);
		timeEndSimple(tmpTimeStart, &t_se);

		tmpTimeStart = timeStartSimple();
		Prev_Dt = Dt;

		/* evolve stars up to new time */
		DMse = 0.0;

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

		tmpTimeStart = timeStartSimple();
		print_results();

		print_snapshot_windows();

		tcount++;

		timeEndSimple(tmpTimeStart, &t_io);

		tmpTimeStart = timeStartSimple();
		//Print out timer file
		if(TIMER)
		{
			rootfprintf(timerfile, "%ld\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", tcount, t_cen_calc, t_timestep, t_dyn, t_se, t_orb, t_tid_str, t_sort, t_postsort_comm, t_pot_cal, t_ener_con3, t_calc_io_vars1, t_calc_io_vars2, t_comp_ener, t_upd_vars, t_io, t_io_ignore, t_oth, t_sort_lsort1, t_sort_splitters, t_sort_a2a, t_sort_lsort2, t_sort_oth, t_sort_lb, t_sort_only);
		}
		timeEndSimple(tmpTimeStart, &t_io_ignore);

		update_tspent(tmsbufref);

		if(CheckCheckpoint())
			save_restart_file();

	} /* End time step iteration loop */

	save_restart_file();

	times(&tmsbuf);
	t_full = MPI_Wtime() - tmpTimeStart_full;

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

	/* free RNG */
	gsl_rng_free(rng);

#ifdef USE_CUDA
	cuCleanUp();
#endif

	/* flush buffers before returning */
	close_buffers();
	free_arrays();

	free(mpiDisp);
	free(mpiLen);
	free(curr_st);

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

	MPI_Finalize();

	return(0);
}
