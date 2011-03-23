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
#include <mpi.h>

#define USE_MPI

#ifdef USE_CUDA
#include "cuda/cmc_cuda.h"
#endif

int main(int argc, char *argv[])
{
	#ifdef USE_MPI
	//int myid, procs; //Declared in cmc.h
	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD,&procs);
  	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	#endif

	#ifdef USE_MPI
	FILE *ftest1;
	FILE *ftest2;
	N_STAR_DIM = 2 + clus.N_STAR + 2 * clus.N_BINARY;
	N_BIN_DIM = 2 + clus.N_STAR / 2 + clus.N_BINARY;
	N_STAR_DIM = (long) floor(1.1 * ((double) N_STAR_DIM));
	N_BIN_DIM = (long) floor(1.1 * ((double) N_BIN_DIM));
	//double *r, *phi;
	//r = (double*)malloc(sizeof(double)*N_STAR_DIM);
	//phi = (double*)malloc(sizeof(double)*N_BIN_DIM);
	#endif

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
	trap_sigs(); // captures i/o signals to close

	/* initialize GSL RNG */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);

	#ifdef MPI
	if(myid == 0)
	#endif
	/* parse input */
		parser(argc, argv, rng); //to do parallel i/o

	#ifdef USE_MPI	
	MPI_Barrier( MPI_COMM_WORLD ) ;
	printf("\nI am process ID: %d\n", myid);

	double test1[N_STAR_DIM];
	double test2[N_STAR_DIM];	
	double test3[N_STAR_DIM];	
	int si;
	for (si = 0; si < N_STAR_DIM; si++)
	{	
		test1[si] = 0.0;
		test2[si] = 0.0;
		test3[si] = 0.0;
	}
	#endif

	/* print version information to log file */
	//commenting out print_version temporarily for basic mpi
	//print_version(logfile);

        /* initialize the Search_Grid r_grid */
	//If we use the GPU code, we dont need the SEARCH_GRID. So commenting it out
/*        if (SEARCH_GRID) { //Parallelize Search Grid
          r_grid= search_grid_initialize(SG_POWER_LAW_EXPONENT, \
              SG_MATCH_AT_FRACTION, SG_STARSPERBIN, SG_PARTICLE_FRACTION);
        };
*/
        /* initialize the root finder algorithm */
        q_root = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

#ifdef DEBUGGING
        /* create a new hash table for the star ids */
        star_ids= g_hash_table_new(g_int_hash, g_int_equal);
        /* load a trace list */
        //load_id_table(star_ids, "trace_list");
#endif
	/* Set up initial conditions */
	#ifdef USE_MPI
	if(myid == 0)
	#endif
		//This fn populates the star array with the the data obtained after parsing
		load_fits_file_data(); 


	#ifdef USE_MPI
	if(myid == 0)
	#endif
		//do only on root node
		star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;
	

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

//MPI_Broadcast(); 
/*broadcast the following:
star structure
binary structure
clus structure
*/
	//MPI_Bcast( star, clus.N_STAR, Startype, 0, MPI_COMM_WORLD);
	//MPI_Bcast( binary, clus.N_STAR, Binarytype, 0, MPI_COMM_WORLD);

	#ifdef USE_MPI
	MPI_Bcast(star, N_STAR_DIM * sizeof(star_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(binary,N_BIN_DIM * sizeof(binary_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&clus, sizeof(clus_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);

	//variables set in in load_fits_file_data()
	MPI_Bcast(&units, sizeof(units_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cenma, sizeof(struct CenMa), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&MEGA_YEAR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&SOLAR_MASS_DYN, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&madhoc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Mtotal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&newstarid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	#endif

	reset_rng_t113(IDUM);
	
	/* binary remainders */
	clus.N_MAX = clus.N_STAR;
	N_b = clus.N_BINARY;

	//for first step, do on all nodes. later will need scattering of particles and also some ghost particles data.
	calc_sigma_r(); //requires data from some neighbouring particles. must be handled during data read


	#ifdef USE_MPI
	if(myid == 0)
	#endif
	central_calculate();

	//broadcast central structure
	#ifdef USE_MPI
	MPI_Bcast(&central, sizeof(clus_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core_nb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&core_radius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Trc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_single, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_bin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

/*
        if (WRITE_EXTRA_CORE_INFO) {
          no_remnants= no_remnants_core(6);
        }
*/
	M_b = 0.0;
	E_b = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		j = star[i].binind;
		if (j && binary[j].inuse) {
			M_b += star[i].m;
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}
//will have to do a reduce for M_b abd E_b later on. for the 1st step ignore.

	
	/* print out binary properties to a file */
	//skipping outputs for initial steps of MPI
	//print_initial_binaries();
	
	orbit_r = R_MAX;

	//first step, ignore
	potential_calculate(); //requires the whole r and phi array. Communication needed. have to calculate cum.sum.

	total_bisections= 0;

/*
Skipping search grid for MPI
	if (SEARCH_GRID) 
		search_grid_update(r_grid);
*/	

	/* compute energy initially */
	star[0].E = star[0].J = 0.0;
	//needs cum.sum. skip for 1st step
	ComputeEnergy();

//has been moved after load_fits_file_data() - above
//	star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;

	/* Noting the total initial energy, in order to set termination energy. */
//do on root
	#ifdef USE_MPI
	if (myid==0)
	#endif
		Etotal.ini = Etotal.tot;
//broadcast Etotal	
	#ifdef USE_MPI
	MPI_Bcast(&Etotal, sizeof(Etotal_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif


	Etotal.New = 0.0;
	Eescaped = 0.0;
	Jescaped = 0.0;
	Eintescaped = 0.0;
	Ebescaped = 0.0;
	Eoops = 0.0;

	//ignore for first step	
	comp_mass_percent(); //used for diagnostics. needs neighouring particles to compute cum. sum.
	//ignore for first step	
	comp_multi_mass_percent();

	/* If we don't set it here, new stars created by breaking binaries (BSE) will
	 * end up in the wrong place */
	clus.N_MAX_NEW = clus.N_MAX;
	/* initialize stellar evolution things */
	DMse = 0.0;

	if (STELLAR_EVOLUTION > 0) {
		stellar_evolution_init(); //whole stellar evol. part does not need any data from other particles
	}
	
	update_vars(); //might need communication for bin. index array. needs cum.sum.

	//root node	
	#ifdef USE_MPI
	if (myid==0)
	#endif
		times(&tmsbufref);

	#ifdef USE_MPI
	MPI_Bcast(&tmsbufref, sizeof(struct tms), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif

	/* calculate central quantities */
	#ifdef USE_MPI
	if(myid == 0)
	#endif
	central_calculate();

	//broadcast central structure
	#ifdef USE_MPI
	MPI_Bcast(&central, sizeof(clus_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core_nb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&core_radius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Trc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_single, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_bin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	/* can skip for MPI
        if (WRITE_EXTRA_CORE_INFO) {
          no_remnants= no_remnants_core(6);
        }
	*/

	/* calculate dynamical quantities */

	//done on root node
	#ifdef USE_MPI
	if (myid==0)
	#endif
		clusdyn_calculate(); //parallel reduction
	//broadcast clusdyn struct	
	#ifdef USE_MPI
	MPI_Bcast(&clusdyn, sizeof(clusdyn_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif

	/* Printing Results for initial model */
	//skip outputs for MPI
	//print_results();

	/* print handy script for converting output files to physical units */
	//skip outputs for MPI
	//print_conversion_script();
	
#ifdef USE_CUDA
	cuInitialize();
#endif

	//printf("\n%d\tN_Core=%g\t\n", myid, N_core);

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) {
	//while (timecheck == 0) {
		/* calculate central quantities */


		#ifdef USE_MPI
		if(myid == 0)
		#endif
		central_calculate(); //parallel reduction

		//broadcast central struct or not? Ask Stefan.
		//MPI_Bcast(&central, sizeof(clus_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
		#ifdef USE_MPI
		MPI_Bcast(&central, sizeof(clus_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&N_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&N_core_nb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&rho_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&core_radius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Trc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&rho_core_single, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&rho_core_bin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif

		/* Get new time step */
		//for the next step, only simul_relax needs to be done in parallel, and DTrel needs to be broadcasted. rest can be done on each node.
		Dt = GetTimeStep(rng); //reduction again. Timestep needs to be communicated to all procs.

		/*dprintf ("before tide: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		     star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
	     		star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		/* if tidal mass loss in previous time step is > 5% reduce PREVIOUS timestep by 20% */
		//root node
		#ifdef USE_MPI
		if(myid==0) {
		#endif
		if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
			diaprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
			Dt = Prev_Dt * 0.8;
		} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
			diaprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
			Dt = Prev_Dt * 1.1;
		}
		#ifdef USE_MPI
		}
		#endif

		//broadcast Dt		
		#ifdef USE_MPI
		MPI_Bcast(&Dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif

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
		qsorts(star+1,clus.N_MAX_NEW); //parallel sort


		#ifdef USE_MPI
		for (si = 0; si < N_STAR_DIM; si++)
			test1[si] = star[si].r;
		if(myid==0){
		for (si = 0; si < N_STAR_DIM; si++)
			test2[si] = test1[si];
		}

		MPI_Bcast(&test2, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double check=5;
		double test_red=0;
		if(myid!=0){
			for (si = 0; si < N_STAR_DIM; si++)
			{
				test1[si] = test1[si]-test2[si];
				test_red += test1[si];
			}
		}
	
		MPI_Barrier(MPI_COMM_WORLD);	
		MPI_Reduce(&test_red, &check, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0)
		printf("\n\nCHECK = %.10g\n\n", check);
		#endif	

/*
		if(myid==0)
		{
			ftest1 = fopen("test_mpi1.dat","w");
			for (si = 0; si < N_STAR_DIM; si++)
				fprintf(ftest1, "%ld\t%g\n", si, test1[si]);
				//fprintf(ftest1, "%ld\t%g\t%g\n", si, test1[si], test2[si]);
			fclose(ftest1);
		}

		if(myid==1)
		{
			ftest2 = fopen("test_mpi2.dat","w");
			for (si = 0; si < N_STAR_DIM; si++)
				fprintf(ftest2, "%d\t%g\n", si, test1[si]);
				//fprintf(ftest2, "%d\t%g\t%g\n", si, test1[si], test3[si]);
			fclose(ftest2);
		}
		#endif
*/

/*
		#ifdef USE_MPI
		long si;
		for (si = 0; si < N_STAR_DIM; si++)
		{
			r[si] = star[si].r;
			phi[si] = star[si].phi;
		}
		#endif

		#ifdef USE_MPI
		//MPI: broadcast r array
		MPI_Bcast(r, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
*/
		potential_calculate();
/*		
		#ifdef USE_MPI
		//MPI: broadcast phi array	
		MPI_Bcast(phi, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
*/
		//commenting out for MPI
		/*
                if (SEARCH_GRID)
                  search_grid_update(r_grid);
		*/

		set_velocities3();

		comp_mass_percent();
		comp_multi_mass_percent();

		ComputeEnergy();

		/* reset interacted flag */
		for (i = 1; i <= clus.N_MAX; i++) {
			star[i].interacted = 0;
		}
		
		/* update variables, then print */
		update_vars();
		tcount++;
		
		/* calculate dynamical quantities */
		clusdyn_calculate();

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

	} /* End FOR (time step iteration loop) */

	//root node
	#ifdef USE_MPI
	if (myid == 0)
	#endif
		times(&tmsbuf);

	#ifdef USE_MPI
		MPI_Bcast(&tmsbuf, sizeof(struct tms), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif

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


//compare results(E,J,vr,vt,r,phi,m) of all nodes for first step

	
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

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return(0);
}
