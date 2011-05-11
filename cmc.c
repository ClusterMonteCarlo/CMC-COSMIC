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
//#include "cmc_mpi.h"

//#define USE_MPI

#ifdef USE_CUDA
#include "cuda/cmc_cuda.h"
#endif

int main(int argc, char *argv[])
{
#ifdef USE_MPI
	//int myid, procs; //Declared in cmc.h to make it available across all files
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
	FILE *ftest1;
	FILE *ftest2;
	FILE *ftime;
	int si;

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

	//Currently doing on all nodes to avoid broadcasting of global variables. Later, might have to be split into 2 or more functions, and store all global variables into a structure for easy broadcast.
	/* parse input */
	parser(argc, argv, rng); //to do parallel i/o


#ifdef USE_MPI	
	//MPI_Barrier( MPI_COMM_WORLD ) ;
	//printf("\nI am process ID: %d\n", myid);

	//These have to be declared here because N_STAR_DIM is initialized only in parser().
	double test1[N_STAR_DIM];
	double test2[N_STAR_DIM];	
	double test3[N_STAR_DIM];	
	int mpiBegin, mpiEnd;
   double starttime, endtime, temptime; 
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
	/*	Commenting out for MPI
		q_root = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
	 */

	/* Commenting out for MPI
#ifdef DEBUGGING
	// create a new hash table for the star ids 
	star_ids= g_hash_table_new(g_int_hash, g_int_equal);
	// load a trace list 
	//load_id_table(star_ids, "trace_list");
#endif
	 */

	/* Set up initial conditions */
	//can be done on all nodes? which means no need to broadcast star structure.
	/*
#ifdef USE_MPI
if(myid == 0)
#endif
	 */
	//This fn populates the star array with the the data obtained after parsing
	load_fits_file_data(); 

	/*
#ifdef USE_MPI
if(myid == 0)
#endif
	 */
	//Step 2: do only on root node....not
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

	//Distrubute star array. Broadcast r and phi arrays.
	/*
#ifdef MPI
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
	 */

/*
#ifdef USE_MPI
	long mpiBeginN_STAR, mpiEndN_STAR;
	MPI_Barrier(MPI_COMM_WORLD);
	mpiFindIndices( clus.N_STAR, &mpiBeginN_STAR, &mpiEndN_STAR );
	printf("\n\nProc %d:\tLow=%ld\tHigh=%ld\n\n", myid , mpiBeginN_STAR, mpiEndN_STAR);
#endif
*/
	//Step 2:currently do on all nodes. But all procs will use same random numbers. So, comparing with host code might be a problem
	reset_rng_t113(IDUM);

	/* binary remainders */
	clus.N_MAX = clus.N_STAR;
	N_b = clus.N_BINARY;

#ifdef USE_MPI
	/*MPI2: Initializing and extracting global arrays that will be needed by all processors.*/
	//MPI2: Tested
	star_r = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_m = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_phi = (double *) malloc(N_STAR_DIM * sizeof(double));

	if(myid==0) {
		for(i=0; i<=N_STAR_DIM; i++) {
			star_r[i] = star[i].r;
			star_m[i] = star[i].m;
		}
	}
	MPI_Bcast(star_m, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(star_r, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
/*
if(myid==0){
	ftest1 = fopen("test_mpi_g0.dat","w");
	for (si = 0; si <= N_STAR_DIM; si++)
		fprintf(ftest1, "%ld\t%g\t%g\n", si, star[si].m, star_m[si] );
	fclose(ftest1);
}
if(myid==1){
	ftest1 = fopen("test_mpi_g1.dat","w");
	for (si = 0; si <= N_STAR_DIM; si++)
		fprintf(ftest1, "%g\n", star_m[si]);
	fclose(ftest1);
}
*/
#endif

/*
#ifdef USE_MPI
	long *mpiDisp, *mpiLen;
	mpiDisp = (long *) malloc(procs * sizeof(long));
	mpiLen = (long *) malloc(procs * sizeof(long));
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0)
	{
		mpiFindDispAndLen( clus.N_STAR, mpiDisp, mpiLen );
		for(i=0;i<procs;i++)
			printf("\n\nDisp=%ld\tLen=%ld\n\n", mpiDisp[i], mpiLen[i]);
	}
	MPI_Finalize();
	return(1);
#endif
*/

//MPI2: mpi_calc_sigma_r(): Tests for precision: Biggest errors are ~ 1e-13. Might be because of catastrphic cancellation/rounding errors?. This might be due to the already imprecise but faster serial code . In a sense, the MPI version might be more precise, because for every processor, actual average is performed for the 1st local star which means errors dont carry over from the previous set of stars which are handled by another processor.
//	double testMpi[clus.N_MAX];

#ifdef USE_MPI
	int *mpiDisp, *mpiLen;
	mpiDisp = (int *) malloc(procs * sizeof(int));
	mpiLen = (int *) malloc(procs * sizeof(int));

	ftime = fopen("test_mpi_time.dat","w");
	if(myid==0)
   fprintf(ftime, "function\t\t\t\t\tmpitime\t\tserialtime\n"); 
   starttime = MPI_Wtime(); 
	mpi_calc_sigma_r();

	//if(myid==0)
	mpiFindDispAndLen( clus.N_STAR, mpiDisp, mpiLen );
/*
	for(i=0;i<procs;i++)
	{
		if(myid==0)
		printf("\n\nProc=%d\tDisp=%ld\tLen=%ld\n\n",myid, mpiDisp[i], mpiLen[i]);
	}
*/
	MPI_Gatherv( &sigma_array.sigma[mpiDisp[myid]], mpiLen[myid], MPI_DOUBLE, &sigma_array.sigma[0], mpiLen, mpiDisp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//Testing mpi_calc_sigma_r()
/*
if(myid==0){
	ftest1 = fopen("test_mpi_sigma1.dat","w");
	for (si = 1; si <= clus.N_MAX; si++){
		testMpi[si] = sigma_array.sigma[si];
		fprintf(ftest1, "%ld\t%g\n", si, sigma_array.sigma[si]);
	}
	fclose(ftest1);
}
*/
   endtime = MPI_Wtime(); 
	temptime = endtime-starttime;
	if(myid==0)
   fprintf(ftime, "calc_sigma_r\t\t\t%g",temptime); 
//#else
   starttime   = MPI_Wtime(); 
	//for first step, do on all nodes. later will need scattering of particles and also some ghost particles data.
	if(myid==0)
	calc_sigma_r(); //requires data from some neighbouring particles. must be handled during data read. look inside function for more comments
   endtime   = MPI_Wtime(); 
	if(myid==0)
   fprintf(ftime, "\t%g\t%g\n",endtime-starttime,(endtime-starttime)/temptime); 
#endif
//Testing mpi_calc_sigma_r()
/*
#ifdef USE_MPI
if(myid==0){
	ftest1 = fopen("test_mpi_sigma.dat","w");
	for (si = 1; si <= clus.N_MAX; si++)
		fprintf(ftest1, "%ld\t%g\t%g\t%g\n", si, sigma_array.sigma[si], testMpi[si],(testMpi[si]-sigma_array.sigma[si])/sigma_array.sigma[si] );
	fclose(ftest1);
}
#endif
*/


	//Step 2: stay on root node
	//MPI2: Split into 2 functions: part 1 can be parallelized after making m and r arrays global. Part 2 has to be done on root node.
#ifdef USE_MPI
   starttime = MPI_Wtime(); 
	//MPI2: Tested! Errors of 1e-14.
	mpi_central_calculate1();
	mpi_central_calculate2();
/*
	central_t central_test;
	central_test = central;
*/
   endtime = MPI_Wtime(); 
	temptime = endtime-starttime;
	if(myid==0)
   fprintf(ftime, "central_calculate\t\t%g",temptime); 
//#else
   starttime   = MPI_Wtime(); 
	if(myid==0)
	central_calculate();
   endtime   = MPI_Wtime(); 
	if(myid==0)
   fprintf(ftime, "\t%g\t%g\n",endtime-starttime,(endtime-starttime)/temptime); 
#endif
/*
if(myid==1)
		printf( "----------------------\n\n\t%g\t%g\t%g\n\n--------------------", central.a_ave, central_test.a_ave, (central.a_ave-central_test.a_ave)/central.a_ave);
*/

/*
if(myid==1)
printf( "----------------------\n\n%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n\n--------------------", central.N_sin, central.N_bin, central.v_sin_rms, central.w2_ave, central.R2_ave, central.mR_ave, central.a_ave, central.a2_ave, central.ma_ave);
*/


	//MPI1:broadcast central structure
/*
#ifdef USE_MPI
	MPI_Bcast(&central, sizeof(central_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&v_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_core_nb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&core_radius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Trc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_single, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rho_core_bin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
*/
	
	/*
		if (WRITE_EXTRA_CORE_INFO) {
		no_remnants= no_remnants_core(6);
		}
	 */


	//MPI2: Binaries; Ignore.
	M_b = 0.0;
	E_b = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		j = star[i].binind;
		if (j && binary[j].inuse) {
			M_b += star[i].m;
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}

/*
	//Step 2: MPI reduce
	//will have to do a reduce for M_b abd E_b later on. for the 1st step ignore.
#ifdef USE_MPI
	double M_b_local, E_b_local;
	double M_b_test, E_b_test;
	M_b_local = 0.0;
	E_b_local = 0.0;
	M_b_test = 0.0;
	E_b_test = 0.0;

	//int indLowNSTAR = myid * clus.N_STAR / procs + 1;
	//int indHighNSTAR = (myid+1 == procs) ? clus.N_STAR : ( (myid + 1) * clus.N_STAR / procs ) ;
	//printf("NoFunc\tProc %d:\tLow=%d\tHigh=%d\n",myid, indLowNSTAR, indHighNSTAR);

	int mpiBegin, mpiEnd;
	mpiFindIndices( clus.N_STAR, &mpiBegin, &mpiEnd );
	//printf("Func\tProc %d:\tLow=%ld\tHigh=%ld\n",myid, mpiBegin, mpiEnd);

	for (i=mpiBegin; i<=mpiEnd; i++) {
		j = star[i].binind;
		//MPI Note: Binary index from other proc. used?
		if (j && binary[j].inuse) {
			M_b_local += star[i].m;
			E_b_local += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}
	//Later must check if the sums are done in order
	//printf("Proc %d:\tE_b_test=%g\tM_b_test=%g\n",myid, E_b_local, M_b_local);

	MPI_Reduce(&M_b_local, &M_b_test, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&E_b_local, &E_b_test, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	printf("\n\nM_b = %g\t M_b_test = %g\nE_b = %g\t E_b_test = %g\n", M_b, M_b_test, E_b, E_b_test);
#endif
*/

	/* print out binary properties to a file */
	//skipping outputs for initial steps of MPI
	//print_initial_binaries();

	orbit_r = R_MAX;

	//MPI2: Do on root node with global phi array.
	//MPI2: Tested. Relative errors are 0. Perfect :)
#ifdef USE_MPI
   starttime = MPI_Wtime(); 
if(myid==0) 
	mpi_potential_calculate();

	MPI_Bcast(star_phi, clus.N_MAX+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
/*
if(myid==1) {
	for(i=1; i<=clus.N_MAX; i++) {
		test1[i] = star_phi[i];
	}
}
*/
   endtime = MPI_Wtime(); 
	temptime = endtime-starttime;
	if(myid==0)
   fprintf(ftime, "potential_calculate\t%g",temptime); 
//#else
   starttime   = MPI_Wtime(); 
	if(myid==0)
	potential_calculate();
   endtime   = MPI_Wtime(); 
	if(myid==0)
   fprintf(ftime, "\t%g\t%g\n",endtime-starttime,(endtime-starttime)/temptime); 
/*
	ftest2 = fopen("test_mpi_phi2.dat","w");
	for(i=1; i<=clus.N_MAX; i++) {
			//fprintf( ftest2,  "%d\t%g\t%g\t%g\n",i, star[i].phi, test1[i], (star[i].phi-test1[i])/star[i].phi);
			fprintf( ftest2,  "%d\t%g\n",i, star[i].phi);
	}
	fclose(ftest2);
*/
#endif

/*
#ifdef USE_MPI
if(myid==1)
{
	ftest2 = fopen("test_mpi_phi3.dat","w");
	for(i=1; i<=clus.N_MAX; i++) {
			fprintf( ftest2,  "%d\t%g\t%g\t%g\n",i, star[i].phi, test1[i], (star[i].phi-test1[i])/star[i].phi);
			//fprintf( ftest2,  "%d\t%g\n",i, star[i].phi);
	}
	fclose(ftest2);
}
#endif
*/
	total_bisections= 0;

	/*
		Skipping search grid for MPI
		if (SEARCH_GRID) 
		search_grid_update(r_grid);
	 */	

	/* compute energy initially */
	//Step 2: do on root.
	//Bharath: doesnt have to be done on root since all procs have the star array, and they probably need the 0th star's values.
/*
#ifdef USE_MPI
	if (myid==0)
#endif
*/
		star[0].E = star[0].J = 0.0;

//MPI2: Tested! No errors.
#ifdef USE_MPI
   starttime = MPI_Wtime(); 
	mpi_ComputeEnergy();
/*
	if(myid==2)
		printf("\nParallel =%g\n", si, Etotal.tot);
*/
   endtime = MPI_Wtime(); 
	temptime = endtime-starttime;
	if(myid==0)
   fprintf(ftime, "ComputeEnergy\t\t\t%g",temptime); 
//#else
   starttime   = MPI_Wtime(); 
	//MPI2: see inline for comments
	if(myid==0)
	ComputeEnergy();
   endtime   = MPI_Wtime(); 
	if(myid==0)
   fprintf(ftime, "\t%g\t%g\n",endtime-starttime,(endtime-starttime)/temptime); 
/*
	if(myid==2)
		printf("\nSerial =%g\n", si, Etotal.tot);
*/	
#endif

	//MPI2: The following line has been moved after load_fits_file_data() - above
	//	star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;

	/* Noting the total initial energy, in order to set termination energy. */
	Etotal.ini = Etotal.tot;

	Etotal.New = 0.0;
	Eescaped = 0.0;
	Jescaped = 0.0;
	Eintescaped = 0.0;
	Ebescaped = 0.0;
	Eoops = 0.0;

#ifndef USE_MPI
	//MPI2: ignore
	comp_mass_percent(); //used for diagnostics. needs neighouring particles to compute cum. sum.
	//MPI2: ignore
	comp_multi_mass_percent();
#endif

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

	//root node	
/*
#ifdef USE_MPI
	if (myid==0)
#endif
*/
	times(&tmsbufref);

//MPI2: Is this broadcast needed?
/*
#ifdef USE_MPI
	MPI_Bcast(&tmsbufref, sizeof(struct tms), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
*/

	/* calculate central quantities */
#ifdef USE_MPI
	//MPI2: Tested! Errors of 1e-14.
	mpi_central_calculate1();
	mpi_central_calculate2();
#else
	central_calculate();
#endif

	/* can skip for MPI
		if (WRITE_EXTRA_CORE_INFO) {
		no_remnants= no_remnants_core(6);
		}
	 */

	/* calculate dynamical quantities */

	//Step 2: done on root node
	//MPI2: Tested! No Errors.
#ifdef USE_MPI
   starttime = MPI_Wtime(); 
	if (myid==0)
		mpi_clusdyn_calculate(); //parallel reduction

	//MPI2: broadcast clusdyn struct	
	MPI_Bcast(&clusdyn, sizeof(clusdyn_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);

   endtime = MPI_Wtime(); 
	temptime = endtime-starttime;
	if(myid==0)
   fprintf(ftime, "clusdyn_calculate\t\t%g",temptime); 
//#else
   starttime   = MPI_Wtime(); 
	if(myid==0)
	clusdyn_calculate(); //parallel reduction
   endtime   = MPI_Wtime(); 
	if(myid==0)
   fprintf(ftime, "\t%g\t%g\n",endtime-starttime,(endtime-starttime)/temptime); 
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

	/*******          Starting evolution               ************/
	/******* This is the main loop in the program *****************/
	while (CheckStop(tmsbufref) == 0) {

		/* calculate central quantities */
#ifdef USE_MPI
		//MPI2: Tested! Errors of 1e-14.
		mpi_central_calculate1();
		mpi_central_calculate2();
#else
		central_calculate();
#endif

		/* Get new time step */
		//for the next step, only simul_relax needs to be done in parallel, and DTrel needs to be broadcasted. rest can be done on each node.
		//Step 2: Do simul_relax() on all procs and others only on root. then broadcast Dt
		Dt = GetTimeStep(rng); //reduction again. Timestep needs to be communicated to all procs.

//printf("\ntimestep=%g\n", Dt);

		/*dprintf ("before tide: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		  star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
		  star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		/* if tidal mass loss in previous time step is > 5% reduce PREVIOUS timestep by 20% */
		//MPI2: root node
#ifdef USE_MPI
		if(myid==0) 
#endif
		{
			if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
				diaprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
				Dt = Prev_Dt * 0.8;
			} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
				diaprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
				Dt = Prev_Dt * 1.1;
			}
		}

		//MPI2: broadcast Dt		
#ifdef USE_MPI
		MPI_Bcast(&Dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//printf("proc=%d\ttimestep=%g\n", myid, Dt);
#endif

		/*dprintf ("after tide before dynamics: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		  star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
		  star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		TotalTime += Dt;

		/* set N_MAX_NEW here since if PERTURB=0 it will not be set below in perturb_stars() */
		clus.N_MAX_NEW = clus.N_MAX;

		/* Perturb velocities of all N_MAX stars. 
		 * Using sr[], sv[], get NEW E, J for all stars */

		//MPI2: Tested for outputs: vr, vt, E and J. Tests performed with same seed for rng of all procs. Check done only for proc 0's values as others cant be tested due to rng. Must test after rng is replaced.
		if (PERTURB > 0) {
			dynamics_apply(Dt, rng);
		}

/*
#ifdef USE_MPI
if(myid==1)
#endif
{
	ftest2 = fopen("test_mpi_dyn1.dat","w");
	for(i=1; i<=clus.N_MAX; i++) {
			fprintf( ftest2,  "%d\t%g\n",i, star[i].J);
	}
	fclose(ftest2);
}
*/

		/*dprintf ("after dynamics before SE: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		  star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
		  star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		/* if N_MAX_NEW is not incremented here, then stars created using create_star()
			will disappear! */
		clus.N_MAX_NEW++;

		/* evolve stars up to new time */
		DMse = 0.0;

		//MPI2: Tested for outputs: rad, m E. Check if rng is used at all. Testing done only for proc 0.
		if (STELLAR_EVOLUTION > 0) {
			do_stellar_evolution(rng);
		}
/*
#ifdef USE_MPI
if(myid==0)
#endif
{
	ftest2 = fopen("test_mpi_se0.dat","w");
	for(i=1; i<=clus.N_MAX; i++) {
			fprintf( ftest2,  "%d\t%g\n",i, star[i].rad);
	}
	fclose(ftest2);
}
*/

		/*dprintf ("after SE: id=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
		  star[787].id,star[787].se_k,star[787].se_mass,star[787].se_mt,star[787].se_radius,star[787].se_lum,star[787].se_mc,star[787].se_rc,
		  star[787].se_menv,star[787].se_renv,star[787].se_ospin,star[787].se_epoch,star[787].se_tms,star[787].se_tphys,star[787].phi, star[787].r);*/

		Prev_Dt = Dt;

		/* some numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
#ifdef USE_MPI 
		mpiFindIndices( clus.N_MAX_NEW, &mpiBegin, &mpiEnd );

		for (i=mpiBegin; i<=mpiEnd; i++) {
#else
		for (i = 1; i <= clus.N_MAX_NEW; i++) {
#endif
			/* saving velocities */
			star[i].vtold = star[i].vt;
			star[i].vrold = star[i].vr;

			/* the following will get updated after sorting and
			 * calling potential_calculate(), needs to be saved 
			 * now */  
#ifdef USE_MPI
			star[i].Uoldrold = star_phi[i] + MPI_PHI_S(star_r[i], i);
#else
			star[i].Uoldrold = star[i].phi + PHI_S(star[i].r, i);
#endif

			/* Unewrold will be calculated after 
			 * potential_calculate() using [].rOld
			 * Unewrnew is [].phi after potential_calculate() */
		}

		/*Sourav: checking all stars for their possible extinction from old age*/
		//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
		DMrejuv = 0.0;
		if (STAR_AGING_SCHEME > 0) {
#ifdef USE_MPI 
			//int i, mpiBegin, mpiEnd;
			mpiFindIndices( clus.N_MAX, &mpiBegin, &mpiEnd );

			for (i=mpiBegin; i<=mpiEnd; i++) {
#else
			for (i=1; i<=clus.N_MAX; i++){
#endif
				remove_old_star(TotalTime, i);
			}	
		}

		/* this calls get_positions() */
		tidally_strip_stars();
#ifdef USE_MPI
if(myid==0)
#endif
{
	ftest2 = fopen("test_mpi_tss0.dat","w");
	for(i=1; i<=clus.N_MAX; i++) {
			fprintf( ftest2,  "%d\t%g\n",i, star[i].r);
	}
	fclose(ftest2);
}


#ifdef USE_MPI 
		//int mpiBegin, mpiEnd;
		mpiFindIndices( clus.N_MAX_NEW, &mpiBegin, &mpiEnd );
		for (i=mpiBegin; i<=mpiEnd; i++) {
#else
		/* more numbers necessary to implement Stodolkiewicz's
		 * energy conservation scheme */
		for (i = 1; i <= clus.N_MAX_NEW; i++) {
#endif
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
		
			if(myid==0)
			{
			ftest1 = fopen("test_mpi1.dat","w");
			for (si = 0; si < clus.N_MAX_NEW; si++)
			fprintf(ftest1, "%ld\t%g\n", si, test1[si]);
			fclose(ftest1);
			}
		/*
			if(myid==4)
			{
			ftest2 = fopen("test_mpi2.dat","w");
			for (si = 0; si < N_STAR_DIM; si++)
			fprintf(ftest2, "%d\t%g\n", si, test1[si]);
			fclose(ftest2);
			}
		 */
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

} /* End WHILE (time step iteration loop) */

//root node?
/*
#ifdef USE_MPI
if (myid == 0)
#endif
*/
	times(&tmsbuf);
/*
#ifdef USE_MPI
	MPI_Bcast(&tmsbuf, sizeof(struct tms), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
*/
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


	//MPI1: compare results(E,J,vr,vt,r,phi,m) of all nodes for first step


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
