/* -*- linux-c -*- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <fitsio.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#include "cmc_vars.h"

void printerror( int status) {
	/*****************************************************/
	/* Print out cfitsio error messages and exit program */
	/*****************************************************/
	if (status) {
		fits_report_error(stderr, status); /* print error report */
		exit_cleanly( status );    /* terminate the program,
				      returning error status */
	}
	return;
}

int move_to_hdu(fitsfile *fptr, char *sea_extname){
	/* a routine to find an HDU in a FITS file by name
	 * returns 0 on success, 1 on failure */
	int hdunum, hdutype, numhdu;
	int status;
	char extname[1024];

	status = 0; numhdu = 0;
	fits_get_num_hdus(fptr, &numhdu, &status);
	printerror(status);
	for(hdunum=1; hdunum<=numhdu; hdunum++){
		fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status);
		if (status==KEY_NO_EXIST) status = 0;
		printerror(status);
		if (strcmp(extname, sea_extname)==0) {
			return 0;
		}
	}
	return 1;
}

void read_fits_file_data(char *fits_file_name) {
	fitsfile *fptr;
	int status;

	status = 0;
	/* open file */
	fits_open_file(&fptr, fits_file_name, READONLY, &status);

	printerror(status);

	if (move_to_hdu(fptr, "CLUSTER_STARS") == 0) {
		read_fits_file_data_old(fptr);
	} else if (move_to_hdu(fptr, "META_INFO") == 0) {
		read_fits_file_data_new(fptr);
	} else {
		eprintf("Cannot read the input file\n");
		exit_cleanly(143);
	}
	
}

void read_fits_file_data_old(fitsfile *fptr) {
	int status, hdunum, hdutype, anynull;
	long frow, felem, nelem;
	unsigned long int i, NSTAR;
	double *dbl_arr;
	double avemass, totmass;

	status = 0;
	
	/* move to proper part of file */
	hdunum = 2; 		/* data table is in second HDU */
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);

	fits_read_key(fptr, TULONG, "NSTAR", &NSTAR, NULL, &status);
	printerror(status);
	if (clus.N_STAR != NSTAR){
		eprintf("There is an inconsistency regarding NSTAR\n");
		exit_cleanly(EXIT_FAILURE);
	}
	
	/* read columns into arrays and copy array into structure */
	dbl_arr = malloc((NSTAR+2)*sizeof(double));
	frow = 1; felem = 1; 
	nelem = NSTAR+2;
	fits_read_col(fptr, TDOUBLE, 1, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=NSTAR+1; i++){
		star[i].m = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 2, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=NSTAR+1; i++){
		star[i].r = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 3, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=NSTAR+1; i++){
		star[i].vr = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=NSTAR+1; i++){
		star[i].vt = dbl_arr[i];
	}
	printerror(status);

	/* some assignments so the code won't break */
	star[0].r = ZERO; 
	star[clus.N_STAR + 1].r = SF_INFINITY;
	Mtotal = 1.0;

	/* some normalization that set_imf() used to take care of */
	totmass = 0.0;
	for(i=1; i<=NSTAR; i++){
		totmass += star[i].m;
	}
	avemass = totmass/NSTAR;
	for(i=1; i<=NSTAR; i++){
		star[i].m /= avemass;
	}
	initial_total_mass = totmass;
	
	/* assign radii based on very simple mass--radius relationship
	   this will overridden by Belczynski's stellar evolution if STELLAR_EVOLUTION 
	   in the input file is non-zero */
	for (i=0; i<=NSTAR+1; i++) {
		star[i].rad = r_of_m(star[i].m);
	}

	/* zero each star's internal energy */
	for (i=0; i<=NSTAR+1; i++) {
		star[i].Eint = 0.0;
	}

	fits_close_file(fptr, &status);
	printerror(status);

	free(dbl_arr);
}

void read_fits_file_data_new(fitsfile *fptr) {
	int status, anynull;
	long frow, felem, nelem;
	unsigned long int i, N, NSTAR, NSTARN, NMAX, NMAXN;
	double *dbl_arr;
	long *lng_arr;
	double avemass, totmass;

	status = 0;
	
	/* NMAX and NSTAR are read from EXT: SS_DYN_PARAM */
	move_to_hdu(fptr, "SS_DYN_PARAM");
	fits_read_key(fptr, TULONG, "NMAX", &NMAX, NULL, &status);
	printerror(status);
	if (clus.N_MAX != NMAX){
		eprintf("There is an inconsistency regarding NMAX\n");
		eprintf("NMAX(read) = %ld, NMAX(set) = %ld\n",
				NMAX, clus.N_MAX);
		exit_cleanly(EXIT_FAILURE);
	}
	fits_read_key(fptr, TULONG, "NMAXN", &NMAXN, NULL, &status);
	printerror(status);
	if (clus.N_MAX_NEW != NMAXN){
		eprintf("There is an inconsistency regarding NMAX_NEW\n");
		eprintf("NMAXN(read) = %ld, NMAXN(set) = %ld\n",
				NMAXN, clus.N_MAX_NEW);
		exit_cleanly(EXIT_FAILURE);
	}
	fits_read_key(fptr, TULONG, "NSTAR", &NSTAR, NULL, &status);
	printerror(status);
	if (clus.N_STAR != NSTAR){
		eprintf("There is an inconsistency regarding NSTAR\n");
		exit_cleanly(EXIT_FAILURE);
	}
	fits_read_key(fptr, TULONG, "NSTARN", &NSTARN, NULL, &status);
	printerror(status);
	if (clus.N_STAR_NEW != NSTARN){
		eprintf("There is an inconsistency regarding NSTARNEW\n");
		exit_cleanly(EXIT_FAILURE);
	}
	
	/* read columns into arrays and copy array into structure */
	/* these columns are read from EXT: SS_DYN_PARAM */
	dbl_arr = malloc((NSTAR+2)*sizeof(double));
	lng_arr = malloc((NSTAR+2)*sizeof(long));
	frow = 1; felem = 1; 
	N = NSTARN;
	nelem = N+2;
	fits_read_col(fptr, TDOUBLE, 1, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].m = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 2, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].r = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 3, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].vr = dbl_arr[i];
	}
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].vt = dbl_arr[i];
	}

	fits_read_col(fptr, TDOUBLE, 5, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].r_peri = dbl_arr[i];
	}
	fits_read_col(fptr, TLONG, 6, frow, felem, nelem, NULL, lng_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].binind = lng_arr[i];
	}
	fits_read_col(fptr, TLONG, 7, frow, felem, nelem, NULL, lng_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].interacted = lng_arr[i];
	}
	/* 
        fits_read_col(fptr, TLONG, 8, frow, felem, nelem, NULL, lng_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].oldk = lng_arr[i];
	}
	*/
	fits_read_col(fptr, TDOUBLE, 9, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	for(i=0; i<=N+1; i++){
		star[i].phi = dbl_arr[i];
	}
	printerror(status);

	if(!ReadSnapshot){
		/* some assignments so the code won't break */
		Mtotal = 1.0;

		/* some normalization that set_imf() used to take care of */
		totmass = 0.0;
		for(i=1; i<=NSTAR; i++){
			totmass += star[i].m;
		}
		avemass = totmass/NSTAR;
		for(i=1; i<=NSTAR; i++){
			star[i].m /= avemass;
		}
		initial_total_mass = totmass;
	} else {
		comp_mass_percent();
	}

	/* assign radii based on very simple mass--radius relationship
	   this will overridden by Belczynski's stellar evolution if STELLAR_EVOLUTION 
	   in the input file is non-zero */
	for (i=0; i<=N+1; i++) {
		star[i].rad = r_of_m(star[i].m);
	}

	/* zero each star's internal energy */
	for (i=0; i<=NSTAR+1; i++) {
		star[i].Eint = 0.0;
	}

	fits_close_file(fptr, &status);
	printerror(status);

	free(dbl_arr); free(lng_arr);
}

void read_fits_file_parameters_old(fitsfile *fptr, gsl_rng *rng) {

	int status, hdunum, hdutype;

	status = 0;
	
	/* move to proper part of file */
	hdunum = 2; 		/* data table is in second HDU */
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);

	fits_read_key(fptr, TULONG, "NSTAR", &(clus.N_STAR), NULL, &status);
	printerror(status);
	
	fits_close_file(fptr, &status);
	printerror(status);
}

void read_fits_file_parameters_new(fitsfile *fptr, gsl_rng *rng) {

	long i;
	int status, anynull;
	int IS_SS;
	struct rng_t113_state rng_st;
	size_t rng_size;
	char *rng_st_ptr;
	int *int_arr;
	double *dbl_arr;
	long int *lng_arr;
	long frow, felem, nelem;
	int no_of_doub;
	double *dvar;

	status = 0;
	anynull = 0;
	rng_size = gsl_rng_size(rng);
	no_of_doub = 58;
	dvar = malloc(no_of_doub*sizeof(double));
	
	/* move to META_INFO to read ... meta info :-) */
	move_to_hdu(fptr, "META_INFO");
	IS_SS = 0;
	fits_read_key(fptr, TINT, "IS_SS", &IS_SS, NULL, &status);
	printerror(status);
	if(IS_SS){
		ReadSnapshot = 1;
	} 

	/* FIXME */
	/* need to add MMIN */
	/* FIXME */

	/* move to SF_RESTART_PARAMETERS */
	move_to_hdu(fptr, "SF_RESTART_PARAM");
	/* variables related to progress */
	fits_read_key(fptr, TDOUBLE, "Time", &TotalTime, NULL, &status);
	fits_read_key(fptr, TLONG, "Step", &tcount, NULL, &status);
	printerror(status);
	/* random number state */
	fits_read_key(fptr, TULONG, "RNG_Z1", &(rng_st.z1), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z2", &(rng_st.z2), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z3", &(rng_st.z3), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z4", &(rng_st.z4), NULL, &status);
	set_rng_t113(rng_st);
	printerror(status);
	/* Total Energy in various forms */
	fits_read_key(fptr, TDOUBLE, "Etotal", &(Etotal.tot), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "EtotalN", &(Etotal.New), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "EtotalI", &(Etotal.ini), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "KEtotal", &(Etotal.K), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "PEtotal", &(Etotal.P), NULL, &status);
	printerror(status);
	/* Sub timestep stuff */
	fits_read_key(fptr, TLONG, "S_NMAX", &(sub.N_MAX), NULL, &status);
	fits_read_key(fptr, TLONG, "S_CNT", &(sub.count), NULL, &status);
	fits_read_key(fptr, TLONG, "S_FACT", &(sub.FACTOR), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "S_Ttime", &(sub.totaltime), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "S_rmax", &(sub.rmax), NULL, &status);
	printerror(status);
	/* various N's of cluster */
	fits_read_key(fptr, TLONG, "NSTAR", &(clus.N_STAR), NULL, &status);
	fits_read_key(fptr, TLONG, "NSTARN", &(clus.N_STAR_NEW), NULL, &status);
	fits_read_key(fptr, TLONG, "NMAX", &(clus.N_MAX), NULL, &status);
	fits_read_key(fptr, TLONG, "NMAXN", &(clus.N_MAX_NEW), NULL, &status);
	fits_read_key(fptr, TLONG, "NBINARY", &(clus.N_BINARY), NULL, &status);
	printerror(status);
	/* Input file parameters */
	fits_read_key(fptr, TLONG, "NSTRDIM", &(N_STAR_DIM), NULL, &status);
	fits_read_key(fptr, TLONG, "MPCCNT", &(MASS_PC_COUNT), NULL, &status);
	/* fits_read_key(fptr, TLONG, "NTRY", &(N_TRY), NULL, &status); */
	fits_read_key(fptr, TLONG, "STEVOL", &(STELLAR_EVOLUTION), NULL, &status);
	fits_read_key(fptr, TLONG, "SSCOLL", &(SS_COLLISION), NULL, &status);
	fits_read_key(fptr, TLONG, "DUMPS", &(DUMPS), NULL, &status);
	fits_read_key(fptr, TLONG, "ECONS", &(E_CONS), NULL, &status);
	/* fits_read_key(fptr, TLONG, "DTMOD", &(DT_MODE), NULL, &status); */
	fits_read_key(fptr, TLONG, "MAXIND", &(MAX_INDEX), NULL, &status);
	fits_read_key(fptr, TLONG, "INDUNT", &(INDEX_UNIT), NULL, &status);
	fits_read_key(fptr, TLONG, "PERTURB", &(PERTURB), NULL, &status);
	fits_read_key(fptr, TLONG, "NUMMASS", &(NUM_MASS), NULL, &status);
	fits_read_key(fptr, TLONG, "TOTPAR", &(TOTAL_PARAMS), NULL, &status);
	printerror(status);
	fits_read_key(fptr, TLONG, "NMSRADB", &(NUM_MASS_RADII_BINS),
		     								NULL, &status);
	fits_read_key(fptr, TLONG, "NCNSTR", &(NUM_CENTRAL_STARS), NULL, &status);
	printerror(status);
	fits_read_key(fptr, TDOUBLE, "TPRSTP", &(T_PRINT_STEP), NULL, &status);
	printerror(status);
	fits_read_key(fptr, TDOUBLE, "TMAX", &(T_MAX), NULL, &status);
	printerror(status);
	/* fits_read_key(fptr, TDOUBLE, "TRCFAC", &(TRC_FACTOR), NULL, &status); */
	fits_read_key(fptr, TDOUBLE, "S2BMAX", &(SIN2BETA_MAX), NULL, &status);
	printerror(status);
	fits_read_key(fptr, TDOUBLE, "TENEDI", &(TERMINAL_ENERGY_DISPLACEMENT),
		     								NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RMAX", &(R_MAX), NULL, &status);
	/* fits_read_key(fptr, TDOUBLE, "DTMAX", &(DT_MAX), NULL, &status); */
	/* fits_read_key(fptr, TDOUBLE, "DTMIN", &(DT_MIN), NULL, &status); */
	printerror(status);
	//fits_read_key(fptr, TDOUBLE, "MINLR", &(MIN_LAGRANGIAN_RADIUS),
	//									NULL, &status);
	fits_read_key(fptr, TDOUBLE, "MEGAYR", &(MEGA_YEAR), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "SMDYN", &(SOLAR_MASS_DYN), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "METZ", &(METALLICITY), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "WINFAC", &(WIND_FACTOR), NULL, &status);
	//fits_read_key(fptr, TDOUBLE, "MINR", &(MINIMUM_R), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CENMAS", &(cenma.m), NULL, &status);
	//fits_read_key(fptr, TDOUBLE, "CENMAE", &(cenma.E), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DTFACT", &(DT_FACTOR), NULL, &status);
	fits_read_key(fptr, TINT, "BINSIN", &(BINSINGLE), NULL, &status);
	fits_read_key(fptr, TINT, "BINBIN", &(BINBIN), NULL, &status);
	fits_read_key(fptr, TINT, "BSFEWB", &(BINSINGLE_FEWBODY), NULL, &status);
	fits_read_key(fptr, TINT, "BBFEWB", &(BINBIN_FEWBODY), NULL, &status);
	fits_read_key(fptr, TINT, "OPERST", &(ORIGINAL_PERTURB_STARS), 
									NULL, &status);
	printerror(status);
	/* variables related to tidal truncation */
	fits_read_key(fptr, TDOUBLE, "MAXR", &(max_r), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RTIDAL", &(Rtidal), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "TMLOSS", &(TidalMassLoss), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "PREVDT", &(Prev_Dt), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "ORBR", &(orbit_r), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "OTMLOS", &(OldTidalMassLoss), 
									NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DTMLOS", &(DTidalMassLoss), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "PREVDT", &(Prev_Dt), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "ETIDAL", &(Etidal), NULL, &status);
	printerror(status);
	/* variables related to binaries */
	fits_read_key(fptr, TLONG, "NB", &(N_b), NULL, &status);
	fits_read_key(fptr, TLONG, "NBB", &(N_bb), NULL, &status);
	fits_read_key(fptr, TLONG, "NBS", &(N_bs), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "MB", &(M_b), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "EB", &(E_b), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DEBB", &(DE_bb), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DEBS", &(DE_bs), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DBEBB", &(Delta_BE_bb), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DBEBS", &(Delta_BE_bs), NULL, &status);
	printerror(status);
	/* variables related to core  */
	fits_read_key(fptr, TDOUBLE, "NCORE", &(N_core), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "TRC", &(Trc), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "COREDE", &(rho_core), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "COREV", &(v_core), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CORER", &(core_radius), NULL, &status);
	printerror(status);
	/* variables related to escaped stars  */
	fits_read_key(fptr, TDOUBLE, "EESC", &(Eescaped), NULL, &status);
	/* FIXME_ATO: Need to put Ebescaped (the total binding energy of escaped binaries) here */
	fits_read_key(fptr, TDOUBLE, "JESC", &(Jescaped), NULL, &status);
	printerror(status);
	/* Total mass and co.  */
	fits_read_key(fptr, TDOUBLE, "INMTOT", &(initial_total_mass),
		     							NULL, &status);
	fits_read_key(fptr, TDOUBLE, "MTOTAL", &(Mtotal), NULL, &status);
	printerror(status);
	/* everything else except arrays */
	/* reading outprefix does not make much sense :-) 
	 * fits_read_key(fptr, TSTRING, "OUTPRE", &(outprefix), NULL, &status);*/
	fits_read_key(fptr, TINT, "SEFCNT", &(se_file_counter), NULL, &status);
	/* fits_read_key(fptr, TLONG, "ERRSTA", &(errstat), NULL, &status); */
	fits_read_key(fptr, TLONG, "TCOUNT", &(tcount), NULL, &status);
	fits_read_key(fptr, TLONG, "ECHECK", &(Echeck), NULL, &status);
	fits_read_key(fptr, TLONG, "SNAPNO", &(snap_num), NULL, &status);
	fits_read_key(fptr, TLONG, "STCNT", &(StepCount), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DECORS", &(rho_core_single), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DECORB", &(rho_core_bin), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RHSIN", &(rh_single), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RHBIN", &(rh_binary), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DT", &(Dt), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "S2BETA", &(Sin2Beta), NULL, &status);
	printerror(status);

	nelem = MASS_PC_COUNT;
	if ((MAX_INDEX + 5) > nelem){
		nelem = MAX_INDEX + 5;
	}
	if (rng_size>nelem){
		nelem = rng_size;
	}
	if (no_of_doub>nelem){
		nelem = no_of_doub;
	}
	int_arr = malloc(nelem*sizeof(int));
	dbl_arr = malloc(nelem*sizeof(double));
	lng_arr = malloc(nelem*sizeof(long));

	/* reading and restoring mass_pc[] and related arrays*/
	mass_r = malloc((NUM_MASS_RADII_BINS + 1) 
				* MASS_PC_COUNT * sizeof(double));
	ave_mass_r = malloc((NUM_MASS_RADII_BINS + 1) 
				* MASS_PC_COUNT * sizeof(double));
	no_star_r = malloc((NUM_MASS_RADII_BINS + 1) 
				* MASS_PC_COUNT * sizeof(double));
	densities_r = malloc((NUM_MASS_RADII_BINS + 1) 
				* MASS_PC_COUNT * sizeof(double));
	mass_pc = malloc(MASS_PC_COUNT * sizeof(double));
	frow = 1; felem = 1;  
	fits_read_col(fptr, TDOUBLE, 1, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	printerror(status);
	for(i=0; i<MASS_PC_COUNT; i++){
		mass_pc[i] = dbl_arr[i];
	}
	/* comp_mass_percent() will be called after reading position etc. */
	/* reading and restoring RNG state */
	frow = 1; felem = 1;  
	fits_read_col(fptr, TINT, 3, frow, felem, nelem, NULL, int_arr, 
			&anynull, &status);
	printerror(status);
	rng_st_ptr = rng->state;
	for(i=0; i<rng_size; i++){
		rng_st_ptr[i] = int_arr[i];
	}

	/* reading and restoring double precision variables */
	frow = 1; felem = 1;  
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, NULL, dbl_arr, 
			&anynull, &status);
	printerror(status);
	for(i=0; i<no_of_doub; i++){
		dvar[i] = dbl_arr[i];
	}
	/* variables related to progress */
	TotalTime = dvar[0];
	/* Total Energy in various forms */
	Etotal.tot = dvar[1];
	Etotal.New = dvar[2];
	Etotal.ini = dvar[3];
	Etotal.K = dvar[4];
	Etotal.P = dvar[5];
	/* Sub timestep stuff */
	sub.totaltime = dvar[6];
	sub.rmax = dvar[7];
	/* Input file parameters */
	T_PRINT_STEP = dvar[8];
	T_MAX = dvar[9];
	/* TRC_FACTOR = dvar[10]; */
	SIN2BETA_MAX = dvar[12];
	TERMINAL_ENERGY_DISPLACEMENT = dvar[13];
	R_MAX = dvar[15];
	/* DT_MAX = dvar[16]; */
	/* DT_MIN = dvar[17]; */
	//MIN_LAGRANGIAN_RADIUS = dvar[18];
	MEGA_YEAR = dvar[19];
	SOLAR_MASS_DYN = dvar[20];
	METALLICITY = dvar[21];
	WIND_FACTOR = dvar[22];
	//MINIMUM_R = dvar[23];
	cenma.m = dvar[24];
	DT_FACTOR = dvar[27];
	/* variables related to tidal truncation */
	max_r = dvar[28];
	Rtidal = dvar[29];
	TidalMassLoss = dvar[30];
	Prev_Dt = dvar[31];
	orbit_r = dvar[32];
	OldTidalMassLoss = dvar[33];
	DTidalMassLoss = dvar[34];
	Prev_Dt = dvar[35];
	Etidal = dvar[36];
	/* variables related to binaries */
	M_b = dvar[37];
	E_b = dvar[38];
	DE_bb = dvar[39];
	DE_bs = dvar[40];
	Delta_BE_bb = dvar[41];
	Delta_BE_bs = dvar[42];
	/* variables related to core  */
	N_core = dvar[43];
	Trc = dvar[44];
	rho_core = dvar[45];
	v_core = dvar[46];
	core_radius = dvar[47];
	/* variables related to escaped stars  */
	Eescaped = dvar[48];
	/* FIXME_ATO: Need to put Ebescaped (the total binding energy of escaped binaries) here */
	Jescaped = dvar[49];
	/* Total mass and co.  */
	initial_total_mass = dvar[50];
	Mtotal = dvar[51];
	/* everything else except arrays */
	rho_core_single = dvar[52];
	rho_core_bin = dvar[53];
	rh_single = dvar[54];
	rh_binary = dvar[55];
	Dt = dvar[56];
	Sin2Beta = dvar[57];
	cenma.E = dvar[58];

	/* memory allocation for various arrays */
	star = (star_t *) malloc(N_STAR_DIM * sizeof(star_t));
	binary = (binary_t *) malloc(N_STAR_DIM * sizeof(binary_t));

	fits_close_file(fptr, &status);
	printerror(status);
	free(int_arr); free(dbl_arr); free(lng_arr); free(dvar);
}

void read_fits_file_parameters(char *fits_file_name, gsl_rng *rng) {
	fitsfile *fptr;
	int status;

	status = 0;
	/* open file */
	fits_open_file(&fptr, fits_file_name, READONLY, &status);

	printerror(status);

	if (move_to_hdu(fptr, "CLUSTER_STARS") == 0) {
		read_fits_file_parameters_old(fptr, rng);
	} else if (move_to_hdu(fptr, "META_INFO") == 0) {
		read_fits_file_parameters_new(fptr, rng);
	} else {
		eprintf("Cannot read the input file\n");
		exit_cleanly(143);
	}
	
}
