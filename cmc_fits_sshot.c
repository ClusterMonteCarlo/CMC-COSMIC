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

#define ALLOC_TSTUFF \
	ttype = saf_malloc(tfields*sizeof(char *)); \
	tform = saf_malloc(tfields*sizeof(char *)); \
	tunit = saf_malloc(tfields*sizeof(char *)); \
	for(i=0; i<tfields; i++){ \
		ttype[i] = saf_malloc(1024*sizeof(char)); \
		tform[i] = saf_malloc(1024*sizeof(char)); \
		tunit[i] = saf_malloc(1024*sizeof(char)); \
	} \

#define FREE_TSTUFF \
	for(i=0;i<tfields;i++){ \
		free(ttype[i]); \
		free(tform[i]); \
		free(tunit[i]); \
	} \
	free(ttype); free(tform); free(tunit); 

/* Using DBL_MAX in FITS files is not possible with glibc in GNU/Linux 
 * so we use an artificial large number */
#define LARGE_DISTANCE 1.0e40

/* exit with the error message and the status provided */
static void ext_err(char *mess, int status){
	eprintf("%s\n", mess);
	exit_cleanly(status);
}

/* safe malloc */
static void *saf_malloc(size_t size){
	void *dummy;

	if((dummy=malloc(size)) == NULL){
		ext_err("Memory allocation failed.", 102);
	} 
	return dummy;
}

void write_restart_param(fitsfile *fptr, gsl_rng *rng){
	char extname[1024];	/* extension name */
	int tfields;       	/* number of columns */
	long nrows;	
	long firstrow, firstelem;
	char **ttype, **tform, **tunit;
	int status;
	long i;
	struct rng_t113_state rng_st;
	double *dbl_arr;
	long *lng_arr;
	int *int_arr;
	char *rng_st_ptr;
	size_t rng_size=gsl_rng_size(rng); 
	double *dvar;
	int no_of_doub;


	status = 0;
	sprintf(extname,"SF_RESTART_PARAM");
	tfields = 4; 
	no_of_doub = 59;
	dvar=malloc(no_of_doub*sizeof(double));

	nrows = MASS_PC_COUNT;
	if (rng_size>nrows){
		nrows = rng_size;
	}
	if (no_of_doub>nrows){
		nrows = no_of_doub;
	}
	dbl_arr = saf_malloc(nrows*sizeof(double));
	lng_arr = saf_malloc(nrows*sizeof(long));
	int_arr = saf_malloc(nrows*sizeof(int));
	firstrow  = 1;  /* first row in table to write   */
	firstelem = 1;  /* first element in row          */

	ALLOC_TSTUFF
	sprintf(ttype[0],"Mass_PC");
	sprintf(tform[0],"1D");
	sprintf(tunit[0],"Internal");
	sprintf(ttype[1],"IndTab");
	sprintf(tform[1],"1J");
	sprintf(tunit[1],"Internal");
	sprintf(ttype[2],"RNG State");
	sprintf(tform[2],"1I");
	sprintf(tunit[2],"Internal");
	sprintf(ttype[3],"Variables");
	sprintf(tform[3],"1D");
	sprintf(tunit[3],"Various");

	/* FIXME */
	/* need to add MMIN */
	/* FIXME */
	
	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
              tunit, extname, &status);
	/* variables related to progress */
	fits_write_key(fptr, TDOUBLE, "Time", &TotalTime, 
			"Age of cluster", &status);
	fits_write_key(fptr, TLONG, "Step", &tcount, 
			"Iteration Step", &status);
	/* random number state */
	get_rng_t113(&rng_st);
	fits_write_key(fptr, TULONG, "RNG_Z1", &(rng_st.z1), 
			"RNG STATE Z1", &status);
	fits_write_key(fptr, TULONG, "RNG_Z2", &(rng_st.z2), 
			"RNG STATE Z2", &status);
	fits_write_key(fptr, TULONG, "RNG_Z3", &(rng_st.z3), 
			"RNG STATE Z3", &status);
	fits_write_key(fptr, TULONG, "RNG_Z4", &(rng_st.z4), 
			"RNG STATE Z4", &status);
	/* Total Energy in various forms */
	fits_write_key(fptr, TDOUBLE, "Etotal", &(Etotal.tot), 
			"Total Energy of Cluster", &status);
	fits_write_key(fptr, TDOUBLE, "EtotalN", &(Etotal.New), 
			"Total Energy of Cluster", &status);
	fits_write_key(fptr, TDOUBLE, "EtotalI", &(Etotal.ini), 
			"Initial Total Energy", &status);
	fits_write_key(fptr, TDOUBLE, "KEtotal", &(Etotal.K), 
			"Total Kinetic Energy", &status);
	fits_write_key(fptr, TDOUBLE, "PEtotal", &(Etotal.P), 
			"Total Potential Energy", &status);
	/* Sub timestep stuff */
	fits_write_key(fptr, TLONG, "S_NMAX", &(sub.N_MAX), 
			"Substep NMAX", &status);
	fits_write_key(fptr, TLONG, "S_CNT", &(sub.count), 
			"Substep count", &status);
	fits_write_key(fptr, TLONG, "S_FACT", &(sub.FACTOR), 
			"Substep FACTOR", &status);
	fits_write_key(fptr, TDOUBLE, "S_Ttime", &(sub.totaltime), 
			"Substep Total time", &status);
	fits_write_key(fptr, TDOUBLE, "S_rmax", &(sub.rmax), 
			"Substep rmax", &status);
	/* various N's of cluster */
	fits_write_comment(fptr,
		"It is not clear what the following different N's do", &status);
	fits_write_key(fptr, TLONG, "NMAX", &(clus.N_MAX), 
			"Number of stars in cluster (bound)", &status);
	fits_write_key(fptr, TLONG, "NMAXN", &(clus.N_MAX_NEW), 
			"Number of stars in cluster ", &status);
	fits_write_key(fptr, TLONG, "NSTAR", &(clus.N_STAR), 
			"Number of stars in cluster (initial)", &status);
	fits_write_key(fptr, TLONG, "NSTARN", &(clus.N_STAR_NEW), 
			"Number of stars in cluster (all)", &status);
	fits_write_key(fptr, TLONG, "NBINARY", &(clus.N_BINARY), 
			"Number of binaries in cluster", &status);
	/* Input file parameters */
	/* removed DUMP_ENERGY_DISPLACEMENT, BINSINGLE_PROB_SCALE_FAC, and BINBIN_PROB_SCALE_FAC since they are unused */
	fits_write_comment(fptr,
		"Input file parameters", &status);
	fits_write_key(fptr, TLONG, "NSTRDIM", &(N_STAR_DIM), 
			"N_STAR_DIM", &status);
	fits_write_key(fptr, TLONG, "MPCCNT", &(MASS_PC_COUNT), 
			"MASS_PC_COUNT", &status);
	/* fits_write_key(fptr, TLONG, "NTRY", &(N_TRY), 
	   "N_TRY", &status); */
	fits_write_key(fptr, TLONG, "STEVOL", &(STELLAR_EVOLUTION), 
			"STELLAR_EVOLUTION", &status);
	fits_write_key(fptr, TLONG, "SSCOLL", &(SS_COLLISION), 
			"SS_COLLISION", &status);
	fits_write_key(fptr, TLONG, "DUMPS", &(DUMPS), 
			"DUMPS", &status);
	fits_write_key(fptr, TLONG, "ECONS", &(E_CONS), 
			"E_CONS", &status);
	/* fits_write_key(fptr, TLONG, "DTMOD", &(DT_MODE), 
	   "DT_MODE", &status); */
	fits_write_key(fptr, TLONG, "PERTURB", &(PERTURB), 
			"PERTURB", &status);
	fits_write_key(fptr, TLONG, "NUMMASS", &(NUM_MASS), 
			"NUM_MASS", &status);
	fits_write_key(fptr, TLONG, "TOTPAR", &(TOTAL_PARAMS), 
			"TOTAL_PARAMS", &status);
	fits_write_key(fptr, TLONG, "NMSRADB", &(NUM_MASS_RADII_BINS), 
			"NUM_MASS_RADII_BINS", &status);
	fits_write_key(fptr, TLONG, "NCNSTR", &(NUM_CENTRAL_STARS), 
			"NUM_CENTRAL_STARS", &status);
	fits_write_key(fptr, TDOUBLE, "TPRSTP", &(T_PRINT_STEP), 
			"T_PRINT_STEP", &status);
	fits_write_key(fptr, TDOUBLE, "TMAX", &(T_MAX), 
			"T_MAX", &status);
	/* fits_write_key(fptr, TDOUBLE, "TRCFAC", &(TRC_FACTOR), 
	   "TRC_FACTOR", &status); */
	fits_write_key(fptr, TDOUBLE, "S2BMAX", &(SIN2BETA_MAX), 
			"SIN2BETA_MAX", &status);
	fits_write_key(fptr, TDOUBLE, "TENEDI", &(TERMINAL_ENERGY_DISPLACEMENT), 
			"TERMINAL_ENERGY_DISPLACEMENT", &status);
	fits_write_key(fptr, TDOUBLE, "RMAX", &(R_MAX), 
			"R_MAX", &status);
	/* fits_write_key(fptr, TDOUBLE, "DTMAX", &(DT_MAX), 
	   "DT_MAX", &status);  */
	/* fits_write_key(fptr, TDOUBLE, "DTMIN", &(DT_MIN), 
	   "DT_MIN", &status); */
	fits_write_key(fptr, TDOUBLE, "MINLR", &(MIN_LAGRANGIAN_RADIUS), 
			"MIN_LAGRANGIAN_RADIUS", &status);
	fits_write_key(fptr, TDOUBLE, "MEGAYR", &(MEGA_YEAR), 
			"MEGA_YEAR", &status);
	fits_write_key(fptr, TDOUBLE, "SMDYN", &(SOLAR_MASS_DYN), 
			"SOLAR_MASS_DYN", &status);
	fits_write_key(fptr, TDOUBLE, "METZ", &(METALLICITY), 
			"METALLICITY", &status);
	fits_write_key(fptr, TDOUBLE, "WINFAC", &(WIND_FACTOR), 
			"WIND_FACTOR", &status);
	fits_write_key(fptr, TDOUBLE, "MINR", &(MINIMUM_R), 
			"MINIMUM_R", &status);
	fits_write_key(fptr, TDOUBLE, "CENMAS", &(cenma.m), 
			"CentralMass_mass", &status);
	fits_write_key(fptr, TDOUBLE, "CENMAE", &(cenma.E), 
			"CentralMass_energy", &status);
	fits_write_key(fptr, TDOUBLE, "DTFACT", &(DT_FACTOR), 
			"DT_FACTOR", &status);
	fits_write_key(fptr, TINT, "BINSIN", &(BINSINGLE), 
			"BINSINGLE", &status);
	fits_write_key(fptr, TINT, "BINBIN", &(BINBIN), 
			"BINBIN", &status);
	fits_write_key(fptr, TINT, "BSFEWB", &(BINSINGLE_FEWBODY), 
			"BINSINGLE_FEWBODY", &status);
	fits_write_key(fptr, TINT, "BBFEWB", &(BINBIN_FEWBODY), 
			"BINBIN_FEWBODY", &status);
	fits_write_key(fptr, TINT, "OPERST", &(ORIGINAL_PERTURB_STARS), 
			"ORIGINAL_PERTURB_STARS", &status);
	/* variables related to tidal truncation */
	fits_write_comment(fptr,
		"Tidal Truncation parameters", &status);
	fits_write_key(fptr, TDOUBLE, "MAXR", &max_r, 
			"max_r", &status);
	fits_write_key(fptr, TDOUBLE, "RTIDAL", &Rtidal, 
			"Rtidal", &status);
	fits_write_key(fptr, TDOUBLE, "TMLOSS", &TidalMassLoss, 
			"TidalMassLoss", &status);
	fits_write_key(fptr, TDOUBLE, "PREVDT", &Prev_Dt, 
			"Prev_Dt", &status);
	fits_write_key(fptr, TDOUBLE, "ORBR", &orbit_r, 
			"orbit_r", &status);
	fits_write_key(fptr, TDOUBLE, "OTMLOS", &OldTidalMassLoss, 
			"OldTidalMassLoss", &status);
	fits_write_key(fptr, TDOUBLE, "DTMLOS", &DTidalMassLoss, 
			"DTidalMassLoss", &status);
	fits_write_key(fptr, TDOUBLE, "PREVDT", &Prev_Dt, 
			"Prev_Dt", &status);
	fits_write_key(fptr, TDOUBLE, "ETIDAL", &Etidal, 
			"Etidal", &status);
	/* variables related to binaries */
	fits_write_comment(fptr,
		"Binary parameters", &status);
	fits_write_key(fptr, TLONG, "NB", &N_b, 
			"N_b", &status);
	fits_write_key(fptr, TLONG, "NBB", &N_bb, 
			"N_bb", &status);
	fits_write_key(fptr, TLONG, "NBS", &N_bs, 
			"N_bs", &status);
	fits_write_key(fptr, TDOUBLE, "MB", &M_b, 
			"M_b", &status);
	fits_write_key(fptr, TDOUBLE, "EB", &E_b, 
			"E_b", &status);
	fits_write_key(fptr, TDOUBLE, "DEBB", &DE_bb, 
			"DE_bb", &status);
	fits_write_key(fptr, TDOUBLE, "DEBS", &DE_bs, 
			"DE_bs", &status);
	fits_write_key(fptr, TDOUBLE, "DBEBB", &Delta_BE_bb, 
			"Delta_BE_bb", &status);
	fits_write_key(fptr, TDOUBLE, "DBEBS", &Delta_BE_bs, 
			"Delta_BE_bs", &status);
	/* variables related to core  */
	fits_write_comment(fptr,
		"Core parameters", &status);
	fits_write_key(fptr, TDOUBLE, "NCORE", &N_core, 
			"N_core", &status);
	fits_write_key(fptr, TDOUBLE, "TRC", &Trc, 
			"Trc", &status);
	fits_write_key(fptr, TDOUBLE, "COREDE", &rho_core, 
			"rho_core", &status);
	fits_write_key(fptr, TDOUBLE, "COREV", &v_core, 
			"v_core", &status);
	fits_write_key(fptr, TDOUBLE, "CORER", &core_radius, 
			"core_radius", &status);
	/* variables related to escaped stars  */
	fits_write_comment(fptr,
		"Escaped stars' parameters", &status);
	fits_write_key(fptr, TDOUBLE, "EESC", &Eescaped, 
			"Eescaped", &status);
	/* FIXME_ATO: Need to put Ebescaped (the total binding energy of escaped binaries) here */
	fits_write_key(fptr, TDOUBLE, "JESC", &Jescaped, 
			"Jescaped", &status);
	/* Total mass and co.  */
	fits_write_comment(fptr,
		"Total mass and co.", &status);
	fits_write_key(fptr, TDOUBLE, "INMTOT", &initial_total_mass, 
			"initial_total_mass", &status);
	fits_write_key(fptr, TDOUBLE, "MTOTAL", &Mtotal, 
			"Mtotal", &status);
	/* everything else except arrays */
	fits_write_comment(fptr,
		"Everything else", &status);
	fits_write_key(fptr, TINT, "SEFCNT", &se_file_counter, 
			"se_file_counter", &status);
	/* fits_write_key(fptr, TLONG, "ERRSTA", &errstat, 
	   "errstat", &status); */
	fits_write_key(fptr, TLONG, "TCOUNT", &tcount, 
			"tcount", &status);
	fits_write_key(fptr, TLONG, "ECHECK", &Echeck, 
			"Echeck", &status);
	fits_write_key(fptr, TLONG, "SNAPNO", &snap_num, 
			"snap_num", &status);
	fits_write_key(fptr, TLONG, "STCNT", &StepCount, 
			"StepCount", &status);
	fits_write_key(fptr, TDOUBLE, "DECORS", &rho_core_single, 
			"rho_core_single", &status);
	fits_write_key(fptr, TDOUBLE, "DECORB", &rho_core_bin, 
			"rho_core_bin", &status);
	fits_write_key(fptr, TDOUBLE, "RHSIN", &rh_single, 
			"rh_single", &status);
	fits_write_key(fptr, TDOUBLE, "RHBIN", &rh_binary, 
			"rh_binary", &status);
	fits_write_key(fptr, TDOUBLE, "DT", &Dt, 
			"Dt", &status);
	fits_write_key(fptr, TDOUBLE, "S2BETA", &Sin2Beta, 
			"Sin2Beta", &status);

	/* variables related to progress */
	dvar[0]  = TotalTime;
	/* Total Energy in various forms */
	dvar[1]  = Etotal.tot;
	dvar[2]  = Etotal.New;
	dvar[3]  = Etotal.ini;
	dvar[4]  = Etotal.K;
	dvar[5]  = Etotal.P;
	/* Sub timestep stuff */
	dvar[6]  = sub.totaltime;
	dvar[7]  = sub.rmax;
	/* Input file parameters */
	dvar[8]  = T_PRINT_STEP;
	dvar[9]  = T_MAX;
	/* dvar[10] = TRC_FACTOR; */
	dvar[11] = 0.0; /* this was mistakenly skipped in orig. imple. (ato)*/
	dvar[12] = SIN2BETA_MAX;
	dvar[13] = TERMINAL_ENERGY_DISPLACEMENT;
	dvar[14] = 0.0;
	dvar[15] = R_MAX;
	/* dvar[16] = DT_MAX; */
	/* dvar[17] = DT_MIN; */
	dvar[18] = MIN_LAGRANGIAN_RADIUS;
	dvar[19] = MEGA_YEAR;
	dvar[20] = SOLAR_MASS_DYN;
	dvar[21] = METALLICITY;
	dvar[22] = WIND_FACTOR;
	dvar[23] = MINIMUM_R;
	dvar[24] = cenma.m;
	dvar[25] = 0.0;
	dvar[26] = 0.0;
	dvar[27] = DT_FACTOR;
	/* variables related to tidal truncation */
	dvar[28] = max_r;
	dvar[29] = Rtidal;
	dvar[30] = TidalMassLoss;
	dvar[31] = Prev_Dt;
	dvar[32] = orbit_r;
	dvar[33] = OldTidalMassLoss;
	dvar[34] = DTidalMassLoss;
	dvar[35] = Prev_Dt;
	dvar[36] = Etidal;
	/* variables related to binaries */
	dvar[37] = M_b;
	dvar[38] = E_b;
	dvar[39] = DE_bb;
	dvar[40] = DE_bs;
	dvar[41] = Delta_BE_bb;
	dvar[42] = Delta_BE_bs;
	/* variables related to core  */
	dvar[43] = N_core;
	dvar[44] = Trc;
	dvar[45] = rho_core;
	dvar[46] = v_core;
	dvar[47] = core_radius;
	/* variables related to escaped stars  */
	dvar[48] = Eescaped;
	/* FIXME_ATO: Need to put Ebescaped (the total binding energy of escaped binaries) here */
	dvar[49] = Jescaped;
	/* Total mass and co.  */
	dvar[50] = initial_total_mass;
	dvar[51] = Mtotal;
	/* everything else except arrays */
	dvar[52] = rho_core_single;
	dvar[53] = rho_core_bin;
	dvar[54] = rh_single;
	dvar[55] = rh_binary;
	dvar[56] = Dt;
	dvar[57] = Sin2Beta;
	dvar[58] = cenma.E;

	for(i=0; i<nrows; i++){
		dbl_arr[i] = 0.0;
	}
	for(i=0; i<MASS_PC_COUNT; i++){
		dbl_arr[i] = mass_pc[i];
	}
	fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows,
			dbl_arr, &status);
	printerror(status);
	
	rng_st_ptr = rng->state;
	for(i=0; i<nrows; i++){
		int_arr[i] = 0;
	}
	for(i=0; i<rng_size; i++){
		/* XXX I am a not sure the below casting works !!! (ato) */
		int_arr[i] += (int) rng_st_ptr[i];
	}
	fits_write_col(fptr, TINT, 3, firstrow, firstelem, nrows,
			int_arr, &status);
	printerror(status);
	
	for(i=0; i<nrows; i++){
		dbl_arr[i] = 0.0;
	}
	for(i=0; i<no_of_doub; i++){
		dbl_arr[i] = dvar[i];
	}
	fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows,
			dbl_arr, &status);
	printerror(status);
	FREE_TSTUFF
	free(dbl_arr); free(lng_arr); free(int_arr); free(dvar);

}

void write_ss_dyn_param(fitsfile *fptr){
	char extname[1024];	/* extension name */
	int tfields;       	/* number of columns */
	long nrows;	
	long firstrow, firstelem;
	char **ttype, **tform, **tunit;
	int status;

	long i, N;
	double *dbl_arr;
	long *lng_arr;

	status = 0;

	sprintf(extname,"SS_DYN_PARAM");
     	N = clus.N_STAR_NEW;
	tfields = 9; nrows = N+2;
	firstrow  = 1;  /* first row in table to write   */
	firstelem = 1;  /* first element in row          */
	dbl_arr = saf_malloc((N+2)*sizeof(double));
	lng_arr = saf_malloc((N+2)*sizeof(long));

	ALLOC_TSTUFF
	sprintf(ttype[0],"Mass");
	sprintf(tform[0],"1D");
	sprintf(tunit[0],"Nbody");
	sprintf(ttype[1],"Position");
	sprintf(tform[1],"1D");
	sprintf(tunit[1],"Nbody");
	sprintf(ttype[2],"Vr");
	sprintf(tform[2],"1D");
	sprintf(tunit[2],"Nbody");
	sprintf(ttype[3],"Vt");
	sprintf(tform[3],"1D");
	sprintf(tunit[3],"Nbody");
	sprintf(ttype[4],"Pericenter");
	sprintf(tform[4],"1D");
	sprintf(tunit[4],"Nbody");
	sprintf(ttype[5],"Binary index");
	sprintf(tform[5],"1J");
	sprintf(tunit[5],"Index");
	sprintf(ttype[6],"Interaction Flag");
	sprintf(tform[6],"1J");
	sprintf(tunit[6],"Flag");
	sprintf(ttype[7],"Old k");
	sprintf(tform[7],"1J");
	sprintf(tunit[7],"Index");
	sprintf(ttype[8],"Sphi");
	sprintf(tform[8],"1D");
	sprintf(tunit[8],"Nbody");

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status);
	fits_write_key(fptr, TLONG, "NSTAR", &(clus.N_STAR), 
			"No of Stars (initial)", &status);
	fits_write_key(fptr, TLONG, "NSTARN", &(clus.N_STAR_NEW), 
			"No of Stars (all)", &status);
	fits_write_key(fptr, TLONG, "NMAX", &(clus.N_MAX), 
			"No of Stars (bound)", &status);
	fits_write_key(fptr, TLONG, "NMAXN", &(clus.N_MAX_NEW), 
			"No of Stars (bound, new)", &status);

	/* mass */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].m;
	}
	fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, dbl_arr,
                   &status);
	/* position */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].r;
	}
	fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, dbl_arr,
                   &status);
	/* vr */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].vr;
	}
	fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nrows, dbl_arr,
                   &status);
	/* vt */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].vt;
	}
	fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows, dbl_arr,
                   &status);
	/* r_peri */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].r_peri;
	}
	fits_write_col(fptr, TDOUBLE, 5, firstrow, firstelem, nrows, dbl_arr,
                   &status);

	/* binindex */
	for(i=0;i<=N+1;i++){
		lng_arr[i] = star[i].binind;
	}
	fits_write_col(fptr, TLONG, 6, firstrow, firstelem, nrows, lng_arr,
                   &status);

	/* interacted? */
	for(i=0;i<=N+1;i++){
		lng_arr[i] = star[i].interacted;
	}
	fits_write_col(fptr, TLONG, 7, firstrow, firstelem, nrows, lng_arr,
                   &status);

	/* oldk */
	/* 
	for(i=0;i<=N+1;i++){
		lng_arr[i] = star[i].oldk;
	}
	fits_write_col(fptr, TLONG, 8, firstrow, firstelem, nrows, lng_arr, &status);
	*/
	/* potential */
	for(i=0;i<=N+1;i++){
		dbl_arr[i] = star[i].phi;
	}
	fits_write_col(fptr, TDOUBLE, 9, firstrow, firstelem, nrows, dbl_arr,
                   &status);

	FREE_TSTUFF
	free(dbl_arr); free(lng_arr);

	printerror(status);
}

void write_ss_se_param(fitsfile *fptr){
	char extname[1024];	/* extension name */
	int tfields;       	/* number of columns */
	long nrows;	
	long firstrow, firstelem;
	char **ttype, **tform, **tunit;
	int status;

	long i, N;
	double *dbl_arr;
	int *int_arr;

	status = 0;

	sprintf(extname,"SS_SE_PARAM");
     	N = clus.N_STAR;
	tfields = 3; nrows = N;
	firstrow  = 1;  /* first row in table to write   */
	firstelem = 1;  /* first element in row          */
	dbl_arr = saf_malloc((N+2)*sizeof(double));
	int_arr = saf_malloc((N+2)*sizeof(int));

	ALLOC_TSTUFF
	sprintf(ttype[0],"Radius");
	sprintf(tform[0],"1D");
	sprintf(tunit[0],"Solar");
	sprintf(ttype[1],"Mass");
	sprintf(tform[1],"1D");
	sprintf(tunit[1],"Solar");
	sprintf(ttype[2],"Type");
	sprintf(tform[2],"1I");
	sprintf(tunit[2],"Internal");
	
	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status);

	/* radius */
	dbl_arr[0] = 0.0;
	for(i=1;i<=N;i++){
		dbl_arr[i] = star[i].rad;
	}
	dbl_arr[N+1] = 0.0;
	fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, dbl_arr,
                   &status);
#ifdef SE
	/* mass */
	dbl_arr[0] = 0.0;
	for(i=1;i<=N;i++){
		dbl_arr[i] = star[i].mass;
	}
	dbl_arr[N+1] = 0.0;
	fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, dbl_arr,
                   &status);
	/* type */
	int_arr[0] = 0.0;
	for(i=1;i<=N;i++){
		int_arr[i] = star[i].k;
	}
	int_arr[N+1] = 0.0;
	fits_write_col(fptr, TINT, 3, firstrow, firstelem, nrows, int_arr,
                   &status);
#endif

	FREE_TSTUFF
	free(dbl_arr); free(int_arr);

	printerror(status);
}

void write_basic_info(fitsfile *fptr){
	char extname[1024];	/* extension name */
	int tfields;       	/* number of columns */
	long nrows;	
	int is_snapshot;
	int status;

	status = 0;
	sprintf(extname,"META_INFO");
	tfields = 0; nrows = 0;		/* this extension has table */

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, NULL, NULL,
                NULL, extname, &status);

	is_snapshot = 1;
	fits_write_key(fptr, TINT, "IS_SS", &is_snapshot, 
			"Is this a snapshot", &status);
	printerror(status);
}


/* Initially the extensions will be:
 * 0: This is just headers, 
 * 1: Meta information written here
 * 2: Restart parameters
 * 3: Single star dynamical parameters
 * 4: Single star stellar evolution parameters
 * 5: Binaries' dynamical parameters
 * 6: Binaries' stellar evolution parameters
 */

void sshot_fits(gsl_rng *rng){
	fitsfile *fptr;
	char filename[1024];
	int status;

	if (sub.count != sub.FACTOR-1){ 
		return;
	}
	status = 0;
	sprintf(filename, "!dummy%ld.fit", tcount);
	fits_create_file(&fptr, filename, &status);
	printerror(status);

	/* 1st Extension, Basic info */
	write_basic_info(fptr);
	/* 2nd Extension: Restart Parameters */
	write_restart_param(fptr, rng);
	/* 3rd Extension: Single Star Dynamical Parameters */
	write_ss_dyn_param(fptr);
	/* 4th Extension: Single Star Stellar Evolution Parameters */
	if (STELLAR_EVOLUTION > 0){
		write_ss_se_param(fptr);
	}

	fits_close_file(fptr, &status);
	printerror(status);
}

#undef ALLOC_TSTUFF 
#undef FREE_TSTUFF
#undef LARGE_DISTANCE
