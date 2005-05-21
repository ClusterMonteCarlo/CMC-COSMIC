 
/* Monte Carlo simulation of Globular Cluster evolution 

	- Kriten Joshi (M.I.T.)  (1997-1999)

*/

/*** Uncomment this for Debug mode (verbose printing) */ 
#define DEBUGMODE

/*** Uncomment this to include code for parallel jobs using MPI ***
#define MPIMODE
*/

#define NRANSI
#include <stdlib.h>
#include <fcntl.h> 
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include "nrutil_d.h"
#include "nrutil_d.c"
#if defined(MPIMODE)
	#include <mpi.h>   
#endif


/*************************** Parameters ******************************/

#define INFINITY 1.0e10	/* Large number, but still INFINITY - 1 <> INFINITY */
#define ZERO 1.0e-20	/* Radius of the zeroth star so that 1/sr[0] is finite */       
#define R1 1.0      /* Factor by which to multiply all radii in mass component 
                       that is left-most in input-file specifications for generated 
					   Plummer model. */
#define R2 1.0      /* Factor by which to multiply all radii in other mass components. */

#define SOLAR_MASS 2.0e33
#define G 6.7e-8
#define KPC 3.0e21
#define YEAR 3.15e7

/*********************** Global Variables ****************************/

double pi, Etotal, Etotal_New, Etotal_initial, KEtotal, PEtotal, Eescaped, Jescaped; 
double TotalTime, max, min, range, actual_total_mass, Dt, max_r; 
double sin2beta, max_sin2beta, avg_sin2beta, avg_X, avg_r, avg_dv, avg_DE, avg_DJ ;
long errstat, N_MAX, idum, n_circular, n_pericenter, n_apocenter, tcount ;
long n_orbit, n_extreme, current_star_i ;
int printstatus;

double *sr, *svr, *svt, *sm ;
double *sLifeTime, *sInitialMass, *sRemnantMass, RealTime, Rtidal, StartTime ;
double *sE, *sJ, *sgravity, *sEI, *sJI, StellarMassLoss, TidalMassLoss ;
double OldTidalMassLoss, DTidalMassLoss, Prev_Dt, Etidal, phi_rtidal, phi_zero ;
double *sKE, *sPE, *sMi_ri, *sMi, *sl, orbit_t[10000], orbit_p[10000], orbit_Mg[10000];
double *srescape, *mass_pc, *mass_r;
double *srnew, *svrnew, *svtnew;
double DEavg, DEmax, DEdiff, *srOld, *sX, *sDE, *sDJ, *sSin2Beta, *sr_peri, *sr_apo ;
double *mass_spec, *dist_params, *actual_dist, *Mtotal, Sin2Beta ;
long *sindex, *IndexTable, *trace, *trace_stat, orbit_step_count, *stype, *old_k ;
long count_1, count_2, count_3, count_4, count_5, N_bb=0, N_bs=0 ;
double mass_1, mass_2, mass_3, mass_4, mass_5, turnoff_mass=0, N_core ;
double *trace_r_ini, *trace_tr, Trc, rho_core, v_core, core_radius ;
double ecc_anomaly, mean_anomaly, orbit_r, r_core, rho_core, v_core ;

double *Smooth_r, *Smooth_m;
FILE *out[7], *logfile, *tracefile1, *tracefile2, *tracefile3, *escfile, *snapfile ;
FILE *stellarfile, *binaryfile;
double *dens_veloc_radii, *countmass, *num_within, *dens_within, *temp_within, *velaver;
double *temp_within_core, *dens_within_core, *num_within_core;
long *whichrad;

int myid=1000, numprocs;
char outprefix[100];

/****** Binary parameters *******/
double *ba, *be, *bm1, *bm2, *bini_a, *bini_e, *bini_tformed, *bini_tdestroyed ;
double *bini_m1, *bini_m2, *by, M_binary, sec_binary_mass, Delta_BE_bb, Delta_BE_bs ;
double rho_core_single, rho_core_bin, rh_single, rh_binary, E_binding, DE_bb, DE_bs  ;
long *s_binindex, *b_sindex, N_BINARY, *sinteracted, N_STAR_NEW, N_MAX_NEW ;
long sub_N_MAX, sub_count, sub_FACTOR ;
double sub_totaltime, sub_rmax ;

/******************* Input file parameters *************************/

long N_STAR, N_STAR_DIM, T_MAX_COUNT, MASS_PC_COUNT, N_TRY, MODEL, STELLAR_EVOLUTION ;
long DUMPS=0, E_CONS, DT_MODE, MAX_INDEX, INDEX_UNIT, E_NORMALIZE, RESTART ;
long IDUM=0, PERTURB=1; /* IDUM=0 --> set default to -115; PERTURB=1 (ON by default) */
long NUM_MASS, TOTAL_PARAMS, NUM_MASS_RADII_BINS, B_VAR_SWITCH, NUM_CORE_STARS ;
double T_PRINT_STEP, R_PRINT_STEP, T_MAX, TOTAL_MASS, TRC_FACTOR, SIN2BETA_MAX ;
double TERMINAL_ENERGY_DISPLACEMENT, DUMP_ENERGY_DISPLACEMENT, R_MAX ;
double DT_MAX, DT_MIN, MIN_LAGRANGIAN_RADIUS, ALPHA ;
char MASS_PC[1000], MASS_SPEC[1000], MASS_DIST_PARAMS[1000], INPUT_SNAPSHOT_FILE[1000];
long TRACE_N, TRACE_STARS, ORBIT_ELLIPTICAL ;
double TRACE_MASS1, TRACE_MASS2 ;
double ORBIT_PERI, ORBIT_VC, ORBIT_PERIOD, ORBIT_ECCENTRICITY, ORBIT_PHASE, ORBIT_MAX_DT ;

/* Conversion factors from code units to physical units */
double UT, UL, FAMILY, T_RELAX ;

/********************** Function Declarations ************************/

long plummer(double r_max) ;	/* Assign radii to stars */
long plummer_mass_spec(double r_max);
double plummer_f_r(double r) ;	/* Plummer radial profile */
double potential(double r) ;	/* get potential using sgravity[] */
double phi(double r) ;	/* get potential w/o using sgravity[] */
double mass(double r) ;	/* get mass contained within radius <= r */
long gravity(void) ; /* get potential at star locations in sgravity[] */
double get_positions(void) ;	/* get positions and velocities */
void perturb_stars(double Dt) ;	/* take a time step (perturb E,J) */
void check_potential_convergence(void);
double function_Q(long k, double E, double J) ;	/* Q for get_positions() */
double function_minusAbsQ(long k, double E, double J) ; /* -abs(Q) */
double golden(long ax, long bx, long cx, long *xmin, 
			  double E, double J, double (*Q)(long, double, double)) ;
void indexx(long n, double arr[], long indx[]) ;
long FindZero(double (*func)(long, double, double), long x1, long x2,
			  double E, double J) ;
long FindZero_Q(long x1, long x2, double E, double J) ;
double ran2(long *idum) ;
double gauss(double x, double x0, double sigma) ;
long zbrent_Q(long x1, long x2, double E, double J) ;
void ComputeEnergy(void) ;
void NormalizeEnergy(void) ;
void PrintLogOutput(void) ;
double GetTimeStep(void) ;
long CheckStop(void) ;
void ComputeIntermediateEnergy(void) ;
void RecomputeEnergy(void) ;
double kepler_anomaly(double E) ;
double rtbis(double (*func)(double), double x1, double x2, double xacc) ;
double trapzd(double (*func)(double), double a, double b, long n) ;
double qsimp(double (*func)(double), double a, double b) ;
double log_potential(double rd) ;
double log_m(double rd) ;

void splint(double xa[], double ya[], double y2a[], long n, double x, double *y);
void spline(double x[], double y[], long n, double yp1, double ypn, double y2[]);
void read_stellar_data(void);
void DoStellarEvolution(void) ;
void CheckStellarTimeStep(void);
void print_2Dsnapshot(void) ;
void get_physical_units(void) ;
double midpnt(double (*func)(double), double a, double b, int n) ;
double qromo(double (*func)(double), double a, double b,
	double (*choose)(double(*)(double), double, double, int)) ;
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

/******* Binary routines **********/

void assign_binaries(void) ;
double bin_single_sigma(double y) ;

/******* Restart Routines *********/
void restart(void);
void output_restartdump(void);

/*====== I/O & component calculation routines -- Wesley Watters =======*/
double timestep_profile(double x);
void output_rawdump(void);
void output_density_profile(void);
int parser(int argc, char *argv[]);
int starlab_model_fetch_with_masses(void);
int starlab_model_fetch_without_masses(void);
int starlab_model_spec_fetch(void);
int find_mass_bin(double mass);
void GetInitialModel(void) ;
void ComputeLagrangianRadii(void) ;
void PrintFileOutput(void) ;
void ComputeComponentKE(void) ;

/*======= Mass-Spectrum functions -- Wesley Watters =======*/
double find_min(double *list, int NumMembers);
double find_max(double *list, int NumMembers);
double find_nearest(double m_pick, double *list, int NumMembers);
int find_position(double item, double *list, int NumMembers);
double mass_spec_dist(double argument, double *parameters);
long mass_spec_assign(int type);
void CountComponentMass(void) ;

/* Mass function by Cody */
long initial_mass_powerlaw(double alpha);

/************************* Main function call ***************************/

int main(int argc, char *argv[])

{
	static long i, j, k, imin, imax, p, si_minus_p, si_plus_p, sub_imin, si, zk ;

	static clock_t timecount ; 
	static double m_prev, v_prev, rmin, rmax, t_r, n_r, rho_r, avg_v ;
	static double Ai, Dt_local, m_avg, w2_avg, zr_min, zr_max;
/*	static char cont ;
*/
#if defined(MPIMODE)

   	double          startwtime, mytime;
   	int             namelen;
   	char            processor_name[MPI_MAX_PROCESSOR_NAME];
   	MPI_Status      stat ;

    MPI_Init(&argc,&argv);
   	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   	MPI_Get_processor_name(processor_name, &namelen);

   	fprintf(stderr,"Process %d on %s\n", myid, processor_name);

#endif

	if(parser(argc, argv)==0) return 0;


/******  Get tabulated t(r) & Mg(r) for cluster orbit in the galaxy *******/

	if (ORBIT_ELLIPTICAL == 2) {
		orbit_t[0] = orbit_t[1] = 0.0 ; orbit_p[0] = orbit_p[1] = 1.0 ;
		orbit_Mg[0] = qsimp(&log_m, 0, 1.0) ;

		/* Get t(r) for first half of the orbit from r_peri to r_apo  */
		for (i=2 ; i<=10000 ; i++) {
			if ( (1.0+ORBIT_ECCENTRICITY) + 2.0*log(2.0) \
					- 2.0*log(1.0+(1.0+i*0.001)*(1.0+i*0.001)) \
					- (1.0+ORBIT_ECCENTRICITY)/((1.0+i*0.001)*(1.0+i*0.001)) > 0) {
				orbit_t[i] = orbit_t[i-1] + qsimp(&log_potential,1.0+(i-1)*0.001, 1.0+i*0.001) ;
				orbit_p[i] = 1.0 + i*0.001 ;
				orbit_Mg[i] = qsimp(&log_m, 0, orbit_p[i]) ;	
				if (i==2) {
					orbit_t[1] = orbit_t[2] ;
					orbit_p[1] = (orbit_p[0] + orbit_p[2])/2.0 ;
					orbit_Mg[1] = qsimp(&log_m, 0, orbit_p[1]) ;
					orbit_t[2] += orbit_t[1] ;
				}
			}
			else {
				break ;
			}
		}
		orbit_t[i] = orbit_t[i-1] + (orbit_t[i-1]-orbit_t[i-2]) ;
		orbit_p[i] = orbit_p[i-1] + (orbit_p[i-1]-orbit_p[i-2])/2.0 ;
		orbit_Mg[i] = qsimp(&log_m, 0, orbit_p[i]) ;
		orbit_step_count = i ;

		/* Assign same values for second half of orbit, from r_apo to r_peri */
		for (i=orbit_step_count+1 ; i<=2*orbit_step_count ; i++) {
				orbit_t[i] = 2.0*orbit_t[orbit_step_count] - orbit_t[2*orbit_step_count - i] ;
				orbit_p[i] = orbit_p[2*orbit_step_count - i] ;
				orbit_Mg[i] = orbit_Mg[2*orbit_step_count - i] ;
		}
		orbit_step_count = i-1 ;

		/* Compute orbital period in years, by converting to physical units.
		   r_peri = R_c (for log potential) = G = Mg(r_peri) = 1  
		   This implies that, V_c(r_peri) = sqrt(G*M/r_peri) = 1, and 
		   unit of time is r_peri / V_c(r_peri).
		*/
		ORBIT_PERIOD = orbit_t[orbit_step_count] * (ORBIT_PERI*3.0e21) / (ORBIT_VC*1.0e5) / 3.15e7 ;
		
	}


/* Restarting the evolution from saved out_restart file - no need for initial setup*/
	if(RESTART==1) {
		restart();
		goto START_EVOLUTION ;
	}

/* Set up initial conditions (positions and velocities ONLY) 
*/
	GetInitialModel() ;	

	if (STELLAR_EVOLUTION > 0) {
		read_stellar_data() ;
	}

/* Get conversion factors from code units to physical units -- UT, UL, T_RELAX 
*/
	get_physical_units() ;


/* Assigning special TRACE stars */

if (TRACE_STARS > 0) {

	for (i = 1; i <= TRACE_N ; i++) {
			sm[i] = TRACE_MASS2 ;		/* Mass of trace component */
			stype[i] = 2 ;
	}
	for (i = TRACE_N+1 ; i <= N_STAR ; i++) {
			sm[i] = TRACE_MASS1 ;		/* Mass of background stars */
			stype[i] = 1 ;
	}
	count_1 = N_STAR - TRACE_N ;
	count_2 = TRACE_N ;
	mass_1 = count_1 * TRACE_MASS1 ;
	mass_2 = count_2 * TRACE_MASS2 ;

	count_3 = count_4 = count_5 = 0 ;
	mass_3 = mass_4 = mass_5 = 0.0 ;

}

/* Indexing stars by radius and setting the 0th star 
   at radius 0 and (N_STAR+1)th star at INFINITY 
*/
	indexx(N_STAR, sr, sindex) ;
	sindex[0] = 0 ;
	sindex[N_STAR+1] = N_STAR+1 ;

	N_MAX = N_STAR; /* For the time 0, N_MAX is assumed to equal N_STAR */


/* Assign binaries -- requires Indexing to be completed in order
   to compute core velocity dispersion. 
   
   ****** MUST call indexx() first before calling assign_binaries() ******
*/
	if (N_BINARY > 0) {
		assign_binaries() ;		/* assign binaries */
	}

/* Count mass in each component 
*/
	CountComponentMass() ;	


/*	Computing the gravity at the star locations provided in sr[].
	Results returned in sgravity[]. Returns N_MAX. Also computes
	Rtidal using Mtotal[0] and orbit_r.
*/
	orbit_r = R_MAX ;
	errstat = gravity() ;


/****** Get sindex[] positions & get local t_r for TRACE stars ******/
/* Save initial trace star positions and density in out_trace1 file */

if (TRACE_STARS > 0) {

	for (i = 1 ; i <= TRACE_N ; i++) {
		for (j = 1 ; j <= N_MAX ; j++) {
			if (sindex[j] == i) {
				trace[i] = j ;
				break ;
			}
		}
		/* Compute local t_r */
		m_prev = 0.0 ; v_prev = 0.0 ;
		imin = trace[i] - 50 ; imax = trace[i] + 50 ;
		if (imin < 0) imin = 0 ;
		if (imax > N_MAX) imax = N_MAX ;
		for (j = imin+1 ; j < imax ; j++) {
			m_prev += sm[sindex[j]] / N_STAR ;
			v_prev += (svr[sindex[j]]*svr[sindex[j]] + svt[sindex[j]]*svt[sindex[j]]) ;
		}
		rmin = sr[sindex[imin]] ;
		rmax = sr[sindex[imax]] ;
		
		rho_r = m_prev / \
			(4.0*pi/3.0) / (pow(rmax,3.0) - pow(rmin,3.0)) ;

		n_r = (imax - imin) / \
			(4.0*pi/3.0) / (pow(rmax,3.0) - pow(rmin,3.0)) ;

		avg_v = sqrt(v_prev / (imax - imin)) ;

		t_r = 0.065 * avg_v*avg_v*avg_v / rho_r ;	

/*		t_r = t_r * 4.5e10 ;				 convert to years */
/*		rho_r = rho_r *	1.06e6 / pow(8.15,3.0) ;	 convert to M_sun / pc^3 */
/*		n_r = n_r / pow(8.15,3.0) * 10 ;			 convert to pc^-3 and N = 10^6 */

		fprintf(tracefile1, "%5d   %5.3g   %5.3g   %5.3g  %5.3g \n", \
			i, sr[i], rho_r, n_r, t_r) ;

		trace_stat[i] = 0 ;
		trace_r_ini[i] = sr[i] ;
		trace_tr[i] = t_r ;
	}

	fclose(tracefile1) ;
}
/**********************************************************************/

/* Compute initial lagrangian radii 
	--- innermost lagrangian radius is used in CheckStop() 

	ComputeLagrangianRadii() ;
*/

/*	Calculate E and J for stars (per unit mass), and compute Etotal, KEtotal & PEtotal
	E = Potential + 1/2 * (vr^2 + vt^2)  and  J = r * vt
*/
	ComputeEnergy() ;
	Etotal_New = 0.0 ; Eescaped = 0.0 ; Jescaped = 0.0 ;


/*	Normalize the total energy to -0.25 EXACTLY by scaling KE as required
	Etotal, KEtotal and PEtotal are recomputed.
*/	
	if (E_NORMALIZE > 0) NormalizeEnergy() ;

	Etotal_initial = Etotal;  /* Noting the total initial energy, in order 
								to set termination energy. */

/* Compute KE's for each component, etc. 

   ComputeComponentKE() ;

*/


/*	Iteration of computation of positions and velocities to test for 
    convergence to correct potential. (In fact, it does not converge!) 

	check_potential_convergence();  
*/	
	

/*******  Save trace star distribution *******/

	if (TRACE_STARS > 0) {
		for (i = 1 ; i <= TRACE_N ; i++) {
			fprintf(tracefile2, \
				"%5d  %.8g  %5.3g %5.3g  %5.3g  %5.3g  %5.3g  %5.3g  %5.3g  %5.3g\n", \
				i, TotalTime, sr[i], sr_peri[i], sr_apo[i], sl[i], sE[i], sJ[i], \
				sDE[i], sDJ[i]) ;
		}
	}


START_EVOLUTION: 

/*	Printing Results for initial model
*/
	PrintFileOutput() ;


/******************* INITIAL MODEL CREATED. Starting evolution ****************/


#if defined(DEBUGMODE)
	if (printstatus==0 || printstatus==6) 
		fprintf(stderr,"\nStarting time steps... \n\n") ;
#endif

	timecount = clock() ;

	if (RESTART == 0) {
		TotalTime = 0.0 ;
		tcount = 1 ;
	}
	
	sub_count = 0 ;
	for (tcount = tcount; tcount <= 100000 ; tcount++) {

	/* Check for END of simulation 
	*/
		if (CheckStop() != 0) break ;

	/* Get new time step --- also enforces DT_MIN & ORBIT_MAX_DT for
	   elliptical orbits (i.e., if ORBIT_ELLIPTICAL > 0)  
	*/
		Dt = GetTimeStep() ;

	/* Make sure stellar mass loss is not greater than 1% in the timestep
	*/
		if (STELLAR_EVOLUTION > 0) {
		  CheckStellarTimeStep();
		}

		/* if tidal mass loss in previous time step is > 5% reduce 
		   PREVIOUS timestep by 20%
		*/
		if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
			fprintf(stderr,"Prev TidalMassLoss = %.6G  Reducing previous Dt by 20pc\n", \
				TidalMassLoss - OldTidalMassLoss) ;
			fprintf(logfile,"Prev TidalMassLoss = %.6G  Reducing previous Dt by 20pc\n", \
				TidalMassLoss - OldTidalMassLoss) ;
			Dt = Prev_Dt * 0.8 ;
		}
		else if (Dt > 1.1*Prev_Dt && Prev_Dt > 0  &&  (TidalMassLoss - OldTidalMassLoss) > 0.02) {
			fprintf(stderr,"Increasing Dt by 10pc -- Dt = %.6G\n", Dt) ;
			fprintf(logfile,"Increasing Dt by 10pc -- Dt = %.6G\n", Dt) ;
			Dt = Prev_Dt * 1.1 ;
		}

		TotalTime = TotalTime + Dt ; 


	/***************************** setup sub time step *******************************/

		/***** Timestep obtained from GetTimeStep is based on the relaxation time
			   in the core. Now look for a suitable boundary between core and halo.
		******/

		if (sub_count == 0) {		/** do a full timestep **/

			p=20 ;	
			sub_N_MAX = N_MAX ;	/* default -- if Dt does not change much up to r_h */
			sub_FACTOR = 1 ;
			sub_imin = 2000  ;

			for (si = sub_imin ; si <= N_MAX ; si += N_MAX/100) {
				si_minus_p = si-p ; 
				si_plus_p = si+p+1 ;

				if (si_minus_p < 0) {
					si_minus_p = 0 ;
					si_plus_p = 2*p + 1 ;
				}
				else if (si_plus_p > N_MAX) {
					si_plus_p = N_MAX ;
					si_minus_p = N_MAX - 2*p - 1 ;
				}

				/* calculate Ai for zone
				*/	
				m_avg = 0.0 ; w2_avg = 0.0 ; zk = 2*p+2 ;
				for (i = si_minus_p ; i <= si_plus_p ; i++) {
					k = sindex[i] ;
					m_avg += sm[k] ;
					w2_avg += sm[k] * (svr[k]*svr[k]+svt[k]*svt[k]) ;
				}
				m_avg = m_avg / zk ;
				w2_avg = w2_avg * 2.0 / m_avg / zk ;
				zr_min = sr[sindex[si_minus_p]] ; zr_max = sr[sindex[si_plus_p]];
				Ai = 6.0 * zk * m_avg*m_avg / \
					(zr_max*zr_max*zr_max - zr_min*zr_min*zr_min) / pow(w2_avg, 1.5) ;


				Dt_local = SIN2BETA_MAX / Ai * N_STAR ;
				
				if (Dt_local >= Dt * 500 && sub_FACTOR < 500) {	
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/500.0) {				
						sub_FACTOR = 500 ; sub_N_MAX = si ; break ;
					}
				}
				if (Dt_local >= Dt * 100 && sub_FACTOR < 100) {
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/100.0) {
						sub_FACTOR = 100 ; sub_N_MAX = si ;
					}
				}
				else if (Dt_local >= Dt * 50 && sub_FACTOR < 50) {
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/50.0) {
						sub_FACTOR = 50 ; sub_N_MAX = si ;	
					}
				}
				else if (Dt_local >= Dt * 25 && sub_FACTOR < 25) {
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/25.0) {
						sub_FACTOR = 25 ; sub_N_MAX = si ;	
					}
				}
				else if (Dt_local >= Dt * 10 && sub_FACTOR < 10) {
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/10.0) {
						sub_FACTOR = 10 ; sub_N_MAX = si ;		
					}
				}
				else if (Dt_local >= Dt * 5 && sub_FACTOR < 5) {
					if (si*1.0/N_MAX < 1.0/sub_FACTOR*(sub_N_MAX*1.0/N_MAX*sub_FACTOR + 1.0) - 1.0/5.0) {
						sub_FACTOR = 5 ; sub_N_MAX = si ;		
					}
				}
				else if (Dt_local >= Dt * 2 && sub_FACTOR < 2 && si < N_MAX/2) {
					sub_FACTOR = 2 ; sub_N_MAX = si ;		
				}
				else if (si > N_MAX/2) {	
					break ;
				}
			}

			sub_rmax = sr[sindex[sub_N_MAX]] ;
			sub_count = 1 ;		
			sub_totaltime = Dt ; 

		}

		else {							/** core timestep only **/

			sub_count++ ;			
			sub_totaltime += Dt ;	

			if (sub_count == sub_FACTOR) {	/* last round, so do a FULL time step. */
				sub_count = 0 ;
			}
		}


	/*************** end setup sub time step ***************************/




	/*	Perturb velocities of all N_MAX stars. Using sr[], sv[], get NEW E, J for
		all stars
	*/
		if (PERTURB > 0) {
			perturb_stars(Dt) ;
		}

		/* Set Real (physical) time in yr for stellar evolution and elliptical orbits */

		RealTime = T_RELAX * TotalTime ;


		if (STELLAR_EVOLUTION > 0) {
			DoStellarEvolution() ;
		}

	/*	Computing New positions and velocities for stars from sE[], sJ[]
		and sgravity[]. Returns max radius for all stars.
	*/
		j = 0 ;
		Etidal = 0.0 ;
		OldTidalMassLoss = TidalMassLoss ;
		max_r = get_positions() ;
		DTidalMassLoss = TidalMassLoss - OldTidalMassLoss ;

		fprintf(stderr,"Iteration %d: OldTidalMassLoss = %.6g   DTidalMassLoss = %.6G \n", \
			j, OldTidalMassLoss, DTidalMassLoss) ;
		fprintf(logfile,"Iteration %d: OldTidalMassLoss = %.6g   DTidalMassLoss = %.6G \n", \
			j, OldTidalMassLoss, DTidalMassLoss) ;

	/* Iterate the removal of tidally stripped stars by reducing Rtidal */
		do {
			Rtidal = orbit_r * pow(Mtotal[0]-(TidalMassLoss-OldTidalMassLoss), 1.0/3.0) ;
			phi_rtidal = potential(Rtidal) ;
			phi_zero = potential(0.0) ;
			DTidalMassLoss = 0.0 ;
			for (i = 1 ; i <= N_MAX ; i++) {
				if (sr_apo[sindex[i]] > Rtidal && srnew[sindex[i]] < 1000000) {	
					srnew[sindex[i]] = INFINITY ;		/* tidally stripped star */
					svrnew[sindex[i]] = 0.0 ;
					svtnew[sindex[i]] = 0.0 ;
					Eescaped += sE[sindex[i]] * sm[sindex[i]] / N_STAR ;
					Jescaped += sJ[sindex[i]] * sm[sindex[i]] / N_STAR ;
					DTidalMassLoss += sm[sindex[i]] / N_STAR ;
					Etidal += sE[sindex[i]] * sm[sindex[i]] / N_STAR ;
					fprintf(escfile, \
						"%ld  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g \n", \
						tcount, TotalTime, sm[sindex[i]], sInitialMass[sindex[i]], \
						sr[sindex[i]], svr[sindex[i]], svt[sindex[i]], sr_peri[sindex[i]], \
						sr_apo[sindex[i]], Rtidal, phi_rtidal, phi_zero, \
						sE[sindex[i]], sJ[sindex[i]]) ;

					if (Etotal - Etidal >= 0) break ;
				}
			}
			j++ ;
			TidalMassLoss = TidalMassLoss + DTidalMassLoss ;
			fprintf(stderr,"Iteration %d: TidalMassLoss = %.6g   DTidalMassLoss = %.6G \n", \
				j, TidalMassLoss, DTidalMassLoss) ;
			fprintf(logfile,"Iteration %d: TidalMassLoss = %.6g   DTidalMassLoss = %.6G \n", \
				j, TidalMassLoss, DTidalMassLoss) ;
		}
		while (DTidalMassLoss > 0 && (Etotal - Etidal) < 0) ;

		Prev_Dt = Dt ;

	/* Compute Intermediate Energies of stars. Also transfers new positions and
	   velocities from srnew[], svrnew[], svtnew[] to sr[], svr[], svt[], and saves srOld[] 
	*/
		ComputeIntermediateEnergy() ;


	/*	Indexing stars by radius. The 0th star at radius 0 
		and (N_STAR+1)th star at INFINITY are already set earlier.
	*/
		indexx(N_STAR_NEW, sr, sindex) ;
		sindex[N_STAR_NEW+1] = N_STAR_NEW+1 ;

	/* Count mass in each component 
	*/
		CountComponentMass() ;	


	/*  For elliptical orbit, compute phase, and get new orbit_r
	*/
		if (ORBIT_ELLIPTICAL == 1) {
			mean_anomaly = 2*pi* TotalTime*T_RELAX/ORBIT_PERIOD ;
			/* get mean anomaly modulo 2*pi */
			mean_anomaly -= floor(mean_anomaly/(2*pi))*2*pi ;
			/* Subtract initial phase from mean anomaly 
			mean_anomaly -= (pi/2.0 - ORBIT_ECCENTRICITY) ; */
			ecc_anomaly = rtbis(kepler_anomaly, -pi/2.0, 2*pi, 1.0e-6) ;
			orbit_r = R_MAX / (1.0 - ORBIT_ECCENTRICITY) * (1.0 - ORBIT_ECCENTRICITY*cos(ecc_anomaly)) ;
		}
		else if (ORBIT_ELLIPTICAL == 2) {
			mean_anomaly = RealTime ;
			/* This is the evolution time... reduce it modulo the Orbital Period */
			mean_anomaly -= floor(mean_anomaly/ORBIT_PERIOD)*ORBIT_PERIOD ;
			/* Get time in orbital time units */
			mean_anomaly = mean_anomaly / ORBIT_PERIOD * orbit_t[orbit_step_count] ;
			for (i=0 ; i < orbit_step_count ; i++) {
				if (orbit_t[i] >= mean_anomaly) break ;
			}
			orbit_r = R_MAX * orbit_p[i] / pow(orbit_Mg[i], 1.0/3.0) ;
		}
		else {
			orbit_r = R_MAX ;
		}
	
	/*	Computing the gravity at the star locations provided in sr[].
		Results returned in sgravity[]. Computes Mtotal[0] (total mass).
		Returns N_MAX. 
		This is the	only place where N_MAX can change. (Also in phi(),
		but that function is only used as a check, not used in the regular 
		computation.

		Also computes new Rtidal using Mtotal[0] and orbit_r
	*/
		errstat = gravity() ;


	/* Get new sindex[] positions of TRACE stars 
	*/
		if (TRACE_STARS > 0) {
			for (i = 1 ; i <= TRACE_N ; i++) {
				for (j = 1 ; j <= N_MAX ; j++) {
					if (sindex[j] == i) {
						trace[i] = j ;
						break ;
					}
				}
			}
		}

	/* Compute Lagrangian Radii 
	
		ComputeLagrangianRadii() ;
	*/

	/* Recompute Energy. Uses the specified ECONS_MODE 
	*/
		RecomputeEnergy() ;


	/* Compute KE's for each component, etc. 
	
	   ComputeComponentKE() ;
	*/


/*	-------------- Printing Results --------------- */

	/* save distribution of trace stars */

	if (TRACE_STARS > 0) {
	for (i = 1 ; i <= TRACE_N ; i++) {

	/* output positions of ALL trace stars (trace star snapshot) every 50 timesteps */

		if (tcount % 100 == 0) {
			fprintf(tracefile2, \
				"%5d  %.8g  %5.3g %5.3g  %5.3g  %5.3g  %5.3g  %5.3g  %5.3g  %5.3g\n", \
				i, TotalTime, sr[i], sr_peri[i], sr_apo[i], sl[i], sE[i], sJ[i], \
				sDE[i], sDJ[i]) ;
		}

		/* output t_mass_seg for each trace star, when it enters the core */
		/* Check if trace star orbit is ENTIRELY within the core */

		if (trace_stat[i] == 0 && sr_apo[i] < core_radius) {
			trace_stat[i] = 1 ;
			fprintf(tracefile3, \
				"%5d  %5d  %.8g  %5.3g %5.3g  %5.3g  %5.3g  %5.3g  %5.3g  %5.3g\n", \
				i, tcount, TotalTime, sr[i], sr_peri[i], sr_apo[i], sl[i], \
				trace_r_ini[i], trace_tr[i], core_radius) ;
		}
	}
	}

	/* Output stellar evolution data to out_stellar file */

	fprintf(stellarfile, \
		"%ld  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %ld  %ld  %ld  %ld  %.8g  %ld  %.8g\n", \
		tcount, TotalTime, RealTime, mass_1, mass_2, mass_3, mass_4, \
		count_1, count_2, count_3, count_4, turnoff_mass, count_5, mass_5) ;


	/* Output binary data 
	   Note: N_BINARY counts ALL binaries (including escaped/destroyed ones)
	   whereas count_5 only counts EXISTING BOUND binaries. 
	   M_binary = mass_5 (always).
	*/

	if (N_BINARY > 0) {
		fprintf(binaryfile, \
		"%ld  %.8g  %ld  %ld  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %ld  %ld\n", \
		tcount, TotalTime, N_BINARY, count_5, M_binary, E_binding, rh_single, \
		rh_binary, rho_core_single, rho_core_bin, Delta_BE_bb, Delta_BE_bs, DE_bb, DE_bs, N_bb, N_bs) ;
	}

	PrintLogOutput() ;

	PrintFileOutput() ;

	fflush(NULL) ;

}     /* End FOR (time step iteration loop) */

if(printstatus==0 || printstatus==6)  
	fprintf(stderr,"Total time = %.8G seconds", \
		(double) (clock() - timecount)/CLOCKS_PER_SEC) ;


if(printstatus==0) 
	for(i=1;i<=5;i++) fclose(out[i]);

fclose(logfile) ;

if (TRACE_STARS > 0 || N_BINARY > 0) {
	fclose(tracefile2) ;
	fclose(tracefile3) ;
}

if (N_BINARY > 0) {
	fclose(binaryfile) ;
}

fclose(escfile) ;
fclose(stellarfile) ;

#if defined(MPIMODE)
	MPI_Finalize();     /* Closes MPI communications */ 
#endif

}	


/************************ end main **********************/





/***** Get conversion factors from code units to physical units ***/


void get_physical_units(void) {

/* Input FAMILY (input file), actual_total_mass, N_STAR, Rtidal (inputfile) */

  double cluster_mass, Rg_kpc, Rg_cm, Mg_g, Mg_sm ;
  double rt_pc, rt_rvir, minus_4Eo ;

  rt_rvir = R_MAX ;		/* initial Rtidal = RMAX specified in input file */
  cluster_mass = actual_total_mass ;

  Rg_kpc = FAMILY * log(N_STAR) / (cluster_mass);
  fprintf(stderr,"\nRg is %g kpc\n", Rg_kpc);
  fprintf(logfile,"\nRg is %g kpc\n", Rg_kpc);

  Rg_cm = Rg_kpc * KPC;

  Mg_g=pow(220e5, 2.0)*Rg_cm/G;
  Mg_sm=Mg_g/SOLAR_MASS;
  fprintf(stderr,"M_galaxy in solar masses is %g sm\n", Mg_sm);
  fprintf(logfile,"M_galaxy in solar masses is %g sm\n", Mg_sm);

  rt_pc= Rg_kpc * 1000 * pow((cluster_mass/(3*Mg_sm)) , (1.0/3.0) );
  UL=rt_pc/rt_rvir;  /* unit of distance in pc */
  fprintf(stderr,"Rtidal = %g pc ; [UL] = %g pc\n", rt_pc, UL);
  fprintf(logfile,"Rtidal = %g pc ; [UL] = %g pc\n", rt_pc, UL);
    
  minus_4Eo=pow( (UL / 1000 * KPC / G / pow((cluster_mass * SOLAR_MASS), 2.0)), -1.0);
  fprintf(stderr,"-4Eo = %g\n", minus_4Eo);
  fprintf(logfile,"-4Eo = %g\n", minus_4Eo);

  /* Here assuming gamma = 1 so that T_RELAX = UT * N/ln(N)  */
  UT = G * pow((cluster_mass * SOLAR_MASS), 2.5) * pow(minus_4Eo, -1.5) / YEAR;
  T_RELAX= UT * N_STAR / log(N_STAR);

  fprintf(stderr,"[UT] = %g yr ;  T_RELAX  %g yr\n\n", UT, T_RELAX);
  fprintf(logfile,"[UT] = %g yr ;  T_RELAX  %g yr\n\n", UT, T_RELAX);
  

} 




/*********** Output 2D/3D snapshots **************/

void print_2Dsnapshot(void) {

	double x, y, z, mass, theta, phi, rxy ;
	long i, j ;
	static long snap_num=0 ;
	static char outfile[100] ;

	snap_num++ ;	/* persistent counter for _snapXX snapshot output file */

	/* open file for 2D snapshot */
	sprintf(outfile,"%ssnap%ld", outprefix, snap_num);
	if ((snapfile = fopen(outfile, "w+t"))==NULL){
	  printf("\n\nCannot create 2D snapshot file %s\n", outfile);
	}

	for (i=1; i<=N_MAX ; i++) {
		j = sindex[i] ;
		phi = ran2(&idum) * 2.0*pi ;
		theta = acos(ran2(&idum)*2.0 - 1.0) ;
		x = sr[j] * sin(theta) * cos(phi) ;
		y = sr[j] * sin(theta) * sin(phi) ;
		z = sr[j] * cos(theta) ;
		rxy = sqrt(x*x + y*y) ;

		/* mass is output in SOLAR MASSES */
		mass = sm[j] * actual_total_mass / N_STAR ;

		fprintf(snapfile, \
			"%ld  %ld  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %ld\n", \
			tcount, i, sr[j], rxy, x, y, z, mass, svr[j], svt[j], stype[j]) ;

	}
	
	fprintf(logfile, "\n\n2D snapshot stored... t = %.8G\n", TotalTime) ;
	fprintf(stderr, "\n\n2D snapshot stored... t = %.8G\n", TotalTime) ;


	fclose(snapfile) ;
}



/****************** Binary routines *********************/

void assign_binaries(void) {

	long i, j ;
	double min_bindingenergy, max_bindingenergy, ba_min, ba_max ;


/**** assign the first N_BINARY stars (NOT sorted by radius- hence 
	  distributed according to density profile) as binaries, and select
	  a, e, m1, m2, for them. Also, update sm = m1 + m2

	Adding mass due to secondary stars causes the TOTAL cluster mass to 
	increase. Hence the average mass also changes. Hence ALL masses must be
	changed to the new unit of average mass.
*/

	/*	First calculate velocity dispersion in the core to get 
		Min binding energy G*m1*m2/2a = m v_core^2 
		Note: v_core is a global variable updated every timestep in
		GetTimeStep().
	*/  

    v_core=0.0;
	for(i=1; i <= NUM_CORE_STARS; i++) {
		v_core += pow(svr[sindex[i]],2.0) + pow(svt[sindex[i]],2.0); 
	}
	v_core = sqrt(v_core / NUM_CORE_STARS);

	/* [BE_min] = [<m> * v_c^2] = v_c^2/N_STAR, with <m>=1 in avg mass units */
	min_bindingenergy = 1.0 * v_core*v_core/N_STAR ;		
	max_bindingenergy = 133 * v_core*v_core/N_STAR ;

	/* Since m1=m2 = <m> (approx) and [G]=1, [BE] = [G*m1*m2/2a] = 1/(2a*N_STAR^2) 
	    => a_max = 0.5 / min_bindingenergy 
	*/
	ba_max = (0.5 / min_bindingenergy)/N_STAR/N_STAR ;
	ba_max = log10(ba_max) ;
	ba_min = (0.5 / max_bindingenergy)/N_STAR/N_STAR ;
	ba_min = log10(ba_min) ;

	/* Begin assigning binary parameters */

	M_binary = 0.0 ;			/* Total mass in binaries */
	sec_binary_mass = 0.0 ;		/* Total mass due to secondary stars */

	for (i=1 ; i<=N_BINARY ; i++) {
		
		bm1[i] = sm[i] * actual_total_mass / N_STAR ;	/* primary mass in SOLAR MASSES */
		bm2[i] = sm[i] * actual_total_mass / N_STAR ;	/* secondary mass in SOLAR MASSES */

		for (j=1; j<=N_TRY; j++) {			/* eccentricity distribution P(e) = 2*e */
			be[i] = ran2(&idum) ;
			if (2.1*ran2(&idum) < 2.0*be[i]) break ;
		}

		/* select log(a) uniform between ba_min and ba_max */
		ba[i] = ran2(&idum)*(ba_max-ba_min) + ba_min ;	
		ba[i] = pow(10.0, ba[i]) ;			


	/* Save all INITIAL  binary parameters */

		bini_m1[i] = bm1[i] ;
		bini_m2[i] = bm2[i] ;
		bini_a[i] = ba[i] ;
		bini_e[i] = be[i] ;
		bini_tformed[i] = TotalTime ;
		bini_tdestroyed[i] = 0.0 ;


		by[i] = 0.0 ;			/* change in binding energy delta_E_bin */

		b_sindex[i] = i ;		/* index from binary into list of stars */
		s_binindex[i] = i ;		/* index from star into list of binaries */
		stype[i] = 5 ;			/* Set star type = 5 (binary) */

		sm[i] = bm1[i] + bm2[i] ;		/* update the total binary mass */
		M_binary += sm[i] ;				/* total mass of binaries in SOLAR MASSES */
		sec_binary_mass += bm2[i] ;		/* total sec mass added due to secondaries */

	}

	/* convert masses of remaining single stars to NEW UNIT of AVERAGE MASS */

	for (i=N_BINARY+1 ; i<=N_STAR ; i++) {
		sm[i] *= actual_total_mass / (actual_total_mass + sec_binary_mass) ;
		s_binindex[i] = 0 ;		/* NOT a binary */
		stype[i] = 1 ;			/* MS star */
	}

	actual_total_mass += sec_binary_mass ;	/* NEW actual total mass in SOLAR MASSES */

	
	/* convert all binary masses from solar masses to NEW AVG MASS units */
	
	E_binding = 0.0 ; DE_bb = DE_bs = 0.0 ;

	for (i=1 ; i<=N_BINARY ; i++) {

		sm[i] = sm[i] / actual_total_mass * N_STAR ;	/* convert mass to NEW AVERAGE MASS units */
		bm1[i] = bm1[i] / actual_total_mass * N_STAR ;	/* primary mass in SOLAR MASSES */
		bm2[i] = bm2[i] / actual_total_mass * N_STAR ;	/* secondary mass in SOLAR MASSES */

		bini_m1[i] = bm1[i] ;
		bini_m2[i] = bm2[i] ;
		
		E_binding += bm1[i]*bm2[i]/N_STAR/N_STAR/2.0/ba[i] ;
	}

	M_binary = M_binary / actual_total_mass * N_STAR ;
	sec_binary_mass = sec_binary_mass / actual_total_mass * N_STAR ;

	/* TRACE_STARS = 1 ;		 treat binaries as trace stars 
	*/
	TRACE_N = N_BINARY ;
	count_1 = N_STAR - TRACE_N ;
	count_5 = TRACE_N ;
	mass_1 = N_STAR - M_binary ;	/* note: total mass is N_STAR in avg mass units */
	mass_5 = M_binary ;

}





/**************** Integration by Quadrature **************/

double log_potential(double rd) {

	return ( 1.0 / sqrt( (1.0+ORBIT_ECCENTRICITY) + 2.0*log(2.0) - 2.0*log(1.0+rd*rd) \
			- (1.0+ORBIT_ECCENTRICITY)/(rd*rd) ) ) ;

}

double log_m(double rd) {

	return ( 2.0* rd*rd * (3.0 + rd*rd) / pow(1.0 + rd*rd, 2.0) ) ;

}


#define EPS 1.0e-6
#define JMAX 40

double qsimp(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, long n);
	void nrerror(char error_text[]);
	long j;
	double s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}
#undef EPS
#undef JMAX

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, long n)
{
	double x,tnm,sum,del;
	static double s;
	long it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

/************************************************************/


void RecomputeEnergy(void) {


static double dtemp, dtemp2 ;
static long k, i ;

/* Recalculating Energies */
	Etotal = 0.0 ;	
	KEtotal = 0.0 ;
	PEtotal = 0.0 ;
	n_extreme = 0 ; n_orbit = 0 ; dtemp2 = 0 ; 
	avg_X = 0 ; avg_r = 0 ; avg_dv = 0 ; 

	if (E_CONS == 0) {		/* recompute new sE[] and sJ[], using the new potential */	

		for (i = 1 ; i <= N_MAX ; i++) {
			k = sindex[i] ;
			sE[k] = sgravity[k] + 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) ;
	
		/*	sE[k] = (-sMi[k]/sr[k] - sMi_ri[k]) + 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) ;
		*/	
			sJ[k] = sr[k] * svt[k] ;

			KEtotal += 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) * sm[k] / N_STAR ;

		/*	Compute PE using Henon method using sgravity[]   */
			PEtotal += sgravity[k] * sm[k] / N_STAR ;

		/*	Compute PE using sMi_Ri[] and sMi[]
			PEtotal += (-sMi[k]/sr[k] - sMi_ri[k]) * sm[k] / N_STAR ;		
		*/
		/*	--- Compute Etotal by summing sE[] to ensure E conservation
			Etotal += sE[k] * sm[k] / N_STAR;
		*/	
			sinteracted[k] = 0 ;	/* reset sinteracted[] flag to 0 */
		}
	}
	else {	/* Try to conserve energy by using intermediate potential */

		for (i = 1 ; i <= N_MAX ; i++) {
			k = sindex[i] ;
		/*	if (potential(srOld[k]) == 1) {
				fprintf(stderr,"Error: Escaped star: i = %5ld  k = %5ld  r = %.8G\n",i,k,srOld[k]) ;
				if(printstatus==6) gets(&cont) ;
			}
		*/

			/* Note: svt[] = J/r_new is already computed in get_positions() */
			/* ignore stars near pericenter, and those with strong interactions */

			if (sX[k] > 0.05 && sinteracted[k] == 1) {
				dtemp2 = potential(srOld[k]) ;
				dtemp = sEI[k] - sgravity[k] + dtemp2  ;

				/* svt[k] = sJI[k]/sr[k] ; ---- Already computed in get_positions() */

				if (dtemp - svr[k]*svr[k] > 0) {
					svt[k] = sqrt(dtemp - svr[k]*svr[k]) ;  /* preserve svr[k] and change svt[k] */

				/* PRESERVE svt[k], and change svr[k] ONLY while preserving its sign 
					if (svr[k] > 0) 
						svr[k] = sqrt(dtemp - svt[k]*svt[k]) ;	
					else 
						svr[k] = -1.0 * sqrt(dtemp - svt[k]*svt[k]) ;
				*/
				}
				else {
					if (dtemp > 0) {	
						avg_dv += dtemp - (svt[k]*svt[k] + svr[k]*svr[k]) ;

						if (dtemp > svt[k]*svt[k]) {
							svr[k] = sqrt(dtemp - svt[k]*svt[k]) ;
						}
						else {
							svt[k] = sqrt(dtemp/2.0) ;
							svr[k] = sqrt(dtemp/2.0) ;
						}
						n_orbit++ ;
						avg_r += sr[k] ;
						avg_X += sqrt((sX[k] - 0.5)*(sX[k] - 0.5)) ;
						avg_DE += fabs(sDE[k]) ;
						avg_DJ += fabs(sDJ[k]) ;

					/*	svt[k] = sqrt(dtemp) ;						
						n_circular++ ;
						printf("error: i = %5d n_circular = %6d v_new^2 = %.9g   J = %.9g  svt = %.9g\n", \
							i, n_circular, dtemp, sJI[k], svt[k]) ;
					*/	
					/*	 reduce the energy of the next star to compensate 
						if (i < N_MAX) 
							sEI[sindex[i+1]] += dtemp - svt[k]*svt[k] ;
					*/	
					}
					else {
						avg_dv += dtemp - (svt[k]*svt[k] + svr[k]*svr[k]) ;
						
						/* velocity dtemp and v^2 are already small, so keep unchanged 
						svr[k] = 0 ;
						*/
						n_extreme++ ;
						avg_r += sr[k] ;
						avg_X += sqrt((sX[k] - 0.5)*(sX[k] - 0.5)) ;
						avg_DE += fabs(sDE[k]) ;
						avg_DJ += fabs(sDJ[k]) ;

					/*	printf("error: i = %5d n_circular = %6d v_new^2 = %.9g   J = %.9g  svt = %.9g\n", \
							i, n_circular, dtemp, sJI[k], svt[k]) ;
					*/	
					/* reduce the energy of the next star to compensate 
					*/	if (i < N_MAX) 
							sEI[sindex[i+1]] += dtemp - (svt[k]*svt[k] + svr[k]*svr[k]);
					}
				}
			
			}

		/* recompute new sE[] and sJ[], using the new potential
		*/	 
			sE[k] = sgravity[k] + 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) ;
			sJ[k] = sr[k] * svt[k] ;

			KEtotal += 0.5*(svr[k]*svr[k] + svt[k]*svt[k]) * sm[k] / N_STAR ;
			PEtotal += sgravity[k] * sm[k] / N_STAR ;

		/*	if ((printstatus==6 || printstatus==0) && tcount == 102 && (i == 6201 || i==6202))
				fprintf(stderr, \
				"ERROR: i = %5d  v_new^2 = %.6g  E = %.6G  J = %.6g  r = %.6G  svr = %.6g  svt = %.6g\n", \
				i, dtemp, sE[k], sJ[k], sr[k], svr[k], svt[k]) ;
		*/

			sinteracted[k] = 0 ;	/* reset sinteracted[] flag to 0 */

		}

	}	/* End ELSE - conserve energy */


	/*	Etotal *= 1.0/3.0 ;	  
	*/
		PEtotal = PEtotal / 2.0 ;
	/*	--- Computing Etotal using gravity (PE+KE)  */
		Etotal = KEtotal + PEtotal ;

		avg_r = avg_r / (n_extreme + n_orbit) ; /* <r> for n_orbit and n_circular cases */
		avg_X = avg_X / (n_extreme + n_orbit) ;
		avg_DE = avg_DE / (n_extreme + n_orbit) ;	
		avg_dv = avg_dv / (n_extreme + n_orbit) ;
		avg_DJ = avg_DJ / (n_extreme + n_orbit) ;

	/* renormalize the total energy

		Dnorm = fabs(Etotal_New / Etotal) ;
		for (i = 1 ; i <= N_MAX ; i++) {
			sE[sindex[i]] *= Dnorm ;
			 v = sqrt(2.0*(sE[sindex[i]] - sgravity[sindex[i]])) \
				/ sqrt(svr[sindex[i]]*svr[sindex[i]] + svt[sindex[i]]*svt[sindex[i]]) ;
			svr[sindex[i]] *= v ;
			svt[sindex[i]] *= v ;
		}
		Etotal = Etotal_New ;
	*/

}


void ComputeIntermediateEnergy(void) {

	static long i, k, j ;

	/* compute intermediate energies for stars due to change in potential 
	*/	for (i = 1 ; i <= N_MAX_NEW ; i++) {

			j = k = sindex[i] ;
			/* new star from destroyed binary -- refer back to original binary for old params */
			if (i > N_MAX) k = sindex[old_k[j]] ;
			
			if (srnew[j] < 1000000) {	/* do only for NON-Escaped stars */
				sEI[j] = (svr[k]*svr[k]+svt[k]*svt[k] + sgravity[k] - potential(srnew[j]));
				sJI[j] = svt[k] * sr[k] ;
			}

		}

	/*	Transferring new positions to sr[], svr[] and svt[] arrays 
		from srnew[], svrnew[] and svtnew[] arrays
	*/	
		for (i = 1 ; i <= N_MAX_NEW ; i++) {
			j = k = sindex[i] ;
			/* new star from destroyed binary -- refer back to original binary for old params */
			if (i > N_MAX) k = sindex[old_k[j]] ;
			srOld[j] = sr[k] ;
			sr[j] = srnew[j] ;
			svr[j] = svrnew[k] ;
			svt[j] = svtnew[k] ;
		}

}

long CheckStop(void) {

	static long Echeck ;

	if (tcount >= T_MAX_COUNT) { 
	   if(DUMPS==1) output_rawdump(); 
	   if(DUMPS==1) print_2Dsnapshot(); 
		fprintf(stderr,"\n\nNo. of timesteps > T_MAX_COUNT ... Terminating.\n\n"); 
		fprintf(logfile,"\n\nNo. of timesteps > T_MAX_COUNT ... Terminating.\n\n"); 
	   return (1) ;
	}

	if (TotalTime >= T_MAX) { 
	   if(DUMPS==1) output_rawdump(); 
	   if(DUMPS==1) print_2Dsnapshot(); 
		fprintf(stderr,"\n\nTotalTime > T_MAX ... Terminating.\n\n"); 
		fprintf(logfile,"\n\nTotalTime > T_MAX ... Terminating.\n\n"); 
	   return (1) ;
	}

	/* Stop if cluster is disrupted -- N_MAX is too small */
	if (N_MAX < (0.02*N_STAR)) {
	  if(DUMPS==1) output_rawdump();
	  if(DUMPS==1) print_2Dsnapshot(); 
		fprintf(stderr,"\n\nN_MAX < 0.02 * N_STAR ... Terminating.\n\n"); 
		fprintf(logfile,"\n\nN_MAX < 0.02 * N_STAR ... Terminating.\n\n"); 
	  return (1);
	}

	/* Stop of Etotal > 0 */
	if (Etotal > 0.0) {
	  if(DUMPS==1) output_rawdump();
	  if(DUMPS==1) print_2Dsnapshot(); 
		fprintf(stderr,"\n\nEtotal > 0 ... Terminating.\n\n"); 
		fprintf(logfile,"\n\nEtotal > 0 ... Terminating.\n\n"); 
	  return(1);
	}


/*	If inner-most Lagrangian radius is too small, then stop: 
*/
		if ( mass_r[0] < MIN_LAGRANGIAN_RADIUS) {
			if(DUMPS==1) output_rawdump();
			if(DUMPS==1) print_2Dsnapshot(); 
			fprintf(stderr,\
				"\n\nMin Lagrange radius < %.6G ... Terminating.\n\n", \
				MIN_LAGRANGIAN_RADIUS); 
			fprintf(logfile, \
				"\n\nMin Lagrange radius < %.6G ... Terminating.\n\n", \
				MIN_LAGRANGIAN_RADIUS); 
			return (1) ;
		} 

/*      Output some snapshots near core collapse (if core density is high enough)
*/
		if (rho_core > 50.0 && Echeck==0) { 
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
        }
        else if (rho_core > 1.0e2 && Echeck==1) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }
        else if (rho_core > 5.0e2 && Echeck==2) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 1.0e3 && Echeck==3) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 5.0e3 && Echeck==4) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 1.0e4 && Echeck==5) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 5.0e4 && Echeck==6) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 1.0e5 && Echeck==7) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 5.0e5 && Echeck==8) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          
        else if (rho_core > 1.0e6 && Echeck==9) {
			if(DUMPS==1) {
/*				output_rawdump();
*/				print_2Dsnapshot();
			}
			Echeck++;
	    }          


	/* If total Energy has drifted too far, then output_rawdump() once. 

		if( Etotal < Etotal_initial - DUMP_ENERGY_DISPLACEMENT && Echeck!=1) {
		  if(DUMPS==1) output_rawdump();
		  if(DUMPS==1) print_2Dsnapshot(); 
		  Echeck=1;
		} 
	*/	
        
	/* If total Energy has diminished by TERMINAL_ENERGY_DISPLACEMENT, then stop */

		if( Etotal < Etotal_initial - TERMINAL_ENERGY_DISPLACEMENT ) {
		  if(DUMPS==1) output_rawdump();
		  if(DUMPS==1) print_2Dsnapshot(); 
		  fprintf(stderr,"\n\nTerminal Energy reached... Terminating.\n\n"); 
		  fprintf(logfile,"\n\nTerminal Energy reached... Terminating.\n\n"); 
		  return (1) ;  
		}

		return (0) ; /* NOT stopping time yet */

}


double GetTimeStep(void) {

	static double m_avg, w2_avg, Ai, zr_min, zr_max ;
	static double mtemp, m_single, m_bin ;
	static long i, zk, k ;

/*	Calculation of Relaxation time in the core in order to compute Dt    */  

    mtemp = m_single = m_bin = 0.0;
    v_core=0.0; rho_core_single = rho_core_bin = 0.0 ;

	for(i=1; i<=NUM_CORE_STARS; i++) {
		k = sindex[i] ;
		mtemp += sm[k];
		if (stype[k] != 5) m_single += sm[k] ;
		else m_bin += sm[k] ;
		v_core += pow(svr[sindex[i]],2.0) + pow(svt[sindex[i]],2.0); 
	}
	rho_core = mtemp/ N_STAR / (4.0/3.0*pi*(pow(sr[sindex[NUM_CORE_STARS+1]],3.0)));
	v_core = sqrt(v_core / NUM_CORE_STARS);

	/* compute binary and single star densities in the core */

	rho_core_single = m_single/ N_STAR / (4.0/3.0*pi*(pow(sr[sindex[NUM_CORE_STARS+1]],3.0)));
	rho_core_bin = m_bin/ N_STAR / (4.0/3.0*pi*(pow(sr[sindex[NUM_CORE_STARS+1]],3.0)));

	Trc = (0.065 * pow(v_core, 3.0) / rho_core) / (mtemp/NUM_CORE_STARS) / (log(N_MAX)/log(N_STAR)); 


	/********* Computing the core radius ********/

	core_radius = sqrt(3.0 * v_core*v_core / (4*pi * rho_core)) ;


	/* calculate Ai for innermost zone
	*/	
	m_avg = 0.0 ; w2_avg = 0.0 ; zk = NUM_CORE_STARS ;
	for (i = 1 ; i <= zk ; i++) {
		k = sindex[i] ;
		m_avg += sm[k] ;
		w2_avg += sm[k] * (svr[k]*svr[k]+svt[k]*svt[k]) ;
	}
	m_avg = m_avg / zk ;
	w2_avg = w2_avg * 2.0 / m_avg / zk ;
	zr_min = sr[sindex[1]] ; zr_max = sr[sindex[zk]];
	Ai = 6.0 * zk * m_avg*m_avg / \
		(zr_max*zr_max*zr_max - zr_min*zr_min*zr_min) / pow(w2_avg, 1.5) ;


	if (DT_MODE == 0) {			/* Compute Dt using henon's method */
		Dt = TRC_FACTOR * (double) N_MAX / (double) N_STAR * pow(Mtotal[0], 2.5) / \
			pow(fabs(Etotal), 1.5) ;
	}
	else if (DT_MODE == 1) {	/* Compute Dt = TRC_FACTOR * T_rc */
		Dt = TRC_FACTOR * Trc; 
	}	
	else if (DT_MODE == 2) {	/* Compute Dt s.t. <sin2beta> < SIN2BETA_MAX in the core */
		Dt = DT_MAX ;
		Sin2Beta = Ai * Dt / N_STAR ;
		while (Sin2Beta > SIN2BETA_MAX) {
			Dt = Dt * 0.8 ;
			Sin2Beta = Ai * Dt / N_STAR ;
		}
	}
	Sin2Beta = Ai * Dt / N_STAR ;

/*	Enforce Minimum allowed size of Dt --- DT_MIN
*/
	if ( Dt < DT_MIN ) Dt = DT_MIN;   

	
/*  For an elliptical orbit, check that the timestep is at 
	most ORBIT_MAX_DT * ORBIT_PERIOD
*/
	if (ORBIT_ELLIPTICAL > 0) {
		if (Dt * T_RELAX > ORBIT_PERIOD * ORBIT_MAX_DT) {
			Dt = ORBIT_PERIOD * ORBIT_MAX_DT / T_RELAX ;
		}
	}

	return (Dt) ;


}

void CheckStellarTimeStep(void) {
  double tempRealTime;  
  long i=0, k ;
  double m_real, m_final;
  
  double projectedmassloss = 1.0;

  /* so can write only one line at end of while loop */
  Dt/=0.8;
  
  /* Do NOT change masses of escaped stars... need original mass for RESTART */

  while(projectedmassloss > 0.01) {

    projectedmassloss=0.0;
    tempRealTime = T_RELAX * (TotalTime+Dt) ;

    fprintf(stderr,"N_MAX is %ld and tempRealTime is %g\n", N_MAX, tempRealTime);
    fprintf(logfile,"N_MAX is %ld and tempRealTime is %g\n", N_MAX, tempRealTime);

    for (k=1; k<=N_MAX ; k++) {
		i = sindex[k] ;
		m_real = sm[i] * actual_total_mass / N_STAR ;
		m_final = sm[i];

		if(STELLAR_EVOLUTION==1) {
			if ((tempRealTime > sLifeTime[i]) && ((sInitialMass[i]-sm[i]) < 1.0e-10)) {
			  if (m_real < 4.7) {
				m_final = (0.58 + 0.22 * (m_real - 1.0)) * (N_STAR/actual_total_mass) ;
			  }
			  else if ((m_real >= 4.7) && (m_real < 8.0)) {
				m_final = 0.0 ;
			  }
			  else if (m_real <= 15.0) {
				m_final = 1.4 * (N_STAR/actual_total_mass) ;
			  }
			}
		}
		else if(STELLAR_EVOLUTION==2) {
			if ((tempRealTime > sLifeTime[i]) && ((sInitialMass[i]-sm[i]) < 1.0e-10)) {
			  m_final=sRemnantMass[i] * (N_STAR/actual_total_mass);
			}
		}
		else
			fprintf(stderr, "ERROR: Wrong input for STELLAR_EVOLUTION variable\n");
			  projectedmassloss += ( (sm[i]-m_final) / N_STAR);
		}

		printf("projected mass loss is %g \n", projectedmassloss);
		Dt*=0.8;
	}

}



void PrintLogOutput(void) {

	char cont ;
	double mh, rh, trh, conc_param, m_single, m_binary ;
	long ih, rhcount, k ;

/* Computing half-mass radius, and relaxation time */

	mh = rh = trh = 0.0 ;
	rh_binary = rh_single = m_binary = m_single = 0.0 ;
	rhcount = 0 ;

	if (count_5 > 0) {		/* if there are any binaries */

		for(ih=1; ih < N_MAX; ih++) {
			k = sindex[ih] ;
			mh+=sm[k]/N_STAR ;
			if (stype[k] != 5) 
				m_single += sm[k]/N_STAR ;
			else  
				m_binary += sm[k]/N_STAR ;
			if (rh < ZERO && (mh/Mtotal[0]) > 0.5) {
				rhcount++ ;
				rh = sr[k] ;
			}
			if (rh_single < ZERO && m_single/(Mtotal[0]-(mass_5/N_STAR)) > 0.5) {
				rhcount++ ;
				rh_single = sr[k] ;
			}
			if (rh_binary < ZERO && m_binary/mass_5*N_STAR > 0.5) {
				rhcount++ ;
				rh_binary = sr[k] ;
			}

			if (rhcount >= 3) break ;
		}
	}
	else {				/* No binaries ... so only compute rh */

		for(ih=1; ih < N_MAX; ih++) {
			k = sindex[ih] ;
			mh+=sm[k]/N_STAR ;
			if ((mh/Mtotal[0]) > 0.5) {
				rh = sr[k] ;
				break ;
			}
		}

	}

	trh=((0.138 * N_MAX)/log(N_MAX)) * sqrt((rh*rh*rh)/Mtotal[0])*log(N_MAX)/N_MAX;



/* Concentration parameter --- note that max_r is the max radius of all bound
   stars, returned by get_positions(). When a finite R_MAX (i.e. initial Rtidal)
   is specified, max_r should be approximately equal to Rtidal. But when 
   Rtidal is very large (isolated cluster), r_max still gives the "actual"
   size of the cluster, indicating the maximum radius.

   The conc parameter is actually defined as LOG_10(rtidal/rcore). But for
   easier reading, it is calculated here only as rtidal/rcore.

*/
	conc_param = (max_r/core_radius) ;  /* max_r & core_radius are output in out_3 file */



	if(printstatus==0 || printstatus==6) {
		fprintf(stderr,"\n>>>sub_count = %ld  sub_FACTOR = %ld  sub_N_MAX = %ld sub_rmax = %.6g \n", \
			sub_count, sub_FACTOR, sub_N_MAX, sub_rmax) ;
		fprintf(stderr,"\n>>>Tcount = %ld  Time = %.6G  Dt = %.6G  Trc_factor = %.6G \n", \
			tcount,TotalTime, Dt, TRC_FACTOR) ;
		fprintf(stderr,"Etotal = %.6G Rmax = %.6G N_MAX = %5ld  Rtidal = %.6G \n", \
			Etotal, max_r, errstat, Rtidal) ;
		fprintf(stderr,"Total Mass = %.6G PEtotal = %.6G KEtotal = %.4G VRatio = %.6G \n", \
			Mtotal[0], PEtotal, KEtotal, -2.0*KEtotal/PEtotal) ;
		fprintf(stderr,"Real Time = %.6G Stellar Mass Loss = %.6G Tidal Mass Loss = %.6G \n", \
			RealTime, StellarMassLoss, TidalMassLoss) ;
		fprintf(stderr,
			"Core (%ld stars): core radius = %.6G  density = %.6G  <v> = %.6G  t_rc = %.6G conc = %.6G  N_core = %.6G\n", \
			NUM_CORE_STARS, core_radius, rho_core, v_core, Trc, conc_param, N_core);
		fprintf(stderr,\
			"DEavg = %.8g  DEmax = %.8g  DEdiff = %.8g  Sin2Beta = %.6G  <sin2beta> = %.8g  max(sin2beta) = %.8g\n", \
			DEavg, DEmax, DEdiff, Sin2Beta, avg_sin2beta, max_sin2beta) ;
		fprintf(stderr, \
			"n_orbit = %ld  n_extreme = %ld n_pericenter = %ld  n_apocenter = %ld\n", \
			n_orbit, n_extreme, n_pericenter, n_apocenter) ;
		fprintf(stderr, \
			"Half-mass: index = %ld  M_h = %.6G  rh = %.6G  t_rh = %.6G\n",ih, mh, rh, trh);
		fprintf(stderr, "\n") ;

/***** Repeat log output for out_0 file ********/

		fprintf(logfile,"\n>>>sub_count = %ld  sub_FACTOR = %ld  sub_N_MAX = %ld sub_rmax = %.6g \n", \
			sub_count, sub_FACTOR, sub_N_MAX, sub_rmax) ;
		fprintf(logfile,"\n>>>Tcount = %ld  Time = %.6G  Dt = %.6G  Trc_factor = %.6G \n", \
			tcount,TotalTime, Dt, TRC_FACTOR) ;
		fprintf(logfile,"Etotal = %.6G Rmax = %.6G N_MAX = %5ld  Rtidal = %.6G \n", \
			Etotal, max_r, errstat, Rtidal) ;
		fprintf(logfile,"Total Mass = %.6G PEtotal = %.6G KEtotal = %.4G VRatio = %.6G \n", \
			Mtotal[0], PEtotal, KEtotal, -2.0*KEtotal/PEtotal) ;
		fprintf(logfile,"Real Time = %.6G Stellar Mass Loss = %.6G Tidal Mass Loss = %.6G \n", \
			RealTime, StellarMassLoss, TidalMassLoss) ;
		fprintf(logfile,
			"Core (%ld stars): core radius = %.6G  density = %.6G  <v> = %.6G  t_rc = %.6G conc = %.6G  N_core = %.6G\n", \
			NUM_CORE_STARS, core_radius, rho_core, v_core, Trc, conc_param, N_core);
		fprintf(logfile,\
			"DEavg = %.8g  DEmax = %.8g  DEdiff = %.8g  Sin2Beta = %.6G  <sin2beta> = %.8g  max(sin2beta) = %.8g\n", \
			DEavg, DEmax, DEdiff, Sin2Beta, avg_sin2beta, max_sin2beta) ;
		fprintf(logfile, \
			"n_orbit = %ld  n_extreme = %ld n_pericenter = %ld  n_apocenter = %ld\n", \
			n_orbit, n_extreme, n_pericenter, n_apocenter) ;
		fprintf(logfile, \
			"Half-mass: index = %ld  M_h = %.6G  rh = %.6G  t_rh = %.6G\n",ih, mh, rh, trh);
		fprintf(logfile, "\n") ;

		if(printstatus==6) gets(&cont) ;

	}

}


void NormalizeEnergy(void) {

	long i, k ;
	double Dnorm ;


	Dnorm = sqrt(fabs((-0.25 - PEtotal) / KEtotal)) ;  /* normalization needed for velocity */
	for (i = 1 ; i <= N_MAX ; i++) {
		svr[sindex[i]] *= Dnorm ;
		svt[sindex[i]] *= Dnorm ;
	}

/* recompute total E -- should be -0.25 EXACTLY 
*/
	Etotal = 0.0 ;
	KEtotal = 0.0 ; PEtotal = 0.0 ;

	for (i = 1 ; i <= N_MAX ; i++) {
		k = sindex[i] ;
		sE[k] = sgravity[k] 
			+ 1.0/2.0*(svr[k]*svr[k] + \
			svt[k]*svt[k]) ;

		sJ[k] = sr[k] * svt[k] ;

		KEtotal = KEtotal + 1.0/2.0*(svr[k]*svr[k] + \
			svt[k]*svt[k]) * sm[k] / N_STAR;

		PEtotal = PEtotal + sgravity[k] * sm[k] / N_STAR;
		Etotal_New = Etotal_New + sE[k] * sm[k] / N_STAR;
	}
	/* PEtotal = 1/2 * Integral(phi(r)*rho(r) d3r) */
	PEtotal = PEtotal / 2.0 ;	
	Etotal = KEtotal + PEtotal ;

	if(printstatus==0 || printstatus==6) {	
		fprintf(stderr,"Energy Normalized to -0.25\n") ;
		fprintf(stderr,"Time = %.8G   Tcount = %ld\n", TotalTime, tcount) ;
		fprintf(stderr,"N = %ld, Total E = %.8G, Total Mass = %.8G, Virial ratio = %.8G\n", \
				N_MAX, Etotal, Mtotal[0], -2.0*KEtotal/PEtotal);
		fprintf(stderr,"Total KE = %.8G, Total PE = %.8G\n", KEtotal, PEtotal);
	}
}

void ComputeEnergy(void) {

	long k, i ;

	Etotal = 0.0 ;
	KEtotal = 0.0 ; PEtotal = 0.0 ;

	for (i = 1 ; i <= N_MAX ; i++) {
		k = sindex[i] ;
		sE[k] = sgravity[k] + 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) ;

		sJ[k] = sr[k] * svt[k] ;

		KEtotal += 1.0/2.0*(svr[k]*svr[k] + svt[k]*svt[k]) * sm[k] / N_STAR;

		PEtotal += sgravity[k] * sm[k] / N_STAR;
	}
	/* PEtotal = 1/2 * Integral(phi(r)*rho(r) d3r) */
	PEtotal = PEtotal / 2.0 ;	
	Etotal = KEtotal + PEtotal ;

	if(printstatus==0 || printstatus==6) {	
		fprintf(stderr,"Time = %.8G   Tcount = %ld\n", TotalTime, tcount) ;
		fprintf(stderr,"N = %ld, Total E = %.8G, Total Mass = %.8G, Virial ratio = %.8G\n", \
				N_MAX, Etotal, Mtotal[0], -2.0*KEtotal/PEtotal);
		fprintf(stderr,"Total KE = %.8G, Total PE = %.8G\n", KEtotal, PEtotal);
	}

}



/* binary - single interaction cross section */

double bin_single_sigma(double y) {

return (1.0/(pow(1.0+y, 4.0) * sqrt(y)) ) ;

}


/*	Peturbing the velocities of stars. Interacting star 1 with 3,
	2 with 4, and so on. Energy is exactly conserved in this step.
	If N_MAX is not divisible by 4, some stars are not perturbed 
	at the end.
*/


void perturb_stars(double Dt) {

static double v[4], vp[4], w[4], W, phi, psi, l, beta, wp, w1[4], w2[4] ;
static double v_new[4], vp_new[4], vr_new, vt_new, vrp_new, vtp_new ;
static double r, rp, rm, Dr, m, mp, w_new[4] ;
static double DeltaE, DeltaEp, vr, vt, vrp, vtp, Eold, Epold, Ai, SaveDt ;
static long p = 20, si, i, sip, k, kp, j, si_minus_p, si_plus_p, RepeatEncounter ;
static double min_dist, p_max, P_enc, n_local, n_bin_local, del_v, a, b, c, theta ;
static double BE, DBE, DBE1, DBE2, DBE3, rpp, vrpp, vtpp, mpp, vrpp_new, vtpp_new ;
static double Eppold, DeltaEpp, DY, BE1, BE2 ;
static long interactiondone, bi, bcount, i_nearest, kpp, N_LIMIT ;

max_sin2beta = 0.0 ;
avg_sin2beta = 0.0 ;
DEavg = 0.0 ;	/* Average fractional energy change <deltaE/E> */
DEmax = 0.0 ;	/* MAX fractional energy change Max{deltaE/E} */
DEdiff = 0.0 ;	/* Total Energy change -- should remain ZERO */
SaveDt = Dt ;	/* Original Dt provided, saved for repeated encounters, where Dt is changed */
Delta_BE_bb = 0.0 ; /* Total energy generated by binary-binary interactions */
Delta_BE_bs = 0.0 ; /* Total energy generated by binary-single interactions */
N_MAX_NEW = N_MAX ; /* keeps track of stars from disrupted binaries */

/* First loop over ALL binaries, and check for binary-binary interactions */

if (sub_count == 0) {
	N_LIMIT = N_MAX ;		/* do FULL time step */
}
else {	
	N_LIMIT = sub_N_MAX ;	/* do sub step only */
}

for (si = 1 ; si <= N_LIMIT ; si++) {

	if (si > sub_N_MAX)	{	/* part of outer halo in a FULL time step */
		Dt = sub_totaltime ;
	}
	else {
		Dt = SaveDt ;
	}

	current_star_i = si ;
	k = sindex[si] ;  
	if (s_binindex[k] > 0 && sinteracted[k]==0) {	
		/* Uninteracted binary -- so check for bin-bin interaction */

		/* Compute local number density of binaries -- find nearest 5 binaries
		   (or less), by looking out to the nearest 100 stars. 
		   Also, find nearest UNINTERACTED binary to interact with.
		*/
		bcount = 0 ;
		si_minus_p = si_plus_p = si ;
		i_nearest = si ;
		for (i=1 ; i<=100 ; i++){
			if (si+i <= N_LIMIT && s_binindex[sindex[si+i]] > 0) {
				bcount++ ;
				si_plus_p = si+i ;
				if (i_nearest == si && sinteracted[sindex[si+i]]==0) {
						i_nearest = si+i ;
				}
			}
			if (si-i > 0 && s_binindex[sindex[si-i]] > 0) {
				bcount++ ;
				si_minus_p = si-i ;
				if (i_nearest == si && sinteracted[sindex[si-i]]==0) {
						i_nearest = si-i ;
				}
			}
			if (bcount > 10) break ;	/* stop when 10 nearest binaries are found */
		}

		bcount++ ;		/* To account for the binary at si */	
		kp = sindex[i_nearest] ;

		if (bcount == 1) {		/* No other binary found nearby - so set density = 0 */
			n_bin_local = 0.0 ;
			P_enc = 0.0 ;		/* set collision probability = 0 */
		}
		if (bcount == 2) {		/* only one other nearby binary found */
			n_bin_local =  \
				(bcount-1)*3.0/4.0/pi/(pow(sr[sindex[si_plus_p]],3.0) \
				-pow(sr[sindex[si_minus_p]],3.0)) ;		
		}
		else {  	/* in computing the density, DON'T count the boundary stars */		
			n_bin_local =  \
				(bcount-2)*3.0/4.0/pi/(pow(sr[sindex[si_plus_p]],3.0) \
				-pow(sr[sindex[si_minus_p]],3.0)) ;
		}
		

		/* check for bin-bin interaction between star si & i_nearest */
		
		if (i_nearest != si && n_bin_local > 0) {

			r = sr[k] ; rp = sr[kp] ; rm = (r+rp)/2.0 ;
			vr = svr[k] ; vt = svt[k] ; vrp = svr[kp] ; vtp = svt[kp] ;	
			Dr = fabs(rp - r) ;
			m = sm[k] ; mp = sm[kp] ;

			/* Compute relative speed W */
			phi = ran2(&idum) * 2*pi ;
			v[1] = vt ; v[2] = 0.0 ; v[3] = vr ;
			vp[1] = vtp*cos(phi) ; vp[2] = vtp*sin(phi) ; vp[3] = vrp ;
			
			for (j=1;j<=3;j++) w[j] = vp[j] - v[j] ;
			W = sqrt(w[1]*w[1] + w[2]*w[2] + w[3]*w[3]) ;

			if (W == 0.0) {
				fprintf(logfile, "ERROR: Star %ld : W = 0 in perturb_stars() - 1\n",si) ;
				fprintf(logfile, \
				"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
				si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
				fprintf(stderr, "ERROR: Star %ld : W = 0 in perturb_stars() - 1\n",si) ;
				fprintf(stderr, \
				"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
				si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
			}

			/* set distance of closest approach for collision based on softer binary */
			
			if (ba[s_binindex[k]] > ba[s_binindex[kp]]) 
				min_dist = 1.7 * (ba[s_binindex[k]]) ;
			 else
				min_dist = 1.7 * (ba[s_binindex[kp]]) ;
			

			/* min_dist = ba[s_binindex[k]] + ba[s_binindex[kp]] ;
			*/

			/* spitzer eq. 6-15 for max impact parameter for min_dist 
			p_max = sqrt(min_dist*min_dist*(1.0 + 2.0*(m+mp)/N_STAR/min_dist/(W*W))) ;
			*/

			p_max = sqrt(6.0*m/N_STAR*min_dist/(W*W)) ;

			/* Probability of bin-bin encounter P = pi*p_max^2*W*n*Dt 
			P_enc = pi * p_max*p_max * W * Dt*N_STAR/log(N_STAR) * n_bin_local/1.0 ;
			*/
			
			/* compute N_core = 2/3*pi*rcore^3/<m>
			   = number of binaries or single stars in the core (since m_bin = 2 <m>)
			*/
			N_core = 2.0/3.0*pi*core_radius*core_radius*core_radius*(1.0*N_STAR) * rho_core ;

			P_enc = pi * p_max*p_max * W * Dt * (1.0*N_STAR)/log(0.4*N_core) * n_bin_local/1.0 ;

		}

		if (ran2(&idum) < P_enc) {		/* check if bin-bin interaction is DUE */

			/* DO interaction */				
			fprintf(stderr, "\n\nbinary - binary interaction DUE \n\n") ;
			fprintf(logfile, "\n\nbinary - binary interaction DUE \n\n") ;
			N_bb++ ;

			/* Compute energy generated from binary -- Delta BE  */
			
			BE1 = bm1[s_binindex[k]]*bm2[s_binindex[k]]/N_STAR/N_STAR/ \
				2.0/ba[s_binindex[k]] ;
			BE2 = bm1[s_binindex[kp]]*bm2[s_binindex[kp]]/N_STAR/N_STAR/ \
			    2.0/ba[s_binindex[kp]] ;
			BE = BE1+BE2 ;

			for (j=1; j<=N_TRY ; j++) {
				DY = ran2(&idum) * 1.0 ;
				if (ran2(&idum)*1.9 < 12.25*DY/pow((1.0+3.5*DY*DY),2.75) )break ;
			}
			if (j==N_TRY+1) {
				fprintf(stderr, "\n\nERROR - N_TRY exceeded in bin-bin cross section\n\n") ;
				fprintf(logfile,"\n\nERROR - N_TRY exceeded in bin-bin cross section\n\n") ;
			}
			DBE = DY * BE ;
			DBE1 = DBE / 4.0 ;		/* recoil energy of binary */
			DBE2 = (DBE - DBE1)*ran2(&idum) ;	/* recoil energy of first star */
			DBE3 = DBE - DBE1 - DBE2 ;			/* recoil energy of second star */

			/* check which binary is harder */
			if (BE1 < BE2) {			/* binary kp is harder than k */
				bi=k ; k=kp ; kp=bi ;	/* so interchange binary k & kp */
				BE = BE1; BE1 = BE2 ; BE2 = BE ; BE = BE1+BE2 ;
				r = sr[k] ; rp = sr[kp] ; rm = (r+rp)/2.0 ;
				vr = svr[k] ; vt = svt[k] ; vrp = svr[kp] ; vtp = svt[kp] ;	
				Dr = fabs(rp - r) ;
				m = sm[k] ; mp = sm[kp] ;
			}
			
			/* Set new (smaller) semi-major axis for binary k (harder binary) 
			   New BE = original(BE1+BE2)+DBE = 2*DBE + DBE = 3*DBE 
			*/
			ba[s_binindex[k]] = bm1[s_binindex[k]]*bm2[s_binindex[k]]/N_STAR/N_STAR/2.0/(BE1+BE2+DBE) ;

			/* binding energy of soft binary is lost since it is destroyed.
			   But it is compensated by harder binary getting harder. So total BE
			   does not decrease due to the destruction of the binary. In addition,
			   more translational energy 0.5*(BE1+BE2) is released by the harder binary.
			*/
			E_binding += DBE ;		/* total BE in binaries */
			Delta_BE_bb += DBE ;	/* total bin-bin energy generated in the timestep */
			DE_bb += DBE ;			/* ongoing total bin-bin energy generated */

		/* An additional second star is also given energy */
			N_STAR_NEW++ ; N_MAX_NEW++ ;
			bi = N_MAX_NEW ;				/* new star */
			kpp = sindex[bi] = N_STAR_NEW ;	/* new star is placed at the end for now */
			stype[kpp] = 1 ;	/* asume new star is single MS star */
			s_binindex[kpp] = 0 ;
			sm[kpp]=bm2[s_binindex[kp]] ;

			/* use same positions and velocities as star kp */
			/* NOTE: DO NOT ASSIGN position sr[kpp] yet, since routines
			   potential() and get_positions() require that sr[N_MAX+1] = INFINITY.
			   Changing this will cause problems in potential().
			   Simply remember that position of kpp is same as kp = sindex[old_k[kpp]].
			   We also go ahead and assign sr[sindex[N_MAX+1]] = INFINITY
			*/
			sr[kpp] = INFINITY ;  /* needed for normal operation until N_MAX is updated */
			sgravity[kpp] = 0.0 ;
			sMi_ri[kpp] = 0.0 ;
			sMi[kpp] = Mtotal[0] ;
			svr[kpp] = svr[kp] ; svt[kpp] = svt[kp] ; /* OK to assign velocities though */

			rpp = sr[kp] ; vrpp = svr[kpp] ; vtpp = svt[kpp] ; mpp = sm[kpp] ;			

			/* save indexed position for new star for use in get_positions() */
			if (kp == sindex[i_nearest])	/* kp is sindex[i_nearest] */
				old_k[kpp] = i_nearest ;
			else							/* kp is sindex[si] */
				old_k[kpp] = si ;


			/* binary kp is destroyed */
			stype[kp] = 1 ; 
			count_5-- ; mass_5 -= sm[kp] ;  			
			count_1++ ; mass_1 += sm[kp] ;
			mp = sm[kp] = bm1[s_binindex[kp]] ;	/* keep first star of the binary as single */
			bini_tdestroyed[s_binindex[kp]] = TotalTime ;
			/* label binary kp as destroyed */
			s_binindex[kp] = 0 ;

/*		for (j=si+1 ; j<= N_MAX ; j++) {
				bi = sindex[j] ;
				if (stype[sindex[j]] != 5 && sinteracted[sindex[j]] == 0 \
					&& bi != k && bi != kp) break ;
			}
			kpp = bi ;
*/



			
			/*********** compute recoil velocities - binary-binary int **********/

			theta = acos(ran2(&idum)*2.0 - 1.0) ;
			phi = 2*pi*ran2(&idum) ;

			/* recoil of the binary */
			a = 1 ; b = 2*(vr*cos(theta)+vt*sin(theta)*cos(phi)) ; c = -2*DBE1/m*N_STAR ;
			del_v = (-b + sqrt(b*b - 4.0*a*c))/2.0/a ;

			vr_new = vr + del_v*cos(theta) ; 
			vt_new = sqrt( pow(vt + del_v*sin(theta)*cos(phi), 2.0) + \
						   pow (del_v*sin(theta)*sin(phi), 2.0) ) ;

			/* recoil of the star */
			a = 1 ; b = 2*(vrp*cos(theta)+vtp*sin(theta)*cos(phi)) ; c = -2*DBE2/mp*N_STAR ;
			del_v = (-b + sqrt(b*b - 4.0*a*c))/2.0/a ;

			vrp_new = vrp + del_v*cos(theta) ; 
			vtp_new = sqrt( pow(vtp + del_v*sin(theta)*cos(phi), 2.0) + \
						   pow(del_v*sin(theta)*sin(phi), 2.0) ) ;

			/* recoil of the second star*/
			a = 1 ; b = 2*(vrpp*cos(theta)+vtpp*sin(theta)*cos(phi)) ; c = -2*DBE3/mpp*N_STAR ;
			del_v = (-b + sqrt(b*b - 4.0*a*c))/2.0/a ;

			vrpp_new = vrpp + del_v*cos(theta) ; 
			vtpp_new = sqrt( pow(vtpp + del_v*sin(theta)*cos(phi), 2.0) + \
						   pow(del_v*sin(theta)*sin(phi), 2.0) ) ;


	/* recompute energies of binary and star */

			DeltaE =  0.5 * (vr_new*vr_new + vt_new*vt_new - vr*vr - vt*vt) ;
			DeltaEp = 0.5 * (vrp_new*vrp_new + vtp_new*vtp_new - vrp*vrp - vtp*vtp) ;
			DeltaEpp = 0.5 * (vrpp_new*vrpp_new + vtpp_new*vtpp_new - vrpp*vrpp - vtpp*vtpp) ;

			Eold = sE[k] ; Epold = sE[kp] ; Eppold = sE[kpp] ;

			/* MAX fractional change in energy over ALL stars */
			if (fabs(DeltaE/Eold) > DEmax) DEmax = fabs(DeltaE/Eold) ;
			if (fabs(DeltaEp/Epold) > DEmax) DEmax = fabs(DeltaEp/Epold) ;
			

	/*	Calculate new energies by recomputing E = PE + KE using new velocity
	*/		sE[k] = sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
			sJ[k] = r * vt_new ;
			sE[kp] = sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
			sJ[kp] = rp * vtp_new ;
			/* sgravity[] is not yet computed for kpp, so use position of kp */
			sE[kpp] = sgravity[kp] + 0.5*(vrpp_new*vrpp_new + vtpp_new*vtpp_new) ;
			sJ[kpp] = rpp * vtpp_new ;
	/*		sE[k] = -sMi[k]/sr[k] - sMi_ri[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
			sE[kp] = -sMi[kp]/sr[kp] - sMi_ri[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
	*/		

	/* set new velocities for both stars 
	*/
			sDE[k] = DeltaE/Eold ;
			sDE[kp] = DeltaEp/Epold ;	
			sDJ[k] = (sJ[k] - r*vt) / (r*vt) ;
			sDJ[kp] = (sJ[kp] - rp*vtp) / (rp*vtp) ;
			sDE[kpp] = DeltaEpp/Eppold ;
			sDJ[kpp] = (sJ[kpp] - rpp*vtpp) / (rpp*vtpp) ;
			sSin2Beta[k] = 0 ;
			sSin2Beta[kp] = 0 ;
			sSin2Beta[kpp] = 0 ;

			svr[k] = vr_new ; 
			svt[k] = vt_new ;
			svr[kp] = vrp_new; 
			svt[kp] = vtp_new ;
			svr[kpp] = vrpp_new; 
			svt[kpp] = vtpp_new ;

			DEdiff += (Eold - (sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new))) ;
			DEdiff += (Epold - (sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new))) ;
			DEdiff += (Eppold - (sgravity[kpp] + 0.5*(vrpp_new*vrpp_new + vtpp_new*vtpp_new))) ;

		/* binary-binary interaction DONE */
			sinteracted[k] = sinteracted[kp] = sinteracted[kpp] = 3 ;

		}
		else {
		/* NO binary-binary interaction due... so label binaries uninteracted.
		   It is possible for k==kp if i_nearest = si, i.e., no nearby bin was found.
		*/
			sinteracted[k] = sinteracted[kp] = 0 ;		
		}


	}	/* End IF (uninteracted binary) */

}		/* End For loop (over all stars to look for uninteracted binaries) */




/**** Compute binary-single and single-single (2-body relaxation) interactions ****/

for (si = 1 ; si <= N_LIMIT ; si += 2) {

	if (si > sub_N_MAX)	{	/* part of outer halo in a FULL time step */
		Dt = sub_totaltime ;
	}
	else {
		Dt = SaveDt ;
	}

	while (si < N_LIMIT && sinteracted[sindex[si]] > 0) {
		si++ ;
	}
	bi=si+1 ;
	if (bi > N_LIMIT) break ;	/* ignore remaining stars if not paired */
	while (sinteracted[sindex[bi]] > 0) {
		bi++ ;	
		if (bi > N_LIMIT) break ;	/* ignore remaining stars if not paired */
	}
	if (bi > N_LIMIT) break ;	/* ignore remaining stars if not paired */

	interactiondone = 0 ;
	k = sindex[si] ; kp = sindex[bi] ; 
	
	/* ensure that k & kp represent binary/single and single, respectively */
	if (s_binindex[k] > 0 || s_binindex[kp] > 0)  {
		if (s_binindex[k] > 0) {
			/* k and kp already represent bin & star/binary, respectively. */
		}
		else {
			sip = kp; kp = k ; k = sip ; bi = si+1 ;	/* interchange k & kp */
		}
	}

	/* if one of the objects is a binary, then check for binary-single interaction */
	if (s_binindex[k] > 0 || s_binindex[kp] > 0)  {	

		r = sr[k] ; rp = sr[kp] ; rm = (r+rp)/2.0 ;
		vr = svr[k] ; vt = svt[k] ; vrp = svr[kp] ; vtp = svt[kp] ;	
		Dr = fabs(rp - r) ;
		m = sm[k] ; mp = sm[kp] ;

		/* Compute relative speed W */
		phi = ran2(&idum) * 2*pi ;
		v[1] = vt ; v[2] = 0.0 ; v[3] = vr ;
		vp[1] = vtp*cos(phi) ; vp[2] = vtp*sin(phi) ; vp[3] = vrp ;
		
		for (j=1;j<=3;j++) w[j] = vp[j] - v[j] ;
		W = sqrt(w[1]*w[1] + w[2]*w[2] + w[3]*w[3]) ;

		if (W == 0.0) {
			fprintf(logfile, "ERROR: Star %ld : W = 0 in perturb_stars() - 2\n",si) ;
			fprintf(logfile, \
			"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
			si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
			fprintf(stderr, "ERROR: Star %ld : W = 0 in perturb_stars() - 2\n",si) ;
			fprintf(stderr, \
			"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
			si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
		}
		
		/* check for binary-single interaction */

		min_dist = 3.5 * ba[s_binindex[k]] ;

		/* spitzer eq. 6-15 for max impact parameter for min_dist */
		p_max = sqrt(min_dist*min_dist*(1.0 + 2.0*(m+mp)/N_STAR/min_dist/(W*W))) ;
		
		/* Computing local density */
		si_minus_p = si-p ; 
		si_plus_p = si+p+1 ;

		if (si_minus_p < 0) {
			si_minus_p = 0 ;
			si_plus_p = 2*p + 1 ;
		}
		else if (si_plus_p > N_LIMIT) {
			si_plus_p = N_LIMIT ;
			si_minus_p = N_LIMIT - 2*p - 1 ;
		}

		n_local =  \
		2.0*p*3.0/4.0/pi/(pow(sr[sindex[si_plus_p]],3.0)-pow(sr[sindex[si_minus_p]],3.0)) ;

		/* compute number of stars/binaries in the core 
		*/
		N_core = 2.0/3.0*pi*core_radius*core_radius*core_radius*(1.0*N_STAR) * rho_core ;

		/* Probability of binary-single collision P = pi*p_max^2*W*n*Dt 
		*/
		P_enc = pi * p_max*p_max * W * Dt * (1.0*N_STAR)/log(0.4*N_core) * n_local / 1.0;

		if (ran2(&idum) < P_enc) {	/* binary-single Collision DUE */

			/* DO collision */				
			
			fprintf(stderr, "\n\nbinary - single interaction DUE \n\n") ;
			fprintf(logfile, "\n\nbinary - single interaction DUE \n\n") ;
			N_bs++ ;

			/* Compute energy generated from binary -- Delta BE  */
			
			BE = bm1[s_binindex[k]]*bm2[s_binindex[k]]/N_STAR/N_STAR/2.0/ba[s_binindex[k]] ;
			
			for (j=1; j<=N_TRY ; j++) {
				DY = ran2(&idum)*1.999 + 0.001 ;
				if (ran2(&idum)*1.1*bin_single_sigma(0.001) < bin_single_sigma(DY)) break ;
			}
			if (j==N_TRY+1) {
				fprintf(stderr, "\n\nERROR - N_TRY exceeded in bin-single cross section\n\n") ;
				fprintf(logfile,"\n\nERROR - N_TRY exceeded in bin-single cross section\n\n") ;
			}

			DBE = DY * BE ;
/*			DBE = 0.2 * BE ; 
*/
			DBE1 = DBE / 3 ;		/* recoil energy of binary */
			DBE2 = DBE - DBE1 ;		/* recoil energy of star */

			/* Set new (smaller) semi-major axis for binary */
			ba[s_binindex[k]] = bm1[s_binindex[k]]*bm2[s_binindex[k]]/N_STAR/N_STAR/2.0/(BE+DBE) ;

			E_binding += DBE ;		/* total BE in binaries */
			Delta_BE_bs += DBE ;	/* total bin-single energy generated in the timestep */
			DE_bs += DBE ;			/* ongoing total bin-single energy generated */

			theta = acos(ran2(&idum)*2.0 - 1.0) ;
			phi = 2*pi*ran2(&idum) ;

			/* recoil of the binary */
			a = 1 ; b = 2*(vr*cos(theta)+vt*sin(theta)*cos(phi)) ; c = -2*DBE1/m*N_STAR ;
			del_v = (-b + sqrt(b*b - 4.0*a*c))/2.0/a ;

			vr_new = vr + del_v*cos(theta) ; 
			vt_new = sqrt( pow(vt + del_v*sin(theta)*cos(phi), 2.0) + \
						   pow (del_v*sin(theta)*sin(phi), 2.0) ) ;

			/* recoil of the star/binary */
			a = 1 ; b = 2*(vrp*cos(theta)+vtp*sin(theta)*cos(phi)) ; c = -2*DBE2/mp*N_STAR ;
			del_v = (-b + sqrt(b*b - 4.0*a*c))/2.0/a ;

			vrp_new = vrp + del_v*cos(theta) ; 
			vtp_new = sqrt( pow(vtp + del_v*sin(theta)*cos(phi), 2.0) + \
						   pow (del_v*sin(theta)*sin(phi), 2.0) ) ;


	/* recompute energies of binary and star */

			DeltaE =  0.5 * (vr_new*vr_new + vt_new*vt_new - vr*vr - vt*vt) ;
			DeltaEp = 0.5 * (vrp_new*vrp_new + vtp_new*vtp_new - vrp*vrp - vtp*vtp) ;

			Eold = sE[k] ; Epold = sE[kp] ;

			/* fractional change in energy -- total ALL stars -- divide by N_MAX in MAIN */
			DEavg += fabs(DeltaE/Eold) + fabs(DeltaEp/Epold) ;

			/* MAX fractional change in energy over ALL stars */
			if (fabs(DeltaE/Eold) > DEmax) DEmax = fabs(DeltaE/Eold) ;
			if (fabs(DeltaEp/Epold) > DEmax) DEmax = fabs(DeltaEp/Epold) ;
			

	/*	Calculate new energies by recomputing E = PE + KE using new velocity
	*/		sE[k] = sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
			sJ[k] = r * vt_new ;
			sE[kp] = sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
			sJ[kp] = rp * vtp_new ;
	/*		sE[k] = -sMi[k]/sr[k] - sMi_ri[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
			sE[kp] = -sMi[kp]/sr[kp] - sMi_ri[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
	*/		

	/* set new velocities for both stars 
	*/
			sDE[k] = DeltaE/Eold ;
			sDE[kp] = DeltaEp/Epold ;	
			sDJ[k] = (sJ[k] - r*vt) / (r*vt) ;
			sDJ[kp] = (sJ[kp] - rp*vtp) / (rp*vtp) ;
			sSin2Beta[k] = sin2beta ;
			sSin2Beta[kp] = sin2beta ;

			svr[k] = vr_new ; 
			svt[k] = vt_new ;
			svr[kp] = vrp_new; 
			svt[kp] = vtp_new ;


			DEdiff += (Eold - (sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new))) ;
			DEdiff += (Epold - (sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new))) ;

		/* binary-single interaction DONE */
			interactiondone = 1 ;
			sinteracted[k] = sinteracted[kp] = 2 ;

		}
		else {
			/* NO binary interaction due -- continue with relaxation step */
			interactiondone = 0 ;
			sinteracted[k] = sinteracted[kp] = 0 ;
		}
	}


/* If not binary-binary or binary-single interaction, do regular relaxation step.
*/

	if (interactiondone == 0) {
	  
DoEncounterAgain:

		r = sr[k] ; rp = sr[kp] ; rm = (r+rp)/2.0 ;
		vr = svr[k] ; vt = svt[k] ; vrp = svr[kp] ; vtp = svt[kp] ;	
		Dr = fabs(rp - r) ;
		m = sm[k] ; mp = sm[kp] ;

		phi = ran2(&idum) * 2*pi ;
		v[1] = vt ; v[2] = 0.0 ; v[3] = vr ;
		vp[1] = vtp*cos(phi) ; vp[2] = vtp*sin(phi) ; vp[3] = vrp ;
		
		for (j=1;j<=3;j++) w[j] = vp[j] - v[j] ;
		W = sqrt(w[1]*w[1] + w[2]*w[2] + w[3]*w[3]) ;

		if (W == 0.0) {
			fprintf(logfile, "ERROR: Star %ld : W = 0 in perturb_stars() - 3\n",si) ;
			fprintf(logfile, \
			"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
			si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
			fprintf(stderr, "ERROR: Star %ld : W = 0 in perturb_stars() - 3\n",si) ;
			fprintf(stderr, \
			"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
			si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
		}

/* henon calculation of beta -- usually using p = 2 
*/		l = sqrt((2.0*rm*rm * Dr * N_STAR) / ((p-1.0)* W * Dt))  ;
		beta = 2.0*atan( (m+mp)/(W*W*l) ) ;


/* Stodolkiewicz calculation of beta */

		si_minus_p = si-p ; 
		si_plus_p = si+p+1 ;

		if (si_minus_p < 0) {
			si_minus_p = 0 ;
			si_plus_p = 2*p + 1 ;
		}
		else if (si_plus_p > N_LIMIT) {
			si_plus_p = N_LIMIT ;
			si_minus_p = N_LIMIT - 2*p - 1 ;
		}

		Ai = 3.0 * p / (pow(sr[sindex[si_plus_p]], 3.0) - pow(sr[sindex[si_minus_p]], 3.0)) \
			 * pow((m+mp), 2.0) / pow(W, 3.0) ;
		
		sin2beta = Ai * Dt / N_STAR;

/*		if (si == 9991)
			fprintf(stderr, \
			"si = %5d sin2beta = %.6G Ai = %.6G  svr = %.6G svt = %.6G sr = %.6G W = %.6G\n", \
			si, sin2beta, Ai, svr[k], svt[k], sr[k], W) ;
*/		

		/* sometimes, sin2beta becomes > 1, so need to handle that differently */

		if (sin2beta <=1) {		/* normal case -- only one encounter required */
			beta = 2.0 * asin(sqrt(sin2beta)) ;
			if (beta > pi/2) beta = pi - beta ;
			RepeatEncounter = 0 ;
			Dt = SaveDt ;		/* reset Dt to original value */
		}
		else {					/* abnormal case -- compute multiple encounters */
			beta = pi ;	
			sin2beta = 1.0 ;	/* use beta = pi for the first encounter */
			Dt = Dt / Ai ;		/* reduce Dt by Ai for the next encounter */
			RepeatEncounter = 1 ;	
		}
/*		sin2beta = pow(sin(beta/2.0), 2.0) ;
*/
		if (RepeatEncounter == 0) {
			avg_sin2beta += sin2beta ;
			if (sin2beta > max_sin2beta) max_sin2beta = sin2beta ;
		}

		wp = sqrt(w[1]*w[1] + w[2]*w[2]) ;
		if (wp == 0.0) {
			fprintf(stderr, "ERROR: Star %ld : wp = 0 in perturb_stars()\n",si) ;
			fprintf(stderr, \
			"ERROR: si = %5ld  r = %.6g  vr = %.6G  vt = %.6g  rp = %.6G  vrp = %.6g  vtp = %.6g\n", \
			si, sr[k], svr[k], svt[k], sr[kp], svr[kp], svt[kp]) ;
		}

		w1[1] = w[2]*W/wp ; w1[2] = -w[1]*W/wp ; w1[3] = 0.0 ;
		w2[1] = -w[1]*w[3]/wp ; w2[2] = -w[2]*w[3]/wp ; w2[3] = wp ;

		psi = ran2(&idum) * 2*pi ;
		for (j=1;j<=3;j++) {
			w_new[j] = w[j]*cos(beta)+w1[j]*sin(beta)*cos(psi)+w2[j]*sin(beta)*sin(psi) ;
		}

		for (j=1;j<=3;j++) {
			v_new[j] = v[j] - mp/(m+mp)*(w_new[j] - w[j]) ;
			vp_new[j] = vp[j] + m/(m+mp)*(w_new[j] - w[j]) ;
		}

/* New velocities found 		
*/		vr_new = v_new[3] ; 
		vt_new = sqrt(v_new[1]*v_new[1] + v_new[2]*v_new[2]) ;
		vrp_new = vp_new[3] ; 
		vtp_new = sqrt(vp_new[1]*vp_new[1] + vp_new[2]*vp_new[2]) ;

/*  Use random directions for new velocities 
		psi = acos(ran2(&idum)*2.0 - 1.0) ;	
		V = sqrt(v_new[1]*v_new[1] + v_new[2]*v_new[2] + v_new[3]*v_new[3]) ;
		vr_new = V * cos(psi) ;
		vt_new = V * sin(psi) ;
		psi = acos(ran2(&idum)*2.0 - 1.0) ;	
		V = sqrt(vp_new[1]*vp_new[1] + vp_new[2]*vp_new[2] + vp_new[3]*vp_new[3]) ;
		vrp_new = V * cos(psi) ;
		vtp_new = V * sin(psi) ;
*/

		/* DEBUG */
		if (0) {
		/* when sin2beta > 1, encounter is repeated until sin2beta <= 1 */
		if (RepeatEncounter > 0) {	
			vr = vr_new ;	/* assign new velocities for use in repeated encounter */
			vt = vt_new ;
			vrp = vrp_new ;
			vtp = vtp_new ;
			goto DoEncounterAgain ;	
		}


/******** Check to see if a STRONG interaction is due for these stars   ******
*/		
		/* save value of beta due to 2-body relaxation */
/*		beta_relax = beta ;	
		local_density = 2*p / (4*pi/3) / (pow(sr[sindex[si_plus_p]], 3.0) - pow(sr[sindex[si_minus_p]], 3.0))

		lambda = local_density * W * Dt ;	// number of interactions per unit area

		for (j=1 ; j<=N_TRY ; j++) {
			
			prob_b = 2*pi*lambda* b_close * exp(-lambda*pi*b_close*b_close) ;
		}
		sin2beta = Ai * Dt / N_STAR;
*/
/*****************************************************************************/



		DeltaE =  0.5 * (vr_new*vr_new + vt_new*vt_new - vr*vr - vt*vt) ;
		DeltaEp = 0.5 * (vrp_new*vrp_new + vtp_new*vtp_new - vrp*vrp - vtp*vtp) ;

		Eold = sE[k] ; Epold = sE[kp] ;

		/* fractional change in energy -- total ALL stars -- divide by N_MAX in MAIN */
		DEavg += fabs(DeltaE/Eold) + fabs(DeltaEp/Epold) ;

		/* MAX fractional change in energy over ALL stars */
		if (fabs(DeltaE/Eold) > DEmax) DEmax = fabs(DeltaE/Eold) ;
		if (fabs(DeltaEp/Epold) > DEmax) DEmax = fabs(DeltaEp/Epold) ;
		
/*	Calculate new energies by just adding the change in KE (due to new velocity) 
		sE[k] = sE[k] + DeltaE ;
		sJ[k] = r * vt_new ;
		sE[kp] = sE[kp] + DeltaEp ;
		sJ[kp] = rp * vtp_new ;
*/
/*	Calculate new energies by recomputing E = PE + KE using new velocity
*/		sE[k] = sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
		sJ[k] = r * vt_new ;
		sE[kp] = sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
		sJ[kp] = rp * vtp_new ;
/*		sE[k] = -sMi[k]/sr[k] - sMi_ri[k] + 0.5*(vr_new*vr_new + vt_new*vt_new) ;
		sE[kp] = -sMi[kp]/sr[kp] - sMi_ri[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new) ;
*/		

/* set new velocities for both stars 
*/
		sDE[k] = DeltaE/Eold ;
		sDE[kp] = DeltaEp/Epold ;	
		sDJ[k] = (sJ[k] - r*vt) / (r*vt) ;
		sDJ[kp] = (sJ[kp] - rp*vtp) / (rp*vtp) ;
		sSin2Beta[k] = sin2beta ;
		sSin2Beta[kp] = sin2beta ;

		svr[k] = vr_new ; 
		svt[k] = vt_new ;
		svr[kp] = vrp_new; 
		svt[kp] = vtp_new ;


		DEdiff += (Eold - (sgravity[k] + 0.5*(vr_new*vr_new + vt_new*vt_new))) ;
		DEdiff += (Epold - (sgravity[kp] + 0.5*(vrp_new*vrp_new + vtp_new*vtp_new))) ;
		}
		/* DEBUG */
		sinteracted[k] = sinteracted[kp] = 1 ;

/*		check to make sure DE + DEp = 0 ... conservation of energy.
		if (DeltaE + DeltaEp > 1.0e-15 || DEdiff > 1.0e-12) {
            if(printstatus==0 || X==0) fprintf(stderr,"ERROR: si = %8ld ; DE + DEp = %.8G \n", \
				si, DeltaE + DeltaEp) ;
			if(printstatus==6) gets(&cont) ;
		}
*/

	}   /** End IF interaction NOT done (relaxation step) **/


	
}	/* End For-- next star */


	avg_sin2beta /= N_LIMIT ;		/* get AVERAGE of sin^2(beta/2) */
	DEavg /= N_LIMIT  ;			/* get AVERAGE fractional change in energy <deltaE/E> */

#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr,"Velocity perturbation complete for N_LIMIT = %ld\n", N_LIMIT) ;
#endif

}



/****** Rhomberg Integration for Improper Integrals ********/


#define FUNC(x) ((*func)(x))

double midpnt(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC


#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

double qromo(double (*func)(double), double a, double b,
	double (*choose)(double(*)(double), double, double, int))
{
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K


#define NRANSI

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
#undef NRANSI





/********* Find root using bisection method ***********/

#define JMAX 100
double rtbis(double (*func)(double), double x1, double x2, double xacc)
{
	void nrerror(char error_text[]);
	int j;

	double dx,f,fmid,xmid,rtb;
	f=(*func)(x1);
	fmid=(*func)(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}

#undef JMAX


/* function f(E) = M - E + e * sin(E) to be solved using rtbis()
   to get eccentric anomaly E, for a Keplerian orbit of eccentricity e
*/

double kepler_anomaly(double E)  {

	return (mean_anomaly - E + ORBIT_ECCENTRICITY * sin(E)) ;
}


/* used by FindZero() to get the largest k such that sr[sindex[k]] < r */

double function_r(long k, double r, double dummy) {

	return (sr[sindex[k]] - r) ;

}

double function_Q(long k, double E, double J) {

return (2.0*E - 2.0*sgravity[sindex[k]] - J*J/(sr[sindex[k]]*sr[sindex[k]])) ;

}

double function_minusAbsQ(long k, double E, double J) {

return(-fabs(2.0*E - 2.0*sgravity[sindex[k]] - J*J/(sr[sindex[k]]*sr[sindex[k]]))) ;

}



/***************************** Find Zero ***************************
	
	Finding the root of the function defined in 'func' by bisection. 
	Requires the root to be bracketed between x1 and x2.
*/
  
#define JMAX 500

long FindZero(double (*func)(long, double, double), long x1, long x2,
			  double E, double J)
{
	void nrerror(char error_text[]);
	int j;

	double f,fmid ;
	long xmid,rtb,dx ;

	f=(*func)(x1,E,J);
	fmid=(*func)(x2,E,J);
	if (f*fmid >= 0.0) 
		return (-1) ;
/*		nrerror("Root must be bracketed for bisection in rtbis");
*/
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+ (dx=(dx*0.5)),E,J);
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) == 0 || fmid == 0.0) break ;
	}
	if (j == JMAX + 1){
		nrerror("Too many bisections in rtbis");
		return 0 ;
	}
	f=(*func)(x1,E,J);
	if (f < 0.0) 
		dx = 1 ;
	else
		dx = -1;

	if ((*func)(rtb,E,J) > 0.0) 
		fprintf(stderr, "ERROR: f > 0 in FindZero\n") ;

	while ((*func)(rtb,E,J) < 0.0) {	/* look for change of sign */
			rtb = rtb + dx ;
	}
	if (dx == 1) rtb-- ;	/* if going from left to right, get lower point */
	return (rtb) ;			/* otherwise return upper point */
}

#undef JMAX

 
/***************************** Find Zero OF Q ***************************
	
	Finding the root of the function defined in 'func' by bisection. 
	Requires the root to be bracketed between x1 and x2.
*/
  
#define JMAX 500

long FindZero_Q(long x1, long x2, double E, double J)
{
	void nrerror(char error_text[]);
	int j;

	double f,fmid ;
	long xmid,rtb,dx ;

	f=(2.0*E - 2.0*sgravity[sindex[x1]] - J*J/(sr[sindex[x1]]*sr[sindex[x1]]));
	fmid=(2.0*E - 2.0*sgravity[sindex[x2]] - J*J/(sr[sindex[x2]]*sr[sindex[x2]]));
	if (f*fmid >= 0.0) 
		return (-1) ;
/*		nrerror("Root must be bracketed for bisection in rtbis");
*/
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		xmid=rtb+(dx = (dx*0.5)) ;
		fmid=(2.0*E - 2.0*sgravity[sindex[xmid]] - J*J/(sr[sindex[xmid]]*sr[sindex[xmid]]));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) == 0 || fmid == 0.0) break ;
	}
	if (j == JMAX + 1){
		nrerror("Too many bisections in rtbis");
		return 0 ;
	}
	f=(2.0*E - 2.0*sgravity[sindex[x1]] - J*J/(sr[sindex[x1]]*sr[sindex[x1]]));
	if (f < 0.0) 
		dx = 1 ;
	else
		dx = -1;

	if ((2.0*E - 2.0*sgravity[sindex[rtb]] - J*J/(sr[sindex[rtb]]*sr[sindex[rtb]])) > 0.0) 
		fprintf(stderr, "ERROR: f > 0 in FindZero\n") ;

	/* look for change of sign */
	while ((2.0*E - 2.0*sgravity[sindex[rtb]] - J*J/(sr[sindex[rtb]]*sr[sindex[rtb]])) < 0.0) {	
			rtb = rtb + dx ;
	}
	if (dx == 1) rtb-- ;	/* if going from left to right, get lower point */
	return (rtb) ;			/* otherwise return upper point */
}

#undef JMAX

 


/**************** Get Positions and Velocities ***********************/

/*	Requires indexed (sorted in increasing r) stars with potential
	computed in sgravity[] and N_MAX set. Uses sE[], sJ[] and sr[]
	from previous iteration. Returns positions and velocities in
	srnew[], svrnew[], and svtnew[]. Returns Max r for all stars, or
	-1 on error.
*/

double get_positions(void) {

  static long k, kmin, kmax, i, i1, si, ktemp, kk, j;	
  static double r, vr, vt, rk, rk1, a, b, rmin, rmax, E, J, Uk, Uk1 ;
  static double g1, g2, F, s0, g0, dQdr_min, dQdr_max, drds, max_rad, pot;
  static double Qtemp, X, Q ;	/* Max value of Q found. Should be > 0 */
  static long N_LIMIT, upperbound, lowerbound ;

max_rad = 0.0 ;	/* max radius for all stars, returned on success */
n_pericenter = 0 ;
n_apocenter = 0 ;
n_circular = 0 ;

phi_rtidal = potential(Rtidal) ;
phi_zero = potential(0.0) ;


if (sub_count == 0) {
	N_LIMIT = N_MAX ;		/* do FULL time step */
}
else {	
	N_LIMIT = sub_N_MAX ;	/* do sub step only */
}

for (si = 1 ; si <= N_MAX_NEW ; si++) {	/*	Repeat for all stars */

	current_star_i = si ;
	j = sindex[si] ;

	E = sE[j] ;
	J = sJ[j] ;

/* Remove stars with positive Energy, or those with mass ZERO (< 1.0e-20) due to
   stellar evolution. Note that energy lost due to stellar evolution is subtracted
   at the time of mass loss in DoStellarEvolution
*/
	if (si > N_LIMIT && si <= N_MAX && sr_peri[j] > sub_rmax) {
		
		/* do nothing, since star this is a sub timestep and star is in the halo */
	}
	else if (E >= 0.0 || sm[j] < ZERO) {	
		srnew[j] = INFINITY ;		/*escaped star */
		Eescaped += E * sm[j] / N_STAR ;
		Jescaped += J * sm[j] / N_STAR ;
		Etidal += E * sm[j] / N_STAR ;
		TidalMassLoss += sm[j] / N_STAR ;
		fprintf(escfile, \
			"%ld  %.8g  %.8g  %.8g  %.8g  %.8g %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g %ld\n", \
			tcount, TotalTime, sm[j], sInitialMass[j], sr[j], svr[j], svt[j], \
			0.0, 0.0, Rtidal, phi_rtidal, phi_zero, E, J, stype[j]) ; 

		if (stype[j] == 1) {
			mass_1 -= sm[j] ;	/* MS star lost */
			count_1-- ;
		}
		else if (stype[j] == 2) {
			mass_2 -= sm[j] ;	/* WD star lost */
			count_2-- ;
		}
		else if (stype[j] == 3) {
			mass_3 -= sm[j] ;	/* NS star lost */
			count_3-- ;
		}
		else if (stype[j] == 4) {
			mass_4 -= sm[j] ;	/* BH star lost */
			count_4-- ;
		}
		else if (stype[j] == 5) {
			mass_5 -= sm[j] ;	/* BINARY lost */
			count_5-- ;			
			M_binary -= sm[j] ;
		}
	}
	else {

/*	if((printstatus==6 || printstatus==0) && (tcount==163)) fprintf(stderr,"%ld \n",si) ;
*/

	/*	For rmin- Look for rk, rk1 such that Q(rk) < 0 < Q(rk1) */

/*		golden(1, N_MAX/3, N_MAX, &ktemp, E, J, &function_Q) ;
*/
		ktemp = si ; /* Q(si) must be +ve (selected that way in the last time step!) */

/*		for new stars due to binary disruptions, position si 
		is not ordered, so use k_old[]. 
*/
		if (si > N_MAX)	
			ktemp = old_k[j] ; 

		Qtemp = function_Q(ktemp, E, J) ;

        kk = ktemp;
		if (Qtemp <= 0) {

          lowerbound = 0;
          upperbound = N_MAX;

		  while(function_Q(kk+1,E,J) > function_Q(kk,E,J) && function_Q(kk,E,J) < 0 && kk <= N_MAX) {
		       lowerbound = kk;
                       kk++;
		  }
		  while(function_Q(kk-1,E,J) > function_Q(kk,E,J) && function_Q(kk,E,J) < 0 && kk >= 1 ) {
		    upperbound = kk;
		    kk--;	  
		  }

		  ktemp = kk;

          if ( function_Q(ktemp,E,J) < 0 ) {
		    
			/* star is on an almost circular orbit... so rmin = rmax (approx) */
			srnew[j] = sr[j] ;	
			svrnew[j] = svr[j] ;
			svtnew[j] = svt[j] ;
            fprintf(stderr, \
				"Circular orbit found: si = %ld, sr = %G, svr = %G, svt = %G, J = %G, E = %G\n", \
                si, sr[j], svr[j], svt[j], sJ[j], sE[j]);
		  }
		}
	
		else {		/*	if (Qtemp > 0)	 */

		kmin = FindZero_Q(0,ktemp,E,J);

/*		knew = zbrent_Q(0,ktemp,E,J); 
		if (abs(kmin-knew) > 1)
			fprintf(stderr,"ERROR: si = %ld   kmin = %ld  knew = %ld\n",si,kmin,knew) ;
*/
        while(function_Q(kmin, E, J) > 0 && kmin > 0) kmin--;
        while(function_Q(kmin+1,E,J) < 0 && kmin+1 < ktemp) kmin++;    

		i = sindex[kmin] ;
		i1 = sindex[kmin+1] ;
		rk = sr[i] ; 
		rk1 = sr[i1] ;
		Uk = sgravity[i] ;
		Uk1 = sgravity[i1] ;
		Q = 2.0*E - 2.0*Uk1 - J*J/(rk1*rk1) ;

		a = (Uk1 - Uk)/(1/rk1 - 1/rk) ;
		b = (Uk/rk1 - Uk1/rk)/(1/rk1 - 1/rk) ;

		rmin = J*J/(-a + sqrt(a*a - 2.0*J*J*(b - E))) ;
		dQdr_min = 2.0*J*J/(rmin*rmin*rmin) + 2.0*a/(rmin*rmin) ;

	/*	For rmax- Look for rk, rk1 such that Q(rk) > 0 > Q(rk1) */

		kmax = FindZero_Q(ktemp, N_MAX+1, E, J) ;

/*		knew = zbrent_Q(ktemp,N_MAX+1,E,J);
		if (abs(kmax-knew) > 1)
			fprintf(stderr,"ERROR: si = %ld   kmax = %ld  knew = %ld\n",si,kmax,knew) ;
*/
		while(function_Q(kmax, E, J) < 0 && kmax > ktemp) kmax--;
		while(function_Q(kmax+1,E,J) > 0 && kmax+1 < N_MAX) kmax++;    

		i = sindex[kmax] ;
		i1 = sindex[kmax+1] ;
		rk = sr[i] ; 
		rk1 = sr[i1] ;
		Uk = sgravity[i] ;
		Uk1 = sgravity[i1] ;
		Q = 2.0*E - 2.0*Uk1 - J*J/(rk1*rk1) ;

		a = (Uk1 - Uk)/(1/rk1 - 1/rk) ;
		b = (Uk/rk1 - Uk1/rk)/(1/rk1 - 1/rk) ;

		rmax = (-a + sqrt(a*a - 2.0*J*J*(b - E)))/(2.0*(b - E)) ;
		dQdr_max = 2.0*J*J/(rmax*rmax*rmax) + 2.0*a/(rmax*rmax) ;

        if( rmin > rmax ) {

			srnew[j] = sr[j] ;	
			svrnew[j] = svr[j] ;
			svtnew[j] = svt[j] ;
            fprintf(stderr, \
				"\n ERROR: rmin > rmax:  si = %ld, sr = %G, svr = %G, svt = %G, J = %G, E = %G\n", \
                si, sr[j], svr[j], svt[j], sJ[j], sE[j]);
            fprintf(stderr,"rmin = %G, rmax =  %G, kmin = %ld, kmax = %ld\n", rmin, rmax, kmin, kmax);
		} 
        else {       

		/* Check for rmax > R_MAX (tidal radius) */

		if (rmax >= Rtidal) {	
			srnew[j] = INFINITY ;		/* tidally stripped star */
			svrnew[j] = 0.0 ;
			svtnew[j] = 0.0 ;
			Eescaped += E * sm[j] / N_STAR ;
			Jescaped += J * sm[j] / N_STAR ;
			TidalMassLoss += sm[j] / N_STAR ;
			Etidal += E * sm[j] / N_STAR ;
			fprintf(escfile, \
				"%ld  %.8g  %.8g  %.8g  %.8g  %.8g %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g %ld\n", \
				tcount, TotalTime, sm[j], sInitialMass[j], sr[j], svr[j], svt[j], \
				rmin, rmax, Rtidal, phi_rtidal, phi_zero, E, J, stype[j]) ; 

			if (stype[j] == 1) {
				mass_1 -= sm[j] ;	/* MS star lost */
				count_1-- ;
			}
			else if (stype[j] == 2) {
				mass_2 -= sm[j] ;	/* WD star lost */
				count_2-- ;
			}
			else if (stype[j] == 3) {
				mass_3 -= sm[j] ;	/* NS star lost */
				count_3-- ;
			}
			else if (stype[j] == 4) {
				mass_4 -= sm[j] ;	/* BH star lost */
				count_4-- ;
			}
			else if (stype[j] == 5) {
				mass_5 -= sm[j] ;	/* BINARY lost */
				count_5-- ;
				M_binary -= sm[j] ;
			}
		}
		else {   /* NOT a tidally stripped star */ 

			g1 = sqrt(3.0*(rmax - rmin)/dQdr_min) ; /* g(-1) */
			g2 = sqrt(-3.0*(rmax - rmin)/dQdr_max) ; /* g(+1) */

			if (g2 > g1) 
				F = 1.2 * g2 ;		/* F = 1.2 * Max(g(-1),g(+1)) */
			else
				F = 1.2 * g1 ;

	/*	if (tcount==163 && si == 3353)
			fprintf(stderr, \
			"ERROR: si = %5d r = %.6g rmin = %.6G rmax = %.6g vr = %.6G X = %.6g E = %.6G J = %.6g pot = %.6G\n", \
			si, r, rmin, rmax, vr, X, E, J, pot) ;
	*/
			for (k = 1 ; k <= N_TRY ; k++) {
				X = ran2(&idum) ;

	/* for innermost 1000 stars, avoid pericenter & apocenter */
	/*			if (si < 1000) {
					if (X < 0.05) 
						X = 0.03 ;
					else if (X > 0.95)
						X = 0.97 ;
				}
	*/
				s0 = 2.0 * X - 1.0 ;	/* random -1 < s0 < 1 */
				g0 = F * ran2(&idum) ;			/* random 0 < g0 < F */

				r = 0.5*(rmin+rmax) + 0.25*(rmax-rmin)*(3.0*s0-s0*s0*s0) ;
				
			/*	if (tcount==163 && si == 3353)
					fprintf(stderr, \
					"ERROR: si = %5d r = %.6g rmin = %.6G rmax = %.6g vr = %.6G X = %.6g E = %.6G J = %.6g pot = %.6G\n", \
					si, r, rmin, rmax, vr, X, E, J, pot) ;
			*/
				pot = potential(r) ;

				drds = 0.25*(rmax-rmin)*(3.0 - 3.0*s0*s0) ;
				pot = 2.0*E - 2.0*pot - J*J/r/r ;

				if (pot >= 0.0) {
					vr = sqrt(pot) ;
				}
				else {
					fprintf(stderr, \
					"\nCircular orbit: vr^2 < 0 : setting vr = 0 in get_positions \n si = %5ld r = %.6g rmin = %.6G rmax = %.6g vr^2 = %.6G X = %.6g E = %.6G J = %.6g\n", \
					si, r, rmin, rmax, pot, X, E, J) ;
					fprintf(logfile, \
					"\nCircular orbit: vr^2 < 0 : setting vr = 0 in get_positions \n si = %5ld r = %.6g rmin = %.6G rmax = %.6g vr^2 = %.6G X = %.6g E = %.6G J = %.6g\n", \
					si, r, rmin, rmax, pot, X, E, J) ;
					vr = 0 ; 				
				}
				if (g0 < 1.0/vr * drds)		/* if g0 < g(s0) then success! */
					break ;
			}
			if ( k == N_TRY + 1) {
				fprintf(stderr,"Error: N_TRY trials exceeded in get_positions") ;
				return(-1) ;
			}

	/*------ check to see if pericenter or apocenter */
			if (X < 0.05) {		
				n_pericenter++ ;
	/*			printf("Si = %6d  Pericenter: X = %.9g  n_peri = %6d \n",si,X,n_pericenter) ;
	*/		}
			else if (X > 0.95) {
				n_apocenter++ ;			
	/*			printf("Si = %6d  Apocenter: X = %.9g  n_apo = %6d \n",si,X,n_apocenter) ;
	*/		}
			sX[j] = X ;
			sr_peri[j] = rmin ;
			sr_apo[j] = rmax ;

	/*------ pick random sign for vr */	

			if (ran2(&idum) < 0.5)		/* pick random sign for vr */
				vr = -vr ;

			vt = J / r ;

			srnew[j] = r ;
			svrnew[j] = vr ;
			svtnew[j] = vt ;

			if (r > max_rad) max_rad = r ;

		}	/* End Else ( Not a tidally stripped star ) */

	}     /* End Else ( Not the case that rmin > rmax) */

    }     /* End Else ( Not a circular orbit) */

	}	/* End Else (not escaped star) */

}	/* Next si */


#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr, 
		"New Positions Computed: Max Radius = %.8G for N_LIMIT = %ld\n", 
		max_rad, N_LIMIT) ;	
#endif

return (max_rad);

}




/* ********************  gravity  *********************

	Computing the gravity at each star indexed in sindex[] by increasing 
	radius. Units: G = 1  and  Mass is in units of total INITIAL mass.
	Total mass is computed by SUMMING over all stars that have NOT ESCAPED 
	i.e., over all stars upto N_MAX <= N_STAR. N_MAX is computed in this 
	routine by counting all stars with radius < INFINITY and Radius of the 
	(N_MAX+1)th star is set to infinity i.e., sr[sindex[N_MAX+1]] = 
	INFINITY. Also setting gravity[sindex[N_MAX+1]] = 0. Assuming 
	sr[sindex[0]] = 0. gravity[] is also indexed i.e. gravity[sindex[k]] 
	is the gravity of the kth star at radius sr[sindex[k]] 
	NOTE: Assming here that NO two stars are at the SAME RADIUS upto 
	double precision. Returns N_MAX. Gravity returned in sgravity[]
*/	

long gravity(void) {

	static long k, ii, mcount ;
	static double mprev, rtemp ;
    
	mprev = 0.0 ;

	ii = 0 ;
	for (k = 1 ; k <= N_STAR_NEW ; k++) {      /* Recompute total Mass and N_MAX */
		rtemp = sr[sindex[k]] ;
		if (rtemp < 1000000 - 1)	{
			mprev = mprev + sm[sindex[k]] ;
		}
		else {
			break ;
		}
		while (ii <= MAX_INDEX && rtemp >= (double) ii/INDEX_UNIT) {
				IndexTable[ii] = k ;
				ii++ ;
		}
	}
	N_MAX = k - 1 ;		/* New N_MAX */
	Mtotal[0] = mprev / N_STAR ;	/* New total Mass; This IS correct for multiple components */
	
	for (k=ii ; k<=MAX_INDEX ; k++) IndexTable[k] = N_MAX+1 ;

	/* Compute new tidal radius using new Mtotal */

	Rtidal = orbit_r * pow(Mtotal[0], 1.0/3.0) ; 

	sr[sindex[N_MAX+1]] = INFINITY ;
	sgravity[sindex[N_MAX+1]] = 0.0 ;
	sMi_ri[sindex[N_MAX+1]] = 0.0 ;
	sMi[sindex[N_MAX+1]] = Mtotal[0] ;
	sMi_ri[sindex[N_MAX]] = 0.0 ;
	sMi[sindex[N_MAX]] = Mtotal[0] ;
	mprev = Mtotal[0] ;
	for (k = N_MAX ; k >= 1 ; k--) {		/* Recompute gravity at star locations */
		sgravity[sindex[k]] = sgravity[sindex[k+1]] - mprev * \
			(1.0/sr[sindex[k]] - 1.0/sr[sindex[k+1]]) ;
		mprev = mprev - sm[sindex[k]] / N_STAR ;
		sMi_ri[sindex[k-1]] = sMi_ri[sindex[k]] + sm[sindex[k]]/sr[sindex[k]] / N_STAR ;
		sMi[sindex[k-1]] = sMi[sindex[k]] - sm[sindex[k]] / N_STAR ;
	}
	sgravity[sindex[0]] = sgravity[sindex[1]] ;	/* Gravity at center is the same as at first star */

#if defined(DEBUGMODE)
	if (printstatus==0 || printstatus==6) 
		fprintf(stderr,"Potential computed at stars: N_MAX = %ld \n", N_MAX) ;
#endif

/* Computing radii containing mass_pc[i] % of the mass
*/
	mprev = 0.0 ;
	mcount = 0 ;
	for (k = 1 ; k <= N_MAX ; k++) {	/* Only need to count up to N_MAX */
		mprev += sm[sindex[k]]/N_STAR ;
		if (mprev/Mtotal[0] > mass_pc[mcount]){
			mass_r[mcount] = sr[sindex[k]] ;
			mcount++ ;
			if (mcount == MASS_PC_COUNT-1) break ;
		}
	}

#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr,"Lagrange radii computed \n") ;
#endif



	return (N_MAX) ;
}





/*********************** potential **********************

  The potential computed using the sgravity[] computed at the star 
  locations in sr[] indexed in sindex[] by increasing r.
*/

double potential(double r) {
	
	static long i, k ;
	static double henon ;
	static char cont ;

/* VERY inefficient search!! replace with bisection method!
	for (i = 1 ; i <= N_MAX ; i++) {
		if (sr[sindex[i]] >= r) break ; 
	}
	kold = i - 1 ;
*/	
/*	use bisection method to find the r_k and r_k+1 that bracket r 

	i = FindZero(&function_r, 0, N_MAX+1, r, 0) ;
	if (i < 0) {
		fprintf(stderr,"Error finding zero: i = %5ld   r = %.5G\n", i, r) ;
		if(printstatus==6) gets(&cont) ;
		return (1) ;
	}
	while (sr[sindex[i]] > r){
		i-- ;
	}
	while (sr[sindex[i]] < r){
		i++ ;
	}
	kold = i - 1 ;
*/

/* New method using indexed values of sr[] & bisection */

	if (r < sr[sindex[1]]) return (sgravity[sindex[1]]) ;

	if (r*INDEX_UNIT <= MAX_INDEX) {		/* use tabulated value of r */
		i = r*INDEX_UNIT ;
		i = IndexTable[i] ;
	}
	else if (IndexTable[MAX_INDEX] == N_MAX+1){	/* ALL stars are WITHIN MAX_INDEX*1000 */
		i = N_MAX ;
	}
	else {							/* beyond tabulated values of r */
		i = FindZero(&function_r, IndexTable[MAX_INDEX]-1, N_MAX+1, r, 0) ;
		if (i < 0) {
			fprintf(logfile,"Potential: Error finding zero: i = %5ld   r = %.5G\n", i, r) ;
			fprintf(stderr,"Potential: Error finding zero: i = %5ld   r = %.5G\n", i, r) ;
			if(printstatus==6) gets(&cont) ;
			return (1) ;
		}
	}
	while (sr[sindex[i]] > r){
		i-- ;
	}
	while (sr[sindex[i]] < r){
		i++ ;
	}
	k = i - 1 ;

/* check if new bisection method gives same k or not 
	if (k != kold) {
		fprintf(stderr,"ERROR finding k in potential(): r = %.6G  k = %5ld  kold = %5ld\n",r, k, kold) ;
		if(printstatus==0) gets(&cont) ;
	}
*/
/*	Henon's method of computing the potential using sgravity[]
*/	if (k == 0) 
		henon = (sgravity[sindex[1]]) ;
	else
		henon = (sgravity[sindex[k]] + (sgravity[sindex[k+1]] \
			- sgravity[sindex[k]]) * (1.0/sr[sindex[k]] - 1.0/r) / \
			(1.0/sr[sindex[k]] - 1.0/sr[sindex[k+1]])) ;

/*	Computing potential using Phi = - Sum{i=1 to k, m_i} / r - Sum{i=k+1 to N_MAX, m_i/r_i} 

	if (k == 0) 
		kj = (- sMi_ri[sindex[k]]) ;
	else
		kj = ( -sMi[sindex[k]]/r - sMi_ri[sindex[k]] ) ;

	if (fabs(kj - henon) > 1.0e10)
		if(printstatus==6) gets(&cont) ;
	else
*/		return (henon) ;

}



/***************************** Output Density Profile **************************/
/* N.B. No longer in use (6/23/98)                                             */

void output_density_profile(void)

{
static double  r_prev, r_mean, rho_r, rho_plum, m_prev,  rmin, rmax;
static long i, X, imin, imax;

	/* else if (PLOT == 2) {*/  X=2;	  /* output density profile */

/*		r_count = 1 ;
		m_prev = 0.0 ;
		for (i=1 ; i <= N_MAX ; i++) {
*/			/* Arrays to store smoothed m,r */

/*			if (sr[sindex[i]] > r_count * R_PRINT_STEP) {
				Smooth_r[r_count] = r_count * R_PRINT_STEP ;	
				Smooth_m[r_count] = \ 
					m_prev + (r_count*R_PRINT_STEP - r_prev)/ \
					(sr[sindex[i]] - r_prev) * sm[sindex[i]] ;
				r_count++ ;	
				i-- ;	*/	/* recheck for same interval again */
/*			}
			else {
				m_prev += sm[sindex[i]] ;
				r_prev = sr[sindex[i]] ;
			}
		
*/

/* Average density over blocks of 50 stars, and also over first 10 stars */
		m_prev = 0.0 ;
		r_prev = 0.0 ; imin = 0 ; imax = 0 ;
		for (i=1 ; i <= 10 ;i++) {
			m_prev += sm[sindex[i]] / N_STAR;
			r_prev += sr[sindex[i]] ;
		}
		rmin = sr[sindex[i - 10]] ;
		rmax = sr[sindex[i-1]] ;
		r_mean = r_prev / 10.0 ;
		
		rho_r = m_prev / \
			(4.0*pi/3.0) / (pow(rmax,3.0) - pow(rmin,3.0)) ;

		rho_plum = (3.0/4.0/pi)/pow((3.0*pi/16.0),3.0) * \
			pow((1.0 + pow(r_mean/(3.0*pi/16.0),2.0)),-2.5);

		if(printstatus==X) X=0; 
                if(printstatus==0 || X==0) fprintf(out[X],"%.6G   %.6G   %.6G   %.6G \n", \
			TotalTime, r_mean, rho_r, rho_plum) ;
		m_prev = 0.0 ;
		r_prev = 0.0 ; imin = 0 ; imax = 0 ;
		for (i=1 ; i <= N_MAX ; i++) {
			if ((i % 50 == 0) || (i == N_MAX)) {
				m_prev += sm[sindex[i]] / N_STAR ;
				r_prev += sr[sindex[i]] ;
				rmin = sr[sindex[i - 49]] ;
				rmax = sr[sindex[i]] ;
				r_mean = r_prev / 50.0 ;
				rho_r = m_prev / \
					(4.0*pi/3.0) / (pow(rmax,3.0) - pow(rmin,3.0)) ;
				rho_plum = (3.0/4.0/pi)/pow((3.0*pi/16.0),3.0) * \
					pow((1.0 + pow(r_mean/(3.0*pi/16.0),2.0)),-2.5);

				if(printstatus==X) X=0; 
				if(printstatus==0 || X==0) fprintf(out[X],"%.6G   %.6G   %.6G   %.6G \n", \
					TotalTime, r_mean, rho_r, rho_plum) ;
				m_prev = 0.0 ;
				r_prev = 0.0 ;
			}
			else {
				m_prev += sm[sindex[i]] / N_STAR ;
				r_prev += sr[sindex[i]] ;
			}

		} 


/*		for (i=1 ; i <= N_MAX-2 ; i++) {
			r_mean = (sr[sindex[i+3]] - sr[sindex[i]]) ;
			rho_r = sm[sindex[i]]*(3.0 - 1.0)/(4.0*pi*pow(sr[sindex[i]], 2.0)) / r_mean ;
			rho_plum = (3.0/4.0/pi)/pow((3.0*pi/16.0),3.0) * \
				pow((1.0 + pow(sr[sindex[i]]/(3.0*pi/16.0),2.0)),-2.5);
			if(printstatus==X) X=0; 
			if(printstatus==0 || X==0) fprintf(out[X],"T = %.6G   r = %.6G    rho = %.6G   rho plum = %.6G \n", \
				TotalTime, sr[sindex[i]], rho_r, rho_plum) ;
		}
*/

}  



/*	Brent's method for finding the zero's of the function Q 
	Function Q has been HARD CODED in to the routine.	

*/

#define ITMAX 100
#define EPS 3.0e-8

long zbrent_Q(long x1, long x2, double E, double J)
{
	long iter;
	long a=x1,b=x2,c=x2,d,e,min1,min2,xm;
	double fa,fb,fc,p,q,r,s,tol1;
	
	fa = (2.0*E - 2.0*sgravity[sindex[a]] - J*J/(sr[sindex[a]]*sr[sindex[a]]));
	fb = (2.0*E - 2.0*sgravity[sindex[b]] - J*J/(sr[sindex[b]]*sr[sindex[b]]));

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5;
		xm= (0.5*(c-b));
		if (fabs(xm) < 1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=(3.0*xm*q-fabs(tol1*q));
			min2=(fabs(e*q));
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=(p/q);
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (SIGN(tol1,xm));
		fb=(2.0*E - 2.0*sgravity[sindex[b]] - J*J/(sr[sindex[b]]*sr[sindex[b]]));
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0 ;
}
#undef ITMAX
#undef EPS





/************************* mass ************************

	Computing the mass as a function of radius. Should be used 
	after gravity is calculated for new positions so that 
	sr[sindex[N_MAX+1]] = INFINITY.
	Returns mass contained within (and including) radius r for 
	r <= R_MAX.
	Returns the total (unescaped) mass for R_MAX <= r < INFINITY. 
	Returns 1 for r >= INFINITY. 
*/

double mass(double r) {

	static double masstemp ;
	static long i ;

	masstemp = 0.0 ;

	for (i=1 ; i <= N_STAR; i++) {
		if (sr[sindex[i]] > r) break ;
		masstemp = masstemp + sm[sindex[i]] ;
	}

	return (masstemp / N_STAR) ;
}

/************************* phi **************************

	Computing the gravity directly using the masses and positions, 
	without using sgravity[] computed at the star locations
	Assumes N_MAX is set correctly and that stars are indexed by radius.
*/

double phi(double r) {

	static long i, k ;
	static double U ;

	for (i = 1 ; i <= N_MAX ; i++) {
		if (sr[sindex[i]] >= r) break ; 
	}

	k = i - 1 ;
	U = 0 ;

	for (i = 1 ; i <= k ; i++) {
		U = U - sm[sindex[i]] / r ;
	}

	for (i = k+1 ; i <= N_MAX ; i++) {
		U = U - sm[sindex[i]] / sr[sindex[i]] ;
	}

	return (U / N_STAR) ;

}		


/************************** plummer *************************

	Assigning radii to stars between 0 and r_max (excluding endpoints) 
	according to a 	plummer profile. Returns error (-1) if N_STAR 
	values cannot be generated in N_TRY trials. 
	(Enables single mass-component only --- NO LONGER IN USE 2/28/99)
*/

long plummer(double r_max) {

	static double g, r, v ;
	static long i, count ;
	static double q, theta ;
	
	sr[0] = ZERO ;		/* tiny number so that 1/sr[0] is still finite */
	sr[N_STAR+1] = INFINITY ;


	for (count = 1 ; count <= N_STAR ; count++) {	

		r = 1.0 / sqrt(pow(ran2(&idum),-2.0/3.0) - 1.0) ;	
		sr[count] = r ;
	
		for (i=1 ; i <= N_TRY ; i++) {

			q = ran2(&idum) ;
			g = q*q * pow((1.0 - q*q), 7.0/2.0) ;

			if (ran2(&idum) * 0.1 < g) {	
				v = q * sqrt(2.0)*pow((1.0+r*r),(-1.0/4.0)) ;
				theta = acos(ran2(&idum)*2.0 - 1.0) ;
				svr[count] = v * cos(theta) ;	
				svt[count] = v * sin(theta) ;
				break ;
			}
		}
		if (i == N_TRY+1) return(-1) ;
	}

/*	Convert Units of r from  R  to  R / (3*Pi/16) -- ONLY FOR PLUMMER MODEL !! 
	i.e. Multiply sr[] by 3*Pi/16
*/
	  for (i = 1 ; i <= N_STAR ; i++) {
	    sr[i] = sr[i] * 3.0*pi/16.0 ;
	    svr[i] = svr[i] * sqrt(16.0/3.0/pi) ;
	    svt[i] = svt[i] * sqrt(16.0/3.0/pi) ;
	  }

	return(N_STAR) ;
}

/********************* Plummer model w/mass components *********************

Plummer models with velocities scaled according to the mass of stars
in each component, in order to ensure initial equipartition.
Corresponds to the MODEL = 0 option.

*/ 

long plummer_mass_spec(double r_max) 

{
	static double g, r, v ;
	static long i, j, count, curr_count, limit;
	static double q, theta ;

 	sr[0] = ZERO ;	/* tiny number so that 1/sr[0] is still finite */
	sr[N_STAR+1] = INFINITY ;

    curr_count = 1;
 
 	for( j=0; j<NUM_MASS; j++) {

	  limit = curr_count + ((actual_total_mass * actual_dist[j])/mass_spec[j]);
             	  
	  for (count = curr_count ; count <= limit ; count++) {	

	    r = 1.0 / sqrt(pow(ran2(&idum),-2.0/3.0) - 1.0) ;	
	    if ( j == 0 ) sr[count] = r * R1;
	    else sr[count] = r * R2;
 
	    for (i=1 ; i <= N_TRY ; i++) {

	      q = ran2(&idum) ;
	      g = q*q * pow((1.0 - q*q), 7.0/2.0) ;
	      
	      if (ran2(&idum) * 0.1 < g) {	
			v = q * sqrt(2.0)*pow((1.0+r*r),(-1.0/4.0)) ;
			theta = acos(ran2(&idum)*2.0 - 1.0) ;
			svr[count] = v * cos(theta) / sqrt(mass_spec[j] / (actual_total_mass / N_STAR));   
			svt[count] = v * sin(theta) / sqrt(mass_spec[j] / (actual_total_mass / N_STAR));
			break ;
	      }
	    }
	    if (i == N_TRY+1) return(-1) ;
	  }

          curr_count = count;
	}

/*	Convert Units of r from  R  to  R / (3*Pi/16) -- ONLY FOR PLUMMER MODEL !! 
	i.e. Multiply sr[] by 3*Pi/16
*/
	  for (i = 1 ; i <= N_STAR ; i++) {
	    sr[i] = sr[i] * 3.0*pi/16.0 ;
	    svr[i] = svr[i] * sqrt(16.0/3.0/pi) ;
	    svt[i] = svt[i] * sqrt(16.0/3.0/pi) ;
	  }

	return(N_STAR) ;
}


/************************ plummer_f_r *************************

	Plummer model radial profile. Returns Probability of having a star at 
	radius r. 	For r <= 0, returns an upper bound for the distribution 
	function (0.6 in this case).
*/

double plummer_f_r(double r) {

	r = r / (3.0*pi/16.0) ;	/* Converting from E units to R units */

	if (r > 0) {
		return (3 * r*r * 1/pow((1 + r*r), 5.0/2.0));
	}
	else {
		return (0.6) ;
	}
}


/* Golden search routine for maximum Q */

#define RC 0.61803399
#define CC (1.0-RC)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double golden(long ax, long bx, long cx,
	long *xmin, double E, double J, double (*Q)(long, double, double)) 
{
	double Q1,Q2 ;
	long x0,x1,x2,x3;

	x0=ax;
	x3=cx;
	if (labs(cx-bx) > labs(bx-ax)) {
		x1=bx;
		x2=bx+(CC*(cx-bx));
	} else {
		x2=bx;
		x1=bx-(CC*(bx-ax));
	}
	Q1=(*Q)(x1,E,J);
	Q2=(*Q)(x2,E,J);
	while ((x2 != x1) && (x1 != x0)) {
		if (Q2 > Q1) {
			SHFT3(x0,x1,x2,(RC*x1+CC*x3))
			SHFT2(Q1,Q2,(*Q)(x2,E,J))
		} else {
			SHFT3(x3,x2,x1,(RC*x2+CC*x0))
			SHFT2(Q2,Q1,(*Q)(x1,E,J))
		}
	}
	if (Q1 > Q2) {
		Q2 = Q1 ; x2 = x1 ;
	}
	else {
		Q1 = Q2 ; x1 = x2 ;
	}
	while (Q2 < 0) {
		x2++ ;
		Q2 = (*Q)(x2,E,J) ;
		if (Q2 < Q1) {
			x2 = x1 ; Q2 = Q1 ;
			break ;
		}
	}
	while (Q2 < 0) {
		x2-- ;
		Q2 = (*Q)(x2,E,J) ;
		if (Q2 < Q1) {
			x2 = x1 ; Q2 = Q1 ;
			break ;
		}
	}
	*xmin=x2;
	return Q2;

}
#undef CC
#undef RC
#undef SHFT2
#undef SHFT3



/******** Iterate computation of positions and velocities 
          to test for convergence to correct potential.
          (In fact, it does not converge!) ****************
*/

void check_potential_convergence(void)

{
 long i,j;
 double dtemp ;

 j=0;

 printf("%5ld  %8ld  %.8G  %.8G  %.8G  %.8G \n", j, N_MAX, Etotal, KEtotal, PEtotal, -2.0*KEtotal/PEtotal) ;

/* perturb initial potential 
	for (i=1; i<=N_MAX ; i++) {
		sgravity[sindex[i]] *= 1.05 ;
	}
*/

/* iterate computation of positions and velocities and gravity a few times to check for 
   convergence to correct potential*/

	for (j = 1 ; j <= 0 ; j++) {

		/*	Computing New positions and velocities for stars from sE[], sJ[]
			and sgravity[]. Returns max radius for all stars.
		*/
			dtemp = get_positions() ;
		/*	if(printstatus==0 || printstatus==6) 
				fprintf(stderr,"New Positions Computed: Max Radius = %.8G \n", dtemp) ;
		*/

		/*	Transferring new positions to sr[], svr[] and svt[] arrays 
			from srnew[], svrnew[] and svtnew[] arrays
		*/	
			for (i = 1 ; i <= N_MAX ; i++) {
				sr[sindex[i]] = srnew[sindex[i]] ;
				svr[sindex[i]] = svrnew[sindex[i]] ;
				svt[sindex[i]] = svtnew[sindex[i]] ;
			}

		/*	Indexing stars by radius. The 0th star at radius 0 
			and (N_STAR+1)th star at INFINITY are already set earlier.
		*/
			indexx(N_STAR, sr, sindex) ;
		/*	if(printstatus==0 || printstatus==6) 
				fprintf(stderr,"Stars indexed by radius \n") ;
		*/

		/*	Computing the gravity at the star locations provided in sr[].
			Results returned in sgravity[]. Returns N_MAX. 
			This is the	only place where N_MAX can change. (Also in phi(),
			but that function is only used as a check, not used in the regular 
			computation.
		*/
			errstat = gravity() ;
		/*	if(printstatus==0 || printstatus==6) 
				fprintf(stderr,"New Potential Computed: N_MAX = %8ld \n", errstat) ;
		*/

		/* recompute new sE[] and sJ[], using the new potential
		*/
		PEtotal = 0.0 ;
		KEtotal = 0.0 ;
		for (i = 1 ; i <= N_MAX ; i++) {
			 sE[sindex[i]] = sgravity[sindex[i]] 
				+ 1.0/2.0*(svr[sindex[i]]*svr[sindex[i]] + \
				svt[sindex[i]]*svt[sindex[i]]) ;
			/*
			sE[sindex[i]] = (-sMi[sindex[i]]/sr[sindex[i]] - sMi_ri[sindex[i]]) \
				+ 1.0/2.0*(svr[sindex[i]]*svr[sindex[i]] + \
				svt[sindex[i]]*svt[sindex[i]]) ;
			*/
			sJ[sindex[i]] = sr[sindex[i]] * svt[sindex[i]] ;
			KEtotal = KEtotal + 1.0/2.0*(svr[sindex[i]]*svr[sindex[i]] + \
				svt[sindex[i]]*svt[sindex[i]]) * sm[sindex[i]] / N_STAR;

			PEtotal = PEtotal + sgravity[sindex[i]] * sm[sindex[i]] / N_STAR;
			Etotal_New = Etotal_New + sE[sindex[i]] * sm[sindex[i]] / N_STAR;
		}

		PEtotal = PEtotal / 2.0 ;
		Etotal = KEtotal + PEtotal ;
		
	/* print plummer parameters for each iteration 
	*/
                printf("%5ld  %8ld  %.8G  %.8G  %.8G  %.8G \n", \
				j, N_MAX, Etotal, KEtotal, PEtotal, -2.0*KEtotal/PEtotal) ;
		

	}	/* End iteration of position/velocity computation */


#if defined(DEBUGMODE)
	if (printstatus==0 || printstatus==6) 
		fprintf(stderr,"Checked potential convergence \n") ;
#endif


}





/************************** gauss ****************************

	Gaussian distribution function
*/

double gauss(double x, double x1, double sigma) {

return(exp(-(x - x1)*(x - x1)/2.0/sigma/sigma)/2.50663/sigma) ;

}


/************************** ran2 *****************************

	Random number generator (double precision) with a long period > 2e18.
	Returns a value between 0.0 and 1.0 (exclusive of end point). Call 
	with idum a negative integer to initialize; thereafter do not alter 
	idum between successive calls in a sequence (idum can be any 
	positive long integer). Note: function ran2 requires a POINTER to a 
	long integer number as argument, not the number itself! Keeps double 
	precision upto 16 decimal places i.e. EPS = 1.2e-16
*/

#define IM1 2147483563L
#define IM2 2147483399L
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014L
#define IA2 40692L
#define IQ1 53668L
#define IQ2 52774L
#define IR1 12211L
#define IR2 3791L
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-16
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789L;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/***************************** func ****************************

	Function whose root is to be found using rtbis (bisection method) 
	after bracketing
*/

double func(double x)
{

/* find root (x) of  G/k *  m1*2.0e33 * x^3/(1+x)^2 -1 = 0  */

return(6.71e-8 / fabs(2.0)* 1.7 * 2.0e33 * x*x*x / (1.0+x) /  \
	   (1.0+x) - 1.0 ) ;

}




/**************************** indexx ***************************

	Indexing routine used to index the stars by increasing radius
*/

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(long n, double arr[], long indx[])
{
	long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);

#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr,"Stars indexed by radius \n") ;
#endif

}

#undef M
#undef NSTACK
#undef SWAP


/************* Data smoothing routines *******************/

#define NRANSI

void savgol(double c[], long np, long nl, long nr, long ld, long m)
{
	void lubksb(double **a, long n, long *indx, double b[]);
	void ludcmp(double **a, long n, long *indx, double *d);
	long imj,ipj,j,k,kk,mm,*indx;
	double d,fac,sum,**a,*b;

	if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m)
	nrerror("bad args in savgol");
	indx=ivector(1,m+1);
	a=matrix(1,m+1,1,m+1);
	b=vector(1,m+1);
	for (ipj=0;ipj<=(m << 1);ipj++) {
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
		for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
		mm=FMIN(ipj,2*m-ipj);
		for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	}
	ludcmp(a,m+1,indx,&d);
	for (j=1;j<=m+1;j++) b[j]=0.0;
	b[ld+1]=1.0;
	lubksb(a,m+1,indx,b);
	for (kk=1;kk<=np;kk++) c[kk]=0.0;
	for (k = -nl;k<=nr;k++) {
		sum=b[1];
		fac=1.0;
		for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
		kk=((np-k) % np)+1;
		c[kk]=sum;
	}
	free_vector(b,1,m+1);
	free_matrix(a,1,m+1,1,m+1);
	free_ivector(indx,1,m+1);
}
#undef NRANSI



#define NRANSI
#define TINY 1.0e-20;

void ludcmp(double **a, long n, long *indx, double *d)
{

long i,imax,j,k;
double big,dum,sum,temp;
double *vv;
vv=vector(1,n);
*d=1.0;

for (i=1;i<=n;i++) {
	big=0.0;
	for (j=1;j<=n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
	if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
	vv[i]=1.0/big;
}

for (j=1;j<=n;j++) {
	for (i=1;i<j;i++) {
		sum=a[i][j];
		for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
	}
	big=0.0;
	for (i=j;i<=n;i++) {
		sum=a[i][j];
		for (k=1;k<j;k++)
		sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		if ( (dum=vv[i]*fabs(sum)) >= big) {
			big=dum;
			imax=i;
		}
	}
	if (j != imax) {
		for (k=1;k<=n;k++) {
			dum=a[imax][k];
			a[imax][k]=a[j][k];
			a[j][k]=dum;
		}
		*d = -(*d);
		vv[imax]=vv[j];
	}
	indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=TINY;

if (j != n) {
	dum=1.0/(a[j][j]);
	for (i=j+1;i<=n;i++) a[i][j] *= dum;
}

}

free_vector(vv,1,n);

}

#undef TINY
#undef NRANSI

void lubksb(double **a, long n, long *indx, double b[])
{

long i,ii=0,ip,j;
double sum;

for (i=1;i<=n;i++) {
	ip=indx[i];
	sum=b[ip];
	b[ip]=b[i];
	if (ii)
		for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
	else if (sum) ii=i;
	b[i]=sum;
}

for (i=n;i>=1;i--) {
	sum=b[i];
	for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
	b[i]=sum/a[i][i];
}

}


/******************************************************************************/


/*** Componentwise KE, density, lagrangian radii, mass spectrum, and file I/O routines */
/***  --- Wesley Andres Watters Farfan (watters@mit.edu)	 **********/


void ComputeComponentKE(void) {


	long i, jj, whichmass, check_core ;
	double mtemp, vtemp, rhotemp ;

/*=================== Calculation and printing of averaged kinetic energies, 
                      densities, and the number of stars contained within the 
                      total cluster core-radius, and the 25%, 50%, and 75% 
                      lagrange radii for each component =====================================*/

        mtemp=0.0;
        vtemp=0.0;
        check_core = 0;  /* Used to check that core densities and temperature are 
                            calculated only once. */

        /* ---- Calculate core radius for the whole cluster ---- */

        for(i=1; i <= NUM_CORE_STARS; i++) {
            mtemp += sm[sindex[i]]*(actual_total_mass/N_STAR)/N_STAR;
            vtemp += pow(svr[sindex[i]],2.0) + pow(svt[sindex[i]],2.0); 
		}

        rhotemp = mtemp / (4*pi*(pow(sr[sindex[NUM_CORE_STARS+1]],3.0))/3);
        vtemp = sqrt(vtemp / NUM_CORE_STARS);

        r_core = sqrt( (3 * pow(vtemp,2)) / (4*pi*rhotemp) );

        /* --- Calculate three lagrange radii per component, densities, and temperatures  --- */
       
        for (jj = 0; jj < NUM_MASS; jj++) {
            whichrad[jj] = 1;   /* Used to keep track of which radii is being calculated */
            countmass[jj] = 0;
            velaver[jj] = 0;
	}
     
	for(i = 1; i <= N_MAX; i++) {

           if(MODEL==1) whichmass = find_position(sm[sindex[i]]*(TOTAL_MASS/N_STAR),mass_spec,NUM_MASS);
             else  whichmass = find_position(sm[sindex[i]]*(actual_total_mass/N_STAR),mass_spec,NUM_MASS);

           if ( sr[sindex[i]] > r_core && check_core == 0 ) {
 
               for (jj = 0; jj < NUM_MASS; jj++) {

                   if(countmass[jj] != 0) temp_within_core[jj] =  velaver[jj] /  countmass[jj]; 
                     else temp_within_core[jj] = INFINITY; 
                   dens_within_core[jj] =  countmass[jj] / ((4*pi*pow(sr[sindex[i]],3))/3);  

                   num_within_core[jj] = countmass[jj];
               	   }

              check_core = 1;
             
	   }

           if ( countmass[whichmass]*sm[sindex[i]]/N_STAR > whichrad[whichmass]*0.25*Mtotal[whichmass+1]) {

              dens_veloc_radii[3*whichmass + (whichrad[whichmass]-1)] = sr[sindex[i]];
              
              /*  Data structure: component, component of mass radius, index (1-3) of mass radius */

	      for (jj = 0; jj < NUM_MASS; jj++) {

		if (countmass[jj] != 0) temp_within[jj*3*NUM_MASS + 3*whichmass + (whichrad[whichmass]-1)] = velaver[jj] /  countmass[jj]; 
                   else temp_within[jj*3*NUM_MASS + 3*whichmass + (whichrad[whichmass]-1)] = INFINITY; 
		dens_within[jj*3*NUM_MASS + 3*whichmass + (whichrad[whichmass]-1)] = countmass[jj] / ((4*pi*pow(sr[sindex[i]],3))/3);

		num_within[jj*3*NUM_MASS + 3*whichmass + (whichrad[whichmass]-1)] = countmass[jj];
              		       
	      }
              whichrad[whichmass]++;    
	   }

          countmass[whichmass]++;     
          velaver[whichmass] += 0.5 * sm[sindex[i]] * (pow(svr[sindex[i]],2.0) + pow(svt[sindex[i]],2.0)); 
           
	}

#if defined(DEBUGMODE)
	if (printstatus==0 || printstatus==6) 
		fprintf(stderr,"Computed KE for each component \n") ;
#endif

}



void PrintFileOutput(void) {

	static long i, j, X, jj, kk, StepCount ;

	X=1;	/* print lagrangian radii */


/**** KJ modified.... set NUM_MASS_RADII_BINS = 0 to get only one 
      set of radii for the total mass. Other mass bins causing trouble.
*/
	NUM_MASS_RADII_BINS = 0 ;
	
/** Also note that there are only MASS_PC_COUNT-1 radii from 0...MASS_PC_COUNT-1 
    So output only needs to go over i = 0 to MASS_PC_COUNT-2. 
**/
/****************************************/


	for(j=0;j<NUM_MASS_RADII_BINS+1;j++) {

		if(printstatus==X) X=0; 
        if(printstatus==0 || X==0) fprintf(out[X],"%ld  %.6G  ",j,TotalTime);

		for (i=0 ; i < MASS_PC_COUNT-1 ; i++) {
			if(printstatus==X) X=0; 
			if(printstatus==0 || X==0) fprintf(out[X],"%.6G   ", mass_r[j*MASS_PC_COUNT+i]) ;
		}
		if(printstatus==X) X=0; 
		if(printstatus==0 || X==0) fprintf(out[X],"\n") ;
	}

	X=3;	    /* output Time,N_MAX,TotalE,TotalKE,TotalPE,Mtotal[0]*/

	if(printstatus==X) X=0; 
	if(printstatus==0 || X==0) 
		fprintf(out[X], \
		"%.8G  %8ld  %.8G  %.8G  %.8G  %.8G  %.8G  %.8G  %.8G  %8ld  %.8G  %.8G  %.8G  %.8G  %.8G  %.8G\n", \
		TotalTime, N_MAX, Etotal, KEtotal, PEtotal, Mtotal[0], Etotal_New, Eescaped, Jescaped,\
		tcount, DEavg, DEmax, DEdiff, max_r, core_radius, N_core);	



	X=5;		/* print temperatures & densities */

	if(printstatus==X) X=0; 
	if(printstatus==0 || X==0) {
	  
	  for (jj = 0; jj < NUM_MASS; jj++) {
	    
		  for (i = 0; i < NUM_MASS; i++) {
	      
			fprintf(out[X],"%G %ld %ld %G %G %G %G ", \
			  TotalTime, jj+1, i+1, r_core, num_within_core[jj], dens_within_core[jj], \
			  temp_within_core[jj]); 
	      
			for(kk= 0; kk<3; kk++) 
			  fprintf(out[X],"%G %G %G %G ", \
				dens_veloc_radii[3*i+kk], num_within[jj*NUM_MASS*3+3*i+kk], \
				dens_within[jj*NUM_MASS*3+3*i+kk], temp_within[jj*NUM_MASS*3+3*i+kk]);
	      
			fprintf(out[X],"\n");

/* Output format is: 
  Time, mass component A, mass component B used to compute the 3 lagrange radii, 
  number of A within core radius, core density of A, core temperature of A,
  number of A within 25% lagrange radius of B, density of A within 25% of B, temperature for A within 25% of B  
  number of A within 50% lagrange radius of B, density of A within 50% of B, temperature for A within 50% of B
  number of A within 75% lagrange radius of B, density of A within 75% of B, temperature for A within 75% of B
*/ 
		 }
	  }

	}

/* ---------- output SNAPSHOT for computing density & velocity profiles ------------ */

	X=2;	/* output density profile */

/* Set StepCount=1 in the first timestep 
	if (StepCount == 0) StepCount = 1 ;		
*/
	if ((TotalTime-StartTime) >= T_PRINT_STEP * StepCount) {	/* also saves INITIAL snapshot (StepCount=0) */
		if (StepCount > 0 || RESTART == 0) {	/* Don't save initial snapshot if Restarting */
			if(DUMPS==1) output_rawdump();		/* outputs both DENSITY and RESTART snapshots */
		}
		StepCount++ ;	
	}

}



void CountComponentMass(void) {

/*==================== Counting the mass in each component ============================

     Counts the total mass in each component, as a fraction of the 
     initial mass in each component.  Same procedure occurs once later.
     Needs to be binned if you are using the continuous mass spectrum of MODEL 4 */

  long i, X;
  long whichmass;

  /******** Added by Cody for Model 4 ************/
  double Max, step, curr_mass, j;
  long k;

  if(MODEL==4) {
    max=find_max(mass_spec, NUM_MASS);
    Max=max*(N_STAR/actual_total_mass);
    
    step=Max/NUM_MASS_RADII_BINS;

    for(i=1; i<=NUM_MASS_RADII_BINS; i++) Mtotal[i] = 0.0;
    for(i=1; i<=N_MAX; i++) {
      curr_mass=sm[sindex[i]];
      j=0.0;
      k=1;
      /* No equal signs here because ran2 is exclusive of 0 and 1 */
	  while( !( (curr_mass >= j) && (curr_mass < (j+step)))) {
		j+=step;
		k++;
		/*	printf("k=%ld curr_mass is %g\n", k, curr_mass);*/
	  }
      Mtotal[k]+=sm[sindex[i]];
    }
    for(i=1; i<=NUM_MASS_RADII_BINS; i++) Mtotal[i] = Mtotal[i]/N_STAR;

    X = 4;
    if(printstatus==X) X=0; 
    if(printstatus==0 || X==0) {
      for(i = 0; i<= NUM_MASS_RADII_BINS; i++) {
	fprintf(out[X],"%ld %G %G",i,TotalTime,Mtotal[i]); 
	fprintf(out[X],"\n");
      }
    }
    return;  
  }
  /********* End of Cody's additions **************/     
  
  
  for(i = 1; i <= NUM_MASS; i++) Mtotal[i] = 0.0;
  
  for(i = 1 ; i <= N_MAX; i++) {
    
    if(MODEL==1) 
      whichmass = find_position(sm[sindex[i]]*(TOTAL_MASS/N_STAR),mass_spec,NUM_MASS);
    else 
      whichmass = find_position(sm[sindex[i]]*(actual_total_mass/N_STAR),mass_spec,NUM_MASS);
    
    Mtotal[whichmass+1] =  Mtotal[whichmass+1] + sm[sindex[i]];
  }
  
  for(i = 1; i <= NUM_MASS; i++) Mtotal[i] = Mtotal[i]/N_STAR;
  
  X = 4;
  if(printstatus==X) X=0; 
  if(printstatus==0 || X==0) {
    for(i = 0; i<= NUM_MASS; i++) {
      fprintf(out[X],"%ld %G %G",i,TotalTime,Mtotal[i]); 
      fprintf(out[X],"\n");
    }
  }

/*=====================================================================================*/


#if defined(DEBUGMODE)
	if (printstatus==0 || printstatus==6) 
	  printf("CountComponentMass() done\n");
#endif


}


void GetInitialModel(void) {

/*============================ Generating initial models =============================*/
/*
	MODEL = 0 ---> Initial Plummer model with user-defined mass-spectra.
	MODEL = 1 ---> Initial Starlab-format model in INPUT_SNAPSHOT_FILE
         	       (with starlab masses assigned)
	MODEL = 2 ---> Initial Starlab-format superposed (starlab#.dat) models with 
               	       masses assigned according to specifications in input file.    
	MODEL = 3 ---> Starlab ``snapshot'' with single mass component, and mass spectra 
		       assigned internally according to input-file parameters.
	MODEL = 4 ---> Starlab "snapshot" with single mass component, add mass spectra assigned
	               with N ~ m^ALPHA where ALPHA is an input parameter


*/
/*======= Assigning masses to stars according to user-specified parameters ============*/

	if(MODEL == 0) {
	  mass_spec_assign(1);	/* NOTE: sm[0]=0 and sm[N_STAR+1]=0  */
	  errstat = plummer_mass_spec(R_MAX) ; /* NOTE: sr[0]=ZERO and sr[N_STAR+1]=INFINITY */
	  if (printstatus==0 || printstatus==6) {
	    fprintf(stderr,"MODEL = %ld - Created Plummer model...%ld stars\n", MODEL, errstat) ;
	  }
	} 
	else if (MODEL == 1 || MODEL == 3) {
	  errstat = starlab_model_fetch_with_masses(); /* NOTE: sr[0]=ZERO and sr[N_STAR+1]=INFINITY */
	  if (MODEL != 1) {
	    mass_spec_assign(1); /* NOTE: sm[0]=0 and sm[N_STAR+1]=0  */
	  }
	  if (printstatus==0 || printstatus==6) {
	    fprintf(stderr,"MODEL = %ld - Retrieved Starlab model with masses...%ld stars\n",MODEL, errstat) ;
	  }
	}
	else if (MODEL == 2) {
	  errstat = starlab_model_spec_fetch(); /* NOTE: sr[0]=ZERO and sr[N_STAR+1]=INFINITY */
	  mass_spec_assign(1); /* NOTE: sm[0]=0 and sm[N_STAR+1]=0  */
	  if (printstatus==0 || printstatus==6) {
	    fprintf(stderr,"MODEL = %ld - Retrieved Starlab mass-component models...%ld stars\n",MODEL, errstat) ;
	  }
	}
	else if (MODEL==4) {	/* continuous power law m^(alpha) mass distribution */
	  errstat = starlab_model_fetch_with_masses();
	  if (printstatus==0 || printstatus==6) {
	    fprintf(stderr,"MODEL = %ld - Retrieved starlab model with masses...%ld stars\n", MODEL, errstat) ;
	  }
	  initial_mass_powerlaw(ALPHA);
	  if (printstatus==0 || printstatus==6) {
	    fprintf(stderr,"Power law mass spectrum assigned\n") ;
	  }
	}

	/* Outputs masses to outputfile4 --- No longer in use (7/18/99)
	if (printstatus==0) {
	for(i=1; i<=N_STAR;i++) fprintf(out[4],"%G\n",sm[i]);
	fclose(out[4]);
	} 
	*/

}


void ComputeLagrangianRadii(void) {

  static double *mprevv, Min, Max, Range ;
  static double step,j;
  static long i, k, *mcountt, l;     
  
  mprevv=(double *)malloc((NUM_MASS_RADII_BINS+1)*sizeof(double));
  mcountt=(long *)malloc((NUM_MASS_RADII_BINS+1)*sizeof(int));
  
  /* Added by Cody to fix for initial_mass_powerlaw() */
  max=find_max(mass_spec, NUM_MASS);
  min=find_min(mass_spec, NUM_MASS);
  

  Min=min*(N_STAR/actual_total_mass);
  Max=max*(N_STAR/actual_total_mass);

  if(MODEL == 1) {
    max=find_max(mass_spec, NUM_MASS);
    min=find_min(mass_spec, NUM_MASS);
    Min=min*(N_STAR/TOTAL_MASS);
    Max=max*(N_STAR/TOTAL_MASS);        
  }
  
  Range=Max;
  
  /* Computing radii containing mass_pc[i] % of the mass (LAGRANGE RADII calculations)
   */
  step=Range/NUM_MASS_RADII_BINS;
  
  if(MODEL==1) step = Range/NUM_MASS_RADII_BINS + 0.25;  /* Note: 0.25 here is an 
							    adjustable parameter that 
							    enables easy placement of 
							    individual masses in the 
							    mass-spectrum bins for 
							    MODEL = 1. 
							    Adjust of lagrange radii 
							    calculations malfunction. */
  
  for(i=0;i<(NUM_MASS_RADII_BINS+1);i++) {
    mprevv[i] = 0.0 ;
    mcountt[i] = 0 ;
    for(k=0;k<MASS_PC_COUNT;k++) mass_r[i*MASS_PC_COUNT+k]=1.0;
    
  }
  
  l=find_position(Max*(actual_total_mass/N_STAR),mass_spec,NUM_MASS);
  if(MODEL==1) l=find_position(Max*(TOTAL_MASS/N_STAR),mass_spec,NUM_MASS); 
  
  
  for (k = 1 ; k <= N_MAX ; k++) {	/* Only need to count up to N_MAX */
    
    /* The following conditional and for loop finds in which mass
       bin the current mass resides, and calls this i 
    */
    
    j=0;
    for(i=0;j<sm[sindex[k]];i++)j=j+step;
    i--;
    
    
    /* The MASS_PC_COUNT Lagrange radii are assigned to mass_r[m] 
       for m=1,...,NUM_MASS_BINS, and the total (mass-insensitive)
       Lagrange radii are assigned to mass_r[0] 
       */
    
    for(j=0;j<2;j++) {   /* Calculation for individual mass component is 
			    made first, and then repeated for the total mass 
			    (i.e., iteration, j = 1, for which the bin/index 
			    ``i'' is set equal to -1, e.g.,  so that 
			    Mtotal[i+1] = Mtotal[0]) */  
      if(j==1) i=-1; 
      
      /* The following keeps track of mass as fraction of the total mass 
       */
      mprevv[i+1] = mprevv[i+1] + sm[sindex[k]] / N_STAR;   
      
      
      /* In the following, if the fraction of total mass in the current component i
	 exceeds the current lagrange radius mass percentage (mass_pc[mcountt[i+1]), 
	 then... 
	 */
      if (mprevv[i+1]/Mtotal[i+1] > mass_pc[mcountt[i+1]]){  
	
	/* As long as this isn't the last (outermost) lagrange radius for the 
	   current component, then go ahead and choose the radius of the 
	   current star for the current lagrange radius for the current 
	   mass component (mcountt[i+1]) 
	*/
	
	if(mcountt[i+1]!=MASS_PC_COUNT){ 
	  mass_r[(i+1)*MASS_PC_COUNT+mcountt[i+1]] = sr[sindex[k]] ;
	  mcountt[i+1]++ ;
	}
	
	/* The following lines ensure that the radius of the current
	   star, if it contains a mass percentage larger than that of 
	   subsequent lagrange radii, will be assigned to the value
	   of those radii. 
	   */
	
	while(mcountt[i+1]<MASS_PC_COUNT-1) {
	  if (mprevv[i+1]/Mtotal[i+1] > mass_pc[mcountt[i+1]]) {
	    if(mcountt[i+1]!=MASS_PC_COUNT){ 
	      mass_r[(i+1)*MASS_PC_COUNT+mcountt[i+1]] = sr[sindex[k]] ;  
	      mcountt[i+1]++ ;
	    }
	  } else break;
	}
	
      }
    }  
  }

#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr,"Lagrangian radii computed \n") ;
#endif

}


/***************** Parsing of Input Parameters / Memory allocation / File I/O **************/


int parser(int argc, char *argv[])

{

static  double normalization, difference;
static  char inputfile[100], outfile[100], outfilemode[5];
static  char parameter_name[400], values[400];
static  char *curr_mass, *curr_param;
static  int i;
static  FILE *in; 

/*======= Checking & storing of command-line arguments =======*/

	if(argc!=4) { 
	  printf("\n\n Usage: cluster [input_file] [output_file_prefix] [print_status]\n\n");
	  printf(" where [print_status] = {0,1,3,4,5}\n");
	  printf(" i.e., [print_status] = 0 --> Batch mode: logs & data to screen + 5 ouput files.\n");
	  printf("       [print_status] = {1,3,4,5} --> Prints to screen the data types 1,3,4,5.\n"); 
          printf("       [print_status] = 6 ---> Interactive mode: logs & data to screen w/pause.\n\n");
          printf(" The data types/files are labeled as follows: \n\n");
          printf("       [1] = Lagrange radii\n");
          printf("       [2] = Binary system snapshots\n");
          printf("       [3] = Total energy, virial ratio\n");
          printf("       [4] = Component-wise fraction of total mass\n");
          printf("       [5] = Temperature & Density inside 7 radii\n\n\n");
	  return 0;
	}

	if (myid != 1000) {
		sprintf(inputfile,"%s_%d",argv[1],myid); 
		sprintf(outprefix,"%s_%d_",argv[2],myid);
		printstatus=argv[3][0]-'0';
	}
	else{
		sprintf(inputfile,"%s",argv[1]); 
		sprintf(outprefix,"%s_",argv[2]);
		printstatus=argv[3][0]-'0';
	}
  
/*======= Opening of input & output files =======*/

	if ((in = fopen(inputfile, "r"))==NULL){      
	  printf("\n\nCannot open input file...\n");
	  return 0;
	}


/*======= Reading of parameter values from input file =======*/

	while(fscanf(in,"%s %s",parameter_name, values)>0) {

	    if(strcmp(parameter_name,"N_STAR")==0) { sscanf(values,"%ld",&N_STAR); N_STAR_NEW=N_STAR; N_STAR_DIM=N_STAR*2;}
	    else if(strcmp(parameter_name,"PERTURB")==0) sscanf(values,"%ld",&PERTURB); 
	    else if(strcmp(parameter_name,"MODEL")==0) sscanf(values,"%ld",&MODEL);
	    else if(strcmp(parameter_name,"N_TRY")==0) sscanf(values,"%ld",&N_TRY);
	    else if(strcmp(parameter_name,"R_MAX")==0) sscanf(values,"%lf",&R_MAX); 
	    else if(strcmp(parameter_name,"T_MAX_COUNT")==0) sscanf(values,"%ld",&T_MAX_COUNT);
	    else if(strcmp(parameter_name,"T_PRINT_STEP")==0) sscanf(values,"%lf",&T_PRINT_STEP);
	    else if(strcmp(parameter_name,"R_PRINT_STEP")==0) sscanf(values,"%lf",&R_PRINT_STEP);
	    else if(strcmp(parameter_name,"T_MAX")==0) sscanf(values,"%lf",&T_MAX);
	    else if(strcmp(parameter_name,"TOTAL_MASS")==0) sscanf(values,"%lf",&TOTAL_MASS);
	    else if(strcmp(parameter_name,"NUM_CORE_STARS")==0) sscanf(values,"%ld",&NUM_CORE_STARS);
	    else if(strcmp(parameter_name,"TERMINAL_ENERGY_DISPLACEMENT")==0) sscanf(values,"%lf",&TERMINAL_ENERGY_DISPLACEMENT);
	    else if(strcmp(parameter_name,"DUMP_ENERGY_DISPLACEMENT")==0) sscanf(values,"%lf",&DUMP_ENERGY_DISPLACEMENT);
	    else if(strcmp(parameter_name,"TRC_FACTOR")==0) sscanf(values,"%lf",&TRC_FACTOR);
	    else if(strcmp(parameter_name,"SIN2BETA_MAX")==0) sscanf(values,"%lf",&SIN2BETA_MAX);
	    else if(strcmp(parameter_name,"DT_MAX")==0) sscanf(values,"%lf",&DT_MAX);
	    else if(strcmp(parameter_name,"DT_MIN")==0) sscanf(values,"%lf",&DT_MIN);
	    else if(strcmp(parameter_name,"DUMPS")==0) sscanf(values,"%ld",&DUMPS);
	    else if(strcmp(parameter_name,"E_CONS")==0) sscanf(values,"%ld",&E_CONS);
	    else if(strcmp(parameter_name,"DT_MODE")==0) sscanf(values,"%ld",&DT_MODE);
	    else if(strcmp(parameter_name,"MAX_INDEX")==0) sscanf(values,"%ld",&MAX_INDEX);
	    else if(strcmp(parameter_name,"INDEX_UNIT")==0) sscanf(values,"%ld",&INDEX_UNIT);
	    else if(strcmp(parameter_name,"INPUT_SNAPSHOT_FILE")==0) sscanf(values,"%s",INPUT_SNAPSHOT_FILE);
	    else if(strcmp(parameter_name,"E_NORMALIZE")==0) sscanf(values,"%ld",&E_NORMALIZE);
	    else if(strcmp(parameter_name,"MIN_LAGRANGIAN_RADIUS")==0) sscanf(values,"%lf",&MIN_LAGRANGIAN_RADIUS);
	    else if(strcmp(parameter_name,"ALPHA")==0) sscanf(values,"%lf",&ALPHA);
	    else if(strcmp(parameter_name,"T_RELAX")==0) sscanf(values,"%lf",&T_RELAX);
	    else if(strcmp(parameter_name,"FAMILY")==0) sscanf(values,"%lf",&FAMILY);
	    else if(strcmp(parameter_name,"RESTART")==0) sscanf(values,"%ld",&RESTART);
	    else if(strcmp(parameter_name,"STELLAR_EVOLUTION")==0) sscanf(values,"%ld",&STELLAR_EVOLUTION);
	    else if(strcmp(parameter_name,"TRACE_N")==0) sscanf(values,"%ld",&TRACE_N);
	    else if(strcmp(parameter_name,"TRACE_MASS1")==0) sscanf(values,"%lf",&TRACE_MASS1);
	    else if(strcmp(parameter_name,"TRACE_MASS2")==0) sscanf(values,"%lf",&TRACE_MASS2);
		else if(strcmp(parameter_name,"TRACE_STARS")==0) sscanf(values,"%ld",&TRACE_STARS);	
	    else if(strcmp(parameter_name,"ORBIT_ELLIPTICAL")==0) sscanf(values,"%ld",&ORBIT_ELLIPTICAL);
	    else if(strcmp(parameter_name,"ORBIT_ECCENTRICITY")==0) sscanf(values,"%lf",&ORBIT_ECCENTRICITY);
	    else if(strcmp(parameter_name,"ORBIT_PERIOD")==0) sscanf(values,"%lf",&ORBIT_PERIOD);
	    else if(strcmp(parameter_name,"ORBIT_PHASE")==0) sscanf(values,"%lf",&ORBIT_PHASE);
	    else if(strcmp(parameter_name,"ORBIT_MAX_DT")==0) sscanf(values,"%lf",&ORBIT_MAX_DT);
	    else if(strcmp(parameter_name,"ORBIT_PERI")==0) sscanf(values,"%lf",&ORBIT_PERI);
	    else if(strcmp(parameter_name,"ORBIT_VC")==0) sscanf(values,"%lf",&ORBIT_VC);
	    else if(strcmp(parameter_name,"N_BINARY")==0) sscanf(values,"%ld",&N_BINARY);
		else if(strcmp(parameter_name,"IDUM")==0) {
                   sscanf(values,"%ld",&IDUM);
                   idum=IDUM;
	    }
        else if(strcmp(parameter_name,"MASS_SPEC")==0) {
               strcpy(MASS_SPEC,values);
               curr_mass=(char *)strtok(values,",; ");
               for(NUM_MASS=1; (curr_mass=(char *)strtok(NULL," ,;"))!=NULL;NUM_MASS++);  
        }
        else if(strcmp(parameter_name,"MASS_PC")==0) {
               strcpy(MASS_PC,values);
               curr_mass=(char *)strtok(values,",; ");
               sscanf(curr_mass,"%ld",&NUM_MASS_RADII_BINS);
               for(MASS_PC_COUNT=1; (curr_mass=(char *)strtok(NULL," ,;"))!=NULL; MASS_PC_COUNT++);  
	    }
        else if(strcmp(parameter_name,"MASS_DIST_PARAMS")==0) {
               strcpy(MASS_DIST_PARAMS,values);
               curr_mass=(char *)strtok(values,",; ");
               for(TOTAL_PARAMS=1; (curr_mass=(char *)strtok(NULL," ,;"))!=NULL; TOTAL_PARAMS++);  
	    }
	}

	fclose(in);
    
	if(IDUM==0) idum=(-115);
        
/*======= Allocation of memory for global variables using parameter values =======*/

 	sr=(double *)malloc(N_STAR_DIM*sizeof(double));
	svr=(double *)malloc(N_STAR_DIM*sizeof(double));
	svt=(double *)malloc(N_STAR_DIM*sizeof(double));
	sm=(double *)malloc(N_STAR_DIM*sizeof(double));
	sLifeTime=(double *)malloc(N_STAR_DIM*sizeof(double));
	sInitialMass=(double *)malloc(N_STAR_DIM*sizeof(double));
	sRemnantMass=(double *)malloc(N_STAR_DIM*sizeof(double));

	sE=(double *)malloc(N_STAR_DIM*sizeof(double));
	sJ=(double *)malloc(N_STAR_DIM*sizeof(double));
	sgravity=(double *)malloc(N_STAR_DIM*sizeof(double));
	sEI=(double *)malloc(N_STAR_DIM*sizeof(double));
	sJI=(double *)malloc(N_STAR_DIM*sizeof(double));

	sKE=(double *)malloc(N_STAR_DIM*sizeof(double));
	sPE=(double *)malloc(N_STAR_DIM*sizeof(double));
	sMi_ri=(double *)malloc(N_STAR_DIM*sizeof(double));
	sMi=(double *)malloc(N_STAR_DIM*sizeof(double));

	srescape=(double *)malloc(N_STAR_DIM*sizeof(double));
	mass_r=(double *)malloc((NUM_MASS_RADII_BINS+1)*MASS_PC_COUNT*sizeof(double));

	srnew=(double *)malloc(N_STAR_DIM*sizeof(double));
	svrnew=(double *)malloc(N_STAR_DIM*sizeof(double));
	svtnew=(double *)malloc(N_STAR_DIM*sizeof(double));

	srOld=(double *)malloc(N_STAR_DIM*sizeof(double));
	sX = (double *)malloc(N_STAR_DIM*sizeof(double));
	sDE = (double *)malloc(N_STAR_DIM*sizeof(double));
	sDJ = (double *)malloc(N_STAR_DIM*sizeof(double));
	sSin2Beta = (double *)malloc(N_STAR_DIM*sizeof(double));
	sr_peri = (double *)malloc(N_STAR_DIM*sizeof(double));
	sr_apo = (double *)malloc(N_STAR_DIM*sizeof(double));
	sl=(double *)malloc(N_STAR_DIM*sizeof(double));

/******** binary parameters *******/
	by=(double *)malloc(N_STAR_DIM*sizeof(double));
	ba=(double *)malloc(N_STAR_DIM*sizeof(double));
	be=(double *)malloc(N_STAR_DIM*sizeof(double));
	bm1=(double *)malloc(N_STAR_DIM*sizeof(double));
	bm2=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_a=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_e=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_m1=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_m2=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_tdestroyed=(double *)malloc(N_STAR_DIM*sizeof(double));
	bini_tformed=(double *)malloc(N_STAR_DIM*sizeof(double));
	s_binindex=(long*)malloc(N_STAR_DIM*sizeof(long));
	b_sindex=(long*)malloc(N_STAR_DIM*sizeof(long));
	sinteracted=(long*)malloc(N_STAR_DIM*sizeof(long));
/***********************************/

	sindex=(long*)malloc(N_STAR_DIM*sizeof(long));
	Smooth_r=(double *)malloc(N_STAR_DIM*sizeof(double));
	Smooth_m=(double *)malloc(N_STAR_DIM*sizeof(double));
	old_k=(long*)malloc(N_STAR_DIM*sizeof(long));

	trace=(long*)malloc((TRACE_N+5)*sizeof(long));
	trace_stat=(long*)malloc((TRACE_N+5)*sizeof(long));
	trace_r_ini=(double *)malloc((TRACE_N+5)*sizeof(double));
	trace_tr=(double *)malloc((TRACE_N+5)*sizeof(double));
	stype=(long*)malloc(N_STAR_DIM*sizeof(long));

	mass_spec=(double *)malloc(NUM_MASS*sizeof(double));
	mass_pc=(double *)malloc(MASS_PC_COUNT*sizeof(double));
	dist_params=(double *)malloc(TOTAL_PARAMS*sizeof(double));
	actual_dist=(double *)malloc(NUM_MASS*sizeof(double));

	if(MODEL==4) {
		Mtotal=(double *)malloc((NUM_MASS_RADII_BINS+1)*sizeof(double));
	}
	else {
		Mtotal=(double *)malloc((NUM_MASS+1)*sizeof(double));
	}

	IndexTable=(long*)malloc((MAX_INDEX+5)*sizeof(long));

	dens_veloc_radii = (double *)malloc(3*NUM_MASS*sizeof(double));
	num_within = (double *)malloc(3*NUM_MASS*NUM_MASS*sizeof(double));
	temp_within = (double *)malloc(3*NUM_MASS*NUM_MASS*sizeof(double));
	dens_within = (double *)malloc(3*NUM_MASS*NUM_MASS*sizeof(double));
	velaver = (double *)malloc(NUM_MASS*sizeof(double));
	countmass = (double *)malloc(NUM_MASS*sizeof(double));
	temp_within_core = (double *)malloc(NUM_MASS*sizeof(double));
	dens_within_core = (double *)malloc(NUM_MASS*sizeof(double));
	num_within_core = (double *)malloc(NUM_MASS*sizeof(double));
	whichrad = (long *)malloc(NUM_MASS*sizeof(long));

	
/*======= Reading of values for the mass spectrum =======*/

curr_mass=(char *)strtok(MASS_SPEC,",; ");
sscanf(curr_mass,"%lf",&mass_spec[0]);

for(i=0; (curr_mass=(char *)strtok(NULL," ,;"))!=NULL; i++)  
   sscanf(curr_mass,"%lf",&mass_spec[i+1]); 
    
/*======= Reading of values for the Lagrange radii =======*/

curr_mass=(char *)strtok(MASS_PC,",; ");

for(i=0; (curr_mass=(char *)strtok(NULL," ,;"))!=NULL; i++)  
   sscanf(curr_mass,"%lf",&mass_pc[i]); 


/*======= Reading of values for the mass distribution parameters =======*/

curr_param=(char *)strtok(MASS_DIST_PARAMS,",; ");
sscanf(curr_param,"%lf",&dist_params[0]);

for(i=0; (curr_param=(char *)strtok(NULL," ,;"))!=NULL; i++)  
   sscanf(curr_param,"%lf",&dist_params[i+1]); 

   if(dist_params[0]==1) {
      if((TOTAL_PARAMS-1)!=NUM_MASS) {
         fprintf(stderr,"\n\n ERROR: Mass distribution function of type 1:\n");
         fprintf(stderr," Some mass-classes have not been assigned a probability-value. \n\n");
         return 0;
      }
      normalization=0;      

      for(i=1;i<TOTAL_PARAMS+1;i++){
         normalization=normalization+dist_params[i];
	   }
      difference=1-normalization;
      if(difference<0) difference=0-difference;  
      if(difference>0.00001){
        fprintf(stderr,"\n\n ERROR: Mass distribution function is not properly normalized\n\n");
        return 0; 
      }      
   }

/*======= Other initial values =======*/
   
	pi = acos(-1.0) ;

	TidalMassLoss = 0.0 ;	/* Reset total mass loss due to stellar evolution & tidal radius */
	Etidal = 0.0 ;			/* Immediate total energy of tidally stripped stars */
	StellarMassLoss = 0.0 ;
	StartTime = 0.0 ;		/* Time where the simulation begins... useful for Restarts */

/*======= Opening of output files =======*/

	if(RESTART==1) 
		sscanf("a", "%s", outfilemode);
	else 
		sscanf("w", "%s", outfilemode);

	if(printstatus==0) {
	  
	  for(i=0; i<5;i++){
	    
	    sprintf(outfile,"%s%d", outprefix,i+1);
	    
	    if ((out[i+1] = fopen(outfile, outfilemode))==NULL){
	      printf("\n\nCannot create output file %s\n", outfile);
	      return 0;
	    }
	  }
	}
	
	out[0]=stdout; /* Used to redirect output from files
			 to screen for the case: printstatus={1,2,3,4,5} */    


	sprintf(outfile,"%s0", outprefix);
	if ((logfile = fopen(outfile, outfilemode))==NULL){
	  printf("\n\nCannot create log output file %s\n", outfile);
	  return 0;
	}
	
	/* output files for TRACE stars -- binaries are also output as trace stars */

	if (TRACE_STARS > 0 || N_BINARY > 0) {
		sprintf(outfile,"%strace1", outprefix);
		if ((tracefile1 = fopen(outfile, outfilemode))==NULL){
		  printf("\n\nCannot create trace output file %s\n", outfile);
		  return 0;
		}
		sprintf(outfile,"%strace2", outprefix);
		if ((tracefile2 = fopen(outfile, outfilemode))==NULL){
		  printf("\n\nCannot create trace output file %s\n", outfile);
		  return 0;
		}
		sprintf(outfile,"%strace3", outprefix);
		if ((tracefile3 = fopen(outfile, outfilemode))==NULL){
		  printf("\n\nCannot create trace output file %s\n", outfile);
		  return 0;
		}
	}

	/* File for parameters of escaping stars */
	sprintf(outfile,"%sesc", outprefix);
	if ((escfile = fopen(outfile, outfilemode))==NULL){
	  printf("\n\nCannot create escapers file %s\n", outfile);
	  return 0;
	}

	/* File for parameters of escaping stars */
	sprintf(outfile,"%sstellar", outprefix);
	if ((stellarfile = fopen(outfile, outfilemode))==NULL){
	  printf("\n\nCannot create stellar file %s\n", outfile);
	  return 0;
	}

	if (N_BINARY > 0) {
		/* File for parameters of escaping stars */
		sprintf(outfile,"%sbinary", outprefix);
		if ((binaryfile = fopen(outfile, outfilemode))==NULL){
		  printf("\n\nCannot create binary file %s\n", outfile);
		  return 0;
		}
	}

return (1);

}


/**************** Output of radii, energies, angular momenta, masses, and velocities *********/

void output_rawdump(void)

{

static long i ;
double line_to_store[7];


	for(i=1;i<=N_MAX;i++) {

		line_to_store[0]=TotalTime;
		line_to_store[1]=sr[sindex[i]];
		line_to_store[2]=sm[sindex[i]];
		line_to_store[3]=sE[sindex[i]];
		line_to_store[4]=sJ[sindex[i]];
		line_to_store[5]=svr[sindex[i]];
		line_to_store[6]=svt[sindex[i]];   

		if(printstatus==2) 
			fprintf(out[0],"%.6G %.6G %.6G %.6G %.6G %.6G %.6G \n ",line_to_store[0], \
			line_to_store[1], line_to_store[2], line_to_store[3], line_to_store[4], \
			line_to_store[5], line_to_store[6]);
		else if(printstatus==0) 
			fwrite(&line_to_store,sizeof(double),7,out[2]);
	}

	fprintf(logfile, "\n\nRaw snapshot stored... t = %.8G\n", TotalTime) ;
	fprintf(stderr, "\n\nRaw snapshot stored... t = %.8G\n", TotalTime) ;

/* Dump snapshot for restart also... */

	output_restartdump() ;
}


/************************* Mass-Spectrum Functions ***************************/

double find_nearest(double m_pick, double *list, int NumMembers)
{
static int i;
static double distance, new_distance, nearest;

distance=100;

 for(i=0;i<NumMembers;i++) {
    new_distance = list[i]-m_pick;
    if(new_distance<0) new_distance=(0-new_distance);
    if(new_distance < distance) {
      distance=new_distance;
      nearest=list[i];
    }
 }

 return nearest; 
}

double find_max(double *list, int NumMembers)
{
static int i;
static double max=0;

     for(i=0;i<NumMembers;i++) if(list[i]>max) max=list[i];
     
return max;
}

double find_min(double *list, int NumMembers)
{
static int i;
static double min=50;

     for(i=0;i<NumMembers;i++) if(list[i]<min) min=list[i];
     
return min;
}

int find_position(double item, double *list, int NumMembers)
{
static int i;
for(i=0; i<NumMembers; i++) {
   if(item < list[i]+0.05 && item > list[i]-0.05) return i;
   }
return 0;
}

double mass_spec_dist(double argument, double *parameters)

{
  /* Uses a properly-normalized probability distribution to calculate 
     the discrete probability for the mass-value "argument".  The 
     distribution is specified according to the number (#) that 
     precedes a list of parameters (parem1,parem2...paremN) in the 
     input-file expression:  
 
     MASS_DIST_PARAMS #;parem1,parem2,...,paremN 

     - Distribution 0 is the equal-mass distribution.  The switch for 
       distribution 0 is intercepted in mass_spec_assign.  
 
     - Syntax for the other distributions are given below. */

  double output, b,a, normalization, TotalNormalization, missing;
  long location, switch_value, n, T, min_location;
  long i; 
  double l,q;
  
  a=min; b=max;
  T=0;
  normalization = 0;
  switch_value =  parameters[0];

  switch(switch_value) {
  case 1: /*=== User-Defined Distribution ===*/

          /* User specifies a probability for each value in 
             in the user-defined mass-spectrum (MASS_SPEC).
             Must be normalized to 1.  (Else, an earlier check
             will halt the program.)  Each parameter specifies
             the percentage of the total mass that is made up
             from the corresponding mass species.  In the 
             inputfile the user must use the following syntax: 
             e.g., for a mass-spectrum with three members:

             MASS_DIST_PARAMS 1;0.3,0.2,0.5     

             N.B. It is possible to assign equal masses (of value 1.0)
             to all stars in the following way:

             MASS_SPEC 1.0
             MASS_DIST_PARAMS 1:1        */
 
 
          location=find_position(argument,mass_spec,NUM_MASS);

          min_location=find_position(a,mass_spec, NUM_MASS);    

          for(i=0;i<NUM_MASS;i++) {
             if(i!=min_location) normalization=normalization+(TOTAL_MASS/N_STAR)*parameters[i+1]/mass_spec[i];
             else missing=(TOTAL_MASS/N_STAR)*parameters[i+1]/mass_spec[i];
	  }

          TotalNormalization = normalization+missing;

          if(TotalNormalization>1.005) {
            fprintf(stderr,"\n\nThe user-defined mass distribution is not normalizable.\n\n");
            exit(0);
		    }

          if(argument!=a) output=(TOTAL_MASS/N_STAR)*parameters[location+1]/mass_spec[location];
            else output=1-normalization;  

          break; 

  case 2: /*=== Uniform Distribution ===*/ 

          output=1;
          break;

  case 3: /*=== Descending Step-Function Distribution ===*/
       
          /* The number of steps is given by the first and only
             number in the parameter list.  e.g., the descending
             step-function with four steps is specified by:

             MASS_DIST_PAREMS 3;4 */

          n=parameters[1];
          for(i=0;i<n;i++)T=T+(i+1);
          l=n/(T*(b-a));
          q=(argument-0.001)*n*pow(b,-1);
          for(i=0;i<q;i++);
          output=(n-(i-1))*l;
          break; 

  case 4: /*=== 1/x Distribution ===*/

          output=pow(log(b*pow(a,-1)),-1)*pow(argument,-1);
          break;

  case 5: /*=== Half-sine Distribution ===*/

          output=pow(((a-b)/pi)*(cos((pi*b)/(b-a))-cos((pi*a)/(b-a))),-1)*sin(pi*pow(b-a,-1));
          break;
  }

return output; 
}

long mass_spec_assign(int type) 
{
	static double m, m_pick, y_pick, y, AverageMass;
	static long i,j, count, curr_count, limit;
	static int mass_pos;
          
        /* Assigning mass=0.0 to 0th star and N_STAR+1st star and Mtotal[0]=1.0 */
         
	  sm[0]=0.0;
	  sm[N_STAR+1] = 0.0 ;
	  Mtotal[0] = 1.0 ;
	  actual_total_mass= 0.0;
	  for(i=0;i<NUM_MASS;i++) actual_dist[i]=0.0;

	  max=find_max(mass_spec, NUM_MASS);
	  min=find_min(mass_spec, NUM_MASS);

	  /* An attempt to diminish edge-discrepancy (called EDGE): 
	     num_mass=NUM_MASS; scale_factor=pow(num_mass,-1); */       

	  range=(max-min);  /* (EDGE): +(max-min)*(2*scale_factor); */

	if( type == 0) { 

	  /* The TYPE=0 option assigns mass values (in units of average 
	     mass) according to the discrete distribution assigned in the 
             input file. Uses the Neuman rejection technique. */

	  for (count = 1; count <= N_STAR ; count++) {	
	    
	    for (i=1 ; i <= N_TRY ; i++) {
         
	      if(dist_params[0]==1) {
			m_pick=NUM_MASS*ran2(&idum);
			for(j=0;j<m_pick;j++);
			m=mass_spec[j-1];
	      }
	      else {
			m_pick = range*ran2(&idum)+min;   /*(EDGE): +(min-(max-min)*scale_factor);*/   
			m=find_nearest(m_pick, mass_spec, NUM_MASS);	
	      }
	      
	      y_pick = 1.005*ran2(&idum) ;
	      y = mass_spec_dist(m, dist_params);
	      
	      if (y_pick < y) {	
			sm[count]=m;
			mass_pos=find_position(m,mass_spec,NUM_MASS);
			actual_total_mass=actual_total_mass+m;                  
			actual_dist[mass_pos]++;            
			break ;
	      }
	    }
	    if (i == N_TRY+1) return(-1) ;
	  }
	  
	}

	else if (type == 1) {
	  
/*	     The TYPE=1 option assigns mass-values (in terms of average-mass) to 
	     stars in each component, in order. That is, the first component is 
	     filled first of all (in increasing index of sm[]), and the others
	     follow, in order.  Once radii are assigned and the stars re-indexed
             according to radius, this order is, of course, obliterated.
*/

	  curr_count = 1; 
 
	  for(j=1; j <= NUM_MASS; j++) {
 
        y = mass_spec_dist(mass_spec[j-1], dist_params);  /* Just to test normalizability */

        if(j !=NUM_MASS) 
			limit = curr_count + (dist_params[j]*TOTAL_MASS/mass_spec[j-1]);  
   	    else 
			limit = N_STAR; 	    
	    
	    for(count = curr_count; count <= limit; count++) {
                sm[count] = mass_spec[j-1];
                actual_total_mass+=mass_spec[j-1];
                actual_dist[j-1]++;
		}
		curr_count = count;

	  }
	}

	for(i=0;i<NUM_MASS;i++) actual_dist[i]=(mass_spec[i]*actual_dist[i])/actual_total_mass;
	  
	AverageMass = actual_total_mass/N_STAR;
	for(i=1;i<=N_STAR;i++) {
	  sm[i]=sm[i]/AverageMass;
	}
        
	return(N_STAR);

}




int find_mass_bin(double Mass)

{
 static int i;
 static double step, j;
 static double MAX, MIN;

 MAX = max * (N_STAR / actual_total_mass);
 MIN = min * (N_STAR / actual_total_mass);

 step = (MAX-MIN)/(NUM_MASS_RADII_BINS); 

 j=MIN;  
     if(Mass<MIN+0.001 && Mass> MIN-0.001) i=0;
	  	  else {
                      for(i=0;j<Mass-0.001;i++)j=j+step;
		      i--;
		  }
  return i;
}



/*************************** Cody's Changes ******************************/

/* Uses initial mass power law of the form P(M) ~ M^alpha */

long initial_mass_powerlaw(double alpha)
{
  int i,j, mass_pos;
  double max, min, current_mass, current_value, current_test_output, output_range, AverageMass;


/* Assigning mass=0.0 to 0th star and N_STAR+1st star and Mtotal[0]=1.0 */
     
  sm[0]=0.0;
  sm[N_STAR+1] = 0.0 ;
  Mtotal[0] = 1.0 ;
  actual_total_mass= 0.0;
  for(i=0;i<NUM_MASS;i++) actual_dist[i]=0.0;


  max=find_max(mass_spec, NUM_MASS);
  min=find_min(mass_spec, NUM_MASS);
  output_range= pow(min, alpha);

  if (printstatus==0 || printstatus==6) {
	fprintf(stderr,"mass range: min = %g  max = %g  alpha = %g mass_function range = %g\n", \
	  min,max,alpha,output_range);
	fprintf(logfile,"mass range: min = %g  max = %g  alpha = %g mass_function range = %g\n", \
	  min,max,alpha,output_range);
  }

  for(i=1; i<=N_STAR; i++) {

    for (j=1 ; j <= N_TRY ; j++) {    
      current_mass = ran2(&idum) * (max-min) + min;
      current_value = pow(current_mass, alpha);
      current_test_output=ran2(&idum) * 1.005 * output_range;
      if(current_test_output < current_value) {
		sm[i]=current_mass;
		mass_pos=find_position(current_mass,mass_spec,NUM_MASS);
		actual_total_mass=actual_total_mass+current_mass;                  
		actual_dist[mass_pos]++;            
		break;
      }
    }
    if (j==N_TRY+1) return(-1);
  }

  for(i=0;i<NUM_MASS;i++) actual_dist[i]=(mass_spec[i]*actual_dist[i])/actual_total_mass;
	  
  AverageMass = actual_total_mass/N_STAR;
  for(i=1;i<=N_STAR;i++) {
    sm[i]=sm[i]/AverageMass;
  }
  
  if (printstatus==0 || printstatus==6) {
	fprintf(stderr,"Total mass  = %g   Average Mass = %g  (solar masses)\n", \
		actual_total_mass, AverageMass);
	fprintf(logfile,"Total mass  = %g   Average Mass = %g  (solar masses)\n", \
		actual_total_mass, AverageMass);
  }

  return (0) ;
}



/*======= End of Mass Spectrum Functions =======*/



/****** Starlab Initial Model Fetch With Masses -- Modified 08/03/99 -- KJ *****

  Retrieves starlab-model-assigned masses, radii and velocities, in a file 
  specified in the input file. Correspondes to the MODEL = 1 option.
*/

int starlab_model_fetch_with_masses(void)

{
	  double rad1, rad2, rad3, vel1, vel2, vel3 ;
	  double v_dot_r, angle, mass; 
	  double absol_v, absol_r, AverageMass;
	  int i;

/*	  double cm_vx, cm_vy, cm_vz, cm_x, cm_y, cm_z ;
*/

  static  char heading_name[100];
  FILE *into;

  sr[0]=ZERO;
  sr[N_STAR+1]=INFINITY;
  sm[0]=0.0;
  sm[N_STAR+1] = 0.0 ;

/*======= Opening of input & output files =======*/

	if ((into = fopen(INPUT_SNAPSHOT_FILE, "r"))==NULL){      
	  printf("\nERROR: Cannot open input snapshot file...\n");
	  return 0;
	}

/*===============================================*/

i = 0;

actual_total_mass = 0.0;
/* cm_x = 0.0 ; cm_y = 0.0 ; cm_z = 0.0 ;
   cm_vx = 0.0 ; cm_vy = 0.0 ; cm_vz = 0.0 ;
*/

while(fscanf(into,"%s",heading_name)>0) {

if(strcmp(heading_name,"(Dynamics")==0) {   

  fscanf(into, "%*s %*s %lf\n",&mass); 
  fscanf(into, "%*s %*s %lf %lf %lf\n", &rad1, &rad2, &rad3);
  fscanf(into, "%*s %*s %lf %lf %lf\n", &vel1, &vel2, &vel3);

/* Convert position and velocity to CM coordinates */
/*  ----- Not used, since coordinates are ALREADY relative to CM ---- 
	rad1 -= cm_x ; rad2 -= cm_y ; rad3 -= cm_z ;
	vel1 -= cm_vx ; vel2 -= cm_vy ; vel3 -= cm_vz ;
*/

  absol_v = sqrt( vel1*vel1 + vel2*vel2 + vel3*vel3);
  absol_r = sqrt( rad1*rad1 + rad2*rad2 + rad3*rad3);

  v_dot_r = vel1*rad1 + vel2*rad2 + vel3*rad3;

  angle  = acos ( v_dot_r / (absol_v * absol_r) );

/* Ignore first entry -- gives Total Mass & System center of mass params */

  if(i == 0) {    
/*	cm_x = rad1 ; cm_y = rad2 ; cm_z = rad3 ;
	cm_vx = vel1 ; cm_vy = vel2 ; cm_vz = vel3 ;
*/
  }
  else {
    sm[i] = mass;
    actual_total_mass += mass;
    svt[i] = absol_v * sin(angle);
    svr[i] = absol_v * cos(angle); 

    sr[i] = absol_r;
  }
 
  i++;

}
}
	AverageMass = actual_total_mass/N_STAR;
	for(i=1;i<=N_STAR;i++) {
	  sm[i]=sm[i]/AverageMass;
	}
        
  fclose(into);
  return i-1;

}




/*********** Superposed Starlab "Snapshot" Models Fetched **************

  Starlab-generaged models with mass components in the files starlab#.dat where 
  # = 1, 2, 3,..., must contain exactly the number of stars that is implied for 
  each component by the input file parameters TOTAL_MASS, MASS_SPEC,
  MAST_DIST_PARAMS, and N_STAR.  Assigns velocities and radii from the starlab
  file, and masses from the user-defined input-file mass spectrum parameters.
  Corresponds to the MODEL = 2 option.
*/

int starlab_model_spec_fetch(void)

{
  static double rad1, rad2, rad3;
  static double vel1, vel2, vel3;
  static double v_dot_r, angle; 
  static double absol_v, absol_r;
  static long start;

  static long count, curr_count;
  static int j;

  static char heading_name[100];
  static char file_name[20];
  FILE *into;

  sr[0] = ZERO ;	/* tiny number so that 1/sr[0] is still finite */
  sr[N_STAR+1] = INFINITY ;

  curr_count = 1;
 
  for( j=0; j<NUM_MASS; j++) {

	  /*======= Opening of input & output files =======*/
	  sprintf(file_name,"starlab%d.dat",j+1); 
	  
	  if ((into = fopen(file_name, "r"))==NULL){      
	    printf("\n\nCannot open input file...\n");
	    return 0;
	  }
	    
	  /*===============================================*/             	  


	  start = 1;
	  count = curr_count;

	  while(fscanf(into,"%s",heading_name)>0) {
	      
	    if(strcmp(heading_name,"(Dynamics")==0) {   
		
	      fscanf(into, "%*s %*s %*s\n"); 

	      fscanf(into, "%*s %*s %lf %lf %lf\n", &rad1, &rad2, &rad3);
	      fscanf(into, "%*s %*s %lf %lf %lf\n", &vel1, &vel2, &vel3);
	      
	      absol_v = sqrt( vel1*vel1 + vel2*vel2 + vel3*vel3);
	      absol_r = sqrt( rad1*rad1 + rad2*rad2 + rad3*rad3);
	      
	      v_dot_r = vel1*rad1 + vel2*rad2 + vel3*rad3;

	      angle  = acos ( v_dot_r / (absol_v * absol_r) );

	      if(start != 1) {    /* Ignore first entry (suspicious) in starlab.dat */
		svt[count] = absol_v * sin(angle) / sqrt(mass_spec[j] / (actual_total_mass / N_STAR));
		svr[count] = absol_v * cos(angle) / sqrt(mass_spec[j] / (actual_total_mass / N_STAR));
		sr[count] = absol_r;
                count++;	      
               
	      }
	      start = 0;
	    }    
	  }
	   
	  fclose(into);
	  curr_count = count;
  }
  
  return(N_STAR) ;
}



/**************** Starlab Initial Model Fetch Without Masses **********

  Retrieves starlab model-assigned radii and velocities, in the input snapshot file.
  Masses are assigned according to the user-defined mass-spectrum in the input file. 
  Corresponds to the MODEL = 3 option.
  */

int starlab_model_fetch_without_masses(void)

{

  double rad1, rad2, rad3;
  double vel1, vel2, vel3;
  double v_dot_r, angle; 
  double absol_v, absol_r;
  int i;

  static  char heading_name[100];
  FILE *into;

  sr[0]=ZERO;
  sr[N_STAR+1]=INFINITY;

/*======= Opening of input & output files =======*/

	if ((into = fopen(INPUT_SNAPSHOT_FILE, "r"))==NULL){      
	  printf("\n\nCannot open input file...\n");
	  return 0;
	}

/*===============================================*/

i = 1;

while(fscanf(into,"%s",heading_name)>0) {

if(strcmp(heading_name,"(Dynamics")==0) {   

  fscanf(into, "%*s %*s %*s\n"); 

  fscanf(into, "%*s %*s %lf %lf %lf\n", &rad1, &rad2, &rad3); 
  fscanf(into, "%*s %*s %lf %lf %lf\n", &rad1, &rad2, &rad3);
  fscanf(into, "%*s %*s %lf %lf %lf\n", &vel1, &vel2, &vel3);

  absol_v = sqrt( vel1*vel1 + vel2*vel2 + vel3*vel3);
  absol_r = sqrt( rad1*rad1 + rad2*rad2 + rad3*rad3);

  v_dot_r = vel1*rad1 + vel2*rad2 + vel3*rad3;

  angle  = acos ( v_dot_r / (absol_v * absol_r) );

  if(i != 1) {    /* Ignore first entry (suspicious) in starlab.dat */
    svt[i-1] = absol_v * sin(angle) / sqrt(sm[i-1]);
    svr[i-1] = absol_v * cos(angle) / sqrt(sm[i-1]);

    sr[i-1] = absol_r;
  }

  i++;

}
}
  fclose(into);
  return i-2;

}



double timestep_profile(double x)

{
 double y;

 y = 0.9*pow(1-x,0.1)+0.1;

 return y;
}




#define NRANSI

void splint(double xa[], double ya[], double y2a[], long n, double x, double *y)
{
	void nrerror(char error_text[]);
	long klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void spline(double x[], double y[], long n, double yp1, double ypn, double y2[])
{
	long i,k;
	double p,qn,sig,un,*u;

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}
#undef NRANSI



/******** Read in stellar evolution data -- stellar masses and lifetimes *******/

void read_stellar_data(void)
{
static  double initialmass[100], remnantmass[100], lifetime[100] ;
static  double  y2[100], y3[100], ltime, realmass, remmass;
static  long j, i;
static  FILE *in; 


/*======= Opening of input & output files =======*/

	if ((in = fopen("stellardata.dat", "r"))==NULL){      
	  printf("\n\nCannot open stellardata.dat file...\n");
	  return;
	}

/*======= Reading of parameter values from input file =======*/

	i=1;

	initialmass[0]=0.0; lifetime[0]=100.0; remnantmass[0] = 0.0 ;

	while(fscanf(in,"%lf %lf %lf",&initialmass[i],&lifetime[i],&remnantmass[i])>0) {
	  i++;}

	if (printstatus==0 || printstatus==6) {
		fprintf(stderr,"Stellar evolution data loaded\n") ;
	}

	fclose(in);

	spline(initialmass, lifetime, i-1, 1.0e50, 1.0e50, y2);
	spline(initialmass, remnantmass, i-1, 1.0e50, 1.0e50, y3);

	turnoff_mass = 0.0 ;	/* MAX mass of main sequence star */

	for (j=1; j<=N_STAR; j++) {
	  
		realmass=sm[j]*actual_total_mass/N_STAR;  /* convert mass into SOLAR MASSES */

		/* get lifetime from spline interpolation of stellar data */
		splint(initialmass, lifetime, y2, i-1, sm[j]*actual_total_mass/N_STAR, &ltime) ;
		sLifeTime[j] = pow(10.0,ltime) ;

		/* get remnant mass from spline interpolation of stellar data */
		splint(initialmass, remnantmass, y3, i-1, realmass, &remmass);
		sRemnantMass[j]=remmass;	/* remnant mass in SOLAR MASSES */

		sInitialMass[j] = sm[j] ;	/* save initial mass in CODE units */

		/***** assign TYPE to each star --  initially ALL main-sequence stars.
		1 - Main-sequence 
		2 - White dwarf
		3 - Neutron star
		4 - Black hole
		*/

		stype[j] = 1 ;	/* Main sequence star */

		if (sm[j] > turnoff_mass) turnoff_mass = sm[j] ;	/* keep track of max mass */
	
	}

	/* Set total counts and masses of all mass components */

	count_1 = N_STAR ; count_2 = 0 ; count_3 = 0 ; count_4 = 0 ; count_5 = 0 ;

	/* Masses sm[] are in units of average mass, so total mass 
	   in units of average mass = actual_total_mass / average mass = N_STAR.
	   Here mass_1 .. mass_4 are in units of average mass.
	*/
	mass_1 = N_STAR ; mass_2 = 0.0 ; mass_3 = 0.0 ; mass_4 = 0.0 ; mass_5 = 0.0;

	
	if (printstatus==0 || printstatus==6) {
		fprintf(stderr,"Spline interpolation of stellar data complete\n") ;
	}

}



/********** Compute mass loss due to stellar evolution *******************/

void DoStellarEvolution(void) {

  long i, k ;
  double m_real ;

  /* Do NOT change masses of escaped stars... need original mass for RESTART 
		-- hence only check bound N_MAX stars.
  */

  turnoff_mass = 0.0 ;	/* recalculate turnoff_mass */

  for (k=1 ; k<=N_MAX ; k++) {

	i = sindex[k] ;
    m_real = sm[i] * actual_total_mass / N_STAR ;

	/*  Chernoff & Weinberg 1990 stellar evolution recipe */
    if(STELLAR_EVOLUTION==1) {
		if ((stype[i] == 1) && (RealTime > sLifeTime[i])) {
			if (m_real < 4.7) {
				sm[i] = (0.58 + 0.22 * (m_real - 1.0)) * (N_STAR/actual_total_mass) ;
				count_2++ ;	/* New white dwarf created */
				count_1-- ; /* One MS star lost */
				mass_2 += sm[i] ;	/* Add mass to WD bin */
				mass_1 -= sInitialMass[i] ; /* Remove mass from MS bin */	
				stype[i] = 2 ;
			}
			else if ((m_real >= 4.7) && (m_real < 8.0)) {
				sm[i] = 0.0 ;	/* Star COMPLETELY destroyed according to CW90 */
				count_2++ ;	/* Treat destroyed star as new WD created */
				count_1-- ;	/* One MS star lost */
				mass_2 += sm[i] ;	/* Add mass to WD bin */
				mass_1 -= sInitialMass[i] ; /* Remove mass from MS bin */
				stype[i] = 2 ;
			/* Destroyed star is later removed (as escaped star) in get_positions() */

			}
			else if (m_real <= 15.0) {
				sm[i] = 1.4 * (N_STAR/actual_total_mass) ;
				count_3++ ;	/* New NS created */
				count_1-- ;	/* One MS star lost */
				mass_3 += sm[i] ;	/* Add mass to NS bin */
				mass_1 -= sInitialMass[i] ; /* Remove mass from MS bin */
				stype[i] = 3 ;
			}
			
			Eescaped += (sInitialMass[i]-sm[i]) / N_STAR * sE[i] ;
			Jescaped += (sInitialMass[i]-sm[i]) / N_STAR * sJ[i] ;
			StellarMassLoss += (sInitialMass[i] - sm[i]) / N_STAR ;
		}

	}

	/* Starlab stellar evolution recipe */
	else if(STELLAR_EVOLUTION==2) {
		if ((stype[i] == 1) && (RealTime > sLifeTime[i])) {
			if (m_real > 8.0) {
				sm[i] = 1.4 * (N_STAR/actual_total_mass) ;
				count_3++ ;	/* New NS created */
				count_1-- ;	/* One MS star lost */
				mass_3 += sm[i] ;	/* Add mass to NS bin */
				mass_1 -= sInitialMass[i] ; /* Remove mass from MS bin */
				stype[i] = 3 ;
			}
			else {
				/* convert remnant mass to CODE units */
				sm[i]=sRemnantMass[i] * (N_STAR/actual_total_mass);
				count_2++ ;	/* New white dwarf created */
				count_1-- ; /* One MS star lost */
				mass_2 += sm[i] ;	/* Add mass to WD bin */
				mass_1 -= sInitialMass[i] ; /* Remove mass from MS bin */	
				stype[i] = 2 ;				
			}

			Eescaped += (sInitialMass[i]-sm[i]) / N_STAR * sE[i] ;
			Jescaped += (sInitialMass[i]-sm[i]) / N_STAR * sJ[i] ;
			StellarMassLoss += (sInitialMass[i] - sm[i]) / N_STAR ;
		}
	}

	/* Unknown stellar recipe --- ERROR */
    else {
      fprintf(stderr, "ERROR: Wrong input for STELLAR_EVOLUTION variable\n");
    }

	/* recompute turnoff_mass by checking ONLY remaining MS stars */
	if ((stype[i] == 1) && (sm[i] > turnoff_mass)) turnoff_mass = sm[i] ;

  }

#if defined(DEBUGMODE)
	if(printstatus==0 || printstatus==6) 
		fprintf(stderr,"Stellar evolution complete\n") ;
#endif	

}


/************ Read saved complete system snapshot and RESTART ****************/

void restart()
{
  long i ;
  double line_to_read[29] ;
  char binaryinputfile[100];
  FILE *into;

  if (printstatus==0 || printstatus==6) {
	printf("Restarting... loading snapshot ...\n" );
  }

  sprintf(binaryinputfile,"%srestart", outprefix);

  /*============ Opening Input File ============*/
  if ((into = fopen(binaryinputfile, "r+b"))==NULL){      
    printf("\nERROR: Cannot open binary input file...\n");
    return;
  }

/* first line specifies N_STAR 
*/	
  fread(&N_STAR,sizeof(long),1,into);	


/* read in N_STAR + 1 lines containing star data */

for(i=0;i<=N_STAR+1;i++) {

    fread(&line_to_read, sizeof(double), 29, into);

    TotalTime=line_to_read[0];
    sr[i]=line_to_read[1];
    svr[i]=line_to_read[2];
    svt[i]=line_to_read[3];
    sm[i]=line_to_read[4];
    sLifeTime[i]=line_to_read[5];
    sInitialMass[i]=line_to_read[6];
    sE[i]=line_to_read[7];
    sJ[i]=line_to_read[8];
    sgravity[i]=line_to_read[9];
    sEI[i]=line_to_read[10];
    sJI[i]=line_to_read[11];
    sKE[i]=line_to_read[12];
    sPE[i]=line_to_read[13];
    sMi_ri[i]=line_to_read[14];
    sMi[i]=line_to_read[15];
    srescape[i]=line_to_read[16];
    srnew[i]=line_to_read[17];
    svrnew[i]=line_to_read[18];
    svtnew[i]=line_to_read[19];
    srOld[i]=line_to_read[20];
    sX[i]=line_to_read[21];
    sDE[i]=line_to_read[22];
    sDJ[i]=line_to_read[23];
    sSin2Beta[i]=line_to_read[24];
    sr_peri[i]=line_to_read[25];
    sr_apo[i]=line_to_read[26];
    sindex[i]=line_to_read[27];
    sRemnantMass[i]=line_to_read[28];
}

/* read in IndexTable[] --- MAX_INDEX + 5 entries 
*/  
  fread(&MAX_INDEX,sizeof(long),1,into);
  fread(IndexTable,sizeof(long),MAX_INDEX+5,into);

/* read mass spectrum related ARRAYS -- NUM_MASS specifies size of each array */

  fread(&NUM_MASS,sizeof(long),1,into);
  fread(dens_veloc_radii,sizeof(double),3*NUM_MASS,into);
  fread(num_within,sizeof(double),3*NUM_MASS*NUM_MASS,into);
  fread(temp_within,sizeof(double),3*NUM_MASS*NUM_MASS,into);
  fread(dens_within,sizeof(double),3*NUM_MASS*NUM_MASS,into);
  fread(velaver,sizeof(double),NUM_MASS,into);
  fread(countmass,sizeof(double),NUM_MASS,into);
  fread(temp_within_core,sizeof(double),NUM_MASS,into);
  fread(dens_within_core,sizeof(double),NUM_MASS,into);
  fread(num_within_core,sizeof(double),NUM_MASS,into);
  fread(whichrad,sizeof(long),NUM_MASS,into);
  fread(mass_r,sizeof(double),(NUM_MASS_RADII_BINS+1)*MASS_PC_COUNT,into);

/* MODEL = 4 specifies a CONTINUOUS mass spectrum, which requires
   NUM_MASS_RADII_BINS different components in Mtotal[].
   Otherwise, for discrete mass components, Mtotal[] has NUM_MASS components 
*/
  if(MODEL==4) 
	  fread(Mtotal,sizeof(double),(NUM_MASS_RADII_BINS+1),into);
  else 
	  fread(Mtotal,sizeof(double),(NUM_MASS+1),into);

  fread(actual_dist,sizeof(double),NUM_MASS,into);


/* Read in all other evolution parameters -- Read in TotalTime as StartTime 
*/
  fread(&tcount,sizeof(long),1,into);
  fread(&StartTime,sizeof(double),1,into);	
  fread(&N_MAX,sizeof(long),1,into);
  fread(&Eescaped,sizeof(double),1,into);
  fread(&Jescaped,sizeof(double),1,into);
  fread(&Etotal,sizeof(double),1,into);
  fread(&Etotal_New,sizeof(double),1,into);
  fread(&Etotal_initial,sizeof(double),1,into);
  fread(&KEtotal,sizeof(double),1,into);
  fread(&PEtotal,sizeof(double),1,into);
  fread(&RealTime,sizeof(double),1,into);
  fread(&Rtidal,sizeof(double),1,into);

  fread(&max,sizeof(double),1,into);
  fread(&min,sizeof(double),1,into);
  fread(&range,sizeof(double),1,into);
  fread(&actual_total_mass,sizeof(double),1,into);
  fread(&Dt,sizeof(double),1,into);
  fread(&max_r,sizeof(double),1,into);

  fread(&sin2beta,sizeof(double),1,into);
  fread(&max_sin2beta,sizeof(double),1,into);
  fread(&avg_sin2beta,sizeof(double),1,into);
  fread(&avg_X,sizeof(double),1,into);
  fread(&avg_r,sizeof(double),1,into);
  fread(&avg_dv,sizeof(double),1,into);
  fread(&avg_DE,sizeof(double),1,into);
  fread(&avg_DJ,sizeof(double),1,into);

  fread(&errstat,sizeof(long),1,into);
  fread(&idum,sizeof(long),1,into);
  fread(&n_circular,sizeof(long),1,into);
  fread(&n_pericenter,sizeof(long),1,into);
  fread(&n_apocenter,sizeof(long),1,into);

  fread(&n_orbit,sizeof(long),1,into);
  fread(&n_extreme,sizeof(long),1,into);

  fread(&DEavg,sizeof(double),1,into);
  fread(&DEmax,sizeof(double),1,into);
  fread(&DEdiff,sizeof(double),1,into);
  fread(&Sin2Beta,sizeof(double),1,into);
  fread(&Trc,sizeof(double),1,into);
  fread(&rho_core,sizeof(double),1,into);
  fread(&v_core,sizeof(double),1,into);
  fread(&r_core,sizeof(double),1,into);
  fread(&StellarMassLoss,sizeof(double),1,into);
  fread(&TidalMassLoss,sizeof(double),1,into);

  fclose(into);

  if(printstatus==0 || printstatus==6) {	
	fprintf(stderr,"Time = %.8G   Tcount = %ld\n", TotalTime, tcount) ;
	fprintf(stderr,"N = %ld, Total E = %.8G, Total Mass = %.8G, Virial ratio = %.8G\n", \
			N_MAX, Etotal, Mtotal[0], -2.0*KEtotal/PEtotal);
	fprintf(stderr,"Total KE = %.8G, Total PE = %.8G\n", KEtotal, PEtotal);
	fprintf(stderr,"Snapshot loaded.\n") ;
  }

}


/**** Output complete system snapshot (all global variables) to allow restarting ***/

void output_restartdump(void)
{

  long i ;
  double line_to_store[29];
  FILE *restartout;
  char restartfile[100];

  sprintf(restartfile,"%srestart", outprefix);

  if ((restartout = fopen(restartfile, "w+b"))==NULL){      
    printf("\nERROR: Cannot open restart output file...\n");
    return;
  }


  fwrite(&N_STAR,sizeof(long),1,restartout);

  for(i=0;i<=N_STAR+1;i++) {

    line_to_store[0]=TotalTime;
    line_to_store[1]=sr[i];
    line_to_store[2]=svr[i];
    line_to_store[3]=svt[i];
    line_to_store[4]=sm[i];
    line_to_store[5]=sLifeTime[i];
    line_to_store[6]=sInitialMass[i];
    line_to_store[7]=sE[i];
    line_to_store[8]=sJ[i];
    line_to_store[9]=sgravity[i];
    line_to_store[10]=sEI[i];
    line_to_store[11]=sJI[i];
    line_to_store[12]=sKE[i];
    line_to_store[13]=sPE[i];
    line_to_store[14]=sMi_ri[i];
    line_to_store[15]=sMi[i];
    line_to_store[16]=srescape[i];
    line_to_store[17]=srnew[i];
    line_to_store[18]=svrnew[i];
    line_to_store[19]=svtnew[i];
    line_to_store[20]=srOld[i];
    line_to_store[21]=sX[i];
    line_to_store[22]=sDE[i];
    line_to_store[23]=sDJ[i];
    line_to_store[24]=sSin2Beta[i];
    line_to_store[25]=sr_peri[i];
    line_to_store[26]=sr_apo[i];
    line_to_store[27]=sindex[i];
    line_to_store[28]=sRemnantMass[i];

    fwrite(line_to_store,sizeof(double),29,restartout);
  }


  fwrite(&MAX_INDEX,sizeof(long),1,restartout);
  fwrite(IndexTable,sizeof(long),MAX_INDEX+5,restartout);


  fwrite(&NUM_MASS,sizeof(long),1,restartout);
  fwrite(dens_veloc_radii,sizeof(double),3*NUM_MASS,restartout);
  fwrite(num_within,sizeof(double),3*NUM_MASS*NUM_MASS,restartout);
  fwrite(temp_within,sizeof(double),3*NUM_MASS*NUM_MASS,restartout);
  fwrite(dens_within,sizeof(double),3*NUM_MASS*NUM_MASS,restartout);
  fwrite(velaver,sizeof(double),NUM_MASS,restartout);
  fwrite(countmass,sizeof(double),NUM_MASS,restartout);
  fwrite(temp_within_core,sizeof(double),NUM_MASS,restartout);
  fwrite(dens_within_core,sizeof(double),NUM_MASS,restartout);
  fwrite(num_within_core,sizeof(double),NUM_MASS,restartout);
  fwrite(whichrad,sizeof(long),NUM_MASS,restartout);
  fwrite(mass_r,sizeof(double),(NUM_MASS_RADII_BINS+1)*MASS_PC_COUNT,restartout);

  if(MODEL==4) 
	  fwrite(Mtotal,sizeof(double),(NUM_MASS_RADII_BINS+1),restartout);
  else 
	  fwrite(Mtotal,sizeof(double),(NUM_MASS+1),restartout);

  fwrite(actual_dist,sizeof(double),NUM_MASS,restartout);

 
  fwrite(&tcount,sizeof(long),1,restartout);
  fwrite(&TotalTime,sizeof(double),1,restartout);
  fwrite(&N_MAX,sizeof(long),1,restartout);
  fwrite(&Eescaped,sizeof(double),1,restartout);
  fwrite(&Jescaped,sizeof(double),1,restartout);
  fwrite(&Etotal,sizeof(double),1,restartout);
  fwrite(&Etotal_New,sizeof(double),1,restartout);
  fwrite(&Etotal_initial,sizeof(double),1,restartout);
  fwrite(&KEtotal,sizeof(double),1,restartout);
  fwrite(&PEtotal,sizeof(double),1,restartout);
  fwrite(&RealTime,sizeof(double),1,restartout);
  fwrite(&Rtidal,sizeof(double),1,restartout);

  fwrite(&max,sizeof(double),1,restartout);
  fwrite(&min,sizeof(double),1,restartout);
  fwrite(&range,sizeof(double),1,restartout);
  fwrite(&actual_total_mass,sizeof(double),1,restartout);
  fwrite(&Dt,sizeof(double),1,restartout);
  fwrite(&max_r,sizeof(double),1,restartout);

  fwrite(&sin2beta,sizeof(double),1,restartout);
  fwrite(&max_sin2beta,sizeof(double),1,restartout);
  fwrite(&avg_sin2beta,sizeof(double),1,restartout);
  fwrite(&avg_X,sizeof(double),1,restartout);
  fwrite(&avg_r,sizeof(double),1,restartout);
  fwrite(&avg_dv,sizeof(double),1,restartout);
  fwrite(&avg_DE,sizeof(double),1,restartout);
  fwrite(&avg_DJ,sizeof(double),1,restartout);

  fwrite(&errstat,sizeof(long),1,restartout);
  fwrite(&idum,sizeof(long),1,restartout);
  fwrite(&n_circular,sizeof(long),1,restartout);
  fwrite(&n_pericenter,sizeof(long),1,restartout);
  fwrite(&n_apocenter,sizeof(long),1,restartout);
  fwrite(&n_orbit,sizeof(long),1,restartout);
  fwrite(&n_extreme,sizeof(long),1,restartout);

  fwrite(&DEavg,sizeof(double),1,restartout);
  fwrite(&DEmax,sizeof(double),1,restartout);
  fwrite(&DEdiff,sizeof(double),1,restartout);
  fwrite(&Sin2Beta,sizeof(double),1,restartout);
  fwrite(&Trc,sizeof(double),1,restartout);
  fwrite(&rho_core,sizeof(double),1,restartout);
  fwrite(&v_core,sizeof(double),1,restartout);
  fwrite(&r_core,sizeof(double),1,restartout);
  fwrite(&StellarMassLoss,sizeof(double),1,restartout);
  fwrite(&TidalMassLoss,sizeof(double),1,restartout);

  fclose(restartout);

  if (printstatus==0 || printstatus==6) {
	printf("\nRestart snapshot stored - TotalTime = %g  Tcount = %ld \n", \
		TotalTime, tcount) ;
  }

}



