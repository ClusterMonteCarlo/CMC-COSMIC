/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

/* a fast square function */
inline double sqr(double x)
{
        return(x*x);
}

/* a fast cube function */
inline double cub(double x)
{
        return(x*x*x);
}

/* toggle debugging */
void toggle_debugging(int signal)
{
	if (debug) {
		fprintf(stderr, "toggle_debugging(): turning debugging off on signal %d\n", signal);
		debug = 0;
	} else {
		fprintf(stderr, "toggle_debugging(): turning debugging on on signal %d\n", signal);
		debug = 1;
	}
}

/* close buffers, then exit */
void exit_cleanly(int signal)
{
	close_buffers();
	free_arrays();

	exit(signal);
}

void free_arrays(void){
	free(mass_pc); free(densities_r); free(no_star_r); 
	free(ave_mass_r); free(mass_r);
	free(star); free(binary);
	free(IndexTable);
}

/* GSL error handler */
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno)
{
	fprintf(stderr, "gsl: %s:%d: ERROR: %s\n", file, line, reason);
	exit_cleanly(gsl_errno);
}

void setup_sub_time_step(void){
	long p, i, k, si_minus_p, si_plus_p, sub_imin, si, zk;
	double Ai, Dt_local, m_avg, w2_avg, zr_min, zr_max;
	
	/********** setup sub time step *************************/
	/* Timestep obtained from GetTimeStep is based on the 
	 * relaxation time in the core. Now look for a suitable 
	 * boundary between core and halo. */

	if (sub.count == 0) {		/** do a full timestep **/
		p = 20;
		sub.N_MAX = clus.N_MAX;	/* default -- if Dt does not 
					   change much up to r_h */
		sub.FACTOR = 1;
		sub_imin = 2000;

		for (si = sub_imin; si <= clus.N_MAX; si += clus.N_MAX / 100) {
			si_minus_p = si - p;
			si_plus_p = si + p + 1;

			if (si_minus_p < 0) {
				si_minus_p = 0;
				si_plus_p = 2 * p + 1;
			} else if (si_plus_p > clus.N_MAX) {
				si_plus_p = clus.N_MAX;
				si_minus_p = clus.N_MAX - 2 * p - 1;
			}

			/* calculate Ai for zone */
			m_avg = 0.0;
			w2_avg = 0.0;
			zk = 2 * p + 2;
			for (i = si_minus_p; i <= si_plus_p; i++) {
				k = i;
				m_avg += star[k].m;
				w2_avg += star[k].m * (star[k].vr * star[k].vr + star[k].vt * star[k].vt);
			}
			m_avg = m_avg / zk;
			w2_avg = w2_avg * 2.0 / m_avg / zk;
			zr_min = star[si_minus_p].r;
			zr_max = star[si_plus_p].r;
			Ai = 6.0 * zk * m_avg * m_avg /
			    (zr_max * zr_max * zr_max - zr_min * zr_min * zr_min) / pow(w2_avg, 1.5);

			Dt_local = SIN2BETA_MAX / Ai * clus.N_STAR;
			Dt_local /= DT_FACTOR;
			/* DEBUG */
			Dt_local *= Mtotal;
			/* DEBUG */

			/* DEBUG */
			/* do away with sub timestep business for now */
			/* break; */
			/* DEBUG */

//			if (Dt_local >= Dt * 500 && sub.FACTOR < 500) {
//				if (si * 1.0 / clus.N_MAX <
//				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 500.0) {
//					sub.FACTOR = 500;
//					sub.N_MAX = si;
//					break;
//				}
//			}
//			if (Dt_local >= Dt * 100 && sub.FACTOR < 100) {
//				if (si * 1.0 / clus.N_MAX <
//				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 100.0) {
//					sub.FACTOR = 100;
//					sub.N_MAX = si;
//				}
//			} else if (Dt_local >= Dt * 50 && sub.FACTOR < 50) {
//				if (si * 1.0 / clus.N_MAX <
//				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 50.0) {
//					sub.FACTOR = 50;
//					sub.N_MAX = si;
//				}
//			} else 
			if (Dt_local >= Dt * 25 && sub.FACTOR < 25) {
				if (si * 1.0 / clus.N_MAX <
				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 25.0) {
					sub.FACTOR = 25;
					sub.N_MAX = si;
				}
			} else 
			if (Dt_local >= Dt * 10 && sub.FACTOR < 10) {
				if (si * 1.0 / clus.N_MAX <
				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 10.0) {
					sub.FACTOR = 10;
					sub.N_MAX = si;
				}
			} else 
			if (Dt_local >= Dt * 5 && sub.FACTOR < 5) {
				if (si * 1.0 / clus.N_MAX <
				    1.0 / sub.FACTOR * (sub.N_MAX * 1.0 / clus.N_MAX * sub.FACTOR + 1.0) - 1.0 / 5.0) {
					sub.FACTOR = 5;
					sub.N_MAX = si;
				}
			} else 
			if (Dt_local >= Dt * 2 && sub.FACTOR < 2 && si < clus.N_MAX / 2) {
				sub.FACTOR = 2;
				sub.N_MAX = si;
			} else if (si > clus.N_MAX / 2) {
				break;
			}
		}

		sub.rmax = star[sub.N_MAX].r;
		sub.count = 1;
		/* DEBUG */
		/* I'm pretty sure this should be 0 to turn off the sub timestep stuff */
		/* sub.count = 0; */
		/* DEBUG */
		sub.totaltime = Dt;

	} else { /** core timestep only **/
		sub.count++;
		sub.totaltime += Dt;

		if (sub.count == sub.FACTOR) {	/* last round, so do a FULL time step. */
			sub.count = 0;
		}
	} /******* end setup sub time step **************/

}

void RecomputeEnergy(void) {
	double dtemp;
	long k, i;

	/* Recalculating Energies */
	Etotal.tot = 0.0;
	Etotal.K = 0.0;
	Etotal.P = 0.0;
	Etotal.Eint = 0.0;
	dtemp = 0;

	if (E_CONS == 0) { /* recompute new sE[] and sJ[], using the new potential */
		for (i = 1; i <= clus.N_MAX; i++) {
			k = i;
			star[k].E = star[k].phi + 0.5 * (SQR(star[k].vr) + SQR(star[k].vt));
			star[k].J = star[k].r * star[k].vt;

			Etotal.K += 0.5 * (SQR(star[k].vr) + SQR(star[k].vt)) * star[k].m / clus.N_STAR;

			/* Compute PE using Henon method using star[].phi */
			Etotal.P += star[k].phi * star[k].m / clus.N_STAR;

			/* add up internal energies */
			Etotal.Eint += star[k].Eint;

			/* reset star[].interacted flag to 0 */
			star[k].interacted = 0;
		}
	} else { /* Try to conserve energy by using intermediate potential */
		for (i = 1; i <= clus.N_MAX; i++) {
			k = i;
			
			/* Note: svt[] = J/r_new is already computed in get_positions() */
			/* ignore stars near pericenter, and those with strong interactions */
			if (star[k].X > 0.05 && star[k].interacted == 1) {
				dtemp = star[k].EI - star[k].phi + potential(star[k].rOld);

				if (dtemp - star[k].vr * star[k].vr > 0) {
					/* preserve star[k].vr and change star[k].vt */
					star[k].vt = sqrt(dtemp - star[k].vr * star[k].vr);
				} else {
					if (dtemp > 0) {
						if (dtemp > star[k].vt * star[k].vt) {
							star[k].vr = sqrt(dtemp - star[k].vt * star[k].vt);
						} else {
							star[k].vt = sqrt(dtemp / 2.0);
							star[k].vr = sqrt(dtemp / 2.0);
						}
					} else {
						/* reduce the energy of the next star to compensate */ 
						if (i < clus.N_MAX)
							star[i + 1].EI += dtemp - (star[k].vt * star[k].vt + star[k].vr * star[k].vr);
					}
				}
			}

			/* recompute new sE[] and sJ[], using the new potential */
			star[k].E = star[k].phi + 0.5 * (star[k].vr * star[k].vr + star[k].vt * star[k].vt);
			star[k].J = star[k].r * star[k].vt;

			Etotal.K += 0.5 * (star[k].vr * star[k].vr + star[k].vt * star[k].vt) * star[k].m / clus.N_STAR;
			Etotal.P += star[k].phi * star[k].m / clus.N_STAR;

			/* add up internal energies */
			Etotal.Eint += star[k].Eint;

			/* reset star[].interacted flag to 0 */
			star[k].interacted = 0;
		}
	}

	Etotal.P *= 0.5;
	Etotal.tot = Etotal.K + Etotal.P + Etotal.Eint + cenma.E/clus.N_STAR;
}

/* computes intermediate energies, and transfers "new" dynamical params to the standard variables */
void ComputeIntermediateEnergy(void)
{
	long j;

	/* compute intermediate energies for stars due to change in pot */ 
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
		/* but do only for NON-Escaped stars */
		if (star[j].rnew < 1.0e6) {
			star[j].EI = sqr(star[j].vr) + sqr(star[j].vt) + star[j].phi - potential(star[j].rnew);
		}
	}
	
	/* Transferring new positions to .r, .vr, and .vt from .rnew, .vrnew, and .vtnew */
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
		star[j].rOld = star[j].r;
		star[j].r = star[j].rnew;
		star[j].vr = star[j].vrnew;
		star[j].vt = star[j].vtnew;
	}
}

long CheckStop(void) {

	if (tcount >= T_MAX_COUNT) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("No. of timesteps > T_MAX_COUNT ... Terminating.\n");
		return (1);
	}

	if (TotalTime >= T_MAX) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("TotalTime > T_MAX ... Terminating.\n");
		return (1);
	}

	/* Stop if cluster is disrupted -- N_MAX is too small */
	/* if (clus.N_MAX < (0.02 * clus.N_STAR)) { */
	if (clus.N_MAX < (0.005 * clus.N_STAR)) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("N_MAX < 0.005 * N_STAR ... Terminating.\n");
		return (1);
	}

	/* Stop if Etotal > 0 */
	if (Etotal.K + Etotal.P > 0.0) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("Etotal > 0 ... Terminating.\n");
		return (1);
	}


	/* If inner-most Lagrangian radius is too small, then stop: */
	if (mass_r[0] < MIN_LAGRANGIAN_RADIUS) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("Min Lagrange radius < %.6G ... Terminating.\n", MIN_LAGRANGIAN_RADIUS);
		return (1);
	}

	/* Output some snapshots near core collapse 
	 * (if core density is high enough) */
	if (DUMPS==1){
		if (rho_core > 50.0 && Echeck == 0) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 1.0e2 && Echeck == 1) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 5.0e2 && Echeck == 2) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 1.0e3 && Echeck == 3) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 5.0e3 && Echeck == 4) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 1.0e4 && Echeck == 5) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 5.0e4 && Echeck == 6) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 1.0e5 && Echeck == 7) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 5.0e5 && Echeck == 8) {
			print_2Dsnapshot();
			Echeck++;
		} else if (rho_core > 1.0e6 && Echeck == 9) {
			print_2Dsnapshot();
			Echeck++;
		}

		/* added by ato 
		 * to try to take snapshots for core bounce as well. 
		 * idea is if we reduced core density by 10 percent the
		 * last time we took snapshot, take another one and adjust
		 * parameters to take further snapshots if further collapse
		 * occurs */
		if (rho_core < 0.9e6 && Echeck == 10){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9*5e5 && Echeck == 9){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9e5 && Echeck == 8){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9*5e4 && Echeck == 7){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9e4 && Echeck == 6){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9*5e3 && Echeck == 5){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9e3 && Echeck == 4){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9*5e2 && Echeck == 3){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9e2 && Echeck == 2){
			print_2Dsnapshot();
			Echeck--;
		} else if (rho_core < 0.9*50 && Echeck == 1){
			print_2Dsnapshot();
			Echeck--;
		} 
	}

	/* If total Energy has diminished by TERMINAL_ENERGY_DISPLACEMENT, then stop */
	if (Etotal.tot < Etotal.ini - TERMINAL_ENERGY_DISPLACEMENT) {
		if (DUMPS == 1)
			print_2Dsnapshot();
		dprintf("Terminal Energy reached... Terminating.\n");
		return (1);
	}
	return (0); /* NOT stopping time yet */
}

/* energy calculation function that is called for a restart (so it's not necessary to 
   re-set all global energy variables */
void ComputeEnergy2(void)
{
	long k, i;

	for (i = 1; i <= clus.N_MAX; i++) {
		k = i;
		star[k].E = star[k].phi + 0.5 * (star[k].vr * star[k].vr + star[k].vt * star[k].vt);

		star[k].J = star[k].r * star[k].vt;
	}
	
	fprintf(stdout, "Time = %.8G   Tcount = %ld\n", TotalTime, tcount);
	fprintf(stdout, "N = %ld, Total E = %.8G, Total Mass = %.8G, Virial ratio = %.8G\n",
		clus.N_MAX, Etotal.tot, Mtotal, -2.0 * Etotal.K / Etotal.P);
	fprintf(stdout, "Total KE = %.8G, Total PE = %.8G\n", Etotal.K, Etotal.P);
}

void ComputeEnergy(void)
{
	long k, i;

	Etotal.tot = 0.0;
	Etotal.K = 0.0;
	Etotal.P = 0.0;
	Etotal.Eint = 0.0;

	star[0].E = star[0].J = 0.0;
	for (i = 1; i <= clus.N_MAX; i++) {
		k = i;
		star[k].E = star[k].phi + 0.5 * (star[k].vr * star[k].vr + star[k].vt * star[k].vt);

		star[k].J = star[k].r * star[k].vt;

		Etotal.K += 0.5 * (star[k].vr * star[k].vr + star[k].vt * star[k].vt) * star[k].m / clus.N_STAR;

		Etotal.P += star[k].phi * star[k].m / clus.N_STAR;

		Etotal.Eint += star[k].Eint;
	}
	star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;
	
	Etotal.P *= 0.5;
	Etotal.tot = Etotal.K + Etotal.P + Etotal.Eint + cenma.E/clus.N_STAR;

	fprintf(stdout, "Time = %.8G   Tcount = %ld\n", TotalTime, tcount);
	fprintf(stdout, "N = %ld, Total E = %.8G, Total Mass = %.8G, Virial ratio = %.8G\n",
		clus.N_MAX, Etotal.tot, Mtotal, -2.0 * Etotal.K / Etotal.P);
	fprintf(stdout, "Total KE = %.8G, Total PE = %.8G\n", Etotal.K, Etotal.P);
}

/* Computing the potential at each star sorted by increasing 
   radius. Units: G = 1  and  Mass is in units of total INITIAL mass.
   Total mass is computed by SUMMING over all stars that have NOT ESCAPED 
   i.e., over all stars upto N_MAX <= N_STAR. N_MAX is computed in this 
   routine by counting all stars with radius < SF_INFINITY and Radius of the 
   (N_MAX+1)th star is set to infinity i.e., star[N_MAX+1].r = 
   SF_INFINITY. Also setting star[N_MAX+1].phi = 0. Assuming 
   star[0].r = 0. star[].phi is also indexed i.e. it
   is the value of the potential at radius star[k].r 
   NOTE: Assming here that NO two stars are at the SAME RADIUS upto 
   double precision. Returns N_MAX. Potential given in star[].phi
*/
long potential_calculate(void) {
	long k, ii;
	double mprev, rtemp;

	mprev = 0.0;

	ii = 0;
	for (k = 1; k <= clus.N_STAR_NEW; k++) { /* Recompute total Mass and N_MAX */
		rtemp = star[k].r;
		if (rtemp < 1000000 - 1) {
			mprev = mprev + star[k].m;
			if(isnan(mprev)){
				eprintf("NaN (2) detected\n");
				exit_cleanly(-1);
			}
		} else {
			break;

		}
		while (ii <= MAX_INDEX && rtemp >= (double) ii / INDEX_UNIT) {
			IndexTable[ii] = k;
			ii++;
		}
	}
	clus.N_MAX = k - 1;		/* New N_MAX */
	/* New total Mass; This IS correct for multiple components */
	Mtotal = mprev/clus.N_STAR + cenma.m/clus.N_STAR;	

	for (k = ii; k <= MAX_INDEX; k++)
		IndexTable[k] = clus.N_MAX + 1;

	/* Compute new tidal radius using new Mtotal */

	Rtidal = orbit_r * pow(Mtotal, 1.0 / 3.0);

	star[clus.N_MAX + 1].r = SF_INFINITY;
	star[clus.N_MAX + 1].phi = 0.0;
	mprev = Mtotal;
	for (k = clus.N_MAX; k >= 1; k--) {/* Recompute potential at each r */
		star[k].phi = star[k + 1].phi 
			- mprev * (1.0 / star[k].r 
					- 1.0 / star[k + 1].r);
		mprev -= star[k].m / clus.N_STAR;
	}

	for (k = 1; k <= clus.N_MAX; k++){
		star[k].phi -= cenma.m / clus.N_STAR
			/ star[k].r;
		if(isnan(star[k].phi)){
			eprintf("NaN detected\n");
			exit_cleanly(-1);
		}
	}
	star[0].phi = star[1].phi; /* U(r=0) is U_1 */

	return (clus.N_MAX);
}

void comp_mass_percent(){
	double mprev;
	long int k, mcount;

	/* Computing radii containing mass_pc[i] % of the mass */
	mprev = cenma.m/clus.N_STAR;
	for(mcount=0; mcount<MASS_PC_COUNT; mcount++){
		if ( mprev/Mtotal > mass_pc[mcount] ) {
			mass_r[mcount] = MINIMUM_R;
			ave_mass_r[mcount] = 0.0;
			no_star_r[mcount] = 0;
			densities_r[mcount] = 0.0;
		} else {
			break;
		}
	}
	for (k = 1; k <= clus.N_MAX; k++) {	/* Only need to count up to N_MAX */
		mprev += star[k].m / clus.N_STAR;
		if (mprev / Mtotal > mass_pc[mcount]) {
			mass_r[mcount] = star[k].r;
			ave_mass_r[mcount] = mprev/Mtotal/k*initial_total_mass;
			no_star_r[mcount] = k;
			densities_r[mcount] = mprev*clus.N_STAR/
				(4/3*3.1416*pow(star[k].r,3));
			mcount++;
			if (mcount == MASS_PC_COUNT - 1)
				break;
		}
	}
}

/* The potential computed using the star[].phi computed at the star 
   locations in star[].r sorted by increasing r. */
double potential(double r)
{
	long i, k;
	double henon;

	/* root finding using indexed values of sr[] & bisection */
	if (r < star[1].r)
		return (star[1].phi);

	if (r * INDEX_UNIT <= MAX_INDEX) {	/* use tabulated value of r */
		i = r * INDEX_UNIT;
		i = IndexTable[i];
	} else if (IndexTable[MAX_INDEX] == clus.N_MAX + 1) {	/* ALL stars are WITHIN MAX_INDEX*1000 */
		i = clus.N_MAX;
	} else {		/* beyond tabulated values of r */
		i = FindZero_r(IndexTable[MAX_INDEX] - 1, clus.N_MAX + 1, r);
		if (i < 0) {
			eprintf("Error finding zero: i = %5ld   r = %.5G\n", i, r);
			exit_cleanly(-1);
		}
	}
	while (star[i].r > r) {
		i--;
	}
	while (star[i].r < r) {
		i++;
	}
	k = i - 1;


	/* Henon's method of computing the potential using star[].phi */ 
	if (k == 0)
		henon = (star[1].phi);
	else
		henon = (star[k].phi + (star[k + 1].phi
		       - star[k].phi) 
			* (1.0/star[k].r - 1.0/r) /
			 (1.0/star[k].r - 1.0/star[k + 1].r));
	
	return (henon);
}

/*****************************************/
/* Unmodified Numerical Recipes Routines */
/*****************************************/
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit_cleanly(1);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

#undef NR_END
#undef FREE_ARG

/* update some important global variables */
void update_vars(void)
{
	long i, j, k;
	
	/* update total number, mass, and binding energy of binaries in cluster */
	N_b = 0;
	M_b = 0.0;
	E_b = 0.0;
	for (i=1; i<=clus.N_MAX; i++) {
		j = i;
		k = star[j].binind;
		if (k != 0) {
			N_b++;
			M_b += star[j].m;
			E_b += (binary[k].m1/clus.N_STAR) * (binary[k].m2/clus.N_STAR) / (2.0 * binary[k].a);
		}
	}
	
	/* ejected binaries */
	Ebescaped = 0.0;
	for (i=clus.N_MAX+1; i<=clus.N_STAR; i++) {
		k = star[i].binind;
		if (k != 0) {
			Ebescaped += (binary[k].m1/clus.N_STAR) * (binary[k].m2/clus.N_STAR) / (2.0 * binary[k].a);
		}
	}
}

void mini_sshot(){
	FILE *mss;
	char *mss_fname;
	int fname_len;
	int i;

	fname_len = strlen(outprefix);
	fname_len += strlen("miniss.");
	fname_len += 10;
	mss_fname = malloc(fname_len*sizeof(char));
	sprintf(mss_fname,"%sminiss.%05ld", outprefix, tcount);
	mss = fopen(mss_fname, "w+");
	for(i=0; i<1000; i++){
		fprintf(mss, "%8ld %.16e %.16e %.16e %.16e ", 
				star[i].id, star[i].m, star[i].r_peri, 
				star[i].r, star[i].r_apo);
		fprintf(mss, "%.16e %.16e %.16e %.16e\n", 
				star[i].vr, star[i].vt, star[i].E, star[i].phi);
	}
	fclose(mss);
}

/* set the units */
void units_set(void)
{
	/* define (N-body) units (in CGS here): U_l = U_t^(2/3) G^(1/3) U_m^(1/3) */
	units.t = log(GAMMA * clus.N_STAR)/(clus.N_STAR * MEGA_YEAR) * 1.0e6 * YEAR;
	units.m = clus.N_STAR * initial_total_mass / SOLAR_MASS_DYN * MSUN;
	units.l = pow(units.t, 2.0/3.0) * pow(G, 1.0/3.0) * pow(units.m, 1.0/3.0);
}
