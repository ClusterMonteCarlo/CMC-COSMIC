/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "cmc.h"
#include "cmc_vars.h"

/* a fast square function */
double sqr(double x)
{
        return(x*x);
}

/* a fast cube function */
double cub(double x)
{
        return(x*x*x);
}

/* The potential computed using the star[].phi computed at the star 
   locations in star[].r sorted by increasing r. */
double potential_serial(double r) {
	long i, kmax, kmin;
	double henon;
	struct Interval star_interval;

	/* root finding using indexed values of sr[] & bisection */
	if (r== 0.0) {
		return (star[0].phi);
	};

	if (r < star[1].r) {
	     return (star[0].phi-(star[0].phi-star[1].phi)*star[1].r/r);
	};
   //i = check_if_r_around_last_index(last_index, r);
   i=-1;
   if (i== -1) {
           if (SEARCH_GRID) {
             star_interval= search_grid_get_interval(r_grid, r);
             kmax= star_interval.max;
             kmin= star_interval.min;
           } else {
             kmax= clus.N_MAX+1;
             kmin= 1;
           };
           if (kmin==kmax-1) {
             i= kmin;
           } else {
             i =  FindZero_r_serial(kmin, kmax, r);
           };
	   last_index= i;
   };

	if(star[i].r > r || star[i+1].r < r){
		eprintf("binary search (FindZero_r) failed!!\n");
		eprintf("pars: i=%ld, star[i].r = %e, star[i+1].r = %e, star[i+2].r = %e, star[i+3].r = %e, r = %e\n",
				i, star[i].r, star[i+1].r, star[i+2].r, star[i+3].r, r);
		eprintf("pars: star[i].m=%g star[i+1].m=%g star[i+2].m=%g star[i+3].m=%g\n",
			star[i].m, star[i+1].m, star[i+2].m, star[i+3].m);
		exit_cleanly(-2);
	}

	/* Henon's method of computing the potential using star[].phi */ 
	if (i == 0){ /* I think this is impossible, due to early return earlier,
			    but I am keeping it. -- ato 23:17,  3 Jan 2005 (UTC) */
		henon = (star[1].phi);
	} else {
		henon = (star[i].phi + (star[i + 1].phi - star[i].phi) 
			 * (1.0/star[i].r - 1.0/r) /
			 (1.0/star[i].r - 1.0/star[i + 1].r));
	}
	
	return (henon);
}

/* The potential computed using the star[].phi computed at the star 
   locations in star[].r sorted by increasing r. */
double potential(double r) {
	long i, kmax, kmin;
	double henon;
	struct Interval star_interval;

	/* root finding using indexed values of sr[] & bisection */
	if (r== 0.0) {
#ifdef USE_MPI
		return (star_phi[0]);
#else
		return (star[0].phi);
#endif
	};

#ifdef USE_MPI
	if (r < star_r[1]) {
	     return (star_phi[0]-(star_phi[0]-star_phi[1])*star_r[1]/r);
	};

#else
	if (r < star[1].r) {
	     return (star[0].phi-(star[0].phi-star[1].phi)*star[1].r/r);
	};
#endif
   //i = check_if_r_around_last_index(last_index, r);
   i=-1;
   if (i== -1) {
           if (SEARCH_GRID) {
             star_interval= search_grid_get_interval(r_grid, r);
             kmax= star_interval.max;
             kmin= star_interval.min;
           } else {
             kmax= clus.N_MAX+1;
             kmin= 1;
           };
           if (kmin==kmax-1) {
             i= kmin;
           } else {
             i =  FindZero_r(kmin, kmax, r);
           };
	   last_index= i;
   };

#ifdef USE_MPI
	if(star_r[i] > r || star_r[i+1] < r){
#else
	if(star[i].r > r || star[i+1].r < r){
#endif

		//MPI2: Need to be changed to global arrays later. Ignoring for now since this is only for debugging.
		eprintf("binary search (FindZero_r) failed!!\n");
		eprintf("pars: i=%ld, star[i].r = %e, star[i+1].r = %e, star[i+2].r = %e, star[i+3].r = %e, r = %e\n",
				i, star[i].r, star[i+1].r, star[i+2].r, star[i+3].r, r);
		eprintf("pars: star[i].m=%g star[i+1].m=%g star[i+2].m=%g star[i+3].m=%g\n",
			star[i].m, star[i+1].m, star[i+2].m, star[i+3].m);
		exit_cleanly(-2);
	}

	/* Henon's method of computing the potential using star[].phi */ 
	if (i == 0){ /* I think this is impossible, due to early return earlier,
			    but I am keeping it. -- ato 23:17,  3 Jan 2005 (UTC) */
#ifdef USE_MPI
		henon = (star_phi[1]);
#else
		henon = (star[1].phi);
#endif
	} else {
#ifdef USE_MPI
		henon = (star_phi[i] + (star_phi[i + 1] - star_phi[i]) 
			 * (1.0/star_r[i] - 1.0/r) /
			 (1.0/star_r[i] - 1.0/star_r[i + 1]));
#else
		henon = (star[i].phi + (star[i + 1].phi - star[i].phi) 
			 * (1.0/star[i].r - 1.0/r) /
			 (1.0/star[i].r - 1.0/star[i + 1].r));
#endif
	}
	
	return (henon);
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
#ifdef USE_MPI
	MPI_Finalize();
#endif

	close_buffers();
	free_arrays();

	exit(signal);
}

void free_arrays(void){
	free(mass_pc); free(densities_r); free(no_star_r); 
	free(ke_rad_r); free(ke_tan_r); free(v2_rad_r); free(v2_tan_r);
	free(ave_mass_r); free(mass_r);
	free(star); free(binary);

	/* MPI Stuff */
#ifdef USE_MPI
	free(star_r); free(star_m); free(star_phi);
#endif
}

/* GSL error handler */
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno)
{
	fprintf(stderr, "gsl: %s:%d: ERROR: %s\n", file, line, reason);
	exit_cleanly(gsl_errno);
}

void set_velocities3(void){
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	/* set velocities a la Stodolkiewicz to be able to conserve energy */
	double vold2, vnew2, Unewrold, Unewrnew;
	double Eexcess, exc_ratio;
	double q=0.5; /* q=0.5 -> Stodolkiewicz, q=0 -> Delta U all from rnew */
	double alpha;
	long i;

	Eexcess = 0.0;
#ifdef USE_MPI
	if(myid==0)
#endif
	{
		for (i = 1; i <= clus.N_MAX; i++) {
			/* modify velocities of stars that have only undergone relaxation */
			if (star[i].interacted == 0) {
				Unewrold = potential_serial(star[i].rOld) + PHI_S(star[i].rOld, i);
				Unewrnew = star[i].phi + PHI_S(star[i].r, i);
				vold2 = star[i].vtold*star[i].vtold + 
					star[i].vrold*star[i].vrold;
				/* predict new velocity */
				vnew2 = vold2 + 2.0*(1.0-q)*(star[i].Uoldrold - star[i].Uoldrnew)
					+ 2.0*q*(Unewrold - Unewrnew);
				/* new velocity can be unphysical, so just use value predicted by old potential
					(this is already set in .vr and .vt) */
				if (vnew2 <= 0.0) {
					Eexcess += 0.5*(sqr(star[i].vr)+sqr(star[i].vt)-vnew2)*star[i].m;
				} else {
					/* scale velocity, preserving v_t/v_r */
					alpha = sqrt(vnew2/(sqr(star[i].vr)+sqr(star[i].vt)));
					star[i].vr *= alpha;
					star[i].vt *= alpha;
					/* if there is excess energy added, try to remove at 
						least part of it from this star */

					if(Eexcess > 0 && Eexcess < 0.5*(sqr(star[i].vt)+sqr(star[i].vr))*star[i].m){
						exc_ratio = 
							sqrt( (sqr(star[i].vt)+sqr(star[i].vr)-2*Eexcess/star[i].m)/
									(sqr(star[i].vt)+sqr(star[i].vr)) );
						star[i].vr *= exc_ratio;
						star[i].vt *= exc_ratio;
						Eexcess = 0.0;
					}
				}
			}
		}

		/* keep track of the energy that's vanishing due to our negligence */
		Eoops += -Eexcess * madhoc;
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

/* computes intermediate energies, and transfers "new" dynamical params to the standard variables */
void ComputeIntermediateEnergy(void)
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	long j;

#ifdef USE_MPI 
	//MPI2: Only running till N_MAX for now since no new stars are created, later new loop has to be introduced from N_MAX+1 to N_MAX_NEW.
	for (j=mpiBegin; j<=mpiEnd; j++) {
#else
	/* compute intermediate energies for stars due to change in pot */ 
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
#endif
		/* but do only for NON-Escaped stars */
		if (star[j].rnew < 1.0e6) {
#ifdef USE_MPI
			star[j].EI = sqr(star[j].vr) + sqr(star[j].vt) + star_phi[j] - potential(star[j].rnew);
#else
			star[j].EI = sqr(star[j].vr) + sqr(star[j].vt) + star[j].phi - potential(star[j].rnew);
#endif
		}
	}
	
#ifdef USE_MPI
	for (j=mpiBegin; j<=mpiEnd; j++) {
#else
	/* Transferring new positions to .r, .vr, and .vt from .rnew, .vrnew, and .vtnew */
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
#endif
#ifdef USE_MPI
		star[j].rOld = star_r[j];
		star_r[j] = star[j].rnew;
#else
		star[j].rOld = star[j].r;
		star[j].r = star[j].rnew;
#endif
		star[j].vr = star[j].vrnew;
		star[j].vt = star[j].vtnew;
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

long CheckStop(struct tms tmsbufref) {
	struct tms tmsbuf;
	long tspent;

	times(&tmsbuf);
	tspent  = tmsbuf.tms_utime-tmsbufref.tms_utime;
	tspent += tmsbuf.tms_stime-tmsbufref.tms_stime;
	tspent += tmsbuf.tms_cutime-tmsbufref.tms_cutime;
	tspent += tmsbuf.tms_cstime-tmsbufref.tms_cstime;
	tspent /= sysconf(_SC_CLK_TCK);
	tspent /= 60;

	if (tspent >= MAX_WCLOCK_TIME) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("MAX_WCLOCK_TIME exceeded ... Terminating.\n");
		return (1);
	}
	
	if (tcount >= T_MAX_COUNT) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("No. of timesteps > T_MAX_COUNT ... Terminating.\n");
		return (1);
	}

	if (TotalTime >= T_MAX) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("TotalTime > T_MAX ... Terminating.\n");
		return (1);
	}

	if (TotalTime / (1.0e3*MEGA_YEAR) >= T_MAX_PHYS) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("TotalTime > T_MAX_PHYS ... Terminating.\n");
		return (1);
	}

	/* Stop if cluster is disrupted -- N_MAX is too small */
	/* if (clus.N_MAX < (0.02 * clus.N_STAR)) { */
	if (clus.N_MAX < (0.005 * clus.N_STAR)) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("N_MAX < 0.005 * N_STAR ... Terminating.\n");
		return (1);
	}

	/* Stop if Etotal > 0 */
	if (Etotal.K + Etotal.P > 0.0) {
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("Etotal > 0 ... Terminating.\n");
		return (1);
	}

	/* Stop at core collapse, if requested.  Core collapsed is detected by looking at the
	   number of core stars.  Anything less than ~100 represents core collapse, independent of 
	   the number of initial cluster stars (see Heggie & Hut's "Gravitational Million Body
	   Problem").  The number can get arbitrarily small if there is no three-body
	   binary formation. */
	if (STOPATCORECOLLAPSE) {
		if (N_core <= 100.0) {
			print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
			diaprintf("N_core < 100.0; terminating.\n");
			return (1);
		}
	}

	/* Output some snapshots near core collapse 
	 * (if core density is high enough) */
	if (SNAPSHOT_CORE_COLLAPSE) { 
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
	}

	/* added by ato 
	 * to try to take snapshots for core bounce as well. 
	 * idea is if we reduced core density by 10 percent the
	 * last time we took snapshot, take another one and adjust
	 * parameters to take further snapshots if further collapse
	 * occurs */
	if (SNAPSHOT_CORE_BOUNCE) {
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
		print_2Dsnapshot();
#ifdef USE_MPI
		if(myid==0)
#endif
		diaprintf("Terminal Energy reached... Terminating.\n");
		return (1);
	}
	return (0); /* NOT stopping time yet */
}


#ifdef USE_MPI
void mpi_ComputeEnergy(void)
{
	strcpy(funcName, __FUNCTION__);
	double Etotal_tot = 0.0;
	double Etotal_K = 0.0;
	double Etotal_P = 0.0;
	double Etotal_Eint = 0.0;
	double Etotal_Eb = 0.0;
	double phi0 = 0.0;

	int i;
	for (i=mpiBegin; i<=mpiEnd; i++) {
		star[i].E = star_phi[i] + 0.5 * (sqr(star[i].vr) + sqr(star[i].vt));
		star[i].J = star_r[i] * star[i].vt;		
	}

	phi0 = star_phi[0];

	for (i=mpiBegin; i<=mpiEnd; i++) {
		Etotal_K += 0.5 * (sqr(star[i].vr) + sqr(star[i].vt)) * star_m[i] / clus.N_STAR;
		Etotal_P += star_phi[i] * star_m[i] / clus.N_STAR;
		//change this so that star[0] is not accesses by all procs. do on root node
		Etotal_P += phi0 * cenma.m*madhoc/ clus.N_STAR;

		//MPI2: Binary indices are involved. Ignore.
		if (star[i].binind == 0) {
			Etotal_Eint += star[i].Eint;
		} else if (binary[star[i].binind].inuse) {
			Etotal_Eb += -(binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
				(2.0 * binary[star[i].binind].a);
			Etotal_Eint += binary[star[i].binind].Eint1 + binary[star[i].binind].Eint2;
		}
	}
	
	Etotal_P *= 0.5;
	Etotal_tot = Etotal_K + Etotal_P + Etotal_Eint + Etotal_Eb;

	MPI_Reduce(&Etotal_K, &Etotal.K, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&Etotal_P, &Etotal.P, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&Etotal_Eint, &Etotal.Eint, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&Etotal_Eb, &Etotal.Eb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&Etotal_tot, &Etotal.tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		

	MPI_Bcast(&Etotal.K, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Etotal.P, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Etotal.Eint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Etotal.Eb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Etotal.tot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	Etotal.tot += cenma.E + Eescaped + Ebescaped + Eintescaped;
}
#endif

void ComputeEnergy(void)
{
	strcpy(funcName, __FUNCTION__);
	long i;
	
	Etotal.tot = 0.0;
	Etotal.K = 0.0;
	Etotal.P = 0.0;
	Etotal.Eint = 0.0;
	Etotal.Eb = 0.0;

	//split into 2 parts. till next comment can be done on all nodes
	for (i=1; i<=clus.N_MAX; i++) {
		star[i].E = star[i].phi + 0.5 * (sqr(star[i].vr) + sqr(star[i].vt));
		star[i].J = star[i].r * star[i].vt;
		
	//MPI2:this is a reduce - write as new function
		Etotal.K += 0.5 * (sqr(star[i].vr) + sqr(star[i].vt)) * star[i].m / clus.N_STAR;
		Etotal.P += star[i].phi * star[i].m / clus.N_STAR;
		//change this so that star[0] is not accesses by all procs. do on root node
		Etotal.P += star[0].phi * cenma.m*madhoc/ clus.N_STAR;
		
		if (star[i].binind == 0) {
			Etotal.Eint += star[i].Eint;
		} else if (binary[star[i].binind].inuse) {
			Etotal.Eb += -(binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
				(2.0 * binary[star[i].binind].a);
			Etotal.Eint += binary[star[i].binind].Eint1 + binary[star[i].binind].Eint2;
		}
	}
	
	Etotal.P *= 0.5;
	Etotal.tot = Etotal.K + Etotal.P + Etotal.Eint + Etotal.Eb + cenma.E + Eescaped + Ebescaped + Eintescaped;
}

#ifdef USE_MPI
long mpi_potential_calculate(void) {
	strcpy(funcName, __FUNCTION__);
	long k;
	double mprev;

	/* count up all the mass and set N_MAX */
	k = 1;
	mprev = 0.0;
	
	//while (star_r[k] < SF_INFINITY && k <= clus.N_STAR_NEW) {
	//MPI2: Temporarily changing this to N_MAX.
	while (star_r[k] < SF_INFINITY && k <= clus.N_MAX) {
		mprev += star_m[k];
		/* I guess NaNs do happen... */
		if(isnan(mprev)){
			eprintf("NaN (2) detected\n");
			exit_cleanly(-1);
		}
		k++;
	}


	/* New N_MAX */
	clus.N_MAX = k - 1;

	/* update central BH mass */
	cenma.m= cenma.m_new;

	/* New total Mass; This IS correct for multiple components */
	Mtotal = mprev * madhoc + cenma.m * madhoc;	
        dprintf("Mtotal is %lf, cenma.m is %lf, madhoc is %lg, mprev is %lf\n", Mtotal, cenma.m, madhoc, mprev);

	/* Compute new tidal radius using new Mtotal */

	Rtidal = orbit_r * pow(Mtotal, 1.0 / 3.0);

	/* zero boundary star first for safety */
	zero_star(clus.N_MAX + 1);

	star_r[clus.N_MAX + 1] = SF_INFINITY;
	star_phi[clus.N_MAX + 1] = 0.0;

	mprev = Mtotal;
	for (k = clus.N_MAX; k >= 1; k--) {/* Recompute potential at each r */
		star_phi[k] = star_phi[k + 1] - mprev * (1.0 / star_r[k] - 1.0 / star_r[k + 1]);
		mprev -= star_m[k] / clus.N_STAR;
		if (isnan(star_phi[k])) {
		  eprintf("NaN in phi[%li] detected\n", k);
		  eprintf("phi[k+1]=%g mprev=%g, r[k]=%g, r[k+1]=%g, m[k]=%g, clus.N_STAR=%li\n", 
		  	star_phi[k + 1], mprev, star_r[k], star_r[k + 1], star_m[k], clus.N_STAR);
		  //exit_cleanly(-1);
		}
	}

	star_phi[0] = star_phi[1]+ cenma.m*madhoc/star_r[1]; /* U(r=0) is U_1 */
	if (isnan(star_phi[0])) {
		eprintf("NaN in phi[0] detected\n");
		exit_cleanly(-1);
	}

	return (clus.N_MAX);
}
#endif

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
	strcpy(funcName, __FUNCTION__);
	long k;
	double mprev;

	/* count up all the mass and set N_MAX */
	k = 1;
	mprev = 0.0;

	while (star[k].r < SF_INFINITY && k <= clus.N_STAR_NEW) {
		mprev += star[k].m;
		/* I guess NaNs do happen... */
		if(isnan(mprev)){
			eprintf("NaN (2) detected\n");
			exit_cleanly(-1);
		}
		k++;
	}

	/* New N_MAX */
	clus.N_MAX = k - 1;

	/* update central BH mass */
	cenma.m= cenma.m_new;

	/* New total Mass; This IS correct for multiple components */
	Mtotal = mprev * madhoc + cenma.m * madhoc;	
	dprintf("Mtotal is %lf, cenma.m is %lf, madhoc is %lg, mprev is %lf\n", Mtotal, cenma.m, madhoc, mprev);

	/* Compute new tidal radius using new Mtotal */

	Rtidal = orbit_r * pow(Mtotal, 1.0 / 3.0);

	/* zero boundary star first for safety */
	zero_star(clus.N_MAX + 1);

	star[clus.N_MAX + 1].r = SF_INFINITY;
	star[clus.N_MAX + 1].phi = 0.0;

	mprev = Mtotal;
	for (k = clus.N_MAX; k >= 1; k--) {/* Recompute potential at each r */
		star[k].phi = star[k + 1].phi - mprev * (1.0 / star[k].r - 1.0 / star[k + 1].r);
		mprev -= star[k].m / clus.N_STAR;
		if (isnan(star[k].phi)) {
		  eprintf("NaN in phi[%li] detected\n", k);
		  eprintf("phi[k+1]=%g mprev=%g, r[k]=%g, r[k+1]=%g, m[k]=%g, clus.N_STAR=%li\n", 
		  	star[k + 1].phi, mprev, star[k].r, star[k + 1].r, star[k].m, clus.N_STAR);
		  //exit_cleanly(-1);
		}
	}

	/*for (k = 1; k <= clus.N_MAX; k++){
		star[k].phi -= cenma.m * madhoc / star[k].r;
		if(isnan(star[k].phi)){
			eprintf("NaN detected\n");
			exit_cleanly(-1);
		}
	}*/
	
	star[0].phi = star[1].phi+ cenma.m*madhoc/star[1].r; /* U(r=0) is U_1 */
	if (isnan(star[0].phi)) {
		eprintf("NaN in phi[0] detected\n");
		exit_cleanly(-1);
	}

	return (clus.N_MAX);
}

#define GENSEARCH_NAME 				m_binsearch
#define GENSEARCH_TYPE 				double
#define GENSEARCH_KEYTYPE			double
#define GENSEARCH_GETKEY(a)			a
#define GENSEARCH_COMPAREKEYS(k1, k2)	k1 < k2

#include "common/gensearch.h"

int find_stars_mass_bin(double smass){
	/* find the star[i]'s mass bin */
	/* return -1 on failure */
	int bn;

	if ( (smass < mass_bins[0]) || 
	     (smass > mass_bins[NO_MASS_BINS-1]) ) return -1;
	
	bn = m_binsearch(mass_bins, 0, NO_MASS_BINS-1, smass);
	return bn;
}

void comp_multi_mass_percent(){
	/* computing the Lagrange radii for various mass bins */
	/* mass bins are stored in the array mass_bins[NO_MASS_BINS] */
	double *mtotal_inbin; // total mass in each mass bin
	long *number_inbin;   // # of stars in each mass bin
	double *mcount_inbin, *r_inbin;
	long *ncount_inbin;
	int *star_bins;       // array holding which bin each star is in
	double **rs, **percents;
	long i;

	/* GSL interpolation function and accelerators. See:
	 * http://sources.redhat.com/gsl/ref/gsl-ref_26.html#SEC391*/
	gsl_interp_accel **acc;
	gsl_spline **spline;

	if (NO_MASS_BINS <=1) return;

	mtotal_inbin = (double *) calloc(NO_MASS_BINS, sizeof(double));
	number_inbin = (long *) calloc(NO_MASS_BINS, sizeof(long));
	r_inbin = (double *) calloc(NO_MASS_BINS, sizeof(double));
	star_bins = (int *) malloc((clus.N_MAX+2)*sizeof(int));
	for (i = 1; i <= clus.N_MAX; i++) {
		star_bins[i] = find_stars_mass_bin(star[i].m/SOLAR_MASS_DYN);
		if (star_bins[i] == -1) continue; /* -1: star isn't in legal bin */
		mtotal_inbin[star_bins[i]] += star[i].m; /* no unit problem, since 
								        we are interested in
							   	        percentage only. */
		number_inbin[star_bins[i]]++;
		r_inbin[star_bins[i]] = star[i].r;
	}
	/* populate arrays rs[NO_MASS_BINS][j] and percents[][] */
	rs = (double **) malloc(NO_MASS_BINS*sizeof(double *));
	percents = (double **) malloc(NO_MASS_BINS*sizeof(double *));
	for(i=0; i<NO_MASS_BINS; i++){
		/* +1 below is to accomodate rs=0 <-> percents=0 point */
		rs[i] = (double *) malloc((number_inbin[i]+1)*sizeof(double));
		percents[i] = (double *) malloc((number_inbin[i]+1)*sizeof(double));
	}
	for(i=0; i<NO_MASS_BINS; i++){
		/* at r=0 there is 0% of mass */
		rs[i][0] = percents[i][0] = 0.0;
	}
	
	mcount_inbin = (double *) calloc(NO_MASS_BINS, sizeof(double));
	ncount_inbin = (long *) calloc(NO_MASS_BINS, sizeof(long));
	for (i = 1; i <= clus.N_MAX; i++) {
		int sbin = star_bins[i];
		if (sbin == -1) continue;
		mcount_inbin[sbin] += star[i].m;
		ncount_inbin[sbin]++;
		rs[sbin][ncount_inbin[sbin]] = star[i].r;
		percents[sbin][ncount_inbin[sbin]] = 
				mcount_inbin[sbin]/mtotal_inbin[sbin];
	}
	free(mcount_inbin); free(ncount_inbin);
	
	acc = (gsl_interp_accel **) malloc(NO_MASS_BINS*sizeof(gsl_interp_accel));
	spline = (gsl_spline **) malloc(NO_MASS_BINS*sizeof(gsl_spline));
	for(i=0; i<NO_MASS_BINS; i++){
		if((number_inbin[i] == 1) || (number_inbin[i] == 0)) continue;
		acc[i] = gsl_interp_accel_alloc();
		/* change gsl_interp_linear to gsl_interp_cspline below,
		 * for spline interpolation; however, this may result in
		 * negative LR values for small mass percentages! */ 
		/* +1's below are to accomodate rs=0 <-> percents=0 point */
		spline[i] = gsl_spline_alloc(gsl_interp_linear, number_inbin[i]+1);
		gsl_spline_init (spline[i], percents[i], rs[i], number_inbin[i]+1);
	}
	for(i=0; i<NO_MASS_BINS; i++){
		free(rs[i]); free(percents[i]);
	}
	free(rs); free(percents);
	
	/* fill multi_mass_r[][] by calling gsl_spline_eval() */
	for (i = 0; i <NO_MASS_BINS; i++) {
		int mcnt;
		if (number_inbin[i] == 0) {
			continue;
		} else if (number_inbin[i] == 1) {
			for(mcnt=0; mcnt<MASS_PC_COUNT; mcnt++){
				multi_mass_r[i][mcnt] = r_inbin[i];
			}
		} else {
			for(mcnt=0; mcnt<MASS_PC_COUNT; mcnt++){
				multi_mass_r[i][mcnt] = 
					gsl_spline_eval(spline[i], mass_pc[mcnt], acc[i]);
			}
		}
	}
	for(i=0; i<NO_MASS_BINS; i++){
		if((number_inbin[i] == 1) || (number_inbin[i] == 0)) continue;
		gsl_interp_accel_free(acc[i]); 
		gsl_spline_free(spline[i]);
	}
	free(acc); free(spline);
	
	free(mtotal_inbin); free(number_inbin); free(r_inbin); 
	free(star_bins);
}
		
void comp_mass_percent(){
	double mprev, ke_rad_prev=0.0, ke_tan_prev=0.0, v2_rad_prev=0.0, v2_tan_prev=0.0;
	long int k, mcount;

	/* Computing radii containing mass_pc[i] % of the mass */
	if (MASS_PC_BH_INCLUDE) {
		mprev = cenma.m * madhoc;
		for(mcount=0; mcount<MASS_PC_COUNT; mcount++){
			if ( mprev/Mtotal > mass_pc[mcount] ) {
				mass_r[mcount] = MINIMUM_R;
				ave_mass_r[mcount] = 0.0;
				no_star_r[mcount] = 0;
				densities_r[mcount] = 0.0;
				ke_rad_r[mcount] = 0.0;
				ke_tan_r[mcount] = 0.0;
				v2_rad_r[mcount] = 0.0;
				v2_tan_r[mcount] = 0.0;
			} else {
				break;
			}
		}
	} else {
		mprev= 0.;
		mcount=0;
	}

	for (k = 1; k <= clus.N_MAX; k++) {	/* Only need to count up to N_MAX */
		mprev += star[k].m / clus.N_STAR;
		ke_rad_prev += 0.5 * star[k].m * madhoc * star[k].vr * star[k].vr;
		ke_tan_prev += 0.5 * star[k].m * madhoc * star[k].vt * star[k].vt;
		v2_rad_prev += star[k].vr * star[k].vr;
		v2_tan_prev += star[k].vt * star[k].vt;
		if (mprev / Mtotal > mass_pc[mcount]) {
			mass_r[mcount] = star[k].r;
			ave_mass_r[mcount] = mprev/Mtotal/k*initial_total_mass;
			no_star_r[mcount] = k;
			densities_r[mcount] = mprev*clus.N_STAR/
				(4/3*3.1416*pow(star[k].r,3));
			ke_rad_r[mcount] = ke_rad_prev;
			ke_tan_r[mcount] = ke_tan_prev;
			v2_rad_r[mcount] = v2_rad_prev;
			v2_tan_r[mcount] = v2_tan_prev;
			mcount++;
			if (mcount == MASS_PC_COUNT)
				break;
		}
	}
}

/* The potential computed using the star[].phi computed at the star 
   locations in star[].r sorted by increasing r. */
double fastpotential(double r, long kmin, long kmax) {
	long i;
	double henon;

	/* root finding using indexed values of sr[] & bisection */
	if (r < star[1].r)
		return (star[1].phi);

	i =  FindZero_r(kmin, kmax, r);
	
	if(star[i].r > r || star[i+1].r < r){
		eprintf("binary search (FindZero_r) failed!!\n");
		eprintf("pars: i=%ld, star[i].r = %e, star[i+1].r = %e, star[i+2].r = %e, star[i+3].r = %e, r = %e\n",
				i, star[i].r, star[i+1].r, star[i+2].r, star[i+3].r, r);
		eprintf("pars: star[i].m=%g star[i+1].m=%g star[i+2].m=%g star[i+3].m=%g\n",
			star[i].m, star[i+1].m, star[i+2].m, star[i+3].m);
		exit_cleanly(-2);
	}

	/* Henon's method of computing the potential using star[].phi */ 
	if (i == 0){ /* I think this is impossible, due to early return earlier,
			    but I am keeping it. -- ato 23:17,  3 Jan 2005 (UTC) */
		henon = (star[1].phi);
	} else {
		henon = (star[i].phi + (star[i + 1].phi - star[i].phi) 
			 * (1.0/star[i].r - 1.0/r) /
			 (1.0/star[i].r - 1.0/star[i + 1].r));
	}
	
	return (henon);
}

long check_if_r_around_last_index(long last_index, double r) {
   long index_found, i;

   //printf("Entering the routine. Last index is %li", last_index);
   index_found= -1;
   for (i=-1; i<2; i++) {
      if (((last_index+i) >= 0) && ((last_index+i+1) <= clus.N_STAR)) {
         if ((star[last_index+i].r < r) && (star[last_index+i+1].r > r)) {
	 	//printf("found it!");
            index_found= last_index+ i;
            break;
         };
      };
    };

   return (index_found);
};

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
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	long i, j, k;
	
	/* update total number, mass, and binding energy of binaries in cluster */
	N_b = 0;
	M_b = 0.0;
	E_b = 0.0;
#ifdef USE_MPI
	for (i=mpiBegin; i<=mpiEnd; i++) 
#else
	for (i=1; i<=clus.N_MAX; i++)
#endif
	{
		j = i;
		k = star[j].binind;
		if (k != 0) {
			N_b++;
			M_b += star[j].m;
			if (binary[k].inuse){
				E_b += (binary[k].m1/clus.N_STAR) * (binary[k].m2/clus.N_STAR) / (2.0 * binary[k].a);
			}
		}
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

/* set the units */
void units_set(void)
{
	units.m = cfd.Mclus * MSUN;
	units.l = cfd.Rvir * PARSEC;
	units.t = pow(units.l * pow(G * units.m, -1.0/3.0), 3.0/2.0);
	units.E = G * sqr(units.m) / units.l;
	units.mstar = units.m / ((double) clus.N_STAR);
	
	/* derivative quantities */
	MEGA_YEAR = log(GAMMA * clus.N_STAR)/(clus.N_STAR * units.t) * 1.0e6 * YEAR;
	SOLAR_MASS_DYN = clus.N_STAR / units.m * MSUN;

	/* Masses such as star.m and binary.m1 and binary.m2 are not stored in code units, 
	   but rather in code units * clus.N_STAR.  This means that whenever you want to
	   calculate a quantity that involves masses and lengths, or masses and times, you
	   have to divide any masses in the expression by clus.N_STAR.  This is that 
	   factor.
	*/
	madhoc = 1.0/((double) clus.N_STAR);

	/* print out diagnostic information */
#ifdef USE_MPI
	if(myid==0)
#endif
	{
		diaprintf("METALLICITY= %g\n",METALLICITY);
		diaprintf("MEGA_YEAR=%g\n", MEGA_YEAR);
		diaprintf("SOLAR_MASS_DYN=%g\n", SOLAR_MASS_DYN);
		diaprintf("initial_total_mass=%g\n", initial_total_mass);
		diaprintf("units.t=%g YEAR\n", units.t/YEAR);
		diaprintf("units.m=%g MSUN\n", units.m/MSUN);
		diaprintf("units.mstar=%g MSUN\n", units.mstar/MSUN);
		diaprintf("units.l=%g PARSEC\n", units.l/PARSEC);
		diaprintf("units.E=%g erg\n", units.E);
		diaprintf("t_rel=%g YEAR\n", units.t * clus.N_STAR / log(GAMMA * clus.N_STAR) / YEAR);
	}
}

#ifdef USE_MPI
void mpi_central_calculate(void)
{
	mpi_central_calculate1();
	mpi_central_calculate2();
}
#endif

#ifdef USE_MPI
/* calculate central quantities */
void mpi_central_calculate1(void)
{
	double m=0.0, *rhoj, mrho, Vrj, rhojsum, rhoj2sum;
	long J=6, i, j, jmin, jmax, nave;

	/* average over all stars out to half-mass radius */
	if(myid==0) {
		nave = 1;
		while (m < 0.5 * Mtotal) {
		/*MPI2: Using global m array*/
			m += star_m[nave] / clus.N_STAR;
			nave++;
		}
	}
	MPI_Bcast(&nave, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	/* exit if not enough stars */
	if (clus.N_STAR <= 2*J || nave >= clus.N_STAR-6) {
		eprintf("clus.N_STAR <= 2*J || nave >= clus.N_STAR-6\n");
		exit_cleanly(-1);
	}

	/* allocate array for local density calculations */
	rhoj = (double *) malloc((nave+1) * sizeof(double));

	int mpiBeginLocal, mpiEndLocal;
	mpiFindIndicesCustom( nave, 20, myid, &mpiBeginLocal, &mpiEndLocal );
	//mpiFindIndicesSpecial( nave, &mpiBeginLocal, &mpiEndLocal );

	/* calculate rhoj's (Casertano & Hut 1985) */
	for (i=mpiBeginLocal; i<=mpiEndLocal; i++) {
		jmin = MAX(i-J/2, 1);
		jmax = jmin + J;
		mrho = 0.0;
		/* this is equivalent to their J-1 factor for the case of equal masses,
		   and seems like a good generalization for unequal masses */
		for (j=jmin+1; j<=jmax-1; j++) {
			mrho += star_m[j] * madhoc;
		}
		Vrj = 4.0/3.0 * PI * (fb_cub(star_r[jmax]) - fb_cub(star_r[jmin]));
		rhoj[i] = mrho / Vrj;
	}

	/* calculate core quantities using density weighted averages (note that in 
	   Casertano & Hut (1985) only rho and rc are analyzed and tested) */

	double mpi_rhojsum, mpi_rhoj2sum, mpi_rc_nb, mpi_c_rho, mpi_c_vrms, mpi_c_rc, mpi_c_mave;
	mpi_rhojsum = 0.0;
	mpi_rhoj2sum = 0.0;
	mpi_rc_nb = 0.0;
	mpi_c_rho = 0.0;
	mpi_c_vrms = 0.0;
	mpi_c_rc = 0.0;
	mpi_c_mave = 0.0;

	for (i=mpiBeginLocal; i<=mpiEndLocal; i++) {
		mpi_rhojsum += rhoj[i];
		mpi_rhoj2sum += sqr(rhoj[i]);
		mpi_c_rho += sqr(rhoj[i]);
		mpi_c_vrms += rhoj[i] * (sqr(star[i].vr) + sqr(star[i].vt));
		mpi_c_rc += rhoj[i] * star_r[i];
		mpi_rc_nb += sqr(rhoj[i] * star_r[i]);
		mpi_c_mave += rhoj[i] * star_m[i] * madhoc;
	}

	MPI_Reduce(&mpi_rhojsum, &rhojsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&mpi_rhoj2sum, &rhoj2sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&mpi_c_rho, &central.rho, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&mpi_c_vrms, &central.v_rms, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&mpi_c_rc, &central.rc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
	MPI_Reduce(&mpi_c_mave, &central.m_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		

	if(myid==0) {
		central.rho /= rhojsum;
		/* correction for inherent bias in estimator */
		central.rho *= 4.0/5.0;
		/* and now correct for the fact this estimate of density is systematically smaller than the
			theoretical by a factor of about 2 (for a range of King models and for the Plummer model) */
		central.rho *= 2.0;
		central.v_rms /= rhojsum;
		central.v_rms = sqrt(central.v_rms);
		central.rc /= rhojsum;
		rc_nb = sqrt(rc_nb/rhoj2sum);
		central.m_ave /= rhojsum;

		/* quantities derived from averages */
		central.n = central.rho / central.m_ave;
		central.rc_spitzer = sqrt(3.0 * sqr(central.v_rms) / (4.0 * PI * central.rho));
	}
	free(rhoj);
}
#endif

#ifdef USE_MPI
void mpi_central_calculate2(void) {
	double Msincentral, Mbincentral, Vcentral, rcentral;
	long i, Ncentral;

	if(myid==0) {
		/* calculate other quantities using old method */
		Ncentral = 0;
		Msincentral = 0.0;
		Mbincentral = 0.0;
		central.N_sin = 0;
		central.N_bin = 0;
		central.v_sin_rms = 0.0;
		central.v_bin_rms = 0.0;
		central.w2_ave = 0.0;
		central.R2_ave = 0.0;
		central.mR_ave = 0.0;
		central.a_ave = 0.0;
		central.a2_ave = 0.0;
		central.ma_ave = 0.0;

		for (i=1; i<=MIN(NUM_CENTRAL_STARS, clus.N_STAR); i++) {
			/* for (i=1; i<=MIN((long) N_core, clus.N_STAR); i++) { */
			Ncentral++;
			/* use only code units here, so always divide star[].m by clus.N_STAR */
			central.w2_ave += 2.0 * star_m[i] / ((double) clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt));

			if (star[i].binind == 0) {
				central.N_sin++;
				Msincentral += star_m[i] / ((double) clus.N_STAR);
				central.v_sin_rms += sqr(star[i].vr) + sqr(star[i].vt);
				central.R2_ave += sqr(star[i].rad);
				central.mR_ave += star_m[i] / ((double) clus.N_STAR) * star[i].rad;
			} else {
				central.N_bin++;
				Mbincentral += star_m[i] / ((double) clus.N_STAR);
				central.v_bin_rms += sqr(star[i].vr) + sqr(star[i].vt);
				central.a_ave += binary[star[i].binind].a;
				central.a2_ave += sqr(binary[star[i].binind].a);
				central.ma_ave += star_m[i] / ((double) clus.N_STAR) * binary[star[i].binind].a;
			}
		}

		/* object quantities */
		rcentral = star_r[Ncentral + 1];
		Vcentral = 4.0/3.0 * PI * cub(rcentral);
		central.w2_ave /= central.m_ave * ((double) Ncentral);

		/* single star quantities */
		central.n_sin = ((double) central.N_sin) / Vcentral;
		central.rho_sin = Msincentral / Vcentral;
		if (central.N_sin != 0) {
			central.m_sin_ave = Msincentral / ((double) central.N_sin);
			central.v_sin_rms = sqrt(central.v_sin_rms / ((double) central.N_sin));
			central.R2_ave /= ((double) central.N_sin);
			central.mR_ave /= ((double) central.N_sin);
		} else {
			central.m_sin_ave = 0.0;
			central.v_sin_rms = 0.0;
			central.R2_ave = 0.0;
			central.mR_ave = 0.0;
		}

		/* binary star quantities */
		central.n_bin = ((double) central.N_bin) / Vcentral;
		central.rho_bin = Mbincentral / Vcentral;
		if (central.N_bin != 0) {
			central.m_bin_ave = Mbincentral / ((double) central.N_bin);
			central.v_bin_rms = sqrt(central.v_bin_rms / ((double) central.N_bin));
			central.a_ave /= ((double) central.N_bin);
			central.a2_ave /= ((double) central.N_bin);
			central.ma_ave /= ((double) central.N_bin);
		} else {
			central.m_bin_ave = 0.0;
			central.v_bin_rms = 0.0;
			central.a_ave = 0.0;
			central.a2_ave = 0.0;
			central.ma_ave = 0.0;
		}
	}

	MPI_Bcast(&central, sizeof(central_t), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* set global code variables */
	v_core = central.v_rms;
	rho_core = central.rho;
	core_radius = central.rc;
	N_core = 4.0 / 3.0 * PI * cub(core_radius) * (central.n / 2.0);
	//Sourav
	N_core_nb = 4.0 / 3.0 * PI * cub(rc_nb) * (central.n / 2.0);

	/* core relaxation time, Spitzer (1987) eq. (2-62) */
	Trc = 0.065 * cub(central.v_rms) / (central.rho * central.m_ave);

	/* set global variables that are used throughout the code */
	rho_core_single = central.rho_sin;
	rho_core_bin = central.rho_bin;
	
}
#endif

/* calculate central quantities */
void central_calculate(void)
{
	double m=0.0, *rhoj, mrho, Vrj, rhojsum, Msincentral, Mbincentral, Vcentral, rcentral;
	long J=6, i, j, jmin, jmax, nave, Ncentral;
	double rhoj2sum;

	/* average over all stars out to half-mass radius */
	nave = 1;
	while (m < 0.5 * Mtotal) {
		m += star[nave].m / clus.N_STAR;
		nave++;
	}
	
	/* DEBUG */
	/* fprintf(stderr, "nave=%ld\n", nave); */
	/* DEBUG */

	/* exit if not enough stars */
	if (clus.N_STAR <= 2*J || nave >= clus.N_STAR-6) {
		eprintf("clus.N_STAR <= 2*J || nave >= clus.N_STAR-6\n");
		exit_cleanly(-1);
	}

	/* allocate array for local density calculations */
	rhoj = (double *) malloc((nave+1) * sizeof(double));

	/* calculate rhoj's (Casertano & Hut 1985) */
	for (i=1; i<=nave; i++) {
		jmin = MAX(i-J/2, 1);
		jmax = jmin + J;
		mrho = 0.0;
		/* this is equivalent to their J-1 factor for the case of equal masses,
		   and seems like a good generalization for unequal masses */
		for (j=jmin+1; j<=jmax-1; j++) {
			mrho += star[j].m * madhoc;
		}
		Vrj = 4.0/3.0 * PI * (fb_cub(star[jmax].r) - fb_cub(star[jmin].r));
		rhoj[i] = mrho / Vrj;
	}

	/* calculate core quantities using density weighted averages (note that in 
	   Casertano & Hut (1985) only rho and rc are analyzed and tested) */
	rhojsum = 0.0;
	rhoj2sum = 0.0;
	rc_nb = 0.0;
	central.rho = 0.0;
	central.v_rms = 0.0;
	central.rc = 0.0;
	central.m_ave = 0.0;
	for (i=1; i<=nave; i++) {
		rhojsum += rhoj[i];
		rhoj2sum += sqr(rhoj[i]);
		central.rho += sqr(rhoj[i]);
		central.v_rms += rhoj[i] * (sqr(star[i].vr) + sqr(star[i].vt));
		central.rc += rhoj[i] * star[i].r;
		rc_nb += sqr(rhoj[i] * star[i].r);
		central.m_ave += rhoj[i] * star[i].m * madhoc;
	}
	central.rho /= rhojsum;
	/* correction for inherent bias in estimator */
	central.rho *= 4.0/5.0;
	/* and now correct for the fact this estimate of density is systematically smaller than the
	   theoretical by a factor of about 2 (for a range of King models and for the Plummer model) */
	central.rho *= 2.0;
	central.v_rms /= rhojsum;
	central.v_rms = sqrt(central.v_rms);
	central.rc /= rhojsum;
	rc_nb = sqrt(rc_nb/rhoj2sum);
	central.m_ave /= rhojsum;
	
	/* quantities derived from averages */
	central.n = central.rho / central.m_ave;
	central.rc_spitzer = sqrt(3.0 * sqr(central.v_rms) / (4.0 * PI * central.rho));

	/* set global code variables */
	v_core = central.v_rms;
	rho_core = central.rho;
	core_radius = central.rc;
	N_core = 4.0 / 3.0 * PI * cub(core_radius) * (central.n / 2.0);
	//Sourav
	N_core_nb = 4.0 / 3.0 * PI * cub(rc_nb) * (central.n / 2.0);
	/* core relaxation time, Spitzer (1987) eq. (2-62) */
	Trc = 0.065 * cub(central.v_rms) / (central.rho * central.m_ave);

	/* calculate other quantities using old method */
	Ncentral = 0;
	Msincentral = 0.0;
	Mbincentral = 0.0;
	central.N_sin = 0;
	central.N_bin = 0;
	central.v_sin_rms = 0.0;
	central.v_bin_rms = 0.0;
	central.w2_ave = 0.0;
	central.R2_ave = 0.0;
	central.mR_ave = 0.0;
	central.a_ave = 0.0;
	central.a2_ave = 0.0;
	central.ma_ave = 0.0;
	for (i=1; i<=MIN(NUM_CENTRAL_STARS, clus.N_STAR); i++) {
	/* for (i=1; i<=MIN((long) N_core, clus.N_STAR); i++) { */
		Ncentral++;
		/* use only code units here, so always divide star[].m by clus.N_STAR */
		central.w2_ave += 2.0 * star[i].m / ((double) clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt));

		if (star[i].binind == 0) {
			central.N_sin++;
			Msincentral += star[i].m / ((double) clus.N_STAR);
			central.v_sin_rms += sqr(star[i].vr) + sqr(star[i].vt);
			central.R2_ave += sqr(star[i].rad);
			central.mR_ave += star[i].m / ((double) clus.N_STAR) * star[i].rad;
		} else {
			central.N_bin++;
			Mbincentral += star[i].m / ((double) clus.N_STAR);
			central.v_bin_rms += sqr(star[i].vr) + sqr(star[i].vt);
			central.a_ave += binary[star[i].binind].a;
			central.a2_ave += sqr(binary[star[i].binind].a);
			central.ma_ave += star[i].m / ((double) clus.N_STAR) * binary[star[i].binind].a;
		}
	}
	/* object quantities */
	rcentral = star[Ncentral + 1].r;
	Vcentral = 4.0/3.0 * PI * cub(rcentral);
	central.w2_ave /= central.m_ave * ((double) Ncentral);
	
	/* single star quantities */
	central.n_sin = ((double) central.N_sin) / Vcentral;
	central.rho_sin = Msincentral / Vcentral;
	if (central.N_sin != 0) {
		central.m_sin_ave = Msincentral / ((double) central.N_sin);
		central.v_sin_rms = sqrt(central.v_sin_rms / ((double) central.N_sin));
		central.R2_ave /= ((double) central.N_sin);
		central.mR_ave /= ((double) central.N_sin);
	} else {
		central.m_sin_ave = 0.0;
		central.v_sin_rms = 0.0;
		central.R2_ave = 0.0;
		central.mR_ave = 0.0;
	}
	
	/* binary star quantities */
	central.n_bin = ((double) central.N_bin) / Vcentral;
	central.rho_bin = Mbincentral / Vcentral;
	if (central.N_bin != 0) {
		central.m_bin_ave = Mbincentral / ((double) central.N_bin);
		central.v_bin_rms = sqrt(central.v_bin_rms / ((double) central.N_bin));
		central.a_ave /= ((double) central.N_bin);
		central.a2_ave /= ((double) central.N_bin);
		central.ma_ave /= ((double) central.N_bin);
	} else {
		central.m_bin_ave = 0.0;
		central.v_bin_rms = 0.0;
		central.a_ave = 0.0;
		central.a2_ave = 0.0;
		central.ma_ave = 0.0;
	}

	/* set global variables that are used throughout the code */
	rho_core_single = central.rho_sin;
	rho_core_bin = central.rho_bin;
	
	free(rhoj);
}

#ifdef USE_MPI
void mpi_clusdyn_calculate(void)
{
	strcpy(funcName, __FUNCTION__);
	double m=0.0;
	long k=1;
	
	while (m < 0.5 * Mtotal) {
		m += star_m[k] / clus.N_STAR;
		k++;
	}
	clusdyn.rh = star_r[k];
}
#endif

/* calculate cluster dynamical quantities */
void clusdyn_calculate(void)
{
	strcpy(funcName, __FUNCTION__);
	double m=0.0;
	long k=1;
	
	while (m < 0.5 * Mtotal) {
		m += star[k].m / clus.N_STAR;
		k++;
	}
	clusdyn.rh = star[k].r;
}

/* This is probably the nth incarnation of vr^2 (aka Q). It differs in as much
 * as it is not a macro but a real function with input parameters all being
 * long double. This has the effect that, on IA-32 systems, vr^2 has always a
 * predictable precision, even when the 387 fpu is used. The reason is simply
 * that vr^2 is now forced to be calculated with extended double precision
 * matching the precision that is used internally by the 387 fpu. This way,
 * possibly mixed type (double, extended) expressions are avoided when 387
 * instructions are generated by the compiler (This problem does not exist with
 * -msse2 -mfpmath=sse). Previously, the results depended on the order of
 *  evaluation of the subexpressions, which changes with the optimization level
 *  and surrounding code, sometimes varying by much more than DBL_EPSILON. It
 *  is clear that there is little one can do about this problem as long as vr^2
 *  is a macro, so I decided to write vr^2 as a function with all variables
 *  declared long double. I also declared it as inlined (see cmc.h), so there
 *  should be no performance impact.
 *
 * -- Stefan, 10/01/07
 */
inline double function_q(long j, long double r, long double pot, long double E, long double J) { 
  double res; 
  long double Jr, phis;

  Jr= SQR((J)/(r));
  phis= PHI_S(r, j);
  res= (2.0 * ((E) - (pot + phis)) - Jr);
  //if (j==3265 && r>1e6) printf("Jr=%Lg, phis= %Lg, pot=%Lg, r=%Lg, E=%Lg\n", Jr, phis, pot, r, E);
  return (res);
};

void timeStart()
{
/*
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	startTime = MPI_Wtime();
#else
	startTime = clock();
#endif
*/
}

void timeEnd(char* fileName, char *funcName, double *tTime)
{
/*
	double temp;
	static double totTime = 0.0;
	FILE *file;

#ifdef USE_MPI
	endTime = MPI_Wtime();
	temp = endTime - startTime;
	MPI_Reduce(&temp, tTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else
	endTime = clock();
	*tTime = ((double) (endTime - startTime)) / CLOCKS_PER_SEC;
#endif

#ifdef USE_MPI
	if(myid==0)
#endif
	{
		totTime += *tTime;
		file = fopen(fileName,"a");
		fprintf(file, "%-5.8lf\t\t%-25s\t\t\t%5.8lf\n", *tTime, funcName, totTime);
		fclose(file);
	}
*/
}

void timeStart2(double *st)
{
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	*st = MPI_Wtime();
#else
	*st = clock();
#endif
}

void timeEnd2(char* fileName, char *funcName, double *st, double *end, double *tot)
{
	double temp;
	//FILE *file;
#ifdef USE_MPI
	*end = MPI_Wtime();
	temp = *end - *st;
   MPI_Reduce(&temp, tot, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else
	*end = clock();
	*tot = ((double) (*end - *st)) / CLOCKS_PER_SEC;
#endif

#ifdef USE_MPI
	if(myid==0)
#endif
	{
		//file = fopen(fileName,"a");
		//fprintf(file, "%-5.8lf\t\t%-25s\t\t\t%5.8lf\n", *tot, funcName, 0.0 );
		printf("Total time = %g\n", *tot);
		//fclose(file);
	}
}

void create_timing_files()
{
	strcpy(fileTime, "mpi_time.dat");
	FILE *file = fopen(fileTime, "w");
	fprintf(file, "Time(s)\t\t\tFunction\t\t\t\t\t\t\t\t\tTotalTime\n");
	fclose(file);
}

void set_global_vars1()
{
	quiet = 0;
	debug = 0;
	initial_total_mass = 1.0;
	newstarid = 0;
	cenma.E = 0.0;
}

void set_global_vars2()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

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

	timeEnd(fileTime, funcName, &timeTotLoc);
}


void calc_sigma_new()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	//MPI2: mpi_calc_sigma_r(): Tests for precision: Biggest errors are ~ 1e-13. Might be because of catastrphic cancellation/rounding errors?. This might be due to the already imprecise but faster serial code . In a sense, the MPI version might be more precise, because for every processor, actual average is performed for the 1st local star which means errors dont carry over from the previous set of stars which are handled by another processor.
#ifdef USE_MPI
	mpi_calc_sigma_r();
#else
	calc_sigma_r(); //requires data from some neighbouring particles. must be handled during data read. look inside function for more comments
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void calc_central_new()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	//Step 2: stay on root node
	//MPI2: Split into 2 functions: part 1 can be parallelized after making m and r arrays global. Part 2 has to be done on root node.
#ifdef USE_MPI
	//MPI2: Tested! Errors of 1e-14.
	mpi_central_calculate();
#else
	central_calculate();
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void bin_vars_calculate()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	int i, j;
	//MPI2: Binaries; Ignore.
	M_b = 0.0;
	E_b = 0.0;
#ifdef USE_MPI
	for (i=mpiBegin; i<=mpiEnd; i++)
#else
	for (i=1; i<=clus.N_STAR; i++)
#endif
	{
		j = star[i].binind;
		if (j && binary[j].inuse) {
			M_b += star[i].m;
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void calc_potential_new()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	//MPI2: Do on root node with global phi array.
	//MPI2: Tested. Relative errors are 0. Perfect :)
#ifdef USE_MPI
	//MPI2: Does it need to be done on root, or can it be done on all nodes for all stars? Think.
	if(myid==0) 
		mpi_potential_calculate();

	MPI_Bcast(&clus.N_MAX, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(star_phi, clus.N_MAX+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Mtotal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Rtidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cenma.m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//MPI2: Calculating indices which will be used in all loops till beginning of the main loop. The value 20 depends on the p value used in calc_sigma_new()
   mpiFindIndicesCustom( clus.N_MAX, 20, myid, &mpiBegin, &mpiEnd );
#else
	potential_calculate();
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void calc_potential_new2()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	//MPI2: Do on root node with global phi array.
	//MPI2: Tested. Relative errors are 0. Perfect :)
#ifdef USE_MPI
	//MPI2: Does it need to be done on root, or can it be done on all nodes for all stars? Think.
	if(myid==0) 
#endif
		potential_calculate();

#ifdef USE_MPI
	MPI_Bcast(&clus.N_MAX, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Mtotal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Rtidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cenma.m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//MPI2: Calculating new indices which will be used in all loops till end of next timestep (qsorts).
	mpiFindIndicesCustom( clus.N_MAX, 20, myid, &mpiBegin, &mpiEnd );
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void compute_energy_new()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	//MPI2: Tested! No errors.
#ifdef USE_MPI
	mpi_ComputeEnergy();
#else
	//MPI2: see inline for comments
	ComputeEnergy();
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void set_energy_vars()
{
	/* Noting the total initial energy, in order to set termination energy. */
	Etotal.ini = Etotal.tot;

	Etotal.New = 0.0;
	Eescaped = 0.0;
	Jescaped = 0.0;
	Eintescaped = 0.0;
	Ebescaped = 0.0;
	Eoops = 0.0;
}

void reset_interaction_flags()
{
	int i;
#ifdef USE_MPI
		for (i=mpiBegin; i<=mpiEnd; i++)
#else
		for (i = 1; i <= clus.N_MAX; i++) 
#endif
			/* reset interacted flag */
			star[i].interacted = 0;
}

void calc_clusdyn_new()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	/* calculate dynamical quantities */
	//Step 2: done on root node
	//MPI2: Tested! No Errors.
#ifdef USE_MPI
	if(myid==0)
		mpi_clusdyn_calculate(); //parallel reduction

	//MPI2: broadcast clusdyn struct	
	MPI_Bcast(&clusdyn, sizeof(clusdyn_struct_t), MPI_BYTE, 0, MPI_COMM_WORLD);
#else
	clusdyn_calculate(); //parallel reduction
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void calc_timestep(gsl_rng *rng)
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	/* Get new time step */
	//for the next step, only simul_relax needs to be done in parallel, and DTrel needs to be broadcasted. rest can be done on each node.
	//Step 2: Do simul_relax() on all procs and others only on root. then broadcast Dt
	Dt = GetTimeStep(rng); //reduction again. Timestep needs to be communicated to all procs.

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
#endif

	TotalTime += Dt;

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void energy_conservation1()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	/* some numbers necessary to implement Stodolkiewicz's
	 * energy conservation scheme */
	int i;
#ifdef USE_MPI 
	//MPI2: Only running till N_MAX for now since no new stars are created, later new loop has to be introduced from N_MAX+1 to N_MAX_NEW.
	for (i=mpiBegin; i<=mpiEnd; i++)
#else
	for (i = 1; i <= clus.N_MAX_NEW; i++)
#endif
	{
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
	timeEnd(fileTime, funcName, &timeTotLoc);
}

void toy_rejuvenation()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	/*Sourav: checking all stars for their possible extinction from old age*/
	//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
	int i;
	DMrejuv = 0.0;
	if (STAR_AGING_SCHEME > 0) {
#ifdef USE_MPI 
		for (i=mpiBegin; i<=mpiEnd; i++)
#else
		//MPI2: Why is this only till N_MAX?
		for (i=1; i<=clus.N_MAX; i++)
#endif
			remove_old_star(TotalTime, i);
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void tidally_strip_stars1()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

#ifdef USE_MPI
		//MPI2: The 2nd part of tidally_strip_stars() has been moved just before sort.
		OldTidalMassLoss = TidalMassLoss;
		/******************************/
		/* get new particle positions */
		/******************************/
		max_r = get_positions();

		double temp = 0.0;
		MPI_Reduce(&Eescaped, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) Eescaped = temp;
		MPI_Bcast(&Eescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Reduce(&Jescaped, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) Jescaped = temp;
		MPI_Bcast(&Jescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Reduce(&Eintescaped, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) Eintescaped = temp;
		MPI_Bcast(&Eintescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Reduce(&Ebescaped, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) Ebescaped = temp;
		MPI_Bcast(&Ebescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Reduce(&TidalMassLoss, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) TidalMassLoss = temp;
		MPI_Bcast(&TidalMassLoss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Reduce(&Etidal, &temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);		
		if(myid==0) Etidal = temp;
		MPI_Bcast(&Etidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
		/* this calls get_positions() */
		tidally_strip_stars();
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void energy_conservation2()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	int i;
#ifdef USE_MPI 
	//MPI2: Should be changed into 2 loops later to include new stars.
	for (i=mpiBegin; i<=mpiEnd; i++) 
#else
	/* more numbers necessary to implement Stodolkiewicz's
	 * energy conservation scheme */
	for (i = 1; i <= clus.N_MAX_NEW; i++) 
#endif
	/* the following cannot be calculated after sorting 
	 * and calling potential_calculate() */
	{
#ifdef USE_MPI 
		star[i].Uoldrnew = potential(star[i].rnew) + MPI_PHI_S(star[i].rnew, i);
#else
		star[i].Uoldrnew = potential(star[i].rnew) + PHI_S(star[i].rnew, i);
#endif
	}

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void tidally_strip_stars2(void)
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

#ifdef USE_MPI
	if(myid==0)
	{
		double phi_rtidal, phi_zero, gierszalpha;
		long i, j, k;
		k=0;

		j = 0;
		Etidal = 0.0;
		OldTidalMassLoss = TidalMassLoss;
		DTidalMassLoss = TidalMassLoss - OldTidalMassLoss;

		gprintf("tidally_strip_stars(): iteration %ld: OldTidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
				j, OldTidalMassLoss, DTidalMassLoss);
		fprintf(logfile, "tidally_strip_stars(): iteration %ld: OldTidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
				j, OldTidalMassLoss, DTidalMassLoss);


		/* Iterate the removal of tidally stripped stars 
		 * by reducing Rtidal */
		do {
			Rtidal = orbit_r * pow(Mtotal 
					- (TidalMassLoss - OldTidalMassLoss), 1.0/3.0);
			phi_rtidal = potential(Rtidal);
			phi_zero = potential(0.0);
			DTidalMassLoss = 0.0;

			/* XXX maybe we should use clus.N_MAX_NEW below?? */
			//MPI2: Only running till N_MAX for now since no new stars are created, later new loop has to be introduced from N_MAX+1 to N_MAX_NEW.
			for (i = 1; i <= clus.N_MAX; i++) 
			{
				if (TIDAL_TREATMENT == 0){
					/*radial cut off criteria*/

					if (star[i].r_apo > Rtidal && star[i].rnew < 1000000) { 
						dprintf("tidally stripping star with r_apo > Rtidal: i=%ld id=%ld m=%g E=%g binind=%ld\n", i, star[i].id, star[i].m, star[i].E, star[i].binind);
						star[i].rnew = SF_INFINITY;	/* tidally stripped star */
						star[i].vrnew = 0.0;
						star[i].vtnew = 0.0;

						Eescaped += star[i].E * star[i].m / clus.N_STAR;
						Jescaped += star[i].J * star[i].m / clus.N_STAR;

						if (star[i].binind == 0) {
							Eintescaped += star[i].Eint;
						} else {
							Ebescaped += -(binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
								(2.0 * binary[star[i].binind].a);
							Eintescaped += binary[star[i].binind].Eint1 + binary[star[i].binind].Eint2;
						}

						DTidalMassLoss += star[i].m / clus.N_STAR;
						Etidal += star[i].E * star[i].m / clus.N_STAR;

						/* logging */
						fprintf(escfile,
								"%ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld ",
								tcount, TotalTime, star[i].m,
								star[i].r, star[i].vr, star[i].vt, star[i].r_peri,
								star[i].r_apo, Rtidal, phi_rtidal, phi_zero, star[i].E, star[i].J, star[i].id);

						if (star[i].binind) {
							k = star[i].binind;
							fprintf(escfile, "1 %.8g %.8g %ld %ld %.8g %.8g ", 
									binary[k].m1 * (units.m / clus.N_STAR) / MSUN, 
									binary[k].m2 * (units.m / clus.N_STAR) / MSUN, 
									binary[k].id1, binary[k].id2,
									binary[k].a * units.l / AU, binary[k].e);
						} else {
							fprintf(escfile, "0 0 0 0 0 0 0 ");	
						}

						if (star[i].binind == 0) {
							fprintf(escfile, "%d na na ", 
									star[i].se_k);
						} else {
							fprintf(escfile, "na %d %d",
									binary[k].bse_kw[0], binary[k].bse_kw[1]);
						}
						fprintf (escfile, "\n");



						/* perhaps this will fix the problem wherein stars are ejected (and counted)
							multiple times */
						dprintf ("before SE: id=%ld k=%ld kw=%d m=%g mt=%g R=%g L=%g mc=%g rc=%g menv=%g renv=%g ospin=%g epoch=%g tms=%g tphys=%g phi=%g r=%g\n",
								star[i].id,i,star[i].se_k,star[i].se_mass,star[i].se_mt,star[i].se_radius,star[i].se_lum,star[i].se_mc,star[i].se_rc,
								star[i].se_menv,star[i].se_renv,star[i].se_ospin,star[i].se_epoch,star[i].se_tms,star[i].se_tphys,star[i].phi, star[i].r);
						destroy_obj(i);
						if (Etotal.K + Etotal.P - Etidal >= 0)
							break;

					}
				}


				else if (TIDAL_TREATMENT == 1){
					/* DEBUG: Now using Giersz prescription for tidal stripping 
						(Giersz, Heggie, & Hurley 2008; arXiv:0801.3709).
						Note that this alpha factor behaves strangely for small N (N<~10^3) */

					gierszalpha = 1.5 - 3.0 * pow(log(GAMMA * ((double) clus.N_STAR)) / ((double) clus.N_STAR), 0.25);
					if (star[i].E > gierszalpha * phi_rtidal && star[i].rnew < 1000000) {
						dprintf("tidally stripping star with E > phi rtidal: i=%ld id=%ld m=%g E=%g binind=%ld\n", i, star[i].id, star[i].m, star[i].E, star[i].binind); 
						star[i].rnew = SF_INFINITY;	/* tidally stripped star */
						star[i].vrnew = 0.0;
						star[i].vtnew = 0.0;
						Eescaped += star[i].E * star[i].m / clus.N_STAR;
						Jescaped += star[i].J * star[i].m / clus.N_STAR;
						if (star[i].binind == 0) {
							Eintescaped += star[i].Eint;
						} else {
							Ebescaped += -(binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
								(2.0 * binary[star[i].binind].a);
							Eintescaped += binary[star[i].binind].Eint1 + binary[star[i].binind].Eint2;
						}

						DTidalMassLoss += star[i].m / clus.N_STAR;
						Etidal += star[i].E * star[i].m / clus.N_STAR;

						/* logging */
						fprintf(escfile,
								"%ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld ",
								tcount, TotalTime, star[i].m,
								star[i].r, star[i].vr, star[i].vt, star[i].r_peri,
								star[i].r_apo, Rtidal, phi_rtidal, phi_zero, star[i].E, star[i].J, star[i].id);

						if (star[i].binind) {
							k = star[i].binind;
							fprintf(escfile, "1 %.8g %.8g %ld %ld %.8g %.8g ", 
									binary[k].m1 * (units.m / clus.N_STAR) / MSUN, 
									binary[k].m2 * (units.m / clus.N_STAR) / MSUN, 
									binary[k].id1, binary[k].id2,
									binary[k].a * units.l / AU, binary[k].e);
						} else {
							fprintf(escfile, "0 0 0 0 0 0 0 ");	
						}

						if (star[i].binind == 0) {
							fprintf(escfile, "%d na na ", 
									star[i].se_k);
						} else {
							fprintf(escfile, "na %d %d",
									binary[k].bse_kw[0], binary[k].bse_kw[1]);
						}
						fprintf (escfile, "\n");

						/* perhaps this will fix the problem wherein stars are ejected (and counted)
							multiple times */
						destroy_obj(i);

						if (Etotal.K + Etotal.P - Etidal >= 0)
							break;

					}
				}

			}
			j++;
			TidalMassLoss = TidalMassLoss + DTidalMassLoss;
			gprintf("tidally_strip_stars(): iteration %ld: TidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
					j, TidalMassLoss, DTidalMassLoss);
			fprintf(logfile, "tidally_strip_stars(): iteration %ld: TidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
					j, TidalMassLoss, DTidalMassLoss);
		} while (DTidalMassLoss > 0 && (Etotal.K + Etotal.P - Etidal) < 0);

	}
	MPI_Bcast(&Eescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Jescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Eintescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Ebescaped, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&DTidalMassLoss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Etidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Rtidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void pre_sort_comm()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	int i;
#ifdef USE_MPI
	MPI_Status stat;
	//MPI2: Collecting the r and m arrays into the original star structure for sorting.
	//MPI2: Only running till N_MAX for now as no new stars are created,later this has to changed to include the new stars somehow. Note that N_MAX_NEW will be different for different processors.
	mpiFindDispAndLenCustom( clus.N_MAX, 20, mpiDisp, mpiLen );

	for(i=0;i<procs;i++)
		mpiLen[i] *= sizeof(star_t); 

	//MPI2: Only running till N_MAX for now as no new stars are created,later this has to changed to include the new stars somehow. Note that N_MAX_NEW will be different for different processors.
	for(i=mpiBegin; i<=mpiEnd; i++) {
		star[i].r = star_r[i];
		star[i].m = star_m[i];
	}

	//MPI2: To be refactored into separate function later.
	if(myid!=0)
		MPI_Send(&star[mpiDisp[myid]], mpiLen[myid], MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	else
		for(i=1;i<procs;i++)
			MPI_Recv(&star[mpiDisp[i]], mpiLen[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &stat);
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void qsorts_new(void)
{
	strcpy(funcName, "qsorts_new");
	static double timeTotLoc;
	timeStart();

#ifdef USE_MPI
	if(myid==0)
#endif
		/* Sorting stars by radius. The 0th star at radius 0 
			and (N_STAR+1)th star at SF_INFINITY are already set earlier.
		 */
		//MPI2: Only running till N_MAX for now since no new stars are created, later has to be changed to N_MAX_NEW.
		//MPI2: Changin to N_MAX+1, later change to N_MAX_NEW+1
		qsorts(star+1,clus.N_MAX);

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void post_sort_comm()
{
	strcpy(funcName, __FUNCTION__);
	static double timeTotLoc;
	timeStart();

	int i;
#ifdef USE_MPI
	MPI_Status stat;
	mpiFindDispAndLenCustom( clus.N_MAX, 20, mpiDisp, mpiLen );

	for(i=0;i<procs;i++)
		mpiLen[i] *= sizeof(star_t); 

	//MPI2: To be refactored into separate function later.
	if(myid==0)
		for(i=1;i<procs;i++)
			MPI_Send(&star[mpiDisp[i]], mpiLen[i], MPI_BYTE, i, 0, MPI_COMM_WORLD);
	else
		MPI_Recv(&star[mpiDisp[myid]], mpiLen[myid], MPI_BYTE, 0, 0, MPI_COMM_WORLD, &stat);

	if(myid==0)
		for(i=0; i<=clus.N_MAX+1; i++) {
			star_r[i] = star[i].r;
			star_m[i] = star[i].m;
			star_phi[i] = star[i].phi;
		}

	MPI_Bcast(star_m, clus.N_MAX+2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(star_r, clus.N_MAX+2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(star_phi, clus.N_MAX+2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	timeEnd(fileTime, funcName, &timeTotLoc);
}

void findIndices( long N, int blkSize, int i, int* begin, int* end )
{
	long chunkSize =  ( N / procs ) * blkSize;
   if ( i < N % procs )
   {
      *begin = i * chunkSize + i * blkSize + 1; //+1 since for loops go from 1 to N
      *end = *begin + chunkSize + 1 * blkSize - 1;
   } else {
      *begin = i * chunkSize + ( N % procs ) * blkSize + 1; //+1 since for loops go from 1 to N
      *end =  *begin + chunkSize - 1;
   }
}

void findLimits( long N, int blkSize )
{
	int i, blocks, remain;

   blocks = N / blkSize;
   remain = N % blkSize;

	for( i = 0; i < procs; i++ )
		findIndices( blocks, blkSize, i, &Start[i], &End[i] );

	End[procs-1] += remain;
}

int findProcForIndex( int j )
{
	int i;
	for( i = 0; i < procs; i++ )
		if( j >= Start[i] && j <= End[i] )
			break;

	return i;
}

void set_rng_states()
{
	int i;
#ifdef USE_MPI
	curr_st = (struct rng_t113_state*) malloc(sizeof(struct rng_t113_state));
	reset_rng_t113_new(IDUM, curr_st);

	for(i = 0; i < myid; i++)
		*curr_st = rng_t113_jump( *curr_st , JPoly_2_20);
#else
	st = (struct rng_t113_state*) malloc(procs * sizeof(struct rng_t113_state));
	reset_rng_t113_new(IDUM, &st[0]);

	for(i = 1; i < procs; i++)
		st[i] = rng_t113_jump( st[i-1] , JPoly_2_20);
#endif
}
