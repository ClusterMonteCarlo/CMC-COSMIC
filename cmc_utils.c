/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/time.h>
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
		exit_cleanly(-2, __FUNCTION__);
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
             i = FindZero_r(kmin, kmax, r);
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
		exit_cleanly(-2, __FUNCTION__);
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

void exit_cleanly_old(int signal)
{

	close_buffers();
	free_arrays();

#ifdef USE_MPI
	//mpi_files_merge();
	//OPT: Another way would be to use assert at all error checks
	MPI_Abort(MPI_COMM_WORLD, signal);
#endif

	exit(signal);
}
/* close buffers, then exit */
void exit_cleanly(int signal, const char* fn)
{

	close_buffers();
	free_arrays();

#ifdef USE_MPI
	//mpi_files_merge();
	//OPT: Another way would be to use assert at all error checks
	printf("------->MPI_Abort called in %s by proc %d\n", fn, myid);
	MPI_Abort(MPI_COMM_WORLD, signal);
#endif

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
	free(new_size); free(disp); free(len);
#endif
}

/* GSL error handler */
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno)
{
	fprintf(stderr, "gsl: %s:%d: ERROR: %s\n", file, line, reason);
	exit_cleanly(gsl_errno, __FUNCTION__);
}

#ifdef USE_MPI
void mpi_set_velocities3(void){
	/* set velocities a la Stodolkiewicz to be able to conserve energy */
	double vold2, vnew2, Unewrold, Unewrnew;
	double Eexcess, exc_ratio, Eexcess_prev;
	double q=0.5; /* q=0.5 -> Stodolkiewicz, q=0 -> Delta U all from rnew */
	double alpha;
	long i;
	int g_i;
	double m, E_dump, E_dump_capacity, Eexcess_check, v2, v2_new;
	double E_dump_factor = 0.8; //Percentage of excess energy to be dumped on each star.
	char *vnew2_arr = (char *) malloc ((clus.N_MAX_NEW+1) * sizeof(char));
	MPI_Status stat;

	Eexcess = 0.0; Eexcess_prev = 0.0;
	for (i = 1; i <= clus.N_MAX_NEW; i++) {
		g_i = get_global_idx(i);
		m = star_m[g_i];

		/* modify velocities of stars that have only undergone relaxation */
		if (star[i].interacted == 0) {
			Unewrold = potential(star[i].rOld) + MPI_PHI_S(star[i].rOld, g_i);
			Unewrnew = star_phi[g_i] + MPI_PHI_S(star_r[g_i], g_i);

			vold2 = star[i].vtold*star[i].vtold + 
				star[i].vrold*star[i].vrold;

			/* predict new velocity */
			vnew2 = vold2 + 2.0*(1.0-q)*(star[i].Uoldrold - star[i].Uoldrnew)
				+ 2.0*q*(Unewrold - Unewrnew);
			vnew2_arr[i] = (vnew2 > 0.0) ? 1 : 0;

			/* new velocity can be unphysical, so just use value predicted by old potential
				(this is already set in .vr and .vt) */
			if (vnew2 <= 0.0) {
				Eexcess += 0.5*(sqr(star[i].vr)+sqr(star[i].vt)-vnew2)*m;
			} else {
				/* scale velocity, preserving v_t/v_r */
				alpha = sqrt(vnew2/(sqr(star[i].vr)+sqr(star[i].vt)));
				star[i].vr *= alpha;
				star[i].vt *= alpha;
				v2 = sqr(star[i].vr)+sqr(star[i].vt);
				star[i].vrnew = star[i].vr;
				star[i].vtnew = star[i].vt;

				/* if there is excess energy added, try to remove at 
					least part of it from this star */
				E_dump_capacity =  E_dump_factor * 0.5 * m * v2;
				E_dump = MIN( E_dump_capacity, Eexcess );
				if(Eexcess > 0){
					exc_ratio = sqrt( (v2 - 2 * E_dump / m) / v2 );
					star[i].vr *= exc_ratio;
					star[i].vt *= exc_ratio;
					Eexcess -= E_dump;
				}
			}
		}
	}


	if(Eexcess > 0)
		printf("id=%d WARNING!!!!! -------> Excess = %g\n", myid, Eexcess);

	Eexcess_prev = 0.0;
	do {
		if(myid != procs-1)
			MPI_Send(&Eexcess, 1, MPI_DOUBLE, ( myid + 1 ), 0, MPI_COMM_WORLD);
		if(myid != 0)
			MPI_Recv(&Eexcess_prev, 1, MPI_DOUBLE, ( myid - 1), 0, MPI_COMM_WORLD, &stat);

		Eexcess = Eexcess_prev;
		for (i = 1; i <= clus.N_MAX_NEW; i++) {
			g_i = get_global_idx(i);
			m = star_m[g_i];

			v2_new = sqr(star[i].vrnew)+sqr(star[i].vtnew);
			v2 = sqr(star[i].vr)+sqr(star[i].vt);
			if (star[i].interacted == 0) {
				if (vnew2_arr[i] ==1) {
					E_dump_capacity = v2*m*0.5 - E_dump_factor * 0.5 * m * v2_new;
					E_dump = MIN( E_dump_capacity, Eexcess );
					if(Eexcess > 0) {
						exc_ratio = sqrt( (v2 - 2 * E_dump / m) / v2 );
						star[i].vr *= exc_ratio;
						star[i].vt *= exc_ratio;
						Eexcess -= E_dump;
					}
				}
			}
		}
		MPI_Allreduce(&Eexcess, &Eexcess_check, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	} while(Eexcess_check > 0);

	free(vnew2_arr);
	/* keep track of the energy that's vanishing due to our negligence */
	Eoops += -Eexcess * madhoc;
}
#endif

void set_velocities3(void){
	/* set velocities a la Stodolkiewicz to be able to conserve energy */
	double vold2, vnew2, Unewrold, Unewrnew;
	double Eexcess, exc_ratio;
	double q=0.5; /* q=0.5 -> Stodolkiewicz, q=0 -> Delta U all from rnew */
	double alpha;
	long i;
	int g_i;
	double m, v2, E_dump, E_dump_capacity;
	double E_dump_factor = 0.8; //Percentage of excess energy to be dumped on each star.

	Eexcess = 0.0;
	for (i = 1; i <= clus.N_MAX; i++) {
		m = star[i].m;

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
				Eexcess += 0.5*(sqr(star[i].vr)+sqr(star[i].vt)-vnew2)*m;
			} else {
				/* scale velocity, preserving v_t/v_r */
				alpha = sqrt(vnew2/(sqr(star[i].vr)+sqr(star[i].vt)));
				star[i].vr *= alpha;
				star[i].vt *= alpha;
				v2 = sqr(star[i].vr)+sqr(star[i].vt);

				/* if there is excess energy added, try to remove at 
					least part of it from this star */
				// MPI3: Check how the distribution of E_dump is but printing out E_dump to a log file.
				E_dump_capacity =  E_dump_factor * 0.5 * m * v2;
				E_dump = MIN( E_dump_capacity, Eexcess );
				if(Eexcess > 0){
					exc_ratio = sqrt( (v2 - 2 * E_dump / m) / v2 );
					star[i].vr *= exc_ratio;
					star[i].vt *= exc_ratio;
					Eexcess -= E_dump;
				}
			}
		}
	}

	/* keep track of the energy that's vanishing due to our negligence */
	Eoops += -Eexcess * madhoc;
}

void energy_conservation3(void)
{
#ifdef USE_MPI
	mpi_set_velocities3();
#else
	set_velocities3();
#endif
}
/* computes intermediate energies, and transfers "new" dynamical params to the standard variables */
void ComputeIntermediateEnergy(void)
{
	int j = 1, g_j; 
	/* compute intermediate energies for stars due to change in pot */ 
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
		/* but do only for NON-Escaped stars */
		if (star[j].rnew < 1.0e6) {
#ifdef USE_MPI
			g_j = get_global_idx(j);
			star[j].EI = sqr(star[j].vr) + sqr(star[j].vt) + star_phi[g_j] - potential(star[j].rnew);
#else
			star[j].EI = sqr(star[j].vr) + sqr(star[j].vt) + star[j].phi - potential(star[j].rnew);
#endif
		}
	}

	/* Transferring new positions to .r, .vr, and .vt from .rnew, .vrnew, and .vtnew */
	for (j = 1; j <= clus.N_MAX_NEW; j++) {
#ifdef USE_MPI
		g_j = get_global_idx(j);
		star[j].rOld = star_r[g_j];
		star[j].r = star[j].rnew;
		//star_r[j] = star[j].rnew;
		star[j].m = star_m[g_j];
#else
		star[j].rOld = star[j].r;
		star[j].r = star[j].rnew;
#endif
		star[j].vr = star[j].vrnew;
		star[j].vt = star[j].vtnew;
	}
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
/* Calculates E,J for every star. Also, calculates, global energy variabies. */
	double buf_reduce[5], phi0 = 0.0;
	int i, j=0;
	for(i=0; i<5; i++)
		buf_reduce[i] = 0.0;

	for (i=1; i<=mpiEnd-mpiBegin+1; i++) {
		j = get_global_idx(i);
		star[i].E = star_phi[j] + 0.5 * (sqr(star[i].vr) + sqr(star[i].vt));
		star[i].J = star_r[j] * star[i].vt;		
	}

	phi0 = star_phi[0];

	for (i=1; i<=mpiEnd-mpiBegin+1; i++) {
		j = get_global_idx(i);
		buf_reduce[1] += 0.5 * (sqr(star[i].vr) + sqr(star[i].vt)) * star_m[j] / clus.N_STAR;
		buf_reduce[2] += star_phi[j] * star_m[j] / clus.N_STAR;
		buf_reduce[2] += phi0 * cenma.m*madhoc/ clus.N_STAR;

	//MPI3:Ignoring binaries for now.
		if (star[i].binind == 0) {
			buf_reduce[3] += star[i].Eint;
		} else if (binary[star[i].binind].inuse) {
			buf_reduce[4] += -(binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
				(2.0 * binary[star[i].binind].a);
			buf_reduce[3] += binary[star[i].binind].Eint1 + binary[star[i].binind].Eint2;
		}
	}
	
	buf_reduce[2] *= 0.5;
	buf_reduce[0] = buf_reduce[1] + buf_reduce[2] + buf_reduce[3] + buf_reduce[4];

	MPI_Allreduce(MPI_IN_PLACE, buf_reduce, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		

	Etotal.tot = buf_reduce[0];
	Etotal.K = buf_reduce[1];
	Etotal.P = buf_reduce[2];
	Etotal.Eint = buf_reduce[3];
	Etotal.Eb = buf_reduce[4];

	Etotal.tot += cenma.E + Eescaped + Ebescaped + Eintescaped;

	double temp = 0.0;
	MPI_Status stat;
//if(myid==0)
//	Etotal.K = buf_reduce[1];	

	//MPI2: Avoiding reduce to improve accuracy.
/*
	if(myid!=0)
		MPI_Send(&buf_reduce[1], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	else
		for(i=1;i<procs;i++)
		{
			MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
			Etotal.K += temp;
		}
*/
	//if(myid==0)
	//	printf("Etotal.tot=%.18g Etotal.K=%.18g Etotal.P=%.18g Etotal.Eint=%.18g Etotal.Eb=%.18g cenma.E=%.18g Eescaped=%.18g Ebescaped=%.18g Eintescaped=%.18g\n",Etotal.tot, Etotal.K, Etotal.P, Etotal.Eint, Etotal.Eb, cenma.E, Eescaped, Ebescaped, Eintescaped);
}
#endif

void ComputeEnergy(void)
{
	long i;
	
	Etotal.tot = 0.0;
	Etotal.K = 0.0;
	Etotal.P = 0.0;
	Etotal.Eint = 0.0;
	Etotal.Eb = 0.0;

	for (i=1; i<=clus.N_MAX; i++) {
		star[i].E = star[i].phi + 0.5 * (sqr(star[i].vr) + sqr(star[i].vt));
		star[i].J = star[i].r * star[i].vt;
		
		Etotal.K += 0.5 * (sqr(star[i].vr) + sqr(star[i].vt)) * star[i].m / clus.N_STAR;
		Etotal.P += star[i].phi * star[i].m / clus.N_STAR;
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

	//printf("Etotal.tot=%.18g Etotal.K=%.18g Etotal.P=%.18g Etotal.Eint=%.18g Etotal.Eb=%.18g cenma.E=%.18g Eescaped=%.18g Ebescaped=%.18g Eintescaped=%.18g\n", Etotal.tot, Etotal.K, Etotal.P, Etotal.Eint, Etotal.Eb, cenma.E, Eescaped, Ebescaped, Eintescaped);
}

#ifdef USE_MPI
long mpi_potential_calculate(void) {
	long k;
	double mprev;

	/* count up all the mass and set N_MAX */
	k = 1;
	mprev = 0.0;
	
	//while (star_r[k] < SF_INFINITY && k <= clus.N_STAR_NEW) {
	//MPI2: Temporarily changing this to N_MAX.
	//MPI3: Do on all nodes. This can be totally parallelized, reduce k and mprev ath the end.
	while (star_r[k] < SF_INFINITY && k <= clus.N_MAX) {
		mprev += star_m[k];
		/* I guess NaNs do happen... */
		if(isnan(mprev)){
			eprintf("NaN (2) detected\n");
			exit_cleanly(-1, __FUNCTION__);
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

	//MPI3: What will happen to the last sentinel? Ask Stefan.
	/* zero boundary star first for safety */
	//zero_star(clus.N_MAX + 1);

	star_r[clus.N_MAX + 1] = SF_INFINITY;
	star_phi[clus.N_MAX + 1] = 0.0;

	//MPI3: In future consider parallelization. For now, keep on root.
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
		exit_cleanly(-1, __FUNCTION__);
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
			exit_cleanly(-1, __FUNCTION__);
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
		  exit_cleanly(-1, __FUNCTION__);
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
		exit_cleanly(-1, __FUNCTION__);
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
#ifdef USE_MPI
	if (r < star_r[1])
		return (star_phi[1]);
#else
	if (r < star[1].r)
		return (star[1].phi);
#endif

	i = FindZero_r(kmin, kmax, r);
	
#ifdef USE_MPI
	if(star_r[i] > r || star_r[i+1] < r){
#else
	if(star[i].r > r || star[i+1].r < r){
#endif
		eprintf("binary search (FindZero_r) failed!!\n");
		eprintf("pars: i=%ld, star[i].r = %e, star[i+1].r = %e, star[i+2].r = %e, star[i+3].r = %e, r = %e\n",
				i, star[i].r, star[i+1].r, star[i+2].r, star[i+3].r, r);
		eprintf("pars: star[i].m=%g star[i+1].m=%g star[i+2].m=%g star[i+3].m=%g\n",
			star[i].m, star[i+1].m, star[i+2].m, star[i+3].m);
		exit_cleanly(-2, __FUNCTION__);
	}

	/* Henon's method of computing the potential using star[].phi */ 
#ifdef USE_MPI
	if (i == 0){ /* I think this is impossible, due to early return earlier,
			    but I am keeping it. -- ato 23:17,  3 Jan 2005 (UTC) */
		henon = (star_phi[1]);
	} else {
		henon = (star_phi[i] + (star_phi[i + 1] - star_phi[i]) 
			 * (1.0/star_r[i] - 1.0/r) /
			 (1.0/star_r[i] - 1.0/star_r[i + 1]));
	}
#else
	if (i == 0){ /* I think this is impossible, due to early return earlier,
			    but I am keeping it. -- ato 23:17,  3 Jan 2005 (UTC) */
		henon = (star[1].phi);
	} else {
		henon = (star[i].phi + (star[i + 1].phi - star[i].phi) 
			 * (1.0/star[i].r - 1.0/r) /
			 (1.0/star[i].r - 1.0/star[i + 1].r));
	}
#endif

	return (henon);
}

long check_if_r_around_last_index(long last_index, double r) {
   long index_found, i;

   index_found= -1;
   for (i=-1; i<2; i++) {
      if (((last_index+i) >= 0) && ((last_index+i+1) <= clus.N_STAR)) {
         if ((star[last_index+i].r < r) && (star[last_index+i+1].r > r)) {
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
	exit_cleanly(1, __FUNCTION__);
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
	long i, k;
	
	/* update total number, mass, and binding energy of binaries in cluster */
	N_b = 0;
	M_b = 0.0;
	E_b = 0.0;
#ifdef USE_MPI
	for (i=1; i<=mpiEnd-mpiBegin+1; i++) 
#else
	for (i=1; i<=clus.N_MAX; i++)
#endif
	{
		k = star[i].binind;
		if (k != 0) {
			N_b++;
#ifdef USE_MPI
			M_b += star_m[get_global_idx(i)];
#else
			M_b += star[i].m;
#endif
			if (binary[k].inuse){
				E_b += (binary[k].m1/clus.N_STAR) * (binary[k].m2/clus.N_STAR) / (2.0 * binary[k].a);
			}
		}
	}
#ifdef USE_MPI
	MPI_Allreduce(MPI_IN_PLACE, &M_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	MPI_Allreduce(MPI_IN_PLACE, &E_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	MPI_Allreduce(MPI_IN_PLACE, &N_b, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);		
#endif
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

/* calculate central quantities */
void central_calculate(void)
{
	double m=0.0, *rhoj, mrho, Vrj, rhojsum, Msincentral, Mbincentral, Vcentral, rcentral;
	long J=6, i, j, jmin, jmax, nave, Ncentral;
	double rhoj2sum;

	/* average over all stars out to half-mass radius */
	nave = 1;
	while (m < 0.5 * Mtotal) {
#ifdef USE_MPI
		m += star_m[nave] / clus.N_STAR;
#else
		m += star[nave].m / clus.N_STAR;
#endif
		nave++;
	}
	
	/* DEBUG */
	/* fprintf(stderr, "nave=%ld\n", nave); */
	/* DEBUG */

	/* exit if not enough stars */
	if (clus.N_STAR <= 2*J || nave >= clus.N_STAR-6) {
		eprintf("clus.N_STAR <= 2*J || nave >= clus.N_STAR-6\n");
		exit_cleanly(-1, __FUNCTION__);
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
#ifdef USE_MPI
			mrho += star_m[j] * madhoc;
		}
		Vrj = 4.0/3.0 * PI * (fb_cub(star_r[jmax]) - fb_cub(star_r[jmin]));
		rhoj[i] = mrho / Vrj;
#else
			mrho += star[j].m * madhoc;
		}
		Vrj = 4.0/3.0 * PI * (fb_cub(star[jmax].r) - fb_cub(star[jmin].r));
		rhoj[i] = mrho / Vrj;
#endif
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
#ifdef USE_MPI
		if( i >= mpiBegin && i <= mpiEnd )
			central.v_rms += rhoj[i] * (sqr(star[get_local_idx(i)].vr) + sqr(star[get_local_idx(i)].vt));

		central.rc += rhoj[i] * star_r[i];
		rc_nb += sqr(rhoj[i] * star_r[i]);
		central.m_ave += rhoj[i] * star_m[i] * madhoc;
#else
      central.v_rms += rhoj[i] * (sqr(star[i].vr) + sqr(star[i].vt));
		central.rc += rhoj[i] * star[i].r;
		rc_nb += sqr(rhoj[i] * star[i].r);
		central.m_ave += rhoj[i] * star[i].m * madhoc;
#endif
	}

//MPI3: This reduce gives round-off errors which might affect the timestep mildly. Consult Stefan and see if it needs to be changed to be summed up in order.
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &central.v_rms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

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
		//MPI3: In the case when the central stars dont all lie in the root node, this just cant be done only by the root node. But assuming that case will never happen, we'll just throw an error for such cases.
		if(NUM_CENTRAL_STARS > End[0] - Start[0]) 
		{
			eprintf("Central stars dont fit in root node!\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		if(myid==0)
		{
			Ncentral++;

			j = get_global_idx(i);

		/* use only code units here, so always divide star[].m by clus.N_STAR */
#ifdef USE_MPI
			central.w2_ave += 2.0 * star_m[j] / ((double) clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt));
#else
			central.w2_ave += 2.0 * star[i].m / ((double) clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt));
#endif

			if (star[i].binind == 0) {
				central.N_sin++;
#ifdef USE_MPI
				Msincentral += star_m[j] / ((double) clus.N_STAR);
#else
				Msincentral += star[i].m / ((double) clus.N_STAR);
#endif
				central.v_sin_rms += sqr(star[i].vr) + sqr(star[i].vt);
				central.R2_ave += sqr(star[i].rad);
#ifdef USE_MPI
				central.mR_ave += star_m[j] / ((double) clus.N_STAR) * star[i].rad;
#else
				central.mR_ave += star[i].m / ((double) clus.N_STAR) * star[i].rad;
#endif
			} else {
				central.N_bin++;

#ifdef USE_MPI
				Mbincentral += star_m[j] / ((double) clus.N_STAR);
#else
				Mbincentral += star[i].m / ((double) clus.N_STAR);
#endif
				central.v_bin_rms += sqr(star[i].vr) + sqr(star[i].vt);
				central.a_ave += binary[star[i].binind].a;
				central.a2_ave += sqr(binary[star[i].binind].a);
#ifdef USE_MPI
				central.ma_ave += star_m[j] / ((double) clus.N_STAR) * binary[star[i].binind].a;
#else
				central.ma_ave += star[i].m / ((double) clus.N_STAR) * binary[star[i].binind].a;
#endif
			}
		}
	}

#ifdef USE_MPI
	//MPI3: Packing into array to optimize communication.
	double *buf_bcast_dbl = (double*) malloc(10 * sizeof(double));
	buf_bcast_dbl[0] = central.w2_ave;
	buf_bcast_dbl[1] = Msincentral;
	buf_bcast_dbl[2] = central.mR_ave;
	buf_bcast_dbl[3] = Mbincentral;
	buf_bcast_dbl[4] = central.ma_ave;
	buf_bcast_dbl[5] = central.a2_ave;
	buf_bcast_dbl[6] = central.a_ave;
	buf_bcast_dbl[7] = central.v_sin_rms;
	buf_bcast_dbl[8] = central.v_bin_rms;
	buf_bcast_dbl[9] = central.R2_ave;

	MPI_Bcast( buf_bcast_dbl, 10, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	central.w2_ave = buf_bcast_dbl[0];
	Msincentral = buf_bcast_dbl[1];
	central.mR_ave = buf_bcast_dbl[2];
	Mbincentral = buf_bcast_dbl[3];
	central.ma_ave = buf_bcast_dbl[4];
	central.a2_ave = buf_bcast_dbl[5];
	central.a_ave = buf_bcast_dbl[6];
	central.v_sin_rms = buf_bcast_dbl[7];
	central.v_bin_rms = buf_bcast_dbl[8];
	central.R2_ave = buf_bcast_dbl[9];
	free(buf_bcast_dbl);
	
	long *buf_bcast_long = (long*) malloc(3 * sizeof(long));
	buf_bcast_long[0] = Ncentral;
	buf_bcast_long[1] = central.N_sin;
	buf_bcast_long[2] = central.N_bin;

	MPI_Bcast( buf_bcast_long, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	Ncentral = buf_bcast_long[0];
	central.N_sin = buf_bcast_long[1];
	central.N_bin = buf_bcast_long[2];

	free(buf_bcast_long);
#endif

	/* object quantities */
#ifdef USE_MPI
	rcentral = star_r[Ncentral + 1];
#else
	rcentral = star[Ncentral + 1].r;
#endif
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

//printf("\n%d rho=%g v_rms=%g rc=%g m_ave=%g n=%g rc_s=%g N_sin=%ld N_bin=%ld n_sin=%g n_bin=%g\n rho_sin=%g rho_bin=%g m_sin_ave=%g m_bin_ave=%g v_sin=%g v_bin=%g w2=%g R2=%g mR=%g a=%g a2=%g ma=%g\n\n", myid, central.rho, central.v_rms, central.rc, central.m_ave, central.n, central.rc_spitzer, central.N_sin, central.N_bin, central.n_sin, central.n_bin, central.rho_sin, central.rho_bin, central.m_sin_ave, central.m_bin_ave, central.v_sin_rms, central.v_bin_rms, central.w2_ave, central.R2_ave, central.mR_ave, central.a_ave, central.a2_ave, central.ma_ave);

	/* set global variables that are used throughout the code */
	rho_core_single = central.rho_sin;
	rho_core_bin = central.rho_bin;
	
	free(rhoj);
}

double local_kT(long si, int p) {
  int simin, simax, j;
  double mave, kT;

  simin = si - p;
  simax = simin + (2 * p - 1);
  if (simin < 1) {
    simin = 1;
    simax = simin + (2 * p - 1);
  }
  mave= 0.;
  for (j=simin; j< simax; j++) {
    mave+= star[j].m *madhoc;
  }
  mave/= (double)(2 * p);

  kT= (1.0/3.0) * mave * sigma_array.sigma[si] * sigma_array.sigma[si];

  return(kT);
}

double core_kT(int Ncore, int p) {
  int i;
  double kT;

  kT= 0.;
  for (i=1; i<Ncore+1; i++) {
    kT+= local_kT(i, p);
  }
  kT/= Ncore;

  return(kT);
}

central_t
central_hard_binary(double ktmin, central_t old_cent) {
  int i, Ncentral; 
  double Rcentral, Vcentral, kTcore, Mbincentral;
  central_t cent;
  
  cent= old_cent;

  cent.a_ave=0.;
  cent.a2_ave=0.;
  cent.v_bin_rms=0.;
  cent.ma_ave=0.;
  cent.m_bin_ave=0.;
  cent.n_bin=0.;
  cent.rho_bin=0.;
  cent.N_bin=0;

  Ncentral= MIN(NUM_CENTRAL_STARS, clus.N_MAX);
  Rcentral= star[Ncentral+1].r;
  Vcentral= 4./3.*PI*cub(Rcentral);
  kTcore= core_kT(Ncentral, AVEKERNEL);
  printf("kTcore is %g\n", kTcore);

  Mbincentral=0.;
  for (i=1; i<Ncentral+1; i++) {
    double Eb, a, m1, m2;

    if (!(star[i].binind>0)) continue;
    
    a= binary[star[i].binind].a;
    m1= binary[star[i].binind].m1 * madhoc;
    m2= binary[star[i].binind].m2 * madhoc;
    Eb= m1*m2/2./a;
    Eb/= kTcore;

    if (Eb< ktmin) continue;

    cent.N_bin++;
    Mbincentral += star[i].m * madhoc;
    cent.v_bin_rms += sqr(star[i].vr) + sqr(star[i].vt);
    cent.a_ave += a;
    cent.a2_ave += a*a;
    cent.ma_ave += star[i].m *madhoc * a;
  }
  if (cent.N_bin>0) {
    cent.a_ave/= cent.N_bin;
    cent.a2_ave/= cent.N_bin;
    cent.v_bin_rms/= cent.N_bin;
    cent.ma_ave/= cent.N_bin;
    cent.m_bin_ave= Mbincentral/cent.N_bin;
    cent.n_bin= cent.N_bin/Vcentral;
    cent.rho_bin= Mbincentral/Vcentral;
  }
  
  return(cent);
}

#ifdef USE_MPI
void mpi_clusdyn_calculate(void)
{
	double m=0.0;
	long k=1;

	//MPI3: Potential parallelization possibolity for large N.
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
#ifdef USE_MPI
  phis= MPI_PHI_S(r, j);
#else
  phis= PHI_S(r, j);
#endif
  res= (2.0 * ((E) - (pot + phis)) - Jr);
  return (res);
};

/* gettimeofday Example:
 */
/*
	timeval tim;
	gettimeofday(&tim, NULL);
	double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	do_something_long();
	gettimeofday(&tim, NULL);
	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	printf("%.6lf seconds elapsed\n", t2-t1);
 */

double timeStartSimple()
{
double timeStart=0;
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	timeStart = MPI_Wtime();
#else
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timeStart=tim.tv_sec+(tim.tv_usec/1000000.0);
#endif
return timeStart;
}

void timeEndSimple(double timeStart, double *timeAccum)
{
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	double timeEnd = MPI_Wtime();
	*timeAccum += timeEnd - timeStart;
#else
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timeEnd=tim.tv_sec+(tim.tv_usec/1000000.0);
	*timeAccum += timeEnd - timeStart;
#endif
}

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
	double totTime = 0.0;
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

#ifdef USE_MPI
	Eescaped_old = 0.0;
	Jescaped_old = 0.0;
	Eintescaped_old = 0.0;
	Ebescaped_old = 0.0;
	TidalMassLoss_old = 0.0;
	Etidal_old = 0.0;
#endif
}

void set_global_vars2()
{

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
}


void calc_sigma_new()
{

	//MPI2: mpi_calc_sigma_r(): Tests for precision: Biggest errors are ~ 1e-13. Might be because of catastrphic cancellation/rounding errors?. This might be due to the already imprecise but faster serial code . In a sense, the MPI version might be more precise, because for every processor, actual average is performed for the 1st local star which means errors dont carry over from the previous set of stars which are handled by another processor.
	//MPI3: HASNT BEEN TESTED YET!
#ifdef USE_MPI
	mpi_calc_sigma_r();
#else
	calc_sigma_r(); //requires data from some neighbouring particles. must be handled during data read. look inside function for more comments
#endif
}

void bin_vars_calculate()
{
	/* Calculates some global binary variables - total binary mass,and E. */
	int i, j;
	M_b = 0.0;
	E_b = 0.0;
#ifdef USE_MPI
	for (i=1; i<=mpiEnd-mpiBegin+1; i++)
#else
	for (i=1; i<=clus.N_STAR; i++)
#endif
	{
		j = star[i].binind;
		if (j && binary[j].inuse) {
#ifdef USE_MPI
			M_b += star_m[i];
#else
			M_b += star[i].m;
#endif
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}
}

void calc_potential_new()
{
#ifdef USE_MPI
	mpi_potential_calculate();

	//MPI2: Calculating indices which will be used in all loops till beginning of the main loop. The value 20 depends on the p value used in calc_sigma_new()
	mpiFindIndicesCustom( clus.N_MAX, 20, myid, &mpiBegin, &mpiEnd );
#else
	potential_calculate();
#endif
}

void compute_energy_new()
{
#ifdef USE_MPI

	Eescaped += Eescaped_old;
	Jescaped += Jescaped_old;
	Eintescaped += Eintescaped_old;
	Ebescaped += Ebescaped_old;
	//TidalMassLoss += TidalMassLoss_old;
	Etidal += Etidal_old;

	//MPI2: Tested! No errors.
	mpi_ComputeEnergy();
#else
	//MPI2: see inline for comments
	ComputeEnergy();
#endif
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
		for (i = 1; i <= clus.N_MAX_NEW; i++)
#else
		for (i = 1; i <= clus.N_MAX; i++) 
#endif
			/* reset interacted flag */
			star[i].interacted = 0;
}

void calc_clusdyn_new()
{
	/* calculate dynamical quantities */
#ifdef USE_MPI
	mpi_clusdyn_calculate();
#else
	clusdyn_calculate();
#endif
}

void calc_timestep(gsl_rng *rng)
{
	/* Get new time step */
	//MPI3: Do on all nodes, and remove Bcast
	//for the next step, only simul_relax needs to be done in parallel, and DTrel needs to be broadcasted. rest can be done on each node.
	//Step 2: Do simul_relax() on all procs and others only on root. then broadcast Dt

	/* if tidal mass loss in previous time step is > 5% reduce PREVIOUS timestep by 20% */
	//MPI2: root node
	if ((TidalMassLoss - OldTidalMassLoss) > 0.01) {
		diaprintf("prev TidalMassLoss=%g: reducing Dt by 20%%\n", TidalMassLoss - OldTidalMassLoss);
		Dt = Prev_Dt * 0.8;
	} else if (Dt > 1.1 * Prev_Dt && Prev_Dt > 0 && (TidalMassLoss - OldTidalMassLoss) > 0.02) {
		diaprintf("Dt=%g: increasing Dt by 10%%\n", Dt);
		Dt = Prev_Dt * 1.1;
	}

	TotalTime += Dt;
}

void energy_conservation1()
{
	/* some numbers necessary to implement Stodolkiewicz's
	 * energy conservation scheme */
	int i, g_i;
	for (i = 1; i <= clus.N_MAX_NEW; i++)
	{
		/* saving velocities */
		star[i].vtold = star[i].vt;
		star[i].vrold = star[i].vr;

		/* the following will get updated after sorting and
		 * calling potential_calculate(), needs to be saved 
		 * now */  
#ifdef USE_MPI
		g_i = get_global_idx(i);
		star[i].Uoldrold = star_phi[g_i] + MPI_PHI_S(star_r[g_i], g_i);
#else
		star[i].Uoldrold = star[i].phi + PHI_S(star[i].r, i);
#endif

		/* Unewrold will be calculated after 
		 * potential_calculate() using [].rOld
		 * Unewrnew is [].phi after potential_calculate() */
	}
}

void toy_rejuvenation()
{
	/*Sourav: checking all stars for their possible extinction from old age*/
	//Sourav: toy rejuvenation: DMrejuv storing amount of mass loss per time step
	int i;
	DMrejuv = 0.0;
	if (STAR_AGING_SCHEME > 0) {
		//MPI2: Why is this only till N_MAX?
		//MPI3: This hasnt been changed to take global/local indices. For later.
		for (i=1; i<=clus.N_MAX; i++)
			remove_old_star(TotalTime, i);
	}
}

void new_orbits_calculate()
{
#ifdef USE_MPI
		//MPI2: The 2nd part of tidally_strip_stars() has been moved just before sort.
		OldTidalMassLoss = TidalMassLoss_old;

		/******************************/
		/* get new particle positions */
		/******************************/
		max_r = get_positions();

		//MPI2: Since max_r will be used only be the root node, hopefully the Bcast wont be reqd.
		//MPI3: Should it be AllReduce?
		MPI_Allreduce(MPI_IN_PLACE, &max_r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);		
		MPI_Allreduce(MPI_IN_PLACE, &TidalMassLoss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		

		TidalMassLoss += TidalMassLoss_old;
#else
		OldTidalMassLoss = TidalMassLoss;
		max_r = get_positions();
#endif
}

void energy_conservation2()
{
	int i;
	/* more numbers necessary to implement Stodolkiewicz's
	 * energy conservation scheme */
	for (i = 1; i <= clus.N_MAX_NEW; i++) 
	{
	/* the following cannot be calculated after sorting 
	 * and calling potential_calculate() */
#ifdef USE_MPI 
		star[i].Uoldrnew = potential(star[i].rnew) + MPI_PHI_S(star[i].rnew, get_global_idx(i));
#else
		star[i].Uoldrnew = potential(star[i].rnew) + PHI_S(star[i].rnew, i);
#endif
	}
}

void pre_sort_comm()
{
#ifdef USE_MPI
	int i;

	mpiFindDispAndLenCustom( clus.N_MAX, 20, mpiDisp, mpiLen );

	MPI_Datatype startype;
	MPI_Type_contiguous( sizeof(star_t), MPI_BYTE, &startype );
	MPI_Type_commit( &startype );
	
	//MPI2: Collecting the r and m arrays into the original star structure for sorting.
/*
	for(i=mpiBegin; i<=mpiEnd; i++) {
		star[i].r = star_r[i];
		star[i].m = star_m[i];
	}
	for(i=clus.N_MAX+2; i<=clus.N_MAX_NEW; i++) {
		star[i].r = star_r[i];
		star[i].m = star_m[i];
	}
*/

	//MPI2: Collection of old stars.
	if(myid == 0)
		MPI_Gatherv(MPI_IN_PLACE, mpiLen[myid], startype, star, mpiLen, mpiDisp , startype , 0, MPI_COMM_WORLD);
	else
		MPI_Gatherv(&star[mpiDisp[myid]], mpiLen[myid], startype, star, mpiLen, mpiDisp , startype , 0, MPI_COMM_WORLD);


	//MPI2: Collection of new stars.
	MPI_Gather( &(clus.N_MAX_NEW), 1, MPI_INT, disp, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myid==0)
	{
		disp[0] = clus.N_MAX + 2;
		len[0] = clus.N_MAX_NEW - clus.N_MAX - 1;
		for(i=1; i<procs; i++)
		{
			len[i] = disp[i] - clus.N_MAX - 1;
			disp[i] = disp[i-1] + len[i-1];
		}
		clus.N_MAX_NEW = disp[procs-1] + len[procs-1] - 1;
	}

	if(myid == 0)
		MPI_Gatherv(MPI_IN_PLACE, len[0], startype, star, len, disp , startype , 0, MPI_COMM_WORLD);
	else
		MPI_Gatherv(&star[clus.N_MAX+2], clus.N_MAX_NEW - clus.N_MAX - 1, startype, star, len, disp , startype , 0, MPI_COMM_WORLD);
#endif
}

void qsorts_new(void)
{
#ifdef USE_MPI
	double tmpTimeStart = timeStartSimple();
	MPI_Datatype startype;
	MPI_Type_contiguous( sizeof(star_t), MPI_BYTE, &startype );
	MPI_Type_commit( &startype );

	MPI_Datatype binarytype;
	MPI_Type_contiguous( sizeof(binary_t), MPI_BYTE, &binarytype );
	MPI_Type_commit( &binarytype );
/*
	clus.N_MAX = sample_sort_old(star,
									&clus.N_MAX_NEW,
									startype,
									MPI_COMM_WORLD,
									SAMPLESIZE);
*/
	int temp = (int)clus.N_MAX_NEW; //to avoid incompatible pointer type warning
	clus.N_MAX = sample_sort(	star+1,
                					&temp,
                					startype,
                					binarytype,
                					MPI_COMM_WORLD,
                					SAMPLESIZE );
	clus.N_MAX_NEW = temp;

	//MPI_Type_free(startype);
	timeEndSimple(tmpTimeStart, &t_sort_only);
#else
	/* Sorting stars by radius. The 0th star at radius 0 
		and (N_STAR+1)th star at SF_INFINITY are already set earlier.
	 */
	qsorts(star+1,clus.N_MAX_NEW);
#endif
}

void post_sort_comm()
{
#ifdef USE_MPI
	int i;
	MPI_Status stat;
	mpiFindDispAndLenCustom( clus.N_MAX, 20, mpiDisp, mpiLen );

	double *temp_r = (double *) malloc( ((int)clus.N_MAX_NEW+1) * sizeof(double) );
	double *temp_m = (double *) malloc( ((int)clus.N_MAX_NEW+1) * sizeof(double) );
	
	//MPI3: Can be made more efficient by using MPI datatypes.
	for(i=1; i<=clus.N_MAX_NEW; i++) {
		temp_r[i] = star[i].r;
		temp_m[i] = star[i].m;
		//star_r[i] = star[i].r;
		//star_m[i] = star[i].m;
	}
/*
	printf("%d new=%ld\n",myid, clus.N_MAX_NEW);
	if(myid==0)
		for(i=0; i<procs; i++)
			printf("%d %d\n", mpiDisp[i], mpiLen[i]);
*/
	{
		//MPI3: No idea why this is not working. Consult Wei-keng.
		//MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, star_r, mpiLen, mpiDisp, MPI_DOUBLE, MPI_COMM_WORLD);
		//MPI_Allgatherv(MPI_IN_PLACE, mpiLen[myid], MPI_DOUBLE, star_m, mpiLen, mpiDisp, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(temp_r+1, mpiLen[myid], MPI_DOUBLE, star_r, mpiLen, mpiDisp, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(temp_m+1, mpiLen[myid], MPI_DOUBLE, star_m, mpiLen, mpiDisp, MPI_DOUBLE, MPI_COMM_WORLD);
	}

	//MPI3: Is this required?
	MPI_Barrier(MPI_COMM_WORLD);
	free(temp_r);
	free(temp_m);
#endif
}

//OPT: Comment in detail so can be understandable in future
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
	int i, up_bound;

	if(procs == 1)
		return 0;

	for( i = 0; i < procs; i++ )
		if( j >= Start[i] && j <= End[i] )
			return i;

	up_bound = clus.N_MAX + 1;
#ifndef USE_MPI
	for( i=0; i<procs; i++ )
	{
		if( j <= up_bound + created_star_dyn_node[i] )
			return i;
		up_bound += created_star_dyn_node[i];
	}

	for( i=0; i<procs; i++ )
	{
		if( j <= up_bound + created_star_se_node[i] )
			return i;
		up_bound += created_star_se_node[i];
	}
#endif

	eprintf("Star id out of bounds! for id = %d\n", j);
	exit_cleanly(-2, __FUNCTION__);
	return -1;
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


int get_global_idx(int i)
{
//encountering errors while rewriting calc_sigma_r() and the binning spanned over i=0th star and was returning mass as 0 as get_global_idx(0) was returnung the sentinel. Hope all loops start from i=1, otherwise there is going to be a problem.
//	if(i == 0) return 0;
//MPI3: In calc_sigma_new() this happens, so will get this danger warning once. Ignore.
//if(i==0) printf("---------DANGER!!!!!!!!!!! get_global_idx(i) called with i=0\n");

	if(i > End[myid] - Start[myid] + 1)
	{
//		printf("------>NEW STAR? myid=%d i=%d gi=%d\n", myid, i, ( clus.N_MAX + (i - (End[myid] - Start[myid] + 1)) ));
		return ( 1 + clus.N_MAX + (i - (End[myid] - Start[myid] + 1)) ); 
	}

	if(myid == 0) return i;

 	if(i - Start[myid] > clus.N_MAX)
	{
		eprintf("Index less than 1! proc:%d i=%d\n",myid,i);
		exit_cleanly(-1, __FUNCTION__);
	}
	//return ( End[myid-1] + i );
	return ( Start[myid] + i - 1 );
}

int get_local_idx(int i)
{
	if(i == 0) return 0;

	if(myid == 0) return i;

 	if(i - End[myid-1] < 1)
	{
		eprintf("Index less than 1!\n");
		exit_cleanly(-1, __FUNCTION__);
	}

	return ( i - End[myid-1] );
}
