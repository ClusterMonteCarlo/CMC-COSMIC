/* -*- linux-c -*- */

/* use commandline option -DUSE_THREADS to turn on threading */
/*#define USE_THREADS*/
#ifdef USE_THREADS
#include <pthread.h>
#define NUM_THREADS	2
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

double GetTimeStep(void) {
	long i;
	double Ai, m_min, m_max, DTrel, Xcoll, Tcoll, DTcoll;
	double Xbb, Tbb, DTbb, Xbs, Tbs, DTbs;
	
	/* calculate Ai for central region */
	Ai = 6.0 * ((double) central.N) * sqr(central.m_ave * ((double) clus.N_STAR)) / 
		((cub(star[central.N].r)-cub(star[1].r)) * sqrt(cub(central.w2_ave)));

	/* set DT_FACTOR automatically */
	m_max = central.m_ave * ((double) clus.N_STAR);
	m_min = central.m_ave * ((double) clus.N_STAR);
	for (i=1; i<=MIN(((long) N_core), clus.N_MAX); i++) {
		m_max = MAX(m_max, star[i].m);
		m_min = MIN(m_min, star[i].m);
	}
	DT_FACTOR = m_max/m_min;

	/* calculate the relaxation timestep */
	/* set by the maximum allowed value of sin^2 beta */
	DTrel = SIN2BETA_MAX * ((double) clus.N_STAR) / Ai / DT_FACTOR;
	/* set to be a fraction of the central relaxation time, as is done by Freitag */
	/* DTrel = 1.0e-2 * PI/128.0 * cub(central.v_rms) / (central.n*sqr(central.m_ave)) / ((double) clus.N_STAR) / DT_FACTOR; */
	Dt = DTrel;

	/* calculate DTcoll, using the expression from Freitag & Benz (2002) (their paper II) */
	if (central.N_sin != 0) {
		/* X defines pericenter needed for collision: r_p = X (R_1+R_2) */
		Xcoll = 1.0;
		Tcoll = 1.0 / (16.0 * sqrt(PI) * central.n_sin * sqr(Xcoll) * (central.v_sin_rms/sqrt(3.0)) * central.R2_ave * 
			       (1.0 + central.mR_ave/(2.0*Xcoll*sqr(central.v_sin_rms/sqrt(3.0))*central.R2_ave))) * 
			log(GAMMA * ((double) clus.N_STAR)) / ((double) clus.N_STAR);
	} else {
		Tcoll = GSL_POSINF;
	}
	DTcoll = 1.0e-4 * Tcoll;
	Dt = MIN(Dt, DTcoll);

	/* calculate DTbb, using a generalization of the expression for Tcoll */
	if (central.N_bin != 0) {
		/* X defines pericenter needed for "strong" interaction: r_p = X (a_1+a_2) */	
		Xbb = 3.5;
		Tbb = 1.0 / (16.0 * sqrt(PI) * central.n_bin * sqr(Xbb) * (central.v_bin_rms/sqrt(3.0)) * central.a2_ave * 
			     (1.0 + central.ma_ave/(2.0*Xbb*sqr(central.v_bin_rms/sqrt(3.0))*central.a2_ave))) * 
			log(GAMMA * ((double) clus.N_STAR)) / ((double) clus.N_STAR);
	} else {
		Tbb = GSL_POSINF;
	}
	DTbb = 1.0e-4 * Tbb;
	Dt = MIN(Dt, DTbb);

	/* calculate DTbs, using a generalization of the expression for Tcoll */
	if (central.N_bin != 0 && central.N_sin != 0) {
		/* X defines pericenter needed for "strong" interaction: r_p = X a */
		Xbs = 3.5;
		Tbs = 1.0 / (4.0 * sqrt(PI) * central.n_sin * sqr(Xbs) * (central.v_rms/sqrt(3.0)) * central.a2_ave * 
			     (1.0 + central.m_ave/(Xbs*sqr(central.v_rms/sqrt(3.0))*central.a_ave))) * 
			log(GAMMA * ((double) clus.N_STAR)) / ((double) clus.N_STAR);
	} else {
		Tbs = GSL_POSINF;
	}
	DTbs = 1.0e-4 * Tbs;
	Dt = MIN(Dt, DTbs);

	/* DEBUG */
	dprintf("TotalTime=%g Dt=%g DTrel=%g DTcoll=%g DTbb=%g DTbs=%g\n", TotalTime, Dt, DTrel, DTcoll, DTbb, DTbs);
	/* DEBUG */

	/* this variable is not used except to be printed out */
	Sin2Beta = Ai * Dt / ((double) clus.N_STAR);

	return (Dt);
}

/* removes tidally-stripped stars */
void sniff_stars(void) {
	double phi_rtidal, phi_zero;
	long i, j;
	
	j = 0;
	Etidal = 0.0;
	OldTidalMassLoss = TidalMassLoss;
	/* get new particle positions */
	max_r = get_positions();
	DTidalMassLoss = TidalMassLoss - OldTidalMassLoss;
		
	fprintf(stdout, 
	   "sniff_stars(): iteration %ld: OldTidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
		j, OldTidalMassLoss, DTidalMassLoss);
	fprintf(logfile, 
	   "sniff_stars(): iteration %ld: OldTidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
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
		for (i = 1; i <= clus.N_MAX; i++) {
			if (star[i].r_apo > Rtidal
					&& star[i].rnew < 1000000) {
				star[i].rnew = SF_INFINITY;	/* tidally stripped star */
				star[i].vrnew = 0.0;
				star[i].vtnew = 0.0;
				Eescaped += star[i].E * star[i].m / clus.N_STAR;
				Jescaped += star[i].J * star[i].m / clus.N_STAR;
				DTidalMassLoss += star[i].m / clus.N_STAR;
				Etidal += star[i].E * star[i].m / clus.N_STAR;
				fprintf(escfile,
					"%ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
					tcount, TotalTime, star[i].m,
					star[i].r, star[i].vr, star[i].vt, star[i].r_peri,
					star[i].r_apo, Rtidal, phi_rtidal, phi_zero, star[i].E, star[i].J);

				if (Etotal.K + Etotal.P - Etidal >= 0)
					break;
			}
		}
		j++;
		TidalMassLoss = TidalMassLoss + DTidalMassLoss;
		fprintf(stdout, "sniff_stars(): iteration %ld: TidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
			j, TidalMassLoss, DTidalMassLoss);
		fprintf(logfile, "sniff_stars(): iteration %ld: TidalMassLoss=%.6g DTidalMassLoss=%.6g\n",
			j, TidalMassLoss, DTidalMassLoss);
	} while (DTidalMassLoss > 0 && (Etotal.K + Etotal.P - Etidal) < 0);
}

static void remove_star(long j, double phi_rtidal, double phi_zero) {
	double E, J;

	dprintf("removing star: id=%ld m=%g\n", star[j].id, star[j].m);

	E = star[j].E;
	J = star[j].J;
	star[j].rnew = SF_INFINITY;	/* tidally stripped star */
	star[j].vrnew = 0.0;
	star[j].vtnew = 0.0;
	Eescaped += E * star[j].m / clus.N_STAR;
	Jescaped += J * star[j].m / clus.N_STAR;
	TidalMassLoss += star[j].m / clus.N_STAR;
	Etidal += E * star[j].m / clus.N_STAR;
	fprintf(escfile, "%ld %.8g %.8g ",
		tcount, TotalTime, star[j].m);
	fprintf(escfile, "%.8g %.8g %.8g ",
		star[j].r, star[j].vr, star[j].vt);
	fprintf(escfile, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
		0.0, 0.0, Rtidal, phi_rtidal, phi_zero, E, J);
}

void remove_star_center(long j) {
	dprintf("removing star: id=%ld m=%g\n", star[j].id, star[j].m);
	
	star[j].rnew = SF_INFINITY;	/* send star to infinity         */
	star[j].m = DBL_MIN;		/* set mass to very small number */
	star[j].vrnew = 0.0;		/* setup vr and vt for           */
	star[j].vtnew = 0.0;		/*		future calculations  */
}

/**************** Get Positions and Velocities ***********************/
/*	Requires indexed (sorted in increasing r) stars with potential
	computed in star[].phi and N_MAX set. Uses sE[], sJ[] and sr[]
	from previous iteration. Returns positions and velocities in
	srnew[], svrnew[], and svtnew[]. Returns Max r for all stars, or
	-1 on error. */
#ifdef USE_THREADS
void *get_positions_loop(void *get_pos_dat_arg){
	struct get_pos_str *get_pos_dat = (struct get_pos_str *) get_pos_dat_arg;
	int taskid = get_pos_dat->taskid;
	gsl_rng *thr_rng = get_pos_dat->thr_rng;
#else
void get_positions_loop(struct get_pos_str *get_pos_dat){
#endif
	long i, j, k, kmin, kmax, kk, ktemp, i1, si;
	double r=0.0, vr=0.0, vt, rk, rk1, a, b, rmin, rmax, E, J, Uk, Uk1;
	double g1, g2, F, s0, g0, dQdr_min, dQdr_max, drds, max_rad, pot;
	double Qtemp, X=0.0, Q;	/* Max value of Q found. Should be > 0 */
	long N_LIMIT;
	double phi_rtidal, phi_zero;

	N_LIMIT = get_pos_dat->N_LIMIT;
	phi_rtidal = get_pos_dat->phi_rtidal;
	phi_zero = get_pos_dat->phi_zero;
	max_rad = get_pos_dat->max_rad;

#ifdef USE_THREADS
	for (si = taskid+1; si <= clus.N_MAX_NEW; si+=NUM_THREADS) { /* Repeat for all stars */
#else
	for (si = 1; si <= clus.N_MAX_NEW; si++) { /* Repeat for all stars */
#endif
		j = si;

		E = star[j].E;
		J = star[j].J;
		
		/* ignore halo stars during sub timesteps */
		if (si > N_LIMIT && si <= clus.N_MAX && star[j].r_peri > sub.rmax) {
			continue;
		}

		/* remove unbound or massless stars */
		/* note that energy lost due to stellar evolution is subtracted
		   at the time of mass loss in DoStellarEvolution */
		if (E >= 0.0 || star[j].m < ZERO) {
			remove_star(j, phi_rtidal, phi_zero);
			continue;
		}

		/* Q(si) must be positive (selected that way in last time step!) */
		ktemp = si;

		/* for newly created stars, position si is not ordered, so search */
		if (si > clus.N_MAX) {
			dprintf("doing stupid linear search due to newly created star: si=%ld\n", si);
			ktemp = 0;
			while (ktemp < clus.N_MAX && star[ktemp].r < star[j].rnew) {
				ktemp++;
			}
		}

		/* this is not right for newly created stars; can be negative since it just uses
		   the position at ktemp */
		Qtemp = function_Q(ktemp, E, J);

		kk = ktemp;
		if (Qtemp <= 0) { /* possibly a circular orbit */
			while (function_Q(kk + 1, E, J) > function_Q(kk, E, J) 
			    && function_Q(kk, E, J) < 0 && kk <= clus.N_MAX) {
				kk++;
			}
			while (function_Q(kk - 1, E, J) > function_Q(kk, E, J) 
			    && function_Q(kk, E, J) < 0 && kk >= 1) {
				kk--;
			}

			ktemp = kk;

			if (function_Q(ktemp, E, J) < 0) {
				/* star is on an almost circular orbit... 
				 * so rmin = rmax (approx) */
				star[j].rnew = star[j].r;
				star[j].vrnew = star[j].vr;
				star[j].vtnew = star[j].vt;
				dprintf("circular orbit found: si=%ld sr=%g svr=%g svt=%g J=%g E=%g\n",
					si, star[j].r, star[j].vr, star[j].vt, star[j].J, star[j].E);
			}
			continue;
		}

		kmin = FindZero_Q(0, ktemp, E, J);
		
		while (function_Q(kmin, E, J) > 0 && kmin > 0)
			kmin--;
		while (function_Q(kmin + 1, E, J) < 0 
				&& kmin + 1 < ktemp)
			kmin++;
		
		i = kmin;
		i1 = kmin + 1;
		rk = star[i].r;
		rk1 = star[i1].r;
		Uk = star[i].phi;
		Uk1 = star[i1].phi;
		Q = 2.0 * E - 2.0 * Uk1 - J * J / (rk1 * rk1);

		a = (Uk1 - Uk) / (1 / rk1 - 1 / rk);
		b = (Uk / rk1 - Uk1 / rk) / (1 / rk1 - 1 / rk);

		rmin = J * J / (-a + sqrt(a * a - 2.0 * J * J * (b - E)));
		dQdr_min = 2.0 * J * J / (rmin * rmin * rmin) + 2.0 * a / (rmin * rmin);

		/*  For rmax- Look for rk, rk1 such that 
		 *  Q(rk) > 0 > Q(rk1) */

		kmax = FindZero_Q(ktemp, clus.N_MAX + 1, E, J);

		while (function_Q(kmax, E, J) < 0 && kmax > ktemp)
			kmax--;
		while (kmax + 1 < clus.N_MAX
				&& function_Q(kmax + 1, E, J) > 0 )
			kmax++;

		i = kmax;
		i1 = kmax + 1;
		rk = star[i].r;
		rk1 = star[i1].r;
		Uk = star[i].phi;
		Uk1 = star[i1].phi;
		Q = 2.0 * E - 2.0 * Uk1 - J * J / (rk1 * rk1);

		a = (Uk1 - Uk) / (1 / rk1 - 1 / rk);
		b = (Uk / rk1 - Uk1 / rk) / (1 / rk1 - 1 / rk);

		rmax = (-a + sqrt(a * a - 2.0 * J * J * (b - E))) / (2.0 * (b - E));
		dQdr_max = 2.0 * J * J / (rmax * rmax * rmax) + 2.0 * a / (rmax * rmax);

		if (rmin > rmax) {
			star[j].rnew = star[j].r;
			star[j].vrnew = star[j].vr;
			star[j].vtnew = star[j].vt;
			eprintf("error: rmin>rmax:  si=%ld sr=%g svr=%g svt=%g J=%g E=%g\n",
				si, star[j].r, star[j].vr, star[j].vt, star[j].J, star[j].E);
			eprintf("\trmin=%g rmax=%g kmin=%ld kmax=%ld\n", rmin, rmax, kmin, kmax);
			continue;
		}
		/* Check for rmax > R_MAX (tidal radius) */

		if (rmax >= Rtidal) {
			remove_star(j, phi_rtidal, phi_zero);
			continue;
		}

		g1 = sqrt(3.0 * (rmax - rmin) / dQdr_min);	/* g(-1) */
		g2 = sqrt(-3.0 * (rmax - rmin) / dQdr_max);	/* g(+1) */

		F = 1.2 * MAX(g1, g2);
		
		for (k = 1; k <= N_TRY; k++) {
#ifdef USE_THREADS
			X = gsl_rng_uniform(thr_rng);
#else
			X = rng_t113_dbl();
#endif
			/* Stodolkiewicz's method to prevent very circular orbits. 
			   Note that, 1.2 in F = 1.2 * ...  may still be
			   inadequate, and Marc has a better but more complicated
			   method for that problem. What's below helps with energy
			   conservation more than anything else. */
			if(X>0.95) X=0.97;
			if(X<0.05) X=0.03;

			s0 = 2.0 * X - 1.0;	 /* random -1 < s0 < 1 */
#ifdef USE_THREADS
			g0 = F * gsl_rng_uniform(thr_rng); /* random  0 < g0 < F */
#else
			g0 = F * rng_t113_dbl(); /* random  0 < g0 < F */
#endif

			r = 0.5 * (rmin + rmax) + 0.25 * (rmax - rmin) * (3.0 * s0 - s0 * s0 * s0);

			pot = potential(r);

			drds = 0.25 * (rmax - rmin) * (3.0 - 3.0 * s0 * s0);
			pot = 2.0 * E - 2.0 * pot - J * J / r / r;

			if (pot >= 0.0) {
				vr = sqrt(pot);
			} else {
				dprintf("circular orbit: vr^2<0: setting vr=0: si=%ld r=%g rmin=%g rmax=%g vr^2=%g X=%g E=%g J=%g\n", si, r, rmin, rmax, pot, X, E, J);
				if (isnan(pot)) {
					eprintf("fatal error: pot=vr^2==nan!\n");
					exit_cleanly(-1);
				}
				vr = 0;
			}
			if (g0 < 1.0 / vr * drds)	/* if g0 < g(s0) then success! */
				break;
		}
		if (k == N_TRY + 1) {
			eprintf("N_TRY exceeded\n");
			exit_cleanly(-1);
		}

		/* remove stars if they are too close to center,
		 * ie. r < MINIMUM_R.
		 * Add their mass to CentralMass */
		//if (rmax < MINIMUM_R){
		if (r<MINIMUM_R && rmin>0.3*rmax){
#ifdef USE_THREADS
			get_pos_dat->CMincr.m += star[j].m;
			get_pos_dat->CMincr.E += (star[j].phi 
				+ star[j].vr*star[j].vr	+ star[j].vt*star[j].vt ) 
					/2.0 *star[j].m/clus.N_STAR;
#else
			cenma.m += star[j].m;
			cenma.E += (star[j].phi + star[j].vr*star[j].vr
				+ star[j].vt*star[j].vt ) / 2.0 *star[j].m/clus.N_STAR;
#endif
			remove_star_center(j);
			continue;
		}

		star[j].X = X;
		star[j].r_peri = rmin;
		star[j].r_apo = rmax;
		
		/* pick random sign for v_r */
#ifdef USE_THREADS
		if (gsl_rng_uniform(thr_rng) < 0.5)
#else
		if (rng_t113_dbl() < 0.5)
#endif
			vr = -vr;

		vt = J / r;

		star[j].rnew = r;
		star[j].vrnew = vr;
		star[j].vtnew = vt;

		if (r > max_rad)
			max_rad = r;
	}/* Next si */
	get_pos_dat->max_rad = max_rad;
#ifdef USE_THREADS
	if(get_pos_dat->taskid>0) 
		pthread_exit(NULL);
	return (NULL);
#endif
}

double get_positions(){
#ifdef USE_THREADS
	pthread_t threads[NUM_THREADS];
      pthread_attr_t attr;
	int rc, t;
	struct get_pos_str *get_positions_data_array;
#endif
	double max_rad, phi_rtidal, phi_zero;
	long N_LIMIT;
	struct get_pos_str get_positions_data;

#ifdef USE_THREADS
	/* Initialize and set thread detached attribute */
      pthread_attr_init(&attr);
	get_positions_data_array = malloc(NUM_THREADS*sizeof(struct get_pos_str));
#endif

	max_rad = 0.0;	/* max radius for all stars, returned on success */

	phi_rtidal = potential(Rtidal);
	phi_zero = potential(0.0);

	if (sub.count == 0) {
		N_LIMIT = clus.N_MAX; /* do FULL time step */
	} else {
		N_LIMIT = sub.N_MAX; /* do sub step only */
	}

	get_positions_data.max_rad = max_rad;
	get_positions_data.phi_rtidal = phi_rtidal;
	get_positions_data.phi_zero = phi_zero;
	get_positions_data.N_LIMIT = N_LIMIT;
	get_positions_data.taskid = 0;
	get_positions_data.CMincr.m = 0.0;
	get_positions_data.CMincr.E = 0.0;
	get_positions_data.thr_rng = NULL;
	
#ifdef USE_THREADS
	/* create the thread(s) */
	for(t=1;t<NUM_THREADS;t++) {
		get_positions_data_array[t] = get_positions_data;
		get_positions_data_array[t].taskid = t;
		get_positions_data_array[t].thr_rng = gsl_rng_alloc(gsl_rng_taus2);
		gsl_rng_set(get_positions_data_array[t].thr_rng,
				rng_t113_int());
		rc = pthread_create(&threads[t], NULL,
            	            get_positions_loop, 
					(void *) &get_positions_data_array[t]);
		if (rc) {
                  eprintf("return code from pthread_create() is %d\n", rc);
                  exit(-1);
            }
      }
	/* do work on main thread as well */
	get_positions_data_array[0] = get_positions_data;
	get_positions_data_array[0].taskid = 0;
	get_positions_data_array[0].thr_rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(get_positions_data_array[0].thr_rng,
			rng_t113_int());
	get_positions_loop(&get_positions_data_array[0]);
	/* Free attribute and wait for the other threads */
	/* Update variables when threads join */
      pthread_attr_destroy(&attr);
      for(t=1;t < NUM_THREADS;t++) {
            rc = pthread_join(threads[t], NULL);
            if (rc) {
                  eprintf("return code from pthread_join() is %d\n", rc);
                  exit(-1);
            }
		if (get_positions_data_array[t].max_rad > max_rad) {
			max_rad = get_positions_data_array[t].max_rad;
		}
      }
      for(t=0;t < NUM_THREADS;t++) {
		cenma.m += get_positions_data_array[t].CMincr.m;
		cenma.E += get_positions_data_array[t].CMincr.E;
		gsl_rng_free(get_positions_data_array[t].thr_rng);
	}
	free(get_positions_data_array);
#else
	get_positions_loop(&get_positions_data);
	max_rad = get_positions_data.max_rad;
#endif
	return (max_rad);
}
