/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

/* core of the code: applies relaxation, does single-single collisions and binary interactions */
void dynamics_apply(double dt, gsl_rng *rng)
{
	//MPI2: Tested for outputs: vr, vt, E and J. Tests performed with same seed for rng of all procs. Check done only for proc 0's values as others cant be tested due to rng. Must test after rng is replaced.
	long j, si, p=AVEKERNEL, N_LIMIT, k, kp, ksin, kbin;
	double SaveDt, S, S_tc, S_coll, S_lombardi, S_tmp, W, v[4], vp[4], w[4], psi, beta, wp, w1[4], w2[4];
	double v_new[4], vp_new[4], w_new[4], P_enc, n_local, vcm[4], rcm=0.0, rperi=0;
	double Trel12;
	int i;
	long Nrel=0, Nrelbeta[4]={0,0,0,0};
	double relbeta[4]={PI/2.0,PI/4.0,PI/8.0,PI/16.0}, maverelbeta[4]={0.0,0.0,0.0,0.0}, raverelbeta[4]={0.0,0.0,0.0,0.0};
	double qaverelbeta[4]={0.0,0.0,0.0,0.0};
	FILE *binfp;
	char filename[1024];
	double mass_k, mass_kp; //Bharath: MPI

#ifdef USE_MPI
	mpi_calc_sigma_r();
#else
	/* Calculate and store velocity dispersion profile for use with breaking binaries later.
	   This can't be calculated later since the properties of the cluster members are changing with time. */
	calc_sigma_r();
#endif

#ifdef USE_MPI
	if(myid==0) 
#endif
	{
		/* useful debugging and file headers */
		if (tcount == 1) {
			fprintf(relaxationfile, "# time");
			for (i=0; i<4; i++) {
				fprintf(relaxationfile, " thetase>%g:f,q,<M>,<r>", relbeta[i]);
			}
			fprintf(relaxationfile, "\n");
		}

		// MPI2: This hopefully does not execute, so dont have to worry abt binaries here.
		/* DEBUG: print out binary information every N steps */
		if (0) {
			/* if (tcount%50==0 || tcount==1) { */
			sprintf(filename, "a_e2.%04ld.dat", tcount);
			binfp = fopen(filename, "w");
			for (j=1; j<=clus.N_MAX; j++) {
				if (star[j].binind) {
					fprintf(binfp, "%g %g\n", binary[star[j].binind].a, sqr(binary[star[j].binind].e));
				}
			}
			fclose(binfp);
		}
		/* DEBUG */
	}	
	
	/* Original dt provided, saved for repeated encounters, where dt is changed */
	SaveDt = dt;
	/* keeps track of stars from disrupted binaries */
	/* the +1 was put there so newly created stars wouldn't disappear, but it causes
	   tidally stripped stars to be removed multiple times... */
	/* clus.N_MAX_NEW = clus.N_MAX+1; */
	//clus.N_MAX_NEW = clus.N_MAX;
	/* binding energy information */
	DE_bb = 0.0;
	DE_bs = 0.0;

	N_LIMIT = clus.N_MAX;

#ifdef USE_MPI
	if(myid==0)
#endif
	gprintf("%s(): performing interactions:\n", __FUNCTION__);

	if (!quiet) {
		fflush(stdout);
	}

	fprintf(logfile, "%s(): performing interactions:", __FUNCTION__);

#ifdef USE_MPI	
	for (si=mpiBegin; si<=mpiEnd-mpiEnd%2-1; si+=2) {
#else
	/* the big loop, with limits chosen so that we omit the last star if it is not paired */
	for (si=1; si<=N_LIMIT-N_LIMIT%2-1; si+=2) {
#endif 
		dt = SaveDt;
		
		k = si;
		kp = si + 1;
	
		//MPI2: Involves rng. To be handled later.	
		/* set dynamical params for this pair */
		calc_encounter_dyns(k, kp, v, vp, w, &W, &rcm, vcm, rng, 1);

		//MPI: Makes use of r values of stars outside range. Assuming r array is global, no change needed for MPI version.
		/* Compute local density */
		n_local = calc_n_local(k, p, N_LIMIT);
	
#ifdef USE_MPI
		mass_k = star_m[k];
		mass_kp = star_m[kp];
#else
		mass_k = star[k].m;
		mass_kp = star[kp].m;
#endif
	
		if (star[k].binind > 0 && star[kp].binind > 0) {
			/* binary--binary cross section */
			rperi = XBB * (binary[star[k].binind].a + binary[star[kp].binind].a);

			if (BINBIN) {
				S = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(mass_k+mass_kp)/(rperi*sqr(W)));
			} else {
				S = 0.0;
			}
		} else if (star[k].binind > 0 || star[kp].binind > 0) {
			if (star[k].binind > 0) {
				kbin = k;
				ksin = kp;
			} else {
				kbin = kp;
				ksin = k;
			}

			/* binary--single cross section */
			rperi = XBS * binary[star[kbin].binind].a;

			if (BINSINGLE) {
				S = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(mass_k+mass_kp)/(rperi*sqr(W)));
			} else {
				S = 0.0;
			}
		} else {
			if (SS_COLLISION) {
				if (TIDAL_CAPTURE) {
					/* single--single tidal capture cross section (Kim & Lee 1999);
					   here we treat a compact object (k>=10) as a point mass, a massive MS star (k=1) as an 
					   n=3 polytrope, and everything else (k=0,2-9) as an n=1.5 polytrope. */
					if (star[k].se_k >= 10 && star[kp].se_k >= 10) {
						/* two compact objects, so simply use sticky sphere approximation */
						S_tc = 0.0;
					} else if (star[k].se_k >= 10 && star[kp].se_k == 1) {
						/* compact object plus n=3 polytrope */
						S_tc = sigma_tc_nd(3.0, madhoc * mass_kp, star[kp].rad, madhoc * mass_k, W);
					} else if (star[k].se_k >= 10) {
						/* compact object plus n=1.5 polytrope */
						S_tc = sigma_tc_nd(1.5, madhoc * mass_kp, star[kp].rad, madhoc * mass_k, W);
					} else if (star[k].se_k == 1 && star[kp].se_k >= 10) {
						/* n=3 polytrope plus compact object */
						S_tc = sigma_tc_nd(3.0, madhoc * mass_k, star[k].rad, madhoc * mass_kp, W);
					} else if (star[kp].se_k >= 10) {
						/* n=1.5 polytrope plus compact object */
						S_tc = sigma_tc_nd(1.5, madhoc * mass_k, star[k].rad, madhoc * mass_kp, W);
					} else if (star[k].se_k == 1 && star[kp].se_k == 1) {
						/* n=3 polytrope plus n=3 polytrope */
						S_tc = sigma_tc_nn(3.0, madhoc * mass_k, star[k].rad, 3.0, madhoc * mass_kp, star[kp].rad, W);
					} else if (star[k].se_k == 1) {
						/* n=3 polytrope plus n=1.5 polytrope */
						S_tc = sigma_tc_nn(3.0, madhoc * mass_k, star[k].rad, 1.5, madhoc * mass_kp, star[kp].rad, W);
					} else if (star[kp].se_k == 1) {
						/* n=1.5 polytrope plus n=3 polytrope */
						S_tc = sigma_tc_nn(1.5, madhoc * mass_k, star[k].rad, 3.0, madhoc * mass_kp, star[kp].rad, W);
					} else {
						/* n=1.5 polytrope plus n=1.5 polytrope */
						S_tc = sigma_tc_nn(1.5, madhoc * mass_k, star[k].rad, 1.5, madhoc * mass_kp, star[kp].rad, W);
					}
					
					/* cross section estimate for Lombardi, et al. (2006) */
					if ((star[k].se_k <= 1 || star[k].se_k >= 10) && (star[kp].se_k >= 2 && star[kp].se_k <= 9 && star[kp].se_k != 7)) {
						rperi = 1.3 * star[kp].rad;
						S_lombardi = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(mass_k+mass_kp)/(rperi*sqr(W)));
					} else if ((star[kp].se_k <= 1 || star[kp].se_k >= 10) && (star[k].se_k >= 2 && star[k].se_k <= 9 && star[k].se_k != 7)) {
						rperi = 1.3 * star[k].rad;
						S_lombardi = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(mass_k+mass_kp)/(rperi*sqr(W)));
					} else {
						S_lombardi = 0.0;
					}

					S_tmp = MAX(S_tc, S_lombardi);
				} else {
					S_tc = 0.0;
					S_lombardi = 0.0;
					S_tmp = 0.0;
				}
				
				/* standard sticky sphere collision cross section */
				rperi = star[k].rad + star[kp].rad;
				S_coll = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(mass_k+mass_kp)/(rperi*sqr(W)));
				
				/* take the max of all cross sections; the event type will be chosen by sampling the impact parameter */
				S = MAX(S_coll, S_tmp);
				rperi = madhoc*(mass_k+mass_kp)/sqr(W) * (-1.0+sqrt(1.0+S/FB_CONST_PI*sqr(W*W/(madhoc*mass_k+madhoc*mass_kp))));
			} else {
				S = 0.0;
			}
		}

		/* calculate encounter probability */
		/* should it be n_local here even for binaries? */
		P_enc = n_local * W * S * (dt * ((double) clus.N_STAR)/log(GAMMA*((double) clus.N_STAR)));
		
		/* warn if something went wrong with the calculation of Dt */
		if (P_enc >= 1.0) {
			wprintf("P_enc = %g >= 1!\n", P_enc);
		}


#ifndef USE_MPI
		curr_st = &st[findProcForIndex(k)];
#endif
		/* do encounter or two-body relaxation */
		if(rng_t113_dbl_new(curr_st) < P_enc) { 
			/* do encounter */
			if (star[k].binind > 0 && star[kp].binind > 0) {
				/* binary--binary */
				print_interaction_status("BB");

				binint_do(k, kp, rperi, w, W, rcm, vcm, rng);
				/* fprintf(collisionfile, "BB %g %g\n", TotalTime, rcm); */
			} else if (star[k].binind > 0 || star[kp].binind > 0) {
				/* binary--single */
				print_interaction_status("BS");

				binint_do(k, kp, rperi, w, W, rcm, vcm, rng);
				/* fprintf(collisionfile, "BS %g %g\n", TotalTime, rcm); */
			} else {
				/* single--single */
				print_interaction_status("SS");

				/* do collision */
				sscollision_do(k, kp, rperi, w, W, rcm, vcm, rng);
				/* fprintf(collisionfile, "SS %g %g\n", TotalTime, rcm); */
			}
		} else if (RELAXATION) {
			/* do two-body relaxation */
			Trel12 = (PI/32.0) * cub(W) / ( ((double) clus.N_STAR) * n_local * sqr((mass_k+mass_kp)*madhoc) ) ;
			beta = (PI/2.0) * sqrt(dt/Trel12);

			/* record statistics on scattering angles */
			Nrel++;
			for (i=0; i<4; i++) {
				if (beta > relbeta[i]) {
					Nrelbeta[i]++;
					qaverelbeta[i] += MAX(mass_k, mass_kp)/MIN(mass_k, mass_kp);
					maverelbeta[i] += (mass_k + mass_kp)/2.0 * units.mstar/MSUN;
					raverelbeta[i] += rcm;
				}
			}

			/* clamp beta at max value */
			if (beta > PI/2.0) {
				beta = PI/2.0;
			}
			
			/* set up coordinate system */
			wp = sqrt(sqr(w[1]) + sqr(w[2]));
			if (wp == 0.0) {
				eprintf("wp=0 \n");
				exit_cleanly(1);
			}
			
			/* You'll notice here that the sign on w1 is opposite that of what's shown in Kris Joshi's
			   paper.  The sign has now been fixed so that (w1, w2, w) define a right-handed coordinate
			   system, as such: \^w1 x \^w2 = \^w */
			w1[1] = -w[2] * W / wp;
			w1[2] = w[1] * W / wp;
			w1[3] = 0.0;
			w2[1] = -w[1] * w[3] / wp;
			w2[2] = -w[2] * w[3] / wp;
			w2[3] = wp;
			
			psi = rng_t113_dbl_new(curr_st) * 2 * PI;

			for (j = 1; j <= 3; j++) {
				w_new[j] = w[j] * cos(beta) + w1[j] * sin(beta) * cos(psi) + w2[j] * sin(beta) * sin(psi);
			}
			
			for (j = 1; j <= 3; j++) {
				v_new[j] = v[j] - mass_kp / (mass_k + mass_kp) * (w_new[j] - w[j]);
				vp_new[j] = vp[j] + mass_k / (mass_k + mass_kp) * (w_new[j] - w[j]);
			}
			
			/* check to see whether stars should be eaten by central BH */
			if (cenma.m > 0.0 && BH_LOSS_CONE) {
				if (star[k].E < 0.0) {
					bh_rand_walk(k, v_new, vcm, beta, dt);
				}
				if (star[kp].E < 0.0) {
					bh_rand_walk(kp, vp_new, vcm, beta, dt);
				}
			}
			
			/* set new velocities for both stars */
			star[k].vr = v_new[3];
			star[k].vt = sqrt(sqr(v_new[1]) + sqr(v_new[2]));
			star[kp].vr = vp_new[3];
			star[kp].vt = sqrt(sqr(vp_new[1]) + sqr(vp_new[2]));
			if(star[k].vr == 0 || star[k].vt == 0 || star[kp].vr == 0 || star[kp].vt == 0)
			
			/* Calculate new energies by recomputing E = PE + KE using new velocity*/ 
			set_star_EJ(k);
			set_star_EJ(kp);
		}
	}

	/* print relaxation information */
	fprintf(relaxationfile, "%g", TotalTime);
	for (i=0; i<4; i++) {
		fprintf(relaxationfile, " %g %g %g %g", 
				((double) Nrelbeta[i])/((double) Nrel), 
				qaverelbeta[i]/((double) Nrelbeta[i]),
				maverelbeta[i]/((double) Nrelbeta[i]),
				raverelbeta[i]/((double) Nrelbeta[i]));
	}
	fprintf(relaxationfile, "\n");

	/* put newline on "...performing interactions..." line */
#ifdef USE_MPI
	if(myid==0)
#endif
	gprintf("\n");
	fprintf(logfile, "\n");

	//MPI2: Binaries, ignoring for now.
	/* break pathologically wide binaries */
	break_wide_binaries();
}
