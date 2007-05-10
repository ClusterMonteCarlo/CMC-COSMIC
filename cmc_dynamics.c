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
	long j, si, p=AVEKERNEL, N_LIMIT, k, kp, ksin, kbin;
	double SaveDt, S, W, v[4], vp[4], w[4], psi, beta, wp, w1[4], w2[4];
	double v_new[4], vp_new[4], w_new[4], P_enc, n_local, vcm[4], rcm=0.0, rperi;
	double Trel12;
	int i;
	long Nrel=0, Nrelbeta[4]={0,0,0,0};
	double relbeta[4]={PI/2.0,PI/4.0,PI/8.0,PI/16.0}, maverelbeta[4]={0.0,0.0,0.0,0.0}, raverelbeta[4]={0.0,0.0,0.0,0.0};
	double qaverelbeta[4]={0.0,0.0,0.0,0.0};
	FILE *binfp;
	char filename[1024];

	/* Calculate and store velocity dispersion profile for use with breaking binaries later.
	   This can't be calculated later since the properties of the cluster members are changing with time. */
	calc_sigma_r();

	/* useful debugging and file headers */
	if (tcount == 1) {
		fprintf(relaxationfile, "# time");
		for (i=0; i<4; i++) {
			fprintf(relaxationfile, " thetase>%g:f,q,<M>,<r>", relbeta[i]);
		}
		fprintf(relaxationfile, "\n");
	}

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
	
	/* Original dt provided, saved for repeated encounters, where dt is changed */
	SaveDt = dt;
	/* keeps track of stars from disrupted binaries */
	/* the +1 was put there so newly created stars wouldn't disappear, but it causes
	   tidally stripped stars to be removed multiple times... */
	/* clus.N_MAX_NEW = clus.N_MAX+1; */
	clus.N_MAX_NEW = clus.N_MAX;
	/* binding energy information */
	DE_bb = 0.0;
	DE_bs = 0.0;

	N_LIMIT = clus.N_MAX;
	
	gprintf("%s(): performing interactions:", __FUNCTION__);
	if (!quiet) {
		fflush(stdout);
	}
	fprintf(logfile, "%s(): performing interactions:", __FUNCTION__);
	
	/* the big loop, with limits chosen so that we omit the last star if it is not paired */
	for (si=1; si<=N_LIMIT-N_LIMIT%2-1; si+=2) {
		dt = SaveDt;
		
		k = si;
		kp = si + 1;
		
		/* set dynamical params for this pair */
		calc_encounter_dyns(k, kp, v, vp, w, &W, &rcm, vcm, rng, 1);

		/* Compute local density */
		n_local = calc_n_local(k, p, N_LIMIT);
		
		if (star[k].binind > 0 && star[kp].binind > 0) {
			/* binary--binary cross section */
			rperi = XBB * (binary[star[k].binind].a + binary[star[kp].binind].a);

			if (BINBIN) {
				S = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
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
				S = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
			} else {
				S = 0.0;
			}
		} else {
			/* single--single star physical collision cross section */
			rperi = XCOLL * (star[k].rad + star[kp].rad);

			if (SS_COLLISION) {
				S = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
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

		/* do encounter or two-body relaxation */
		if (rng_t113_dbl() < P_enc) {
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
				sscollision_do(k, kp, rcm, vcm);
				/* fprintf(collisionfile, "SS %g %g\n", TotalTime, rcm); */
			}
		} else if (RELAXATION) {
			/* do two-body relaxation */
			Trel12 = (PI/32.0) * cub(W) / ( ((double) clus.N_STAR) * n_local * sqr((star[k].m+star[kp].m)*madhoc) ) ;
			beta = (PI/2.0) * sqrt(dt/Trel12);
			
			/* record statistics on scattering angles */
			Nrel++;
			for (i=0; i<4; i++) {
				if (beta > relbeta[i]) {
					Nrelbeta[i]++;
					qaverelbeta[i] += MAX(star[k].m, star[kp].m)/MIN(star[k].m, star[kp].m);
					maverelbeta[i] += (star[k].m + star[kp].m)/2.0 * units.mstar/MSUN;
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
			
			psi = rng_t113_dbl() * 2 * PI;
			for (j = 1; j <= 3; j++) {
				w_new[j] = w[j] * cos(beta) + w1[j] * sin(beta) * cos(psi) + w2[j] * sin(beta) * sin(psi);
			}
			
			for (j = 1; j <= 3; j++) {
				v_new[j] = v[j] - star[kp].m / (star[k].m + star[kp].m) * (w_new[j] - w[j]);
				vp_new[j] = vp[j] + star[k].m / (star[k].m + star[kp].m) * (w_new[j] - w[j]);
			}
			
			/* set new velocities for both stars */
			star[k].vr = v_new[3];
			star[k].vt = sqrt(sqr(v_new[1]) + sqr(v_new[2]));
			star[kp].vr = vp_new[3];
			star[kp].vt = sqrt(sqr(vp_new[1]) + sqr(vp_new[2]));
			
			/* Calculate new energies by recomputing E = PE + KE using new velocity*/ 
			set_star_EJ(k);
			set_star_EJ(kp);

			/* check to see whether stars should be eaten by central BH */
			if (cenma.m > 0.0 && BH_LOSS_CONE) {
				if (star[k].E < 0.0) {
					bh_rand_walk(k, beta, dt);
				}
				if (star[kp].E < 0.0) {
					bh_rand_walk(kp, beta, dt);
				}
			}
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
	gprintf("\n");
	fprintf(logfile, "\n");

	/* break pathologically wide binaries */
	break_wide_binaries();
}
