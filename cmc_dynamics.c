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
	double SaveDt, S, S_tc, S_coll, S_lombardi, S_tmp, W, v[4], vp[4], w[4], psi, beta, wp, w1[4], w2[4];
	double v_new[4], vp_new[4], w_new[4], P_enc, n_local, vcm[4];
	double vel1[4], vel2[4], vel3[4], vel1a[4], vel2a[4], vel1b[4], vel3b[4], rcm=0.0, rperi;
	double Trel12;
	int i;
	long Nrel=0, Nrelbeta[4]={0,0,0,0};
	double relbeta[4]={PI/2.0,PI/4.0,PI/8.0,PI/16.0}, maverelbeta[4]={0.0,0.0,0.0,0.0}, raverelbeta[4]={0.0,0.0,0.0,0.0};
	double qaverelbeta[4]={0.0,0.0,0.0,0.0};
	FILE *binfp;
	char filename[1024];
/* Meagan: added these variables for three-body binary formation */
	long sq, k1, k2, k3, binary_companion, binary_formed, knew;
	double sigma_local, v1[4], v2[4], v3[4], vrel12[4], vrel13[4], vcm12[4], vcm13[4], vrel3[4], vrel2[4], v3_rel_mag, v2_rel_mag; 
	double angle1, angle2, angle3, angle4; 
	double prefactor_3bb, rate_12, rate_13, P_12, P_13, P_3bb;
	double Y, Z;
	double n_threshold, eta_min=1.0, eta, randno, r_p, ecc_max, ecc, binary_a;
	long THREEBODYBINARIES=1;	


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




/* Meagan 3/21/09  */
	/* Three-body binary formation */
	/* variables
		n_threshold = density threshold: below threshold, don't do 3bb formation
		THREEBODYBINARIES = flag to turn on 3bb formation (set =1 in .cmc file)
		sq = iteration variable for looping through objects, three at a time
		k, kp, kq = indices for the three neighboring stars, possibly interacting to form binary

	 Check probability of star 1 forming a binary with stars 2 and 3; depends on masses, rel. velocity 
		Do Monte Carlo to decide if either binary will form, given probability, based on 

	*/
	eta_min = 1.0; /* TODO: later set this as a global variable */
/*	eta = 1.0;  TODO: remove this line: only did this bc it was telling me that eta was not defined; but I define it later in function get_eta(), in cmc_dynamics_helper.c. Must solve this problem */

	/* TODO: set n_threshold to reasonable value */
	n_threshold = 0;  /* density threshold: below threshold, don't even bother doing 3BB formation */	
	if (THREEBODYBINARIES) { /* level 1 */
		for (sq=1; sq<=N_LIMIT-N_LIMIT%3-2; sq+=3) { /* level 2 */
			dt = SaveDt;

			k1 = sq;
			k2 = sq + 1;
			k3 = sq + 2;
			
		 	fprintf(threebbfile, "three object indices:  %d  %d  %d\n", k1, k2, k3);	
			n_local = calc_n_local(k1, p, N_LIMIT);
				/* n_local about star 1 */
			if (n_local > n_threshold) { /* level 3 */
				if (star[k1].binind == 0 && star[k2].binind == 0 && star[k3].binind == 0) { /* level 4 */
					/* three single stars: try 3bb formation */
					fprintf(threebbfile, "*****************************************\n");
					fprintf(threebbfile, "Above density threshold, and we have three single stars. Checking for three-body binary formation\n");
					/* calculate P12, probability that 1 and 2 form a binary */
		
					/* prefactor for rate expression for star 1 */
					sigma_local = sigma_array.sigma[k1];	
						/* CHECK: make sure sigma_array.sigma[k] is what I think it is */
					prefactor_3bb = PI * sqr(n_local) * pow(madhoc, 5.0) / pow(sigma_local, 9.0);
						/* units are 1/T_cross  */

				// Computes quantities needed to calculate probabilities:
					// generates randomly oriented velocities for three stars
					// computes relative velocities between candidate binary pairs
					// as well as relative veloc. of the single star wrt the binary
					// and the COM veloc. of the binary
					calc_3bb_encounter_dyns(k1, k2, k3, angle1, angle2, v1, v2, v3, vrel12, vrel13, vrel3, vrel2, vcm12, vcm13, rng);

				// Three body binary formation probabilities
					/* Calc P_12, probability that star 1 forms binary with star 2 */
					/* set dynamical params for this pair */	
					

					/* calculate rate and probability that 1 and 2 form a binary */
/*					rate_12 = prefactor_3bb * (v12 / sigma_local) * pow(((star[k1].m * star[k2].m)/eta_min), 5.0)* (1.0 + 2.0*eta_min)* (1.0 + (v3*(pow((2.0*star[k1].m*star[k2].m), 0.5))));
*/

/*					rate_12 = prefactor_3bb * (v12 / sigma_local) * pow(((star[k1].m * star[k2].m)/eta_min), 5.0) * (1.0 + 2.0*eta_min) * (1.0 + (v3_rel_mag/sigma_local)*pow((2.0*star[k1].m*star[k2].m/((star[k1].m+star[k2].m)*eta_min)), 0.5));
*/
					rate_12=0.5;				

					P_12 = rate_12 * (dt * ((double) clus.N_STAR)/log(GAMMA*((double) clus.N_STAR)));

					/* calculate rate and probability that 1 and 3 form a binary */
/*					rate_13 = prefactor_3bb * (v13 / sigma_local) * pow(((star[k1].m * star[k3].m)/eta_min), 5.0) * (1.0 + 2.0*eta_min) * (1.0 + (v2_rel_mag/sigma_local)*pow((2.0*star[k1].m*star[k3].m/((star[k1].m+star[k3].m)*eta_min)), 0.5));
*/
					rate_13=0.5;

					P_13 = rate_13 * (dt * (
(double) clus.N_STAR)/log(GAMMA*((double) clus.N_STAR)));
					P_3bb = (P_12 + P_13)/2.0;
					
					fprintf(threebbfile, "P_12=%g  P_13=%g  P_3bb=%g\n", P_12, P_13, P_3bb);
					Y = exp(-P_3bb); /* Y falls in range (0, 1) [(likely, unlikely)] for P_3bb in range (+inf, 0) [(likely, unlikely)]; so choose random #, if Y<rand#, then binary will be formed*/
					fprintf(threebbfile, "Y=exp(-P)=%g\n", Y);
					
					if (Y < rng_t113_dbl()) { 
						/* Star 1 WILL form a binary; choose binary companion */
						fprintf(threebbfile, "binary to be formed with stars %d  %d  %d\n", k1, k2, k3);
						binary_formed = 1;
						Z = rng_t113_dbl() * (P_12 + P_13);
						if ( Z < P_12) {   /* 1 & 2 will form binary */
							fprintf(threebbfile, "binary companion is star2\n");
							/* force more massive star to be the primary - represented by index 1 */								
							if (star[k1].m > star[k2].m) {
								// primary=star1, secondary=star2, single=star3
								fprintf(threebbfile, "m1 > m2\n");
								for (j=1; j<=3; j++) {
									vel1[j] = vel1a[j];
									vel2[j] = vel2a[j];
									vel3[j] = vel3b[j];	
								}
							} else { // single=star3 already
								fprintf(threebbfile, "m1 < m2\n");
								k1 = sq + 1; // primary=star2
								k2 = sq; // secondary=star1
								for (j=1; j<=3; j++) {
									vel1[j] = vel2a[j];
									vel2[j] = vel1a[j];
									vel3[j] = vel3b[j];	
								}
							}
						} else {   /* 1 & 3 will form binary */
							/* change indices so that binary components indexed by 1,2, and the single by 3 */
							if (star[k1].m > star[k3].m) {
								k1 = sq; // primary=star1
								k2 = sq + 2; //secondary=star3
								k3 = sq + 1; //single=star2
								for (j=1; j<=3; j++) {
									vel1[j] = vel1b[j];
									vel2[j] = vel3b[j];
									vel3[j] = vel2a[j];		
								}
							} else {
								k1 = sq + 2; // primary=star3
								k2 = sq; // secondary=star1
								k3 = sq + 1; // single=star2
								for (j=1; j<=3; j++) {
									vel1[j] = vel3b[j];
									vel2[j] = vel1b[j];
									vel3[j] = vel2a[j];		
								}
							}
						fprintf(threebbfile, "IDs of three stars: %ld  %ld  %ld\n", star[k1].id, star[k2].id, star[k3].id);
						}
					/* NOW we have the primary, secondary, and single represented by indices 1, 2, 3. The absolute velocities are vel1, vel2, and vel3; FIXME:the potentials are phi1, phi2, and phi3 */	
					} else {
						/* No binary formed */
						binary_formed = 0;
						fprintf(threebbfile, "no binary formed:  %d  %d  %d\n", k1, k2, k3); 
						/* do nothing; check 3bb formation for next 3 stars */
					}
					if (binary_formed==1) {
						/* new function: make_threebodybinary(k1, k2, k3) */
						knew = create_binary(); /* returns an unused binary id */
						fprintf(threebbfile, "new binary created: id 'knew' = %ld\n", knew);
						make_threebodybinary(k1, k2, k3, knew, sigma_local, v1, v2, v3, rng);
						fprintf(threebbfile, "make_threebodybinary() function called\n");
							// k1, k2, k3, knew: indices of primary and secondary, 
							// remaining single star, and new binary index 
					}
				} /* level 4 */
			} /* level 3 */
		} /* level 2 */
	} /* level 1 */
			

			//shouldn't need this anymore!********************
			/*		calc_encounter_dyns(k1, k2, v1, v2, w, &W, &rcm, vcm, rng, 1);
					//v12 = W; // store relative speed of stars 1,2
					// store CM veloc. of
					for (j=1; j<=3; j++) {
						// store rel. vel. of 1 & 2 for calculating P(3bb)
						vrel_12[j] = w[j];
						// store cm velocity and cm position
						vcm_12[j] = vcm[j];
						rcm_12[j] = rcm[j];
 						// store abs. veloc vectors for stars 1,2 for later
						vel1a[j] = v1[j];
						vel2a[j] = v2[j];
					}
			
					calc_encounter_dyns(k1, k3, v1, v3, w, &W, &rcm, vcm, rng, 1);
					v13 = W; // store rel. vel. of 1 & 3 for calculating P(3bb)
					for (j=1; j<=3; j++) {
						// store rel. vel. of 1 & 3 for calculating P(3bb)
						vrel_13[j] = w[j];
						// store cm velocity and cm position
						vcm_13[j] = vcm[j];
						rcm_13[j] = rcm[j];
 						// store abs. veloc vectors for stars 1,3 for later
						vel1b[j] = v1[j];
						vel3b[j] = v2[j];
					}
					// find velocities of candidate single stars rel. to cm veloc of candidate binary pair; this is v3 in the rate formula 
					for (j=1; j<=3; j++) {
						v3_rel[j] = v3[j] - vcm_12[j]; // v3 rel to cm v of 1 & 2
						v2_rel[j] = v2[j] - vcm_13[j]; // v2 rel to cm v of 1 & 3 
					}
					// calc relative speed of candidate single stars to candidate binary pair
					// This is the v3 that appears in the probability formula
					v3_rel_mag = sqrt(sqr(v3_rel[1]) + sqr(v3_rel[2]) + sqr(v3_rel[3]));
					v2_rel_mag = sqrt(sqr(v2_rel[1]) + sqr(v2_rel[2]) + sqr(v2_rel[3]));

//					fprintf(threebbfile, "v3=%g\n",v3);
			// ************************************************
			*/

/***********************************************/	


	/* the big loop, with limits chosen so that we omit the last star if it is not paired */

	while (si<=N_LIMIT%2-1) {	
//	for (si=1; si<=N_LIMIT-N_LIMIT%2-1; si+=2) {
		dt = SaveDt;
		
		while (star[si].interacted = 1) {
			si += 1; // iterate until non-interacted object found
		}
		k = si;  // object 1 for interaction
		si += 1;
		while (star[si].interacted = 1) {
			si += 1; // iterate until non-interacted object found
		}	
		kp = si; // object 2 for interaction
		si += 1;
		
		fprintf(threebbfile, "two-body loop: k = %i, kp = %i\n", k, kp);

//		k = si;
//		kp = si + 1;
		
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
						S_tc = sigma_tc_nd(3.0, madhoc * star[kp].m, star[kp].rad, madhoc * star[k].m, W);
					} else if (star[k].se_k >= 10) {
						/* compact object plus n=1.5 polytrope */
						S_tc = sigma_tc_nd(1.5, madhoc * star[kp].m, star[kp].rad, madhoc * star[k].m, W);
					} else if (star[k].se_k == 1 && star[kp].se_k >= 10) {
						/* n=3 polytrope plus compact object */
						S_tc = sigma_tc_nd(3.0, madhoc * star[k].m, star[k].rad, madhoc * star[kp].m, W);
					} else if (star[kp].se_k >= 10) {
						/* n=1.5 polytrope plus compact object */
						S_tc = sigma_tc_nd(1.5, madhoc * star[k].m, star[k].rad, madhoc * star[kp].m, W);
					} else if (star[k].se_k == 1 && star[kp].se_k == 1) {
						/* n=3 polytrope plus n=3 polytrope */
						S_tc = sigma_tc_nn(3.0, madhoc * star[k].m, star[k].rad, 3.0, madhoc * star[kp].m, star[kp].rad, W);
					} else if (star[k].se_k == 1) {
						/* n=3 polytrope plus n=1.5 polytrope */
						S_tc = sigma_tc_nn(3.0, madhoc * star[k].m, star[k].rad, 1.5, madhoc * star[kp].m, star[kp].rad, W);
					} else if (star[kp].se_k == 1) {
						/* n=1.5 polytrope plus n=3 polytrope */
						S_tc = sigma_tc_nn(1.5, madhoc * star[k].m, star[k].rad, 3.0, madhoc * star[kp].m, star[kp].rad, W);
					} else {
						/* n=1.5 polytrope plus n=1.5 polytrope */
						S_tc = sigma_tc_nn(1.5, madhoc * star[k].m, star[k].rad, 1.5, madhoc * star[kp].m, star[kp].rad, W);
					}
					
					/* cross section estimate for Lombardi, et al. (2006) */
					if ((star[k].se_k <= 1 || star[k].se_k >= 10) && (star[kp].se_k >= 2 && star[kp].se_k <= 9 && star[kp].se_k != 7)) {
						rperi = 1.3 * star[kp].rad;
						S_lombardi = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
					} else if ((star[kp].se_k <= 1 || star[kp].se_k >= 10) && (star[k].se_k >= 2 && star[k].se_k <= 9 && star[k].se_k != 7)) {
						rperi = 1.3 * star[k].rad;
						S_lombardi = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
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
				S_coll = PI * sqr(rperi) * (1.0 + 2.0*madhoc*(star[k].m+star[kp].m)/(rperi*sqr(W)));
				
				/* take the max of all cross sections; the event type will be chosen by sampling the impact parameter */
				S = MAX(S_coll, S_tmp);
				rperi = madhoc*(star[k].m+star[kp].m)/sqr(W) * (-1.0+sqrt(1.0+S/FB_CONST_PI*sqr(W*W/(madhoc*star[k].m+madhoc*star[kp].m))));
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
				sscollision_do(k, kp, rperi, w, W, rcm, vcm, rng);
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
	gprintf("\n");
	fprintf(logfile, "\n");

	/* break pathologically wide binaries */
	break_wide_binaries();
}
