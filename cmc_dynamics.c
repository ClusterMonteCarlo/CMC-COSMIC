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
	long sq, k1, k2, k3, form_binary;
	double n_threshold, triplet_count, num_triplets_averaged=200;
	double ave_local_mass, sigma_local, vrel_ave, v1[4], v2[4], v3[4], vrel12[4], vrel3[4]; 
	double angle1, angle2, angle3, angle4; 
	double eta_min=MIN_BINARY_HARDNESS, Y1, rate_3bb, rate_ave, P_3bb, P_ave;

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




/* Added by Meagan 3/10/11  */
	//================================
	//  Three-body binary formation  |
	//================================

	/* TODO: set n_threshold to reasonable value */
	n_threshold = 0;  /* density threshold: below threshold, don't even bother doing 3BB formation */	
	triplet_count = 0; //  # of potential binary-forming triplets (three single stars)
			   //  Used for outputting binary formation rate, averaged over innermost objects


	//=====================================================================================
	//  Loop through objects, 3 at a time, and check whether all three are single object  |
	//  If three singles, check whether a binary should be formed with stars 1 and 2      |
	//=====================================================================================

	if (THREEBODYBINARIES) { // Flag for turning on three-body binary formation
		for (sq=1; sq<=N_LIMIT-N_LIMIT%3-2; sq+=3) { // loop through objects, 3 at a time
			dt = SaveDt;
			form_binary = 0; // reset this to zero; later we decide whether to form a binary, and if so, set form_binary=1
			// Local density about star k1, nearest 20 stars (10 inside, 10 outside)
			n_local = calc_n_local(k1, 10, N_LIMIT);
			// Sort stars by mass (k1 is most massive)
			sort_three_masses(sq, &k1, &k2, &k3);

			// If density above threshold, check for 3bb formation
			if (n_local > n_threshold) {
				// Are all stars singles? If not, exit loop - don't do binary formation
				if (star[k1].binind == 0 && star[k2].binind == 0 && star[k3].binind == 0) {
					triplet_count ++;
					// Calc local velocity dispersion, nearest 20 stars	
					calc_sigma_local(k1, 10, N_LIMIT, &ave_local_mass, &sigma_local);
					// Average relative speed for a Maxwellian, from Binney & Tremaine
					vrel_ave = 4.0 * sigma_local / sqrt(3.0 * PI);
					// Quantities needed for encounter
					calc_3bb_encounter_dyns(k1, k2, k3, angle1, angle2, v1, v2, v3, &vrel12, &vrel3, rng);

					//====================================================================
					//  Calculate probability of binary formation between stars 1 and 2  |
					//====================================================================
					/* units in rate are 1/T_cross  */

					// Calculate RATE of binary formation

					// Below is rate_3bb with all the velocity terms vrel_3 and vrel_12 replaced with the averaged local relative velocity, vrel_ave. We did this because we were finding that when we used the actual relative velocities, vrel12, if too large, the 3bb rate would be extremely low and binaries would not form (since it depends strongly on v: v^-9). When we replaced vrel12 with the average relative velocity (over 20 stars), the 3bb formation rate was high enough that binaries would form. We decided to use the average relative velocity for all relative velocity terms, vrel_3 and vrel_12.

					// *Note* Factor of 0.5 in front of rate_3bb ensures that our sampling method produces the correct overall analytic 3bb rate
					rate_3bb = 0.5 * sqrt(2) * sqr(PI) * sqr(n_local) *  pow(vrel_ave, -9) * pow(((star[k1].m + star[k2].m) * madhoc), 5.0) * pow(eta_min, -5.5) * (1.0 + 2.0*eta_min) * (1.0 + 2.0 * ((star[k1].m + star[k2].m + star[k3].m) / (star[k1].m + star[k2].m)) * eta_min);

					// Calculate PROBABILITY of binary formation
					P_3bb = rate_3bb * (dt * ((double) clus.N_STAR)/log(GAMMA*((double) clus.N_STAR)));

					// print info on probability calculation, for innermost 200 triplets, in each timestep
					if (triplet_count <= num_triplets_averaged) {
						P_ave += P_3bb;
						rate_ave += rate_3bb;
					}

					// For each timestep, output average rate and probability of 3bb formation, 
					// for innermost triplets (actual number is num_triplets_averaged)
					if (triplet_count == num_triplets_averaged) {						P_ave = P_ave/num_triplets_averaged;
						rate_ave = rate_ave/num_triplets_averaged;
						fprintf(threebbprobabilityfile, "%g %g %g %g %g %g\n", TotalTime, dt, dt * ((double) clus.N_STAR)/log(GAMMA*((double) clus.N_STAR)), rate_ave, P_ave, star[k3].r);
					}

					//=======================================================
					//  Monte Carlo - To form binary or not to form binary  |
					//=======================================================

					Y1 = rng_t113_dbl();
					if (P_3bb > Y1) { // Binary should be formed
						//  TODO: should really check if a three-body induced collision would happen - simply check rp to see if stars would be in contact - if so, make them collide instead.

						/* For now, can choose whether to allow any star types 
							to be involved in 3bb formation, or only BHs */
						if (ONLY_FORM_BH_THREEBODYBINARIES) { // let only BHs be involved in 3bb
							if (star[k1].se_k==14 && star[k2].se_k==14 && star[k3].se_k==14) { // ALL objects BHs - form binary
								form_binary = 1; // we will actually form it later
								N3bbformed ++;
							} else { // at least one object is not a BH; keep track of this 'light collision'
								form_binary = 0; // BUT DO NOT form a binary
							}

						/* Allow all star types to be involved with 3bb formation. 
							This is fine in point mass approximation, where we 
							don't have to worry about physical collisions */
						} else {  
							form_binary = 1;
							N3bbformed ++;
						}
						/* Here is where new binary properties are calculated.
						if form_binary=1, the new binary will be created, but if form_binary=0,
						the interaction will be logged in lightcollision log  */
						make_threebodybinary(P_3bb, k1, k2, k3, form_binary, eta_min, ave_local_mass, n_local, sigma_local, v1, v2, v3, vrel12, vrel3, delta_E_3bb, rng);
					} else { // Probability of 3bb formation too low ==> No binary formed
						form_binary = 0;
						/* do nothing; check 3bb formation for next 3 stars */
					}
				} 
			}
		} 
	} 
			


/***********************************************/	


	/* the big loop, with limits chosen so that we omit the last star if it is not paired */

	/* Meagan: 2/10/12  ** Switched 'for' loop to 'while' loop
	Had to change the structure of this loop from one that simply loops over pairs of stars
	and lets them interact, to one where it is possible to skip certain stars, if they have 
	already interacted in the 3bb formation loop above. The reasoning is that, if stars have 
	been involved in strong 3-body interaction, then they've already undergone a relaxation
	interaction, and should not be relaxed again. 
	****  To handle this, switched from 'for' loop to 'while' loop   */

	si = 1;

	// previously: for (si=1; si<=N_LIMIT-N_LIMIT%2-1; si+=2) {
	while (si<=N_LIMIT-N_LIMIT%2-1) {
		/* si is used to iterate over objects, and k, kp are the objects that will interact
		   NOTE: objects k, kp will not always be nearest neighbors, since some stars
			are skipped if they already interacted in 3bb loop!  */
		dt = SaveDt;
		// only let those stars that did not participate in 3bb formation interact/relax
		while (star[si].threebb_interacted == 1) {
			si += 1; // iterate until non-interacted object found
		}

		k = si;  // object 1 for interaction
//		kp = si + 1;
		si += 1;
		while (star[si].threebb_interacted == 1) {
		//	fprintf(threebbfile, "star marked as threebb_interacted!\n");
			si += 1; // iterate until non-interacted object found
		}	
		kp = si; // object 2 for interaction
		si += 1; // iterate for the next interaction

		/* The indices for the 2 stars that will interact are k and kp  */

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
