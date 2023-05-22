/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap.h"

/**
* @brief Does single single collision
*
* @param k index of first star
* @param kp index of second star
* @param rperimax ?
* @param w[4] ?
* @param W ?
* @param rcm ?
* @param vcm[4] ?
* @param rng gsl rng
*/
void sscollision_do(long k, long kp, double rperimax, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng)
{
	//int ST_tide; //PDK addition for hrdiag, will eventually make it an input from cmc_stelar_evolution.c
	long knew;
	double vs[20], bmax, b, rperi, Eorbnew, acoll, ecoll, ace, ece, anew, enew, efinal, afinal;
	double clight5;
	double aj, aj_k, aj_kp, tm, tn, tscls[20], lums[10], GB[10], k2;
	double Einit;
	double mass_k, mass_kp, phi_k, phi_kp, r_k, r_kp;
	double collisions_multiple; //used with BHNS_TDE parameter
        double collisions_multiple_hold;
        int MS_vanish_flag = 0; //Shi: use with giant collision

	int g_knew;
	int g_k = get_global_idx(k);
	int g_kp = get_global_idx(kp);
	mass_k = star_m[g_k];
	mass_kp = star_m[g_kp];
	phi_k = star_phi[g_k];
	phi_kp = star_phi[g_kp];
	r_k = star_r[g_k];
	r_kp = star_r[g_kp];

/* PDK addition for hrdiag... For now these are placed in by hand here and must match values given in evolv2.f.
   Will update later to be inputs within cmc_stellar_evolution.c. This will mean they must become a COMMON. */
        /*ST_tide = 0;*/

	bmax = rperimax * sqrt(1.0 + 2.0 * ((mass_k + mass_kp) * madhoc) / (rperimax * sqr(W)));
	//b = sqrt(rng_t113_dbl()) * bmax;
	b = sqrt(rng_t113_dbl_new(curr_st)) * bmax;

	rperi = madhoc*(mass_k+mass_kp)/sqr(W) * (-1.0+sqrt(1.0+sqr(b*W*W/(madhoc*mass_k+madhoc*mass_kp))));

	/* fprintf(stderr, "\n *** sscollision: rperimax=%g (%g RSUN) bmax=%g (%g RSUN) b=%g (%g RSUN) rperi=%g (%g RSUN)\n", 
	   rperimax, rperimax * units.l / RSUN, bmax, bmax * units.l / RSUN, b, b*units.l/RSUN, rperi, rperi*units.l/RSUN); */

	if (BHNS_TDE) {
		if (star[kp].se_k >= 13 && (star[k].se_k <= 1 || star[k].se_k == 7) && mass_kp >= mass_k) {
			if (mass_k * units.mstar / FB_CONST_MSUN < 0.001) {
				collisions_multiple = pow(mass_kp/(0.001*FB_CONST_MSUN/units.mstar),1./3.);
			} else {
				collisions_multiple = pow(mass_kp/mass_k,1./3.);
			}
		} else if (star[k].se_k >= 13 && (star[kp].se_k <= 1 || star[kp].se_k == 7) && mass_k >= mass_kp) {
			if (mass_kp * units.mstar / FB_CONST_MSUN < 0.001) {
				collisions_multiple = pow(mass_k/(0.001*FB_CONST_MSUN/units.mstar),1./3.);
			} else {
				collisions_multiple = pow(mass_k/mass_kp, 1./3.);
			}
		} else {
			collisions_multiple = COLL_FACTOR;  /* this is just the direct collision limit*/
		}
	} else {
		collisions_multiple = COLL_FACTOR;
	}

	if (TIDAL_CAPTURE && (star[k].se_k <= 1 || star[k].se_k >= 10 || star[k].se_k == 7) && (star[kp].se_k >= 2 && star[kp].se_k <= 9 && star[kp].se_k != 7) && (rperi <= 1.3*star[kp].rad && star[kp].rad > 2*star[k].rad)) {
                collisions_multiple_hold = rperi/star[kp].rad;
		/* log stuff */
                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_kp*madhoc, mass_k*madhoc, star[kp].se_mc*MSUN/units.m, star[kp].se_radius*RSUN/units.l, W);
		parafprintf(tidalcapturefile, "%.6g SS_COLL_Giant %s+%s+(%g,%g,%g)+%g+%g+(-100,%g,-100,%g)->", TotalTime, 
			sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), 
                        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN,
                        W*units.l/units.t/1.e5, rcm*units.l/PARSEC, star[kp].se_mc, star[kp].se_rc);

		/* instead of a merger, form a CV, WD-WD binary, or UCXB from the Ivanova & Lombardi collision mechanism */
                /* Shi: The prescription is referred to the equations in Ivanova et al. 2006, eqn. (5) and (6)*/
		ecoll = 0.88 - rperi/(3.0*star[kp].rad);
		acoll = rperi/(3.3*(1.0-sqr(ecoll)));
		
		ace = coll_CE(mass_kp*madhoc, mass_k*madhoc, star[kp].se_mc*MSUN/units.m, star[kp].se_radius*RSUN/units.l, W);
		ece = 0.0;

		if (ace < acoll) {
			afinal = ace;
			efinal = ece;
		} else {
			afinal = acoll;
			efinal = ecoll;
		}

                /* Decide if there is Roche Lobe Overflow between the star and the giant core. 
                   If not, eject envelope and form binary. 
                   If yes, merge the giant and the star (keep envelope) */
                if (bse_rl(mass_k/(star[kp].se_mc* MSUN / units.mstar))*afinal > star[k].rad && bse_rl((star[kp].se_mc* MSUN / units.mstar)/mass_k)*afinal > star[kp].se_rc*RSUN/units.l) { 
		
		        /* strip envelope of RG */
		        Einit = 0.5 * mass_k * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
		                0.5 * mass_kp * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) + 
			        0.5 * mass_k * madhoc * phi_k + 0.5 * mass_kp * madhoc * phi_kp;
		        DMse += mass_kp * madhoc;

		        aj = star[kp].se_tphys - star[kp].se_epoch;
		        star[kp].se_mt = star[kp].se_mc;
		        bse_star(&(star[kp].se_k), &(star[kp].se_mass), &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars);
		        bse_hrdiag(&(star[kp].se_mass), &aj, &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars,
			           &(star[kp].se_radius), &(star[kp].se_lum), &(star[kp].se_k), &(star[kp].se_mc), 
                                   &(star[kp].se_rc), &(star[kp].se_menv), &(star[kp].se_renv), &k2, &(star[kp].se_bhspin));
		        star[kp].se_epoch = star[kp].se_tphys - aj;
		        star[kp].rad = star[kp].se_radius * RSUN / units.l;
		        mass_kp = star[kp].se_mt * MSUN / units.mstar;
		        DMse -= mass_kp * madhoc;

		        /* check to see if MS star overfills RL at pericenter, and destroy if this is the case */
                        /* Shi: this is turned off by the MS_vanish_flag */
		        if (MS_vanish_flag && star[k].se_k <= 1 && star[k].se_radius*RSUN/units.l >= bse_rl(mass_k/mass_kp)*(1.0-efinal)*afinal) {
			/* keep only RG core as single star */
                    star_r[g_kp] = rcm;
			        star[kp].vr = vcm[3];
			        star[kp].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                    star_phi[g_kp] = potential(star_r[g_kp]);
                    phi_kp = star_phi[g_kp];
			        set_star_EJ(kp);
			        set_star_news(kp);
			        set_star_olds(kp);
			
			        /* mark stars as interacted so they don't undergo E_CONS mode stuff */
			        star[kp].interacted = 1;
			        star[kp].Eint = star[k].Eint + star[kp].Eint + Einit 
				        - 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt)) 
				        - 0.5 * mass_kp * madhoc * phi_kp;
		
			        /* log stuff */
			        parafprintf(tidalcapturefile, "%s\n", sprint_star_dyn(kp, dummystring));
                                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_kp*madhoc, mass_k*madhoc, star[kp].se_mc*MSUN/units.m, star[kp].se_radius*RSUN/units.l, W);

			        destroy_obj(k);
		        } else {
			        /* form compact binary with stripped RG core*/
			        /* put new binary together and destroy original stars */
			        knew = create_binary(k, 0);

			        g_knew = get_global_idx(knew);
			        star_m[g_knew] = mass_k + mass_kp;
			        star_r[g_knew] = rcm;
			        star_phi[g_knew] = potential(star_r[g_knew]);

			        star[knew].vr = vcm[3];
			        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			        set_star_EJ(knew);
			        set_star_news(knew);
			        set_star_olds(knew);

			        star[knew].interacted = 1;
			
			        binary[star[knew].binind].a = afinal;
			        binary[star[knew].binind].e = efinal;
			        binary[star[knew].binind].m1 = mass_k;
			        binary[star[knew].binind].m2 = mass_kp;
			        binary[star[knew].binind].rad1 = star[k].rad;
			        binary[star[knew].binind].rad2 = star[kp].rad;
			        binary[star[knew].binind].Eint1 = star[k].Eint;
			        binary[star[knew].binind].Eint2 = star[kp].Eint;

			        /* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
			        binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
			        binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);



			        binary[star[knew].binind].id1 = star[k].id;
			        binary[star[knew].binind].id2 = star[kp].id;
			        cp_SEvars_to_newbinary(k, -1, knew, 0);
			        cp_SEvars_to_newbinary(kp, -1, knew, 1);
			        binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			        compress_binary(&star[knew], &binary[star[knew].binind]);
			
			        /* log stuff */
			        parafprintf(tidalcapturefile, "%s+(%g,%g)\n", sprint_bin_dyn(knew, dummystring), 
                                        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN);
                        
                                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_kp*madhoc, mass_k*madhoc, star[kp].se_mc*MSUN/units.m, star[kp].se_radius*RSUN/units.l, W);

			        destroy_obj(k);
			        destroy_obj(kp);
                        }
                
                } else if (bse_rl(mass_k/(star[kp].se_mc* MSUN / units.mstar))*afinal <= star[k].rad || bse_rl((star[kp].se_mc* MSUN / units.mstar)/mass_k)*afinal <= star[kp].se_rc*RSUN/units.l) {

                        /* log stuff */
                        parafprintf(tidalcapturefile, "merged\n");

                        /* create new star */
                        knew = create_star(k, 0);

                        /* merge parent stars, setting mass, stellar radius, and SE params */
                        //MPI: Since we pass the star pointer itself into the merging routine, we need to copy the duplicated array values back into the star element before passing it in.
                        copy_globals_to_locals(k);
                        copy_globals_to_locals(kp);
                        merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs, curr_st);
                                
                        g_knew = get_global_idx(knew);
                        star_r[g_knew] = rcm;
                        star_m[g_knew] = star[knew].m;
                        star[knew].vr = vcm[3];
                        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                        star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                        vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
                                
                        star_phi[g_knew] = potential(star_r[g_knew]);
                        set_star_EJ(knew);
                        set_star_news(knew);
                        set_star_olds(knew);

                        /* mark stars as interacted so they don't undergo E_CONS mode stuff */
                        //star[knew].id = star_get_id_new();

                        star[knew].id = star_get_merger_id_new(star[k].id, star[kp].id);
                        star[knew].interacted = 1;

                        star[knew].Eint = star[k].Eint + star[kp].Eint
                                + 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt))
                                + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
                                - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                + 0.5 * mass_k * madhoc * phi_k
                                + 0.5 * mass_kp * madhoc * phi_kp
                                - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew];


                        /* log collision */
			/* Elena: changing format of this output*/

                        parafprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d rad1[RSUN]=%g rad2[RSUN]=%g b[RSUN]=%g vinf[km/s]=%g rperi=%g coll_mult=%g\n",
                                TotalTime,
                                star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN,
                                star[k].id, mass_k * units.mstar / FB_CONST_MSUN,
                                star[kp].id, mass_kp * units.mstar / FB_CONST_MSUN,
                                star_r[get_global_idx(knew)], star[knew].se_k, star[k].se_k, star[kp].se_k,
				star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,
                    b*units.l/RSUN,W*units.l/units.t/1.e5, rperi*units.l/RSUN, collisions_multiple_hold);
                    
			/* units should be okay already */
			double rho0_c = (star[k].se_mc  ) / ((4/3)* PI * pow((star[k].se_rc  ),3));
			double rho1_c = (star[kp].se_mc ) / ((4/3)* PI * pow((star[kp].se_rc ),3));
			double rho0_env = (star[k].se_menv  ) / ((4/3)* PI * pow((star[k].se_renv ),3));
			double rho1_env = (star[kp].se_menv  ) / ((4/3)* PI *pow((star[kp].se_renv  ),3));
			double rhor_c = (star[knew].se_mc ) / ((4/3)* PI *pow((star[knew].se_rc ),3));
			double rhor_env = (star[knew].se_menv ) / ((4/3)* PI * pow((star[knew].se_renv ),3));

			if(isnan(rho0_c)){rho0_c = -100;}
			if(isnan(rho1_c)){rho1_c = -100;}
			if(isnan(rhor_c)){rhor_c = -100;}
			if(isnan(rho0_env)){rho0_env = -100;}
			if(isnan(rho1_env)){rho1_env = -100;}
			if(isnan(rhor_env)){rhor_env = -100;}
			
			// Elena: For some stars, COSMIC assigns default renv and menv values of e-10, which makes my densities exactly 3.1831e-19. I 				will change these vales to output a -100 intead, since it is not physical //
	
			if(rho0_env >= 1.0e19){rho0_env = -100;}
			if(rho1_env >= 1.0e19){rho1_env = -100;}
			if(rhor_env >= 1.0e19){rhor_env = -100;}
					
			/*Elena: Creating a file with additional collision information */
			parafprintf(morecollfile, "%g single-single %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, star[k].id, star[kp].id, mass_k * units.mstar / FB_CONST_MSUN, mass_kp * units.mstar / FB_CONST_MSUN,
				    star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,rho0_c,rho1_c,rho0_env, rho1_env, star[k].se_k, star[kp].se_k,
				    star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN, star[knew].rad*units.l/RSUN, rhor_c, 				             rhor_env, star[knew].se_k, W*units.l/units.t/1.e5, rperi*units.l/RSUN);
                        /* destroy two progenitors */
                        destroy_obj(k);
                        destroy_obj(kp);
                }   

	} else if (TIDAL_CAPTURE && (star[kp].se_k <= 1 || star[kp].se_k >= 10 || star[kp].se_k == 7) && (star[k].se_k >= 2 && star[k].se_k <= 9 && star[k].se_k != 7) && (rperi <= 1.3 * star[k].rad && star[k].rad > 2*star[kp].rad)) {
                collisions_multiple_hold = rperi/star[k].rad;
		/* log stuff */
                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_k*madhoc, mass_kp*madhoc, star[k].se_mc*MSUN/units.m, star[k].se_radius*RSUN/units.l, W);
		parafprintf(tidalcapturefile, "%.6g SS_COLL_Giant %s+%s+(%g,%g,%g)+%g+%g+(%g,-100,%g,-100)->", TotalTime, 
				sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), 
                                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN, 
                                W*units.l/units.t/1.e5, rcm*units.l/PARSEC, star[k].se_mc, star[k].se_rc);

		/* instead of a merger, form a CV, WD-WD binary, or UCXB from the Ivanova & Lombardi collision mechanism */
		ecoll = 0.88 - rperi/(3.0*star[k].rad);
		acoll = rperi/(3.3*(1.0-sqr(ecoll)));
		
		ace = coll_CE(mass_k*madhoc, mass_kp*madhoc, star[k].se_mc*MSUN/units.m, star[k].se_radius*RSUN/units.l, W);
		ece = 0.0;

		if (ace < acoll) {
			afinal = ace;
			efinal = ece;
		} else {
			afinal = acoll;
			efinal = ecoll;
		}
		
                /* Form binary */
                if (bse_rl(mass_kp/(star[k].se_mc* MSUN / units.mstar))*afinal > star[kp].rad && bse_rl((star[k].se_mc* MSUN / units.mstar)/mass_kp)*afinal > star[k].se_rc*RSUN/units.l) {

		        /* strip envelope of RG */
		        Einit = 0.5 * mass_k * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
			        0.5 * mass_kp * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) + 
			        0.5 * mass_k * madhoc * phi_k + 0.5 * mass_kp * madhoc * phi_kp;
		        DMse += mass_k * madhoc;
		        aj = star[k].se_tphys - star[k].se_epoch;
		        star[k].se_mt = star[k].se_mc;
		        bse_star(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars);
		        bse_hrdiag(&(star[k].se_mass), &aj, &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars,
			           &(star[k].se_radius), &(star[k].se_lum), &(star[k].se_k), &(star[k].se_mc), 
                                   &(star[k].se_rc), &(star[k].se_menv), &(star[k].se_renv), &k2, &(star[k].se_bhspin));
		        star[k].se_epoch = star[k].se_tphys - aj;
		        star[k].rad = star[k].se_radius * RSUN / units.l;
		        mass_k = star[k].se_mt * MSUN / units.mstar;
		        DMse -= mass_k * madhoc;

		        /* check to see if MS star overfills RL at pericenter, and destroy if this is the case */
                        /* This is turned off by the MS_vanish_flag */
		        if (MS_vanish_flag && star[kp].se_k <= 1 && star[kp].se_radius*RSUN/units.l >= bse_rl(mass_kp/mass_k)*(1.0-efinal)*afinal) {
			        /* keep only RG core as single star */
                                star_r[g_k] = rcm;
                                star[k].vr = vcm[3];
                                star[k].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                                star_phi[g_k] = potential(star_r[g_k]);
                                phi_k = star_phi[g_k];
			        set_star_EJ(k);
			        set_star_news(k);
			        set_star_olds(k);

			
			        /* mark stars as interacted so they don't undergo E_CONS mode stuff */
			        star[k].interacted = 1;
			        star[k].Eint = star[k].Eint + star[kp].Eint + Einit 
				        - 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				        - 0.5 * mass_k * madhoc * phi_k;
	
			        /* log stuff */
			        parafprintf(tidalcapturefile, "%s\n", sprint_star_dyn(k, dummystring));
                                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_k*madhoc, mass_kp*madhoc, star[k].se_mc*MSUN/units.m, star[k].se_radius*RSUN/units.l, W);

			        destroy_obj(kp);
		        } else {
			        /* form compact binary with stripped RG core*/
			        /* put new binary together and destroy original stars */
			        knew = create_binary(k, 0);

			        g_knew = get_global_idx(knew);
			        star_m[g_knew] = mass_k + mass_kp;
			        star_r[g_knew] = rcm;
			        star_phi[g_knew] = potential(star_r[g_knew]);

			        star[knew].vr = vcm[3];
			        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			        set_star_EJ(knew);
			        set_star_news(knew);
			        set_star_olds(knew);
			        star[knew].interacted = 1;
			
			        binary[star[knew].binind].a = afinal;
			        binary[star[knew].binind].e = efinal;
			        binary[star[knew].binind].m1 = mass_k;
			        binary[star[knew].binind].m2 = mass_kp;
			        binary[star[knew].binind].rad1 = star[k].rad;
			        binary[star[knew].binind].rad2 = star[kp].rad;
			        binary[star[knew].binind].Eint1 = star[k].Eint;
			        binary[star[knew].binind].Eint2 = star[kp].Eint;

			        /* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */

			        binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
			        binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);


			        binary[star[knew].binind].id1 = star[k].id;
			        binary[star[knew].binind].id2 = star[kp].id;
			        cp_SEvars_to_newbinary(k, -1, knew, 0);
			        cp_SEvars_to_newbinary(kp, -1, knew, 1);
			        binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			        compress_binary(&star[knew], &binary[star[knew].binind]);
			
			        /* log stuff */
			        parafprintf(tidalcapturefile, "%s+(%g,%g)\n", sprint_bin_dyn(knew, dummystring), 
                                        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN);
                                //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_k*madhoc, mass_kp*madhoc, star[k].se_mc*MSUN/units.m, star[k].se_radius*RSUN/units.l, W);

			        destroy_obj(k);
			        destroy_obj(kp);
                        }
                } else if (bse_rl(mass_kp/(star[k].se_mc* MSUN / units.mstar))*afinal <= star[kp].rad || bse_rl((star[k].se_mc* MSUN / units.mstar)/mass_kp)*afinal <= star[k].se_rc*RSUN/units.l) {        
                        /* log stuff */
                        parafprintf(tidalcapturefile, "merged\n");

                        /* create new star */
                        knew = create_star(k, 0);

                        /* merge parent stars, setting mass, stellar radius, and SE params */
                        //MPI: Since we pass the star pointer itself into the merging routine, we need to copy the duplicated array values back into the star element before passing it in.
                        copy_globals_to_locals(k);
                        copy_globals_to_locals(kp);
                        merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs, curr_st);

                        g_knew = get_global_idx(knew);
                        star_r[g_knew] = rcm;
                        star_m[g_knew] = star[knew].m;
                        star[knew].vr = vcm[3];
                        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                        star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                        vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);

                        star_phi[g_knew] = potential(star_r[g_knew]);
                        set_star_EJ(knew);
                        set_star_news(knew);
                        set_star_olds(knew);

                        /* mark stars as interacted so they don't undergo E_CONS mode stuff */
                        //star[knew].id = star_get_id_new();

                        star[knew].id = star_get_merger_id_new(star[k].id, star[kp].id);
                        star[knew].interacted = 1;

                        star[knew].Eint = star[k].Eint + star[kp].Eint
                                + 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt))
                                + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
                                - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                + 0.5 * mass_k * madhoc * phi_k
                                + 0.5 * mass_kp * madhoc * phi_kp
                                - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew];


                        /* log collision */
                        parafprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d rad1[RSUN]=%g rad2[RSUN]=%g b[RSUN]=%g vinf[km/s]=%g rperi=%g coll_mult=%g\n",
                                TotalTime,
                                star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN,
                                star[k].id, mass_k * units.mstar / FB_CONST_MSUN,
                                star[kp].id, mass_kp * units.mstar / FB_CONST_MSUN,
                                star_r[get_global_idx(knew)], star[knew].se_k, star[k].se_k, star[kp].se_k,
			        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,
                    b*units.l/RSUN,W*units.l/units.t/1.e5,  rperi*units.l/RSUN, collisions_multiple_hold);

			/*Elena: Creating a file with additional collision information */
			double rho0_c = (star[k].se_mc ) / ((4/3)* PI * pow((star[k].se_rc ),3));
			double rho1_c = (star[kp].se_mc ) / ((4/3)* PI * pow((star[kp].se_rc ),3));
			double rho0_env = (star[k].se_menv ) / ((4/3)* PI * pow((star[k].se_renv ),3));
			double rho1_env = (star[kp].se_menv ) / ((4/3)* PI *pow((star[kp].se_renv ),3));
			double rhor_c = (star[knew].se_mc ) / ((4/3)* PI *pow((star[knew].se_rc ),3));
			double rhor_env =(star[knew].se_menv ) /((4/3)* PI * pow((star[knew].se_renv ),3));
			
			
			if(isnan(rho0_c)){rho0_c = -100;}
			if(isnan(rho1_c)){rho1_c = -100;}
			if(isnan(rhor_c)){rhor_c = -100;}
			if(isnan(rho0_env)){rho0_env = -100;}
			if(isnan(rho1_env)){rho1_env = -100;}
			if(isnan(rhor_env)){rhor_env = -100;}
			
			// Elena: For some stars, COSMIC assigns default renv and menv values of of e-10, which makes my densities exactly 3.1831e-19. I 				will change these vales to output a -100 intead, since it is not physical //
	

			if(rho0_env >= 1.0e19){rho0_env = -100;}
			if(rho1_env >= 1.0e19){rho1_env = -100;}
			if(rhor_env >= 1.0e19){rhor_env = -100;}
	
			parafprintf(morecollfile, "%g single-single %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, star[k].id, star[kp].id, mass_k * units.mstar / FB_CONST_MSUN, mass_kp * units.mstar / FB_CONST_MSUN,
				    star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,rho0_c,rho1_c,rho0_env, rho1_env, star[k].se_k, star[kp].se_k,
				    star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN, star[knew].rad*units.l/RSUN, rhor_c, 					    rhor_env, star[knew].se_k, W*units.l/units.t/1.e5, rperi*units.l/RSUN);

                        /* destroy two progenitors */
                        destroy_obj(k);
                        destroy_obj(kp);
		}

        /* Shi: Special treatment for giant-giant collision. */
        } else if (TIDAL_CAPTURE && (2 <= star[k].se_k && star[k].se_k <= 9 && star[k].se_k != 7) && (2 <= star[kp].se_k && star[kp].se_k <= 9 && star[kp].se_k != 7) && rperi <= 1.3*(star[k].rad+star[kp].rad)){
                collisions_multiple_hold = rperi/(star[k].rad+star[kp].rad);
                /*log stuff*/
                parafprintf(tidalcapturefile, "%.6g SS_COLL_Giant %s+%s+(%g,%g,%g)+%g+%g+(%g,%g,%g,%g)->", TotalTime,
                                sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), 
                                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN, 
                                W*units.l/units.t/1.e5, rcm*units.l/PARSEC, star[k].se_mc, star[kp].se_mc, star[k].se_rc, star[kp].se_rc);

                /* instead of a merger, form a WD-WD binary */
                afinal = coll_CE_twogiant(mass_k*madhoc, mass_kp*madhoc, star[k].se_mc*MSUN/units.m, star[kp].se_mc*MSUN/units.m, star[k].se_radius*RSUN/units.l, star[kp].se_radius*RSUN/units.l, W);
                efinal = 0.0;

                if (bse_rl(star[kp].se_mc/star[k].se_mc)*afinal > star[kp].se_rc*RSUN/units.l && bse_rl(star[k].se_mc/star[kp].se_mc)*afinal > star[k].se_rc*RSUN/units.l) {
                        /* strip envelope of both giants */
                        Einit = 0.5 * mass_k * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
                                0.5 * mass_kp * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) +
                                0.5 * mass_k * madhoc * phi_k + 0.5 * mass_kp * madhoc * phi_kp;

                        DMse += mass_k * madhoc;
                        DMse += mass_kp * madhoc;
                        aj_k = star[k].se_tphys - star[k].se_epoch;
                        star[k].se_mt = star[k].se_mc;
                        bse_star(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars);
                        bse_hrdiag(&(star[k].se_mass), &aj_k, &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars,
                                   &(star[k].se_radius), &(star[k].se_lum), &(star[k].se_k), &(star[k].se_mc), &(star[k].se_rc),
                                   &(star[k].se_menv), &(star[k].se_renv), &k2, &(star[k].se_bhspin));
                        star[k].se_epoch = star[k].se_tphys - aj_k;
                        star[k].rad = star[k].se_radius * RSUN / units.l;
                        mass_k = star[k].se_mt * MSUN / units.mstar;

                        aj_kp = star[kp].se_tphys - star[kp].se_epoch;
                        star[kp].se_mt = star[kp].se_mc;
                        bse_star(&(star[kp].se_k), &(star[kp].se_mass), &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars);
                        bse_hrdiag(&(star[kp].se_mass), &aj_kp, &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars,
                                   &(star[kp].se_radius), &(star[kp].se_lum), &(star[kp].se_k), &(star[kp].se_mc), &(star[kp].se_rc),
                                   &(star[kp].se_menv), &(star[kp].se_renv), &k2, &(star[kp].se_bhspin));
                        star[kp].se_epoch = star[kp].se_tphys - aj_kp;
                        star[kp].rad = star[kp].se_radius * RSUN / units.l;
                        mass_kp = star[kp].se_mt * MSUN / units.mstar;

                        DMse -= mass_k * madhoc;
                        DMse -= mass_kp * madhoc;

                        /* form compact binary with stripped RG cores*/
                        /* put new binary together and destroy original stars */
                        knew = create_binary(k, 0);

                        g_knew = get_global_idx(knew);
                        star_m[g_knew] = mass_k + mass_kp;
                        star_r[g_knew] = rcm;
                        star_phi[g_knew] = potential(star_r[g_knew]);
                        star[knew].vr = vcm[3];
                        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                        set_star_EJ(knew);
                        set_star_news(knew);
                        set_star_olds(knew);

                        star[knew].interacted = 1;

                        binary[star[knew].binind].a = afinal;
                        binary[star[knew].binind].e = efinal;
                        binary[star[knew].binind].m1 = mass_k;
                        binary[star[knew].binind].m2 = mass_kp;
                        binary[star[knew].binind].rad1 = star[k].rad;
                        binary[star[knew].binind].rad2 = star[kp].rad;
                        binary[star[knew].binind].Eint1 = star[k].Eint;
                        binary[star[knew].binind].Eint2 = star[kp].Eint;

                        /* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
                        binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
                                 - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                 - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
                                 + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
                        binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
                                 - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                 - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
                                 + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);

                        binary[star[knew].binind].id1 = star[k].id;
                        binary[star[knew].binind].id2 = star[kp].id;
                        cp_SEvars_to_newbinary(k, -1, knew, 0);
                        cp_SEvars_to_newbinary(kp, -1, knew, 1);
                        binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
                        compress_binary(&star[knew], &binary[star[knew].binind]);

                        /* log stuff */
                        parafprintf(tidalcapturefile, "%s+(%g,%g)\n", sprint_bin_dyn(knew, dummystring),
                                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN);

                        //parafprintf(tidalcapturefile, "coll_CE_debug %.3g %.3g %.3g %.3g %.3g\n", mass_kp*madhoc, mass_k*madhoc, star[kp].se_mc*MSUN/units.m, star[kp].se_radius*RSUN/units.l, W);
                
                        destroy_obj(k);
                        destroy_obj(kp);

                } else if (bse_rl(star[kp].se_mc/star[k].se_mc)*afinal <= star[kp].se_rc*RSUN/units.l || bse_rl(star[k].se_mc/star[kp].se_mc)*afinal <= star[k].se_rc*RSUN/units.l) {
                        /* log stuff */
                        parafprintf(tidalcapturefile, "merged\n");

                        /* create new star */
                        knew = create_star(k, 0);

                        /* merge parent stars, setting mass, stellar radius, and SE params */
                        //MPI: Since we pass the star pointer itself into the merging routine, we need to copy the duplicated array values back into the star element before passing it in.

                        copy_globals_to_locals(k);
                        copy_globals_to_locals(kp);
                        merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs, curr_st);

                        g_knew = get_global_idx(knew);
                        star_r[g_knew] = rcm;
                        star_m[g_knew] = star[knew].m;
                        star[knew].vr = vcm[3];
                        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                        star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                        vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);

                        star_phi[g_knew] = potential(star_r[g_knew]);
                        set_star_EJ(knew);
                        set_star_news(knew);
                        set_star_olds(knew);

                        /* mark stars as interacted so they don't undergo E_CONS mode stuff */
                        //star[knew].id = star_get_id_new();

                        star[knew].id = star_get_merger_id_new(star[k].id, star[kp].id);
                        star[knew].interacted = 1;

                        star[knew].Eint = star[k].Eint + star[kp].Eint
                                + 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt))
                                + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
                                - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                + 0.5 * mass_k * madhoc * phi_k
                                + 0.5 * mass_kp * madhoc * phi_kp
                                - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew];


                        /* log collision */
                        parafprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d rad1[RSUN]=%g rad2[RSUN]=%g b[RSUN]=%g vinf[km/s]=%g rperi=%g coll_mult=%g\n",
                                TotalTime,
                                star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN,
                                star[k].id, mass_k * units.mstar / FB_CONST_MSUN,
                                star[kp].id, mass_kp * units.mstar / FB_CONST_MSUN,
                                star_r[get_global_idx(knew)], star[knew].se_k, star[k].se_k, star[kp].se_k,
				star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,
                    b*units.l/RSUN,W*units.l/units.t/1.e5, rperi*units.l/RSUN, collisions_multiple_hold);

			/*Elena: Creating a file with additional collision information */
			double rho0_c = (star[k].se_mc ) / ((4/3)* PI * pow((star[k].se_rc),3));
			double rho1_c = (star[kp].se_mc ) / ((4/3)* PI * pow((star[kp].se_rc ),3));
			double rho0_env = (star[k].se_menv ) / ((4/3)* PI * pow((star[k].se_renv ),3));
			double rho1_env = (star[kp].se_menv ) / ((4/3)* PI *pow((star[kp].se_renv ),3));
			double rhor_c = (star[knew].se_mc ) / ((4/3)* PI *pow((star[knew].se_rc ),3));
			double rhor_env = (star[knew].se_menv ) / ((4/3)* PI * pow((star[knew].se_renv),3));

			if(isnan(rho0_c)){rho0_c = -100;}
			if(isnan(rho1_c)){rho1_c = -100;}
			if(isnan(rhor_c)){rhor_c = -100;}
			if(isnan(rho0_env)){rho0_env = -100;}
			if(isnan(rho1_env)){rho1_env = -100;}
			if(isnan(rhor_env)){rhor_env = -100;}
			
			// Elena: For some stars, COSMIC assigns default renv and menv values of of e-10, which makes my densities exactly 3.1831e-19. I 				will change these vales to output a -100 intead, since it is not physical //
			

			if(rho0_env >= 1.0e19){rho0_env = -100;}
			if(rho1_env >= 1.0e19){rho1_env = -100;}
			if(rhor_env >= 1.0e19){rhor_env = -100;}
			
			parafprintf(morecollfile, "%g single-single %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, star[k].id, star[kp].id, mass_k * units.mstar / FB_CONST_MSUN, mass_kp * units.mstar / FB_CONST_MSUN,
				    star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,rho0_c,rho1_c,rho0_env, rho1_env, star[k].se_k, star[kp].se_k,
				    star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN, star[knew].rad*units.l/RSUN, rhor_c, rhor_env, star[knew].se_k, W*units.l/units.t/1.e5, rperi*units.l/RSUN);

                        /* destroy two progenitors */
                        destroy_obj(k);
                        destroy_obj(kp);
                }

	} else if (BH_CAPTURE && star[k].se_k == 14 && star[kp].se_k == 14){
                
                Einit = 0.5 * mass_k * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
                        0.5 * mass_kp * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) +
                        0.5 * mass_k * madhoc * phi_k + 0.5 * mass_kp * madhoc * phi_kp;

		/* form compact binary from gravitational-wave bremsstrahlung*/
		/* put new binary together and destroy original black holes*/

		/* log stuff */
		parafprintf(tidalcapturefile, "%.6g SS_COLL_GW %s+%s+%g->", TotalTime, 
				sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), W*units.l/units.t/1.e5);
		

		/* Compute the total energy at infinity, then subtract the total 
		 * amount that is lost to GW during pericenter passage (from 
		 * Quinlan and Shapiro 1987) */
		/* Cabrera 220419: updated to use m1 != m2 expression from Quinlan and Shapiro 1989 */
                /* (85 * pi) / (12 * sqrt(2)) \approx 15.735210 */
		clight5 = pow(2.9979e10 / (units.l/units.t) ,5);

		Eorbnew = 0.5*madhoc*(mass_k*mass_kp)/(mass_k+mass_kp)*sqr(W);
                Eorbnew -= 15.735210 * pow(madhoc, 4.5) * sqr(mass_k * mass_kp) * sqrt(mass_k + mass_kp) / clight5 / pow(rperi,3.5);
        /* 85*pi/12 = 22.252948 */

		afinal = mass_k*mass_kp*madhoc*madhoc / 2. / fabs(Eorbnew); 
		efinal = 1. - rperi/afinal;

		knew = create_binary(k, 0);

		g_knew = get_global_idx(knew);
		star_m[g_knew] = mass_k + mass_kp;
		star_r[g_knew] = rcm;
		star_phi[g_knew] = potential(star_r[g_knew]);
		star[knew].vr = vcm[3];
		star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
		set_star_EJ(knew);
		set_star_news(knew);
		set_star_olds(knew);

		star[knew].interacted = 1;
		
		binary[star[knew].binind].a = afinal;
		binary[star[knew].binind].e = efinal;
		binary[star[knew].binind].m1 = mass_k;
		binary[star[knew].binind].m2 = mass_kp;
		binary[star[knew].binind].rad1 = star[k].rad;
		binary[star[knew].binind].rad2 = star[kp].rad;
		binary[star[knew].binind].Eint1 = star[k].Eint;
		binary[star[knew].binind].Eint2 = star[kp].Eint;

		/* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
		binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
			 - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
			 - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
			 + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
		binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
			 - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
			 - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
			 + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);

		binary[star[knew].binind].id1 = star[k].id;
		binary[star[knew].binind].id2 = star[kp].id;
		cp_SEvars_to_newbinary(k, -1, knew, 0);
		cp_SEvars_to_newbinary(kp, -1, knew, 1);
		binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
		compress_binary(&star[knew], &binary[star[knew].binind]);
		
		/* log stuff */
		parafprintf(tidalcapturefile, "%s\n", sprint_bin_dyn(knew, dummystring));

		destroy_obj(k);
		destroy_obj(kp);

        } else if (TC_FACTOR > 1 && (star[k].se_k != 14 || star[kp].se_k != 14) && (star[k].se_k <= 1 || star[k].se_k == 7 || star[k].se_k >= 10) && (star[kp].se_k <= 1 || star[kp].se_k == 7 || star[kp].se_k >= 10) && rperi > collisions_multiple * (star[k].rad + star[kp].rad) && rperi <= TC_FACTOR * (star[k].rad + star[kp].rad)){
                
                Einit = 0.5 * mass_k * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
                        0.5 * mass_kp * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) +
                        0.5 * mass_k * madhoc * phi_k + 0.5 * mass_kp * madhoc * phi_kp;

                /* Shi: single-single tidal capture for all stars with a radius, except giants. */
                /* put new binary together and destroy original stars */
                /* log stuff */
                parafprintf(tidalcapturefile, "%.6g SS_COLL_TC %s+%s+(%g,%g,%g)+%g+%g->", TotalTime,
                                sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), 
                                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN, W*units.l/units.t/1.e5, rcm*units.l/PARSEC);
                
                /* form a binary with zero eccentricity and semi-major axis = 2*rperi */
                afinal = 2 * rperi;
                efinal = 0.0;
               
                knew = create_binary(k, 0);

                g_knew = get_global_idx(knew);
                star_m[g_knew] = mass_k + mass_kp;
                star_r[g_knew] = rcm;
                star_phi[g_knew] = potential(star_r[g_knew]);

                star[knew].vr = vcm[3];
                star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                set_star_EJ(knew);
                set_star_news(knew);
                set_star_olds(knew);

                star[knew].interacted = 1;

                binary[star[knew].binind].a = afinal;
                binary[star[knew].binind].e = efinal;
                binary[star[knew].binind].m1 = mass_k;
                binary[star[knew].binind].m2 = mass_kp;
                binary[star[knew].binind].rad2 = star[k].rad;
                binary[star[knew].binind].rad2 = star[kp].rad;
                binary[star[knew].binind].Eint1 = star[k].Eint;
                binary[star[knew].binind].Eint2 = star[kp].Eint;
                
                /* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
                binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
                         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
                         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
                binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
                         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
                         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);

                binary[star[knew].binind].id1 = star[k].id;
                binary[star[knew].binind].id2 = star[kp].id;
                cp_SEvars_to_newbinary(k, -1, knew, 0);
                cp_SEvars_to_newbinary(kp, -1, knew, 1);
                binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
                compress_binary(&star[knew], &binary[star[knew].binind]);

                /* log stuff */
                parafprintf(tidalcapturefile, "%s+(%g,%g)\n", sprint_bin_dyn(knew, dummystring), 
                                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN);

                destroy_obj(k);
                destroy_obj(kp); 

	} else if (rperi <= collisions_multiple * (star[k].rad + star[kp].rad)) {
                /* perform standard sticky-sphere merger */
                /* If tidal capture is turned off, the cross section is just large enough to enter this clause,
                   so the next clause should never be entered. */

                /* create new star */
                knew = create_star(k, 0);

                /* merge parent stars, setting mass, stellar radius, and SE params */
        //MPI: Since we pass the star pointer itself into the merging routine, we need to copy the duplicated array values back into the star element before passing it in.
        copy_globals_to_locals(k);
        copy_globals_to_locals(kp);
                merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs, curr_st);

                g_knew = get_global_idx(knew);
                star_r[g_knew] = rcm;
                star_m[g_knew] = star[knew].m;
                star[knew].vr = vcm[3];
                star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
                
                star_phi[g_knew] = potential(star_r[g_knew]);
                set_star_EJ(knew);
                set_star_news(knew);
                set_star_olds(knew);

                /* mark stars as interacted so they don't undergo E_CONS mode stuff */
                //star[knew].id = star_get_id_new();
                star[knew].id = star_get_merger_id_new(star[k].id, star[kp].id);
                star[knew].interacted = 1;
                
                star[knew].Eint = star[k].Eint + star[kp].Eint
                        + 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt))
                        + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
                        - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                        + 0.5 * mass_k * madhoc * phi_k
                        + 0.5 * mass_kp * madhoc * phi_kp
                        - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew];


                /* log collision */
                parafprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d rad1[RSUN]=%g rad2[RSUN]=%g b[RSUN]=%g vinf[km/s]=%g rperi=%g coll_mult=%g\n",
                        TotalTime,
                        star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN,
                        star[k].id, mass_k * units.mstar / FB_CONST_MSUN,
                        star[kp].id, mass_kp * units.mstar / FB_CONST_MSUN,
                        star_r[get_global_idx(knew)], star[knew].se_k, star[k].se_k, star[kp].se_k,
	                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,
            b*units.l/RSUN,W*units.l/units.t/1.e5, rperi*units.l/RSUN, collisions_multiple);

		/*Elena: Creating a file with additional collision information */
			double rho0_c   = (star[k].se_mc)     / ((4/3)* PI * pow((star[k].se_rc),3));
			double rho1_c   = (star[kp].se_mc)    / ((4/3)* PI * pow((star[kp].se_rc),3));
			double rho0_env = (star[k].se_menv)   / ((4/3)* PI * pow((star[k].se_renv),3));
			double rho1_env = (star[kp].se_menv)  / ((4/3)* PI * pow((star[kp].se_renv),3));
			double rhor_c   = (star[knew].se_mc)  / ((4/3)* PI * pow((star[knew].se_rc),3));
			double rhor_env = (star[knew].se_menv)/ ((4/3)* PI * pow((star[knew].se_renv),3));

			if(isnan(rho0_c)){rho0_c = -100;}
			if(isnan(rho1_c)){rho1_c = -100;}
			if(isnan(rhor_c)){rhor_c = -100;}
			if(isnan(rho0_env)){rho0_env = -100;}
			if(isnan(rho1_env)){rho1_env = -100;}
			if(isnan(rhor_env)){rhor_env = -100;}
			
			// Elena: For some stars, COSMIC assigns default renv and menv values of of e-10, which makes my densities exactly 3.1831e-19. I 				will change these vales to output a -100 intead, since it is not physical //
	

			if(rho0_env >= 1.0e19){rho0_env = -100;}
			if(rho1_env >= 1.0e19){rho1_env = -100;}
			if(rhor_env >= 1.0e19){rhor_env = -100;}
				
			parafprintf(morecollfile, "%g single-single %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, star[k].id, star[kp].id, mass_k * units.mstar / FB_CONST_MSUN, mass_kp * units.mstar / FB_CONST_MSUN,
				    star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,rho0_c,rho1_c,rho0_env, rho1_env, star[k].se_k, star[kp].se_k,
				    star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN, star[knew].rad*units.l/RSUN, rhor_c, 					    rhor_env, star[knew].se_k, W*units.l/units.t/1.e5, rperi*units.l/RSUN);

                /* destroy two progenitors */
                destroy_obj(k);
                destroy_obj(kp);

	} else if (TC_POLYTROPE && (star[k].se_k != 14 || star[kp].se_k != 14) && (star[k].se_k <= 1 || star[k].se_k == 7 || star[k].se_k >= 10) && (star[kp].se_k <= 1 || star[kp].se_k == 7 || star[kp].se_k >= 10))  {
		/* apply tidal capture / common envelope test */
                /* Shi: The prescription is referring to the fitting fomulae in Kim & Lee 1999 */
                /* Shi: Don't turn this on at the same time as the TC_FACTOR! */
		Eorbnew = 0.5*madhoc*mass_k*mass_kp/(mass_k+mass_kp)*sqr(W);

		if (star[k].se_k == 1) {
			Eorbnew -= Etide(rperi, mass_k*madhoc, star[k].rad, 3.0, mass_kp*madhoc);
		} else if (star[k].se_k < 10) {
			Eorbnew -= Etide(rperi, mass_k*madhoc, star[k].rad, 1.5, mass_kp*madhoc);
		}

		if (star[kp].se_k == 1) {
			Eorbnew -= Etide(rperi, mass_kp*madhoc, star[kp].rad, 3.0, mass_k*madhoc);
		} else if (star[kp].se_k < 10) {
			Eorbnew -= Etide(rperi, mass_kp*madhoc, star[kp].rad, 1.5, mass_k*madhoc);
		}

		if (Eorbnew < 0.0) {
			/* bound system; don't worry about RL overflow here, since BSE will take care of that automatically */
			//anew = madhoc * mass_k * madhoc * mass_kp / (2.0 * fabs(Eorbnew));
			//enew = MIN(0.0, 1.0-rperi/anew); //Why min?
                        //enew = MAX(0.0, 1.0-rperi/anew);

			/* apply rapid tidal circularization of orbit, assuming no angular momentum is transferred to 
			   stars' internal rotation for simplicity; but only if there is no Roche-lobe overflow at 
			   pericenter */
			//if (bse_rl(mass_k/mass_kp)*anew*(1.0-enew) > star[k].rad && 
			//    bse_rl(mass_kp/mass_k)*anew*(1.0-enew) > star[kp].rad) {
			//	efinal = 0.0;
			//	afinal = anew * (1.0 - enew*enew);
			//} else {
			//	afinal = anew;
			//	efinal = enew;
			//}
                       
                        /*Manually set and circularize the tidal capture binary*/
                        afinal = 2 * rperi;
                        efinal = 0.0;

                        /* Check if the binary is filling its Roche Lobe */
                        /* If yes, treat it as TDE and merge the binary immediately */
                        if (bse_rl(mass_k/mass_kp)*afinal < star[k].rad || bse_rl(mass_kp/mass_k)*afinal < star[kp].rad) {
                                /* log stuff */
                                parafprintf(tidalcapturefile, "%.6g SS_COLL_TC_P %s+%s+(%g,%g,%g)+%g+%g->merged\n", 
                                        TotalTime,
                                        sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2),
                                        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN,
                                        W*units.l/units.t/1.e5, rcm*units.l/PARSEC);

                                /* create new star */
		                knew = create_star(k, 0);
                              
                                /* merge parent stars, setting mass, stellar radius, and SE params */
                        //MPI: Since we pass the star pointer itself into the merging routine, we need to copy the duplicated array values back into the star element before passing it in.
                         copy_globals_to_locals(k);
                         copy_globals_to_locals(kp);
                                merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs, curr_st);
  
                                g_knew = get_global_idx(knew);
                                star_r[g_knew] = rcm;
                                star_m[g_knew] = star[knew].m;
                                star[knew].vr = vcm[3];
                                star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
                                star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                                vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                                //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);         
        
                                star_phi[g_knew] = potential(star_r[g_knew]);
                                set_star_EJ(knew);
                                set_star_news(knew);
                                set_star_olds(knew);

                                /* mark stars as interacted so they don't undergo E_CONS mode stuff */
                                //star[knew].id = star_get_id_new();

                                star[knew].id = star_get_merger_id_new(star[k].id, star[kp].id);
                                star[knew].interacted = 1;

                                star[knew].Eint = star[k].Eint + star[kp].Eint
                                        + 0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt))
                                        + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
                                        - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
                                        + 0.5 * mass_k * madhoc * phi_k
                                        + 0.5 * mass_kp * madhoc * phi_kp
                                        - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew];

                                /* log collision */
                                parafprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d rad1[RSUN]=%g rad2[RSUN]=%g b[RSUN]=%g vinf[km/s]=%g rperi=%g coll_mult=%g\n",
                                        TotalTime,
                                        star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN,
                                        star[k].id, mass_k * units.mstar / FB_CONST_MSUN,
                                        star[kp].id, mass_kp * units.mstar / FB_CONST_MSUN,
                                        star_r[get_global_idx(knew)], star[knew].se_k, star[k].se_k, star[kp].se_k,
			                star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,
                            b*units.l/RSUN,W*units.l/units.t/1.e5, rperi*units.l/RSUN, rperi/(star[k].rad+star[kp].rad));

				/*Elena: Creating a file with additional collision information */
				
			double rho0_c   = (star[k].se_mc)      / ((4/3)* PI * pow((star[k].se_rc),3));
			double rho1_c   = (star[kp].se_mc)     / ((4/3)* PI * pow((star[kp].se_rc),3));
			double rho0_env = (star[k].se_menv)    / ((4/3)* PI * pow((star[k].se_renv),3));
			double rho1_env = (star[kp].se_menv)   / ((4/3)* PI * pow((star[kp].se_renv),3));
			double rhor_c   = (star[knew].se_mc)   / ((4/3)* PI * pow((star[knew].se_rc),3));
			double rhor_env = (star[knew].se_menv) / ((4/3)* PI * pow((star[knew].se_renv),3));
	
			if(isnan(rho0_c)){rho0_c = -100;}
			if(isnan(rho1_c)){rho1_c = -100;}
			if(isnan(rhor_c)){rhor_c = -100;}
			if(isnan(rho0_env)){rho0_env = -100;}
			if(isnan(rho1_env)){rho1_env = -100;}
			if(isnan(rhor_env)){rhor_env = -100;}
			
			// Elena: For some stars, COSMIC assigns default renv and menv values of of e-10, which makes my densities exactly 3.1831e-19. I will change these vales to output a -100 intead, since it is not physical //
	
			if(rho0_env >= 1.0e19){rho0_env = -100;}
			if(rho1_env >= 1.0e19){rho1_env = -100;}
			if(rhor_env >= 1.0e19){rhor_env = -100;}
			
			parafprintf(morecollfile, "%g single-single %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, star[k].id, star[kp].id, mass_k * units.mstar / FB_CONST_MSUN, mass_kp * units.mstar / FB_CONST_MSUN,
				    star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN,rho0_c,rho1_c,rho0_env, rho1_env, star[k].se_k, star[kp].se_k,
				    star[knew].id, star_m[get_global_idx(knew)] * units.mstar / FB_CONST_MSUN, star[knew].rad*units.l/RSUN, rhor_c, 					    rhor_env, star[knew].se_k, W*units.l/units.t/1.e5, rperi*units.l/RSUN);

                                /* destroy two progenitors */
                                destroy_obj(k);
                                destroy_obj(kp);

			} else {
                        
                                /* put new binary together and destroy original stars */
			        knew = create_binary(k, 0);

			        g_knew = get_global_idx(knew);
			        star_m[g_knew] = mass_k + mass_kp;
			        star_r[g_knew] = rcm;
			        star_phi[g_knew] = potential(star_r[g_knew]);
			        star[knew].vr = vcm[3];
			        star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			        set_star_EJ(knew);
			        set_star_news(knew);
			        set_star_olds(knew);
			        star[knew].interacted = 1;
			
			        binary[star[knew].binind].a = afinal;
			        binary[star[knew].binind].e = efinal;
			        binary[star[knew].binind].m1 = mass_k;
			        binary[star[knew].binind].m2 = mass_kp;
			        binary[star[knew].binind].rad1 = star[k].rad;
			        binary[star[knew].binind].rad2 = star[kp].rad;
			        binary[star[knew].binind].Eint1 = star[k].Eint;
			        binary[star[knew].binind].Eint2 = star[kp].Eint;

			        /* put tidal energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */

			        binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * 
				        (0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				         + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         + 0.5 * mass_k * madhoc * phi_k
				         + 0.5 * mass_kp * madhoc * phi_kp
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);
			        binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * 
				        (0.5 * mass_k * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				         + 0.5 * mass_kp * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
				         - 0.5 * star_m[g_knew] * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				         + 0.5 * mass_k * madhoc * phi_k
				         + 0.5 * mass_kp * madhoc * phi_kp
				         - 0.5 * star_m[g_knew] * madhoc * star_phi[g_knew]
				         + 0.5 * mass_k * madhoc * mass_kp * madhoc / binary[star[knew].binind].a);


			        binary[star[knew].binind].id1 = star[k].id;
			        binary[star[knew].binind].id2 = star[kp].id;
			        cp_SEvars_to_newbinary(k, -1, knew, 0);
			        cp_SEvars_to_newbinary(kp, -1, knew, 1);
			        binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			        compress_binary(&star[knew], &binary[star[knew].binind]);

			        /* log stuff */
			        parafprintf(tidalcapturefile, "%.6g SS_COLL_TC_P %s+%s+(%g,%g,%g)+%g+%g->%s+(%g,%g)\n", TotalTime, 
				        sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2), 
                                        star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN, rperi*units.l/RSUN, 
                                        W*units.l/units.t/1.e5, rcm*units.l/PARSEC, 
                                        sprint_bin_dyn(knew, dummystring3), star[k].rad*units.l/RSUN, star[kp].rad*units.l/RSUN);
			
			        destroy_obj(k);
			        destroy_obj(kp);
		//} else {
			//parafprintf(tidalcapturefile, "%.6g SS_COLL_TC_P_FAILED %s+%s->%s+%s\n", TotalTime, 
			//	sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring2),
			//	sprint_star_dyn(k, dummystring3), sprint_star_dyn(kp, dummystring4));
			}
		}
	}
}

/**
* @brief merge two stars using stellar evolution if it's enabled
*
* @param star1 pointer to star 1
* @param star2 pointer to star 2
* @param merged_star pointer to merged star
* @param vs ?
* @param s rng state
*/
void merge_two_stars(star_t *star1, star_t *star2, star_t *merged_star, double *vs, struct rng_t113_state* s) {
	double tphysf, dtp, age, temp_rad, bseaj[2];
	//double vsaddl[12], 
	double tm, tn, tscls[20], lums[10], GB[10], k2, lamb_val;
	binary_t tempbinary, tbcopy;
	int tbi=-1, j, ktry, i, fb, icase, convert;
	fb = 1;

	for(i=0;i<20;i++) {
	  vs[i] = 0.0;
	}

  bse_set_taus113state(*curr_st, 0);

	if (STELLAR_EVOLUTION && !STAR_AGING_SCHEME) {
		/* evolve for just a year for merger */
		tphysf = star1->se_tphys + 1.0e-6;
		
		/* Create temporary binary in which the two stars touch by a factor of 10 at pericenter, 
		   but underfill their Roche lobes by a factor of 10 at R=a.  This ensures that 
		   RL overflow is not invoked in BSE, and that the merger routine is called instead. */
		tempbinary.id1 = star1->id;
		tempbinary.id2 = star2->id;
		tempbinary.rad1 = star1->rad;
		tempbinary.rad2 = star2->rad;
		//tempbinary.m1 = star1->m;
		//tempbinary.m2 = star2->m;
		//MPI: This was missed before in previous version. I suppose this has to be done
		tempbinary.m1 = star1->m;
		tempbinary.m2 = star2->m;

		tempbinary.Eint1 = star1->Eint;
		tempbinary.Eint2 = star2->Eint;
		//tempbinary.a = BSE_WRAP_MAX(10.0*star1->rad/bse_rl(star1->m/star2->m), 
		//			    10.0*star2->rad/bse_rl(star2->m/star1->m));
		//tempbinary.e = 1.0 - 0.1*(star1->rad+star2->rad)/tempbinary.a;
		tempbinary.a = star1->rad + star2->rad; //We are directly merging or CEing, where a merger MUST occur. 
		//So use this for now... But be prepared to revert to previous values.
		tempbinary.e = 0.0;
		tempbinary.inuse = 1;
		tempbinary.bse_mass0[0] = star1->se_mass;
		tempbinary.bse_mass0[1] = star2->se_mass;
		tempbinary.bse_kw[0] = star1->se_k;
		tempbinary.bse_kw[1] = star2->se_k;
		tempbinary.bse_mass[0] = star1->se_mt;
		tempbinary.bse_mass[1] = star2->se_mt;
		tempbinary.bse_ospin[0] = star1->se_ospin;
		tempbinary.bse_ospin[1] = star2->se_ospin;
                tempbinary.bse_B_0[0] = star1->se_B_0; /* PK add pulsar stuff here */
                tempbinary.bse_B_0[1] = star2->se_B_0;
                tempbinary.bse_bacc[0] = star1->se_bacc;
                tempbinary.bse_bacc[1] = star2->se_bacc;
                tempbinary.bse_tacc[0] = star1->se_tacc;
                tempbinary.bse_tacc[1] = star2->se_tacc;
                tempbinary.bse_bcm_B[0] = star1->se_scm_B;
                tempbinary.bse_bcm_B[1] = star2->se_scm_B;
                tempbinary.bse_bcm_formation[0] = star1->se_scm_formation;
                tempbinary.bse_bcm_formation[1] = star2->se_scm_formation;
		tempbinary.bse_epoch[0] = star1->se_epoch;
		tempbinary.bse_epoch[1] = star2->se_epoch;
		tempbinary.bse_tphys = star1->se_tphys;
		tempbinary.bse_radius[0] = star1->se_radius;
		tempbinary.bse_radius[1] = star2->se_radius;
		tempbinary.bse_lum[0] = star1->se_lum;
		tempbinary.bse_lum[1] = star2->se_lum;
		tempbinary.bse_massc[0] = star1->se_mc;
		tempbinary.bse_massc[1] = star2->se_mc;
		tempbinary.bse_radc[0] = star1->se_rc;
		tempbinary.bse_radc[1] = star2->se_rc;
		tempbinary.bse_menv[0] = star1->se_menv;
		tempbinary.bse_menv[1] = star2->se_menv;
		tempbinary.bse_renv[0] = star1->se_renv;
		tempbinary.bse_renv[1] = star2->se_renv;
		tempbinary.bse_tms[0] = star1->se_tms;
		tempbinary.bse_tms[1] = star2->se_tms;
		tempbinary.bse_bhspin[0] = star1->se_bhspin;
		tempbinary.bse_bhspin[1] = star2->se_bhspin;
		
		tempbinary.bse_tb = sqrt(cub(tempbinary.a * units.l / AU)/(tempbinary.bse_mass[0]+tempbinary.bse_mass[1]))*365.25;
		
		tbcopy = tempbinary;

		if(tempbinary.bse_mass[0]!=0.0 && tempbinary.bse_mass[1]!=0.0) { //else system isn't a binary!!!
		/* dtp = tphysf - tempbinary.bse_tphys; */
		/* Since the evolution time is so short in this routine, we can simply set dtp=0.0
		   without worrying about the bcm arrays filling up. */
		  dtp = 0.0;
		  bse_set_id1_pass(tempbinary.id1);
		  bse_set_id2_pass(tempbinary.id2);

		  icase = icase_get(tempbinary.bse_kw[0],tempbinary.bse_kw[1]);

// make ktype callable from cmc then if icase > 100 dont use mix.f
		  if (icase <= 100 ) {
		    //
		    // Merging stars via mixing the stars
		    //
		    dprintf("Calling mix from merge_two_stars 1: id1=%ld id2=%ld m1=%g m2=%g kw1=%d kw2=%d \n",tempbinary.id1,tempbinary.id2,tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		    //aj = star[kp].se_tphys - star[kp].se_epoch;
		    bseaj[0] = tempbinary.bse_tphys - tempbinary.bse_epoch[0];
		    bseaj[1] = tempbinary.bse_tphys - tempbinary.bse_epoch[1];
		    bse_mix(&(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), &(bseaj[0]), &(tempbinary.bse_kw[0]), zpars, &(tempbinary.bse_bhspin[0]));
		    
		    bse_star(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars);
		    
		    bse_hrdiag(&(tempbinary.bse_mass0[0]), &(bseaj[0]), &(tempbinary.bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars,
			       &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_kw[0]), &(tempbinary.bse_massc[0]), &(tempbinary.bse_radc[0]), 
			       &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), &k2, &(tempbinary.bse_bhspin[0]));
		    tempbinary.bse_epoch[0] = tempbinary.bse_tphys - bseaj[0];
		    tempbinary.bse_mass0[1] = 0.0;
		    //tempbinary.bse_mass[1] = 0.0;
		    tempbinary.bse_radius[1] = 0.0;
		    dprintf("Called mix from merge_two_stars 1: id1=%ld id2=%ld m1=%g m2=%g kw1=%d kw2=%d \n",tempbinary.id1,tempbinary.id2,tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		    
		    //tempbinary.bse_epoch[1] = tempbinary.bse_tphys - bseaj[1];
		    //if NS and not an MSP reset? Also, there is no bcm array here, yet below it searches it - must do a check of icase before searching bcm...
		  } else {
		    //
		    // Merging stars via common envelope - yet forcing the stars to merge.
		    //
		  
		  /* Force merger if necessary by resetting pericenter to nearly 0;  This may be needed in some cases
		     because BSE doesn't strictly conserve angular momentum in binaries during common envelope evolution. */
		    ktry = 1;
		    while (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0 && ktry < 10) {
		      dprintf("Attempting to force merger in CE possibly by repeating evolution with a change of sep and lambda.\n");
		      if(ktry>0){
			lamb_val = -0.0001/((float)ktry);
			bse_set_lambdaf(lamb_val); //perhaps should do this in the second attempt? 
		      }
		      tempbinary.a = tbcopy.a/ktry;
		      tempbinary.bse_tb = sqrt(cub(tempbinary.a * units.l / AU)/(tempbinary.bse_mass[0]+tempbinary.bse_mass[1]))*365.25;
		      dtp = 0.0;
		      //tphysf += 1.0e-6;
		      // PDK added...			
		      tempbinary.bse_mass0[0] = star1->se_mass;
		      tempbinary.bse_mass0[1] = star2->se_mass;
		      tempbinary.bse_kw[0] = star1->se_k;
		      tempbinary.bse_kw[1] = star2->se_k;
		      tempbinary.bse_mass[0] = star1->se_mt;
		      tempbinary.bse_mass[1] = star2->se_mt;
		      tempbinary.bse_ospin[0] = star1->se_ospin;
		      tempbinary.bse_ospin[1] = star2->se_ospin;
			tempbinary.bse_B_0[0] = star1->se_B_0; // PK add pulsar stuff here
			tempbinary.bse_B_0[1] = star2->se_B_0;
			tempbinary.bse_bacc[0] = star1->se_bacc;
			tempbinary.bse_bacc[1] = star2->se_bacc;
			tempbinary.bse_tacc[0] = star1->se_tacc;
			tempbinary.bse_tacc[1] = star2->se_tacc;
			tempbinary.bse_bcm_B[0] = star1->se_scm_B;
			tempbinary.bse_bcm_B[1] = star2->se_scm_B;
			tempbinary.bse_bcm_formation[0] = star1->se_scm_formation;
			tempbinary.bse_bcm_formation[1] = star2->se_scm_formation;
		      tempbinary.bse_epoch[0] = star1->se_epoch;
		      tempbinary.bse_epoch[1] = star2->se_epoch;
		      tempbinary.bse_tphys = star1->se_tphys;
		      tempbinary.bse_radius[0] = star1->se_radius;
		      tempbinary.bse_radius[1] = star2->se_radius;
		      tempbinary.bse_lum[0] = star1->se_lum;
		      tempbinary.bse_lum[1] = star2->se_lum;
		      tempbinary.bse_massc[0] = star1->se_mc;
		      tempbinary.bse_massc[1] = star2->se_mc;
		      tempbinary.bse_radc[0] = star1->se_rc;
		      tempbinary.bse_radc[1] = star2->se_rc;
		      tempbinary.bse_menv[0] = star1->se_menv;
		      tempbinary.bse_menv[1] = star2->se_menv;
		      tempbinary.bse_renv[0] = star1->se_renv;
		      tempbinary.bse_renv[1] = star2->se_renv;
		      tempbinary.bse_tms[0] = star1->se_tms;
		      tempbinary.bse_tms[1] = star2->se_tms; 
		      tempbinary.bse_bhspin[0] = star1->se_bhspin;
		      tempbinary.bse_bhspin[1] = star2->se_bhspin; 
		      dprintf("before ce merge1: tb=%g a=%g m1=%g m2=%g e=%g kw1=%d kw2=%d r1=%g r2=%g\n",tempbinary.bse_tb,tempbinary.a,tempbinary.bse_mass[0],tempbinary.bse_mass[1],tempbinary.e,tempbinary.bse_kw[0],tempbinary.bse_kw[1], tempbinary.bse_radius[0],tempbinary.bse_radius[1]);
		      //tphysf += 1.0e-4; //PDK ADD(was1.0e-6): perhaps should add 100 years here, just to be a little safer in evolv2...
		      if(tempbinary.bse_tb==0.0){
			dprintf("tb0ev2_2: %g %g %g %g %g\n",tempbinary.bse_tb,tempbinary.a,tempbinary.bse_mass[0],tempbinary.bse_mass[1],tempbinary.e);
		      }
		      //CALL comenv routine...
		      cmc_bse_comenv(&(tempbinary), units.l, RSUN, zpars, 
				     vs, &fb);
		      dprintf("after ce merge1: tb=%g a=%g m1=%g m2=%g e=%g kw1=%d kw2=%d r1=%g r2=%g\n",tempbinary.bse_tb,tempbinary.a,tempbinary.bse_mass[0],tempbinary.bse_mass[1],tempbinary.e,tempbinary.bse_kw[0],tempbinary.bse_kw[1],tempbinary.bse_radius[0],tempbinary.bse_radius[1]);
		      ktry++;
		    }
		    bse_set_lambdaf(BSE_LAMBDAF);
		    bse_set_merger(-1.0);
		  }  //end of merger attempt, either from mix (non-CE) or evolv2 (for CE) with appropriate separation selection.

		  
		  /* make sure outcome was as expected */
		  if (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0) {
		    eprintf("Artificial stellar evolution of eccentric binary failed: both stars have non-zero mass.\n");
		    eprintf("k1=%d k2=%d kb1=%d kb2=%d\n", star1->se_k, star2->se_k, 
			    tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		    eprintf("m1=%g m2=%g mb1=%g mb2=%g\n", star1->se_mt, star2->se_mt,
			    tempbinary.bse_mass[0], tempbinary.bse_mass[1]);
		    eprintf("tphys=%g tphysf=%g (tphysf-tphys)/tphys=%g\n", tempbinary.bse_tphys, tphysf,
			    (tphysf-tempbinary.bse_tphys)/tempbinary.bse_tphys);
		    eprintf("tb=%g a=%g e=%g\n", tempbinary.bse_tb, tempbinary.a * units.l / AU, tempbinary.e);
		    eprintf("kin1=%d kin2=%d kout1=%d kout2=%d\n", tbcopy.bse_kw[0], tbcopy.bse_kw[1],
			    tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		    eprintf("e=%.18g\n", tbcopy.e);
		    eprintf("bse_tb=%.18g\n", tbcopy.bse_tb);
		    eprintf("bse_mass0=%.18g %.18g\n", tbcopy.bse_mass0[0], tbcopy.bse_mass0[1]);
		    eprintf("bse_kw=%d %d\n", tbcopy.bse_kw[0], tbcopy.bse_kw[1]);
		    eprintf("bse_mass=%.18g %.18g\n", tbcopy.bse_mass[0], tbcopy.bse_mass[1]);
		    eprintf("bse_ospin=%.18g %.18g\n", tbcopy.bse_ospin[0], tbcopy.bse_ospin[1]);
		    eprintf("bse_epoch=%.18g %.18g\n", tbcopy.bse_epoch[0], tbcopy.bse_epoch[1]);
		    eprintf("bse_tphys=%.18g\n", tbcopy.bse_tphys);
		    eprintf("bse_radius=%.18g %.18g\n", tbcopy.bse_radius[0], tbcopy.bse_radius[1]);
		    eprintf("bse_lum=%.18g %.18g\n", tbcopy.bse_lum[0], tbcopy.bse_lum[1]);
		    eprintf("bse_massc=%.18g %.18g\n", tbcopy.bse_massc[0], tbcopy.bse_massc[1]);
		    eprintf("bse_radc=%.18g %.18g\n", tbcopy.bse_radc[0], tbcopy.bse_radc[1]);
		    eprintf("bse_menv=%.18g %.18g\n", tbcopy.bse_menv[0], tbcopy.bse_menv[1]);
		    eprintf("bse_renv=%.18g %.18g\n", tbcopy.bse_renv[0], tbcopy.bse_renv[1]);
		    eprintf("bse_tms=%.18g %.18g\n", tbcopy.bse_tms[0], tbcopy.bse_tms[1]);
		    eprintf("bse_bhspin=%.18g %.18g\n", tbcopy.bse_bhspin[0], tbcopy.bse_bhspin[1]);
		    j = 1;
		    while (bse_get_bpp(j, 1) >= 0.0) {
		      fprintf(stderr, "time=%g m1=%g m2=%g k1=%d k2=%d sep=%g ecc=%g r1/rol1=%g r2/rol2=%g type=%s\n",
			      bse_get_bpp(j, 1), bse_get_bpp(j, 2), bse_get_bpp(j, 3), (int) bse_get_bpp(j, 4), 
			      (int) bse_get_bpp(j, 5), bse_get_bpp(j, 6), bse_get_bpp(j, 8), bse_get_bpp(j, 9), 
			      bse_get_bpp(j, 10), bse_get_bselabel((int) bse_get_bpp(j, 11)));
		      fflush(NULL);
		      j++;
		    }
		    icase = icase_get(tempbinary.bse_kw[0],tempbinary.bse_kw[1]);
		    // make ktype callable from cmc then if icase > 100 dont use mix.f
		    if (icase <= 100 ) {
		      //
		      // Merging stars via mixing the stars
		      //
		      dprintf("Calling mix after merger failure: id1=%ld id2=%ld m1=%g m2=%g kw1=%d kw2=%d \n",tempbinary.id1,tempbinary.id2,tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		      //aj = star[kp].se_tphys - star[kp].se_epoch;
		      bseaj[0] = tempbinary.bse_tphys - tempbinary.bse_epoch[0];
		      bseaj[1] = tempbinary.bse_tphys - tempbinary.bse_epoch[1];
		      bse_mix(&(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), &(bseaj[0]), &(tempbinary.bse_kw[0]), zpars, &(tempbinary.bse_bhspin[0]));
		      
		      bse_star(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars);
		      
		      bse_hrdiag(&(tempbinary.bse_mass0[0]), &(bseaj[0]), &(tempbinary.bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars,
				 &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_kw[0]), &(tempbinary.bse_massc[0]), &(tempbinary.bse_radc[0]), 
				 &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), &k2, &(tempbinary.bse_bhspin[0]));
		      tempbinary.bse_epoch[0] = tempbinary.bse_tphys - bseaj[0];
		      tempbinary.bse_mass0[1] = 0.0;
		      //tempbinary.bse_mass[1] = 0.0;
		      tempbinary.bse_radius[1] = 0.0;
		      dprintf("Called mix after merger failure: id1=%ld id2=%ld m1=%g m2=%g kw1=%d kw2=%d \n",tempbinary.id1,tempbinary.id2,tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		    } else {
		      exit_cleanly(1, __FUNCTION__);
		    }
		  }
		  //Make sure all info that isn't needed is zeroed out...
		  if (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0) {
		    exit_cleanly(1, __FUNCTION__);
		  } else if (tempbinary.bse_mass[0] != 0) {
		    tbi = 0;
		    tempbinary.bse_mass0[1] = 0.0;
		    tempbinary.bse_kw[1] = 0;
		    tempbinary.bse_mass[1] = 0.0;
		    tempbinary.bse_ospin[1] = 0.0;
		    tempbinary.bse_B_0[1] = 0.0;
		    tempbinary.bse_bacc[1] = 0.0;
		    tempbinary.bse_tacc[1] = 0.0;
		    tempbinary.bse_bcm_B[1] = 0.0;
		    tempbinary.bse_bcm_formation[1] = 0;
		    tempbinary.bse_epoch[1] = 0.0;
		    tempbinary.bse_radius[1] = 0.0;
		    tempbinary.bse_lum[1] = 0.0;
		    tempbinary.bse_massc[1] = 0.0;
		    tempbinary.bse_radc[1] = 0.0;
		    tempbinary.bse_menv[1] = 0.0;
		    tempbinary.bse_renv[1] = 0.0;
		    tempbinary.bse_tms[1] = 0.0;
		    tempbinary.bse_bhspin[1] = 0.0;
		  } else if (tempbinary.bse_mass[1] != 0) {
		    tbi = 1;
		    tempbinary.bse_mass0[0] = 0.0;
		    tempbinary.bse_kw[0] = 0;
		    tempbinary.bse_mass[0] = 0.0;
		    tempbinary.bse_ospin[0] = 0.0;
		    tempbinary.bse_B_0[0] = 0.0;
		    tempbinary.bse_bacc[0] = 0.0;
		    tempbinary.bse_tacc[0] = 0.0;
		    tempbinary.bse_bcm_B[0] = 0.0;
		    tempbinary.bse_bcm_formation[0] = 0;
		    tempbinary.bse_epoch[0] = 0.0;
		    tempbinary.bse_radius[0] = 0.0;
		    tempbinary.bse_lum[0] = 0.0;
		    tempbinary.bse_massc[0] = 0.0;
		    tempbinary.bse_radc[0] = 0.0;
		    tempbinary.bse_menv[0] = 0.0;
		    tempbinary.bse_renv[0] = 0.0;
		    tempbinary.bse_tms[0] = 0.0;
		    tempbinary.bse_bhspin[0] = 0.0;
		  } else {
		    if (tempbinary.bse_kw[0] == 15) {
		      tbi = 0;
		      tempbinary.bse_mass0[1] = 0.0;
		      tempbinary.bse_kw[1] = 0;
		      tempbinary.bse_mass[1] = 0.0;
		      tempbinary.bse_ospin[1] = 0.0;
		      tempbinary.bse_B_0[1] = 0.0;
		      tempbinary.bse_bacc[1] = 0.0;
		      tempbinary.bse_tacc[1] = 0.0;
		      tempbinary.bse_bcm_B[1] = 0.0;
		      tempbinary.bse_bcm_formation[1] = 0;
		      tempbinary.bse_epoch[1] = 0.0;
		      tempbinary.bse_radius[1] = 0.0;
		      tempbinary.bse_lum[1] = 0.0;
		      tempbinary.bse_massc[1] = 0.0;
		      tempbinary.bse_radc[1] = 0.0;
		      tempbinary.bse_menv[1] = 0.0;
		      tempbinary.bse_renv[1] = 0.0;
		      tempbinary.bse_tms[1] = 0.0;
		      tempbinary.bse_bhspin[1] = 0.0;
		    } else if (tempbinary.bse_kw[1] == 15) {
		      tbi = 1;
		      tempbinary.bse_mass0[0] = 0.0;
		      tempbinary.bse_kw[0] = 0;
		      tempbinary.bse_mass[0] = 0.0;
		      tempbinary.bse_ospin[0] = 0.0;
		      tempbinary.bse_B_0[0] = 0.0;
		      tempbinary.bse_bacc[0] = 0.0;
		      tempbinary.bse_tacc[0] = 0.0;
		      tempbinary.bse_bcm_B[0] = 0.0;
		      tempbinary.bse_bcm_formation[0] = 0;
		      tempbinary.bse_epoch[0] = 0.0;
		      tempbinary.bse_radius[0] = 0.0;
		      tempbinary.bse_lum[0] = 0.0;
		      tempbinary.bse_massc[0] = 0.0;
		      tempbinary.bse_radc[0] = 0.0;
		      tempbinary.bse_menv[0] = 0.0;
		      tempbinary.bse_renv[0] = 0.0;
		      tempbinary.bse_tms[0] = 0.0;
		      tempbinary.bse_bhspin[0] = 0.0;
		    } else {
		      eprintf("Artificial stellar evolution of eccentric binary failed: both stars have zero mass.\n");
		      eprintf("kin1=%d kin2=%d kout1=%d kout2=%d\n", tbcopy.bse_kw[0], tbcopy.bse_kw[1],
			      tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
		      eprintf("e=%.18g\n", tbcopy.e);
		      eprintf("bse_tb=%.18g\n", tbcopy.bse_tb);
		      eprintf("bse_mass0=%.18g %.18g\n", tbcopy.bse_mass0[0], tbcopy.bse_mass0[1]);
		      eprintf("bse_kw=%d %d\n", tbcopy.bse_kw[0], tbcopy.bse_kw[1]);
		      eprintf("bse_mass=%.18g %.18g\n", tbcopy.bse_mass[0], tbcopy.bse_mass[1]);
		      eprintf("bse_ospin=%.18g %.18g\n", tbcopy.bse_ospin[0], tbcopy.bse_ospin[1]);
		      eprintf("bse_epoch=%.18g %.18g\n", tbcopy.bse_epoch[0], tbcopy.bse_epoch[1]);
		      eprintf("bse_tphys=%.18g\n", tbcopy.bse_tphys);
		      eprintf("bse_radius=%.18g %.18g\n", tbcopy.bse_radius[0], tbcopy.bse_radius[1]);
		      eprintf("bse_lum=%.18g %.18g\n", tbcopy.bse_lum[0], tbcopy.bse_lum[1]);
		      eprintf("bse_massc=%.18g %.18g\n", tbcopy.bse_massc[0], tbcopy.bse_massc[1]);
		      eprintf("bse_radc=%.18g %.18g\n", tbcopy.bse_radc[0], tbcopy.bse_radc[1]);
		      eprintf("bse_menv=%.18g %.18g\n", tbcopy.bse_menv[0], tbcopy.bse_menv[1]);
		      eprintf("bse_renv=%.18g %.18g\n", tbcopy.bse_renv[0], tbcopy.bse_renv[1]);
		      eprintf("bse_tms=%.18g %.18g\n", tbcopy.bse_tms[0], tbcopy.bse_tms[1]);
		      eprintf("bse_bhspin=%.18g %.18g\n", tbcopy.bse_bhspin[0], tbcopy.bse_bhspin[1]);
		      exit_cleanly(1, __FUNCTION__);
		    }
		  }


		  /* Extract info from bcm array. PK.
		     In particular any interesting system that is formed from merger. 
		     At present this includes: NS.                                    */
		  //     ***********Handle any possibility of kicks***************    //
		  if(icase <= 100 && (tempbinary.bse_kw[0]==13 || tempbinary.bse_kw[1]==13)){
		    if(tbcopy.bse_kw[0]==13 || tbcopy.bse_kw[1]==13){
//			    tempbinary.bse_epoch[0] = tempbinary.bse_tphys;
//                            bse_set_merger(BSE_SIGMA);
//                            fprintf(stderr, "merge=265.0");
//                            if(tbcopy.bse_kw[0]==13){//if previously it was a pulsar...
//			      bse_set_merger(-1.0);
//                              fprintf(stderr, "merge=-1.0");
//                              if(((2.0*3.14159*31557600.0)/tempbinary.bse_ospin[0])<=0.030){
//                                bse_set_merger(-2.0);
//                                fprintf(stderr, "merge=-2.0");
//                              }
		      bse_set_merger(-1.0); //don't give a kick cos its already had one...
		    } else if (tbcopy.bse_kw[0]>= 10 && tbcopy.bse_kw[1]>=10 && tbcopy.bse_kw[1]<=12){
		      bse_set_merger(-20.0); //Merger induced collapse (like accretion induced collapse) probably results from electron capture onto Mg etc. 
		      fprintf(stderr, "merge=-20.0");
		      if(tempbinary.bse_kw[0]==13){
			tempbinary.bse_epoch[0] = tempbinary.bse_tphys;
		      } else {
			tempbinary.bse_epoch[1] = tempbinary.bse_tphys;
		      }
		    } else {
		      bse_set_merger(BSE_SIGMA);
		      if(tempbinary.bse_kw[0]==13){
			tempbinary.bse_epoch[0] = tempbinary.bse_tphys;
		      } else {
			tempbinary.bse_epoch[1] = tempbinary.bse_tphys;
		      }
		    }
		    tempbinary.bse_tb = 0.0;
		    dtp = 0.0;
// if msp set merger value to something other than -1.0 then it will trigger possibly a setting of NS particulars, a reset of NS particulars to msp-like values (-2.0), else reset to 'standard' birth properties or magnatar(see Maxim V. Barkov &Serguei S. Komissarov 2011?). Plus, of course we give the NSs kicks etc.
		    fprintf(stderr, "a ns is born? cmc_sscollision.c\n"); 
		    fprintf(stderr, "B4: tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g\n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1], tempbinary.bse_tb);
		    bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
				      &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
				      &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
				      &(tempbinary.bse_ospin[0]), &(tempbinary.bse_B_0[0]), &(tempbinary.bse_bacc[0]), &(tempbinary.bse_tacc[0]),
				      &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
				      &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
				      &(tempbinary.bse_tb), &(tempbinary.e), vs, &(tempbinary.bse_bhspin[0]));
		    fprintf(stderr, "a ns is born? cmc_sscollision.c\n"); 
		    fprintf(stderr, "RftR: tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g\n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1], tempbinary.bse_tb);
		    fprintf(stderr, "vk_y=%g should be>0...\n", vs[2]);
		  } else if (icase <= 100 && (tempbinary.bse_kw[0]==14 || tempbinary.bse_kw[1]==14)){
		    bse_set_merger(BSE_SIGMA);
		    fprintf(stderr, "BH merge=265.0\n");
		    if(tbcopy.bse_kw[0]==14 || tbcopy.bse_kw[1]==14){
		      bse_set_merger(-1.0);
		      fprintf(stderr, "BH merge=-1.0\n");
		    }
		    tempbinary.bse_tb = 0.0;
		    dtp = 0.0;
		    fprintf(stderr, "B4: tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g bhspin1=%g bhspin2=%g\n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1], tempbinary.bse_tb, tempbinary.bse_bhspin[0],tempbinary.bse_bhspin[1]);
		    bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
				      &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
				      &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
				      &(tempbinary.bse_ospin[0]), &(tempbinary.bse_B_0[0]), &(tempbinary.bse_bacc[0]), &(tempbinary.bse_tacc[0]), 
				      &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
				      &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
				      &(tempbinary.bse_tb), &(tempbinary.e), vs, &(tempbinary.bse_bhspin[0]));
		    fprintf(stderr, "RftR: tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g bhspin1=%g bhspin2=%g\n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1], tempbinary.bse_tb,tempbinary.bse_bhspin[0],tempbinary.bse_bhspin[1]);
		    fprintf(stderr, "BH vk_y=%g should be>0...\n", vs[2]);
		  }
		  bse_set_merger(-1.0);
                  j = 1;
                  while (bse_get_bcm(j,1) >= 0.0 && j<50000) {
                    j++;
                  }
                        //if(j>1){
                        //  tempbinary.bse_bcm_B[tbi] == bse_get_bcm(j-1,33+tbi);
                        //}
                  merged_star->se_scm_B = tempbinary.bse_bcm_B[tbi];
// check if the merger of a NS occured and if so make sure we copy the formation value to correct star (NS might have been star1 but merge put it in star2).
                  convert = 0;
                  if((tempbinary.bse_mass[0] == 0.0 || tempbinary.bse_mass[1] == 0.0) && (star1->se_scm_formation!=0 || star2->se_scm_formation!=0)){
                    if(tempbinary.bse_mass[0]!=0.0 && star2->se_scm_formation!=0 && star1->se_scm_formation==0){
                       convert = 1;
                    } else if (tempbinary.bse_mass[1]!=0.0 && star1->se_scm_formation!=0 && star2->se_scm_formation==0){
                       convert = 2;
                    }
                  }

// Check convert to see if formation was updated incorrectly. If so then correct it. PDK added...
                  if(convert>0){
                     if(convert==1){
                        merged_star->se_scm_formation = star2->se_scm_formation;
                        dprintf("sscoll, pulsar bin. convert1: %g %g %g %ld %ld \n",tempbinary.bse_B_0[tbi],tempbinary.bse_ospin[tbi],merged_star->se_scm_formation,tempbinary.id1,tempbinary.id2);

                     } else if (convert==2){
                        merged_star->se_scm_formation = star1->se_scm_formation;
                        dprintf("sscoll, pulsar bin. convert2: %g %g %g %ld %ld \n",tempbinary.bse_B_0[tbi],tempbinary.bse_ospin[tbi],merged_star->se_scm_formation,tempbinary.id1,tempbinary.id2);

                     }
                  }
/*                if(tbi == 0){
                   merged_star->se_scm_B = bse_get_bcm(j,33);
                } else {
                   merged_star->se_scm_B = bse_get_bcm(j,34);
                   } */
                  if(j > 1){
                     merged_star->se_scm_B = bse_get_bcm(j,33+tbi);
                     if(tempbinary.bse_kw[tbi] == 13) {
                        dprintf("sscoll, pulsar bin. merge5: %g %g %g \n",tempbinary.bse_B_0[tbi],tempbinary.bse_ospin[tbi],merged_star->se_scm_B);
                     }
                  }
		  
		  merged_star->se_mass = tempbinary.bse_mass0[tbi];
		  merged_star->se_k = tempbinary.bse_kw[tbi];
		  merged_star->se_mt = tempbinary.bse_mass[tbi];
		  merged_star->se_ospin = tempbinary.bse_ospin[tbi];
                  merged_star->se_B_0 = tempbinary.bse_B_0[tbi]; /* PK add in NS formation search (and updates)*/
                                                         /* i.e. add extraction of bcm array. */
                  merged_star->se_bacc = tempbinary.bse_bacc[tbi];
                  merged_star->se_tacc = tempbinary.bse_tacc[tbi];
		  merged_star->se_epoch = tempbinary.bse_epoch[tbi];
		  merged_star->se_tphys = tempbinary.bse_tphys;
		  merged_star->se_radius = tempbinary.bse_radius[tbi];
		  merged_star->se_lum = tempbinary.bse_lum[tbi];
		  merged_star->se_mc = tempbinary.bse_massc[tbi];
		  merged_star->se_rc = tempbinary.bse_radc[tbi];
		  merged_star->se_menv = tempbinary.bse_menv[tbi];
		  merged_star->se_renv = tempbinary.bse_renv[tbi];
		  merged_star->se_tms = tempbinary.bse_tms[tbi];
		  merged_star->se_bhspin = tempbinary.bse_bhspin[tbi];
          merged_star->Eint = tempbinary.Eint1 + tempbinary.Eint2;
		  
		  merged_star->rad = merged_star->se_radius * RSUN / units.l;
		  merged_star->m = merged_star->se_mt * MSUN / units.mstar;
		  
		  /* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
		  dtp = tphysf - merged_star->se_tphys;
		  dtp = 0.0;
		  /* Update star id for pass through. */
		  if(tbi==1){
		    bse_set_id1_pass(tempbinary.id1);
		    bse_set_id2_pass(0);
		  } else {
		    bse_set_id1_pass(tempbinary.id2);
		    bse_set_id2_pass(0);
		  }
		  temp_rad = tempbinary.bse_radius[0];
		  if(isnan(temp_rad)){
		    fprintf(stderr, "An isnan occured for r1 in sscollision\n");
		    fprintf(stderr, "tphys=%.18g tphysf=%.18g kstar1=%d kstar2=%d m1=%.18g m2=%.18g r1=%.18g r2=%.18g l1=%.18g l2=%.18g \n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1]);
		    exit(1);
		  }
		  temp_rad = tempbinary.bse_radius[1];
		  if(isnan(temp_rad)){
		    fprintf(stderr, "An isnan occured for r2 in sscollision\n");
		    fprintf(stderr, "tphys=%.18g tphysf=%.18g kstar1=%d kstar2=%d m1=%.18g m2=%.18g r1=%.18g r2=%.18g l1=%.18g l2=%.18g \n", tempbinary.bse_tphys, tphysf, tempbinary.bse_kw[0], tempbinary.bse_kw[1], tempbinary.bse_mass[0], tempbinary.bse_mass[1], tempbinary.bse_radius[0], tempbinary.bse_radius[1], tempbinary.bse_lum[0], tempbinary.bse_lum[1]);
		    exit(1);
		  }
		  /* REMOVED EVOLV1 CALL
		  bse_evolv1_safely(&(merged_star->se_k), &(merged_star->se_mass), &(merged_star->se_mt), &(merged_star->se_radius), 
				    &(merged_star->se_lum), &(merged_star->se_mc), &(merged_star->se_rc), &(merged_star->se_menv), 
				    &(merged_star->se_renv), &(merged_star->se_ospin), &(merged_star->se_epoch), &(merged_star->se_tms), 
				    &(merged_star->se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vsaddl);
		  */
		  
		  /*		vs[0] += vsaddl[0];
		    vs[1] += vsaddl[1];
		    vs[2] += vsaddl[2]; */
		  /* PDK vs array has now been increased. Although only 1st 3 would be used here anyways. */
		  /* DONT GIVE MULTIPLE KICKS WHEN NOT USING EVOLV1
		  for (i=0; i<=11; i++) {
                    vs[i] += vsaddl[i];
		  }
		  */
		  merged_star->rad = merged_star->se_radius * RSUN / units.l;
		  merged_star->m = merged_star->se_mt * MSUN / units.mstar;
		} else {
		  // remove the damn thing, how is this getting here anyway?
		  dprintf("In merge_two_stars: attempted to merge something not a binary! So skipping merging\n");
		  dprintf("Some info on skipped system: id1=%ld id2=%ld mass1=%g mass2=%g kw1=%d kw2=%d", tempbinary.id1,tempbinary.id2,tempbinary.bse_mass[0],tempbinary.bse_mass[1],tempbinary.bse_kw[0],tempbinary.bse_kw[1]);
		  //call destroy_obj(i) here (where i = k, star[k]... )? Or call zero out? Probably should update the system somehow, then we can check if it should be removed where the call to merge_two_stars was made from...
		}
	} else if (STAR_AGING_SCHEME && !STELLAR_EVOLUTION){
	  /* Sourav: toy rejuvenation version of stellar mergers */
	        merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);
	  /*		vs[0] = 0.0;
	    vs[1] = 0.0;
	    vs[2] = 0.0; */
	  /* */
	        for (i=0; i<=11; i++) {
	                vs[i] = 0.0;
                }		
		if (STAR_AGING_SCHEME==1){
			merged_star->lifetime = pow(10.0,9.921)*pow(merged_star->m*units.mstar/MSUN,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			age = (merged_star->lifetime/merged_star->m)*
				((star1->m*(TotalTime-star1->createtime)/star1->lifetime)
					+ star2->m*(TotalTime-star2->createtime)/star2->lifetime);
			merged_star->createtime = TotalTime-age;
		}
		else if (STAR_AGING_SCHEME==2){
		/*Sourav: zero lifetime of the collision product all else has infinite lifetime*/
			merged_star->lifetime = 0.0;
			merged_star->createtime = TotalTime;
		}
		else if (STAR_AGING_SCHEME==3) {
		/*Sourav: this is arbitrary lifetime for the collision products.  They are 
		uniformly chosen between 10^6 years to 10^8 years.  Reference: Sills et.al. 2007, 
		Sills et.al. 1997; all other stars have infinite lifetime*/
			//merged_star->lifetime = 1.0e6*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR + rng_t113_dbl() * (1.0e8-1.0e6)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			merged_star->lifetime = 1.0e6*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR + rng_t113_dbl_new(curr_st) * (1.0e8-1.0e6)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			merged_star->createtime = TotalTime;
		}
	} else {
		merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);

/*		vs[0] = 0.0;
		vs[1] = 0.0;
		vs[2] = 0.0; */
                for (i=0; i<=11; i++) {
                   vs[i] = 0.0;
                }
	}
  *curr_st= bse_get_taus113state();
}

/**
* @brief calculate resulting semi-major axis from collisional common envelope event
*
* @param Mrg mass of red giant
* @param Mint mass of intruder (NS, BH, MS, etc.)
* @param Mwd mass of RG core that will become WD
* @param Rrg radius of RG
* @param vinf relative velocity at infinity between RG and intruder
*
* @return resulting semi-major axis from collisional common envelope event
*/
double coll_CE(double Mrg, double Mint, double Mwd, double Rrg, double vinf)
{
	double alpha, lambda;
	
	alpha = bse_get_alpha1();
	lambda = bse_get_lambdaf();

	return(1.0/(2.0*Mrg*(Mrg-Mwd)/(Mwd*Mint*alpha*lambda*Rrg)-(Mrg+Mint)/(Mwd*Mint*alpha)*vinf*vinf));

}

/**
* @brief calculate resulting semi-major axis from collisional common envelope event for two giant stars
*
* @param M1 mass of giant
* @param M2 mass of giant
* @param Mc1 mass of giant core that will become WD
* @param Mc2 mass of giant core that will become WD
* @param R1 radius of giant
* @param R2 radius of giant
* @param vinf relative velocity at infinity between the two giant stars
*
* @return resulting semi-major axis from collisional common envelope event
*/
double coll_CE_twogiant(double M1, double M2, double Mc1, double Mc2, double R1, double R2, double vinf)
{
        double alpha, lambda;

        alpha = bse_get_alpha1();
        lambda = bse_get_lambdaf();

        return(Mc1*Mc2*alpha/(2*(M1*(M1-Mc1)/R1+M2*(M2-Mc2)/R2)/lambda-(M1+M2)*vinf*vinf));
}

/**
* @brief
* the mass--radius relationship, from Freitag, et al. (2004)
* M is expected to be input in the same units as star[].m
* radius is output in code (N-body) units
*
* This is ONLY used if stellar evolution is off (for a toy rejuvination problem)
*
* @param M ?
*
* @return radius in code (N-body) units
*/
double r_of_m(double M)
{
	if (M == 0.0) {
		/* special case of massless particles */
		return(0.0);
	} else if (M / clus.N_STAR * units.m / MSUN < 0.1) {
		return(0.1 * RSUN / units.l);
	} else if (M / clus.N_STAR * units.m / MSUN < 1.0) {
		return(RSUN / units.l * (M / clus.N_STAR * units.m / MSUN));
	} else if (M / clus.N_STAR * units.m / MSUN < 120.0) {
		return(RSUN / units.l * pow(M / clus.N_STAR * units.m / MSUN, 0.57));
	} else {
		return(1.6 * RSUN / units.l * pow(M / clus.N_STAR * units.m / MSUN, 0.47));
	}
}
