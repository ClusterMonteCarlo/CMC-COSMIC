/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap/bse_wrap.h"

void sscollision_do(long k, long kp, double rcm, double vcm[4])
{
	long knew;
	
	/* create new star */
	knew = create_star();

	/* merge parent stars, setting mass, stellar radius, and SE params */
	merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]));
	
	star[knew].r = rcm;
	star[knew].vr = vcm[3];
	star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
	star[knew].phi = potential(star[knew].r);
	set_star_EJ(knew);
	set_star_news(knew);
	set_star_olds(knew);
	
	/* mark stars as interacted so they don't undergo E_CONS mode stuff */
	star[knew].id = star_get_id_new();
	star[knew].interacted = 1;
	star[knew].Eint = star[k].Eint + star[kp].Eint 
		+ 0.5 * star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
		+ 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
		- 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
		+ 0.5 * star[k].m * madhoc * star[k].phi
		+ 0.5 * star[kp].m * madhoc * star[kp].phi
		- 0.5 * star[knew].m * madhoc * star[knew].phi;
	
	
	/* log collision */
	fprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g)\n", 
		TotalTime, 
		star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN, 
		star[k].id, star[k].m * units.mstar / FB_CONST_MSUN, 
		star[kp].id, star[kp].m * units.mstar / FB_CONST_MSUN,
		star[knew].r);
	
	/* destroy two progenitors */
	destroy_obj(k);
	destroy_obj(kp);
}

/* merge two stars using stellar evolution if it's enabled */
void merge_two_stars(star_t *star1, star_t *star2, star_t *merged_star) {
	double tphysf, dtp;
	binary_t tempbinary;
	int tbi=-1;
	
	if (STELLAR_EVOLUTION) {
		/* create temporary tight, eccentric binary that will merge immediately */
		tempbinary.id1 = star1->id;
		tempbinary.id2 = star2->id;
		tempbinary.rad1 = star1->rad;
		tempbinary.rad2 = star2->rad;
		tempbinary.m1 = star1->m;
		tempbinary.m2 = star2->m;
		tempbinary.Eint1 = star1->Eint;
		tempbinary.Eint2 = star2->Eint;
		tempbinary.a = 1.0e-12 * AU / units.l;
		tempbinary.e = 0.999;
		tempbinary.inuse = 1;
		tempbinary.bse_mass0[0] = star1->se_mass;
		tempbinary.bse_mass0[1] = star2->se_mass;
		tempbinary.bse_kw[0] = star1->se_k;
		tempbinary.bse_kw[1] = star2->se_k;
		tempbinary.bse_mass[0] = star1->se_mt;
		tempbinary.bse_mass[1] = star2->se_mt;
		tempbinary.bse_ospin[0] = star1->se_ospin;
		tempbinary.bse_ospin[1] = star2->se_ospin;
		tempbinary.bse_epoch[0] = star1->se_epoch;
		tempbinary.bse_epoch[1] = star2->se_epoch;
		tempbinary.bse_tphys = star1->se_tphys; /* tphys should be the same for both input stars so this should be OK */
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
		
		/* evolve for just a year for merger */
		tphysf = tempbinary.bse_tphys + 1.0e-6;
		dtp = tphysf;
		tempbinary.bse_tb = sqrt(cub(tempbinary.a * units.l / AU)/(tempbinary.bse_mass[0]+tempbinary.bse_mass[1]))*365.25;
		bse_evolv2(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), &(tempbinary.bse_radius[0]), 
			   &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), 
			   &(tempbinary.bse_renv[0]), &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
			   &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
			   &(tempbinary.bse_tb), &(tempbinary.e));
		
		/* make sure outcome was as expected */
		if (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0) {
			eprintf("Artificial stellar evolution of eccentric binary failed: both stars have non-zero mass.\n");
			exit_cleanly(1);
		} else if (tempbinary.bse_mass[0] != 0) {
			tbi = 0;
		} else if (tempbinary.bse_mass[1] != 0) {
			tbi = 1;
		} else {
			if (tempbinary.bse_kw[0] == 15) {
				tbi = 0;
			} else if (tempbinary.bse_kw[1] == 15) {
				tbi = 1;
			} else {
				eprintf("Artificial stellar evolution of eccentric binary failed: both stars have zero mass.\n");
				eprintf("k1=%d k2=%d kb1=%d kb2=%d\n", star1->se_k, star2->se_k, 
					tempbinary.bse_kw[0], tempbinary.bse_kw[1]);
				exit_cleanly(1);
			}
		}
		
		/* debug stuff */
		dprintf("\ninstar1: se_k=%d se_mass=%g se_epoch=%g se_mc=%g se_menv=%g se_tms=%g\n", 
			star1->se_k, star1->se_mass, star1->se_epoch, star1->se_mc, star1->se_menv, star1->se_tms);
		dprintf("instar2: se_k=%d se_mass=%g se_epoch=%g se_mc=%g se_menv=%g se_tms=%g\n", 
			star2->se_k, star2->se_mass, star2->se_epoch, star2->se_mc, star2->se_menv, star2->se_tms);
		dprintf("outcome: se_k=%d se_mass=%g se_epoch=%g se_mc=%g se_menv=%g se_tms=%g\n",
			tempbinary.bse_kw[tbi], tempbinary.bse_mass[tbi], tempbinary.bse_epoch[tbi], tempbinary.bse_massc[tbi], 
			tempbinary.bse_menv[tbi], tempbinary.bse_tms[tbi]);
		
		merged_star->m = tempbinary.bse_mass[tbi] * MSUN / units.mstar;
		merged_star->rad = tempbinary.bse_radius[tbi] * RSUN / units.l;
		merged_star->se_mass = tempbinary.bse_mass0[tbi]; /* initial mass (at curent epoch?) */
		merged_star->se_k = tempbinary.bse_kw[tbi];
		merged_star->se_mt = tempbinary.bse_mass[tbi]; /* current mass */
		merged_star->se_ospin = tempbinary.bse_ospin[tbi];
		merged_star->se_epoch = tempbinary.bse_epoch[tbi];
		merged_star->se_tphys = tempbinary.bse_tphys;
		merged_star->se_radius = tempbinary.bse_radius[tbi];
		merged_star->se_lum = tempbinary.bse_lum[tbi];
		merged_star->se_mc = tempbinary.bse_massc[tbi];
		merged_star->se_rc = tempbinary.bse_radc[tbi];
		merged_star->se_menv = tempbinary.bse_menv[tbi];
		merged_star->se_renv = tempbinary.bse_renv[tbi];
		merged_star->se_tms = tempbinary.bse_tms[tbi];
	} else {
		merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);
	}
}
