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
	double vs[3];
	
	/* create new star */
	knew = create_star();

	/* merge parent stars, setting mass, stellar radius, and SE params */
	merge_two_stars(&(star[k]), &(star[kp]), &(star[knew]), vs);
	
	star[knew].r = rcm;
	star[knew].vr = vcm[3];
	star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
	star[knew].vr += vs[2] * 1.0e5 / (units.l/units.t);
	star[knew].vt += sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
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
void merge_two_stars(star_t *star1, star_t *star2, star_t *merged_star, double *vs) {
	double tphysf, dtp, vsaddl[3];
	binary_t tempbinary, tbcopy;
	int tbi=-1, j;
	
	if (STELLAR_EVOLUTION) {
		/* evolve for just a year for merger */
		tphysf = star1->se_tphys + 1.0e-6;
		
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
		
		tempbinary.bse_tb = sqrt(cub(tempbinary.a * units.l / AU)/(tempbinary.bse_mass[0]+tempbinary.bse_mass[1]))*365.25;
		
		tbcopy = tempbinary;
		
		dtp = tphysf - tempbinary.bse_tphys;
		bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
			   &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
			   &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
			   &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
			   &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
			   &(tempbinary.bse_tb), &(tempbinary.e), vs);
		
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
			j = 1;
			while (bse_get_bpp(j, 1) >= 0.0) {
				fprintf(stderr, "time=%g m1=%g m2=%g k1=%d k2=%d sep=%g ecc=%g r1/rol1=%g r2/rol2=%g type=%s\n",
					bse_get_bpp(j, 1), bse_get_bpp(j, 2), bse_get_bpp(j, 3), (int) bse_get_bpp(j, 4), 
					(int) bse_get_bpp(j, 5), bse_get_bpp(j, 6), bse_get_bpp(j, 7), bse_get_bpp(j, 8), 
					bse_get_bpp(j, 9), bse_get_bselabel((int) bse_get_bpp(j, 10)));
				fflush(NULL);
				j++;
			}
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
				exit_cleanly(1);
			}
		}
	
		merged_star->se_mass = tempbinary.bse_mass0[tbi];
		merged_star->se_k = tempbinary.bse_kw[tbi];
		merged_star->se_mt = tempbinary.bse_mass[tbi];
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

		merged_star->rad = merged_star->se_radius * RSUN / units.l;
		merged_star->m = merged_star->se_mt * MSUN / units.mstar;
		
		/* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
		dtp = tphysf - merged_star->se_tphys;
		bse_evolv1_safely(&(merged_star->se_k), &(merged_star->se_mass), &(merged_star->se_mt), &(merged_star->se_radius), 
				  &(merged_star->se_lum), &(merged_star->se_mc), &(merged_star->se_rc), &(merged_star->se_menv), 
				  &(merged_star->se_renv), &(merged_star->se_ospin), &(merged_star->se_epoch), &(merged_star->se_tms), 
				  &(merged_star->se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vsaddl);
		
		vs[0] += vsaddl[0];
		vs[1] += vsaddl[1];
		vs[2] += vsaddl[2];

		merged_star->rad = merged_star->se_radius * RSUN / units.l;
		merged_star->m = merged_star->se_mt * MSUN / units.mstar;
	} else {
		merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);

		vs[0] = 0.0;
		vs[1] = 0.0;
		vs[2] = 0.0;
	}
}
