/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap/bse_wrap.h"

void sscollision_do(long k, long kp, double rperimax, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng)
{
	long knew;
	double vs[3], bmax, b, rperi, Eorbnew, acoll, ecoll, ace, ece, anew, enew, efinal, afinal;
	double aj, tm, tn, tscls[20], lums[10], GB[10], k2;
	double Einit;

	bmax = rperimax * sqrt(1.0 + 2.0 * ((star[k].m + star[kp].m) * madhoc) / (rperimax * sqr(W)));
	b = sqrt(rng_t113_dbl()) * bmax;
	rperi = madhoc*(star[k].m+star[kp].m)/sqr(W) * (-1.0+sqrt(1.0+sqr(b*W*W/(madhoc*star[k].m+star[kp].m))));

	if (TIDAL_CAPTURE && (star[k].se_k <= 1 || star[k].se_k >= 10) && (star[kp].se_k >= 2 && star[kp].se_k <= 9 && star[kp].se_k != 7) && 
	    rperi <= 1.3 * star[kp].rad) {
		/* log stuff */
		fprintf(tidalcapturefile, "%.3g SS_TC %s+%s->", TotalTime, 
			sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring));

		/* instead of a merger, form a CV, WD-WD binary, or UCXB from the Ivanova & Lombardi collision mechanism */
		ecoll = 0.88 - rperi/(3.0*star[kp].rad);
		acoll = rperi/(3.3*(1.0-sqr(ecoll)));
		
		ace = coll_CE(star[kp].m*madhoc, star[k].m*madhoc, star[kp].se_mc*MSUN/units.mstar, star[kp].se_radius*RSUN/units.l, W);
		ece = 0.0;

		if (ace < acoll) {
			afinal = ace;
			efinal = ece;
		} else {
			afinal = acoll;
			efinal = ecoll;
		}
		
		/* strip envelope of RG */
		Einit = star[k].Eint + 0.5 * star[k].m * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
			star[kp].Eint + 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) + 
			0.5 * star[k].m * madhoc * star[k].phi + 0.5 * star[kp].m * madhoc * star[kp].phi;
		DMse += star[kp].m * madhoc;
		aj = star[kp].se_tphys - star[kp].se_epoch;
		star[kp].se_mt = star[kp].se_mc;
		bse_star(&(star[kp].se_k), &(star[kp].se_mass), &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars);
		bse_hrdiag(&(star[kp].se_mass), &aj, &(star[kp].se_mt), &tm, &tn, tscls, lums, GB, zpars,
			   &(star[kp].se_radius), &(star[kp].se_lum), &(star[kp].se_k), &(star[kp].se_mc), &(star[kp].se_rc), 
			   &(star[kp].se_menv), &(star[kp].se_renv), &k2);
		star[kp].se_epoch = star[kp].se_tphys - aj;
		star[kp].rad = star[kp].se_radius * RSUN / units.l;
		star[kp].m = star[kp].se_mt * MSUN / units.mstar;
		DMse -= star[kp].m * madhoc;

		/* check to see if MS star overfills RL at pericenter, and destroy if this is the case */
		if (star[k].se_k <= 1 && star[k].se_radius*RSUN/units.l >= bse_rl(star[k].m/star[kp].m)*(1.0-efinal)*afinal) {
			/* keep only RG core as single star */
			star[kp].r = rcm;
			star[kp].vr = vcm[3];
			star[kp].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			star[kp].phi = potential(star[kp].r);
			set_star_EJ(kp);
			set_star_news(kp);
			set_star_olds(kp);
			
			/* mark stars as interacted so they don't undergo E_CONS mode stuff */
			star[kp].interacted = 1;
			star[kp].Eint = star[k].Eint + star[kp].Eint + Einit 
				- 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt)) 
				- 0.5 * star[kp].m * madhoc * star[kp].phi;
			
			/* log stuff */
			fprintf(tidalcapturefile, "%s\n", sprint_star_dyn(kp, dummystring));

			destroy_obj(k);
		} else {
			/* form compact binary with stripped RG core*/
			/* put new binary together and destroy original stars */
			knew = create_binary();
			star[knew].m = star[k].m + star[kp].m;
			star[knew].r = rcm;
			star[knew].vr = vcm[3];
			star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			star[knew].phi = potential(star[knew].r);
			set_star_EJ(knew);
			set_star_news(knew);
			set_star_olds(knew);
			star[knew].interacted = 1;
			
			binary[star[knew].binind].a = afinal;
			binary[star[knew].binind].e = efinal;
			binary[star[knew].binind].m1 = star[k].m;
			binary[star[knew].binind].m2 = star[kp].m;
			binary[star[knew].binind].rad1 = star[k].rad;
			binary[star[knew].binind].rad2 = star[kp].rad;
			binary[star[knew].binind].Eint1 = star[k].Eint;
			binary[star[knew].binind].Eint2 = star[kp].Eint;

			/* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
			binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);
			binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);

			binary[star[knew].binind].id1 = star[k].id;
			binary[star[knew].binind].id2 = star[kp].id;
			cp_SEvars_to_newbinary(k, -1, knew, 0);
			cp_SEvars_to_newbinary(kp, -1, knew, 1);
			binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			compress_binary(&star[knew], &binary[star[knew].binind]);
			
			/* log stuff */
			fprintf(tidalcapturefile, "%s\n", sprint_bin_dyn(knew, dummystring));

			destroy_obj(k);
			destroy_obj(kp);
		}
	} else if (TIDAL_CAPTURE && (star[kp].se_k <= 1 || star[kp].se_k >= 10) && (star[k].se_k >= 2 && star[k].se_k <= 9 && star[k].se_k != 7) && 
		   rperi <= 1.3 * star[k].rad) {
		/* log stuff */
		fprintf(tidalcapturefile, "%.3g SS_TC %s+%s->", TotalTime, 
			sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring));

		/* instead of a merger, form a CV, WD-WD binary, or UCXB from the Ivanova & Lombardi collision mechanism */
		ecoll = 0.88 - rperi/(3.0*star[k].rad);
		acoll = rperi/(3.3*(1.0-sqr(ecoll)));
		
		ace = coll_CE(star[k].m*madhoc, star[kp].m*madhoc, star[k].se_mc*MSUN/units.mstar, star[k].se_radius*RSUN/units.l, W);
		ece = 0.0;

		if (ace < acoll) {
			afinal = ace;
			efinal = ece;
		} else {
			afinal = acoll;
			efinal = ecoll;
		}
		
		/* strip envelope of RG */
		Einit = star[k].Eint + 0.5 * star[k].m * madhoc * (sqr(star[k].vr)+sqr(star[k].vt)) +
			star[kp].Eint + 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr)+sqr(star[kp].vt)) + 
			0.5 * star[k].m * madhoc * star[k].phi + 0.5 * star[kp].m * madhoc * star[kp].phi;
		DMse += star[k].m * madhoc;
		aj = star[k].se_tphys - star[k].se_epoch;
		star[k].se_mt = star[k].se_mc;
		bse_star(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars);
		bse_hrdiag(&(star[k].se_mass), &aj, &(star[k].se_mt), &tm, &tn, tscls, lums, GB, zpars,
			   &(star[k].se_radius), &(star[k].se_lum), &(star[k].se_k), &(star[k].se_mc), &(star[k].se_rc), 
			   &(star[k].se_menv), &(star[k].se_renv), &k2);
		star[k].se_epoch = star[k].se_tphys - aj;
		star[k].rad = star[k].se_radius * RSUN / units.l;
		star[k].m = star[k].se_mt * MSUN / units.mstar;
		DMse -= star[k].m * madhoc;

		/* check to see if MS star overfills RL at pericenter, and destroy if this is the case */
		if (star[kp].se_k <= 1 && star[kp].se_radius*RSUN/units.l >= bse_rl(star[kp].m/star[k].m)*(1.0-efinal)*afinal) {
			/* keep only RG core as single star */
			star[k].r = rcm;
			star[k].vr = vcm[3];
			star[k].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			star[k].phi = potential(star[k].r);
			set_star_EJ(k);
			set_star_news(k);
			set_star_olds(k);
			
			/* mark stars as interacted so they don't undergo E_CONS mode stuff */
			star[k].interacted = 1;
			star[k].Eint = star[k].Eint + star[kp].Eint + Einit 
				- 0.5 * star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				- 0.5 * star[k].m * madhoc * star[k].phi;
			
			/* log stuff */
			fprintf(tidalcapturefile, "%s\n", sprint_star_dyn(k, dummystring));

			destroy_obj(kp);
		} else {
			/* form compact binary with stripped RG core*/
			/* put new binary together and destroy original stars */
			knew = create_binary();
			star[knew].m = star[k].m + star[kp].m;
			star[knew].r = rcm;
			star[knew].vr = vcm[3];
			star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			star[knew].phi = potential(star[knew].r);
			set_star_EJ(knew);
			set_star_news(knew);
			set_star_olds(knew);
			star[knew].interacted = 1;
			
			binary[star[knew].binind].a = afinal;
			binary[star[knew].binind].e = efinal;
			binary[star[knew].binind].m1 = star[k].m;
			binary[star[knew].binind].m2 = star[kp].m;
			binary[star[knew].binind].rad1 = star[k].rad;
			binary[star[knew].binind].rad2 = star[kp].rad;
			binary[star[knew].binind].Eint1 = star[k].Eint;
			binary[star[knew].binind].Eint2 = star[kp].Eint;

			/* put lost energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
			binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * (Einit
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);
			binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * (Einit
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);

			binary[star[knew].binind].id1 = star[k].id;
			binary[star[knew].binind].id2 = star[kp].id;
			cp_SEvars_to_newbinary(k, -1, knew, 0);
			cp_SEvars_to_newbinary(kp, -1, knew, 1);
			binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			compress_binary(&star[knew], &binary[star[knew].binind]);
			
			/* log stuff */
			fprintf(tidalcapturefile, "%s\n", sprint_bin_dyn(knew, dummystring));

			destroy_obj(k);
			destroy_obj(kp);
		}
	} else if (rperi <= star[k].rad + star[kp].rad) {
		/* perform standard sticky-sphere merger */
		/* If tidal capture is turned off, the cross section is just large enough to enter this clause, 
		   so the next clause should never be entered. */

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
		fprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typem=%d type1=%d type2=%d\n", 
			TotalTime, 
			star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN, 
			star[k].id, star[k].m * units.mstar / FB_CONST_MSUN, 
			star[kp].id, star[kp].m * units.mstar / FB_CONST_MSUN,
			star[knew].r, star[knew].se_k, star[k].se_k, star[kp].se_k);

		/* destroy two progenitors */
		destroy_obj(k);
		destroy_obj(kp);
	} else if (TIDAL_CAPTURE)  {
		/* apply tidal capture / common envelope test */
		Eorbnew = 0.5*madhoc*star[k].m*star[kp].m/(star[k].m+star[kp].m)*sqr(W);

		if (star[k].se_k == 1) {
			Eorbnew -= Etide(rperi, star[k].m*madhoc, star[k].rad, 3.0, star[kp].m*madhoc);
		} else if (star[k].se_k < 10) {
			Eorbnew -= Etide(rperi, star[k].m*madhoc, star[k].rad, 1.5, star[kp].m*madhoc);
		}

		if (star[kp].se_k == 1) {
			Eorbnew -= Etide(rperi, star[kp].m*madhoc, star[kp].rad, 3.0, star[k].m*madhoc);
		} else if (star[kp].se_k < 10) {
			Eorbnew -= Etide(rperi, star[kp].m*madhoc, star[kp].rad, 1.5, star[k].m*madhoc);
		}

		if (Eorbnew < 0.0) {
			/* bound system; don't worry about RL overflow here, since BSE will take care of that automatically */
			anew = madhoc * star[k].m * madhoc * star[kp].m / (2.0 * fabs(Eorbnew));
			enew = MIN(0.0, 1.0-rperi/anew);
			
			/* apply rapid tidal circularization of orbit, assuming no angular momentum is transferred to 
			   stars' internal rotation for simplicity; but only if there is no Roche-lobe overflow at 
			   pericenter */
			if (bse_rl(star[k].m/star[kp].m)*anew*(1.0-enew) > star[k].rad && 
			    bse_rl(star[kp].m/star[k].m)*anew*(1.0-enew) > star[kp].rad) {
				efinal = 0.0;
				afinal = anew * (1.0 + enew);
			} else {
				afinal = anew;
				efinal = enew;
			}

			/* put new binary together and destroy original stars */
			knew = create_binary();
			star[knew].m = star[k].m + star[kp].m;
			star[knew].r = rcm;
			star[knew].vr = vcm[3];
			star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
			star[knew].phi = potential(star[knew].r);
			set_star_EJ(knew);
			set_star_news(knew);
			set_star_olds(knew);
			star[knew].interacted = 1;
			
			binary[star[knew].binind].a = afinal;
			binary[star[knew].binind].e = efinal;
			binary[star[knew].binind].m1 = star[k].m;
			binary[star[knew].binind].m2 = star[kp].m;
			binary[star[knew].binind].rad1 = star[k].rad;
			binary[star[knew].binind].rad2 = star[kp].rad;
			binary[star[knew].binind].Eint1 = star[k].Eint;
			binary[star[knew].binind].Eint2 = star[kp].Eint;

			/* put tidal energy into Eint of each star, divided equally (it's just for bookkeeping anyway) */
			binary[star[knew].binind].Eint1 = star[k].Eint + 0.5 * 
				(0.5 * star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				 + 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 + 0.5 * star[k].m * madhoc * star[k].phi
				 + 0.5 * star[kp].m * madhoc * star[kp].phi
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);
			binary[star[knew].binind].Eint2 = star[kp].Eint + 0.5 * 
				(0.5 * star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
				 + 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
				 - 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
				 + 0.5 * star[k].m * madhoc * star[k].phi
				 + 0.5 * star[kp].m * madhoc * star[kp].phi
				 - 0.5 * star[knew].m * madhoc * star[knew].phi
				 + 0.5 * star[k].m * madhoc * star[kp].m * madhoc / binary[star[knew].binind].a);

			binary[star[knew].binind].id1 = star[k].id;
			binary[star[knew].binind].id2 = star[kp].id;
			cp_SEvars_to_newbinary(k, -1, knew, 0);
			cp_SEvars_to_newbinary(kp, -1, knew, 1);
			binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
			compress_binary(&star[knew], &binary[star[knew].binind]);

			/* log stuff */
			fprintf(tidalcapturefile, "%.3g SS_TC %s+%s->%s\n", TotalTime, 
				sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring), sprint_bin_dyn(knew, dummystring));
			
			destroy_obj(k);
			destroy_obj(kp);
		} else {
			fprintf(tidalcapturefile, "%.3g SS_TC_FAILED %s+%s->%s+%s\n", TotalTime, 
				sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring),
				sprint_star_dyn(k, dummystring), sprint_star_dyn(kp, dummystring));
			/* this shouldn't happen if the cross section is accurate, so just ignore */
		}
	}
}

/* merge two stars using stellar evolution if it's enabled */
void merge_two_stars(star_t *star1, star_t *star2, star_t *merged_star, double *vs) {
	double tphysf, dtp, vsaddl[3], age;
	binary_t tempbinary, tbcopy;
	int tbi=-1, j, ktry;
	
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
		tempbinary.m1 = star1->m;
		tempbinary.m2 = star2->m;
		tempbinary.Eint1 = star1->Eint;
		tempbinary.Eint2 = star2->Eint;
		tempbinary.a = BSE_WRAP_MAX(10.0*star1->rad/bse_rl(star1->m/star2->m), 
					    10.0*star2->rad/bse_rl(star2->m/star1->m));
		tempbinary.e = 1.0 - 0.1*(star1->rad+star2->rad)/tempbinary.a;
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
		
		/* dtp = tphysf - tempbinary.bse_tphys; */
		/* Since the evolution time is so short in this routine, we can simply set dtp=0.0
		   without worrying about the bcm arrays filling up. */
		dtp = 0.0;
		bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
			   &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
			   &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
			   &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
			   &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
			   &(tempbinary.bse_tb), &(tempbinary.e), vs);
		
		/* Force merger if necessary by resetting pericenter to nearly 0;  This may be needed in some cases
		   because BSE doesn't strictly conserve angular momentum in binaries during common envelope evolution. */
		ktry = 0;
		while (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0 && ktry < 10) {
			dprintf("Attempting to force merger in BSE by repeating evolution with tiny pericenter.\n");
			tempbinary.a = 1.0e-12 * AU / units.l;
			tempbinary.e = 0.999;
			tempbinary.bse_tb = sqrt(cub(tempbinary.a * units.l / AU)/(tempbinary.bse_mass[0]+tempbinary.bse_mass[1]))*365.25;
			dtp = 0.0;
			tphysf += 1.0e-6;
			bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
					  &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
					  &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
					  &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
					  &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
					  &(tempbinary.bse_tb), &(tempbinary.e), vs);
			ktry++;
		}

		/*Sourav:debug if the collision takes more time to happen?*/
		/* if (tempbinary.bse_mass[0] != 0.0 && tempbinary.bse_mass[1] != 0.0) { */
/* 		  tphysf+=1.; */
/* 		  dtp= 0.; */
/* 		  bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]),  */
/* 			   &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]),  */
/* 			   &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]),  */
/* 			   &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]),  */
/* 			   &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars,  */
/* 			   &(tempbinary.bse_tb), &(tempbinary.e), vs); */
/* 		  fprintf (stderr, "\n*******bugfix2 for sscollision********\n"); */
/* 		  j=1; */
/* 		  while (bse_get_bpp(j,1)>=0.0) { */
/* 			  if (bse_get_bpp(j,4)==15 || bse_get_bpp(j,5)==15){ */
/* 				tphysf=bse_get_bpp(j,1)+1.0e-06; */
/* 				fprintf (stderr, "k1=bse_get_bpp(j,4)= %d,k2=bse_get_bpp(j,5)= %d\n",bse_get_bpp(j,4),bse_get_bpp(j,5)); */
/* 				break; */
/* 			  } */
/* 			  if(j>80){ */
/* 				  eprintf ("no 15 found"); */
/* 				  exit_cleanly(1); */
/* 			  } */
/* 			  j++; */
/* 		  } */
/* 		  bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]),  */
/* 			   &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]),  */
/* 			   &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]),  */
/* 			   &(tempbinary.bse_ospin[0]), &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]),  */
/* 			   &(tempbinary.bse_tphys), &tphysf, &dtp, &METALLICITY, zpars,  */
/* 			   &(tempbinary.bse_tb), &(tempbinary.e), vs); */
/* 		} */
	
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
	} else if (STAR_AGING_SCHEME && !STELLAR_EVOLUTION){
		/* Sourav: toy rejuvenation version of stellar mergers */
		merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);
		vs[0] = 0.0;
		vs[1] = 0.0;
		vs[2] = 0.0;
		
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
			merged_star->lifetime = 1.0e6*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR + rng_t113_dbl() * (1.0e8-1.0e6)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			merged_star->createtime = TotalTime;
		}
	} else {
		merged_star->m = star1->m + star2->m;
		merged_star->rad = r_of_m(merged_star->m);

		vs[0] = 0.0;
		vs[1] = 0.0;
		vs[2] = 0.0;
	}
}

/* calculate resulting semi-major axis from collisional common envelope event */
double coll_CE(double Mrg, double Mint, double Mwd, double Rrg, double vinf)
/* 
   Mrg = mass of red giant
   Mint = mass of intruder (NS, BH, MS, etc.)
   Mwd = mass of RG core that will become WD
   Rrg = radius of RG
   vinf = relative velocity at infinity between RG and intruder
*/
{
	double alpha, lambda;
	
	alpha = bse_get_alpha1();
	lambda = bse_get_lambda();

	return(1.0/(2.0*Mrg*(Mrg-Mwd)/(Mwd*Mint*alpha*lambda*Rrg)-(Mrg+Mint)/(Mwd*Mint)*vinf*vinf));

}
