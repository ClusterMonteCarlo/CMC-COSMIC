#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap/bse_wrap.h"

void stellar_evolution_init(void){  
  double tphysf, dtp;
  long k;

  /* SSE */
  /* bse_set_hewind(0.5); */

  /* BSE */
  bse_set_neta(0.5);
  bse_set_bwind(0.0);
  bse_set_hewind(1.0);
  bse_set_alpha1(3.0);
  bse_set_lambda(0.5);
  bse_set_ceflag(0);
  bse_set_tflag(1);
  bse_set_ifflag(0);
  bse_set_wdflag(1);
  bse_set_bhflag(0);
  bse_set_nsflag(1);
  bse_set_mxns(3.0);
  bse_set_idum(29769);
  bse_set_pts1(0.05);
  bse_set_pts2(0.01);
  bse_set_pts3(0.02);
  bse_set_sigma(190.0);
  bse_set_beta(0.125);
  bse_set_xi(1.0);
  bse_set_acc2(1.5);
  bse_set_epsnov(0.001);
  bse_set_eddfac(10.0);
  bse_set_gamma(-1.0);

  /* set parameters relating to metallicity */
  zpars = (double *) malloc(20 * sizeof(double));
  bse_zcnsts(&METALLICITY, zpars);

  /* set collisions matrix */
  bse_instar();

  /* set initial properties of stars */
  for (k=1; k<=clus.N_MAX; k++){
    star[k].se_mass = star[k].m * units.mstar / MSUN;
    /* setting the type */
    if(star[k].se_mass <= 0.7){
      star[k].se_k = 0;
    } else {
      star[k].se_k = 1;
    }
    star[k].se_mt = star[k].se_mass;
    star[k].se_ospin = 0.0;
    star[k].se_epoch = 0.0;
    star[k].se_tphys = 0.0;
  }

  /* evolve slightly (1 year) for initial radii */
  tphysf = 1.0e-6;
  dtp = tphysf;
  for(k=1; k<=clus.N_MAX; k++){
    bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
	       &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
	       &(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
	       &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars);
    star[k].rad = star[k].se_radius * RSUN / units.l;
    star[k].m = star[k].se_mass * MSUN / units.mstar;
  }
}

/* note that this routine is called after perturb_stars() and get_positions() */
void do_stellar_evolution(gsl_rng *rng){
  long k;
  int kprev;
  double dtp, tphysf, theta, vk;
  
  for(k=1; k<=clus.N_MAX; k++){
    tphysf = TotalTime / MEGA_YEAR;
    dtp = tphysf;
    kprev = star[k].se_k;

    bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
	       &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
	       &(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
	       &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars);

    star[k].rad = star[k].se_radius * RSUN / units.l;
    star[k].m = star[k].se_mass * MSUN / units.mstar;

    /* impose compact object birth kick, and add speed to systemic speed */
    if ((star[k].se_k == 13 || star[k].se_k == 14) && star[k].se_k != kprev) {
      /* convert speed in km/s to code units */
      vk = bse_kick_speed(&(star[k].se_k)) * 1.0e5 / (units.l/units.t);
      theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0);
      star[k].vr += cos(theta) * vk;
      star[k].vt += sin(theta) * vk;
      set_star_EJ(k);
    }
  }
}

void write_stellar_data(void){
	long k;
	FILE *stel_file;
	char filename[1024];

	se_file_counter++;
	sprintf(filename, "%s_stellar_info.%05d", outprefix, se_file_counter);
	stel_file = fopen(filename, "w");
	if (stel_file==NULL){
		fprintf(stderr,
			"file cannot be opened to write stellar info\n");
		return;
	}
	fprintf(stel_file, "# time (Myr): %e\n",
			TotalTime/MEGA_YEAR);
	fprintf(stel_file, "# time (FP):  %e\n", TotalTime);
	fprintf(stel_file,
	       "#  id        mass        radius     luminosity  type\n");
	fprintf(stel_file,
	       "#======= ============ ============ ============ ====\n");
	for(k=1; k<=clus.N_MAX; k++){
		fprintf(stel_file, "%08ld ", k);
		fprintf(stel_file, "%e ", star[k].se_mass);
		fprintf(stel_file, "%e ", star[k].se_radius);
		fprintf(stel_file, "%e ", star[k].se_lum);
		fprintf(stel_file, "%d ", star[k].se_k);
		fprintf(stel_file, "\n");
	}
	fclose(stel_file);
}
