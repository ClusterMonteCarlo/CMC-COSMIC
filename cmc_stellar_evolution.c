#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap/bse_wrap.h"

void stellar_evolution_init(void){  
  double tphysf, dtp;
  int i;
  long k, kb;

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
  for (k=1; k<=clus.N_MAX; k++) {
    if (star[k].binind == 0) { /* single star */
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
      
      /* evolve slightly (1 year) for initial radii */
      tphysf = 1.0e-6;
      dtp = tphysf;
      bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
		 &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
		 &(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
		 &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars);
      star[k].rad = star[k].se_radius * RSUN / units.l;
      star[k].m = star[k].se_mass * MSUN / units.mstar;
    } else if (star[k].binind > 0) { /* binary */
      star[k].se_k = NOT_A_STAR; /* just for safety */
      kb = star[k].binind;
      binary[kb].bse_mass[0] = binary[kb].m1 * units.mstar / MSUN;
      binary[kb].bse_mass[1] = binary[kb].m2 * units.mstar / MSUN;
      for (i=0; i<=1; i++) {
	if(binary[kb].bse_mass[i] <= 0.7){
	  binary[kb].bse_kw[i] = 0;
	} else {
	  binary[kb].bse_kw[i] = 1;
	}
	binary[kb].bse_mass0[i] = binary[kb].bse_mass[i];
	binary[kb].bse_ospin[i] = 0.0;
	binary[kb].bse_epoch[i] = 0.0;
      }
      binary[kb].bse_tphys = 0.0;
      
      /* set binary orbital period (in days) from a */
      binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;
      
      /* evolve slightly (1 year) for initial radii */
      tphysf = 1.0e-6;
      dtp = tphysf;
      bse_evolv2(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
		 &(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
		 &(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]), &(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
		 &(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
		 &(binary[kb].bse_tb), &(binary[kb].e));
      handle_bse_outcome(k, kb);
    } else {
      eprintf("totally confused!\n");
      exit_cleanly(-1);
    }
  }
}

/* note that this routine is called after perturb_stars() and get_positions() */
void do_stellar_evolution(gsl_rng *rng){
  long k, kb;
  int kprev;
  double dtp, tphysf, theta, vk;
  
  for(k=1; k<=clus.N_MAX; k++){
    if (star[k].binind == 0) { /* single star */
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
    } else { /* binary */
	tphysf = TotalTime / MEGA_YEAR;
	dtp = tphysf;
	kb = star[k].binind;

	/* set binary orbital period (in days) from a */
	binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;
	
	bse_evolv2(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
		   &(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
		   &(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]), &(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
		   &(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
		   &(binary[kb].bse_tb), &(binary[kb].e));
	handle_bse_outcome(k, kb);
    }
  }
}
void write_stellar_data(void){
  long k, kb;
  FILE *stel_file;
  char filename[1024];
  
  se_file_counter++;

  /* single star info */
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

  /* binary star info */
  sprintf(filename, "%s_binary_stellar_info.%05d", outprefix, se_file_counter);
  stel_file = fopen(filename, "w");
  if (stel_file==NULL){
    fprintf(stderr,
	    "file cannot be opened to write binary stellar info\n");
    return;
  }
  fprintf(stel_file, "# time (Myr): %e\n",
	  TotalTime/MEGA_YEAR);
  fprintf(stel_file, "# time (FP):  %e\n", TotalTime);
  fprintf(stel_file, "#1:id1 #2:id2 #3:m1 #4:m2 #5:R1 #6:R2 #7:kw1 #8:kw2 #9:tb #10:e\n");
  for(k=1; k<=clus.N_MAX; k++){
    if (star[k].binind) {
      kb = star[k].binind;
      fprintf(stel_file, "%08ld %08ld %g %g %g %g %d %d %g %g\n", 
	      binary[kb].id1, binary[kb].id2,
	      binary[kb].bse_mass[0], binary[kb].bse_mass[1],
	      binary[kb].bse_radius[0], binary[kb].bse_radius[1],
	      binary[kb].bse_kw[0], binary[kb].bse_kw[1],
	      binary[kb].bse_tb, binary[kb].e);
    }
  }
  fclose(stel_file);
}

void handle_bse_outcome(long k, long kb)
{
  long knew, knewp;

  if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0 && binary[kb].bse_tb > 0.0) {
    /* normal evolution */
    binary[kb].rad1 = binary[kb].bse_radius[0] * RSUN / units.l;
    binary[kb].rad2 = binary[kb].bse_radius[1] * RSUN / units.l;
    binary[kb].m1 = binary[kb].bse_mass[0] * MSUN / units.mstar;
    binary[kb].m2 = binary[kb].bse_mass[1] * MSUN / units.mstar;
    star[k].m = binary[kb].m1 + binary[kb].m2;
    binary[kb].a = pow((binary[kb].bse_mass[0]+binary[kb].bse_mass[1])*sqr(binary[kb].bse_tb/365.25), 1.0/3.0)
      * AU / units.l;
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0 && binary[kb].bse_tb <= 0.0) {
    /* disruption with both stars "intact" */
    dprintf("binary disrupted via BSE with both stars intact\n");
    knew = create_star();
    knewp = create_star();
    cp_binmemb_to_star(k, 0, knew);
    cp_binmemb_to_star(k, 1, knew);
    destroy_obj(k);
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] == 0.0 && binary[kb].bse_tb <= 0.0) {
    /* secondary star gone */
    dprintf("binary disrupted via BSE with first star intact\n");
    knew = create_star();
    cp_binmemb_to_star(k, 0, knew);
    destroy_obj(k);
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] != 0.0 && binary[kb].bse_tb <= 0.0) {
    /* primary star gone */
    dprintf("binary disrupted via BSE with second star intact\n");
    knew = create_star();
    cp_binmemb_to_star(k, 1, knew);
    destroy_obj(k);
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] == 0.0 && binary[kb].bse_tb <= 0.0) {
    /* both stars gone */
    dprintf("binary disrupted via BSE with no stars intact\n");
    destroy_obj(k);
  }
}

/* k=star index, kbi=0 or 1, knew=index of new star */
void cp_binmemb_to_star(long k, int kbi, long knew)
{
  long kb;

  kb = star[k].binind;
  /* and set the stars' dynamical properties */
  star[knew].r = star[k].r;
  star[knew].vr = star[k].vr;
  star[knew].vt = star[k].vt;
  star[knew].m = binary[kb].bse_mass[kbi] * MSUN / units.mstar;
  star[knew].phi = star[k].phi;
  set_star_EJ(knew);
  set_star_news(knew);
  set_star_olds(knew);
  /* mark stars as interacted so they don't undergo E_CONS mode stuff */
  star[knew].interacted = 1;
  if (kbi == 0) {
    star[knew].Eint = binary[kb].Eint1;
    star[knew].id = binary[kb].id1;
  } else {
    star[knew].Eint = binary[kb].Eint2;
    star[knew].id = binary[kb].id2;
  }
  star[knew].rad = binary[kb].bse_radius[kbi] * RSUN / units.l;
  star[knew].se_mass = binary[kb].bse_mass[kbi];
  star[knew].se_k = binary[kb].bse_kw[kbi];
  star[knew].se_mt = binary[kb].bse_mass0[kbi]; /* FIXME: is this right? */
  star[knew].se_ospin = binary[kb].bse_ospin[kbi];
  star[knew].se_epoch = binary[kb].bse_epoch[kbi];
  star[knew].se_tphys = binary[kb].bse_tphys;
  star[knew].se_radius = binary[kb].bse_radius[kbi];
  star[knew].se_lum = binary[kb].bse_lum[kbi];
  star[knew].se_mc = binary[kb].bse_massc[kbi];
  star[knew].se_rc = binary[kb].bse_radc[kbi];
  star[knew].se_menv = binary[kb].bse_menv[kbi];
  star[knew].se_renv = binary[kb].bse_renv[kbi];
  star[knew].se_tms = binary[kb].bse_tms[kbi];
}
