#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "bse_wrap.h"
#include "popsynth.h"

void printoutput(int *filectr, binary_t *binarray, double time);

int main(void) {
  int filectr=0;
  unsigned long int seed=234UL;
  binary_t *binarray;
  int i, j, idctr;
  double X, norm, mmin, mmax, index, mass;
  gsl_rng *rng;
  const gsl_rng_type *rng_type=gsl_rng_mt19937;
  double z=0.03, *zpars, vs[12];
  double amin, amax, targettphysf, tphysf, dtp, rlprimovera;
  
  /* initialize GSL rng */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(rng, seed);
  
  binarray = (binary_t *) malloc(NBIN * sizeof(binary_t));
  
  /* initialize stellar evolution */
  zpars = (double *) malloc(20 * sizeof(double));

  /* evolution parameters */
  bse_set_neta(0.5);
  bse_set_bwind(0.0);
  bse_set_hewind(1.0);
  bse_set_alpha1(1.0);
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
  bse_set_sigma(265.0);
  bse_set_beta(0.125);
  bse_set_xi(1.0);
  bse_set_acc2(1.5);
  bse_set_epsnov(0.001);
  bse_set_eddfac(10.0);
  bse_set_gamma(-1.0);
  
  bse_zcnsts(&z, zpars);
  bse_instar();

  idctr = 0;
  for (i=0; i<NBIN; i++) {
    /* set arbitrary but sensible IDs */
    binarray[i].id1 = idctr++;
    binarray[i].id2 = idctr++;

    /* zero out some vars */
    binarray[i].bse_ospin[0] = 0.0;
    binarray[i].bse_ospin[1] = 0.0;
    binarray[i].bse_epoch[0] = 0.0;
    binarray[i].bse_epoch[1] = 0.0;
    binarray[i].bse_tphys = 0.0;

    /* set primary mass from power law */
    X = gsl_rng_uniform(rng);
    mmin = 6.0;
    mmax = 150.0;
    index = -2.35;
    norm = pow(mmax/mmin, index+1.0) - 1.0;
    mass = mmin*pow(norm*X+1, 1.0/(index+1.0));

    /* set secondary mass from uniform distribution in mass ratio */
    binarray[i].bse_mass0[0] = mass;
    binarray[i].bse_mass0[1] = gsl_rng_uniform(rng) * (mass - 0.08) + 0.08;
    
    binarray[i].bse_mass[0] = binarray[i].bse_mass0[0];
    binarray[i].bse_mass[1] = binarray[i].bse_mass0[1];
    
    binarray[i].bse_kw[0] = 1;
    binarray[i].bse_kw[1] = 1;

    /* evolve each star slightly to set radius */
    dtp = 0.0;
    binarray[i].bse_tphys = 0.0;
    tphysf = 1.0;
    bse_evolv1(&(binarray[i].bse_kw[0]), &(binarray[i].bse_mass0[0]), 
	       &(binarray[i].bse_mass[0]), &(binarray[i].bse_radius[0]), &(binarray[i].bse_lum[0]), 
	       &(binarray[i].bse_massc[0]), &(binarray[i].bse_radc[0]), 
	       &(binarray[i].bse_menv[0]), &(binarray[i].bse_renv[0]), 
	       &(binarray[i].bse_ospin[0]), &(binarray[i].bse_epoch[0]), 
	       &(binarray[i].bse_tms[0]), 
	       &binarray[i].bse_tphys, 
	       &tphysf, &dtp, &z, zpars, vs);

    dtp = 0.0;
    binarray[i].bse_tphys = 0.0;
    tphysf = 1.0;
    bse_evolv1(&(binarray[i].bse_kw[1]), &(binarray[i].bse_mass0[1]), 
	       &(binarray[i].bse_mass[1]), &(binarray[i].bse_radius[1]), &(binarray[i].bse_lum[1]), 
	       &(binarray[i].bse_massc[1]), &(binarray[i].bse_radc[1]), 
	       &(binarray[i].bse_menv[1]), &(binarray[i].bse_renv[1]), 
	       &(binarray[i].bse_ospin[1]), &(binarray[i].bse_epoch[1]), 
	       &(binarray[i].bse_tms[1]), 
	       &binarray[i].bse_tphys, 
	       &tphysf, &dtp, &z, zpars, vs);
    
    /* set semimajor axis uniform in log from primary radius at half Roche lobe to 10^5 R_sun, 
       eccentricity thermal */
    rlprimovera = 0.46224 * pow(binarray[i].bse_mass[0]/(binarray[i].bse_mass[0]+binarray[i].bse_mass[1]), 1.0/3.0);
    amin = 2.0 * binarray[i].bse_radius[0] / rlprimovera * RSUN / AU;
    amax = 1.0e5 * RSUN / AU;
    binarray[i].a = pow(10.0, gsl_rng_uniform(rng)*(log10(amax)-log10(amin))+log10(amin));
    binarray[i].bse_tb = 365.25 * sqrt(pow(binarray[i].a, 3.0)/(binarray[i].bse_mass0[0]+binarray[i].bse_mass0[1]));
    binarray[i].e = sqrt(gsl_rng_uniform(rng));
  }

  /* loop over time */
  for (targettphysf=9.0e3; targettphysf<=10.0e3; targettphysf+=1.0) {
    for (i=0; i<NBIN; i++) { 
      tphysf = targettphysf;
      dtp = 0.0;
      bse_evolv2(&(binarray[i].bse_kw[0]), &(binarray[i].bse_mass0[0]), &(binarray[i].bse_mass[0]), 
		 &(binarray[i].bse_radius[0]), &(binarray[i].bse_lum[0]), 
		 &(binarray[i].bse_massc[0]), &(binarray[i].bse_radc[0]), 
		 &(binarray[i].bse_menv[0]), &(binarray[i].bse_renv[0]), 
		 &(binarray[i].bse_ospin[0]), &(binarray[i].bse_epoch[0]), 
		 &(binarray[i].bse_tms[0]), 
		 &binarray[i].bse_tphys, &tphysf, &dtp, 
		 &z, zpars, 
		 &binarray[i].bse_tb, &binarray[i].e, vs);
      
      /* extract some binary info from BSE's bcm array */
      j = 1;
      while (bse_get_bcm(j, 1) >= 0.0) {
	j++;
      }
      j--;
      if (j >= 1) {
	if (fabs((binarray[i].bse_tphys - bse_get_bcm(j,1))/binarray[i].bse_tphys) >= 1.0e-6) {
	  fprintf(stderr, "binarray[kb].bse_tphys=%g bcmtime=%g\n", binarray[i].bse_tphys, bse_get_bcm(j,1));
	  /* exit_cleanly(-1); */
	}
	binarray[i].bse_bcm_dmdt[0] = bse_get_bcm(j, 14);
	binarray[i].bse_bcm_dmdt[1] = bse_get_bcm(j, 28);
	binarray[i].bse_bcm_radrol[0] = bse_get_bcm(j, 15);
	binarray[i].bse_bcm_radrol[1] = bse_get_bcm(j, 29);
      } else {
	fprintf(stderr, "Could not extract BSE bcm info!  Input dtp not exactly equal to tphysf-tphys?\n");
	exit(1);
      }
    }

    printoutput(&filectr, binarray, tphysf);
  }

  /* free GSL stuff */
  gsl_rng_free(rng);

  free(binarray);

  return(0);
}

void printoutput(int *filectr, binary_t *binarray, double time)
{
  int i;
  FILE *ofp;
  char ofpname[1024];

  sprintf(ofpname, "popsynth.%05d.dat.gz", *filectr);
  if ((ofp = (FILE *) gzopen(ofpname, "wb")) == NULL) {
    fprintf(stderr, "cannot create output file %s\n", ofpname);
    exit(1);
  }
  
  gzprintf(ofp, "#t=%g\n", time);
  gzprintf(ofp, "#1:id1 #2:id2 #3:m1 #4:m2 #5:k1 #6:k2 #7:tb #8:e #9:rad1 #10:rad2 #11:lum1 #12:lum2 #13:massc1 #14:massc2 #15:radc1 #16:radc2 #17:menv1 #18:menv2 #19:renv1 #20:renv2 #21:ospin1 #22:ospin2 #23:tms1 #24:tms2 #25:dmdt1 #26:dmdt2 #27:radrol1 #28:radrol2\n");
  
  for (i=0; i<NBIN; i++) {
    if ( binarray[i].bse_tb > 0.0 && 
	 ( (binarray[i].bse_kw[0] >= 13 && binarray[i].bse_kw[0] <= 14) || 
	   (binarray[i].bse_kw[1] >= 13 && binarray[i].bse_kw[1] <= 14)) &&
	 (binarray[i].bse_bcm_dmdt[0] != 0.0 || binarray[i].bse_bcm_dmdt[1] != 0.0) ) {
      gzprintf(ofp, "%ld %ld %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
	       binarray[i].id1, binarray[i].id2,
	       binarray[i].bse_mass[0], binarray[i].bse_mass[1],
	       binarray[i].bse_kw[0], binarray[i].bse_kw[1],
	       binarray[i].bse_tb, binarray[i].e,
	       binarray[i].bse_radius[0], binarray[i].bse_radius[1],
	       binarray[i].bse_lum[0], binarray[i].bse_lum[1],
	       binarray[i].bse_massc[0], binarray[i].bse_massc[1],
	       binarray[i].bse_radc[0], binarray[i].bse_radc[1],
	       binarray[i].bse_menv[0], binarray[i].bse_menv[1],
	       binarray[i].bse_renv[0], binarray[i].bse_renv[1],
	       binarray[i].bse_ospin[0], binarray[i].bse_ospin[1],
	       binarray[i].bse_tms[0], binarray[i].bse_tms[1],
	       binarray[i].bse_bcm_dmdt[0], binarray[i].bse_bcm_dmdt[1],
	       binarray[i].bse_bcm_radrol[0], binarray[i].bse_bcm_radrol[1]);
    }
  }

  gzclose(ofp);
  
  (*filectr)++;
}
