#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bse_wrap.h"

int main(void)
{
  int i, j, kw[2];
  double mass[2], mt[2], r[2], lum[2], mc[2], rc[2];
  double menv[2], renv[2], ospin[2], epoch[2], tms[2], tphys[2];
  double tphysf=12000.0, dtp=12000.0, z=0.02, *zpars;

  zpars = (double *) malloc(20 * sizeof(double));

  bse_set_neta(0.5);
  bse_set_bwind(0.0);
  bse_set_hewind(0.5);
  bse_set_sigma(190.0);
  bse_set_ifflag(0);
  bse_set_wdflag(1);
  bse_set_bhflag(0);
  bse_set_nsflag(1);
  bse_set_mxns(3.0);
  bse_set_pts1(0.05);
  bse_set_pts2(0.01);
  bse_set_pts3(0.02);
  bse_set_idum(0);

  bse_zcnsts(&z, zpars);
  
  mass[0] = 30.0;
  mass[1] = 1.0;
  kw[0] = 1;
  kw[1] = 1;
  mt[0] = mass[0];
  mt[1] = mass[1];
  ospin[0] = 0.0;
  ospin[1] = 0.0;
  epoch[0] = 0.0;
  epoch[1] = 0.0;
  tphys[0] = 0.0;
  tphys[1] = 0.0;

  bse_evolv1(&(kw[0]), &(mass[0]), &(mt[0]), &(r[0]), &(lum[0]), &(mc[0]), &(rc[0]), 
	     &(menv[0]), &(renv[0]), &(ospin[0]), &(epoch[0]), &(tms[0]), &(tphys[0]), 
	     &tphysf, &dtp, &z, zpars);
  
  i = 1;
  while (bse_get_spp(i, 1) >= 0.0) {
    fprintf(stdout, "star 0: type=%25s time=%f mass=%f radius=%f\n", 
	    bse_get_sselabel((int) bse_get_spp(i, 2)), 
	    bse_get_spp(i, 1), 
	    bse_get_spp(i, 3),
	    r[0]);
    i++;
  }
  
  bse_evolv1(&(kw[1]), &(mass[1]), &(mt[1]), &(r[1]), &(lum[1]), &(mc[1]), &(rc[1]), 
	     &(menv[1]), &(renv[1]), &(ospin[1]), &(epoch[1]), &(tms[1]), &(tphys[1]), 
	     &tphysf, &dtp, &z, zpars);
  
  i = 1;
  while (bse_get_spp(i, 1) >= 0.0) {
    fprintf(stdout, "star 1: type=%25s time=%f mass=%f radius=%f\n", 
	    bse_get_sselabel((int) bse_get_spp(i, 2)), 
	    bse_get_spp(i, 1), 
	    bse_get_spp(i, 3),
	    r[1]);
    i++;
  }

  free(zpars);

  return(0);
}
