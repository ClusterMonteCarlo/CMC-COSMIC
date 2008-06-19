#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bse_wrap.h"

int main(void)
{
  int i, kw, counter;
  double mass, mt, rad, lum, mc, rc;
  double menv, renv, ospin, epoch, tms, tphys;
  double tphysf, dtp, z=0.001, *zpars;
  double vs[3];

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
  
  mass = 14.0;
  mt = mass;
  kw = 1;
  ospin = 0.0;
  epoch = 0.0;
  tphys = 0.0;

  counter = 0;
  tphysf = 0.0;
  while (kw <= 6 && tphys < 15000.0) {
    counter ++;
    tphysf += ((double) counter) * 0.1;
    dtp = 0.0;

    bse_evolv1(&kw, &mass, &mt, &rad, &lum, &mc, &rc, 
	       &menv, &renv, &ospin, &epoch, &tms, &tphys, 
	       &tphysf, &dtp, &z, zpars, vs);
    
    fprintf(stdout, "mass=%f mt=%f vs=%f\n", mass, mt, sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));

    i = 1;
    while (bse_get_spp(i, 1) >= 0.0) {
      fprintf(stdout, "type=%25s time=%f mass=%f radius=%f\n", 
	      bse_get_sselabel((int) bse_get_spp(i, 2)), 
	      bse_get_spp(i, 1), 
	      bse_get_spp(i, 3),
	      rad);
      i++;
    }
  }
  
  free(zpars);

  return(0);
}
