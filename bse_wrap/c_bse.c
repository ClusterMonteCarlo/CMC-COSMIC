#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bse_wrap.h"

int main(void)
{
  int i, j, kw[2];
  double mass0[2], mass[2], tb, ecc;
  double menv[2], renv[2], ospin[2], epoch[2], tms[2], tphys;
  double rad[2], lum[2], massc[2], radc[2];
  double tphysf, dtp, z=0.02, *zpars;
  double vs[3];

  zpars = (double *) malloc(20 * sizeof(double));

  /* evolution parameters */
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
  
  bse_zcnsts(&z, zpars);
  bse_instar();
  mass0[0] = 120.0;
  mass0[1] = 75.0;
  mass[0] = mass0[0];
  mass[1] = mass0[1];
  kw[0] = 1;
  kw[1] = 1;
  tb = 1.0e-19;
  ecc = 0.9999;
  ospin[0] = 0.0;
  ospin[1] = 0.0;
  epoch[0] = 0.0;
  epoch[1] = 0.0;
  tphys = 0.0;
  tphysf = 1.0;

  /* a particularly troublesome binary merger */
  /* z = 0.02; */
/*   bse_zcnsts(&z, zpars); */
/*   bse_instar(); */
/*   tphys = 6.84714285976709736; */
/*   tphysf = tphys + 1.0e-6; */
/*   dtp = 0.0; */
/*   mass0[0] = 58.9425994842380163; */
/*   mass0[1] = 100.281509857595665; */
/*   kw[0] = 14; */
/*   kw[1] = 1; */
/*   mass[0] = 58.9425994842380163; */
/*   mass[1] = 100.281509857595665; */
/*   ospin[0] = 200000000; */
/*   ospin[1] = 3095.38442913478411; */
/*   epoch[0] = 6.74714285976709771; */
/*   epoch[1] = 1.32731978583074994; */
/*   rad[0] = 0.000249916621813169189; */
/*   rad[1] = 15.8986730963320575; */
/*   lum[0] = 1.00000000000000004e-10; */
/*   lum[1] = 1897127.10106106009; */
/*   massc[0] = 58.9425994842380163; */
/*   massc[1] = 0.0; */
/*   radc[0] = 0.000249916621813169189; */
/*   radc[1] = 0.0; */
/*   menv[0] = 1.00000000000000004e-10; */
/*   menv[1] = 1.00000000000000004e-10; */
/*   renv[0] = 1.00000000000000004e-10; */
/*   renv[1] = 1.00000000000000004e-10; */
/*   tms[0] = 10000000000; */
/*   tms[1] = 3.50096074226224108; */
/*   tb = 2.89458168924826573e-17; */
/*   ecc = 0.998999999999999999; */
  //tb = 1.0e6;
  //ecc = 0.0;

  bse_evolv2(&(kw[0]), &(mass0[0]), &(mass[0]), &(rad[0]), &(lum[0]), &(massc[0]), &(radc[0]), 
	     &(menv[0]), &(renv[0]), &(ospin[0]), &(epoch[0]), &(tms[0]), 
	     &tphys, &tphysf, &dtp, &z, zpars, &tb, &ecc, vs);

  fprintf(stdout, "star 0: mass0=%f mass=%f tms=%g epoch=%g massc=%g rad=%g\n", mass0[0], mass[0], tms[0], epoch[0], massc[0], rad[0]);
  fprintf(stdout, "star 1: mass0=%f mass=%f tms=%g epoch=%g massc=%g rad=%g\n", mass0[1], mass[1], tms[1], epoch[1], massc[1], rad[1]);

  j = 1;
  while (bse_get_bpp(j, 1) >= 0.0) {
    fprintf(stdout, "time=%g m1=%g m2=%g k1=%d k2=%d sep=%g ecc=%g r1/rol1=%g r2/rol2=%g type=%s\n",
	    bse_get_bpp(j, 1), bse_get_bpp(j, 2), bse_get_bpp(j, 3), (int) bse_get_bpp(j, 4), (int) bse_get_bpp(j, 5),
	    bse_get_bpp(j, 6), bse_get_bpp(j, 7), bse_get_bpp(j, 8), bse_get_bpp(j, 9),
	    bse_get_bselabel((int) bse_get_bpp(j, 10)));

    j++;
  }
  
  fprintf(stdout, "m1=%f m2=%f tb=%f e=%f\n", mass[0], mass[1], tb, ecc);
  
  fprintf(stdout, "vs=%g\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));

  free(zpars);

  return(0);
}
