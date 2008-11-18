/* c_bse.c

   Copyright (C) 2008 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

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
  mass0[0] = 21.0;
  mass0[1] = 3.0;
  mass[0] = mass0[0];
  mass[1] = mass0[1];
  kw[0] = 1;
  kw[1] = 1;
  tb = 1.0e3;
  ecc = 0.1;
  ospin[0] = 0.0;
  ospin[1] = 0.0;
  epoch[0] = 0.0;
  epoch[1] = 0.0;
  tphys = 0.0;
  tphysf = 1.0e4;
  dtp = 10.0;

  bse_evolv2(&(kw[0]), &(mass0[0]), &(mass[0]), &(rad[0]), &(lum[0]), &(massc[0]), &(radc[0]), 
	     &(menv[0]), &(renv[0]), &(ospin[0]), &(epoch[0]), &(tms[0]), 
	     &tphys, &tphysf, &dtp, &z, zpars, &tb, &ecc, vs);

/*   tphysf = 100.0; */
/*   dtp = 0.0; */

/*   bse_evolv2(&(kw[0]), &(mass0[0]), &(mass[0]), &(rad[0]), &(lum[0]), &(massc[0]), &(radc[0]),  */
/* 	     &(menv[0]), &(renv[0]), &(ospin[0]), &(epoch[0]), &(tms[0]),  */
/* 	     &tphys, &tphysf, &dtp, &z, zpars, &tb, &ecc, vs); */

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

  /* print additional info stored only in bcm array */
  j = 1;
  while (bse_get_bcm(j, 1) >= 0.0) {
    /* fprintf(stdout, "j=%d t=%g dm/dt=%g %g rad/rol=%g %g\n", j, bse_get_bcm(j, 1),
	    bse_get_bcm(j,14), bse_get_bcm(j,28), 
	    bse_get_bcm(j,15), bse_get_bcm(j,29)); */
    j++;
  }
  j--;
  if (j >= 1) {
    fprintf(stdout, "j=%d dm/dt=%g %g rad/rol=%g %g\n", j, 
	    bse_get_bcm(j,14), bse_get_bcm(j,28), 
	    bse_get_bcm(j,15), bse_get_bcm(j,29));
  }
  
  fprintf(stdout, "m1=%f m2=%f tb=%f e=%f\n", mass[0], mass[1], tb, ecc);
  
  fprintf(stdout, "vs=%g\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));

  free(zpars);

  return(0);
}
