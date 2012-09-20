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
  double menv[2], renv[2], ospin[2], B_0[2], bacc[2], tacc[2], epoch[2], tms[2], tphys;
  double rad[2], lum[2], massc[2], radc[2];
  double tphysf, dtp, z=0.02, *zpars;
  double vs[12];
  double aj, tm, tn, tscls[20], lums[10], GB[10], k2;
  double a;

  zpars = (double *) malloc(20 * sizeof(double));

  /* evolution parameters */
  bse_set_neta(0.5);
  bse_set_bwind(0.0);
  bse_set_hewind(1.0);
  bse_set_windflag(1);
  bse_set_alpha1(3.0);
  bse_set_lambda(0.5);
  bse_set_ceflag(0);
  bse_set_tflag(1);
  bse_set_ifflag(0);
  bse_set_wdflag(1);
  bse_set_bhflag(1);
  bse_set_nsflag(1);
  bse_set_mxns(3.0);
  bse_set_bconst(-3000.0);
  bse_set_CK(-1000.0);
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
  mass0[0] = 1.3;
  mass0[1] = 3.0;
  mass[0] = mass0[0];
  mass[1] = mass0[1];
  kw[0] = 1;
  kw[1] = 1;
  tb = 1.0e5;
  ecc = 0.0;
  ospin[0] = 0.0;
  ospin[1] = 0.0;
  B_0[0] = 0.0;
  B_0[1] = 0.0;
  bacc[0] = 0.0;
  bacc[1] = 0.0;
  tacc[0] = 0.0;
  tacc[1] = 0.0;
  epoch[0] = 0.0;
  epoch[1] = 0.0;
  tphys = 0.0;
  tphysf = 100.0;
  dtp = 10.0;

  while (kw[0]<=1) {
    bse_evolv2(&(kw[0]), &(mass0[0]), &(mass[0]), &(rad[0]), &(lum[0]), &(massc[0]), &(radc[0]), 
	       &(menv[0]), &(renv[0]), &(ospin[0]), &(B_0[0]), &(bacc[0]), &(tacc[0]), &(epoch[0]), &(tms[0]), 
	       &tphys, &tphysf, &dtp, &z, zpars, &tb, &ecc, vs);
    
    j = 1;
    while (bse_get_bpp(j, 1) >= 0.0) {
      fprintf(stdout, "time=%g m1=%g m2=%g k1=%d k2=%d sep=%g ecc=%g r1/rol1=%g r2/rol2=%g type=%s\n",
	      bse_get_bpp(j, 1), bse_get_bpp(j, 2), bse_get_bpp(j, 3), (int) bse_get_bpp(j, 4), (int) bse_get_bpp(j, 5),
	      bse_get_bpp(j, 6), bse_get_bpp(j, 7), bse_get_bpp(j, 8), bse_get_bpp(j, 9),
	      bse_get_bselabel((int) bse_get_bpp(j, 10)));
      
      j++;
    }
    tphysf += 100.0;
  }

  /* fprintf(stdout, "mc=%g rc=%g\n", massc[1], radc[1]); */

  /* strip envelope of giant */
  /* fprintf(stdout, "stripping giant's envelope...\n"); */
/*   aj = tphys - epoch[1]; */
/*   mass[1] = massc[1]; */
/*   bse_star(&(kw[1]), &(mass0[1]), &(mass[1]), &tm, &tn, tscls, lums, GB, zpars); */
/*   bse_hrdiag(&(mass0[1]), &aj, &(mass[1]), &tm, &tn, tscls, lums, GB, zpars, */
/* 	     &(rad[1]), &(lum[1]), &(kw[1]), &(massc[1]), &(radc[1]), &(menv[1]), &(renv[1]), &k2); */
/*   epoch[1] = tphys - aj; */
/*   tphysf = 12.0e3; */

/*   fprintf(stdout, "m=%g r=%g\n", mass[1], rad[1]); */

  /* shrink orbit */
  /* fprintf(stdout, "shrinking orbit...\n"); */
  tphysf = 12.0e3;

  fprintf(stdout, "rad0=%g, rad1=%g\n", rad[0], rad[1]);

  a = BSE_WRAP_MAX(10.0*rad[0]/bse_rl(mass[0]/mass[1]), 10.0*rad[1]/bse_rl(mass[1]/mass[0]));
  /* convert to AU */
  a /= 214.9456;
  
  fprintf(stdout, "new a=%g AU RL1/a=%g RL2/a=%g\n", a, bse_rl(mass[0]/mass[1]), bse_rl(mass[1]/mass[0]));

  tb = 365.25 * sqrt(a*a*a/(mass[0]+mass[1]));
  ecc = 1.0 - 0.1*(rad[0]+rad[1])/214.9456/a;

  bse_evolv2(&(kw[0]), &(mass0[0]), &(mass[0]), &(rad[0]), &(lum[0]), &(massc[0]), &(radc[0]), 
	     &(menv[0]), &(renv[0]), &(ospin[0]), &(B_0[0]), &(bacc[0]), &(tacc[0]), &(epoch[0]), &(tms[0]), 
	     &tphys, &tphysf, &dtp, &z, zpars, &tb, &ecc, vs);

  j = 1;
  while (bse_get_bpp(j, 1) >= 0.0) {
    fprintf(stdout, "time=%g m1=%g m2=%g k1=%d k2=%d sep=%g ecc=%g r1/rol1=%g r2/rol2=%g type=%s\n",
	    bse_get_bpp(j, 1), bse_get_bpp(j, 2), bse_get_bpp(j, 3), (int) bse_get_bpp(j, 4), (int) bse_get_bpp(j, 5),
	    bse_get_bpp(j, 6), bse_get_bpp(j, 7), bse_get_bpp(j, 8), bse_get_bpp(j, 9),
	    bse_get_bselabel((int) bse_get_bpp(j, 10)));

    j++;
  }

  /* fprintf(stdout, "m1=%f m2=%f tb=%f e=%f\n", mass[0], mass[1], tb, ecc); */
  
  /* fprintf(stdout, "vs=%g\n", sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3])); */

  free(zpars);

  return(0);
}
