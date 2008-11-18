/* bse_wrap.h

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

/* prototypes for fortran BSE functions */
void zcnsts_(double *z, double *zpars);
void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum,
	     double *mc, double *rc, double *menv, double *renv, double *ospin,
	     double *epoch, double *tms, double *tphys, double *tphysf, 
	     double *dtp, double *z, double *zpars, double *vs);
void evolv2_(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
	     double *massc, double *radc, double *menv, double *renv, double *ospin,
	     double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
	     double *z, double *zpars, double *tb, double *ecc, double *vs);
void instar_(void);
float ran3_(int *idum);

/* wrapped BSE functions */
void bse_zcnsts(double *z, double *zpars);
void bse_evolv1(int *kw, double *mass, double *mt, double *r, double *lum,
		double *mc, double *rc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf, 
		double *dtp, double *z, double *zpars, double *vs);
void bse_evolv1_safely(int *kw, double *mass, double *mt, double *r, double *lum,
		       double *mc, double *rc, double *menv, double *renv, double *ospin,
		       double *epoch, double *tms, double *tphys, double *tphysf, 
		       double *dtp, double *z, double *zpars, double *vs);
void bse_evolv2(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		double *massc, double *radc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		double *z, double *zpars, double *tb, double *ecc, double *vs);
void bse_evolv2_safely(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		       double *massc, double *radc, double *menv, double *renv, double *ospin,
		       double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		       double *z, double *zpars, double *tb, double *ecc, double *vs);
void bse_instar(void);

/* structs to access BSE common blocks */
/* note the index swap between fortran and C: i,j->j,i */
extern struct { int idum; } value3_;
extern struct { int idum2, iy, ir[32]; } rand3_;
extern struct { int ktype[15][15]; } types_;
extern struct { int ceflag, tflag, ifflag, nsflag, wdflag; } flags_;
extern struct { double neta, bwind, hewind, mxns; } value1_;
extern struct { double alpha1, lambda; } value2_;
extern struct { double sigma; int bhflag; } value4_;
extern struct { double beta, xi, acc2, epsnov, eddfac, gamma; } value5_;
extern struct { double pts1, pts2, pts3; } points_;
extern struct { double dmmax, drmax; } tstepc_;
extern struct { float scm[14][50000], spp[3][20]; } single_;
extern struct { float bcm[34][50000], bpp[10][80]; } binary_;

/* setters */
void bse_set_idum(int idum); /* RNG seed (for NS birth kicks) */
void bse_set_neta(double neta); /* Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally) */
void bse_set_bwind(double bwind); /* binary enhanced mass loss parameter (inactive for single) */
void bse_set_hewind(double hewind); /* helium star mass loss factor (1.0 normally) */
void bse_set_sigma(double sigma); /* dispersion in the Maxwellian for the SN kick speed (190 km/s) */
void bse_set_ifflag(int ifflag); /* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0) */
void bse_set_wdflag(int wdflag); /* wdflag > 0 uses modified-Mestel cooling for WDs (0) */
void bse_set_bhflag(int bhflag); /* bhflag > 0 allows velocity kick at BH formation (0) */
void bse_set_nsflag(int nsflag); /* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_mxns(double mxns); /* maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1) */
void bse_set_pts1(double pts1); /* timestep taken in MS phase (0.05) */
void bse_set_pts2(double pts2); /* timestep taken in GB, CHeB, AGB, HeGB phases (0.01) */
void bse_set_pts3(double pts3); /* timestep taken in HG, HeMS phases (0.02) */
void bse_set_alpha1(double alpha1); /* common-envelope efficiency parameter (1.0) */
void bse_set_lambda(double lambda); /* binding energy factor for common envelope evolution (0.5) */
void bse_set_ceflag(int ceflag); /* ceflag > 0 activates spin-energy correction in common-envelope (0); ceflag = 3 activates de Kool common-envelope model (0) */
void bse_set_tflag(int tflag); /* tflag > 0 activates tidal circularisation */
void bse_set_beta(double beta); /* wind velocity factor: proportional to vwind**2 (1/8) */
void bse_set_xi(double xi); /* wind accretion efficiency factor (1.0) */
void bse_set_acc2(double acc2); /* Bondi-Hoyle wind accretion factor (3/2) */
void bse_set_epsnov(double epsnov); /* fraction of accreted matter retained in nova eruption (0.001) */
void bse_set_eddfac(double eddfac); /* Eddington limit factor for mass transfer (1.0) */
void bse_set_gamma(double gamma); /* angular momentum factor for mass lost during Roche (-1.0) */

/* getters */
float bse_get_spp(int i, int j); /* stellar evolution log */
float bse_get_scm(int i, int j); /* stored stellar parameters at interval dtp */
float bse_get_bpp(int i, int j); /* binary evolution log */
float bse_get_bcm(int i, int j); /* stored binary parameters at interval dtp */
char *bse_get_sselabel(int kw); /* converts stellar type number to text label */
char *bse_get_bselabel(int kw); /* converts binary type number to text label */

/* copied functions */
double bse_kick_speed(int *startype); /* routine for generating birth kick speed from distribution */

/* useful macros */
#define BSE_WRAP_MAX(a, b) ((a)>=(b)?(a):(b))
#define BSE_WRAP_MIN(a, b) ((a)<=(b)?(a):(b))

