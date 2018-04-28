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
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/taus113-v2.h"
#include "../config.h"

/**
* @brief A structure used to pass binary information along to bse_wrap.c
*/
typedef struct{
/**
* @brief semimajor axis
*/
	double a;
/**
* @brief eccentricity
*/
	double e;
/**
* @brief star types
*/
	int bse_kw[2];
/**
* @brief initial masses
*/
	double bse_mass0[2];
/**
* @brief masses
*/
	double bse_mass[2];
/**
* @brief radii
*/
	double bse_radius[2];
/**
* @brief luminosity
*/
	double bse_lum[2];
/**
* @brief ?
*/
	double bse_massc[2];
/**
* @brief ?
*/
	double bse_radc[2];
/**
* @brief ?
*/
	double bse_menv[2];
/**
* @brief ?
*/
	double bse_renv[2];
/**
* @brief original spin
*/
	double bse_ospin[2];
/**
* @brief Pulsar magnetic field
*/
    double bse_B_0[2];
/**
* @brief Amount of mass pulsar has accreted
*/
    double bse_bacc[2];
/**
* @brief Amount of time pulsar has spent accreting
*/
    double bse_tacc[2];
/**
* @brief ?
*/
	double bse_epoch[2];
/**
* @brief ?
*/
	double bse_tms[2];
/**
* @brief physical time
*/
	double bse_tphys;
/**
* @brief binary orbital period
*/
	double bse_tb;
/**
* @brief mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)]
*/
	double bse_bcm_dmdt[2];
/**
* @brief the black hole spins, taken from the CO core mass 
*/
	double bse_bhspin[2];
/**
* @brief radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)]
*/
	double bse_bcm_radrol[2];
/**
* @brief Pulsar magnetic field strength at surface
*/
    double bse_bcm_B[2];
/**
* @brief provides formation pathway of NS
*/
    int bse_bcm_formation[2];
} bse_binary;

/* prototypes for fortran BSE functions */
void zcnsts_(double *z, double *zpars);
void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum,
	     double *mc, double *rc, double *menv, double *renv, double *ospin,
	     double *epoch, double *tms, double *tphys, double *tphysf, 
	     double *dtp, double *z, double *zpars, double *vs);
void evolv2_(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
	     double *massc, double *radc, double *menv, double *renv, double *ospin,
             double *B_0, double *bacc, double *tacc,
	     double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
	     double *z, double *zpars, double *tb, double *ecc, double *vs, double* bhspin);
void instar_(void);
float ran3_(int *idum);
void star_(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, 
	   double *lums, double *GB, double *zpars);
void hrdiag_(double *mass, double *aj, double *mt, double *tm, double *tn, double *tscls, 
	     double *lums, double *GB, double *zpars, double *r, double *lum, int *kw, 
	     double *mc, double *rc, double *menv, double *renv, double *k2, int *ST_tide, double *ecsnp, double *ecsn_mlow, double *bhspin);
void kick_(int *kw, double *m1, double *m1n, double *m2, double *ecc, double *sep, 
	   double *jorb, double *vk, int *snstar, double *r2, double *fallback, double *vs);
void mix_(double *mass, double *mt, double *aj, int *kw, double *zpars, double *ecsnp, double *bhspin);
// note: these function names only work if in lowercase here, even though FORTRAN versions in uppercase.
void comenv_(double *M01, double *M1, double *MC1, double *AJ1, double *JSPIN1, int *KW1, double *M02, double *M2, double *MC2, double *AJ2, double *JSPIN2, int *KW2, double *ZPARS, double *ECC, double *SEP, double *JORB, int *COEL, int *star1, int *star2, double *vk, int *fb, double *bkick, double *ecsnp, double *ecsn_mlow, int *formation1, int *formation2, int *ST_tide, double *bhspin1, double *bhspin2);

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
void bse_evolve_single(int *kw, double *mass, double *mt, double *r, double *lum,
		double *mc, double *rc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf,
		double *dtp, double *z, double *zpars, double *vs, double *bhspin);
void bse_evolv2(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		double *massc, double *radc, double *menv, double *renv, double *ospin,
                double *B_0, double *bacc, double *tacc,
		double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		double *z, double *zpars, double *tb, double *ecc, double *vs, double *bhspin);
void bse_evolv2_safely(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		       double *massc, double *radc, double *menv, double *renv, double *ospin,
                       double *B_0, double *bacc, double *tacc,
		       double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		       double *z, double *zpars, double *tb, double *ecc, double *vs, double *bhspin);
void bse_instar(void);
void bse_star(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, 
	      double *lums, double *GB, double *zpars);
void bse_hrdiag(double *mass, double *aj, double *mt, double *tm, double *tn, double *tscls, 
		double *lums, double *GB, double *zpars, double *r, double *lum, int *kw, 
		double *mc, double *rc, double *menv, double *renv, double *k2, int *ST_tide, double *ecsnp, double *ecsn_mlow, double *bhspin);
void bse_kick(int *kw, double *m1, double *m1n, double *m2, double *ecc, double *sep, 
	      double *jorb, double *vk, int *snstar, double *r2, double *fallback, double *vs);
void bse_mix(double *mass, double *mt, double *aj, int *kw, double *zpars, double *ecsnp, double *bhspin);
void bse_comenv(bse_binary *binary, double *zpars,
                double *vs, int *fb, double *ecsnp, double *ecsn_mlow, int *ST_tide);

/* structs to access BSE common blocks */
/* note the index swap between fortran and C: i,j->j,i */
extern struct { int idum; } value3_;
extern struct { int idum2, iy, ir[32]; } rand3_;
#ifdef USE_TAUS
extern struct { long long int state[4]; int first;} taus113state_;
#endif
extern struct { int ktype[15][15]; } types_;
extern struct { int ceflag, tflag, ifflag, nsflag, wdflag; } flags_;
extern struct { double neta, bwind, hewind, mxns; int windflag; int bhspinflag; double bhspinmag; int ppsn; } value1_;
extern struct { double alpha1, lambda; } value2_;
extern struct { double sigma; double bhsigmafrac; double bconst; double CK; int bhflag; int opening_angle; } value4_;
extern struct { double beta, xi, acc2, epsnov, eddfac, gamma; } value5_;
extern struct { double pts1, pts2, pts3; } points_;
extern struct { double dmmax, drmax; } tstepc_;
extern struct { float scm[14][50000], spp[3][20]; } single_;
extern struct { float bcm[36][50000], bpp[10][80]; } binary_;
extern struct { double merger; long int id1_pass; long int id2_pass; } cmcpass_;

/* setters */
void bse_set_idum(int idum); /* RNG seed (for NS birth kicks) */
void bse_set_neta(double neta); /* Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally) */
void bse_set_bwind(double bwind); /* binary enhanced mass loss parameter (inactive for single) */
void bse_set_hewind(double hewind); /* helium star mass loss factor (1.0 normally) */
void bse_set_windflag(int windflag); /* Sets wind prescription (0=BSE, 1=StarTrack, 2=Vink; 0) */
void bse_set_ppsn(int ppsn); /* Sets Pair-instability pulsations and supernoa */ 
void bse_set_sigma(double sigma); /* dispersion in the Maxwellian for the SN kick speed (190 km/s) */
void bse_set_bhsigmafrac(double bhsigmafrac); /* Ad hoc factor to change BH SN kick speed relative to NS SN kick sigma (1) */
void bse_set_opening_angle(int opening_angle); /* Switch to set the allowed opening angle of SN kicks.  Defaults to 180 degrees*/
void bse_set_ifflag(int ifflag); /* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0) */
void bse_set_wdflag(int wdflag); /* wdflag > 0 uses modified-Mestel cooling for WDs (0) */
void bse_set_bhflag(int bhflag); /* bhflag > 0 allows velocity kick at BH formation (0) */
void bse_set_nsflag(int nsflag); /* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_bhspinflag(int bhflag);/* bhspinflag (0=[bhspinmag], 1=Uniform(0-1)*[bhspinmag], 2=Belczynski2017)*/
void bse_set_bhspinmag(double bhspinmag);/* value of BH spins (default=0.0) */ 
void bse_set_mxns(double mxns); /* maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1) */
void bse_set_bconst(double bconst); /* isolated pulsar field decay timescale */
void bse_set_CK(double CK); /* Pulsar mass accretion field decay factor */
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
void bse_set_merger(double merger); /* pass through a value signifying a merger (>0.d0), evolv2.f will then do the appropriate kick and/or initial spin*/
void bse_set_id1_pass(long int id1_pass); /* pass through cmc star id into bse to help in debugging, this is used for iso star and star 1 in binary */
void bse_set_id2_pass(long int id2_pass); /* pass through cmc star id into bse to help in debugging, this is used for star 2 in binary */
void bse_set_taus113state(struct rng_t113_state state, int first);

/* getters */
double bse_get_alpha1(void); /* get CE alpha */
double bse_get_lambda(void); /* get CE lambda */
float bse_get_spp(int i, int j); /* stellar evolution log */
float bse_get_scm(int i, int j); /* stored stellar parameters at interval dtp */
float bse_get_bpp(int i, int j); /* binary evolution log */
float bse_get_bcm(int i, int j); /* stored binary parameters at interval dtp */
char *bse_get_sselabel(int kw); /* converts stellar type number to text label */
char *bse_get_bselabel(int kw); /* converts binary type number to text label */
struct rng_t113_state bse_get_taus113state(void);
int icase_get(int i, int j); /* use to get mixed type from ktype table */

/* copied functions */
double bse_kick_speed(int *startype); /* routine for generating birth kick speed from distribution */
double bse_rl(double q); /* Roche lobe formula */

/* useful macros */
#define BSE_WRAP_MAX(a, b) ((a)>=(b)?(a):(b))
#define BSE_WRAP_MIN(a, b) ((a)<=(b)?(a):(b))
