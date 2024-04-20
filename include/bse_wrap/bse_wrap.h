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
void evolv2_(int *kstar, double *mass, double *tb, double *ecc, double *z, 
	     double *tphysf, double *dtp, double *mass0, double *rad, double *lum,
             double *massc, double *radc, double *menv, double *renv,
	     double *ospin, double *B_0, double *bacc, double *tacc, double *epoch,
	     double *tms, double *bhspin, double *tphys, double *zpars, double *vs, double *kick_info);
void instar_(void);
float ran3_(int *idum);
void star_(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, 
	   double *lums, double *GB, double *zpars);
void hrdiag_(double *mass, double *aj, double *mt, double *tm, double *tn, double *tscls, 
	     double *lums, double *GB, double *zpars, double *r, double *lum, int *kw, 
	     double *mc, double *rc, double *menv, double *renv, double *k2,  double *bhspin, int *kidx);
void kick_(int *kw, double *m1, double *m1n, double *m2, double *ecc, double *sep, 
	   double *jorb, double *vk, int *snstar, double *r2, double *fallback, double *sigmahold, double *kick_info, int *disrupt, double *vs);
void mix_(double *mass, double *mt, double *aj, int *kw, double *zpars, double *bhspin);
// note: these function names only work if in lowercase here, even though FORTRAN versions in uppercase.
void comenv_(double *M01, double *M1, double *MC1, double *AJ1, double *JSPIN1, int *KW1, double *M02, double *M2, double *MC2, double *AJ2, double *JSPIN2, int *KW2, double *ZPARS, double *ECC, double *SEP, double *JORB, int *COEL, int *star1, int *star2, double *vk, double *kick_info, int *formation1, int *formation2, double *sigmahold, double *bhspin1, double *bhspin2, int *binstate, int *mergertype, int *jp, double *tphys,int *swtichedCE, double *rad, double *tms, double *evolve_type, int *disrupt, double * lumin, double * B_0, double * bacc, double * tacc, double * epoch, double * menv_bpp, double * renv_bpp, double *bkick);


/* wrapped BSE functions */
void bse_zcnsts(double *z, double *zpars);
void bse_evolv2(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		double *massc, double *radc, double *menv, double *renv, double *ospin,
                double *B_0, double *bacc, double *tacc,
		double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		double *z, double *zpars, double *tb, double *ecc, double *vs, double *bhspin);
void bse_evolve_single(int *kw, double *mass, double *mt, double *r, double *lum,
                double *mc, double *rc, double *menv, double *renv, double *ospin,
                double *epoch, double *tms, double *tphys, double *tphysf,
                double *dtp, double *z, double *zpars, double *vs, double *bhspin);
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
		double *mc, double *rc, double *menv, double *renv, double *k2, double *bhspin);
void bse_kick(int *kw, double *m1, double *m1n, double *m2, double *ecc, double *sep, 
	      double *jorb, double *vk, int *snstar, double *r2, double *fallback, double *vs);
void bse_mix(double *mass, double *mt, double *aj, int *kw, double *zpars, double *bhspin);
void bse_comenv(bse_binary *binary, double *zpars,
                double *vs, int *fb);

/* structs to access BSE common blocks */
/* note the index swap between fortran and C: i,j->j,i */
extern struct { int idum1; } rand1_;
extern struct { int idum2, iy, ir[32]; } rand2_;
extern struct { long long int state[4]; int first;} taus113state_;
extern struct { int ktype[15][15]; } types_;
extern struct { int  tflag, ifflag, remnantflag, wdflag, bhflag, windflag,  qcflag, eddlimflag, bhspinflag, aic, rejuvflag,  htpmb, st_cr, st_tide, bdecayfac, grflag, bhms_coll_flag, wd_mass_lim; } flags_;
extern struct { int ceflag,cekickflag,cemergeflag,cehestarflag,ussn; } ceflags_;
extern struct { int pisn_track[2]; } trackers_;
extern struct { double zsun; } metvars_;
extern struct { double don_lim, acc_lim; } mtvars_; 

extern struct { double neta, bwind, hewind, beta, xi, acc2, epsnov, eddfac, gamma; } windvars_;
extern struct { double qcrit_array[16], alpha1, lambdaf; } cevars_;
extern struct { double bconst, ck; } magvars_;
extern struct { double rejuv_fac; } mixvars_;
extern struct { double natal_kick_array[5][2], sigma, sigmadiv, bhsigmafrac, polar_kick_angle, pisn, ecsn, ecsn_mlow, bhspinmag, mxns, rembar_massloss; int kickflag;} snvars_;
extern struct { double fprimc_array[16]; } tidalvars_;
extern struct { double pts1, pts2, pts3; } points_;
extern struct { double dmmax, drmax; } tstepc_;
extern struct { double scm[14][50000], spp[3][20]; } single_;
extern struct { double bcm[38][50000], bpp[43][1000]; } binary_;
extern struct { double merger; long int id1_pass, id2_pass; long int using_cmc; } cmcpass_;

/* setters */
void bse_set_idum(int idum); /* RNG seed (for NS birth kicks) */

void bse_set_neta(double neta); /* Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally) */
void bse_set_bwind(double bwind); /* binary enhanced mass loss parameter (inactive for single) */
void bse_set_hewind(double hewind); /* helium star mass loss factor (1.0 normally) */
void bse_set_windflag(int windflag); /* Sets wind prescription (0=BSE, 1=StarTrack, 2=Vink; 0) */
void bse_set_eddlimflag(int eddlimflag); /* Sets wind prescription (0=BSE, 1=StarTrack, 2=Vink; 0) */
// rejuvflag toggles between the original BSE prescription for MS mixing and 
// lifetimes of stars based on the mass of the MS stars (equation 80) or a
// prescription that uses the ratio of helium core mass of the pre-merger stars
// at the base of the first ascent of the giant branch to determine relative to the
// helium core mass of the merger product at the base of the giant branch
void bse_set_rejuvflag(int rejuvflag);
void bse_set_grflag(int grflag);
void bse_set_bhms_coll_flag(int bhms_coll_flag);
void bse_set_wd_mass_lim(int wd_mass_lim);
void bse_set_kickflag(int kickflag);
void bse_set_using_cmc(void);
void bse_set_pisn(double pisn); /* Sets Pair-instability pulsations and supernoa */ 
void bse_set_aic(int aic); /* Sets Pair-instability pulsations and supernoa */
void bse_set_bdecayfac(int bdecayfac);
void bse_set_htpmb(int htpmb); /* Sets Pair-instability pulsations and supernoa */
void bse_set_st_tide(int st_tide); /* Sets Pair-instability pulsations and supernoa */
void bse_set_st_cr(int st_cr); /* Sets Pair-instability pulsations and supernoa */
void bse_set_ussn(int ussn); /* Sets Pair-instability pulsations and supernoa */
void bse_set_don_lim(double don_lim); /* Set donor limits for RLO mass loss */
void bse_set_acc_lim(double acc_lim); /* Set accretor limits for RLO mass accretion */

void bse_set_zsun(double zsun);

void bse_set_rembar_massloss(double rembar_massloss);
void bse_set_ecsn(double ecsn);
void bse_set_ecsn_mlow(double ecsn_mlow);
void bse_set_sigma(double sigma); /* dispersion in the Maxwellian for the SN kick speed (190 km/s) */
void bse_set_sigmadiv(double sigmadiv); /* dispersion in the Maxwellian for the SN kick speed (190 km/s) */
void bse_set_bhsigmafrac(double bhsigmafrac); /* Ad hoc factor to change BH SN kick speed relative to NS SN kick sigma (1) */
void bse_set_polar_kick_angle(int polar_kick_angle); /* Switch to set the allowed opening angle of SN kicks.  Defaults to 180 degrees*/
void bse_set_ifflag(int ifflag); /* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0) */
void bse_set_wdflag(int wdflag); /* wdflag > 0 uses modified-Mestel cooling for WDs (0) */
void bse_set_bhflag(int bhflag); /* bhflag > 0 allows velocity kick at BH formation (0) */
void bse_set_remnantflag(int remnantflag); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_qcflag(int qcflag); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_qcrit_array(double qcrit_array[16], long len); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_fprimc_array(double fprimc_array[16], long len); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_natal_kick_array(double natal_kick_array[10], long len); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
void bse_set_bhspinflag(int bhflag);/* bhspinflag (0=[bhspinmag], 1=Uniform(0-1)*[bhspinmag], 2=Belczynski2017)*/
void bse_set_bhspinmag(double bhspinmag);/* value of BH spins (default=0.0) */ 
void bse_set_mxns(double mxns); /* maximum NS mass (1.8, remnantflag=0; 3.0, remnantflag=1) */
void bse_set_bconst(double bconst); /* isolated pulsar field decay timescale */
void bse_set_CK(double CK); /* Pulsar mass accretion field decay factor */
void bse_set_rejuv_fac(double rejuv_fac);
void bse_set_pts1(double pts1); /* timestep taken in MS phase (0.05) */
void bse_set_pts2(double pts2); /* timestep taken in GB, CHeB, AGB, HeGB phases (0.01) */
void bse_set_pts3(double pts3); /* timestep taken in HG, HeMS phases (0.02) */
void bse_set_alpha1(double alpha1); /* common-envelope efficiency parameter (1.0) */
void bse_set_lambdaf(double lambdaf); /* binding energy factor for common envelope evolution (0.5) */
void bse_set_ceflag(int ceflag); /* ceflag > 0 activates spin-energy correction in common-envelope (0); ceflag = 3 activates de Kool common-envelope model (0) */
void bse_set_cehestarflag(int cehestarflag); /* cehestarflag > 0 activates spin-energy correction in common-envelope (0); cehestarflag = 3 activates de Kool common-envelope model (0) */
void bse_set_cemergeflag(int cemergeflag); /* cemergeflag > 0 activates spin-energy correction in common-envelope (0); cemergeflag = 3 activates de Kool common-envelope model (0) */
void bse_set_cekickflag(int cekickflag); /* cekickflag > 0 activates spin-energy correction in common-envelope (0); cekickflag = 3 activates de Kool common-envelope model (0) */
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
double bse_get_lambdaf(void); /* get CE lambda */
double bse_get_spp(int i, int j); /* stellar evolution log */
double bse_get_scm(int i, int j); /* stored stellar parameters at interval dtp */
double bse_get_bpp(int i, int j); /* binary evolution log */
double bse_get_bcm(int i, int j); /* stored binary parameters at interval dtp */
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
