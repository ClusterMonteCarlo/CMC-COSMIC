/* bse_wrap.c

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
#include "bse_wrap.h"

/**
* @brief calculate metallicity constants
*
* @param z ?
* @param zpars ?
*/
void bse_zcnsts(double *z, double *zpars)
{
  zcnsts_(z, zpars);
}

void bse_evolve_single(int *kw, double *mass, double *mt, double *r, double *lum,
		double *mc, double *rc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf,
		double *dtp, double *z, double *zpars, double *vs, double *bhspin) {
  bse_binary tempbinary;

  tempbinary.bse_mass0[0] = *mass;
  tempbinary.bse_mass0[1] = 0.0;
  tempbinary.bse_kw[0] = *kw;
  tempbinary.bse_kw[1] = 15;
  tempbinary.bse_mass[0] = *mt;
  tempbinary.bse_mass[1] = 0.0;
  tempbinary.bse_radius[0] = *r;
  tempbinary.bse_radius[1] = 0.0;
  tempbinary.bse_lum[0] = *lum;
  tempbinary.bse_lum[1] = 0.0;
  tempbinary.bse_massc[0] = *mc;
  tempbinary.bse_massc[1] = 0.0;
  tempbinary.bse_radc[0] = *rc;
  tempbinary.bse_radc[1] = 0.0;
  tempbinary.bse_menv[0] = *menv;
  tempbinary.bse_menv[1] = 0.0;
  tempbinary.bse_renv[0] = *renv;
  tempbinary.bse_renv[1] = 0.0;
  tempbinary.bse_ospin[0] = *ospin;
  tempbinary.bse_ospin[1] = 0.0;
  tempbinary.bse_bhspin[0] = *bhspin;
  tempbinary.bse_bhspin[1] = 0.0; 
  tempbinary.bse_B_0[0] = 0.0;
  tempbinary.bse_B_0[1] = 0.0;
  tempbinary.bse_bacc[0] = 0.0;
  tempbinary.bse_bacc[1] = 0.0;
  tempbinary.bse_tacc[0] = 0.0;
  tempbinary.bse_tacc[1] = 0.0;
  tempbinary.bse_epoch[0] = *epoch;
  tempbinary.bse_epoch[1] = 0.0;
  tempbinary.bse_tms[0] = *tms;
  tempbinary.bse_tms[1] = 0.0;
  tempbinary.bse_tb = 0.0;
  tempbinary.e = 0.0;

  bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]),
      &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]),
      &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]),
      &(tempbinary.bse_ospin[0]), &(tempbinary.bse_B_0[0]), &(tempbinary.bse_bacc[0]), &(tempbinary.bse_tacc[0]),
      &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]),
      tphys, tphysf, dtp, z, zpars, &(tempbinary.bse_tb), &(tempbinary.e), vs, &(tempbinary.bse_bhspin[0]));

  *mass = tempbinary.bse_mass0[0];
  *kw = tempbinary.bse_kw[0];
  *mt = tempbinary.bse_mass[0];
  *r = tempbinary.bse_radius[0];
  *lum = tempbinary.bse_lum[0];
  *mc = tempbinary.bse_massc[0];
  *rc = tempbinary.bse_radc[0];
  *menv = tempbinary.bse_menv[0];
  *renv = tempbinary.bse_renv[0];
  *ospin = tempbinary.bse_ospin[0];
  *epoch = tempbinary.bse_epoch[0];
  *tms = tempbinary.bse_tms[0];
  *bhspin = tempbinary.bse_bhspin[0];
}

/**
* @brief evolve a binary
*
* @param kstar ?
* @param mass0 ?
* @param mass ?
* @param rad ?
* @param lum ?
* @param massc ?
* @param radc ?
* @param menv ?
* @param renv ?
* @param ospin ?
* @param B_0 ?
* @param bacc ?
* @param tacc ?
* @param epoch ?
* @param tms ?
* @param tphys ?
* @param tphysf ?
* @param dtp ?
* @param z ?
* @param zpars ?
* @param tb ?
* @param ecc ?
* @param vs ?
*/
void bse_evolv2(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		double *massc, double *radc, double *menv, double *renv, double *ospin,
                double *B_0, double *bacc, double *tacc,
		double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		double *z, double *zpars, double *tb, double *ecc, double *vs, double *bhspin)
{
 int i;
  /* must null out vs, since SSE/BSE is not designed to return it and hence doesn't null it out */
/*  vs[0] = 0.0;
  vs[1] = 0.0;
  vs[2] = 0.0; */
  for(i=0;i<20;i++) {
      vs[i] = 0.0;
  }

  /* used by COSMIC, but not needed here */
   //double bppout[23][1000], bcmout[42][50000];

      evolv2_(kstar,mass,tb,ecc,z,tphysf,dtp,mass0,rad,lum,massc,radc, menv,renv,ospin,B_0,bacc,tacc,epoch,tms,bhspin,tphys,zpars,vs);

}

/**
* @brief evolve a binary star safely: in some cases, a merger can have non self-consistent properties, leading to crazy things like NaN radii---this is the easiest way to get around that problem
*
* @param kstar ?
* @param mass0 ?
* @param mass ?
* @param rad ?
* @param lum ?
* @param massc ?
* @param radc ?
* @param menv ?
* @param renv ?
* @param ospin ?
* @param B_0 ?
* @param bacc ?
* @param tacc ?
* @param epoch ?
* @param tms ?
* @param tphys ?
* @param tphysf ?
* @param dtp ?
* @param z ?
* @param zpars ?
* @param tb ?
* @param ecc ?
* @param vs ?
*/
void bse_evolv2_safely(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		       double *massc, double *radc, double *menv, double *renv, double *ospin,
                       double *B_0, double *bacc, double *tacc,
		       double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		       double *z, double *zpars, double *tb, double *ecc, double *vs, double *bhspin)
{
  int mykstar[2], mykstarprev[2], kattempt=-1, j, i;
  double mymass0[2], mymass[2], myrad[2], mylum[2], mymassc[2], myradc[2], mymenv[2], myrenv[2], myospin[2], myB_0[2], mybacc[2], mytacc[2], myepoch[2];
  double mytms[2], mytphys, mytphysf, mydtp, tphystried, mytb, myecc, myvs[20], mybhspin[2];

  //  do {
    kattempt++;
    for (j=0; j<2; j++) {
      mykstar[j] = kstar[j];
      //Sourav: testing additional controlable kick routine for blackholes
      //Save the initial startypes since formation kick must be added only when the SN happens
      //mykstarprev[j] = kstar[j];
      mymass0[j] = mass0[j];
      mymass[j] = mass[j];
      myrad[j] = rad[j];
      mylum[j] = lum[j];
      mymassc[j] = massc[j];
      myradc[j] = radc[j];
      mymenv[j] = menv[j];
      myrenv[j] = renv[j];
      myospin[j] = ospin[j];
      myB_0[j] = B_0[j]; /* PK */
      mybacc[j] = bacc[j];
      mytacc[j] = tacc[j];
      myepoch[j] = epoch[j];
      mytms[j] = tms[j];
      mybhspin[j] = bhspin[j];
    }
    mytphys = BSE_WRAP_MAX(*tphys - ((float) kattempt) * pow(1.2, kattempt), 0.0);
    tphystried = mytphys;
    mytphysf = *tphysf;
    /* try to interpret what is meant by the input dtp */
    if (*dtp == *tphysf - *tphys) {
      mydtp = mytphysf - mytphys;
    } else {
      mydtp = *dtp;
    }
    mytb = *tb;
    myecc = *ecc;
    bse_evolv2(mykstar, mymass0, mymass, myrad, mylum, mymassc, myradc, mymenv, myrenv, myospin, myB_0, mybacc, mytacc,
	       myepoch, mytms, &mytphys, &mytphysf, &mydtp, z, zpars, &mytb, &myecc, myvs, mybhspin);
    /*
  } while ((isnan(myrad[0]) || mymassc[0] < 0.0 || mymass[0] < 0.0 || mymass0[0] < 0.0 || mylum[0] < 0.0 || 
	    isnan(myrad[1]) || mymassc[1] < 0.0 || mymass[1] < 0.0 || mymass0[1] < 0.0 || mylum[1] < 0.0) && tphystried > 0.0);
  
  if (tphystried == 0.0) {
    fprintf(stderr, "bse_evolv2_safely(): Artifical age reduction failed.\n");
    exit(1);
  } else if (kattempt > 1) {
    fprintf(stderr, "bse_evolv2_safely(): Artifical age reduction succeeded.\n");
    fprintf(stderr, "bse_evolv2_safely(): kattempt=%d age_reduction=%g Myr.\n", kattempt, *tphys-tphystried);
  }
    */
  if(isnan(myrad[0])){
	fprintf(stderr, "An isnan occured for r1\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1]);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  if(isnan(myrad[1])){
	fprintf(stderr, "An isnan occured for r2\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g\n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1], *tb);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  if(mymass[0]<0.0){
	fprintf(stderr, "Mymass1<0.0\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1]);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  if(mymass[1]<0.0){
	fprintf(stderr, "Mymass2<0.0\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1]);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  if(mylum[0]<0.0){
	fprintf(stderr, "Mylum1<0.0\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1]);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  if(mylum[1]<0.0){
	fprintf(stderr, "Mylum2<0.0\n");
	fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", *tphys, *tphysf, kstar[0], kstar[1], mass[0], mass[1], rad[0], rad[1], lum[0], lum[1]);
	fprintf(stderr, "mytphys=%g mytphysf=%g mykw1=%d mykw2=%d mym1=%g mym2=%g myr1=%g myr2=%g myl1=%g myl2=%g \n", mytphys, mytphysf, mykstar[0], mykstar[1], mymass[0], mymass[1], myrad[0], myrad[1], mylum[0], mylum[1]);
	//exit(1);
  }
  
  for (j=0; j<2; j++) {
      kstar[j] = mykstar[j];
      mass0[j] = mymass0[j];
      mass[j] = mymass[j];
      rad[j] = myrad[j];
      lum[j] = mylum[j];
      massc[j] = mymassc[j];
      radc[j] = myradc[j];
      menv[j] = mymenv[j];
      renv[j] = myrenv[j];
      ospin[j] = myospin[j];
      B_0[j] = myB_0[j];
      bacc[j] = mybacc[j];
      tacc[j] = mytacc[j];
      epoch[j] = myepoch[j];
      tms[j] = mytms[j];
	  bhspin[j] = mybhspin[j];
    }
  *tphys = mytphys;
  *tphysf = mytphysf;
  *dtp = mydtp;
  *tb = mytb;
  *ecc = myecc;
/*  vs[0] = myvs[0];
  vs[1] = myvs[1];
  vs[2] = myvs[2];*/
  for(i=0;i<20;i++) {
      vs[i] = myvs[i];
  }
}


/**
* @brief set collision matrix
*/
void bse_instar(void)
{
  instar_();
}

/**
* @brief star routine; shouldn't need to be used often outside of BSE
*
* @param kw stellar type
* @param mass initial mass
* @param mt current mass
* @param tm main sequence timescale
* @param tn nuclear timescale
* @param tscls timescales array (20)
* @param lums luminosities array (10)
* @param GB giant branch parameters (10)
* @param zpars evolution parameters based on metallicity (20)
*/
void bse_star(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, 
	      double *lums, double *GB, double *zpars)
{
  /* INPUTS

     kw = stellar type
     mass = initial mass
     mt = current mass
     zpars = evolution parameters based on metallicity (20)

     OUTPUTS
     
     tm = main sequence timescale
     tn = nuclear timescale
     tscls = timescales array (20)
     lums = luminosities array (10)
     GB = giant branch parameters (10)
   */

  star_(kw, mass, mt, tm, tn, tscls, lums, GB, zpars);
}

/* hrdiag routine; shouldn't need to be used often outside of BSE */
void bse_hrdiag(double *mass, double *aj, double *mt, double *tm, double *tn, double *tscls, 
		double *lums, double *GB, double *zpars, double *r, double *lum, int *kw, 
		double *mc, double *rc, double *menv, double *renv, double *k2, double *bhspin)
{
  /* INPUTS
     
     mass = mass (old value)
     aj = current age
     mt = current mass
     tm = main sequence timescale
     tn = nuclear timescale
     tscls = timescales array (20)
     lums = luminosities array (10)
     GB = giant branch parameters (10)
     zpars = evolution parameters based on metallicity (20)
     kw = stellar type
     
     OUTPUTS
     
     mass = mass (new value)
     aj = curret age (new value)
     mt = current mass (new value)
     r = radius
     lum = luminosity
     kw = stellar type (new value)
     mc = core mass
     rc = core radius
     menv = envelope mass
     renv = envelope radius
     k2 = radius of gyration of envelope
   */
  int kidx=0;

  hrdiag_(mass, aj, mt, tm, tn, tscls, lums, GB, zpars, r, lum, kw, mc, rc, menv, renv, k2, bhspin, &kidx);
}

/**
* @brief kick routine; shouldn't need to be used often outside of BSE
*
* @param kw stellar type
* @param m1 mass of star 1
* @param m1n new mass of star 1
* @param m2 mass of star 2
* @param ecc orbit eccentrity
* @param sep orbit semimajor axis
* @param jorb orbital angular momentum
* @param vk kick magnitude, can be used to help set initial pulsar properties if desired
* @param snstar which star (primary 1 or secondary 2) that goes SN
* @param r2 radius of companion, stars may collide if lucky kick
* @param fallback fraction of pre-SN mass that may fall back onto remnant. Can be 0. Kick strength is limited by such fall back of material
* @param vs
     old -> vs = velocity (3) of center of mass if bound, relative velocity at infinity if unbound
     new -> vs = three possible sets of velocities (3). Is an array of 12, correctly accounts for
*/
void bse_kick(int *kw, double *m1, double *m1n, double *m2, double *ecc, double *sep, 
	      double *jorb, double *vk, int *snstar, double *r2, double *fallback, double *vs)
{
  /* INPUTS
     kw = stellar type
     m1 = mass of star 1
     m1n = new mass of star 1
     m2 = mass of star 2
     ecc = orbit eccentrity
     sep = orbit semimajor axis
     snstar = which star (primary 1 or secondary 2) that goes SN
     r2 = radius of companion, stars may collide if lucky kick
     fallback = fraction of pre-SN mass that may fall back onto remnant. Can be 0.
                Kick strength is limited by such fall back of material
     
     OUTPUTS
     kw = stellar type (new value if input kw<0)
     ecc = eccentricity (new value)
     sep = orbit semimajor axis (new value)
     jorb = orbital angular momentum
     vk = kick magnitude, can be used to help set initial pulsar properties if desired.
     old -> vs = velocity (3) of center of mass if bound, relative velocity at infinity if unbound
     new -> vs = three possible sets of velocities (3). Is an array of 12, correctly accounts for
     run away velocities and outputs them in the original COM frame. 
     Array of 12 for two kicks that can occur - one kick may disrupt the system so each star has 
     output into one of the 3 velocity slots.
     Another 3 values within the array show which star (1 or 2) went SN for that kick.
     This helps in differentiating which kick goes where.
   */
  /* LOGICAL used by COSMIC, but not needed here */
  int disrupt=0;
  kick_(kw, m1, m1n, m2, ecc, sep, jorb, vk, snstar, r2, fallback, vs, &disrupt);
}

/**
* @brief Mix routine, call from cmc_sscollision.c within the merge_two_stars routine. Helps merge troublesome systems that also wont pass through a common envelope.
*
* @param mass ?
* @param mt ?
* @param aj ?
* @param kw ?
* @param zpars ?
* @param ecsnp ?
*/
void bse_mix(double *mass, double *mt, double *aj, int *kw, double *zpars, double *bhspin)
{
  mix_(mass, mt, aj, kw, zpars, bhspin);
}

/**
* @brief ?
*
* @param tempbinary ?
* @param zpars ?
* @param vs ?
* @param fb ?
* @param ecsnp ?
* @param ecsn_mlow ?
*/
void bse_comenv(bse_binary *tempbinary, double *zpars, double *vs, int *fb)
		//double *M0, double *M, double *MC, double *AJ, double *OSPIN, int *KW, 
		//                double *M02, double *M2, double *MC2, double *AJ2, double *JSPIN2, int *KW2,
                //double *ZPARS, double *ECC, double *SEP, double *PORB,  
{
  int i, kcomp1,kcomp2, star1, star2, COEL;
  double vk,OORB,JORB,mce[2],AJ[2], M0ce[2], PI, JSPIN1, JSPIN2;
  double tm, tn, tscls[20], lums[10], GB[10], k2, k3;
  double bhspin[2];
  int binstate=0, mergertype=0;
  int jp=0,switchedCE=0, disrupt=0;
  double tphys=0, evolve_type=0;
  double tms[2], rad[2];
  k3 = 0.21;
  PI = acos(-1.0);
  //
  //  INPUTS
  //  Initial masses [Solar], current masses [Solar], current core masses [Solar], ages, spin vel., stellar types,
  //  zpars (interesting values), eccentricity, separation [Rsun], orbital period [years, convert to days], star values help in bse for memory of star array order,
  //  fb tells code if we want fallback or not during SN, ecsnp and ecsn_mlow provide electron capture kick option and are critical mass limits,
  //  ST_tide tells code if StarTrack-like tides are assumed in code.
  //
  //  OUTPUTS
  //  Final masses etc, ecc, sep, Porb will be modified in CE evolution, vs is an array that might be filled if a NS is formed, formations tell
  //  the type of SN that occured (if one occured).
  //
  //
  // Note that because we are making tempbinary particulars pointers we can't 
  // simple dot to the structure values but instead need to point to them with either b->x or (*b).x
  //
  (*tempbinary).bse_tb = (*tempbinary).bse_tb/365.25; //convert to years
  //  tempbinary.a = tempbinary.a * cmc_l_unit/RSUN; //separation in solar radii as used in comenv.f
  AJ[0] = (*tempbinary).bse_tphys - (*tempbinary).bse_epoch[0];
  AJ[1] = (*tempbinary).bse_tphys - (*tempbinary).bse_epoch[1];
  star1 = 1;
  star2 = 2;
  COEL = 0;
  mce[0] = (*tempbinary).bse_mass[0];
  mce[1] = (*tempbinary).bse_mass[1];
  M0ce[0] = (*tempbinary).bse_mass0[0];
  M0ce[1] = (*tempbinary).bse_mass0[1];
  kcomp1 = (*tempbinary).bse_kw[0];
  kcomp2 = (*tempbinary).bse_kw[1];
  OORB = (2.0*PI)/((*tempbinary).bse_tb); //units: 1/year
  JORB = (((*tempbinary).bse_mass[0])*((*tempbinary).bse_mass[1]))/(((*tempbinary).bse_mass[0])+((*tempbinary).bse_mass[1]))*sqrt(1.0-((*tempbinary).e)*((*tempbinary).e))*((*tempbinary).a)*((*tempbinary).a)*OORB; //solar mass * Rsun^2 / year
  //
  // Call hrdiag.f to receive stellar spin structure (k2); requires output from star.f.
  bse_star(&((*tempbinary).bse_kw[0]), &((*tempbinary).bse_mass0[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars);
  //
  bse_hrdiag(&((*tempbinary).bse_mass0[0]), &(AJ[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[0]), &((*tempbinary).bse_lum[0]), &((*tempbinary).bse_kw[0]), &((*tempbinary).bse_massc[0]), &((*tempbinary).bse_radc[0]), 
	     &((*tempbinary).bse_menv[0]), &((*tempbinary).bse_renv[0]), &k2, &((*tempbinary).bse_bhspin[0]));
  ////
  JSPIN1 = (*tempbinary).bse_ospin[0]*(k2*(*tempbinary).bse_radius[0]*(*tempbinary).bse_radius[0]*(*tempbinary).bse_menv[0] + 
				       k3*(*tempbinary).bse_radc[0]*(*tempbinary).bse_radc[0]*(*tempbinary).bse_massc[0]);
  ////
  ////
  bse_star(&((*tempbinary).bse_kw[1]), &((*tempbinary).bse_mass0[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars);
  //
  bse_hrdiag(&((*tempbinary).bse_mass0[1]), &(AJ[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[1]), &((*tempbinary).bse_lum[1]), &((*tempbinary).bse_kw[1]), &((*tempbinary).bse_massc[1]), &((*tempbinary).bse_radc[1]), 
	     &((*tempbinary).bse_menv[1]), &((*tempbinary).bse_renv[1]), &k2, &((*tempbinary).bse_bhspin[1]));
  ////
  JSPIN2 = (*tempbinary).bse_ospin[1]*(k2*(*tempbinary).bse_radius[1]*(*tempbinary).bse_radius[1]*(*tempbinary).bse_menv[1] + 
				       k3*(*tempbinary).bse_radc[1]*(*tempbinary).bse_radc[1]*(*tempbinary).bse_massc[1]);
  //bse_hrdiag(&(M0ce), &(AJ), &M, double *tm, double *tn, double *tscls, 
  //     double *lums, double *GB, double *zpars, double *r, double *lum, int *kw, 
  //     double *mc, double *rc, double *menv, double *renv, double *k2, int *ST_tide, double *ecsnp, double *ecsn_mlow)
  //
  vk = 0.0;
  for(i=0;i<20;i++) {
      vs[i] = 0.0;
  }
  ////
  //// Why we are here:
  ////
  //  COMENV_((tempbinary.bse_mass0[0]), (tempbinary.bse_mass[0]), (tempbinary.bse_massc[0]), &(AJ[0]), &JSPIN1, (tempbinary.bse_kw[0]),
  //	  (tempbinary.bse_mass0[1]), (tempbinary.bse_mass[1]), (tempbinary.bse_massc[1]), &(AJ[1]), &JSPIN2, (tempbinary.bse_kw[1]),
  //	  zpars, (tempbinary.e), (tempbinary.a), &(JORB), &COEL, &star1, &star2, &vk, fb, vs, ecsnp, ecsn_mlow, &(tempbinary.bse_bcm_formation[0]), &(tempbinary.bse_bcm_formation[1]), ST_tide);
  //printf("bse_wrap ZPARS: \n");
  //for(iii=0;iii<20; iii++){
  //  printf("%g ", zpars[iii]);
  //}
  //printf(" kw1i=%d kw2i=%d m1i=%g m2i=%g r1i=%g r2i=%g epoch1=%g epoch2=%g ", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1], (*tempbinary).bse_epoch[0], (*tempbinary).bse_epoch[1]);
  comenv_(&((*tempbinary).bse_mass0[0]), &((*tempbinary).bse_mass[0]), &((*tempbinary).bse_massc[0]), &(AJ[0]), &JSPIN1, &((*tempbinary).bse_kw[0]),
	  &((*tempbinary).bse_mass0[1]), &((*tempbinary).bse_mass[1]), &((*tempbinary).bse_massc[1]), &(AJ[1]), &JSPIN2, &((*tempbinary).bse_kw[1]),
	  zpars, &((*tempbinary).e), &((*tempbinary).a), &(JORB), &COEL, &star1, &star2, &vk, vs, &((*tempbinary).bse_bcm_formation[0]), &((*tempbinary).bse_bcm_formation[1]), &((*tempbinary).bse_bhspin[0]),&((*tempbinary).bse_bhspin[1]),&binstate,&mergertype,&jp,&tphys,&switchedCE,rad,tms,&evolve_type,&disrupt);
  //printf("kw1i=%d kw2i=%d m1f=%g m2f=%g r1f=%g r2f=%g ", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1]);
  //printf("\n");
  ////
  //// Update things to do with system components...
  ////
  (*tempbinary).bse_epoch[0] = (*tempbinary).bse_tphys - AJ[0];
  (*tempbinary).bse_epoch[1] = (*tempbinary).bse_tphys - AJ[1];
  if(kcomp1==13 && (*tempbinary).bse_kw[0]==15 && (*tempbinary).bse_kw[1]==13){
    (*tempbinary).bse_bcm_formation[1] = (*tempbinary).bse_bcm_formation[0];
  }
  if(kcomp2==13 && (*tempbinary).bse_kw[1]==15 && (*tempbinary).bse_kw[0]==13){
    (*tempbinary).bse_bcm_formation[0] = (*tempbinary).bse_bcm_formation[1];
  }
  if((*tempbinary).a>0.0 && (*tempbinary).e<1.0 && (*tempbinary).e>=0.0 && (*tempbinary).bse_mass[0] > 0.0 && (*tempbinary).bse_mass[1] > 0.0){
    bse_star(&((*tempbinary).bse_kw[0]), &((*tempbinary).bse_mass0[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars);
    //
    bse_hrdiag(&((*tempbinary).bse_mass0[0]), &(AJ[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[0]), &((*tempbinary).bse_lum[0]), &((*tempbinary).bse_kw[0]), &((*tempbinary).bse_massc[0]), &((*tempbinary).bse_radc[0]), 
	     &((*tempbinary).bse_menv[0]), &((*tempbinary).bse_renv[0]), &k2, &((*tempbinary).bse_bhspin[0]));
    bse_star(&((*tempbinary).bse_kw[1]), &((*tempbinary).bse_mass0[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars);
    //
    bse_hrdiag(&((*tempbinary).bse_mass0[1]), &(AJ[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[1]), &((*tempbinary).bse_lum[1]), &((*tempbinary).bse_kw[1]), &((*tempbinary).bse_massc[1]), &((*tempbinary).bse_radc[1]), 
	     &((*tempbinary).bse_menv[1]), &((*tempbinary).bse_renv[1]), &k2, &((*tempbinary).bse_bhspin[1]));
    (*tempbinary).bse_epoch[0] = (*tempbinary).bse_tphys - AJ[0];
    (*tempbinary).bse_epoch[1] = (*tempbinary).bse_tphys - AJ[1];
    OORB = JORB/(((*tempbinary).bse_mass[0]*(*tempbinary).bse_mass[1])/((*tempbinary).bse_mass[0]+(*tempbinary).bse_mass[1])*sqrt(1.0-((*tempbinary).e*(*tempbinary).e))*(*tempbinary).a*(*tempbinary).a);///214.9/214.9); //turn sep in Rsun into sep in AU by dividing by ~215 to find orb velocity freq in 1/years
    (*tempbinary).bse_tb = (2.0*PI*365.25)/OORB; //convert from years to days.
    //    tempbinary.a = tempbinary.a * RSUN/cmc_l_unit;
    printf("after hrdiag: kw1i=%d kw2i=%d m1f=%g m2f=%g r1f=%g r2f=%g epoch1=%g epoch2=%g ", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1], (*tempbinary).bse_epoch[0], (*tempbinary).bse_epoch[1]);
    printf("\n");
  } else { //either disrupted or merged. If disrupted hopefully we will be able to merge it later (even if we must use mix on the post CE stars!).
    (*tempbinary).a = 0.0;
    (*tempbinary).e = -1.0;
    (*tempbinary).bse_tb = 0.0;
    if((*tempbinary).bse_mass[0]==0.0){
      bse_star(&((*tempbinary).bse_kw[1]), &((*tempbinary).bse_mass0[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars);
      //
      bse_hrdiag(&((*tempbinary).bse_mass0[1]), &(AJ[1]), &((*tempbinary).bse_mass[1]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[1]), &((*tempbinary).bse_lum[1]), &((*tempbinary).bse_kw[1]), &((*tempbinary).bse_massc[1]), &((*tempbinary).bse_radc[1]), 
	     &((*tempbinary).bse_menv[1]), &((*tempbinary).bse_renv[1]), &k2, &((*tempbinary).bse_bhspin[1]));
      (*tempbinary).bse_epoch[1] = (*tempbinary).bse_tphys - AJ[1];
      printf("after hrdiag: kw1i=%d kw2i=%d m1f=%g m2f=%g r1f=%g r2f=%g epoch1=%g epoch2=%g ", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1], (*tempbinary).bse_epoch[0], (*tempbinary).bse_epoch[1]);
      printf("\n");
    } else {
      bse_star(&((*tempbinary).bse_kw[0]), &((*tempbinary).bse_mass0[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars);
      //
      bse_hrdiag(&((*tempbinary).bse_mass0[0]), &(AJ[0]), &((*tempbinary).bse_mass[0]), &tm, &tn, tscls, lums, GB, zpars,
	     &((*tempbinary).bse_radius[0]), &((*tempbinary).bse_lum[0]), &((*tempbinary).bse_kw[0]), &((*tempbinary).bse_massc[0]), &((*tempbinary).bse_radc[0]), 
	     &((*tempbinary).bse_menv[0]), &((*tempbinary).bse_renv[0]), &k2, &((*tempbinary).bse_bhspin[0]));
      (*tempbinary).bse_epoch[0] = (*tempbinary).bse_tphys - AJ[0];
      printf("after hrdiag: kw1i=%d kw2i=%d m1f=%g m2f=%g r1f=%g r2f=%g epoch1=%g epoch2=%g ", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1], (*tempbinary).bse_epoch[0], (*tempbinary).bse_epoch[1]);
      printf("\n");
    }
  }

}

/* setters */
//OPT: Use inline before void -> makes it macro
void bse_set_using_cmc(void) {cmcpass_.using_cmc = 1; }
void bse_set_idum(int idum) { rand1_.idum1 = idum; }
void bse_set_eddlimflag(int eddlimflag) { flags_.eddlimflag = eddlimflag; } 
void bse_set_sigmadiv(double sigmadiv) { snvars_.sigmadiv = sigmadiv; } 

void bse_set_natal_kick_array(double natal_kick_array[6], long len_kick) {int i; for(i=0; i<len_kick; i++) snvars_.natal_kick_array[i] = natal_kick_array[i]; }
void bse_set_fprimc_array(double fprimc_array[16], long len_fprimc) {int i; for(i=0;i<len_fprimc;i++) tidalvars_.fprimc_array[i] = fprimc_array[i]; }
void bse_set_qcrit_array(double qcrit_array[16], long len_qcrit) {int i; for(i=0;i<len_qcrit;i++) cevars_.qcrit_array[i] = qcrit_array[i]; }

void bse_set_aic(int aic) { flags_.aic = aic; }
void bse_set_bdecayfac(int bdecayfac) { flags_.bdecayfac = bdecayfac; }
void bse_set_st_cr(int st_cr) { flags_.st_cr = st_cr; }
void bse_set_st_tide(int st_tide) { flags_.st_tide = st_tide; }
void bse_set_htpmb(int htpmb) { flags_.htpmb = htpmb; }
void bse_set_rejuvflag(int rejuvflag) { flags_.rejuvflag = rejuvflag; }
void bse_set_qcflag(int qcflag) { flags_.qcflag = qcflag; }
void bse_set_ussn(int ussn) { ceflags_.ussn = ussn; }
void bse_set_neta(double neta) { windvars_.neta = neta; }
void bse_set_bwind(double bwind) { windvars_.bwind = bwind; }
void bse_set_hewind(double hewind) { windvars_.hewind = hewind; }
void bse_set_windflag(int windflag) { flags_.windflag = windflag; }
void bse_set_pisn(double pisn) { snvars_.pisn = pisn; }
void bse_set_ecsn(double ecsn) { snvars_.ecsn = ecsn; }
void bse_set_ecsn_mlow(double ecsn_mlow) { snvars_.ecsn_mlow = ecsn_mlow; }
void bse_set_sigma(double sigma) { snvars_.sigma = sigma; }
void bse_set_bhsigmafrac(double bhsigmafrac) {snvars_.bhsigmafrac = bhsigmafrac; }
void bse_set_polar_kick_angle(int polar_kick_angle) {snvars_.polar_kick_angle = polar_kick_angle; }
void bse_set_ifflag(int ifflag) { flags_.ifflag = ifflag; }
void bse_set_wdflag(int wdflag) { flags_.wdflag = wdflag; }
void bse_set_bhflag(int bhflag) { flags_.bhflag = bhflag; }
void bse_set_bhspinflag(int bhspinflag) { flags_.bhspinflag = bhspinflag; }
void bse_set_bhspinmag(double bhspinmag) { snvars_.bhspinmag = bhspinmag; }
void bse_set_nsflag(int nsflag) { flags_.nsflag = nsflag; }
void bse_set_mxns(double mxns) { windvars_.mxns = mxns;} 
void bse_set_bconst(double bconst) { magvars_.bconst = bconst; }
void bse_set_CK(double CK) {magvars_.ck = CK;}
void bse_set_rejuv_fac(double rejuv_fac) {mixvars_.rejuv_fac = rejuv_fac;}
void bse_set_pts1(double pts1) { points_.pts1 = pts1; }
void bse_set_pts2(double pts2) { points_.pts2 = pts2; }
void bse_set_pts3(double pts3) { points_.pts3 = pts3; }
void bse_set_alpha1(double alpha1) { cevars_.alpha1 = alpha1; }
void bse_set_lambdaf(double lambdaf) { cevars_.lambdaf = lambdaf; }
void bse_set_ceflag(int ceflag) { ceflags_.ceflag = ceflag; }
void bse_set_cemergeflag(int cemergeflag) { ceflags_.cemergeflag = cemergeflag; }
void bse_set_cekickflag(int cekickflag) { ceflags_.cekickflag = cekickflag; }
void bse_set_cehestarflag(int cehestarflag) { ceflags_.cehestarflag = cehestarflag; }
void bse_set_tflag(int tflag) { flags_.tflag = tflag; }
void bse_set_beta(double beta) { windvars_.beta = beta; }
void bse_set_xi(double xi) { windvars_.xi = xi; }
void bse_set_acc2(double acc2) { windvars_.acc2 = acc2; }
void bse_set_epsnov(double epsnov) { windvars_.epsnov = epsnov; }
void bse_set_eddfac(double eddfac) { windvars_.eddfac = eddfac; }
void bse_set_gamma(double gamma) { windvars_.gamma = gamma; }
void bse_set_merger(double merger) {cmcpass_.merger = merger; }
void bse_set_id1_pass(long int id1_pass) { cmcpass_.id1_pass = id1_pass; }
void bse_set_id2_pass(long int id2_pass) { cmcpass_.id2_pass = id2_pass; }
#ifndef USE_TAUS
/* need to define this here, as it is not defined in the Fortran part
 * when we don't use the Tausworthe generator for BSE kicks */
static struct { long long int state[4]; int first;} taus113state_;
#endif

/**
* @brief copies the C tausworthe rng state variables to the Fortran states
*
* @param state C rng state
* @param first ?
*/
void bse_set_taus113state(struct rng_t113_state state, int first) {

  taus113state_.state[0]= state.z[0];
  taus113state_.state[1]= state.z[1];
  taus113state_.state[2]= state.z[2];
  taus113state_.state[3]= state.z[3];

  if (first>=0) {
    taus113state_.first= first;
  }
}

/* getters */
/* note the index flip and decrement so the matrices are accessed
   as they would be in fortran */
double bse_get_alpha1(void) { return(cevars_.alpha1); }
double bse_get_lambdaf(void) { return(cevars_.lambdaf); }
double bse_get_spp(int i, int j) { return(single_.spp[j-1][i-1]); }
double bse_get_scm(int i, int j) { return(single_.scm[j-1][i-1]); }
double bse_get_bpp(int i, int j) { return(binary_.bpp[j-1][i-1]); }
double bse_get_bcm(int i, int j) { return(binary_.bcm[j-1][i-1]); }

/**
* @brief copies back the Fortran tausworthe rng state variables to the C state
*
* @return state with variables copied from the Fortran rng
*/
struct rng_t113_state bse_get_taus113state(void) {
  struct rng_t113_state state;

  state.z[0]= taus113state_.state[0];
  state.z[1]= taus113state_.state[1];
  state.z[2]= taus113state_.state[2];
  state.z[3]= taus113state_.state[3];

  return (state);
}

/**
* @brief ?
*
* @param i ?
* @param j ?
*
* @return ?
*/
int icase_get(int i, int j) { return(types_.ktype[j][i]); }

/**
* @brief ?
*
* @param kw ?
*
* @return ?
*/
char *bse_get_sselabel(int kw)
{
  if (kw == 0) {
    return("Low Mass MS Star");
  } else if (kw == 1) {
    return("Main sequence Star");
  } else if (kw == 2) {
    return("Hertzsprung Gap");
  } else if (kw == 3) {
    return("Giant Branch");
  } else if (kw == 4) {
    return("Core Helium Burning");
  } else if (kw == 5) {
    return("First AGB");
  } else if (kw == 6) {
    return("Second AGB");
  } else if (kw == 7) {
    return("Naked Helium MS");
  } else if (kw == 8) {
    return("Naked Helium HG");
  } else if (kw == 9) {
    return("Naked Helium GB");
  } else if (kw == 10) {
    return("Helium WD");
  } else if (kw == 11) {
    return("Carbon/Oxygen WD");
  } else if (kw == 12) {
    return("Oxygen/Neon WD");
  } else if (kw == 13) {
    return("Neutron Star");
  } else if (kw == 14) {
    return("Black Hole");
  } else if (kw == 15) {
    return("Massless Supernova");
  } else {
    return("Unkown Stellar Type!");
  }
}

/**
* @brief ?
*
* @param kw ?
*
* @return ?
*/
char *bse_get_bselabel(int kw)
{
  if (kw == 1) {
    return("INITIAL");
  } else if (kw == 2) {
    return("KW CHNGE");
  } else if (kw == 3) {
    return("BEG RCHE");
  } else if (kw == 4) {
    return("END RCHE");
  } else if (kw == 5) {
    return("CONTACT");
  } else if (kw == 6) {
    return("COELESCE");
  } else if (kw == 7) {
    return("COMENV");
  } else if (kw == 8) {
    return("GNTAGE");
  } else if (kw == 9) {
    return("NO REMNT");
  } else if (kw == 10) {
    return("MAX TIME");
  } else if (kw == 11) {
    return("DISRUPT");
  } else if (kw == 12) {
    return("BEG SYMB");
  } else if (kw == 13) {
    return("END SYMB");
  } else if (kw == 14) {
    return("BEG BSS");
  } else {
    return("Unkown Stellar Type!");
  }
}

/**
* @brief kick speed from distribution, taken directly from BSE code.
* @details bse_evolv1() and bse_evolv2() return the kick speed, but the kick function is included here just in case you want to use it. may change startype in certain cases
*
* @param startype star type
*
* @return ?
*/
double bse_kick_speed(int *startype)
{
  int k;
  double u1, u2, s, theta, vk2, vk, v[4]; /* yes, v is supposed to be 4-D */

  for (k=1; k<=2; k++) {
    u1 = ran3_(&(rand1_.idum1));
    u2 = ran3_(&(rand1_.idum1));
    s = snvars_.sigma * sqrt(-2.0*log(1.0-u1));
    theta = 2.0 * M_PI * u2;
    v[2*k-1-1] = s*cos(theta);
    v[2*k-1] = s*sin(theta);
  }
  vk2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  vk = sqrt(vk2);

  if (((*startype) == 14 && flags_.bhflag == 0) || (*startype) < 0) {
    vk2 = 0.0;
    vk = 0.0;
    if ((*startype) < 0) {
      (*startype) = 13;
    }
  }

  return(vk);
}

/**
* @brief Roche lobe formula, taken directly from BSE code.
*
* @param q ?
*
* @return Returns R_L1/a, where q=m_1/m_2.
*/
double bse_rl(double q)
{
  double p;

  p = pow(q, 1.0/3.0);

  return(0.49*p*p/(0.6*p*p+log(1.0+p)));
}
