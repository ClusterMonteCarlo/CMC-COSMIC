#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bse_wrap.h"

/* calculate metallicity constants */
void bse_zcnsts(double *z, double *zpars)
{
  zcnsts_(z, zpars);
}

/* evolve a single star */
void bse_evolv1(int *kw, double *mass, double *mt, double *r, double *lum,
		double *mc, double *rc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf, 
		double *dtp, double *z, double *zpars)
{
  evolv1_(kw, mass, mt, r, lum, mc, rc, menv, renv, ospin,
	  epoch, tms, tphys, tphysf, dtp, z, zpars);
}

/* evolve a binary */
void bse_evolv2(int *kstar, double *mass0, double *mass, double *rad, double *lum, 
		double *massc, double *radc, double *menv, double *renv, double *ospin,
		double *epoch, double *tms, double *tphys, double *tphysf, double *dtp,
		double *z, double *zpars, double *tb, double *ecc)
{
  evolv2_(kstar, mass0, mass, rad, lum, massc, radc, menv, renv, ospin,
	  epoch, tms, tphys, tphysf, dtp, z, zpars, tb, ecc);
}

/* set collision matrix */
void bse_instar(void)
{
  instar_();
}

/* setters */
void bse_set_idum(int idum) { value3_.idum = idum; }
void bse_set_neta(double neta) { value1_.neta = neta; }
void bse_set_bwind(double bwind) { value1_.bwind = bwind; }
void bse_set_hewind(double hewind) { value1_.hewind = hewind; }
void bse_set_sigma(double sigma) { value4_.sigma = sigma; }
void bse_set_ifflag(int ifflag) { flags_.ifflag = ifflag; }
void bse_set_wdflag(int wdflag) { flags_.wdflag = wdflag; }
void bse_set_bhflag(int bhflag) { value4_.bhflag = bhflag; }
void bse_set_nsflag(int nsflag) { flags_.nsflag = nsflag; }
void bse_set_mxns(double mxns) { value1_.mxns = mxns; }
void bse_set_pts1(double pts1) { points_.pts1 = pts1; }
void bse_set_pts2(double pts2) { points_.pts2 = pts2; }
void bse_set_pts3(double pts3) { points_.pts3 = pts3; }
void bse_set_alpha1(double alpha1) { value2_.alpha1 = alpha1; }
void bse_set_lambda(double lambda) { value2_.lambda = lambda; }
void bse_set_ceflag(int ceflag) { flags_.ceflag = ceflag; }
void bse_set_tflag(int tflag) { flags_.tflag = tflag; }
void bse_set_beta(double beta) { value5_.beta = beta; }
void bse_set_xi(double xi) { value5_.xi = xi; }
void bse_set_acc2(double acc2) { value5_.acc2 = acc2; }
void bse_set_epsnov(double epsnov) { value5_.epsnov = epsnov; }
void bse_set_eddfac(double eddfac) { value5_.eddfac = eddfac; }
void bse_set_gamma(double gamma) { value5_.gamma = gamma; }

/* getters */
/* note the index flip and decrement so the matrices are accessed
   as they would be in fortran */
float bse_get_spp(int i, int j) { return(single_.spp[j-1][i-1]); }
float bse_get_scm(int i, int j) { return(single_.scm[j-1][i-1]); }
float bse_get_bpp(int i, int j) { return(binary_.bpp[j-1][i-1]); }
float bse_get_bcm(int i, int j) { return(binary_.bcm[j-1][i-1]); }

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
