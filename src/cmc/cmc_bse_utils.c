/* vi: set filetype=c.doxygen: */

#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap.h" 

/**
* @brief ?
*
* @param bvars ?
* @param svars ?
* @param bmember ?
*/
void update_bse_from_sse(bse_t *bvars, sse_t *svars, int bmember) {
  bvars->bse_mass0[bmember] = svars->se_mass;
  bvars->bse_kw[bmember] = svars->se_k;
  bvars->bse_mass[bmember] = svars->se_mt;
  bvars->bse_bhspin[bmember] = svars->se_bhspin;
  bvars->bse_ospin[bmember] = svars->se_ospin;
  bvars->bse_B_0[bmember] = svars->se_B_0;
  bvars->bse_bacc[bmember] = svars->se_bacc;
  bvars->bse_tacc[bmember] = svars->se_tacc;
  bvars->bse_epoch[bmember] = svars->se_epoch;
  bvars->bse_tphys = svars->se_tphys; /* tphys should be the same for both input stars so this should be OK */
  bvars->bse_radius[bmember] = svars->se_radius;
  bvars->bse_lum[bmember] = svars->se_lum;
  bvars->bse_massc[bmember] = svars->se_mc;
  bvars->bse_radc[bmember] = svars->se_rc;
  bvars->bse_menv[bmember] = svars->se_menv;
  bvars->bse_renv[bmember] = svars->se_renv;
  bvars->bse_tms[bmember] = svars->se_tms;
  bvars->bse_bcm_B[bmember] = svars->se_scm_B;
  bvars->bse_bcm_formation[bmember] = svars->se_scm_formation;
}

/**
* @brief ?
*
* @param bvars ?
* @param svars ?
* @param bmember ?
*/
void update_sse_from_bse(bse_t *bvars, sse_t *svars, int bmember) {
  svars->se_mass = bvars->bse_mass0[bmember];
  svars->se_k = bvars->bse_kw[bmember];
  svars->se_mt =  bvars->bse_mass[bmember];
  svars->se_bhspin =  bvars->bse_bhspin[bmember];
  svars->se_ospin = bvars->bse_ospin[bmember];
  svars->se_B_0 = bvars->bse_B_0[bmember];
  svars->se_bacc = bvars->bse_bacc[bmember];
  svars->se_tacc = bvars->bse_tacc[bmember];
  svars->se_epoch = bvars->bse_epoch[bmember];
  svars->se_tphys = bvars->bse_tphys; /* tphys should be the same for both input stars so this should be OK */
  svars->se_radius = bvars->bse_radius[bmember];
  svars->se_lum = bvars->bse_lum[bmember];
  svars->se_mc = bvars->bse_massc[bmember];
  svars->se_rc = bvars->bse_radc[bmember];
  svars->se_menv = bvars->bse_menv[bmember];
  svars->se_renv = bvars->bse_renv[bmember];
  svars->se_tms = bvars->bse_tms[bmember];
  svars->se_scm_B = bvars->bse_bcm_B[bmember];
  svars->se_scm_formation = bvars->bse_bcm_formation[bmember];
}

/**
* @brief ?
*
* @param svars ?
* @param star ?
*/
void get_sse_from_star(sse_t *svars, star_t *star) {
  svars-> se_mass   = star->se_mass;
  svars-> se_k      = star->se_k;
  svars-> se_mt     = star->se_mt;
  svars-> se_bhspin = star->se_bhspin;
  svars-> se_ospin  = star->se_ospin;
  svars-> se_B_0  = star->se_B_0;
  svars-> se_bacc  = star->se_bacc;
  svars-> se_tacc  = star->se_tacc;
  svars-> se_epoch  = star->se_epoch;
  svars-> se_tphys  = star->se_tphys;
  svars-> se_radius = star->se_radius;
  svars-> se_lum    = star->se_lum;
  svars-> se_mc     = star->se_mc;
  svars-> se_rc     = star->se_rc;
  svars-> se_menv   = star->se_menv;
  svars-> se_renv   = star->se_renv;
  svars-> se_tms    = star->se_tms;
  svars-> se_scm_B  = star->se_scm_B;
  svars-> se_scm_formation  = star->se_scm_formation;
}

/**
* @brief ?
*
* @param star ?
* @param svars ?
*/
void update_star_from_sse(star_t *star, sse_t svars) {
  star-> se_mass   = svars.se_mass;
  star-> se_k      = svars.se_k;
  star-> se_mt     = svars.se_mt;
  star-> se_bhspin = svars.se_bhspin;
  star-> se_ospin  = svars.se_ospin;
  star-> se_B_0  = svars.se_B_0;
  star-> se_bacc = svars.se_bacc;
  star-> se_tacc = svars.se_tacc;
  star-> se_epoch  = svars.se_epoch;
  star-> se_tphys  = svars.se_tphys;
  star-> se_radius = svars.se_radius;
  star-> se_lum    = svars.se_lum;
  star-> se_mc     = svars.se_mc;
  star-> se_rc     = svars.se_rc;
  star-> se_menv   = svars.se_menv;
  star-> se_renv   = svars.se_renv;
  star-> se_tms    = svars.se_tms;
  star-> se_scm_B = svars.se_scm_B;
  star-> se_scm_formation = svars.se_scm_formation;
}

/**
* @brief ?
*
* @param bvars ?
* @param star ?
*/
void get_bse_from_binary(bse_t *bvars, binary_t *star) {
  int i;

  for (i=0; i<2; i++) {
    bvars->bse_kw[i]        = star->bse_kw[i];
    bvars->bse_mass0[i]     = star->bse_mass0[i];
    bvars->bse_mass[i]      = star->bse_mass[i];
    bvars->bse_radius[i]    = star->bse_radius[i];
    bvars->bse_lum[i]       = star->bse_lum[i];
    bvars->bse_massc[i]     = star->bse_massc[i];
    bvars->bse_radc[i]      = star->bse_radc[i];
    bvars->bse_menv[i]      = star->bse_menv[i];
    bvars->bse_renv[i]      = star->bse_renv[i];
    bvars->bse_ospin[i]     = star->bse_ospin[i];
    bvars->bse_B_0[i]     = star->bse_B_0[i];
    bvars->bse_bacc[i]     = star->bse_bacc[i];
    bvars->bse_tacc[i]     = star->bse_tacc[i];
    bvars->bse_epoch[i]     = star->bse_epoch[i];
    bvars->bse_tms[i]       = star->bse_tms[i];
    bvars->bse_bhspin[i]       = star->bse_bhspin[i];
    bvars->bse_bcm_dmdt[i]  = star->bse_bcm_dmdt[i];
    bvars->bse_bcm_radrol[i]= star->bse_bcm_radrol[i];
    bvars->bse_bcm_B[i]     = star->bse_bcm_B[i];
    bvars->bse_bcm_formation[i]     = star->bse_bcm_formation[i];
  }
  bvars->bse_tphys        = star->bse_tphys;
  bvars->bse_tb           = star->bse_tb;
}

/**
* @brief ?
*
* @param star ?
* @param bvars ?
*/
void update_binary_from_bse(binary_t *star, bse_t *bvars) {
  int i;

  for (i=0; i<2; i++) {
    star->bse_kw[i]        =bvars->bse_kw[i];
    star->bse_mass0[i]     =bvars->bse_mass0[i];
    star->bse_mass[i]      =bvars->bse_mass[i];
    star->bse_radius[i]    =bvars->bse_radius[i];
    star->bse_lum[i]       =bvars->bse_lum[i];
    star->bse_massc[i]     =bvars->bse_massc[i];
    star->bse_radc[i]      =bvars->bse_radc[i];
    star->bse_menv[i]      =bvars->bse_menv[i];
    star->bse_renv[i]      =bvars->bse_renv[i];
    star->bse_ospin[i]     =bvars->bse_ospin[i];
    star->bse_B_0[i]     =bvars->bse_B_0[i];
    star->bse_bacc[i]     =bvars->bse_bacc[i];
    star->bse_tacc[i]     =bvars->bse_tacc[i];
    star->bse_epoch[i]     =bvars->bse_epoch[i];
    star->bse_tms[i]       =bvars->bse_tms[i];
    star->bse_bhspin[i]       =bvars->bse_bhspin[i];
    star->bse_bcm_dmdt[i]  =bvars->bse_bcm_dmdt[i];
    star->bse_bcm_radrol[i]=bvars->bse_bcm_radrol[i];
    star->bse_bcm_B[i]     =bvars->bse_bcm_B[i];
    star->bse_bcm_formation[i]     =bvars->bse_bcm_formation[i];
  }
  star->bse_tphys        =bvars->bse_tphys;
  star->bse_tb           =bvars->bse_tb;
}

/** 
 * @brief Breaks a binary if it contains a k=15 "star", i.e. a massless remnant.
 *
 * Use this function only if STELLAR_EVOLUTION is turned on (Spares one an 
 * extra if-statement)!
 *
 * @param bincom center-of-mass properties
 * @param bin binary properties
 */
void compress_binary(star_t *bincom, binary_t *bin) {
  bse_t bse_vars;
  sse_t sse_vars;
  long rem_idx, rem_id;
  double rem_lifetime, rem_createtime, rem_Eint;
  
  if(!STELLAR_EVOLUTION)
    return;

  if (bin->bse_mass[0]==0) {
    rem_idx= 1;
    rem_id= bin->id2;
    rem_lifetime= bin->lifetime_m2;
    rem_createtime= bin->createtime_m2;
    rem_Eint= bin->Eint2;
  } else if (bin->bse_mass[1]==0) {
    rem_idx= 0;
    rem_id= bin->id1;
    rem_lifetime= bin->lifetime_m1;
    rem_createtime= bin->createtime_m1;
    rem_Eint= bin->Eint1;
  } else {
    return;
  }

  /* Copy the BSE vars to the bincom SSE vars. */
  get_bse_from_binary(&bse_vars, bin);
  update_sse_from_bse(&bse_vars, &sse_vars, rem_idx);
  update_star_from_sse(bincom, sse_vars);

  /* Copy the remaining vars from surviving bin member. */
  bincom->m= bin->bse_mass[rem_idx] / units.mstar * MSUN;
  bincom->rad= bin->bse_radius[rem_idx] / units.l * RSUN;
  bincom->id= rem_id;
  bincom->Eint= rem_Eint;

  /* Vars from Sourav's toy rejuvenation prescription */
  bincom->lifetime= rem_lifetime;
  bincom->createtime= rem_createtime;

  /* guess we are done ... */
  dprintf("Compressed binary with idx=%li, id1=%li and id2=%li.\n",
      bincom->binind, bin->id1, bin->id2);

  destroy_binary(bincom->binind);
  bincom->binind= 0;
}


/**
* @brief ?
*
* @param tempbinary ?
* @param cmc_l_unit ?
* @param RbloodySUN ?
* @param zpars ?
* @param vs ?
* @param fb ?
* @param ecsnp ?
* @param ecsn_mlow ?
* @param ST_tide ?
*/
void cmc_bse_comenv(binary_t *tempbinary, double cmc_l_unit, double RbloodySUN, double *zpars, double *vs, int *fb)
		//double *M0, double *M, double *MC, double *AJ, double *OSPIN, int *KW, 
		//                double *M02, double *M2, double *MC2, double *AJ2, double *JSPIN2, int *KW2,
                //double *ZPARS, double *ECC, double *SEP, double *PORB,  
{
  int ii,jj;
  bse_binary binary;
  binary_t temphold;

  //Because the FORTRAN routine, comenv.f, takes star one as the initiator of the CE we must ensure
  //that the first star is a giant or the most massive giant.
  if((tempbinary->bse_kw[0])>=10 || (tempbinary->bse_kw[0]) <= 1 || (tempbinary->bse_kw[0] == 7)){
    //switch
    temphold = *tempbinary;
    for(ii=0;ii<=1;ii++){
      jj = abs(ii-1);
      tempbinary->bse_epoch[ii] = temphold.bse_epoch[jj];
      tempbinary->bse_mass0[ii] = temphold.bse_mass0[jj];
      tempbinary->bse_kw[ii] = temphold.bse_kw[jj];
      tempbinary->bse_mass[ii] = temphold.bse_mass[jj];
      tempbinary->bse_ospin[ii] = temphold.bse_ospin[jj];
      tempbinary->bse_B_0[ii] = temphold.bse_B_0[jj];
      tempbinary->bse_bacc[ii] = temphold.bse_bacc[jj];
      tempbinary->bse_tacc[ii] = temphold.bse_tacc[jj];
      tempbinary->bse_bcm_B[ii] = temphold.bse_bcm_B[jj];
      tempbinary->bse_bcm_formation[ii] = temphold.bse_bcm_formation[jj];
      tempbinary->bse_radius[ii] = temphold.bse_radius[jj];
      tempbinary->bse_lum[ii] = temphold.bse_lum[jj];
      tempbinary->bse_massc[ii] = temphold.bse_massc[jj];
      tempbinary->bse_radc[ii] = temphold.bse_radc[jj];
      tempbinary->bse_menv[ii] = temphold.bse_menv[jj];
      tempbinary->bse_renv[ii] = temphold.bse_renv[jj];
      tempbinary->bse_tms[ii] = temphold.bse_tms[jj];
      tempbinary->bse_bhspin[ii] = temphold.bse_bhspin[jj];
      tempbinary->bse_bcm_radrol[ii] = temphold.bse_bcm_radrol[jj];
      tempbinary->bse_bcm_dmdt[ii] = temphold.bse_bcm_dmdt[jj];
    }
  } else if ((tempbinary->bse_kw[1]>=2 && tempbinary->bse_kw[1]<=9 && tempbinary->bse_kw[1]!=7) && (tempbinary->bse_mass[1]>tempbinary->bse_mass[0])) {
    //switch
    temphold = *tempbinary;
    for(ii=0;ii<=1;ii++){
      jj = abs(ii-1);
      tempbinary->bse_epoch[ii] = temphold.bse_epoch[jj];
      tempbinary->bse_mass0[ii] = temphold.bse_mass0[jj];
      tempbinary->bse_kw[ii] = temphold.bse_kw[jj];
      tempbinary->bse_mass[ii] = temphold.bse_mass[jj];
      tempbinary->bse_ospin[ii] = temphold.bse_ospin[jj];
      tempbinary->bse_B_0[ii] = temphold.bse_B_0[jj];
      tempbinary->bse_bacc[ii] = temphold.bse_bacc[jj];
      tempbinary->bse_tacc[ii] = temphold.bse_tacc[jj];
      tempbinary->bse_bcm_B[ii] = temphold.bse_bcm_B[jj];
      tempbinary->bse_bcm_formation[ii] = temphold.bse_bcm_formation[jj];
      tempbinary->bse_radius[ii] = temphold.bse_radius[jj];
      tempbinary->bse_lum[ii] = temphold.bse_lum[jj];
      tempbinary->bse_massc[ii] = temphold.bse_massc[jj];
      tempbinary->bse_radc[ii] = temphold.bse_radc[jj];
      tempbinary->bse_menv[ii] = temphold.bse_menv[jj];
      tempbinary->bse_renv[ii] = temphold.bse_renv[jj];
      tempbinary->bse_tms[ii] = temphold.bse_tms[jj];
      tempbinary->bse_bhspin[ii] = temphold.bse_bhspin[jj];
      tempbinary->bse_bcm_radrol[ii] = temphold.bse_bcm_radrol[jj];
      tempbinary->bse_bcm_dmdt[ii] = temphold.bse_bcm_dmdt[jj];
    }
  }
  
  binary.a = tempbinary->a * cmc_l_unit/RSUN; //separation in solar radii as used in comenv.f
  binary.bse_tb = tempbinary->bse_tb;
  binary.e = tempbinary->e;
  binary.bse_tphys = tempbinary->bse_tphys;
  binary.bse_epoch[0] = tempbinary->bse_epoch[0];
  binary.bse_epoch[1] = tempbinary->bse_epoch[1];
  binary.bse_mass[0] = tempbinary->bse_mass[0];
  binary.bse_mass[1] = tempbinary->bse_mass[1];
  binary.bse_mass0[0] = tempbinary->bse_mass0[0];
  binary.bse_mass0[1] = tempbinary->bse_mass0[1];
  binary.bse_kw[0] = tempbinary->bse_kw[0];
  binary.bse_kw[1] = tempbinary->bse_kw[1];
  binary.bse_bhspin[0] = tempbinary->bse_bhspin[0];
  binary.bse_bhspin[1] = tempbinary->bse_bhspin[1];
  binary.bse_radius[0] = tempbinary->bse_radius[0];
  binary.bse_radius[1] = tempbinary->bse_radius[1];
  binary.bse_lum[0] = tempbinary->bse_lum[0];
  binary.bse_lum[1] = tempbinary->bse_lum[1];
  binary.bse_massc[0] = tempbinary->bse_massc[0];
  binary.bse_massc[1] = tempbinary->bse_massc[1];
  binary.bse_radc[0] = tempbinary->bse_radc[0];
  binary.bse_radc[1] = tempbinary->bse_radc[1];
  binary.bse_menv[0] = tempbinary->bse_menv[0];
  binary.bse_menv[1] = tempbinary->bse_menv[1];
  binary.bse_renv[0] = tempbinary->bse_renv[0];
  binary.bse_renv[1] = tempbinary->bse_renv[1];
  binary.bse_ospin[0] = tempbinary->bse_ospin[0];
  binary.bse_ospin[1] = tempbinary->bse_ospin[1];
  binary.bse_B_0[0] = tempbinary->bse_B_0[0];
  binary.bse_B_0[1] = tempbinary->bse_B_0[1];
  binary.bse_bacc[0] = tempbinary->bse_bacc[0];
  binary.bse_bacc[1] = tempbinary->bse_bacc[1];
  binary.bse_tacc[0] = tempbinary->bse_tacc[0];
  binary.bse_tacc[1] = tempbinary->bse_tacc[1];
  binary.bse_bcm_B[0] = tempbinary->bse_bcm_B[0];
  binary.bse_bcm_B[1] = tempbinary->bse_bcm_B[1];
  binary.bse_bcm_formation[0] = tempbinary->bse_bcm_formation[0];
  binary.bse_bcm_formation[1] = tempbinary->bse_bcm_formation[1];
  //
  bse_comenv(&(binary), zpars, vs, fb);
  //
  tempbinary->a = binary.a * RSUN/cmc_l_unit;
  tempbinary->bse_tb = binary.bse_tb;
  tempbinary->e = binary.e;
  tempbinary->bse_tphys = binary.bse_tphys;
  tempbinary->bse_epoch[0] = binary.bse_epoch[0];
  tempbinary->bse_epoch[1] = binary.bse_epoch[1];
  tempbinary->bse_mass[0] = binary.bse_mass[0];
  tempbinary->bse_mass[1] = binary.bse_mass[1];
  tempbinary->bse_mass0[0] = binary.bse_mass0[0];
  tempbinary->bse_mass0[1] = binary.bse_mass0[1];
  tempbinary->bse_kw[0] = binary.bse_kw[0];
  tempbinary->bse_kw[1] = binary.bse_kw[1];
  tempbinary->bse_bhspin[0] = binary.bse_bhspin[0];
  tempbinary->bse_bhspin[1] = binary.bse_bhspin[1];
  tempbinary->bse_radius[0] = binary.bse_radius[0];
  tempbinary->bse_radius[1] = binary.bse_radius[1];
  tempbinary->bse_lum[0] = binary.bse_lum[0];
  tempbinary->bse_lum[1] = binary.bse_lum[1];
  tempbinary->bse_massc[0] = binary.bse_massc[0];
  tempbinary->bse_massc[1] = binary.bse_massc[1];
  tempbinary->bse_radc[0] = binary.bse_radc[0];
  tempbinary->bse_radc[1] = binary.bse_radc[1];
  tempbinary->bse_menv[0] = binary.bse_menv[0];
  tempbinary->bse_menv[1] = binary.bse_menv[1];
  tempbinary->bse_renv[0] = binary.bse_renv[0];
  tempbinary->bse_renv[1] = binary.bse_renv[1];
  tempbinary->bse_ospin[0] = binary.bse_ospin[0];
  tempbinary->bse_ospin[1] = binary.bse_ospin[1];
  tempbinary->bse_B_0[0] = binary.bse_B_0[0];
  tempbinary->bse_bacc[0] = binary.bse_bacc[0];
  tempbinary->bse_tacc[0] = binary.bse_tacc[0];
  tempbinary->bse_bcm_B[0] = binary.bse_bcm_B[0];
  tempbinary->bse_bcm_formation[0] = binary.bse_bcm_formation[0];
  tempbinary->bse_bcm_formation[1] = binary.bse_bcm_formation[1];
  dprintf("after bse_comenv in cmc_bse_comenv: \n");
  dprintf("after hrdiag: kw1i=%d kw2i=%d m1f=%g m2f=%g r1f=%g r2f=%g epoch1=%g epoch2=%g \n", (*tempbinary).bse_kw[0], (*tempbinary).bse_kw[1], (*tempbinary).bse_mass[0], (*tempbinary).bse_mass[1], (*tempbinary).bse_radius[0], (*tempbinary).bse_radius[1], (*tempbinary).bse_epoch[0], (*tempbinary).bse_epoch[1]);
}
