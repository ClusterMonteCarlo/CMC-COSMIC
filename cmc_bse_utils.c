/* vi: set filetype=c.doxygen: */

#include "cmc.h"
#include "cmc_vars.h"

void update_bse_from_sse(bse_t *bvars, sse_t *svars, int bmember) {
  bvars->bse_mass0[bmember] = svars->se_mass;
  bvars->bse_kw[bmember] = svars->se_k;
  bvars->bse_mass[bmember] = svars->se_mt;
  bvars->bse_ospin[bmember] = svars->se_ospin;
  bvars->bse_epoch[bmember] = svars->se_epoch;
  bvars->bse_tphys = svars->se_tphys; /* tphys should be the same for both input stars so this should be OK */
  bvars->bse_radius[bmember] = svars->se_radius;
  bvars->bse_lum[bmember] = svars->se_lum;
  bvars->bse_massc[bmember] = svars->se_mc;
  bvars->bse_radc[bmember] = svars->se_rc;
  bvars->bse_menv[bmember] = svars->se_menv;
  bvars->bse_renv[bmember] = svars->se_renv;
  bvars->bse_tms[bmember] = svars->se_tms;
}

void update_sse_from_bse(bse_t *bvars, sse_t *svars, int bmember) {
  svars->se_mass = bvars->bse_mass0[bmember];
  svars->se_k = bvars->bse_kw[bmember];
  svars->se_mt =  bvars->bse_mass[bmember];
  svars->se_ospin = bvars->bse_ospin[bmember];
  svars->se_epoch = bvars->bse_epoch[bmember];
  svars->se_tphys = bvars->bse_tphys; /* tphys should be the same for both input stars so this should be OK */
  svars->se_radius = bvars->bse_radius[bmember];
  svars->se_lum = bvars->bse_lum[bmember];
  svars->se_mc = bvars->bse_massc[bmember];
  svars->se_rc = bvars->bse_radc[bmember];
  svars->se_menv = bvars->bse_menv[bmember];
  svars->se_renv = bvars->bse_renv[bmember];
  svars->se_tms = bvars->bse_tms[bmember];
}

void get_sse_from_star(sse_t *svars, star_t *star) {
  svars-> se_mass   = star->se_mass;
  svars-> se_k      = star->se_k;
  svars-> se_mt     = star->se_mt;
  svars-> se_ospin  = star->se_ospin;
  svars-> se_epoch  = star->se_epoch;
  svars-> se_tphys  = star->se_tphys;
  svars-> se_radius = star->se_radius;
  svars-> se_lum    = star->se_lum;
  svars-> se_mc     = star->se_mc;
  svars-> se_rc     = star->se_rc;
  svars-> se_menv   = star->se_menv;
  svars-> se_renv   = star->se_renv;
  svars-> se_tms    = star->se_tms;
}

void update_star_from_sse(star_t *star, sse_t svars) {
  star-> se_mass   = svars.se_mass;
  star-> se_k      = svars.se_k;
  star-> se_mt     = svars.se_mt;
  star-> se_ospin  = svars.se_ospin;
  star-> se_epoch  = svars.se_epoch;
  star-> se_tphys  = svars.se_tphys;
  star-> se_radius = svars.se_radius;
  star-> se_lum    = svars.se_lum;
  star-> se_mc     = svars.se_mc;
  star-> se_rc     = svars.se_rc;
  star-> se_menv   = svars.se_menv;
  star-> se_renv   = svars.se_renv;
  star-> se_tms    = svars.se_tms;
}

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
    bvars->bse_epoch[i]     = star->bse_epoch[i];
    bvars->bse_tms[i]       = star->bse_tms[i];
    bvars->bse_bcm_dmdt[i]  = star->bse_bcm_dmdt[i];
    bvars->bse_bcm_radrol[i]= star->bse_bcm_radrol[i];
  }
  bvars->bse_tphys        = star->bse_tphys;
  bvars->bse_tb           = star->bse_tb;
}

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
    star->bse_epoch[i]     =bvars->bse_epoch[i];
    star->bse_tms[i]       =bvars->bse_tms[i];
    star->bse_bcm_dmdt[i]  =bvars->bse_bcm_dmdt[i];
    star->bse_bcm_radrol[i]=bvars->bse_bcm_radrol[i];
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

