/* vi: set filetype=c.doxygen: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cmc.h"
#include "cmc_vars.h"
#include "bse_wrap.h"

/**
* @brief This is the same as the stellar_evolution_init, except we're only
* setting the global BSE properties (the individual star and binary properties
* were set in the restart binary files)
*/
void restart_stellar_evolution(void){
  bse_set_using_cmc();
  bse_set_neta(BSE_NETA);
  bse_set_bwind(BSE_BWIND);
  bse_set_hewind(BSE_HEWIND);
  bse_set_windflag(BSE_WINDFLAG);
  bse_set_eddlimflag(BSE_EDDLIMFLAG);
  bse_set_pisn(BSE_PISN);
  bse_set_ecsn(BSE_ECSN);
  bse_set_ecsn_mlow(BSE_ECSN_MLOW);
  bse_set_aic(BSE_AIC);
  bse_set_bdecayfac(BSE_BDECAYFAC);
  bse_set_st_cr(BSE_ST_CR);
  bse_set_st_tide(BSE_ST_TIDE);
  bse_set_htpmb(BSE_HTPMB);
  bse_set_rejuvflag(BSE_REJUVFLAG);
  bse_set_ussn(BSE_USSN);
  bse_set_qcrit_array(bse_qcrit_array, NO_BSE_QCRIT_ARRAY); 
  bse_set_fprimc_array(bse_fprimc_array, NO_BSE_FPRIMC_ARRAY);
  bse_set_natal_kick_array(bse_natal_kick_array, NO_BSE_NATAL_KICK_ARRAY); 
  bse_set_sigmadiv(BSE_SIGMADIV);
  bse_set_alpha1(BSE_ALPHA1); /* FIXME: is 3 too high? (normally 1.0) */
  bse_set_lambdaf(BSE_LAMBDAF);
  bse_set_ceflag(BSE_CEFLAG);
  bse_set_cehestarflag(BSE_CEHESTARFLAG);
  bse_set_cemergeflag(BSE_CEMERGEFLAG);
  bse_set_cekickflag(BSE_CEKICKFLAG);
  bse_set_tflag(BSE_TFLAG);
  bse_set_qcflag(BSE_QCFLAG);
  bse_set_ifflag(BSE_IFFLAG);
  bse_set_wdflag(BSE_WDFLAG);
  bse_set_bhflag(BSE_BHFLAG);
  bse_set_grflag(BSE_GRFLAG);
  bse_set_kickflag(BSE_KICKFLAG);
  bse_set_zsun(BSE_ZSUN);
  bse_set_rembar_massloss(BSE_REMBAR_MASSLOSS);
  bse_set_remnantflag(BSE_REMNANTFLAG);
  bse_set_bhspinflag(BSE_BHSPINFLAG);
  bse_set_bhms_coll_flag(BSE_BHMS_COLL_FLAG);
  bse_set_bhspinmag(BSE_BHSPINMAG);
  bse_set_mxns(BSE_MXNS); //3 if remnantflag=1 or 2, 1.8 if remnantflag=0 (see evolv2.f)
  bse_set_bconst(BSE_BCONST);
  bse_set_CK(BSE_CK);
  bse_set_rejuv_fac(BSE_REJUV_FAC);
  /* need to suppress the self-initialization of the Fortran Tausworthe generator */
  if (BSE_IDUM<0) {
    BSE_IDUM= -BSE_IDUM;
  }
  bse_set_idum(BSE_IDUM);
  bse_set_pts1(BSE_PTS1);
  bse_set_pts2(BSE_PTS2);
  bse_set_pts3(BSE_PTS3);
  bse_set_sigma(BSE_SIGMA);
  bse_set_bhsigmafrac(BSE_BHSIGMAFRAC);
  bse_set_polar_kick_angle(BSE_POLAR_KICK_ANGLE);
  bse_set_beta(BSE_BETA); //set -0.125 if variable beta (following startrack), otherwise 0.125 for bse.
  bse_set_xi(BSE_XI);
  bse_set_acc2(BSE_ACC2);
  bse_set_epsnov(BSE_EPSNOV);
  bse_set_eddfac(BSE_EDDFAC); /* (normally 1.0) */
  bse_set_gamma(BSE_GAMMA);
  bse_set_merger(-1.0);
  
  /* set parameters relating to metallicity */
  zpars = (double *) malloc(20 * sizeof(double));
  bse_zcnsts(&METALLICITY, zpars);

  /* set collisions matrix */
  bse_instar();

  /*Set rng to saved rng seed*/
  bse_set_taus113state(*curr_st, 0);
}

/**
* @brief Initializes stellar evolution, both the global BSE properties and the
* individual properties for each star.  Runs BSE for a dt->0 time to set a bunch
* of the parameters for individual stars
*/
void stellar_evolution_init(void){  
  double tphysf, dtp, vs[20];
  int i;
  long k, kb;
  int kprev0=-100;
  int kprev1=-100;
  binary_t tempbinary;

  /* SSE */
  /* bse_set_hewind(0.5); */

  /* BSE */
  bse_set_using_cmc();
  bse_set_neta(BSE_NETA);
  bse_set_bwind(BSE_BWIND);
  bse_set_hewind(BSE_HEWIND);
  bse_set_windflag(BSE_WINDFLAG);
  bse_set_eddlimflag(BSE_EDDLIMFLAG);
  bse_set_pisn(BSE_PISN);
  bse_set_ecsn(BSE_ECSN);
  bse_set_ecsn_mlow(BSE_ECSN_MLOW);
  bse_set_aic(BSE_AIC);
  bse_set_bdecayfac(BSE_BDECAYFAC);
  bse_set_st_cr(BSE_ST_CR);
  bse_set_st_tide(BSE_ST_TIDE);
  bse_set_htpmb(BSE_HTPMB);
  bse_set_rejuvflag(BSE_REJUVFLAG);
  bse_set_ussn(BSE_USSN);
  bse_set_qcrit_array(bse_qcrit_array, NO_BSE_QCRIT_ARRAY); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
  bse_set_fprimc_array(bse_fprimc_array, NO_BSE_FPRIMC_ARRAY); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
  bse_set_natal_kick_array(bse_natal_kick_array, NO_BSE_NATAL_KICK_ARRAY); /* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1) */
  bse_set_sigmadiv(BSE_SIGMADIV);
  bse_set_alpha1(BSE_ALPHA1); /* FIXME: is 3 too high? (normally 1.0) */
  bse_set_lambdaf(BSE_LAMBDAF);
  bse_set_ceflag(BSE_CEFLAG);
  bse_set_cehestarflag(BSE_CEHESTARFLAG);
  bse_set_cemergeflag(BSE_CEMERGEFLAG);
  bse_set_cekickflag(BSE_CEKICKFLAG);
  bse_set_tflag(BSE_TFLAG);
  bse_set_qcflag(BSE_QCFLAG);
  bse_set_ifflag(BSE_IFFLAG);
  bse_set_wdflag(BSE_WDFLAG);
  bse_set_bhflag(BSE_BHFLAG);
  bse_set_grflag(BSE_GRFLAG);
  bse_set_kickflag(BSE_KICKFLAG);
  bse_set_zsun(BSE_ZSUN);
  bse_set_rembar_massloss(BSE_REMBAR_MASSLOSS);
  bse_set_remnantflag(BSE_REMNANTFLAG);
  bse_set_bhspinflag(BSE_BHSPINFLAG);
  bse_set_bhms_coll_flag(BSE_BHMS_COLL_FLAG);
  bse_set_bhspinmag(BSE_BHSPINMAG);
  bse_set_mxns(BSE_MXNS); //3 if remnantflag=1 or 2, 1.8 if remnantflag=0 (see evolv2.f)
  bse_set_bconst(BSE_BCONST);
  bse_set_CK(BSE_CK);
  bse_set_rejuv_fac(BSE_REJUV_FAC);
  /* need to suppress the self-initialization of the Fortran Tausworthe generator */
  if (BSE_IDUM<0) {
    BSE_IDUM= -BSE_IDUM;
  }
  bse_set_idum(BSE_IDUM);
  bse_set_pts1(BSE_PTS1);
  bse_set_pts2(BSE_PTS2);
  bse_set_pts3(BSE_PTS3);
  bse_set_sigma(BSE_SIGMA);
  bse_set_bhsigmafrac(BSE_BHSIGMAFRAC);
  bse_set_polar_kick_angle(BSE_POLAR_KICK_ANGLE);
  bse_set_beta(BSE_BETA); //set -0.125 if variable beta (following startrack), otherwise 0.125 for bse.
  bse_set_xi(BSE_XI);
  bse_set_acc2(BSE_ACC2);
  bse_set_epsnov(BSE_EPSNOV);
  bse_set_eddfac(BSE_EDDFAC); /* (normally 1.0) */
  bse_set_gamma(BSE_GAMMA);
  bse_set_merger(-1.0);
  
  /* set parameters relating to metallicity */
  zpars = (double *) malloc(20 * sizeof(double));
  bse_zcnsts(&METALLICITY, zpars);

  /* set collisions matrix */
  bse_instar();
  dprintf("se_init: %g %g %g %d %g %g %g %d %d %d %d %d %d %g %d %g %g %g %g %g %g\n", BSE_NETA, BSE_BWIND, BSE_HEWIND, BSE_WINDFLAG, BSE_PISN, BSE_ALPHA1, BSE_LAMBDAF, BSE_CEFLAG, BSE_TFLAG, BSE_IFFLAG, BSE_WDFLAG, BSE_BHFLAG, BSE_REMNANTFLAG, BSE_MXNS, BSE_IDUM, BSE_SIGMA, BSE_BHSIGMAFRAC, BSE_BETA, BSE_EDDFAC, BSE_GAMMA, BSE_POLAR_KICK_ANGLE);

  for (k=1; k<=mpiEnd-mpiBegin+1; k++) {
    long g_k = get_global_idx(k);

    if (star[k].binind == 0) { /* single star */
      star[k].se_mass = star_m[g_k] * units.mstar / MSUN;
	  star[k].zams_mass = star[k].se_mass;
      /* setting the type */
      if(star[k].se_mass <= 0.7){
        star[k].se_k = 0;
      } else {
        star[k].se_k = 1;
      }
      star[k].se_mt = star[k].se_mass;
      star[k].se_ospin = 0.0;
      star[k].se_B_0 = 0.0; /* PK */
      star[k].se_bacc = 0.0;
      star[k].se_tacc = 0.0;
      star[k].se_epoch = 0.0;
      star[k].se_tphys = 0.0;

      /* evolve slightly (1 year) for initial radii */
      tphysf = 1.0e-6;
      dtp = tphysf - star[k].se_tphys;
      dtp = 0.0;
      DMse += star_m[g_k] * madhoc;
      /* Update star id for pass through. */
      bse_set_id1_pass(star[k].id);
      bse_set_id2_pass(0);
      tempbinary.bse_mass0[0] = star[k].se_mass;
      tempbinary.bse_mass0[1] = 0.0;
      tempbinary.bse_kw[0] = star[k].se_k;
      tempbinary.bse_kw[1] = 15;
      tempbinary.bse_mass[0] = star[k].se_mt;
      tempbinary.bse_mass[1] = 0.0;
      tempbinary.bse_radius[0] = star[k].se_radius;
      tempbinary.bse_radius[1] = 0.0;
      tempbinary.bse_lum[0] = star[k].se_lum;
      tempbinary.bse_lum[1] = 0.0;
      tempbinary.bse_massc[0] = star[k].se_mc;
      tempbinary.bse_massc[1] = 0.0;
      tempbinary.bse_radc[0] = star[k].se_rc;
      tempbinary.bse_radc[1] = 0.0;
      tempbinary.bse_menv[0] = star[k].se_menv;
      tempbinary.bse_menv[1] = 0.0;
      tempbinary.bse_renv[0] = star[k].se_renv;
      tempbinary.bse_renv[1] = 0.0;
      tempbinary.bse_ospin[0] = star[k].se_ospin;
      tempbinary.bse_ospin[1] = 0.0;
      tempbinary.bse_B_0[0] = star[k].se_B_0;
      tempbinary.bse_B_0[1] = 0.0;
      tempbinary.bse_bacc[0] = star[k].se_bacc;
      tempbinary.bse_bacc[1] = 0.0;
      tempbinary.bse_tacc[0] = star[k].se_tacc;
      tempbinary.bse_tacc[1] = 0.0;
      tempbinary.bse_epoch[0] = star[k].se_epoch;
      tempbinary.bse_epoch[1] = 0.0;
      tempbinary.bse_tms[0] = star[k].se_tms;
      tempbinary.bse_tms[1] = 0.0;
      tempbinary.bse_bhspin[0] = star[k].se_bhspin;
      tempbinary.bse_bhspin[1] = 0.0;
      tempbinary.bse_tb = 0.0;
      tempbinary.e = 0.0;

      /*
        bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
        &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
        &(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
        &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
       */

//MPI: To allow the serial version to mimic the parallel, we draw random numbers from the appropriate streams. And hence, we switch the current random state to appropriate state as in a parallel run.
      bse_set_taus113state(*curr_st, 0);
      bse_evolv2(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
          &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
          &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
          &(tempbinary.bse_ospin[0]), &(tempbinary.bse_B_0[0]), &(tempbinary.bse_bacc[0]), &(tempbinary.bse_tacc[0]),
          &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
          &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
          &(tempbinary.bse_tb), &(tempbinary.e), vs,&(tempbinary.bse_bhspin[0]));
      *curr_st=bse_get_taus113state();

      star[k].se_mass = tempbinary.bse_mass0[0];
      star[k].se_k = tempbinary.bse_kw[0];
      star[k].se_mt = tempbinary.bse_mass[0];
      star[k].se_radius = tempbinary.bse_radius[0];
      star[k].se_lum = tempbinary.bse_lum[0];
      star[k].se_mc = tempbinary.bse_massc[0];
      star[k].se_rc = tempbinary.bse_radc[0];
      star[k].se_menv = tempbinary.bse_menv[0];
      star[k].se_renv = tempbinary.bse_renv[0];
      star[k].se_ospin = tempbinary.bse_ospin[0];
      star[k].se_B_0 = tempbinary.bse_B_0[0];
      star[k].se_bacc = tempbinary.bse_bacc[0];
      star[k].se_tacc = tempbinary.bse_tacc[0];
      star[k].se_epoch = tempbinary.bse_epoch[0];
      star[k].se_tms = tempbinary.bse_tms[0];
	  star[k].se_bhspin = tempbinary.bse_bhspin[0];

      star[k].rad = star[k].se_radius * RSUN / units.l;
      star_m[g_k] = star[k].se_mt * MSUN / units.mstar;
      DMse -= star_m[g_k] * madhoc;
      /* birth kicks */
      if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[2]*vs[2]) != 0.0) {
        //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
      }
      star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);


      vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
      //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
      set_star_EJ(k);
    } else if (star[k].binind > 0) { /* binary */
      star[k].se_k = NOT_A_STAR; /* just for safety */
      kb = star[k].binind;
      binary[kb].bse_mass0[0] = binary[kb].m1 * units.mstar / MSUN;
      binary[kb].bse_mass0[1] = binary[kb].m2 * units.mstar / MSUN;
	  binary[kb].bse_zams_mass[0] = binary[kb].bse_mass0[0];
	  binary[kb].bse_zams_mass[1] = binary[kb].bse_mass0[1];
      for (i=0; i<=1; i++) {
        if(binary[kb].bse_mass0[i] <= 0.7){
          binary[kb].bse_kw[i] = 0;
        } else {
          binary[kb].bse_kw[i] = 1;
        }
        binary[kb].bse_mass[i] = binary[kb].bse_mass0[i];
        binary[kb].bse_ospin[i] = 0.0;
        binary[kb].bse_B_0[i] = 0.0;
        binary[kb].bse_bacc[i] = 0.0;
        binary[kb].bse_tacc[i] = 0.0;
        binary[kb].bse_epoch[i] = 0.0;
        binary[kb].bse_bhspin[i] = 0.0;
      }
      binary[kb].bse_tphys = 0.0;

      /* set binary orbital period (in days) from a */
      binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;

      /* evolve slightly (1 year) for initial radii */
      tphysf = 1.0e-6;
      dtp = tphysf - binary[kb].bse_tphys;
      dtp = 0.0;
      DMse += (binary[kb].m1 + binary[kb].m2) * madhoc;
      /* Update star id for pass through. */
      bse_set_id1_pass(binary[kb].id1);
      bse_set_id2_pass(binary[kb].id2);
      bse_set_taus113state(*curr_st, 0);
      bse_evolv2(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
              &(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
              &(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]), 
              &(binary[kb].bse_B_0[0]), &(binary[kb].bse_bacc[0]), &(binary[kb].bse_tacc[0]),
              &(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
              &(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
              &(binary[kb].bse_tb), &(binary[kb].e), vs,&(binary[kb].bse_bhspin[0]));
      *curr_st=bse_get_taus113state();

      handle_bse_outcome(k, kb, vs, tphysf, kprev0, kprev1);
    } else {
      eprintf("totally confused!\n");
      exit_cleanly(-1, __FUNCTION__);
    }
  }

	double tmpTimeStart = timeStartSimple();
  //if(myid==0)
    MPI_Allreduce(MPI_IN_PLACE, &DMse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //else
  //  MPI_Allreduce(&DMse, &DMse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	timeEndSimple(tmpTimeStart, &t_comm);
}

/* note that this routine is called after perturb_stars() and get_positions() */
/**
* @brief does stellar evolution using sse and bse packages.
*
* @param rng gsl rng
*/
void do_stellar_evolution(gsl_rng *rng)
{
  long k, kb, j, jj;
  int kprev,i, ii;
  int kprev0, kprev1;
  double dtp, tphysf, vs[20], VKO;
  double M_beforeSE, M10_beforeSE, M100_beforeSE, M1000_beforeSE, Mcore_beforeSE;
  double M_afterSE, M10_afterSE, M100_afterSE, M1000_afterSE, Mcore_afterSE;
  double r10_beforeSE, r100_beforeSE, r1000_beforeSE, rcore_beforeSE;
  double dM_dt_SE10, dM_dt_SE100, dM_dt_SE1000, dM_dt_SEcore; 
  struct rng_t113_state temp_state;
  int reduced_timestep=0;
  binary_t tempbinary;
  bse_set_merger(-1.0);
  /* double vk, theta; */

  //MPI: The serial version runs till N_MAX_NEW+1 to account for the sentinel. But in the parallel version, there is no sentinel, so runs only till N_MAX_NEW.
  for(k=1; k<=clus.N_MAX_NEW; k++){ 
    int g_k = get_global_idx(k);
    if (star[k].binind == 0) { /* single star */
      tphysf = TotalTime / MEGA_YEAR;
      dtp = tphysf;
      dtp = 0.0;
      kprev = star[k].se_k;   
      kprev0 = -100; /* set the previous stellar type variable for binary, just so they are initialized) */
      kprev1 = -100;
      if (star_m[get_global_idx(k)]<=DBL_MIN && star[k].vr==0. && star[k].vt==0. && star[k].E==0. && star[k].J==0.){ //ignoring zeroed out stars
        dprintf ("zeroed out star: skipping SE:\n"); 
        dprintf ("k=%ld m=%g r=%g phi=%g vr=%g vt=%g E=%g J=%g\n", k, star_m[g_k], star_r[g_k], star_phi[g_k], star[k].vr, star[k].vt, star[k].E, star[k].J);
      } else {
        DMse += star_m[g_k] * madhoc;
        /* Update star id for pass through. */
        bse_set_id1_pass(star[k].id);
        bse_set_id2_pass(0);
        tempbinary.bse_mass0[0] = star[k].se_mass;
        tempbinary.bse_mass0[1] = 0.0;
        tempbinary.bse_kw[0] = star[k].se_k;
        tempbinary.bse_kw[1] = 15;
        tempbinary.bse_mass[0] = star[k].se_mt;
        tempbinary.bse_mass[1] = 0.0;
        tempbinary.bse_radius[0] = star[k].se_radius;
        tempbinary.bse_radius[1] = 0.0;
        tempbinary.bse_lum[0] = star[k].se_lum;
        tempbinary.bse_lum[1] = 0.0;
        tempbinary.bse_massc[0] = star[k].se_mc;
        tempbinary.bse_massc[1] = 0.0;
        tempbinary.bse_radc[0] = star[k].se_rc;
        tempbinary.bse_radc[1] = 0.0;
        tempbinary.bse_menv[0] = star[k].se_menv;
        tempbinary.bse_menv[1] = 0.0;
        tempbinary.bse_renv[0] = star[k].se_renv;
        tempbinary.bse_renv[1] = 0.0;
        tempbinary.bse_ospin[0] = star[k].se_ospin;
        tempbinary.bse_ospin[1] = 0.0;
        tempbinary.bse_B_0[0] = star[k].se_B_0;
        tempbinary.bse_B_0[1] = 0.0;
        tempbinary.bse_bacc[0] = star[k].se_bacc;
        tempbinary.bse_bacc[1] = 0.0;
        tempbinary.bse_tacc[0] = star[k].se_tacc;
        tempbinary.bse_tacc[1] = 0.0;
        tempbinary.bse_epoch[0] = star[k].se_epoch;
        tempbinary.bse_epoch[1] = 0.0;
        tempbinary.bse_tms[0] = star[k].se_tms;
        tempbinary.bse_tms[1] = 0.0;
        tempbinary.bse_bhspin[0] = star[k].se_bhspin;
        tempbinary.bse_bhspin[1] = 0.0;
        tempbinary.bse_tb = 0.0;
        tempbinary.e = 0.0;

		  /*If we've got a large MS star, we need to reduce the timestep, otherwise
		   * we miss the transition from MS to HG to giant, and won't start applying
		   * winds for massive stars at the right time*/
		  if(star[k].zams_mass > 18){
			  bse_set_pts1(BSE_PTS1/10.);
			  reduced_timestep = 1;
		  }
        /*
          bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
          &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
          &(star[k].se_renv), &(star[k].se_ospin), &(star[k].se_epoch), &(star[k].se_tms), 
          &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
         */
        bse_set_taus113state(*curr_st, 0);
        bse_evolv2_safely(&(tempbinary.bse_kw[0]), &(tempbinary.bse_mass0[0]), &(tempbinary.bse_mass[0]), 
            &(tempbinary.bse_radius[0]), &(tempbinary.bse_lum[0]), &(tempbinary.bse_massc[0]), 
            &(tempbinary.bse_radc[0]), &(tempbinary.bse_menv[0]), &(tempbinary.bse_renv[0]), 
            &(tempbinary.bse_ospin[0]), &(tempbinary.bse_B_0[0]), &(tempbinary.bse_bacc[0]), &(tempbinary.bse_tacc[0]), 
            &(tempbinary.bse_epoch[0]), &(tempbinary.bse_tms[0]), 
            &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
            &(tempbinary.bse_tb), &(tempbinary.e), vs, &(tempbinary.bse_bhspin[0]));
        *curr_st=bse_get_taus113state();

        star[k].se_mass = tempbinary.bse_mass0[0];
        star[k].se_k = tempbinary.bse_kw[0];
        star[k].se_mt = tempbinary.bse_mass[0];
        star[k].se_radius = tempbinary.bse_radius[0];
        star[k].se_lum = tempbinary.bse_lum[0];
        star[k].se_mc = tempbinary.bse_massc[0];
        star[k].se_rc = tempbinary.bse_radc[0];
        star[k].se_menv = tempbinary.bse_menv[0];
        star[k].se_renv = tempbinary.bse_renv[0];
        star[k].se_ospin = tempbinary.bse_ospin[0];
        star[k].se_B_0 = tempbinary.bse_B_0[0];
        star[k].se_bacc = tempbinary.bse_bacc[0];
        star[k].se_tacc = tempbinary.bse_tacc[0];
        star[k].se_epoch = tempbinary.bse_epoch[0];
        star[k].se_tms = tempbinary.bse_tms[0];
	star[k].se_bhspin = tempbinary.bse_bhspin[0];

		  /*Reset the MS timestep once we're done*/
		  if(reduced_timestep == 1)
			  bse_set_pts1(BSE_PTS1);

        star[k].rad = star[k].se_radius * RSUN / units.l;
        star_m[g_k] = star[k].se_mt * MSUN / units.mstar;
        DMse -= star_m[g_k] * madhoc;

        /* extract info from scm array */ /* PK looping over a large number anticipating further changes */
	        i = 1;
	        j = 1;
        	while (bse_get_bcm(i,1)>=0.0 && i < 50000) {
          		if(i > 1) {
            			if(bse_get_bcm(i,2) == 13 && bse_get_bcm(i-1,2) < 13){
              				if(bse_get_bcm(i+1,1) >= 0.0){
                				star[k].se_scm_formation = bse_get_bcm(i+1,39);
                //printf("stel, pulsar sin. 1 : %g %g %g %g %ld \n",star[k].se_B_0,star[k].se_ospin,bse_get_bcm(i+1,33),star[k].se_scm_formation,star[k].id);
              				} else {
                				star[k].se_scm_formation = bse_get_bcm(i,39);
                //printf("stel, pulsar sin. 2 : %g %g %g %g %ld \n",star[k].se_B_0,star[k].se_ospin,bse_get_scm(i,33),star[k].se_scm_formation,star[k].id);
              				}
            			}
          		}
          		i++;
        	}
        	i--;
        	if(i+1 > 50000){
          		i = 0;
        	}
        	if(i>=1) {
          		if(i>1){
            			star[k].se_scm_B = bse_get_bcm(i,33);
            			if(star[k].se_k == 13) {
              //printf("stel, pulsar sin. 3 : %g %g %g \n",star[k].se_B_0,star[k].se_ospin,star[k].se_scm_B);
            			}
          		} else {
            			star[k].se_scm_B = bse_get_bcm(i,33);
            			if(star[k].se_k == 13) {
              //printf("stel, pulsar sin. 4 : %g %g %g \n",star[k].se_B_0,star[k].se_ospin,star[k].se_scm_B);
            			}
          		}
        	} else {
          /*         bse_evolv1(&(star[k].se_k), &(star[k].se_mass), &(star[k].se_mt), &(star[k].se_radius), 
                     &(star[k].se_lum), &(star[k].se_mc), &(star[k].se_rc), &(star[k].se_menv), 
                     &(star[k].se_renv), &(star[k].se_ospin), 
                     &(star[k].se_B_0), &(star[k].se_bacc), &(star[k].se_tacc), 
                     &(star[k].se_epoch), &(star[k].se_tms), 
                     &(star[k].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs); */
          		eprintf("Couldn't extract iso star bse info (looking for pulsar data)...");
          		eprintf("Evolv1 info from non scm extraction: k=%ld, kw=%d mass=%g mt=%g rad=%g lum=%g tphysf=%g dtp=%g ",k,star[k].se_k,star[k].se_mass,star[k].se_mt,star[k].se_radius,star[k].se_lum,tphysf,dtp);
          //           exit_cleanly(-1); //should only enter here if no bcm array entry, 
          //                               and that should only happen to uninteresting systems and/or outcomes.
        	}

        /* birth kicks */
        if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
          //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
        }
        star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);


        vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
        //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
        set_star_EJ(k);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
        /* birth kicks */
        if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
          //   dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
          dprintf("birth kick(iso): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g VKO=%.18g type=%d star_id=%ld Pi=%g\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[k].se_k, star[k].id,star[k].se_ospin);
        }

        /* WD birth kicks, just in case they exist */
        /* if ((star[k].se_k >= 10 && star[k].se_k <= 12) && star[k].se_k != kprev) { */
        /* vk = 2.0e5 / (units.l/units.t); */
        /*  theta = acos(2.0 * gsl_rng_uniform(rng) - 1.0); */
        /*  star[k].vr += cos(theta) * vk; */
        /*  star[k].vt += sin(theta) * vk; */
        /*  set_star_EJ(k); */
        /* } */
		
	/* PDK search for boom stuff and write pulsar data. */
		//bcm_boom_search(k, vs, getCMCvalues);
		if(WRITE_PULSAR_INFO)
		{
			pulsar_write(k, VKO);
		}

		//Shi
		if(tcount%PULSAR_DELTACOUNT==0){
			if (WRITE_MOREPULSAR_INFO){
				write_morepulsar(k);
			}
		}

		if (WRITE_BH_INFO) {
			if (kprev!=14 && star[k].se_k==14) { // newly formed BH
				parafprintf(newbhfile, "%.18g %g 0 %ld %g %g %g %g %g", TotalTime, star_r[g_k], star[k].id,star[k].zams_mass,star[k].se_mass, star[k].se_mt, star[k].se_bhspin, VKO);
				for (ii=0; ii<16; ii++){
					parafprintf (newbhfile, " %g", vs[ii]);
				}
				parafprintf (newbhfile, "\n");
//m_init, m_bh, time, id, kick, r, vr_init, vt_init, vr_final, vt_final, binflag, m0_init, m1_init, m0_final, m1_final, 
			}
		}

      }
    } else { /* binary */
      tphysf = TotalTime / MEGA_YEAR;
      kb = star[k].binind;
      dtp = tphysf - binary[kb].bse_tphys;
      dtp = 0.0;
      if (star_m[g_k]<=DBL_MIN && binary[kb].a==0. && binary[kb].e==0. && binary[kb].m1==0. && binary[kb].m2==0.){ //ignoring zeroed out binaries
          dprintf ("zeroed out star: skipping SE:\n");  
          dprintf ("k=%ld kb=%ld m=%g m1=%g m2=%g a=%g e=%g r=%g\n", k, kb, star_m[g_k], binary[kb].m1, binary[kb].m2, binary[kb].a, binary[kb].e, star_r[g_k]);
      } else {
		/* store previous star types for binary components, before evolving binary */
		kprev0=binary[kb].bse_kw[0];
		kprev1=binary[kb].bse_kw[1];

		/*If we've got a large MS star, we need to reduce the timestep, otherwise
		 * we miss the transition from MS to HG to giant, and won't start applying
		 * winds for massive stars at the right time*/
		if(binary[kb].bse_zams_mass[0] > 18 || binary[kb].bse_zams_mass[1] > 18){
		    bse_set_pts1(BSE_PTS1/10.);
		    reduced_timestep = 1;
		}

        /* set binary orbital period (in days) from a */
        binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;
        DMse += (binary[kb].m1 + binary[kb].m2) * madhoc;
        /* Update star id for pass through. */
        bse_set_id1_pass(binary[kb].id1);
        bse_set_id2_pass(binary[kb].id2);
		/* If this is a binary black hole, skip BSE and explicitly integrate the
		 * Peters equations*/
		if(binary[kb].bse_kw[0] == 14 && binary[kb].bse_kw[1] == 14){
			integrate_a_e_peters_eqn(kb);
			for (ii = 0 ; ii < 16 ; ii++) vs[ii] = 0.;
		} else{
			bse_set_taus113state(*curr_st, 0);
			bse_evolv2_safely(&(binary[kb].bse_kw[0]), &(binary[kb].bse_mass0[0]), &(binary[kb].bse_mass[0]), &(binary[kb].bse_radius[0]), 
				&(binary[kb].bse_lum[0]), &(binary[kb].bse_massc[0]), &(binary[kb].bse_radc[0]), &(binary[kb].bse_menv[0]), 
					&(binary[kb].bse_renv[0]), &(binary[kb].bse_ospin[0]),
						&(binary[kb].bse_B_0[0]), &(binary[kb].bse_bacc[0]), &(binary[kb].bse_tacc[0]),
				&(binary[kb].bse_epoch[0]), &(binary[kb].bse_tms[0]), 
				&(binary[kb].bse_tphys), &tphysf, &dtp, &METALLICITY, zpars, 
				&(binary[kb].bse_tb), &(binary[kb].e), vs, &(binary[kb].bse_bhspin[0]));
			*curr_st=bse_get_taus113state();
		}

		/*Reset the MS timestep once we're done*/
		if(reduced_timestep == 1)
		    bse_set_pts1(BSE_PTS1);

        if(isnan(binary[kb].bse_radius[0])){
		printf("id1=%ld id2=%ld\n",binary[kb].id1,binary[kb].id2);
          fprintf(stderr, "An isnan occured for r1 cmc_stellar_evolution.c\n");
          fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g tb=%g\n", binary[kb].bse_tphys, tphysf, binary[kb].bse_kw[0], binary[kb].bse_kw[1], binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_radius[0], binary[kb].bse_radius[1], binary[kb].bse_lum[0], binary[kb].bse_lum[1], binary[kb].bse_tb);
          fprintf(stderr, "k= %ld kb=%ld star_id=%ld bin_id1=%ld bin_id2=%ld \n", k, kb, star[k].id, binary[kb].id1, binary[kb].id2);
          exit(1);
        } 
        if(isnan(binary[kb].bse_radius[1])){
          fprintf(stderr, "An isnan occured for r2 cmc_stellar_evolution.c\n");
          fprintf(stderr, "tphys=%g tphysf=%g kstar1=%d kstar2=%d m1=%g m2=%g r1=%g r2=%g l1=%g l2=%g \n", binary[kb].bse_tphys, tphysf, binary[kb].bse_kw[0], binary[kb].bse_kw[1], binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_radius[0], binary[kb].bse_radius[1], binary[kb].bse_lum[0], binary[kb].bse_lum[1]);
          fprintf(stderr, "k= %ld kb=%ld star_id=%ld bin_id1=%ld bin_id2=%ld \n", k, kb, star[k].id, binary[kb].id1, binary[kb].id2);
          exit(1);
        }

	handle_bse_outcome(k, kb, vs, tphysf, kprev0, kprev1);

	if (WRITE_BH_INFO) {
		if (kprev0!=14 && binary[kb].bse_kw[0]==14) { // newly formed BH
			parafprintf(newbhfile, "%.18g %g 1 %ld %g %g %g %g %g", TotalTime, star_r[g_k], binary[kb].id1, binary[kb].bse_zams_mass[0], binary[kb].bse_mass0[0], binary[kb].bse_mass[0], binary[kb].bse_bhspin[0], VKO);
			for (ii=0; ii<16; ii++){
				parafprintf (newbhfile, " %g", vs[ii]);
			}
			parafprintf (newbhfile, "\n");
		}
		if (kprev1!=14 && binary[kb].bse_kw[1]==14 && binary[kb].id2 != 0) { // newly formed BH
			parafprintf(newbhfile, "%.18g %g 1 %ld %g %g %g %g %g", TotalTime, star_r[g_k], binary[kb].id2, binary[kb].bse_zams_mass[1],binary[kb].bse_mass0[1], binary[kb].bse_mass[1], binary[kb].bse_bhspin[1],VKO);
			for (ii=0; ii<16; ii++){
				parafprintf (newbhfile, " %g", vs[ii]);
			}
			parafprintf (newbhfile, "\n");
		}
	}
	//handle_bse_outcome(k, kb, vs, tphysf, kprev0, kprev1);
      }
    }
    bh_count(k);
  }

  double tmpTimeStart = timeStartSimple();
  double temp = 0.0;

  MPI_Status stat;
  
  //MPI: If we use the MPI_Reduce call, we have no control over the order in which the summing is done, so, we avoid that by using send/recv calls and sum up in order to avoid round-off errors and comparison with serial version.
  if(myid!=0)
    MPI_Send(&DMse, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  else
    for(i=1;i<procs;i++)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
      DMse += temp;
    }
	timeEndSimple(tmpTimeStart, &t_comm);
}

/**
* @brief ?
*/
void write_stellar_data(void){
  long k, kb;
  char filename[1024];
  
  se_file_counter++;
  
  /* single star info */
  sprintf(filename, "%s_stellar_info.%05d.dat", outprefix, se_file_counter);

  MPI_File mpi_stel_file;
  char mpi_stel_file_buf[10000], mpi_stel_file_wrbuf[10000000];
  long long mpi_stel_file_len=0, mpi_stel_file_ofst_total=0;
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_stel_file);
  MPI_File_set_size(mpi_stel_file, 0);

  pararootfprintf(stel_file, "# time (Myr): %e\n",
    TotalTime/MEGA_YEAR);
  pararootfprintf(stel_file, "# time (FP):  %e\n", TotalTime);
  pararootfprintf(stel_file,
	  "#  id        mass        radius     luminosity  type      ospin         B           formation     bacc      tacc  \n");
  pararootfprintf(stel_file,
	  "#======= ============ ============ ============ ====    ==========  ===========    ===========  ========  ========\n");

  for(k=1; k<=clus.N_MAX_NEW; k++){
    parafprintf(stel_file, "%08d ", get_global_idx(k));
    parafprintf(stel_file, "%e ", star[k].se_mt);
    parafprintf(stel_file, "%e ", star[k].se_radius);
    parafprintf(stel_file, "%e ", star[k].se_lum);
    parafprintf(stel_file, "%d ", star[k].se_k);
    parafprintf(stel_file, "%e ", star[k].se_ospin);
    parafprintf(stel_file, "%e ", star[k].se_scm_B);
    parafprintf(stel_file, "%g ", star[k].se_scm_formation);
    parafprintf(stel_file, "%g ", star[k].se_bacc);
    parafprintf(stel_file, "%g ", star[k].se_tacc);
    parafprintf(stel_file, "\n");
    parafprintf(stel_file, "%08d ", get_global_idx(k));
  }

  mpi_para_file_write(mpi_stel_file_wrbuf, &mpi_stel_file_len, &mpi_stel_file_ofst_total, &mpi_stel_file);
  MPI_File_close(&mpi_stel_file);

  /* binary star info */
  sprintf(filename, "%s_binary_stellar_info.%05d.dat", outprefix, se_file_counter);

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_stel_file);
  MPI_File_set_size(mpi_stel_file, 0);

  pararootfprintf(stel_file, "# time (Myr): %e\n",
    TotalTime/MEGA_YEAR);
  pararootfprintf(stel_file, "# time (FP):  %e\n", TotalTime);
  pararootfprintf(stel_file, "#1:id1 #2:id2 #3:M1[MSUN] #4:M2 #5:R1[RSUN] #6:R2 #7:k1 #8:k2 #9:Porb[day] #10:e #11:L1[LSUN] #12:L2 #13:Mcore1[MSUN] #14:Mcore2 #15:Rcore1[RSUN] #16:Rcore2 #17:Menv1[MSUN] #18:Menv2 #19:Renv1[RSUN] #20:Renv2 #21:Tms1[MYR] #22:Tms2 #23:bhspin1 #24:bhspin2 #25:Mdot1[MSUN/YR] #26:Mdot2 #27:R1/ROL1 #28:R2/ROL2 #29:ospin1 #30:ospin2 #31:B1 #32:B2 #33:formation1 #34:formation2 #35:bacc1 #36:bacc2 #37:tacc1 #38:tacc2\n");

  for(k=1; k<=clus.N_MAX_NEW; k++){
    if (star[k].binind) {
      kb = star[k].binind;
      parafprintf(stel_file, "%08ld %08ld %g %g %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        binary[kb].id1, binary[kb].id2,
        binary[kb].bse_mass[0], binary[kb].bse_mass[1],
        binary[kb].bse_radius[0], binary[kb].bse_radius[1],
        binary[kb].bse_kw[0], binary[kb].bse_kw[1],
        binary[kb].bse_tb, binary[kb].e,
        binary[kb].bse_lum[0], binary[kb].bse_lum[1],
        binary[kb].bse_massc[0], binary[kb].bse_massc[1],
        binary[kb].bse_radc[0], binary[kb].bse_radc[1],
        binary[kb].bse_menv[0], binary[kb].bse_menv[1],
        binary[kb].bse_renv[0], binary[kb].bse_renv[1],
        binary[kb].bse_tms[0], binary[kb].bse_tms[1],
        binary[kb].bse_bhspin[0], binary[kb].bse_bhspin[1],
        binary[kb].bse_bcm_dmdt[0], binary[kb].bse_bcm_dmdt[1],
        binary[kb].bse_bcm_radrol[0], binary[kb].bse_bcm_radrol[1],
        binary[kb].bse_ospin[0], binary[kb].bse_ospin[1],
        binary[kb].bse_bcm_B[0], binary[kb].bse_bcm_B[1],
        binary[kb].bse_bcm_formation[0], binary[kb].bse_bcm_formation[1],
        binary[kb].bse_bacc[0], binary[kb].bse_bacc[1],
        binary[kb].bse_tacc[0], binary[kb].bse_tacc[1]);
    }
  }
  mpi_para_file_write(mpi_stel_file_wrbuf, &mpi_stel_file_len, &mpi_stel_file_ofst_total, &mpi_stel_file);
  MPI_File_close(&mpi_stel_file);
}

/**
* @brief ?
*
* @param k index of star 1
* @param kb index of star 2
* @param vs ?
* @param tphysf ?
*/
void handle_bse_outcome(long k, long kb, double *vs, double tphysf, int kprev0, int kprev1)
{
  int j, jj;
  long knew=0, knewp=0, convert;
  double dtp, VKO;
  
  knew = 0;
  VKO = 0.0;

  /* PK: extract some stellar/binary info from BSE's bcm array */
  /* but don't loop forever, just until maximum array length.  */
  /* 5100000 will eventually be a possible max array length... */
  j = 1;
  jj = 1;
  while (bse_get_bcm(j, 1) >= 0 && j<50000) {
    j++;
  }
  j--;
  if(j+1 > 50000) {
     j = 0;
  }
  if (j >= 1) {
   if ((fabs((binary[kb].bse_tphys - bse_get_bcm(j,1))/binary[kb].bse_tphys) >= 1.0e-6) && !(kprev0 == 14 && kprev1 == 14)) {
        wprintf("binary[kb].bse_tphys=%g bcmtime=%g\n", binary[kb].bse_tphys, bse_get_bcm(j,1));
        /* exit_cleanly(-1); */
    }
    binary[kb].bse_bcm_dmdt[0] = bse_get_bcm(j, 14);
    binary[kb].bse_bcm_dmdt[1] = bse_get_bcm(j, 28);
    binary[kb].bse_bcm_radrol[0] = bse_get_bcm(j, 15);
    binary[kb].bse_bcm_radrol[1] = bse_get_bcm(j, 29);
    binary[kb].bse_bcm_B[0] = bse_get_bcm(j, 33); /* PK (note, minus one could be here - is just being paranoid about bcm arrays output, but sometimes you might need to go back a bit further to get an accurate value... Of course this could be a problem if the pulsar was just formed and one step back there was no magnetic field! */
    binary[kb].bse_bcm_B[1] = bse_get_bcm(j, 34);
//    binary[kb].bse_bcm_formation[0] = bse_get_bcm(j, 35);
//    binary[kb].bse_bcm_formation[1] = bse_get_bcm(j, 36);
// Check if merger of a NS occured. If placed into the other star's position this will cause formation to be updated as a zero.
    convert = 0;
    if((binary[kb].bse_mass[0] == 0.0 || binary[kb].bse_mass[1] == 0.0) && (binary[kb].bse_bcm_formation[0]!=0 || binary[kb].bse_bcm_formation[1]!=0)){
        if(binary[kb].bse_mass[0]!=0.0 && binary[kb].bse_bcm_formation[1]!=0 && binary[kb].bse_bcm_formation[0]==0){
                convert = 1;
        } else if (binary[kb].bse_mass[1]!=0.0 && binary[kb].bse_bcm_formation[0]!=0 && binary[kb].bse_bcm_formation[1]==0){
                convert = 2;
        }
    }
    while (bse_get_bcm(jj,1)>=0.0) {
      if(jj > 1){
        if(bse_get_bcm(jj,2) == 13 && bse_get_bcm(jj-1,2) < 13){
           if(bse_get_bcm(jj+1,1) >= 0.0) {
              binary[kb].bse_bcm_formation[0] = bse_get_bcm(jj+1,39); //again just being paranoid about bcm array, but here moving forward by one.
//              printf("stel, pulsar bin. 1 : %g %g %g %g \n",binary[kb].bse_B_0[0],binary[kb].bse_ospin[0],bse_get_bcm(jj+1,33),binary[kb].bse_bcm_formation[0]);
           } else {
              binary[kb].bse_bcm_formation[0] = bse_get_bcm(jj,39);
//              printf("stel, pulsar bin. 2 : %g %g %g %g \n",binary[kb].bse_B_0[0],binary[kb].bse_ospin[0],bse_get_bcm(jj,33),binary[kb].bse_bcm_formation[0]);
           }
        }
        if(bse_get_bcm(jj,16) == 13 && bse_get_bcm(jj-1,16) < 13){
           if(bse_get_bcm(jj+1,1) >= 0.0) {
              binary[kb].bse_bcm_formation[1] = bse_get_bcm(jj+1,40);
//              printf("stel, pulsar bin. 3 : %g %g %g %g \n",binary[kb].bse_B_0[1],binary[kb].bse_ospin[1],bse_get_bcm(jj+1,34),binary[kb].bse_bcm_formation[1]);
           } else {
              binary[kb].bse_bcm_formation[1] = bse_get_bcm(jj,40);
//              printf("stel, pulsar bin. 4 : %g %g %g %g \n",binary[kb].bse_B_0[1],binary[kb].bse_ospin[1],bse_get_bcm(jj,34),binary[kb].bse_bcm_formation[1]);
           }
        }
      }
      jj++;
    }
// Check convert to see if formation was updated incorrectly. If so then correct it.
    if(convert>0){
        if(convert==1){
                binary[kb].bse_bcm_formation[0] = binary[kb].bse_bcm_formation[1];
        } else if (convert==2){
                binary[kb].bse_bcm_formation[1] = binary[kb].bse_bcm_formation[0];
        }
    }
  } else {
    eprintf("Could not extract BSE bcm info!  Input dtp not exactly equal to tphysf-tphys?");
    exit_cleanly(-1, __FUNCTION__);
  }

  if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0 && binary[kb].bse_tb > 0.0) {
    /* normal evolution */
    binary[kb].rad1 = binary[kb].bse_radius[0] * RSUN / units.l;
    binary[kb].rad2 = binary[kb].bse_radius[1] * RSUN / units.l;
    binary[kb].m1 = binary[kb].bse_mass[0] * MSUN / units.mstar;
    binary[kb].m2 = binary[kb].bse_mass[1] * MSUN / units.mstar;
    int g_k = get_global_idx(k);
    star_m[g_k] = binary[kb].m1 + binary[kb].m2;
    DMse -= star_m[g_k] * madhoc;
    binary[kb].a = pow((binary[kb].bse_mass[0]+binary[kb].bse_mass[1])*sqr(binary[kb].bse_tb/365.25), 1.0/3.0)
      * AU / units.l;
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    if (vs[0] <= 0.0) {
        star[k].vr += 0.0;
        star[k].vt += 0.0;
        set_star_EJ(k);
  VKO = 0.0;
    } else if (vs[4] <= 0.0) {
    /* If one kick */
       star[k].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
       //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(k);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else {
    /* Two kicks */
       star[k].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[k].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[k].vt),vs[5],vs[6], curr_st);//will mean a different angle for each kick, should be ok as this is a randomised process...
       //star[k].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(k);
       VKO = sqrt(vs[1]*vs[1] + vs[2]*vs[2] + vs[3]*vs[3]) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
  dprintf("birth kick(bin): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g VK=%.18g id1=%ld id2=%ld pid1=%ld pid2=%ld type1=%d type2=%d \n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,binary[kb].id1,binary[kb].id2,binary[kb].id1,binary[kb].id2, binary[kb].bse_kw[0], binary[kb].bse_kw[1]);
    }
    /* extract some binary info from BSE's bcm array */
//    j = 1;
//    while (bse_get_bcm(j, 1) >= 0.0) {
//      j++;
//    }
//    j--;
//    if (j >= 1) {
//      if (fabs((binary[kb].bse_tphys - bse_get_bcm(j,1))/binary[kb].bse_tphys) >= 1.0e-6) {
//	wprintf("binary[kb].bse_tphys=%g bcmtime=%g\n", binary[kb].bse_tphys, bse_get_bcm(j,1));
//	/* exit_cleanly(-1); */
//      }
//      binary[kb].bse_bcm_dmdt[0] = bse_get_bcm(j, 14);
//      binary[kb].bse_bcm_dmdt[1] = bse_get_bcm(j, 28);
//      binary[kb].bse_bcm_radrol[0] = bse_get_bcm(j, 15);
//      binary[kb].bse_bcm_radrol[1] = bse_get_bcm(j, 29);
//    } else {
//      eprintf("Could not extract BSE bcm info!  Input dtp not exactly equal to tphysf-tphys?");
//      exit_cleanly(-1, __FUNCTION__);
//    }
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] != 0.0) {
    /* disruption with both stars "intact" */
    //dprintf("binary disrupted via BSE with both stars intact\n");
    knew = create_star(k, 1);
    knewp = create_star(k, 1);
    cp_binmemb_to_star(k, 0, knew);
    cp_binmemb_to_star(k, 1, knewp);
    DMse -= (star_m[get_global_idx(knew)] + star_m[get_global_idx(knewp)]) * madhoc;

    parafprintf(semergedisruptfile, "t=%g disruptboth id1=%ld(m1=%g) id2=%ld(m2=%g) (r=%g) type1=%d type2=%d\n",
      TotalTime,
      star[knew].id, star[knew].se_mt, 
      star[knewp].id, star_m[get_global_idx(knewp)] * units.mstar / FB_CONST_MSUN,
      star_r[get_global_idx(k)], kprev0, kprev1);

    destroy_obj(k);
    /* in this case vs is relative speed between stars at infinity */
/*    star[knew].vr += star[knewp].m/(star[knew].m+star[knewp].m) * vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += star[knewp].m/(star[knew].m+star[knewp].m) * sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    star[knewp].vr += -star[knew].m/(star[knew].m+star[knewp].m) * vs[2] * 1.0e5 / (units.l/units.t);
    star[knewp].vt += -star[knew].m/(star[knew].m+star[knewp].m) * sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knewp);
*/
    /* vs is now the recoil speed of both runaway stars owing to SN disruption */
    if (vs[4]>0.0 && vs[8]<=0.0 ) {
      /* 1 kick occured and it disrupted the system */
        if (vs[0]==1) {
      /* Star knew went SN */
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);

            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            set_star_EJ(knew);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knewp].vr += vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);
            //star[knewp].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t); 
            set_star_EJ(knewp);
        } else {
       /* Star knewp went SN */
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knew].vr += vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            //star[knew].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t); //minus at front?
            set_star_EJ(knew);
        }
    } else if ((vs[4]>0.0 && vs[8]>0.0) && (vs[4] == vs[8])) {
      /* Two SNe and the 2nd one disrupts the system. 
         Thus the primary receives vs[1-3] and vs[9-11] secondary receives vs[1-3] and vs[5-7] */
        if (vs[0]==1) {
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            star[knewp].vt = star[knew].vt; //update, now seperate companion tangent velocity, after first SN when still part of binary here.
            vt_add_kick(&(star[knew].vt),vs[9],vs[10], curr_st);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[9]*vs[9]+vs[10]*vs[10]+vs[11]*vs[11]);
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);//finalise companion vt...

            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        } else {
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            star[knew].vt = star[knewp].vt;
            vt_add_kick(&(star[knewp].vt),vs[9],vs[10], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[9]*vs[9]+vs[10]*vs[10]+vs[11]*vs[11]);
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
        }
    } else if ((vs[4]>0.0 && vs[8]>0.0) && (vs[4] != vs[8])) {
      /* Two SNe and the 1st one disrupts the system.
         Primary feels vs[1-3] and secondary feels vs[5-7] and vs[9-11]. */
        if (vs[0]==1) {
            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
            set_star_EJ(knew);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knewp].vr += vs[7] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[5],vs[6], curr_st);
            vt_add_kick(&(star[knewp].vt),vs[9],vs[10], curr_st);
            //star[knewp].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        } else {
            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
            //star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knewp);
        VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
            star[knew].vr += vs[7] * 1.0e5 / (units.l/units.t) + vs[11] * 1.0e5 / (units.l/units.t);
            vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
            vt_add_kick(&(star[knew].vt),vs[9],vs[10], curr_st);
            //star[knew].vt += sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t) + sqrt(vs[9]*vs[9]+vs[10]*vs[10]) * 1.0e5 / (units.l/units.t);
            set_star_EJ(knew);
        }
    } else {
      /* No kick */
        star[knew].vr += 0.0;
        star[knew].vt += 0.0;
        set_star_EJ(knew);
        VKO = 0.0;
        star[knewp].vr += 0.0;
        star[knewp].vt += 0.0;
        set_star_EJ(knewp);
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
      dprintf("birth kick(iso1.2): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g id1=%ld id2=%ld pid1=%ld pid2=%ld type1=%d type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k,star[knewp].se_k);
    }
  } else if (binary[kb].bse_mass[0] != 0.0 && binary[kb].bse_mass[1] == 0.0) {
    /* secondary star gone */
    knew = create_star(k, 1);
    cp_binmemb_to_star(k, 0, knew);

	/*If this was a BBH merger, special things must be done*/
	if(kprev0 == 14 && kprev1 == 14)
		binary_bh_merger(k, kb, knew, kprev0, kprev1, curr_st);

    parafprintf(semergedisruptfile, "t=%g disrupt1 idr=%ld(mr=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typer=%d type1=%d type2=%d\n",
      TotalTime,
      star[knew].id, star[knew].se_mt,
      binary[kb].id1, binary[kb].m1 * units.mstar / FB_CONST_MSUN,
      binary[kb].id2, binary[kb].m2 * units.mstar / FB_CONST_MSUN,
      star_r[get_global_idx(k)], star[knew].se_k, kprev0, kprev1);
    destroy_obj(k);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }

/*    star[knew].vr += vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew); */
    /* Update velocity owing to kick(s?) - could be overkill here... */
    if (vs[0]>0.0 && vs[4]<=0.0 ) {
       /* One SN occured */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]==vs[4])) {
       /* 1 SN occured disrupted system then one star killed itself */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]!=vs[4])) {
       /* 2 SNe then system mergers, star feels both kicks (i.e. it is COM) */
       /* NOTE: merger remnants of double compact (NS or BH) binaries do not receive kicks in BSE as of yet!
                When it is introduced may have to put new else if statement in here collecting all (different) vs[]'s. */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    } else {
      /* No kick - going to set_star_EJ seems overkill (here and elsewhere). */
       star[knew].vr += 0.0;
       star[knew].vt += 0.0;
       set_star_EJ(knew);
       VKO = 0.0;
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
  dprintf("birth kick(iso2): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g id1=%ld id2=na pid1=%ld pid2=%ld type1=%d type2=15\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k);
    }
    /* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
    dtp = tphysf - star[knew].se_tphys;
    dtp = 0.0;
  /* Update star id for pass through. */
    bse_set_id1_pass(star[knew].id);
    bse_set_id2_pass(0);
    /* HAVE REMOVED SAFE SINGLE EVOLVE BELIEVING THAT NOW SELF-CONSISTENT MERGERS SHOULD ALWAYS RESULT
    bse_evolv1_safely(&(star[knew].se_k), &(star[knew].se_mass), &(star[knew].se_mt), &(star[knew].se_radius), 
          &(star[knew].se_lum), &(star[knew].se_mc), &(star[knew].se_rc), &(star[knew].se_menv), 
          &(star[knew].se_renv), &(star[knew].se_ospin), &(star[knew].se_epoch), &(star[knew].se_tms), 
          &(star[knew].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
    */
    
    star[knew].rad = star[knew].se_radius * RSUN / units.l;
    star_m[get_global_idx(knew)] = star[knew].se_mt * MSUN / units.mstar;
    DMse -= star_m[get_global_idx(knew)] * madhoc;

    /* birth kicks */
    /* ALSO REMOVE BELOW AS WE NOW WONT PRODUCE ANOTHER KICK
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
    vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
    //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
  dprintf("birth kick(iso3): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VKO=%.18g, id=%ld id1=%ld id2=%ld type1=%d type2=15\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].id,binary[kb].id1,binary[kb].id2,star[knew].se_k);
    }
    */
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] != 0.0) {
    /* primary star gone */
    //dprintf("binary disrupted via BSE with second star intact\n");
    knew = create_star(k, 1);
    cp_binmemb_to_star(k, 1, knew);

	/*If this was a BBH merger, special things must be done*/
	if(kprev0 == 14 && kprev1 == 14)
		binary_bh_merger(k, kb, knew, kprev0, kprev1, curr_st);

    parafprintf(semergedisruptfile, "t=%g disrupt2 idr=%ld(mr=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g) typer=%d type1=%d type2=%d\n",
      TotalTime,
      star[knew].id, star[knew].se_mt,
      binary[kb].id1, binary[kb].m1 * units.mstar / FB_CONST_MSUN,
      binary[kb].id2, binary[kb].m2 * units.mstar / FB_CONST_MSUN,
      star_r[get_global_idx(k)], star[knew].se_k, kprev0, kprev1);

    destroy_obj(k);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
/*    star[knew].vr += vs[2] * 1.0e5 / (units.l/units.t);
    star[knew].vt += sqrt(vs[0]*vs[0]+vs[1]*vs[1]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew); */
    /* Update velocity owing to kick(s?) - could be overkill here... */
    if (vs[0]>0.0 && vs[4]<=0.0 ) {
       /* One SN occured */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]==vs[4])) {
       /* 1 SN occured disrupted system then one star killed itself */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    } else if (vs[0]>0.0 && vs[4]>0.0 && vs[8]<=0.0 && (vs[0]!=vs[4])) {
       /* 2 SNe then system mergers, star feels both kicks (i.e. it is COM) */
       /* NOTE: merger remnants of double compact (NS or BH) binaries do not receive kicks in BSE as of yet!
                When it is introduced may have to put new else if statement in here collecting all (different) vs[]'s. */
       star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t) + vs[7] * 1.0e5 / (units.l/units.t);
       vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
       vt_add_kick(&(star[knew].vt),vs[5],vs[6], curr_st);
       //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t) + sqrt(vs[5]*vs[5]+vs[6]*vs[6]) * 1.0e5 / (units.l/units.t);
       set_star_EJ(knew);
       VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]+vs[5]*vs[5]+vs[6]*vs[6]+vs[7]*vs[7]);
    } else {
      /* No kick - going to set_star_EJ seems overkill (here and elsewhere). */
       star[knew].vr += 0.0;
       star[knew].vt += 0.0;
       set_star_EJ(knew);
       VKO = 0.0;
    }
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
  dprintf("birth kick(iso4): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VKO=%.18g id=%ld id1=%ld id2=%ld type1=15 type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO, star[knew].id, binary[kb].id1, binary[kb].id2, star[knew].se_k);
    }
    
    /* here we do a safe single evolve, just in case the remaining star is a non self-consistent merger */
    dtp = tphysf - star[knew].se_tphys;
    dtp = 0.0;
  /* Update star id for pass through. */
    /* REMOVED AS PART OF EVOLV1 PHASE OUT
    bse_set_id1_pass(star[knew].id);
    bse_set_id2_pass(0);
    bse_evolv1_safely(&(star[knew].se_k), &(star[knew].se_mass), &(star[knew].se_mt), &(star[knew].se_radius), 
          &(star[knew].se_lum), &(star[knew].se_mc), &(star[knew].se_rc), &(star[knew].se_menv), 
          &(star[knew].se_renv), &(star[knew].se_ospin), &(star[knew].se_epoch), &(star[knew].se_tms), 
          &(star[knew].se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
    */
    
    star[knew].rad = star[knew].se_radius * RSUN / units.l;
    star_m[get_global_idx(knew)] = star[knew].se_mt * MSUN / units.mstar;
    DMse -= star_m[get_global_idx(knew)] * madhoc;

    /* birth kicks */
    /* REMOVED AGAIN NOT TO ADD KICK AGAIN
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
    }
    star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
    vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
    //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
    set_star_EJ(knew);
    VKO = sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]);
    if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
      //dprintf("birth kick of %f km/s\n", sqrt(vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2]));
  dprintf("birth kick(iso6): TT=%.18g, vs[0]=%.18g, vs[1]=%.18g, vs[2]=%.18g, vs[3]=%.18g, vr=%.18g, vt=%.18g, VK=%.18g type1=15 type2=%d\n",TotalTime,vs[0],vs[1],vs[2],vs[3],star[k].vr,star[k].vt,VKO,star[knew].se_k);
    }
    */
  } else if (binary[kb].bse_mass[0] == 0.0 && binary[kb].bse_mass[1] == 0.0) {
    /* both stars gone */
    //dprintf("binary disrupted via BSE with no stars intact\n");
    destroy_obj(k);
  } else {
    dprintf("unhandled binary outcome!\n");
  /* End by looking for boom stuff to log and log pulsar info */

	 if(WRITE_PULSAR_INFO)
	 {
		 if (knew) {
			 //bcm_boom_search(knew, vs, getCMCvalues);
			 pulsar_write(knew, VKO);
		 } else {
			 //bcm_boom_search(k, vs, getCMCvalues);
			 pulsar_write(k, VKO);
		 }
	 }

	 dprintf("bse_mass0=%g bse_mass1=%g tb=%g\n",
      binary[kb].bse_mass[0], binary[kb].bse_mass[1], binary[kb].bse_tb);
    exit_cleanly(-1, __FUNCTION__);
  }

	/*Shi*/
	if(tcount%PULSAR_DELTACOUNT==0){
		if (WRITE_MOREPULSAR_INFO){
        		if (knew){
                		write_morepulsar(knew);

                        	if (knewp){
                        		write_morepulsar(knewp);
                        	}

                	} else {
                        	write_morepulsar(k);

                	}	
        	}	
	}
}

/**
* @brief Output info of pulsars
*
* @param k star index
* @param kick ?
*/
void pulsar_write(long k, double kick)
{
	long b;
	double phi_r0, phi_rt, spin, twopi=6.283185307179586, yearsc=31557600;
	int g_k;
	double r, phi;

	g_k = get_global_idx(k);
	r = star_r[g_k];
	phi = star_phi[g_k];

	phi_r0 = potential(0.0);
	phi_rt = potential(Rtidal);
	b = star[k].binind;
	if(b>0) {
		//BINARIES
		if(binary[b].bse_kw[0]==13){
			spin = (twopi*yearsc)/binary[b].bse_ospin[0];
	  		parafprintf(pulsarfile,"%ld %.8g %ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %d %.8g %g %g %g %g %g %g %g %g %g %g %g\n", tcount, TotalTime, star[k].id, star[k].r_peri, star[k].r_apo, r, star[k].vr, star[k].vt, phi, phi_r0, phi_rt, kick, binary[b].id1, binary[b].id2, binary[b].bse_kw[1], spin, binary[b].bse_bcm_B[0], binary[b].bse_bcm_formation[0], binary[b].bse_bacc[0], binary[b].bse_tacc[0], binary[b].bse_B_0[0], binary[b].bse_tb, binary[b].bse_mass[1], binary[b].bse_mass[0], binary[b].e, binary[b].bse_bcm_radrol[1], binary[b].bse_bcm_dmdt[0]);
		}
		if(binary[b].bse_kw[1]==13){
		        spin = (twopi*yearsc)/binary[b].bse_ospin[0];
	  		parafprintf(pulsarfile,"%ld %.8g %ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %d %.8g %g %g %g %g %g %g %g %g %g %g %g\n", tcount, TotalTime, star[k].id, star[k].r_peri, star[k].r_apo, r, star[k].vr, star[k].vt, phi, phi_r0, phi_rt, kick, binary[b].id2, binary[b].id1, binary[b].bse_kw[0], spin, binary[b].bse_bcm_B[1], binary[b].bse_bcm_formation[1], binary[b].bse_bacc[1], binary[b].bse_tacc[1], binary[b].bse_B_0[1], binary[b].bse_tb, binary[b].bse_mass[0], binary[b].bse_mass[1], binary[b].e, binary[b].bse_bcm_radrol[0], binary[b].bse_bcm_dmdt[1]);

		}
	} else {
		//ISOLATED
		if(star[k].se_k==13){
		  spin = (twopi*yearsc)/star[k].se_ospin;
		  parafprintf(pulsarfile,"%ld %.8g %ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g  0  0  0 %.8g %g %g %g %g %g  0.0  0.0 %g  0.0  0.0  0.0\n", tcount, TotalTime, star[k].id, star[k].r_peri, star[k].r_apo, r, star[k].vr, star[k].vt, phi, phi_r0, phi_rt, kick, spin, star[k].se_scm_B, star[k].se_scm_formation, star[k].se_bacc, star[k].se_tacc, star[k].se_B_0, star[k].se_mt);

		}
	}
}

/**
 * @brief Output info of morepulsars
 * 
 * @param k star index
 */
void write_morepulsar(long i){      //Shi

        long j;
        double spin, spin0, spin1, twopi=6.283185307179586, yearsc=31557600;
	int g_i;
        double r;

        g_i = get_global_idx(i);
        r = star_r[g_i];

        j=star[i].binind;

        if (j==0){ //Single
                if (star[i].se_k==13){
                        spin = (twopi*yearsc)/star[i].se_ospin;
                        parafprintf(morepulsarfile, "%ld %.8g 0 %ld -100 %.8g -100 %g -100 %g -100 %d -100 -100 -100 -100 -100 -100 -100 %.8g %.8g %.8g -100 -100 -100 -100\n", tcount, TotalTime, star[i].id, star[i].se_mt, star[i].se_scm_B, spin, star[i].se_k, r, star[i].vr, star[i].vt);

                }
        } else { //Binary
                if (binary[j].bse_kw[0]==13 || binary[j].bse_kw[1]==13){
                        spin0 = (twopi*yearsc)/binary[j].bse_ospin[0];
                        spin1 = (twopi*yearsc)/binary[j].bse_ospin[1];
                        parafprintf(morepulsarfile, "%ld %.8g 1 %ld %ld %.8g %.8g %g %g %g %g %d %d %.8g %.8g %g %g %g %g %.8g %.8g %.8g %g %g %g %g\n", tcount, TotalTime, binary[j].id1, binary[j].id2, binary[j].bse_mass[0], binary[j].bse_mass[1], binary[j].bse_bcm_B[0], binary[j].bse_bcm_B[1], spin0, spin1, binary[j].bse_kw[0], binary[j].bse_kw[1], binary[j].a* units.l/AU, binary[j].e, binary[j].bse_bcm_radrol[0], binary[j].bse_bcm_radrol[1], binary[j].bse_bcm_dmdt[0], binary[j].bse_bcm_dmdt[1], r, star[i].vr, star[i].vt, binary[j].bse_bacc[0], binary[j].bse_bacc[1], binary[j].bse_tacc[0], binary[j].bse_tacc[1]);

                }
        }
}

/* outputs boom information */
//void boomoutput(long int s_id, double kick, double remtype, int(*getCMCvalues)(long int s_id, *double TC, *double TT, *long ID, *int )

/**
* @brief ?
*
* @param s_id ?
* @param kick ?
* @param remtype1 ?
* @param remtype2 ?
* @param progtype1 ?
* @param progtype2 ?
* @param formation ?
* @param boomtype ?
* @param boomstar ?
*/
void getCMCvalues(long s_id, double kick, int remtype1, int remtype2, int progtype1, int progtype2, double formation, int boomtype, int boomstar)
{
// here kick should be total kick(s) received by that star in the timestep that the star exploded. It might be that two stars exploded in the one timestep. Then kick (for that star) is magnitude of both kicks (direct and indirect) that star received...
	double phi_r0, phi_rt;
	phi_r0 = potential(0.0);
	phi_rt = potential(Rtidal);
	if(boomtype==1){
// 	  fprintf(boomlog,"%s %ld %.8g %ld %d %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %g %d\n", "NS", tcount, TotalTime, star[s_id].id, star[s_id].se_k, remtype1, remtype2, progtype1, progtype2, star[s_id].r_peri, star[s_id].r_apo, star[s_id].r, star[s_id].vr, star[s_id].vt, kick, star[s_id].phi, phi_r0, phi_rt, binary[star[s_id].binind].id1, binary[star[s_id].binind].id2, formation, boomstar);
	} else if (boomtype==2){
//		fprintf(boomlog,"%s %ld %.8g %ld %d %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %g %d\n", "SDS_1a", tcount, TotalTime, star[s_id].id, star[s_id].se_k, remtype1, remtype2, progtype1, progtype2, star[s_id].r_peri, star[s_id].r_apo, star[s_id].r, star[s_id].vr, star[s_id].vt, kick, star[s_id].phi, phi_r0, phi_rt, binary[star[s_id].binind].id1, binary[star[s_id].binind].id2, formation, boomstar);
	} else if (boomtype==3){
//		fprintf(boomlog,"%s %ld %.8g %ld %d %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %g %d\n", "DDS_1a", tcount, TotalTime, star[s_id].id, star[s_id].se_k, remtype1, remtype2, progtype1, progtype2, star[s_id].r_peri, star[s_id].r_apo, star[s_id].r, star[s_id].vr, star[s_id].vt, kick, star[s_id].phi, phi_r0, phi_rt, binary[star[s_id].binind].id1, binary[star[s_id].binind].id2, formation, boomstar);
	} else if (boomtype==4){
//		fprintf(boomlog,"%s %ld %.8g %ld %d %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %g %d\n", "BH", tcount, TotalTime, star[s_id].id, star[s_id].se_k, remtype1, remtype2, progtype1, progtype2, star[s_id].r_peri, star[s_id].r_apo, star[s_id].r, star[s_id].vr, star[s_id].vt, kick, star[s_id].phi, phi_r0, phi_rt, binary[star[s_id].binind].id1, binary[star[s_id].binind].id2, formation, boomstar);
	} else {
//		fprintf(boomlog,"%s %ld %.8g %ld %d %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %ld %ld %g %d\n", "DBUG", tcount, TotalTime, star[s_id].id, star[s_id].se_k, remtype1, remtype2, progtype1, progtype2, star[s_id].r_peri, star[s_id].r_apo, star[s_id].r, star[s_id].vr, star[s_id].vt, kick, star[s_id].phi, phi_r0, phi_rt, binary[star[s_id].binind].id1, binary[star[s_id].binind].id2, formation, boomstar);
        }
}
/**
* @brief Meagan: count different types of bh-objects; at end of timestep, we'll print these totals
*
* @param k index of star
*/
void bh_count(long k) {
	long b;
	double phi_r0, phi_rt;
	phi_r0 = potential(0.0);
	phi_rt = potential(Rtidal);
	b = star[k].binind;
	if(b>0) {
                //BINARIES
		if(binary[b].bse_kw[0]==14 || binary[b].bse_kw[1]==14) { // binary containing at least one bh
			bhbinary += 1;
			if(binary[b].bse_kw[0]==14 && binary[b].bse_kw[1]==14) { // both are bhs
				bhbh += 1;
			} else { // one of the stars is not a BH
				bhnonbh += 1;
				//if (binary[b].bse_kw[0]==14 && binary[b].bse_kw[1]!=14) { // primary bh, secondary other
				if (binary[b].bse_kw[0]==13 || binary[b].bse_kw[1]==13) { //  BH-NS
					bh13 += 1;			
				} else if (binary[b].bse_kw[0]==10 || binary[b].bse_kw[1]==10) {
					bh10 += 1;
				} else if (binary[b].bse_kw[0]==11 || binary[b].bse_kw[1]==11) {
					bh11 += 1; 
				} else if (binary[b].bse_kw[0]==12 || binary[b].bse_kw[1]==12) { 
					bh12 += 1;
				} else if (binary[b].bse_kw[0]==9 || binary[b].bse_kw[1]==9 ||
					   binary[b].bse_kw[0]==8 || binary[b].bse_kw[1]==8) { 
					bh89 += 1;  // BH-POSTMS_HE binary
				} else if (binary[b].bse_kw[0]==7 || binary[b].bse_kw[1]==7) { 
					bh7 += 1;  // BH-MS_HE binary
				} else if ( (binary[b].bse_kw[0]>=2 && binary[b].bse_kw[0]<=6) ||
					    (binary[b].bse_kw[1]>=2 && binary[b].bse_kw[1]<=6) ) {
					bh26 += 1;  // BH-postMS binary
				} else if ((binary[b].bse_kw[0]>=0 && binary[b].bse_kw[0]<=1) ||
					   (binary[b].bse_kw[1]>=0 && binary[b].bse_kw[1]<=1)) {
					bh01 += 1;  // BH-MS binary
				}
				bhstar = bh01 + bh26 + bh7 + bh89;
				bhwd = bh10 + bh11 + bh12;
			}	

		}

	} else {
		//SINGLE BHs
		if(star[k].se_k==14){
			bhsingle += 1;
		}
    }
}

/**
* @brief copies corresponding binary member variables of given star to the star variables of new star
*
* @param k star index
* @param kbi 0 or 1 based on which binary member
* @param knew index of new star
*/
void cp_binmemb_to_star(long k, int kbi, long knew)
{
  long kb;
  
  kb = star[k].binind;
  long g_k = get_global_idx(k);
  long g_knew = get_global_idx(knew);
  star_r[g_knew] = star_r[g_k];
  if(kbi == 0){
    star_m[g_knew] = binary[kb].m1;//bse_mass[kbi] * MSUN / units.mstar;
  }
  else{
    star_m[g_knew] = binary[kb].m2;//bse_mass[kbi] * MSUN / units.mstar;
  }
  star_phi[g_knew] = star_phi[g_k];
  star[knew].vr = star[k].vr;
  star[knew].vt = star[k].vt;
  set_star_EJ(knew);
  set_star_news(knew);
  set_star_olds(knew);
  /* mark stars as interacted so they don't undergo E_CONS mode stuff */
  star[knew].interacted = 1;
  if (kbi == 0) {
    star[knew].Eint = binary[kb].Eint1;
    star[knew].id = binary[kb].id1;
  } else {
    star[knew].Eint = binary[kb].Eint2;
    star[knew].id = binary[kb].id2;
  }
  star[knew].rad = binary[kb].bse_radius[kbi] * RSUN / units.l;
  star[knew].se_mass = binary[kb].bse_mass0[kbi]; /* initial mass (at curent epoch?) */
  star[knew].se_k = binary[kb].bse_kw[kbi];
  star[knew].se_mt = binary[kb].bse_mass[kbi]; /* current mass */
  star[knew].se_ospin = binary[kb].bse_ospin[kbi];
  star[knew].se_B_0 = binary[kb].bse_B_0[kbi]; /* PK */
  star[knew].se_bacc = binary[kb].bse_bacc[kbi];
  star[knew].se_tacc = binary[kb].bse_tacc[kbi];
  star[knew].se_scm_B = binary[kb].bse_bcm_B[kbi];
  star[knew].se_scm_formation = binary[kb].bse_bcm_formation[kbi];
  star[knew].se_epoch = binary[kb].bse_epoch[kbi];
  star[knew].se_tphys = binary[kb].bse_tphys;
  star[knew].se_radius = binary[kb].bse_radius[kbi];
  star[knew].se_lum = binary[kb].bse_lum[kbi];
  star[knew].se_mc = binary[kb].bse_massc[kbi];
  star[knew].se_rc = binary[kb].bse_radc[kbi];
  star[knew].se_menv = binary[kb].bse_menv[kbi];
  star[knew].se_renv = binary[kb].bse_renv[kbi];
  star[knew].se_tms = binary[kb].bse_tms[kbi];
  star[knew].se_bhspin = binary[kb].bse_bhspin[kbi];
  //Sourav: toy rejuvenation- variables updating
  if (kbi==0){
    star[knew].createtime = binary[kb].createtime_m1;
    star[knew].lifetime = binary[kb].lifetime_m1;
  } else {
    star[knew].createtime = binary[kb].createtime_m2;
    star[knew].lifetime = binary[kb].lifetime_m2;
  }
}

/**
* @brief copies variables of old star to stellar evolution variables of new star
*
* @param oldk old star index
* @param kbi 0 or 1 for binary, -1 for non-binary
* @param knew index of new star
*/
void cp_SEvars_to_newstar(long oldk, int kbi, long knew)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    star[knew].se_mass = star[oldk].se_mass;
    star[knew].se_k = star[oldk].se_k;
    star[knew].se_mt = star[oldk].se_mt;
    star[knew].se_ospin = star[oldk].se_ospin;
    star[knew].se_B_0 = star[oldk].se_B_0; /* PK */
    star[knew].se_bacc = star[oldk].se_bacc;
    star[knew].se_tacc = star[oldk].se_tacc;
    star[knew].se_scm_B = star[oldk].se_scm_B;
    star[knew].se_scm_formation = star[oldk].se_scm_formation;
    star[knew].se_epoch = star[oldk].se_epoch;
    star[knew].se_tphys = star[oldk].se_tphys;
    star[knew].se_radius = star[oldk].se_radius;
    star[knew].se_lum = star[oldk].se_lum;
    star[knew].se_mc = star[oldk].se_mc;
    star[knew].se_rc = star[oldk].se_rc;
    star[knew].se_menv = star[oldk].se_menv;
    star[knew].se_renv = star[oldk].se_renv;
    star[knew].se_tms = star[oldk].se_tms;
	star[knew].se_bhspin = star[oldk].se_bhspin;
    //Sourav: toy rejuvenation- updating the createtime and lifetime
    star[knew].createtime = star[oldk].createtime;
    star[knew].lifetime = star[oldk].lifetime;
    star[knew].rad = star[oldk].rad;
  } else { /* star comes from input binary */
    star[knew].se_mass = binary[kb].bse_mass0[kbi];
    star[knew].se_k = binary[kb].bse_kw[kbi];
    star[knew].se_mt = binary[kb].bse_mass[kbi];
    star[knew].se_ospin = binary[kb].bse_ospin[kbi];
    star[knew].se_B_0 = binary[kb].bse_B_0[kbi]; /* PK */
    star[knew].se_bacc = binary[kb].bse_bacc[kbi];
    star[knew].se_tacc = binary[kb].bse_tacc[kbi];
    star[knew].se_scm_B = binary[kb].bse_bcm_B[kbi];
    star[knew].se_scm_formation = binary[kb].bse_bcm_formation[kbi];
    star[knew].se_epoch = binary[kb].bse_epoch[kbi];
    star[knew].se_tphys = binary[kb].bse_tphys;
    star[knew].se_radius = binary[kb].bse_radius[kbi];
    star[knew].se_lum = binary[kb].bse_lum[kbi];
    star[knew].se_mc = binary[kb].bse_massc[kbi];
    star[knew].se_rc = binary[kb].bse_radc[kbi];
    star[knew].se_menv = binary[kb].bse_menv[kbi];
    star[knew].se_renv = binary[kb].bse_renv[kbi];
    star[knew].se_tms = binary[kb].bse_tms[kbi];
	star[knew].se_bhspin = binary[kb].bse_bhspin[kbi];
    //Sourav: toy rejuvenation- updating the rejuv variables for two cases- mass1 and mass2 
    if (kbi==0){
      star[knew].createtime = binary[kb].createtime_m1;
      star[knew].lifetime = binary[kb].lifetime_m1;
        star[knew].rad = binary[kb].rad1;
    } else {
  star[knew].createtime = binary[kb].createtime_m2;
      star[knew].lifetime = binary[kb].lifetime_m2; 
        star[knew].rad = binary[kb].rad2;
    } 
  }
}

/**
* @brief copies mass from old to new star
*
* @param oldk old star index
* @param kbi 0 or 1 for binar, -1 for non-binary
* @param knew index of new star
*/
void cp_m_to_newstar(long oldk, int kbi, long knew)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    star_m[get_global_idx(knew)] = star_m[get_global_idx(oldk)];
  } else { /* star comes from input binary */
    if (kbi == 0) {
      star_m[get_global_idx(knew)] = binary[kb].m1; //should this be multiplied by MSUN/units.mstar ?
    } else {
      star_m[get_global_idx(knew)] = binary[kb].m2;//should this be multiplied by MSUN/units.mstar ?
    }
  }
}

/**
* @brief same as cpSEvars_to_newstar, only difference in function signature
*
* @param oldk old star index
* @param kbi 0 or 1 for binary, -1 for non-binary
* @param target_star pointer to target star
*/
void cp_SEvars_to_star(long oldk, int kbi, star_t *target_star)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    target_star->rad = star[oldk].rad; //PDK addition
    target_star->se_mass = star[oldk].se_mass;
    target_star->se_k = star[oldk].se_k;
    target_star->se_mt = star[oldk].se_mt;
    target_star->se_ospin = star[oldk].se_ospin;
    target_star->se_B_0 = star[oldk].se_B_0; /* PK */
    target_star->se_bacc = star[oldk].se_bacc;
    target_star->se_tacc = star[oldk].se_tacc;
    target_star->se_scm_B = star[oldk].se_scm_B;
    target_star->se_scm_formation = star[oldk].se_scm_formation;
    target_star->se_epoch = star[oldk].se_epoch;
    target_star->se_tphys = star[oldk].se_tphys;
    target_star->se_radius = star[oldk].se_radius;
    target_star->se_lum = star[oldk].se_lum;
    target_star->se_mc = star[oldk].se_mc;
    target_star->se_rc = star[oldk].se_rc;
    target_star->se_menv = star[oldk].se_menv;
    target_star->se_renv = star[oldk].se_renv;
    target_star->se_tms = star[oldk].se_tms;
	target_star->se_bhspin = star[oldk].se_bhspin;
    target_star->Eint = star[oldk].Eint;
    //Sourav: toy rejuvenation- updating rejuv variables
    target_star->createtime = star[oldk].createtime;
    target_star->lifetime = star[oldk].lifetime;
  } else { /* star comes from input binary */
    target_star->se_mass = binary[kb].bse_mass0[kbi];
    target_star->se_k = binary[kb].bse_kw[kbi];
    target_star->se_mt = binary[kb].bse_mass[kbi];
    target_star->se_ospin = binary[kb].bse_ospin[kbi];
    target_star->se_B_0 = binary[kb].bse_B_0[kbi]; /* PK */
    target_star->se_bacc = binary[kb].bse_bacc[kbi];
    target_star->se_tacc = binary[kb].bse_tacc[kbi];
    target_star->se_scm_B = binary[kb].bse_bcm_B[kbi];
    target_star->se_scm_formation = binary[kb].bse_bcm_formation[kbi];
    target_star->se_epoch = binary[kb].bse_epoch[kbi];
    target_star->se_tphys = binary[kb].bse_tphys;
    target_star->se_radius = binary[kb].bse_radius[kbi];
    target_star->se_lum = binary[kb].bse_lum[kbi];
    target_star->se_mc = binary[kb].bse_massc[kbi];
    target_star->se_rc = binary[kb].bse_radc[kbi];
    target_star->se_menv = binary[kb].bse_menv[kbi];
    target_star->se_renv = binary[kb].bse_renv[kbi];
    target_star->se_tms = binary[kb].bse_tms[kbi];
	target_star->se_bhspin = binary[kb].bse_bhspin[kbi];
    //Sourav: toy rejuvenation- updating rejuv variables for two cases mass1 and mass2
    if (kbi==1){
        target_star->createtime = binary[kb].createtime_m1;
        target_star->lifetime = binary[kb].lifetime_m1;
        target_star->rad = binary[kb].rad1; // PDK addition
        target_star->Eint = binary[kb].Eint1;
    } else {
        target_star->createtime = binary[kb].createtime_m2;
        target_star->lifetime = binary[kb].lifetime_m2;
        target_star->rad = binary[kb].rad2; // PDK addition
        target_star->Eint = binary[kb].Eint2;
    }
  }
}

/**
* @brief same as cp_m_to_newstar, only difference in function signature
*
* @param oldk old star index
* @param kbi 0 or 1 for binary, -1 for non-binary
* @param target_star pointer to target star
*/
void cp_m_to_star(long oldk, int kbi, star_t *target_star)
{
  long kb;
  
  kb = star[oldk].binind;
  
  if (kbi == -1) { /* star comes from input single star */
    target_star->m = star_m[get_global_idx(oldk)];
  } else { /* star comes from input binary */
    if (kbi == 0) {
      target_star->m = binary[kb].m1;//should this be multiplied by MSUN/units.mstar ?
    } else {
      target_star->m = binary[kb].m2;//should this be multiplied by MSUN/units.mstar ?
    }
  }
}

/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
/**
* @brief copies stellar evoltion variables to new binary (set everything except tb)
*
* @param oldk old star index
* @param oldkbi 0 or 1 for binary, -1 for non-binary
* @param knew index of new star
* @param kbinew ?
*/
void cp_SEvars_to_newbinary(long oldk, int oldkbi, long knew, int kbinew)
{
  long kbold, kbnew;

  kbold = star[oldk].binind;
  kbnew = star[knew].binind;

  if (oldkbi == -1) { /* star comes from input single star */
    binary[kbnew].bse_mass0[kbinew] = star[oldk].se_mass;
    binary[kbnew].bse_kw[kbinew] = star[oldk].se_k;
    binary[kbnew].bse_mass[kbinew] = star[oldk].se_mt;
    binary[kbnew].bse_ospin[kbinew] = star[oldk].se_ospin;
    binary[kbnew].bse_B_0[kbinew] = star[oldk].se_B_0; /* PK */
    binary[kbnew].bse_bacc[kbinew] = star[oldk].se_bacc;
    binary[kbnew].bse_tacc[kbinew] = star[oldk].se_tacc;
    binary[kbnew].bse_bcm_B[kbinew] = star[oldk].se_scm_B;
    binary[kbnew].bse_bcm_formation[kbinew] = star[oldk].se_scm_formation;
    binary[kbnew].bse_epoch[kbinew] = star[oldk].se_epoch;
    binary[kbnew].bse_tphys = star[oldk].se_tphys; /* tphys should be the same for both input stars so this should be OK */
    binary[kbnew].bse_radius[kbinew] = star[oldk].se_radius;
    binary[kbnew].bse_lum[kbinew] = star[oldk].se_lum;
    binary[kbnew].bse_massc[kbinew] = star[oldk].se_mc;
    binary[kbnew].bse_radc[kbinew] = star[oldk].se_rc;
    binary[kbnew].bse_menv[kbinew] = star[oldk].se_menv;
    binary[kbnew].bse_renv[kbinew] = star[oldk].se_renv;
    binary[kbnew].bse_tms[kbinew] = star[oldk].se_tms;
	binary[kbnew].bse_bhspin[kbinew] = star[oldk].se_bhspin;
    //Sourav: toy rejuv- updating rejuv variables to the binary member from the single star
    if (kbinew==0){
      binary[kbnew].createtime_m1 = star[oldk].createtime;
      binary[kbnew].lifetime_m1 = star[oldk].lifetime;
      binary[kbnew].rad1 = star[oldk].rad; // PDK addition
    } else {
      binary[kbnew].createtime_m2 = star[oldk].createtime;
      binary[kbnew].lifetime_m2 = star[oldk].lifetime;
      binary[kbnew].rad2 = star[oldk].rad; // PDK addition
    }
  } else { /* star comes from input binary */
    binary[kbnew].bse_mass0[kbinew] = binary[kbold].bse_mass0[oldkbi];
    binary[kbnew].bse_kw[kbinew] = binary[kbold].bse_kw[oldkbi];
    binary[kbnew].bse_mass[kbinew] = binary[kbold].bse_mass[oldkbi];
    binary[kbnew].bse_ospin[kbinew] = binary[kbold].bse_ospin[oldkbi];
    binary[kbnew].bse_B_0[kbinew] = binary[kbold].bse_B_0[oldkbi]; /* PK */
    binary[kbnew].bse_bacc[kbinew] = binary[kbold].bse_bacc[oldkbi];
    binary[kbnew].bse_tacc[kbinew] = binary[kbold].bse_tacc[oldkbi];
    binary[kbnew].bse_bcm_B[kbinew] = binary[kbold].bse_bcm_B[oldkbi];
    binary[kbnew].bse_bcm_formation[kbinew] = binary[kbold].bse_bcm_formation[oldkbi];
    binary[kbnew].bse_epoch[kbinew] = binary[kbold].bse_epoch[oldkbi];
    binary[kbnew].bse_tphys = binary[kbold].bse_tphys;
    binary[kbnew].bse_radius[kbinew] = binary[kbold].bse_radius[oldkbi];
    binary[kbnew].bse_lum[kbinew] = binary[kbold].bse_lum[oldkbi];
    binary[kbnew].bse_massc[kbinew] = binary[kbold].bse_massc[oldkbi];
    binary[kbnew].bse_radc[kbinew] = binary[kbold].bse_radc[oldkbi];
    binary[kbnew].bse_menv[kbinew] = binary[kbold].bse_menv[oldkbi];
    binary[kbnew].bse_renv[kbinew] = binary[kbold].bse_renv[oldkbi];
    binary[kbnew].bse_tms[kbinew] = binary[kbold].bse_tms[oldkbi];
    binary[kbnew].bse_bhspin[kbinew] = binary[kbold].bse_bhspin[oldkbi];
    //Sourav: toy rejuv- updating rejuv variables to binary members from binary members
    //There can be four cases. m1,2(old)->m1,2(new) 
    if(kbinew==0){
      if (oldkbi==0){
        binary[kbnew].createtime_m1 = binary[kbold].createtime_m1;
        binary[kbnew].lifetime_m1 = binary[kbold].lifetime_m1;
        binary[kbnew].rad1 = binary[kbold].rad1; // PDK addition
      } else {
        binary[kbnew].createtime_m1 = binary[kbold].createtime_m2;
        binary[kbnew].lifetime_m1 = binary[kbold].lifetime_m2;
        binary[kbnew].rad1 = binary[kbold].rad2; // PDK addition
      }
    } else {
      if (oldkbi==0){
        binary[kbnew].createtime_m2 = binary[kbold].createtime_m1;
        binary[kbnew].lifetime_m2 = binary[kbold].lifetime_m1;
        binary[kbnew].rad2 = binary[kbold].rad1; // PDK addition
      } else {
        binary[kbnew].createtime_m2 = binary[kbold].createtime_m2;
        binary[kbnew].lifetime_m2 = binary[kbold].lifetime_m2;
        binary[kbnew].rad2 = binary[kbold].rad2; // PDK addition
      }
    }
  }
}


/* olsk=old star index; kbi=0,1 for binary, -1 for non-binary; knew=index of new star */
/* set everything except tb */
/**
* @brief ?
*
* @param instar ?
* @param binindex ?
* @param bid ?
*/
void cp_starSEvars_to_binmember(star_t instar, long binindex, int bid)
{
  binary[binindex].bse_mass0[bid] = instar.se_mass;
  binary[binindex].bse_kw[bid] = instar.se_k;
  binary[binindex].bse_mass[bid] = instar.se_mt;
  binary[binindex].bse_ospin[bid] = instar.se_ospin;
  binary[binindex].bse_B_0[bid] = instar.se_B_0; /* PK */
  binary[binindex].bse_bacc[bid] = instar.se_bacc;
  binary[binindex].bse_tacc[bid] = instar.se_tacc;
  binary[binindex].bse_bcm_B[bid] = instar.se_scm_B;
  binary[binindex].bse_bcm_formation[bid] = instar.se_scm_formation;
  binary[binindex].bse_epoch[bid] = instar.se_epoch;
  binary[binindex].bse_tphys = instar.se_tphys; /* tphys should be the same for both input stars so this should be OK */
  binary[binindex].bse_radius[bid] = instar.se_radius;
  binary[binindex].bse_lum[bid] = instar.se_lum;
  binary[binindex].bse_massc[bid] = instar.se_mc;
  binary[binindex].bse_radc[bid] = instar.se_rc;
  binary[binindex].bse_menv[bid] = instar.se_menv;
  binary[binindex].bse_renv[bid] = instar.se_renv;
  binary[binindex].bse_tms[bid] = instar.se_tms;
  binary[binindex].bse_bhspin[bid] = instar.se_bhspin;
  //Sourav: toy rejuv- updating rejuv variables from a single to a binary member
  if (bid==0){
    binary[binindex].createtime_m1 = instar.createtime;
    binary[binindex].lifetime_m1 = instar.lifetime;
    binary[binindex].rad1 = instar.rad; // PDK addition
    binary[binindex].Eint1 = instar.Eint;
  } else {
    binary[binindex].createtime_m2 = instar.createtime;
    binary[binindex].lifetime_m2 = instar.lifetime;
    binary[binindex].rad2 = instar.rad; // PDK addition
    binary[binindex].Eint2 = instar.Eint;
  }
}

/**
* @brief ?
*
* @param instar ?
* @param binindex ?
* @param bid ?
*/
void cp_starmass_to_binmember(star_t instar, long binindex, int bid)
{
  if (bid == 0) {
    binary[binindex].m1 = instar.m;
  } else {
    binary[binindex].m2 = instar.m;
 }
}

int rhs_peters(double t, const double y[], double f[], void *params){
	/* dadt and dedt from Peters 1964; note that this is entirely in code units! */
	double *mG3c5 = (double *) params;
	f[0] = -((12.8) * (*mG3c5) / (y[0]*y[0]*y[0]*pow(1-y[1]*y[1],3.5)))
			*(1+3.041666666666667*y[1]*y[1] + 0.3854166666666667*y[1]*y[1]*y[1]*y[1]);
	f[1] = -((20.2666666666667) * y[1] * (*mG3c5) / (y[0]*y[0]*y[0]*y[0]*pow(1-y[1]*y[1],2.5)))
			*(1+0.3980263157894737*y[1]*y[1]);

	return GSL_SUCCESS;
}

void integrate_a_e_peters_eqn(long kb){
	double t = 0.;
	double h = 1.e-9;
	double eps_abs = 1.e-10;
	double eps_rel = 1.e-10;
	int collision = 0, status;

	/* First, allocate the GSL integrator and all it's accoutrements*/
	const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc(type_ptr, 2);
	gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
	gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc(2);
	gsl_odeiv2_system my_system;

	/* Then define the system; note that this is in code units*/
	double m1 = binary[kb].m1;
	double m2 = binary[kb].m2;
	double clight = 2.9979e10 / (units.l/units.t); 
	// main prefactor for Peters equations, in code units
	double mG3c5 = m1*m2*(m1+m2)*madhoc*madhoc*madhoc / (clight*clight*clight*clight*clight);
	double y[2] = {binary[kb].a, binary[kb].e};
	
	my_system.function = rhs_peters;
	my_system.dimension = 2;
	my_system.params = &mG3c5;


	/* Finally, integrate for this timestep; remember that Dt is in relaxation times,
     * NOT N-body times */
    double t_final = Dt * ((double) clus.N_STAR)/ log(GAMMA * ((double) clus.N_STAR));
	while (t < t_final){
		status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr,
								&my_system, &t, t_final, &h, y); 

		//fprintf(stderr,"lol = %g %g %g %g %g \n",y[0]*units.l/AU,y[1],m1*units.mstar/MSUN,m2*units.mstar/MSUN,t);
		/* Check for collisions at periapse */
		if(y[0]*(1-y[1]) < BH_RADIUS_MULTIPLYER*(binary[kb].rad1 + binary[kb].rad2)){
			collision = 1;
			break;
		}

		/* Make sure the integrator does what it says it does...*/
		if(status != GSL_SUCCESS){
		    eprintf("GSL Integrator for the Peters equations failed; yell at Carl\n");
		    exit_cleanly(-1, __FUNCTION__);
			break;
		}
	}

	/* If we have a collision, then set the BSE mass to zero, and
	 * handle_bse_outcome will treat it correctly
	 *
	 * Otherwise, update eccentricity and semi-major axis accordingly */
	if(collision){
		binary[kb].bse_mass[0] += binary[kb].bse_mass[1];
		binary[kb].bse_mass[1] = 0.;
		binary[kb].bse_tb = 0.;
	} else {
		binary[kb].a = y[0];
		binary[kb].e = y[1];
		binary[kb].bse_tb = sqrt(cub(binary[kb].a * units.l / AU)/(binary[kb].bse_mass[0]+binary[kb].bse_mass[1]))*365.25;
		/* This is used to set a again in handle_bse_outcome; easier to just set
		 * it here as well*/
	}

	/* all done; free up the gsl_odeiv stuff */
	gsl_odeiv2_evolve_free (evolve_ptr);
	gsl_odeiv2_control_free (control_ptr);
	gsl_odeiv2_step_free (step_ptr);

	return;
}

