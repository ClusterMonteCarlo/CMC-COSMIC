/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief
* These routines are for calculating the orbital perturbations of a pair of
* stars due to relaxation. They still depend on the global 'star' array and
* are thus not fully modular.
*
* @param index index of star
*
* @return ?
*/
struct star_coords
get_star_coords_from_star(long index) {
  struct star_coords pos;

  pos.index= index;
  pos.field_index= index;
  pos.E = star[index].E;
  pos.J = star[index].J;
  pos.r = star[index].r;
  pos.vr= star[index].vr;
  pos.vt= star[index].vt;
  pos.pot= star[index].phi;

  return(pos);
}

/**
* @brief ?
*
* @param pos ?
*/
void set_star_coords_for_star(struct star_coords pos) {
  long index;

  index= pos.field_index;

  star[index].E= pos.E;
  star[index].J= pos.J;
  star[index].r= pos.r;
  star[index].vr= pos.vr;
  star[index].vt= pos.vt;
  star[index].phi= pos.pot;
}

/**
* @brief ?
*
* @param star1 ?
* @param star2 ?
* @param rng ?
*
* @return ?
*/
struct encounter
get_encounter_dyns(struct star_coords star1, struct star_coords star2, gsl_rng *rng) {
  struct encounter enc;
  int j;
  double phi;
  long k, kp;

  /* set random angle between vt's */
  /* first store random variable */
  enc.Y = rng_t113_dbl();
  phi = enc.Y * 2.0 * PI;
  enc.v[1] = star1.vt;
  enc.v[2] = 0.0;
  enc.v[3] = star1.vr;
  enc.vp[1] = star2.vt * cos(phi);
  enc.vp[2] = star2.vt * sin(phi);
  enc.vp[3] = star2.vr;

  for (j=1; j<=3; j++) {
    enc.w[j] = enc.vp[j] - enc.v[j];
  }

  enc.W = sqrt(sqr(enc.w[1]) + sqr(enc.w[2]) + sqr(enc.w[3]));
  if (enc.W == 0.0) {
    eprintf("W = 0!\n");
    exit_cleanly(1, __FUNCTION__);
  }

  /* compute CM quantities */
  k= star1.index;
  kp= star2.index;
  enc.rcm = (star[k].m * star1.r + star[kp].m * star2.r) / (star[k].m + star[kp].m);
  for (j=1; j<=3; j++) {
    enc.vcm[j] = (star[k].m * enc.v[j] + star[kp].m * enc.vp[j]) / (star[k].m + star[kp].m);
  };

  enc.k= k;
  enc.kp= kp;
  enc.r= star1.r;
  enc.rp= star2.r;

  return(enc);
};

/**
* @brief ?
*
* @param dt ?
* @param enc ?
*
* @return ?
*/
struct relaxation_params
get_relaxation_params(double dt, struct encounter enc) {
  struct relaxation_params params;
  double n_local;
  long k, kp;

  k= enc.k;
  kp= enc.kp;
  n_local = calc_n_local(kp, AVEKERNEL, clus.N_MAX);
  params.Trel12 = (PI/32.0) * cub(enc.W) / ( ((double) clus.N_STAR) * n_local * 
      sqr((star[k].m+star[kp].m)*madhoc) ) ;
  params.beta = (PI/2.0) * sqrt(dt/params.Trel12);
  dprintf("kp=%li, beta=%g, fdt=%g\n", kp, params.beta, sqr(params.beta*2.0/PI));

  return(params);
};

/**
* @brief ?
*
* @param enc ?
*
* @return ?
*/
double get_relaxation_time(struct encounter enc) {
  double n_local, Trel12;
  long k, kp;

  k= enc.k;
  kp= enc.kp;
  n_local = calc_n_local(kp, AVEKERNEL, clus.N_MAX);
  Trel12 = (PI/32.0) * cub(enc.W) / ( ((double) clus.N_STAR) * n_local * 
      sqr((star[k].m+star[kp].m)*madhoc) ) ;
  return(Trel12);
};

/**
* @brief ?
*
* @param dt ?
* @param Trel ?
*
* @return ?
*/
double scattering_angle(double dt, double Trel) {
  double beta;

  beta = (PI/2.0) * sqrt(dt/Trel);
  dprintf("beta=%g, fdt=%g\n", beta, sqr(beta*2.0/PI));

  return(beta);
};

/**
* @brief ?
*
* @param enc ?
* @param rparams ?
*
* @return ?
*/
struct perturbation
scatter_relax(struct encounter enc, struct relaxation_params rparams) {
  struct perturbation pert;
  long k, kp;
  double wp, w1[4], w2[4];
  double psi, w_new[4], v_new[4], vp_new[4];
  double beta, Trel12;
  int j;

  k= enc.k;
  kp= enc.kp;

  Trel12= rparams.Trel12;
  beta= rparams.beta;

  /* clamp beta at max value */
  if (beta > PI/2.0) {
    beta = PI/2.0;
  }

  /* set up coordinate system */
  wp = sqrt(sqr(enc.w[1]) + sqr(enc.w[2]));
  if (wp == 0.0) {
    eprintf("wp=0 \n");
    exit_cleanly(1, __FUNCTION__);
  }

  /* You'll notice here that the sign on w1 is opposite that of what's shown in Kris Joshi's
     paper.  The sign has now been fixed so that (w1, w2, w) define a right-handed coordinate
     system, as such: \^w1 x \^w2 = \^w */
  w1[1] = -enc.w[2] * enc.W / wp;
  w1[2] = enc.w[1] * enc.W / wp;
  w1[3] = 0.0;
  w2[1] = -enc.w[1] * enc.w[3] / wp;
  w2[2] = -enc.w[2] * enc.w[3] / wp;
  w2[3] = wp;

  psi = rng_t113_dbl() * 2 * PI;
  for (j = 1; j <= 3; j++) {
    w_new[j] = enc.w[j] * cos(beta) + w1[j] * sin(beta) * cos(psi) + w2[j] * sin(beta) * sin(psi);
  }

  for (j = 1; j <= 3; j++) {
    v_new[j] = enc.v[j] - star[kp].m / (star[k].m + star[kp].m) * (w_new[j] - enc.w[j]);
    vp_new[j] = enc.vp[j] + star[k].m / (star[k].m + star[kp].m) * (w_new[j] - enc.w[j]);
  }

  /* set new velocities for both stars */
  pert.vr[0] = v_new[3];
  pert.vt[0] = sqrt(sqr(v_new[1]) + sqr(v_new[2]));
  pert.vr[1] = vp_new[3];
  pert.vt[1] = sqrt(sqr(vp_new[1]) + sqr(vp_new[2]));

  for (j=1; j<=3; j++) {
    pert.v[j]= v_new[j];
    pert.vp[j]= vp_new[j];
  };

  pert.dE[0]= 0.5*(sqr(v_new[1])+sqr(v_new[2])+sqr(v_new[3]))- 
    0.5*(sqr(enc.v[1])+sqr(enc.v[2])+sqr(enc.v[3]));
  pert.dJ[0]= enc.r*(pert.vt[0]- sqrt(sqr(enc.v[1])+ sqr(enc.v[2])));
  pert.dE[1]= 0.5*(sqr(vp_new[1])+sqr(vp_new[2])+sqr(vp_new[3]))- 
    0.5*(sqr(enc.vp[1])+sqr(enc.vp[2])+sqr(enc.vp[3]));
  pert.dJ[1]= enc.rp*(pert.vt[1]- sqrt(sqr(enc.vp[1])+ sqr(enc.vp[2])) );

  return(pert);
}

/**
* @brief ?
*
* @param pos[2] ?
* @param dt ?
* @param rng ?
*
* @return ?
*/
struct perturbation
scatter_relax_old(struct star_coords pos[2], double dt, gsl_rng *rng) {
  struct perturbation pert;
  struct encounter enc;
  struct relaxation_params rparams;
  long k, kp;
  double wp, w1[4], w2[4];
  double psi, w_new[4], v_new[4], vp_new[4];
  double beta, Trel12;
  int j;

  pert.index[0]= pos[0].index;
  pert.index[1]= pos[1].index;

  k= pos[0].index;
  kp= pos[1].index;

  /* set dynamical params for this pair */
  enc= get_encounter_dyns(pos[0], pos[1], rng);
  rparams= get_relaxation_params(dt, enc);
  Trel12= rparams.Trel12;
  beta= rparams.beta;

  /* clamp beta at max value */
  if (beta > PI/2.0) {
    beta = PI/2.0;
  }

  /* set up coordinate system */
  wp = sqrt(sqr(enc.w[1]) + sqr(enc.w[2]));
  if (wp == 0.0) {
    eprintf("wp=0 \n");
    exit_cleanly(1, __FUNCTION__);
  }

  /* You'll notice here that the sign on w1 is opposite that of what's shown in Kris Joshi's
     paper.  The sign has now been fixed so that (w1, w2, w) define a right-handed coordinate
     system, as such: \^w1 x \^w2 = \^w */
  w1[1] = -enc.w[2] * enc.W / wp;
  w1[2] = enc.w[1] * enc.W / wp;
  w1[3] = 0.0;
  w2[1] = -enc.w[1] * enc.w[3] / wp;
  w2[2] = -enc.w[2] * enc.w[3] / wp;
  w2[3] = wp;

  psi = rng_t113_dbl() * 2 * PI;
  for (j = 1; j <= 3; j++) {
    w_new[j] = enc.w[j] * cos(beta) + w1[j] * sin(beta) * cos(psi) + w2[j] * sin(beta) * sin(psi);
  }

  for (j = 1; j <= 3; j++) {
    v_new[j] = enc.v[j] - star[kp].m / (star[k].m + star[kp].m) * (w_new[j] - enc.w[j]);
    vp_new[j] = enc.vp[j] + star[k].m / (star[k].m + star[kp].m) * (w_new[j] - enc.w[j]);
  }

  /* set new velocities for both stars */
  pert.vr[0] = v_new[3];
  pert.vt[0] = sqrt(sqr(v_new[1]) + sqr(v_new[2]));
  pert.vr[1] = vp_new[3];
  pert.vt[1] = sqrt(sqr(vp_new[1]) + sqr(vp_new[2]));

  for (j=1; j<=3; j++) {
    pert.v[j]= v_new[j];
    pert.vp[j]= vp_new[j];
  };

  /* Calculate new energies by recomputing E = PE + KE using new velocity*/ 
  pert.E[0]= pos[0].pot + 0.5 * (sqr(pert.vr[0]) + sqr(pert.vt[0]));
  pert.J[0]= pos[0].r * pert.vt[0];
  pert.E[1]= star[kp].phi + 0.5 * (sqr(pert.vr[1]) + sqr(pert.vt[1]));
  pert.J[1]= pos[1].r * pert.vt[1];

  pert.dE[0]= pert.E[0]- pos[0].E;
  pert.dJ[0]= pert.J[0]- pos[0].J;
  pert.dE[1]= pert.E[1]- pos[1].E;
  pert.dJ[1]= pert.J[1]- pos[1].J;

  return(pert);
}

