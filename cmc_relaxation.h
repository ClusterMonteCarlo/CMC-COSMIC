/**
* @brief ?
*/
struct perturbation {
/**
* @brief ?
*/
  double E[2], J[2];
/**
* @brief ?
*/
  double dE[2], dJ[2];
/**
* @brief ?
*/
  double vr[2], vt[2];
/**
* @brief ?
*/
  double v[4], vp[4];
/**
* @brief ?
*/
  long index[2];
};

/**
* @brief ?
*/
struct encounter {
/**
* @brief ?
*/
  long k, kp;
/**
* @brief ?
*/
  double r, rp;
/**
* @brief ?
*/
  double v[4], vp[4];
/**
* @brief ?
*/
  double w[4], W;
/**
* @brief ?
*/
  double rcm, vcm[4];
/**
* @brief ?
*/
  double Y;
};

/**
* @brief ?
*/
struct relaxation_params {
/**
* @brief ?
*/
  double Trel12;
/**
* @brief ?
*/
  double beta;
};

struct perturbation
scatter_relax_old(struct star_coords pos[2], double dt, gsl_rng *rng);

struct perturbation
scatter_relax(struct encounter enc, struct relaxation_params rparams);

struct encounter
get_encounter_dyns(struct star_coords star1, struct star_coords star2, gsl_rng *rng);

struct relaxation_params
get_relaxation_params(double dt, struct encounter enc);

double get_relaxation_time(struct encounter enc);

double scattering_angle(double dt, double Trel);

struct star_coords
get_star_coords_from_star(long index);

void set_star_coords_for_star(struct star_coords pos);

