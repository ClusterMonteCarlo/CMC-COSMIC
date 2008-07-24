struct perturbation {
  double E[2], J[2];
  double dE[2], dJ[2];
  double vr[2], vt[2];
  double v[4], vp[4];
  long index[2];
};

struct encounter {
  long k, kp;
  double r, rp;
  double v[4], vp[4];
  double w[4], W;
  double rcm, vcm[4];
  double Y;
};

struct relaxation_params {
  double Trel12;
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

