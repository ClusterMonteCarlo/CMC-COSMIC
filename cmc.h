/* -*- linux-c -*- */

#include <sys/times.h>
#include <zlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <fitsio.h>
#include "common/fitslib.h"
#include "common/taus113-v2.h"
#include "fewbody-0.24/fewbody.h"

// This automatic SVN versioning only updates the version of this file.
//#define CMCVERSION "$Revision$"
//#define CMCDATE "$Date$"
#define CMCNICK "SVN"
#define CMCPRETTYNAME "ClusterMonteCarlo"

/*************************** Parameters ******************************/
/* Large number, but still SF_INFINITY - 1 <> SF_INFINITY */
#define SF_INFINITY 1.0e10
/* Radius of the zeroth star so that 1/star[0].r is finite */
#define ZERO 1.0e-20

#define PI 3.14159265358979323

#define N_TRY 50000

#define AVEKERNEL 20

#define MSUN 1.989e+33
#define SOLAR_MASS MSUN
#define RSUN 6.9599e10
#define G 6.67259e-8
#define KPC 3.0857e+21
#define YEAR 3.155693e+7
#define AU 1.496e+13
#define PARSEC 3.0857e+18

/* defining maximum pericenters for "strong" interactions */
/* rperi = XCOLL * (R_1 + R_2) */
#define XCOLL 1.0
/* rperi = XBS * a */
#define XBS 2.0
/* rperi = XBB * (a_1 + a_2) */
#define XBB 2.0
/* hard soft boundary XHS=vorb/W */
#define XHS 0.7
/* safety factor used in the black hole accretion routine */
/* by Freitag & Benz (2002) */
#define CSAFE 0.2

/* star types, from StarTrack */
#define S_MASS_MS_STAR                    (0)
#define L_MASS_MS_STAR                    (1)
#define HERTZ_GAP_STAR                    (2)
#define FIRST_GIANT_BRANCH_STAR           (3)
#define CORE_He_BURNING_STAR              (4)
#define EARLY_ASYM_GIANT_BRANCH_STAR      (5)
#define THERM_PULS_ASYM_GIANT_BRANCH_STAR (6)
#define MS_NAKED_He_STAR                  (7)
#define HERTZ_GAP_NAKED_He_STAR           (8)
#define GIANT_BRANCH_NAKED_He_STAR        (9)
#define He_WD                             (10)
#define CO_WD                             (11)
#define ONe_WD                            (12)
#define NEUTRON_STAR                      (13)
#define BLACK_HOLE                        (14)
#define MASSLESS_REMNANT                  (15)
#define NAKED_CO_CORE                     (-1)
#define NOT_A_STAR                        (-100)

/* binaries */
typedef struct{
	long id1; /* unique id of star 1 */
	long id2; /* unique id of star 2 */
	double rad1; /* radius of star 1 */
	double rad2; /* radius of star 2 */
	double m1; /* mass of star 1 */
	double m2; /* mass of star 2 */
	double Eint1; /* internal energy of star 1 */
	double Eint2; /* internal energy of star 2 */
	double a; /* semimajor axis */
	double e; /* eccentricity */
	int inuse; /* whether or not binary exists */
} binary_t;

/* single stars/objects */
typedef struct{
	/* dynamical evolution variables */
	double r;  /* radial coordinate */
	double vr; /* radial velocity */
	double vt; /* tangential velocity */
	double m; /* mass */
	double E; /* kinetic plus potential energy per unit mass */
	double J; /* angular momentum per unit mass */
	double EI; /* "intermediate" energy per unit mass */
	double Eint; /* internal energy (due to collisions, e.g.) */
	double rnew; /* new radial coordinate */
	double vrnew; /* new radial velocity */
	double vtnew; /* new tangential velocity */
	double rOld; /* ? */
	double X; /* ? */
	double Y; /* random variable that must be stored */
	double r_peri; /* pericenter distance */
	double r_apo; /* apocenter distance */
	double phi; /* value of potential at position of star (only updated at end of timestep) */
	long   interacted; /* whether or not the star has undergone a strong interaction (i.e., not relaxation) */
	long   binind; /* index to the binary */
	long   id; 	/* the star's unique identifier */
	double rad; /* radius */
	double Uoldrold, Uoldrnew; /* variables for Stodolkiewicz */
	double vtold, vrold;       /* energy conservation scheme  */
	/* stellar evolution variables */
	double mzams;  /* ZAMS mass                                       */
	double m0;     /* initial mass                                    */
	double mass;   /* real mass, this is the one relevant to dynamics */
	double tbeg;   /* time where evolution step begins                */
	double tvir;   /* virtual/internal time                           */
	double tend;   /* time where evolution step is supposed to end    */
	double lum;    /* luminosity                                      */
	double mc;     /* core mass                                       */
	double mcHe;   /* He core mass                                    */
	double mcCO;   /* C/O core mass                                   */
	double dt;     /* ???                                             */
	double mpre;   /* pre SN mass                                     */
	double tstart; /* ???                                             */
	long init_no;  /* for bookkeeping                                 */
	int k;	       /* star type                                       */ 
	int flag;      /* ???                                             */
	int kpre;      /* pre SN type                                     */
} star_t;

struct CenMa{
	double m;
	double E;
};

struct get_pos_str {
	double max_rad;
	double phi_rtidal;
	double phi_zero;
	long N_LIMIT;
	int taskid;
	struct CenMa CMincr;
	gsl_rng *thr_rng;
};

typedef struct{
	int BINBIN;
	int BINSINGLE;
	int CENTRAL_MASS;
	int SNAPSHOTTING;
	int SNAPSHOT_DELTAT;
	int SNAPSHOT_DELTACOUNT;
        int SNAPSHOT_CORE_COLLAPSE;
        int SNAPSHOT_CORE_BOUNCE;
	int IDUM;
	int INPUT_FILE;
	int MASS_PC;
	int MASS_BINS;
	int MINIMUM_R;
	int STOPATCORECOLLAPSE;
	int NUM_CENTRAL_STARS;
	int PERTURB;
	int RELAXATION;
	int THETASEMAX;
	int STELLAR_EVOLUTION;
	int SS_COLLISION;
	int TERMINAL_ENERGY_DISPLACEMENT;
	int T_MAX;
	int T_MAX_COUNT;
	int MAX_WCLOCK_TIME;
	int WIND_FACTOR;
	int GAMMA;
        int SEARCH_GRID;
        int SG_STARSPERBIN;
        int SG_MAXLENGTH;
        int SG_MINLENGTH;
        int SG_POWER_LAW_EXPONENT;
        int SG_MATCH_AT_FRACTION;
        int SG_PARTICLE_FRACTION;
        int BH_LOSS_CONE;
        int BH_R_DISRUPT_NB;
        int FORCE_RLX_STEP; 
#ifdef DEBUGGING
        int BH_LC_FDT;
#endif
        int APSIDES_PRECISION;
        int APSIDES_MAX_ITER;
        int APSIDES_CONVERGENCE;
} parsed_t;

/* a struct containing the units used */
typedef struct{
	double t; /* time */
	double m; /* mass */
	double l; /* length */
	double E; /* energy */
	double mstar; /* stars' masses are kept in different units */
} units_t;

/* a struct for the potential, must be malloc'ed */
typedef struct{
	long n; /* number of elements */
	double *r; /* radius */
	double *phi; /* potential */
} potential_t;

/* a struct for the force, must be malloc'ed */
typedef struct{
	long n; /* number of elements */
	double *r; /* radius */
	double *force; /* force */
} force_t;

/* useful structure for central quantities */
typedef struct{
	double rho; /* central mass density */
	double v_rms; /* rms object velocity */
	double rc; /* core radius */
	double m_ave; /* average object mass */
	double n; /* number density of objects */
	double rc_spitzer; /* Spitzer definition of core radius */
	long N_sin; /* number of objects that are single */
	long N_bin; /* number of objects that are binary */
	double n_sin; /* single star number density */
	double n_bin; /* binary star number density */
	double rho_sin; /* central single star mass density */
	double rho_bin; /* central binary star mass density */
	double m_sin_ave; /*average single star mass */
	double m_bin_ave; /* average binary star mass */
	double v_sin_rms; /* rms single star velocity */
	double v_bin_rms; /* rms binary star velocity */
	double w2_ave; /* average of 2*m*v^2 per average mass for all objects */
	double R2_ave; /* average of R^2 for single stars */
	double mR_ave; /* average of m*R for single stars */
	double a_ave; /* average of a for binaries */
	double a2_ave; /* average of a^2 for binaries */
	double ma_ave; /* average of m*a for binaries */
} central_t;

/* useful structure for core quantities */
typedef struct{
	double N;
	double kT;
} core_t;

/* to store the velocity dispersion profile */
typedef struct{
	long n;
	double *r;
	double *sigma;
} sigma_t;

/* parameters for orbit */
typedef struct{
	double rp;
	double ra;
	double dQdrp;
	double dQdra;
	long kmax;
	long kmin;
	int circular_flag;
} orbit_rs_t;

/* parameters for calc_p_orb function */
typedef struct{
	double E;
	double J;
	long index;
	long kmin;
	long kmax;
	double rp;
	double ra;
} calc_p_orb_params_t;

/* other useful structs */
typedef struct{
	long N_MAX;
	long N_MAX_NEW;
	long N_STAR;
	long N_STAR_NEW;
	long N_BINARY;
} clus_struct_t;

typedef struct{
	double tot; /* total = kinetic + potential + internal + binary binding + central mass energy */
	double New; /* ??? */
	double ini; /* initial total energy */
	double K; /* total kinetic */
	double P; /* total potential */
	double Eint; /* total internal */
	double Eb; /* total binary binding energy */
} Etotal_struct_t;

typedef struct{
	double rh; /* half-mass radius */
} clusdyn_struct_t;

/********************** Function Declarations ************************/
double sqr(double x);
double cub(double x);
void tidally_strip_stars(void);
void remove_star_center(long j);
void print_results(void);
double potential(double r);	       /* get potential using star.phi */
double fastpotential(double r, long kmin, long kmax);
long potential_calculate(void);	/* calculate potential at star locations in star.phi */
void comp_mass_percent(void);
void comp_multi_mass_percent(void);
orbit_rs_t calc_orbit_rs(long si, double E, double J);
double get_positions(void);	/* get positions and velocities */
void perturb_stars(double Dt);	/* take a time step (perturb E,J) */
long FindZero_r(long x1, long x2, double r);
long FindZero_Q(long j, long x1, long x2, double E, double J);
double potentialDifference(int particleIndex);
void ComputeEnergy(void);
void PrintLogOutput(void);
double GetTimeStep(gsl_rng *rng);
long CheckStop(struct tms tmsbuf);
void ComputeIntermediateEnergy(void);
void set_velocities3(void);
double rtbis(double (*func) (double), double x1, double x2, double xacc);
double trapzd(double (*func) (double), double a, double b, long n);
double qsimp(double (*func) (double), double a, double b);

void splint(double xa[], double ya[], double y2a[], long n, double x, double *y);
void spline(double x[], double y[], long n, double yp1, double ypn, double y2[]);
void print_2Dsnapshot(void);
void get_physical_units(void);
void update_vars(void);

void print_version(FILE *stream);
void print_usage(FILE *stream, char *argv[]);
void parser(int argc, char *argv[], gsl_rng *r);
void PrintFileOutput(void);

/* Unmodified Numerical Recipes Routines */
void nrerror(char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);

/* fits stuff */
void load_fits_file_data(void);

/* stellar evolution stuff */
void stellar_evolution_init(void);
void do_stellar_evolution(void);
void write_stellar_data(void);
double r_of_m(double M);

/* Fewbody stuff */
void destroy_obj(long i);
void destroy_binary(long i);
long create_star(void);
long create_binary(void);
void dynamics_apply(double dt, gsl_rng *rng);
void perturb_stars_fewbody(double dt, gsl_rng *rng);
void qsorts(star_t *s, long n);
void units_set(void);
void central_calculate(void);
void clusdyn_calculate(void);
void print_interaction_status(char status_text[]);
void print_interaction_error(void);
long star_get_id_new(void);
double calc_n_local(long k, long p, long N_LIMIT);
double calc_Ai_local(long k, long kp, long p, double W, long N_LIMIT);
void calc_encounter_dyns(long k, long kp, double v[4], double vp[4], double w[4], double *W, double *rcm, double vcm[4], gsl_rng *rng, int setY);
void set_star_EJ(long k);
void set_star_news(long k);
void set_star_olds(long k);
void zero_star(long j);
void zero_binary(long j);

void sscollision_do(long k, long kp, double rcm, double vcm[4]);

void print_initial_binaries(void);

void bs_calcunits(fb_obj_t *obj[2], fb_units_t *bs_units);
fb_ret_t binsingle(double *t, long ksin, long kbin, double W, double bmax, fb_hier_t *hier, gsl_rng *rng);

void bb_calcunits(fb_obj_t *obj[2], fb_units_t *bb_units);
fb_ret_t binbin(double *t, long k, long kp, double W, double bmax, fb_hier_t *hier, gsl_rng *rng);

double binint_get_mass(long k, long kp, long id);
void binint_log_obj(fb_obj_t *obj, fb_units_t units);
void binint_log_status(fb_ret_t retval);
void binint_log_collision(const char interaction_type[], long id, double mass, double r, fb_obj_t obj, long k, long kp);
void binint_do(long k, long kp, double rperi, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng);

double simul_relax(gsl_rng *rng);
void break_wide_binaries(void);
void calc_sigma_r(void);
double sigma_r(double r);

/* signal/GSL error handling stuff */
void toggle_debugging(int signal);
void exit_cleanly(int signal);
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno);
void close_buffers(void);
void trap_sigs(void);
void free_arrays(void);

/* if we are using gcc, we can remove the annoying warning */
#ifdef __GNUC__
static double SQR(double a) __attribute__ ((always_inline,pure));
static inline double SQR(double a){return a*a;}
#else
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#endif

/* black hole accretion stuff (loss cone) */
void get_3d_velocities(double *w, double vr, double vt);
void do_random_step(double *w, double beta, double delta);
double check_angle_w_w_new(double *w, double *w_new, double delta);
double calc_P_orb(long index);
double calc_p_orb_f(double x, void *params);
double calc_p_orb_f2(double x, void *params);
double calc_p_orb_gc(double x, void *params);
void bh_rand_walk(long index, double v[4], double vcm[4], double beta, double dt);

/* potential calculation speed-up*/
long check_if_r_around_last_index(long last_index, double r);

struct Search_Grid {
   /* user modifiable parameters */
   /* starsPerBin is the estimated (anticipated) number of stars per search
    * bin. It is only used to calculate interpol_coeff and grid length.
    */
   long starsPerBin;
   long max_length;
   long min_length;
   /* fraction determines the interpol_coeff such that 
    * grid.radius[fraction*grid->length] is approx. 
    * star[clus.N_MAX*fraction].r. Depending how well the power_law_exponent 
    * fits the distribution of star[i::starsPerBin].r the better this is 
    * fulfilled*/
   double fraction;
   /* Only particle_fraction * clus.N_MAX stars should be used to construct the
    * grid (i.e. we ignore a possible outer halo).*/
   double particle_fraction;
   double power_law_exponent;
   /* You should not change any of the following variables yourself 
    * unless, of course, you know what you are doing.
    */
   long *radius;
   long length;
   double interpol_coeff;
};

struct Interval {
  long min;
  long max;
};

struct Search_Grid *
search_grid_initialize(double power_law_exponent, double fraction, long starsPerBin, double part_frac);

double search_grid_estimate_prop_const(struct Search_Grid *grid);
/* This function does not belong to the public API! */
void search_grid_allocate(struct Search_Grid *grid);
void search_grid_update(struct Search_Grid *grid);

struct Interval
search_grid_get_interval(struct Search_Grid *grid, double r);

long search_grid_get_grid_index(struct Search_Grid *grid, double r);
double search_grid_get_r(struct Search_Grid *grid, long index);
double search_grid_get_grid_indexf(struct Search_Grid *grid, double r);
void search_grid_free(struct Search_Grid *grid);

#ifdef DEBUGGING
#include <glib.h>
void load_id_table(GHashTable* ids, char *filename);
double calc_average_mass_sqr(long index, long N_LIMIT);
#endif
struct Interval get_r_interval(double r);
double calc_vr_within_interval(double r, void *p);
double calc_vr(double r, long index, double E, double J);
double find_root_vr(long index, long k, double E, double J);

/* macros */ 
/* correction to potential due to subtracting star's contribution, and adding self-gravity */
/* #define PHI_S(rad, j) ( ((rad)>=star[(j)].r ? star[(j)].m/(rad) : star[(j)].m/star[(j)].r) * (1.0/clus.N_STAR) - 0.5*star[(j)].m/clus.N_STAR/(rad) ) */
/* correction to potential due to subtracting star's contribution (this is what Marc uses) */
#define PHI_S(rad, j) ( ((rad)>=star[(j)].r ? star[(j)].m/(rad) : star[(j)].m/star[(j)].r) * (1.0/clus.N_STAR) )
#define function_Q(j, k, E, J) (2.0 * ((E) - (star[(k)].phi + PHI_S(star[(k)].r, j))) - SQR((J) / star[(k)].r))
#define gprintf(args...) if (!quiet) {fprintf(stdout, args);}
#define diaprintf(args...) if (!quiet) {fprintf(stdout, "DIAGNOSTIC: %s(): ", __FUNCTION__); fprintf(stdout, args);}
#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#define wprintf(args...) {fprintf(stderr, "WARNING: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#define eprintf(args...) {fprintf(stderr, "ERROR: %s:%d in %s(): ", __FILE__, __LINE__, __FUNCTION__); fprintf(stderr, args);}
#ifdef DEBUGGING
#undef MAX
#undef MIN
#endif
#define MAX(a, b) ((a)>=(b)?(a):(b))
#define MIN(a, b) ((a)<=(b)?(a):(b))
