/* -*- linux-c -*- */

#include <zlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <fitsio.h>
#include "taus113-v2.h"
#include "fewbody.h"

#define VERSION "0.1"
#define NICK "Blah"
#define DATE "Mon Dec 20 14:23:25 CST 2004"

/* toggle stellar evolution */
/* #define SE */

/*************************** Parameters ******************************/
/* Large number, but still SF_INFINITY - 1 <> SF_INFINITY */
#define SF_INFINITY 1.0e10
/* Radius of the zeroth star so that 1/star[0].r is finite */
#define ZERO 1.0e-20

#define PI 3.14159265358979323

#define N_TRY 50000

#define MSUN 1.989e+33
#define SOLAR_MASS MSUN
#define RSUN 6.9599e10
#define G 6.67259e-8
#define KPC 3.0857e+21
#define YEAR 3.155693e+7
#define AU 1.496e+13
#define PARSEC 3.0857e+18

typedef struct{
	double m1; /* mass of star 1 */
	double m2; /* mass of star 2 */
	double a; /* semimajor axis */
	double e; /* eccentricity */
} binary_t;

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
	double r_peri; /* pericenter distance */
	double r_apo; /* apocenter distance */
	double phi; /* value of potential at position of star (only updated at end of timestep) */
	long   interacted; /* whether or not the star has undergone an interaction */
	long   binind; /* index to the binary */
	long id; 	/* the star's unique identifier */
	double rad; /* radius */
	/* stellar evolution variables */
#ifdef SE
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
#endif
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
	int MMIN;
	int BINBIN;
	int BINBIN_FEWBODY;
	int BINSINGLE;
	int BINSINGLE_FEWBODY;
	int CENTRAL_MASS;
	int DT_FACTOR;
	int DUMPS;
	int E_CONS;
	int IDUM;
	int INDEX_UNIT;
	int INPUT_FILE;
	int MASS_PC;
	int MAX_INDEX;
	int MEGA_YEAR;
	int METALLICITY;
	int MINIMUM_R;
	int MIN_LAGRANGIAN_RADIUS;
	int N_BINARY;
	int NUM_CORE_STARS;
	int ORIGINAL_PERTURB_STARS;
	int PERTURB;
	int R_MAX;
	int SIN2BETA_MAX;
	int SOLAR_MASS_DYN;
	int STELLAR_EVOLUTION;
	int TERMINAL_ENERGY_DISPLACEMENT;
	int T_MAX;
	int T_MAX_COUNT;
	int T_PRINT_STEP;
	int WIND_FACTOR;
	int GAMMA;
} parsed_t;

/* a struct containing the units used */
typedef struct{
	double t; /* time */
	double m; /* mass */
	double l; /* length */
} units_t;

/* a struct for the potential, must be malloc'ed */
typedef struct{
	long n; /* number of elements */
	double *r; /* radius */
	double *phi; /* potential */
} potential_t;

/********************** Function Declarations ************************/
inline double sqr(double x);
inline double cub(double x);
void setup_sub_time_step(void);
void sniff_stars(void);
void remove_star_center(long j);
void print_results(void);
double potential(double r);	       /* get potential using star.phi */
long potential_calculate(void);	/* calculate potential at star locations in star.phi */
void comp_mass_percent(void);
double get_positions(void);	/* get positions and velocities */
void perturb_stars(double Dt);	/* take a time step (perturb E,J) */
long FindZero_r(long x1, long x2, double r);
long FindZero_Q(long x1, long x2, double E, double J);
void ComputeEnergy(void);
void ComputeEnergy2(void);
void NormalizeEnergy(void);
void PrintLogOutput(void);
double GetTimeStep(void);
long CheckStop(void);
void ComputeIntermediateEnergy(void);
void RecomputeEnergy(void);
double rtbis(double (*func) (double), double x1, double x2, double xacc);
double trapzd(double (*func) (double), double a, double b, long n);
double qsimp(double (*func) (double), double a, double b);

void splint(double xa[], double ya[], double y2a[], long n, double x, double *y);
void spline(double x[], double y[], long n, double yp1, double ypn, double y2[]);
void print_2Dsnapshot(void);
void get_physical_units(void);
void update_vars(void);

/******* Binary routines **********/
void assign_binaries(void);
double bin_single_sigma(double y);

void print_version(FILE *stream);
void print_usage(FILE *stream, char *argv[]);
int parser(int argc, char *argv[], gsl_rng *r);
void PrintFileOutput(void);

/* Unmodified Numerical Recipes Routines */
void nrerror(char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);

/* fits stuff */
void read_fits_file_data(char *fits_file_name);
void read_fits_file_data_old(fitsfile *fptr);
void read_fits_file_data_new(fitsfile *fptr);
void read_fits_file_parameters_new(fitsfile *fptr, gsl_rng *r);
void read_fits_file_parameters_old(fitsfile *fptr, gsl_rng *r);
void read_fits_file_parameters(char *fits_file_name, gsl_rng *r);
void sshot_fits(gsl_rng *r);
void printerror(int status);

/* stellar evolution stuff */
void stellar_evolution_init(void);
void do_stellar_evolution(void);
void write_stellar_data(void);
double r_of_m(double M);

/* Fewbody stuff */
void bb_calcunits(fb_obj_t *obj[2], fb_units_t *units);
fb_ret_t binbin(double *t, double m00, double m01, double m10, double m11, double r00, double r01, double r10, double r11, double a0, double a1, double e0, double e1, double vinf, double b, fb_units_t *units, fb_hier_t *hier, gsl_rng *rng);
void bs_calcunits(fb_obj_t *obj[2], fb_units_t *units);
fb_ret_t binsingle(double *t, double m0, double m10, double m11, double r0, double r10, double r11, double a1, double e1, double vinf, double b, fb_units_t *units, fb_hier_t *hier, gsl_rng *rng);
long destroy_binary(long kp);
long create_star(void);
void perturb_stars_new(double dt, gsl_rng *rng);
void perturb_stars_fewbody(double dt, gsl_rng *rng);
void qsorts(star_t *s, long n);
void fits_sanity_check(void);
int move_to_hdu(fitsfile *fptr, char *sea_extname);
void mini_sshot(void);
void units_set(void);

/* signal/GSL error handling stuff */
void toggle_debugging(int signal);
void exit_cleanly(int signal);
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno);
void close_buffers(void);
void free_arrays(void);

/* if we are using gcc, we can remove the annoying warning */
#ifdef __GNUC__
static double SQR(double a) __attribute__ ((always_inline,pure));
static inline double SQR(double a){return a*a;}
#else
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#endif

/* macros */
#define function_Q(k, E, J) (2.0 * ((E) - star[(k)].phi) - SQR((J) / star[(k)].r))
#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#define wprintf(args...) {fprintf(stderr, "WARNING: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#define eprintf(args...) {fprintf(stderr, "ERROR: %s:%d in %s(): ", __FILE__, __LINE__, __FUNCTION__); fprintf(stderr, args);}
#define MAX(a, b) ((a)>=(b)?(a):(b))
#define MIN(a, b) ((a)<=(b)?(a):(b))
