/* -*- linux-c -*- */

struct {
	long N_MAX;
	long N_MAX_NEW;
	long N_STAR;
	long N_STAR_NEW;
	long N_BINARY;
} clus;

/* sub timestep stuff */
struct {
	long N_MAX, count, FACTOR;
	double totaltime, rmax;
} sub;

struct {
	double tot; /* total = K + P + Eint + cenma.E */
	double New; /* ??? */
	double ini; /* initial total energy */
	double K; /* total kinetic */
	double P; /* total potential */
	double Eint; /* total internal */
} Etotal;

/* tidal truncation stuff */
double max_r, Rtidal, TidalMassLoss, orbit_r;
double OldTidalMassLoss, DTidalMassLoss, Prev_Dt, Etidal; 
double clus_gal_mass_ratio, dist_gal_center;
/* core variables */
double N_core, Trc, rho_core, v_core, core_radius;
/* escaped stars */
double Eescaped, Ebescaped, Jescaped;
/* total mass & co */
double initial_total_mass, Mtotal;
/******************* Input file parameters *************************/
long N_STAR_DIM, T_MAX_COUNT, MASS_PC_COUNT, STELLAR_EVOLUTION;
long DUMPS, E_CONS, MAX_INDEX, INDEX_UNIT;
long IDUM, PERTURB;
long NUM_MASS, TOTAL_PARAMS, NUM_MASS_RADII_BINS, NUM_CORE_STARS;
double T_PRINT_STEP, T_MAX, SIN2BETA_MAX;
double TERMINAL_ENERGY_DISPLACEMENT, R_MAX;
double MIN_LAGRANGIAN_RADIUS, DT_FACTOR;
char MASS_PC[1000], INPUT_FILE[1000];
double MEGA_YEAR, SOLAR_MASS_DYN, METALLICITY, WIND_FACTOR;
double MMIN;
double MINIMUM_R;
double GAMMA;
struct CenMa cenma;
int BINSINGLE, BINBIN, BINSINGLE_FEWBODY, BINBIN_FEWBODY;
int ORIGINAL_PERTURB_STARS;
/* binary stuff */
long N_b, N_bb, N_bs;
double M_b, Delta_BE_bb, Delta_BE_bs, E_b, DE_bb, DE_bs;
binary_t *binary;
/* file pointers */
FILE *out[3], *logfile, *escfile, *snapfile, *ave_mass_file, *densities_file, *no_star_file, *centmass_file;
FILE *binaryfile, *binsinglefile, *binbinfile;
/* everything else except arrays */
char outprefix[100];
int se_file_counter;
long tcount;
long Echeck;
long snap_num, StepCount;
double rho_core_single, rho_core_bin, rh_single, rh_binary;
double TotalTime, Dt;
double Sin2Beta;
int ReadSnapshot;
/* arrays */
star_t *star;
double *mass_pc, *mass_r, *ave_mass_r, *densities_r, *no_star_r;
long *IndexTable;
/* debugging */
int debug;
/* units */
units_t units;
