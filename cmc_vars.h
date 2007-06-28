/* -*- linux-c -*- */

struct {
	long N_MAX;
	long N_MAX_NEW;
	long N_STAR;
	long N_STAR_NEW;
	long N_BINARY;
} clus;

struct {
	double tot; /* total = kinetic + potential + internal + binary binding + central mass energy */
	double New; /* ??? */
	double ini; /* initial total energy */
	double K; /* total kinetic */
	double P; /* total potential */
	double Eint; /* total internal */
	double Eb; /* total binary binding energy */
} Etotal;

struct {
	double rh; /* half-mass radius */
} clusdyn;

/* tidal truncation stuff */
double max_r, Rtidal, TidalMassLoss, orbit_r;
double OldTidalMassLoss, DTidalMassLoss, Prev_Dt, Etidal; 
double clus_gal_mass_ratio, dist_gal_center;
/* core variables */
double N_core, Trc, rho_core, v_core, core_radius;
/* escaped stars */
double Eescaped, Jescaped, Eintescaped, Ebescaped;
/* total mass & co */
double initial_total_mass, Mtotal;
/******************* Input file parameters *************************/
long N_STAR_DIM, N_BIN_DIM, T_MAX_COUNT, MASS_PC_COUNT, STELLAR_EVOLUTION, SS_COLLISION, NO_MASS_BINS;
long SNAPSHOT_PERIOD, CHECKPOINT_PERIOD, MAX_WCLOCK_TIME, E_CONS;
int SNAPSHOT_CORE_COLLAPSE, SNAPSHOT_CORE_BOUNCE;
long IDUM, PERTURB, RELAXATION, MONITOR_COLL;
long NUM_CENTRAL_STARS;
double T_PRINT_STEP, T_MAX, THETASEMAX;
double TERMINAL_ENERGY_DISPLACEMENT, R_MAX;
double DT_FACTOR;
double MEGA_YEAR, SOLAR_MASS_DYN, METALLICITY, WIND_FACTOR;
double MMIN;
double MINIMUM_R;
double GAMMA, BININITEBMIN, BININITEBMAX;
struct CenMa cenma;
char MASS_PC[1000], MASS_BINS[1000], INPUT_FILE[1000];
int BINSINGLE, BINBIN;
int BININITKT, STOPATCORECOLLAPSE;
/* binary stuff */
long N_b, N_bb, N_bs;
double M_b, E_b;
binary_t *binary;
/* file pointers */
FILE *lagradfile, *dynfile, *logfile, *escfile, *snapfile, *ave_mass_file, *densities_file, *no_star_file, *centmass_file, **mlagradfile;
FILE *ke_rad_file, *ke_tan_file, *v2_rad_file, *v2_tan_file;
FILE *binaryfile, *binintfile, *collisionfile, *relaxationfile;
/* everything else except arrays */
char outprefix[100];
int se_file_counter;
long tcount;
long Echeck;
long snap_num, StepCount;
long newstarid;
double rho_core_single, rho_core_bin, rh_single, rh_binary;
double TotalTime, Dt;
double Sin2Beta;
int ReadSnapshot;
/* arrays */
star_t *star;
double *mass_pc, *mass_r, *ave_mass_r, *densities_r, *no_star_r;
double *ke_rad_r, *ke_tan_r, *v2_rad_r, *v2_tan_r;
double *mass_bins, **multi_mass_r;
int quiet;
/* debugging */
int debug;
/* units */
units_t units;
double madhoc;
/* etc. */
central_t central;
sigma_t sigma_array;
double Eoops; /* energy that has vanished from the system for various and sundry reasons */
double E_bb, E_bs, DE_bb, DE_bs;

/* variables for potential calculation (they are not the only ones, just the ones I added!) */
long last_index;
/* parameters for the Search_Grid */
long SEARCH_GRID;
long SG_STARSPERBIN, SG_MAXLENGTH, SG_MINLENGTH;
double SG_POWER_LAW_EXPONENT, SG_MATCH_AT_FRACTION, SG_PARTICLE_FRACTION;
/* The variable */
struct Search_Grid *r_grid;
/* for testing the effect */
long total_bisections;

/* an optional switch to turn off black hole accretion but having a non-zero central mass */
long BH_LOSS_CONE;
int CENTRAL_MASS_FROM_FILE;
double BH_R_DISRUPT_NB;
/* force to use the relaxation step, even if relaxation is turned off */
int FORCE_RLX_STEP;

#ifdef DEBUGGING
GHashTable *star_ids;
GArray *id_array;
#endif
