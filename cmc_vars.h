/* -*- linux-c -*- */

/* "extern" should be uncommented here when compiling with g++ */
#ifndef _MAIN_
#define _EXTERN_ /* extern */
#else
#define _EXTERN_
#endif

_EXTERN_ clus_struct_t clus;
_EXTERN_ Etotal_struct_t Etotal;
_EXTERN_ clusdyn_struct_t clusdyn;

/* tidal truncation stuff */
_EXTERN_ double max_r, Rtidal, TidalMassLoss, orbit_r;
_EXTERN_ double OldTidalMassLoss, DTidalMassLoss, Prev_Dt, Etidal; 
_EXTERN_ double clus_gal_mass_ratio, dist_gal_center;
/* core variables */
_EXTERN_ double N_core, Trc, rho_core, v_core, core_radius, rc_nb;
/* escaped stars */
_EXTERN_ double Eescaped, Jescaped, Eintescaped, Ebescaped;
/* total mass & co */
_EXTERN_ double initial_total_mass, Mtotal;
/* mass lost due to stellar evolution */
_EXTERN_ double DMse, DMrejuv;
/******************* Input file parameters *************************/
_EXTERN_ long N_STAR_DIM, N_BIN_DIM, T_MAX_COUNT, MASS_PC_COUNT, STELLAR_EVOLUTION, TIDAL_TREATMENT, SS_COLLISION, STAR_AGING_SCHEME, PREAGING, NO_MASS_BINS;
_EXTERN_ long SNAPSHOT_DELTACOUNT, MAX_WCLOCK_TIME;
_EXTERN_ int SNAPSHOTTING, SNAPSHOT_CORE_COLLAPSE, SNAPSHOT_CORE_BOUNCE;
_EXTERN_ long IDUM, PERTURB, RELAXATION;
_EXTERN_ long NUM_CENTRAL_STARS;
_EXTERN_ double SNAPSHOT_DELTAT, T_MAX, T_MAX_PHYS, THETASEMAX;
_EXTERN_ double TERMINAL_ENERGY_DISPLACEMENT, R_MAX;
_EXTERN_ double DT_FACTOR;
_EXTERN_ double MEGA_YEAR, SOLAR_MASS_DYN, METALLICITY, WIND_FACTOR;
_EXTERN_ double MMIN;
_EXTERN_ double MINIMUM_R;
_EXTERN_ double GAMMA, BININITEBMIN, BININITEBMAX;
_EXTERN_ struct CenMa cenma;
_EXTERN_ char MASS_PC[1000], MASS_BINS[1000], INPUT_FILE[1000];
_EXTERN_ int MASS_PC_BH_INCLUDE;
_EXTERN_ int BINSINGLE, BINBIN;
_EXTERN_ int BININITKT, STOPATCORECOLLAPSE;
/* binary stuff */
_EXTERN_ long N_b, N_bb, N_bs;
_EXTERN_ double M_b, E_b;
_EXTERN_ binary_t *binary;
/* file pointers */
_EXTERN_ FILE *lagradfile, *dynfile, *logfile, *escfile, *snapfile, *ave_mass_file, *densities_file, *no_star_file, *centmass_file, **mlagradfile;
_EXTERN_ FILE *ke_rad_file, *ke_tan_file, *v2_rad_file, *v2_tan_file;
_EXTERN_ FILE *binaryfile, *binintfile, *collisionfile, *semergedisruptfile, *removestarfile, *relaxationfile;
_EXTERN_ FILE *corefile;
/* everything else except arrays */
_EXTERN_ char outprefix[100];
_EXTERN_ int se_file_counter;
_EXTERN_ long tcount;
_EXTERN_ long Echeck;
_EXTERN_ long snap_num, StepCount;
_EXTERN_ long newstarid;
_EXTERN_ double rho_core_single, rho_core_bin, rh_single, rh_binary;
_EXTERN_ double TotalTime, Dt;
_EXTERN_ double Sin2Beta;
/* arrays */
_EXTERN_ star_t *star;
_EXTERN_ double *mass_pc, *mass_r, *ave_mass_r, *densities_r, *no_star_r;
_EXTERN_ double *ke_rad_r, *ke_tan_r, *v2_rad_r, *v2_tan_r;
_EXTERN_ double *mass_bins, **multi_mass_r;
_EXTERN_ int quiet;
/* debugging */
_EXTERN_ int debug;
/* units */
_EXTERN_ units_t units;
_EXTERN_ double madhoc;
/* etc. */
_EXTERN_ central_t central;
_EXTERN_ sigma_t sigma_array;
_EXTERN_ double Eoops; /* energy that has vanished from the system for various and sundry reasons */
_EXTERN_ double E_bb, E_bs, DE_bb, DE_bs;
/* FITS stuff */
_EXTERN_ cmc_fits_data_t cfd;
/* variables for potential calculation (they are not the only ones, just the ones I added!) */
_EXTERN_ long last_index;
/* parameters for the Search_Grid */
_EXTERN_ long SEARCH_GRID;
_EXTERN_ long SG_STARSPERBIN, SG_MAXLENGTH, SG_MINLENGTH;
_EXTERN_ double SG_POWER_LAW_EXPONENT, SG_MATCH_AT_FRACTION, SG_PARTICLE_FRACTION;
/* The variable */
_EXTERN_ struct Search_Grid *r_grid;
/* for testing the effect */
_EXTERN_ long total_bisections;

/* an optional switch to turn off black hole accretion but having a non-zero central mass */
_EXTERN_ long BH_LOSS_CONE;
_EXTERN_ double BH_R_DISRUPT_NB;
/* force to use the relaxation step, even if relaxation is turned off */
_EXTERN_ int FORCE_RLX_STEP;
/* A threshold value for the difference between apastron and periastron below
 * which an orbit is considered to be circular. Currently this is only used 
 * for the period calculation in the loss-cone routine, as the integrator 
 * barfs when the difference between the two apsides (which form the integration
 * interval) gets too small.
 */
_EXTERN_ double CIRC_PERIOD_THRESHOLD;
_EXTERN_ int WRITE_STELLAR_INFO;
_EXTERN_ int WRITE_RWALK_INFO;
_EXTERN_ int WRITE_EXTRA_CORE_INFO;
#ifdef EXPERIMENTAL
/* scale the square of the deflection angle used for testing of entry into the loss 
 * cone by Trel*BH_LC_FDT/dt 
 */
_EXTERN_ double BH_LC_FDT;
_EXTERN_ long AVEKERNEL;
#endif
_EXTERN_ long N_Q_TRACE;
#ifdef DEBUGGING
_EXTERN_ GHashTable *star_ids;
_EXTERN_ GArray *id_array;
#endif
/* The root solver to find the roots of Q */
_EXTERN_ gsl_root_fsolver *q_root;
_EXTERN_ double APSIDES_PRECISION;
_EXTERN_ long APSIDES_MAX_ITER;
_EXTERN_ double APSIDES_CONVERGENCE;
_EXTERN_ double *zpars;
_EXTERN_ struct core_t no_remnants;
_EXTERN_ double OVERWRITE_RVIR;
_EXTERN_ double OVERWRITE_RTID;
_EXTERN_ double OVERWRITE_Z;
_EXTERN_ double OVERWRITE_MCLUS;

