/* vi: set filetype=c.doxygen: */

/* "extern" should be uncommented here when compiling with g++ */
#ifndef _MAIN_
#define _EXTERN_ /* extern */
#else
#define _EXTERN_
#endif
#include "hdf5.h"

_EXTERN_ clus_struct_t clus;
_EXTERN_ Etotal_struct_t Etotal;
_EXTERN_ clusdyn_struct_t clusdyn;
_EXTERN_ long RESTART_TCOUNT;
_EXTERN_ long NEW_IDUM;
//_EXTERN_ LightCollision LIGHTCOLLISION_TABLE[5000];
//_EXTERN_ int LIGHTCOLLISION_NRECORDS;
_EXTERN_ herr_t status;
//_EXTERN_ char *light_collision_field_names[NFIELDS_LIGHT_COLLISION];
//_EXTERN_ hid_t light_collision_field_type[NFIELDS_LIGHT_COLLISION];
//
///* Define an array of Particles */
//_EXTERN_ LightCollision light_collision_data[1];
///* Calculate the size and the offsets of our struct members in memory */
//_EXTERN_ size_t light_collision_dst_size;
//_EXTERN_ size_t light_collision_dst_offset[NFIELDS_LIGHT_COLLISION];
//
//_EXTERN_ size_t light_collision_dst_sizes[NFIELDS_LIGHT_COLLISION];


/* tidal truncation stuff */
_EXTERN_ double max_r, Rtidal, TidalMassLoss, orbit_r;
_EXTERN_ double OldTidalMassLoss, DTidalMassLoss, Prev_Dt, Etidal; 
_EXTERN_ double clus_gal_mass_ratio, dist_gal_center;
/* core variables */
_EXTERN_ double N_core, N_core_nb, Trc, rho_core, v_core, core_radius, rc_nb;
/* escaped stars */
_EXTERN_ double Eescaped, Jescaped, Eintescaped, Ebescaped;
/* total mass & co */
_EXTERN_ double initial_total_mass, Mtotal;
/* mass lost due to stellar evolution */
_EXTERN_ double DMse, DMrejuv;
/******************* Input file parameters *************************/
_EXTERN_ long N_STAR_DIM, N_STAR_DIM_OPT, N_BIN_DIM, N_BIN_DIM_OPT, T_MAX_COUNT, MASS_PC_COUNT, STELLAR_EVOLUTION, TIDAL_TREATMENT, SS_COLLISION, TIDAL_CAPTURE, BH_CAPTURE, TC_POLYTROPE, STAR_AGING_SCHEME, PREAGING, NO_MASS_BINS, NO_BSE_NATAL_KICK_ARRAY, NO_BSE_FPRIMC_ARRAY, NO_BSE_QCRIT_ARRAY;
_EXTERN_ long SNAPSHOT_DELTACOUNT, BH_SNAPSHOT_DELTACOUNT, PULSAR_DELTACOUNT, MAX_WCLOCK_TIME, CHECKPOINT_INTERVAL, CHECKPOINTS_TO_KEEP;
_EXTERN_ long headnode_time;
_EXTERN_ int SNAPSHOTTING, BH_SNAPSHOTTING, SNAPSHOT_CORE_COLLAPSE, SNAPSHOT_CORE_BOUNCE;
_EXTERN_ long IDUM, PERTURB, RELAXATION, TIDALLY_STRIP_STARS;
_EXTERN_ long NUM_CENTRAL_STARS;
_EXTERN_ double T_MAX, T_MAX_PHYS, THETASEMAX;
_EXTERN_ double TERMINAL_ENERGY_DISPLACEMENT, R_MAX;
_EXTERN_ double DT_FACTOR;
_EXTERN_ double MEGA_YEAR, SOLAR_MASS_DYN, METALLICITY, WIND_FACTOR;
_EXTERN_ double MMIN;
_EXTERN_ double MINIMUM_R;
_EXTERN_ double GAMMA, BININITEBMIN, BININITEBMAX;
_EXTERN_ double TC_FACTOR, COLL_FACTOR;
_EXTERN_ struct CenMa cenma;
_EXTERN_ char MASS_PC[1000], MASS_BINS[1000], INPUT_FILE[1000], BSE_QCRIT_ARRAY[1000], BSE_NATAL_KICK_ARRAY[1000], BSE_FPRIMC_ARRAY[1000];
_EXTERN_ char *TT_FILE;
_EXTERN_ char *DF_FILE;
_EXTERN_ char *SNAPSHOT_WINDOWS;
_EXTERN_ char *SNAPSHOT_WINDOW_UNITS;
_EXTERN_ int MASS_PC_BH_INCLUDE;
/**
* @brief Variable to store the input parameter which indicates the number of samples per processor to be used for Sample Sort. Defaults to the number of processors if not set.
*/
_EXTERN_ int SAMPLESIZE;
_EXTERN_ int BINSINGLE, BINBIN;
/**
* @brief Variable specified as command line argument or input parameter which indicates the number or streams (or number of processors to mimic) to be used in a serial run in order to mimic a parallel run.
*/
_EXTERN_ int STREAMS;
_EXTERN_ int BININITKT, STOPATCORECOLLAPSE;
/* Meagan - 3bb */
_EXTERN_ int THREEBODYBINARIES, ONLY_FORM_BH_THREEBODYBINARIES;;
_EXTERN_ double MIN_BINARY_HARDNESS;
_EXTERN_ double BINARY_DISTANCE_BREAKING;
_EXTERN_ int BINARY_BREAKING_MIN; 
_EXTERN_ int USE_TT_FILE; 
_EXTERN_ int USE_DF_CUTOFF, DF_INTEGRATED_CRITERION; 
_EXTERN_ double INITIAL_VALUE_DF_INTEGRAND;
_EXTERN_ long N3bbformed;
_EXTERN_ double delta_E_3bb;
// Meagan: extra bh output
_EXTERN_ long bhbinary, bhsingle, bhnonbh, bhbh, bh13, bh10, bh11, bh12, bhwd;
_EXTERN_ long bhstar, bh01, bh26, bh7, bh89;
_EXTERN_ long esc_bhsingle, esc_bhbinary, esc_bhnonbh, esc_bhbh, esc_bh13; 
_EXTERN_ long esc_bh10, esc_bh11, esc_bh12, esc_bhwd, esc_bhstar, esc_bh01, esc_bh26;
_EXTERN_ long esc_bh7, esc_bh89; 
_EXTERN_ long esc_bhsingle_tot, esc_bhbinary_tot, esc_bhnonbh_tot, esc_bhbh_tot;
_EXTERN_ long esc_bh13_tot, esc_bh10_tot, esc_bh11_tot, esc_bh12_tot, esc_bhwd_tot;
_EXTERN_ long esc_bhstar_tot, esc_bh01_tot, esc_bh26_tot, esc_bh7_tot, esc_bh89_tot;
_EXTERN_ double fb_bh, esc_fb_bh, esc_fb_bh_tot;
_EXTERN_ FILE *bhsummaryfile, *escbhsummaryfile, *newbhfile, *bhmergerfile;
/* BSE input file parameters */
_EXTERN_ int BSE_CEMERGEFLAG, BSE_CEKICKFLAG, BSE_CEHESTARFLAG, BSE_CEFLAG, BSE_TFLAG, BSE_IFFLAG, BSE_WDFLAG, BSE_RTMSFLAG, BSE_BHFLAG, BSE_BHSPINFLAG, BSE_REMNANTFLAG, BSE_IDUM, BSE_WINDFLAG, BSE_QCFLAG, BSE_EDDLIMFLAG, BSE_AIC, BSE_BDECAYFAC, BSE_HTPMB, BSE_ST_TIDE, BSE_ST_CR, BSE_REJUVFLAG, BSE_USSN, BSE_KICKFLAG, BSE_GRFLAG, BSE_BHMS_COLL_FLAG;
_EXTERN_ double BSE_POLAR_KICK_ANGLE, BH_RADIUS_MULTIPLYER, BSE_BHSPINMAG;
_EXTERN_ double BSE_PTS1, BSE_PTS2, BSE_PTS3, BSE_PTS1_HIGHMASS_CUTOFF, BSE_NETA, BSE_BWIND, BSE_HEWIND, BSE_ALPHA1, BSE_LAMBDAF, BSE_MXNS, BSE_BCONST, BSE_CK, BSE_REJUV_FAC, BSE_SIGMA, BSE_SIGMADIV, BSE_BHSIGMAFRAC, BSE_BETA, BSE_EDDFAC, BSE_GAMMA, BSE_XI, BSE_ACC2, BSE_PISN, BSE_EPSNOV, BSE_ECSN, BSE_ECSN_MLOW, BSE_REMBAR_MASSLOSS, BSE_ZSUN, BSE_DON_LIM, BSE_ACC_LIM;
/* binary stuff */
_EXTERN_ long N_b, N_bb, N_bs, last_hole;
//Probably not needed anymore
//_EXTERN_ long N_b_OLD, N_b_NEW;
_EXTERN_ double M_b, E_b;
/**
* @brief Array to store data of binaries
*/
_EXTERN_ binary_t *binary;
/* timer */
/**
* @brief Variable to store the input parameter which triggers the functionality to profile/time in detal, individual parts of the code.
*/
_EXTERN_ int TIMER;

/* file pointers */
_EXTERN_ FILE *lagradfile, *dynfile, *lagrad10file, *logfile, *escfile, *snapfile, *ave_mass_file, *densities_file, *no_star_file, *centmass_file, **mlagradfile;
_EXTERN_ FILE *ke_rad_file, *ke_tan_file, *v2_rad_file, *v2_tan_file;
_EXTERN_ FILE *binaryfile, *threebbfile, *threebbprobabilityfile, *lightcollisionfile, *threebbdebugfile, *binintfile, *collisionfile, *pulsarfile, *morepulsarfile, *newnsfile, *morecollfile, *triplefile, *tidalcapturefile, *semergedisruptfile, *removestarfile, *relaxationfile;
_EXTERN_ FILE *corefile;
_EXTERN_ FILE *fp_lagrad, *fp_log, *fp_denprof;
_EXTERN_ FILE *timerfile;
// Meagan: file for tracking potential fluctuations for innermost 1000 stars

/**
* @brief MPI: MPI-IO file pointers corresponding to the C File(pointer)s used in the serial code for files that are needed to be written out using MPI-IO
*/
_EXTERN_ MPI_File mpi_logfile, mpi_binintfile, mpi_escfile, mpi_collisionfile, mpi_pulsarfile, mpi_morepulsarfile, mpi_newnsfile, mpi_morecollfile, mpi_triplefile, mpi_tidalcapturefile, mpi_semergedisruptfile, mpi_removestarfile, mpi_relaxationfile;

/**
* @brief MPI: String buffers to store intermediate data that is finally flush out to files using MPI-IO
*/
_EXTERN_ char mpi_logfile_buf[STR_BUF_LEN], mpi_escfile_buf[STR_BUF_LEN], mpi_binintfile_buf[STR_BUF_LEN], mpi_collisionfile_buf[STR_BUF_LEN], mpi_pulsarfile_buf[STR_BUF_LEN], mpi_morepulsarfile_buf[STR_BUF_LEN], mpi_newnsfile_buf[STR_BUF_LEN], mpi_morecollfile_buf[STR_BUF_LEN], mpi_triplefile_buf[STR_BUF_LEN],mpi_tidalcapturefile_buf[STR_BUF_LEN], mpi_semergedisruptfile_buf[STR_BUF_LEN], mpi_removestarfile_buf[STR_BUF_LEN], mpi_relaxationfile_buf[STR_BUF_LEN];

/**
* @brief MPI: String buffers to store intermediate data that is finally flush out to files using MPI-IO
*/
_EXTERN_ char mpi_logfile_wrbuf[STR_WRBUF_LEN], mpi_escfile_wrbuf[STR_WRBUF_LEN], mpi_binintfile_wrbuf[STR_WRBUF_LEN], mpi_collisionfile_wrbuf[STR_WRBUF_LEN], mpi_pulsarfile_wrbuf[STR_WRBUF_LEN], mpi_morepulsarfile_wrbuf[STR_WRBUF_LEN], mpi_newnsfile_wrbuf[STR_WRBUF_LEN], mpi_morecollfile_wrbuf[STR_WRBUF_LEN], mpi_triplefile_wrbuf[STR_WRBUF_LEN], mpi_tidalcapturefile_wrbuf[STR_WRBUF_LEN], mpi_semergedisruptfile_wrbuf[STR_WRBUF_LEN], mpi_removestarfile_wrbuf[STR_WRBUF_LEN], mpi_relaxationfile_wrbuf[STR_WRBUF_LEN];

/**
* @brief MPI: Variables to maintail the length of the buffers until the next flush
*/
_EXTERN_ long long mpi_logfile_len, mpi_escfile_len, mpi_binintfile_len, mpi_collisionfile_len, mpi_pulsarfile_len, mpi_morepulsarfile_len, mpi_newnsfile_len, mpi_morecollfile_len, mpi_triplefile_len, mpi_tidalcapturefile_len, mpi_semergedisruptfile_len, mpi_removestarfile_len, mpi_relaxationfile_len;

/**
* @brief MPI: Variables to maintain the total offset of the file
*/
_EXTERN_ long long mpi_logfile_ofst_total, mpi_escfile_ofst_total, mpi_binaryfile_ofst_total, mpi_binintfile_ofst_total, mpi_collisionfile_ofst_total, mpi_pulsarfile_ofst_total, mpi_morepulsarfile_ofst_total, mpi_newnsfile_ofst_total, mpi_morecollfile_ofst_total, mpi_triplefile_ofst_total, mpi_tidalcapturefile_ofst_total, mpi_semergedisruptfile_ofst_total, mpi_removestarfile_ofst_total, mpi_relaxationfile_ofst_total;

/* Meagan's 3bb files */
/**
* @brief MPI: MPI-IO file pointers corresponding to the C File(pointer)s used in the serial code for files that are needed to be written out using MPI-IO
*/
_EXTERN_ MPI_File mpi_bhsummaryfile, mpi_escbhsummaryfile, mpi_newbhfile, mpi_bhmergerfile, mpi_threebbfile, mpi_threebbprobabilityfile, mpi_lightcollisionfile, mpi_threebbdebugfile;

/**
* @brief MPI: String buffers to store intermediate data that is finally flush out to files using MPI-IO
*/
_EXTERN_ char mpi_bhsummaryfile_buf[STR_BUF_LEN], mpi_escbhsummaryfile_buf[STR_BUF_LEN], mpi_newbhfile_buf[STR_BUF_LEN], mpi_bhmergerfile_buf[STR_BUF_LEN], mpi_threebbfile_buf[STR_BUF_LEN], mpi_threebbprobabilityfile_buf[STR_BUF_LEN], mpi_lightcollisionfile_buf[STR_BUF_LEN], mpi_threebbdebugfile_buf[STR_BUF_LEN];

/**
* @brief MPI: String buffers to store intermediate data that is finally flush out to files using MPI-IO
*/
_EXTERN_ char mpi_bhsummaryfile_wrbuf[STR_WRBUF_LEN], mpi_escbhsummaryfile_wrbuf[STR_WRBUF_LEN], mpi_newbhfile_wrbuf[STR_WRBUF_LEN], mpi_bhmergerfile_wrbuf[STR_WRBUF_LEN], mpi_threebbfile_wrbuf[STR_WRBUF_LEN], mpi_threebbprobabilityfile_wrbuf[STR_WRBUF_LEN], mpi_lightcollisionfile_wrbuf[STR_WRBUF_LEN], mpi_threebbdebugfile_wrbuf[STR_WRBUF_LEN];

/**
* @brief MPI: Variables to maintail the length of the buffers until the next flush
*/
_EXTERN_ long long mpi_bhsummaryfile_len, mpi_escbhsummaryfile_len, mpi_newbhfile_len, mpi_bhmergerfile_len, mpi_threebbfile_len, mpi_threebbprobabilityfile_len, mpi_lightcollisionfile_len, mpi_threebbdebugfile_len;

/**
* @brief MPI: Variables to maintain the total offset of the file
*/
_EXTERN_ long long mpi_bhsummaryfile_ofst_total, mpi_escbhsummaryfile_ofst_total, mpi_newbhfile_ofst_total, mpi_bhmergerfile_ofst_total, mpi_threebbfile_ofst_total, mpi_threebbprobabilityfile_ofst_total, mpi_lightcollisionfile_ofst_total, mpi_threebbdebugfile_ofst_total;


/* everything else except arrays */
_EXTERN_ char outprefix[100];
_EXTERN_ char oldoutprefix[100];
_EXTERN_ char dummystring[MAX_STRING_LENGTH], dummystring2[MAX_STRING_LENGTH], dummystring3[MAX_STRING_LENGTH], dummystring4[MAX_STRING_LENGTH];
_EXTERN_ int se_file_counter;
_EXTERN_ long tcount;
_EXTERN_ long NEXT_RESTART;
_EXTERN_ long next_restart_t;
_EXTERN_ long Echeck;
_EXTERN_ long snap_num, StepCount, bh_snap_num;
_EXTERN_ long newstarid;
_EXTERN_ double rho_core_single, rho_core_bin, rh_single, rh_binary;
_EXTERN_ double TotalTime, Dt;
_EXTERN_ double Sin2Beta;
/* arrays */
/**
* @brief Array to store data of all stars in the simulated system. Binaries are stored in a separate array name binary. If an element in the star array is a binary, it's binind property/variable is a non-zero value. Moreover, the value of this variable indicates the index of the binary array which holds the properties of this binary.
*/
_EXTERN_ star_t *star;
_EXTERN_ double *mass_pc, *mass_r, *ave_mass_r, *densities_r, *no_star_r;
_EXTERN_ double *ke_rad_r, *ke_tan_r, *v2_rad_r, *v2_tan_r;
_EXTERN_ double *mass_bins, *bse_qcrit_array, *bse_fprimc_array, *bse_natal_kick_array, **multi_mass_r;
_EXTERN_ int quiet;

/**
 * @brief the external arrays to save the diagionalized tidal tensor in.  Note 
 * that this is saved AFTER the Eigenvalues have been computed and the variable 
 * have been changed into code units.  
*/
_EXTERN_ double *TT_times, *TT_l1e;
_EXTERN_ long TT_num, TT_num_max;
/**
 * @brief the needed quantities for computing the dymical friction timescale
*/
_EXTERN_ double *DF_times, *DF_Menc, *DF_prefactor;
_EXTERN_ long DF_num, DF_num_max;
_EXTERN_ double t_df_cum, t_df_prev;

/* mpi parallelization stuff */
/**
* @brief MPI: Array to store the entire list of stellar positions. This is duplicated across processors in the beginning, and is synchronized at the end of each timestep.
*/
_EXTERN_ double *star_r;
/**
* @brief MPI: Array to store the entire list of stellar masses. This is duplicated across processors in the beginning, and is synchronized at the end of each timestep.
*/
_EXTERN_ double *star_m;
/**
* @brief MPI: Array to store the gravitaional potential at the stellar positions. This is duplicated across processors in the beginning, and is synchronized at the end of each timestep.
*/
_EXTERN_ double *star_phi;
/**
* @brief MPI: These variables are used for storing data partitioning related information in the parallel version. These variables store the start and end indices of the (hypothetical) global array that each processor is responsible for processing.
*/
_EXTERN_ int mpiBegin, mpiEnd;
_EXTERN_ int *mpiDisp, *mpiLen;
/**
* @brief MPI: There are some global variables that are updated at various places during a timestep, and towards the end need to be summed up across all processors. So, we store the values of these variables from the previous timestep into corresponding _old variables, and reset the actual variables to zero. At the end of the timestep, we cumulate/reduce the actual variables across processors and finally add them to the _old value i.e. total value of the variable from the previous timestep to obtain the updated values for these variables.
*/
_EXTERN_ double Eescaped_old, Jescaped_old, Eintescaped_old, Ebescaped_old, TidalMassLoss_old, Etidal_old;
/**
* @brief Processor id of THIS processor.
*/
_EXTERN_ int myid;
/**
* @brief Total number of processors or streams (in case of serial mimic).
*/
_EXTERN_ int procs;
_EXTERN_ double timeT, startTime, endTime; 
_EXTERN_ char funcName[64];
_EXTERN_ char fileTime[64];

/**
* @brief These variables are used for storing data partitioning related information in the parallel version. These arrays store the start and end indices in the global array that each processor is responsible for processing. In the serial version, these are used to mimic the parallel version to obtain comparable results.
*/
_EXTERN_ int *Start, *End; 
/**
* @brief This pointer is used in different ways in the serial and parallel versions. In the serial version, this is used to switch between different random states corresponding to the different streams used to mimic the parallel version, and draw an appropriate random number for each star when needed. In the parallel version this is THE random state of the respective processor.
*/
_EXTERN_ struct rng_t113_state *curr_st;
/**
* @brief Array used to mimic the parallel random number generator. This array stores the states of all the processors as in a corresponding parallel run.
*/
_EXTERN_ struct rng_t113_state *st;
/* to handle binaries */
/*
//Probably not needed anymore
_EXTERN_ binary_t *binary_buf;
_EXTERN_ int *num_bin_buf;
_EXTERN_ int size_bin_buf;
*/
_EXTERN_ int N_b_local;

_EXTERN_ FILE* ftest2;
_EXTERN_ char num2[5],filename2[20], tempstr2[20];
//Probably not needed anymore
//_EXTERN_ long *new_size;
//_EXTERN_ int *disp, *len;
//Temp file handle for debugging
FILE *ftest;
//MPI: Some variables to assist debugging
char num[5],filename[20], tempstr[20];
//MPI: Variables for detailed timing of the sorting routine and total communication
_EXTERN_ double t_sort_only, t_sort_lb;
_EXTERN_ double t_sort_lsort1, t_sort_splitters, t_sort_a2a, t_sort_lsort2, t_sort_oth;
_EXTERN_ double t_comm;

/* debugging */
_EXTERN_ int debug;
//Probably not needed anymore
//_EXTERN_ int mpi_debug;
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

/**
* @brief an optional switch to turn off black hole accretion but having a non-zero central mass
*/
_EXTERN_ long BH_LOSS_CONE;
_EXTERN_ double BH_R_DISRUPT_NB;
/**
* @brief force to use the relaxation step, even if relaxation is turned off
*/
_EXTERN_ int FORCE_RLX_STEP;
/**
* @brief calculate the binary interaction time steps by only considering hard binaries (0=off, 1=on)
*/
_EXTERN_ int DT_HARD_BINARIES;
/**
* @brief defines the minimum binary binding energy for a hard binary
*/
_EXTERN_ double HARD_BINARY_KT;
/**
* @brief
* A threshold value for the difference between apastron and periastron below
* which an orbit is considered to be circular. Currently this is only used
* for the period calculation in the loss-cone routine, as the integrator
* barfs when the difference between the two apsides (which form the integration
* interval) gets too small.
*/
_EXTERN_ double CIRC_PERIOD_THRESHOLD;
_EXTERN_ int WRITE_STELLAR_INFO;
_EXTERN_ int WRITE_BH_INFO;
_EXTERN_ int WRITE_RWALK_INFO;
_EXTERN_ int WRITE_EXTRA_CORE_INFO;
_EXTERN_ int WRITE_PULSAR_INFO;
_EXTERN_ int WRITE_MOREPULSAR_INFO;
_EXTERN_ int WRITE_MORECOLL_INFO;
_EXTERN_ int BHNS_TDE;
_EXTERN_ int CALCULATE10;

/**
* @brief
* scale the square of the deflection angle used for testing of entry into the loss
* cone by Trel*BH_LC_FDT/dt
*/
_EXTERN_ double BH_LC_FDT;
_EXTERN_ long AVEKERNEL;
_EXTERN_ long MIN_CHUNK_SIZE;
_EXTERN_ long BH_AVEKERNEL;
_EXTERN_ long N_Q_TRACE;
#ifdef DEBUGGING
_EXTERN_ GHashTable *star_ids;
_EXTERN_ GArray *id_array;
#endif
/**
* @brief The root solver to find the roots of Q
*/
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
_EXTERN_ double *snapshot_windows;
_EXTERN_ int *snapshot_window_counters;
_EXTERN_ int snapshot_window_count;
