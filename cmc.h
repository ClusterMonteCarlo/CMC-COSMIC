/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <sys/times.h>
#include <zlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <fitsio.h>
#include <string.h>
#include "common/fitslib.h"
#include "common/taus113-v2.h"
#include "fewbody-0.24/fewbody.h"
#include "cmc_core.h"
#include "config.h"

#ifdef USE_MPI
#include "cmc_mpi.h"
#endif

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

#ifndef EXPERIMENTAL
#define AVEKERNEL 20
#endif

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
/* We'll set XCOLL to 4 here, as a generous estimate of the tidal capture + merger cross section.  XCOLL is used
   only to estimate the timestep, so it doesn't have to be commpletely accurate.  The real cross section is still
   used when collisions are checked for.  */
#define XCOLLTC 4.0 /* tidal capture */
#define XCOLLSS 1.5 /* simple sticky spheres */
/* rperi = XBS * a */
#define XBS 2.0
/* rperi = XBB * (a_1 + a_2) */
#define XBB 2.0
/* hard soft boundary XHS=vorb/W */
#define XHS 0.7
/* safety factor used in the black hole accretion routine */
/* by Freitag & Benz (2002) */
#define CSAFE 0.2

/* extra star type */
#define NOT_A_STAR (-100)

/*Sourav: required for the preaging when STAR_AGING_SCHEME is emabled in Solar mass*/
#define PREAGING_MASS 5.0

#define MAX_STRING_LENGTH 2048

//MPI: For MPI-IO
#define STR_BUF_LEN 10000
#define STR_WRBUF_LEN 100000000

/*-------------------------------------------------------------c
*
*     STELLAR TYPES - KW
*
*        0 - deeply or fully convective low mass MS star
*        1 - Main Sequence star
*        2 - Hertzsprung Gap
*        3 - First Giant Branch
*        4 - Core Helium Burning
*        5 - First Asymptotic Giant Branch
*        6 - Second Asymptotic Giant Branch
*        7 - Main Sequence Naked Helium star
*        8 - Hertzsprung Gap Naked Helium star
*        9 - Giant Branch Naked Helium star
*       10 - Helium White Dwarf
*       11 - Carbon/Oxygen White Dwarf
*       12 - Oxygen/Neon White Dwarf
*       13 - Neutron Star
*       14 - Black Hole
*       15 - Massless Supernova
*
*-------------------------------------------------------------*/


/**
* @brief Data structure for storing binary properties.
* @details One or more stars in the star data structure can be binary stars. To store all the binaries, an array of this structure is used. A binary in the star data structure has a non-zero value for the variable binind. Moreover, the value of binind indicates the index in the binary array which holds the properties of that binary. For example, if star[213].binind has the value 56, binary[56] contains the properties of that binary.
*/
typedef struct{
/**
* @brief  unique id of star 1
*/
	long id1;
/**
* @brief  unique id of star 2
*/
	long id2;
/**
* @brief  radius of star 1
*/
	double rad1;
/**
* @brief  radius of star 2
*/
	double rad2;
/**
* @brief  mass of star 1
*/
	double m1;
/**
* @brief  mass of star 2
*/
	double m2;
/**
* @brief  internal energy of star 1
*/
	double Eint1;
/**
* @brief  internal energy of star 2
*/
	double Eint2;
/**
* @brief  semimajor axis
*/
	double a;
/**
* @brief  eccentricity
*/
	double e;
/**
* @brief  whether or not binary exists
*/
	int inuse;
/**
* @brief  star types
*/
	int bse_kw[2];
/**
* @brief original (t=0) ZAMS mass 
*/
	double bse_zams_mass[2];
/**
* @brief  initial masses before BSE is run
*/
	double bse_mass0[2];
/**
* @brief  masses
*/
	double bse_mass[2];
/**
* @brief  radii
*/
	double bse_radius[2];
/**
* @brief  luminosity
*/
	double bse_lum[2];
/**
* @brief ?
*/
	double bse_massc[2];
/**
* @brief ?
*/
	double bse_radc[2];
/**
* @brief ?
*/
	double bse_menv[2];
/**
* @brief ?
*/
	double bse_renv[2];
/**
* @brief  original spin
*/
	double bse_ospin[2];
/**
* @brief  Pulsar magnetic field
*/
	double bse_B_0[2];
/**
* @brief  Amount of mass pulsar has accreted
*/
	double bse_bacc[2];
/**
* @brief  Amount of time pulsar has spent accreting
*/
	double bse_tacc[2];
/**
* @brief ?
*/
	double bse_epoch[2];
/**
* @brief ?
*/
	double bse_tms[2];
/**
* @brief  physical time
*/
	double bse_tphys;
/**
* @brief  binary orbital period
*/
	double bse_tb;
/**
* @brief  mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)]
*/
	double bse_bcm_dmdt[2];
/**
* @brief  radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)]
*/
	double bse_bcm_radrol[2];
/**
* @brief  Pulsar magnetic field strength at surface
*/
	double bse_bcm_B[2];
/**
* @brief  provides formation pathway of NS
*/
	double bse_bcm_formation[2];
	//Sourav:toy rejuvenation variables
/**
* @brief  Sourav: lifetime of star1
*/
	double lifetime_m1;
/**
* @brief  Sourav: lifetime of star2
*/
	double lifetime_m2;
/**
* @brief  Sourav: createtime of star1
*/
	double createtime_m1;
/**
* @brief  Sourav: createtime of star2
*/
	double createtime_m2;
	// Meagan: to keep track of three-body binaries
/**
* @brief  whether binary was formed via three-body encounter
*/
	//int threebodybinary;
} binary_t;

/**
* @brief ?
*/
struct star_coords {
  long index;
  long field_index;
  double r;
  double vr;
  double vt;
  double E, J;
  double pot;
};

/**
* @brief Data structure for storing single star properties.
* @details An array of this structure serves as the principal data structure to store the data of the cluster being simulated by CMC. One or more stars in the star array can be binary stars.
*/
typedef struct{
	/* dynamical evolution variables */
/**
* @brief   radial coordinate
*/
	double r;
/**
* @brief  radial velocity
*/
	double vr;
/**
* @brief  tangential velocity
*/
	double vt;
/**
* @brief  mass
*/
	double m;
/**
* @brief  kinetic plus potential energy per unit mass
*/
	double E;
/**
* @brief  angular momentum per unit mass
*/
	double J;
/**
* @brief  "intermediate" energy per unit mass
*/
	double EI;
/**
* @brief  internal energy (due to collisions, e.g.)
*/
	double Eint;
/**
* @brief  new radial coordinate
*/
	double rnew;
/**
* @brief  new radial velocity
*/
	double vrnew;
/**
* @brief  new tangential velocity
*/
	double vtnew;
/**
* @brief  ?
*/
	double rOld;
/**
* @brief  ?
*/
	double X;
/**
* @brief  random variable that must be stored
*/
	double Y;
/**
* @brief  pericenter distance
*/
	double r_peri;
/**
* @brief  apocenter distance
*/
	double r_apo;
/**
* @brief  value of potential at position of star (only updated at end of timestep)
*/
	double phi;
/**
* @brief  whether or not the star has undergone a strong interaction (i.e., not relaxation)
*/
	long   interacted;
/**
* @brief whether or not object was involved in three-body binary formation
*/
	long   threebb_interacted;
/**
* @brief  index to the binary
* @details If the star is a binary, this variable has a non-zero value. Moreover, the value of this variable indicates the index of the binary array which holds the properties of this binary. For example, if star[213].binind has the value 56, binary[56] contains the properties of that binary.
*/
	long   binind;
/**
* @brief  the star's unique identifier
*/
	long   id;
/**
* @brief  radius
*/
	double rad;
/**
* @brief  variables for Stodolkiewicz
*/
	double Uoldrold, Uoldrnew;
/**
* @brief  energy conservation scheme
*/
	double vtold, vrold;
	/* TODO: stellar evolution variables */
/**
* @brief mass at beining of (BSE) integration 
*/
	double se_mass;
/**
* @brief original ZAMS mass (at t=0) 
*/
	double zams_mass;
/**
* @brief stellar types (see bse_wrap/bse/bse.f for the list)
*/
	int se_k;
/**
* @brief ?
*/
	double se_mt;
/**
* @brief ?
*/
	double se_ospin;
/**
* @brief  Pulsar initial magentif field
*/
	double se_B_0;
/**
* @brief ?
*/
	double se_bacc;
/**
* @brief ?
*/
	double se_tacc;
/**
* @brief ?
*/
	double se_epoch;
/**
* @brief ?
*/
	double se_tphys;
/**
* @brief ?
*/
	double se_radius;
/**
* @brief ?
*/
	double se_lum;
/**
* @brief ?
*/
	double se_mc;
/**
* @brief ?
*/
	double se_rc;
/**
* @brief ?
*/
	double se_menv;
/**
* @brief ?
*/
	double se_renv;
/**
* @brief ?
*/
	double se_tms;
/**
* @brief  Pulsar surface magnetic field
*/
	double se_scm_B;
/**
* @brief  formation pathway of NS
*/
	double se_scm_formation;
/**
* @brief Sourav: toy rejuvenation variables
*/
	double createtime, createtimenew, createtimeold;
/**
* @brief ?
*/
	double lifetime, lifetimeold, lifetimenew;
} star_t;

/**
* @brief ?
*/
struct CenMa{
/**
* @brief ?
*/
	double m;
/**
* @brief ?
*/
	double m_new;
/**
* @brief ?
*/
	double E;
};

/**
* @brief ?
*/
struct get_pos_str {
/**
* @brief ?
*/
	double max_rad;
/**
* @brief ?
*/
	double phi_rtidal;
/**
* @brief ?
*/
	double phi_zero;
/**
* @brief ?
*/
	long N_LIMIT;
/**
* @brief ?
*/
	int taskid;
/**
* @brief ?
*/
	struct CenMa CMincr;
/**
* @brief ?
*/
	gsl_rng *thr_rng;
};

// This is a total hack for including parameter documentation
/**
* @brief Struct to store the input parameters parsed from the input .cmc file
*/
typedef struct{
#define PARAMDOC_SAMPLESIZE "no.of sample keys contributed per processor for sample sort--binary interactions. Applicable only for the parallel version"
/**
* @brief no.of sample keys contributed per processor for sample sort. By default (if not specified), equal to the number of processors. Applicable only for the parallel version
* @cite [Pattabiraman et al.(2013)]{2013ApJS..204...15P} Pattabiraman, B., Umbreit, S., Liao, W.-k., et al.\ 2013, \apjs, 204, 15
*/
	int SAMPLESIZE;
#define PARAMDOC_BINBIN "toggles binary--binary interactions (0=off, 1=on)"
/**
* @brief toggles binary--binary interactions (0=off, 1=on)
*/
	int BINBIN;
#define PARAMDOC_BINSINGLE "toggles binary--single interactions (0=off, 1=on)"
/**
* @brief toggles binary--single interactions (0=off, 1=on)
*/
	int BINSINGLE;
#define PARAMDOC_STREAMS "to run the serial version with the given number of random streams - primarily used to mimic the parallel version running with the same no.of processors"
	int STREAMS;
/* Meagan - 3bb */
#define PARAMDOC_THREEBODYBINARIES "toggles three-body binary formation (0=off, 1=on)"
/**
* @brief toggles three-body binary formation (0=off, 1=on)
*/
	int THREEBODYBINARIES;
#define PARAMDOC_MIN_BINARY_HARDNESS "minimum hardness for newly formed three-body binaries"
/**
* @brief minimum hardness for newly formed three-body binaries
*/
	int MIN_BINARY_HARDNESS;
#define PARAMDOC_BINARY_DISTANCE_BREAKING "What fraction of the local interparticle seperation do binaries have to be apart to be broken up (default is 0.1)"
/**
* @brief Fraction of the local interparticule seperation to disrupt binaries (default is 0.1)
*/
	int BINARY_DISTANCE_BREAKING;
#define PARAMDOC_BINARY_BREAKING_MIN "Should I use the interparticle seperation or MIN_BINARY_HARDNESS to break soft binaries (0=interparticle, 1=MIN_BINARY_HARDNESS)"
/**
* @brief criterion to break soft binaries (0=interparticle, 1=MIN_BINARY_HARDNESS) 
*/
	int BINARY_BREAKING_MIN;
#define PARAMDOC_ONLY_FORM_BH_THREEBODYBINARIES "allow only black holes to form binaries via three-body binary formation (1=only black holes, 0=any object types)"
/**
* @brief allow only black holes to form binaries via three-body binary formation (1=only black holes, 0=any object types)
*/
	int ONLY_FORM_BH_THREEBODYBINARIES;
#define PARAMDOC_BH_SNAPSHOTTING "toggles output bh snapshotting (0=off, 1=on)"
/**
* @brief toggles output bh snapshotting (0=off, 1=on)
*/
	int BH_SNAPSHOTTING;
#define PARAMDOC_BH_SNAPSHOT_DELTACOUNT "BH snapshotting interval in time steps"
/**
* @brief BH snapshotting interval in time steps
*/
	int BH_SNAPSHOT_DELTACOUNT;
#define PARAMDOC_SNAPSHOTTING "toggles output snapshotting (0=off, 1=on)"
/**
* @brief toggles output snapshotting (0=off, 1=on)
*/
	int SNAPSHOTTING;
#define PARAMDOC_SNAPSHOT_DELTAT "snapshotting time interval (FP units)"
/**
* @brief snapshotting time interval (FP units)
*/
	int SNAPSHOT_DELTAT;
#define PARAMDOC_SNAPSHOT_DELTACOUNT "snapshotting interval in time steps"
/**
* @brief snapshotting interval in time steps
*/
	int SNAPSHOT_DELTACOUNT;
#define PARAMDOC_SNAPSHOT_CORE_COLLAPSE "output extra snapshotting information during core collapse (0=off, 1=on)"
/**
* @brief output extra snapshotting information during core collapse (0=off, 1=on)
*/
        int SNAPSHOT_CORE_COLLAPSE;
#define PARAMDOC_SNAPSHOT_CORE_BOUNCE "output extra snapshotting information during core bounce (0=off, 1=on)"
/**
* @brief output extra snapshotting information during core bounce (0=off, 1=on)
*/
        int SNAPSHOT_CORE_BOUNCE;
#define PARAMDOC_SNAPSHOT_WINDOWS "Output extra snapshots within time windows. \n#The format is start_w0,step_w0,end_w0;start_w1,step_w1,stop_w1 ... etc." 
/**
* @brief Output extra snapshots within time windows. The format is start_w0,step_w0,end_w0;start_w1,step_w1,stop_w1 ... etc.
*/
        int SNAPSHOT_WINDOWS;
#define PARAMDOC_SNAPSHOT_WINDOW_UNITS "Units used for time window parameters. Possible choices: Gyr, Trel, and Tcr"
/**
* @brief Units used for time window parameters. Possible choices: Gyr, Trel, and Tcr
*/
        int SNAPSHOT_WINDOW_UNITS;
#define PARAMDOC_IDUM "random number generator seed"
/**
* @brief random number generator seed
*/
	int IDUM;
#define PARAMDOC_INPUT_FILE "input FITS file"
/**
* @brief input FITS file
*/
	int INPUT_FILE;
#define PARAMDOC_MASS_PC "mass fractions for Lagrange radii"
/**
* @brief mass fractions for Lagrange radii
*/
	int MASS_PC;
#define PARAMDOC_MASS_PC_BH_INCLUDE "Shall the central black hole be included for the calculation of the Lagrange radii? (1=yes, 0=no)"
/**
* @brief Shall the central black hole be included for the calculation of the Lagrange radii? (1=yes, 0=no)
*/
	int MASS_PC_BH_INCLUDE;
#define PARAMDOC_MASS_BINS "mass ranges for calculating derived quantities"
/**
* @brief mass ranges for calculating derived quantities
*/
	int MASS_BINS;
#define PARAMDOC_MINIMUM_R "radius of central mass"
/**
* @brief radius of central mass
*/
	int MINIMUM_R;
#define PARAMDOC_STOPATCORECOLLAPSE "stop calculation at core collapse (0=no, 1=yes)"
/**
* @brief stop calculation at core collapse (0=no, 1=yes)
*/
	int STOPATCORECOLLAPSE;
#define PARAMDOC_NUM_CENTRAL_STARS "number of central stars used to calculate certain averages"
/**
* @brief number of central stars used to calculate certain averages
*/
	int NUM_CENTRAL_STARS;
#define PARAMDOC_PERTURB "perform dynamical perturbations on objects (0=off, 1=on)"
/**
* @brief perform dynamical perturbations on objects (0=off, 1=on)
*/
	int PERTURB;
#define PARAMDOC_RELAXATION "perform two-body relaxation (0=off, 1=on)"
/**
* @brief perform two-body relaxation (0=off, 1=on)
*/
	int RELAXATION;
#define PARAMDOC_TIDALLY_STRIP_STARS "remove stars that pass beyond the tidal boundary"
/**
* @brief strip stars pass the tidal boundary (0 = off, 1=on[default]) 
*/
	int TIDALLY_STRIP_STARS;
#define PARAMDOC_THETASEMAX "maximum super-encounter scattering angle (radians)"
/**
* @brief maximum super-encounter scattering angle (radians)
*/
	int THETASEMAX;
#define PARAMDOC_STELLAR_EVOLUTION "stellar evolution (0=off, 1=on)"
/**
* @brief stellar evolution (0=off, 1=on)
*/
	int STELLAR_EVOLUTION;
#define PARAMDOC_TIDAL_TREATMENT "choose the tidal cut-off criteria (0=radial criteria, 1=Giersz energy criteria)"
/**
* @brief choose the tidal cut-off criteria (0=radial criteria, 1=Giersz energy criteria)
*/
	int TIDAL_TREATMENT;
#define PARAMDOC_SS_COLLISION "perform physical stellar collisions (0=off, 1=on)"
/**
* @brief perform physical stellar collisions (0=off, 1=on)
*/
	int SS_COLLISION;
#define PARAMDOC_TIDAL_CAPTURE "allow for tidal capture in single-single interactions, including Lombardi, et al. (2006) collisional binary formation mechanism (0=off, 1=on)"
/**
* @brief allow for tidal capture in single-single interactions, including Lombardi, et al. (2006) collisional binary formation mechanism (0=off, 1=on)
*/
	int TIDAL_CAPTURE;
	//Sourav: toy rejuvenation flags
#define PARAMDOC_STAR_AGING_SCHEME "the aging scheme of the stars (0=infinite age of all stars, 1=rejuvenation, 2=zero lifetime of collision stars, 3=arbitrary lifetime)"
/**
* @brief the aging scheme of the stars (0=infinite age of all stars, 1=rejuvenation, 2=zero lifetime of collision stars, 3=arbitrary lifetime)
*/
	int STAR_AGING_SCHEME;
#define PARAMDOC_PREAGING "preage the cluster (0=off, 1=on)"
/**
* @brief preage the cluster (0=off, 1=on)
*/
	int PREAGING;
#define PARAMDOC_TERMINAL_ENERGY_DISPLACEMENT "energy change calculation stopping criterion"
/**
* @brief energy change calculation stopping criterion
*/
	int TERMINAL_ENERGY_DISPLACEMENT;
#define PARAMDOC_T_MAX "maximum integration time (FP units)"
/**
* @brief maximum integration time (FP units)
*/
	int T_MAX;
#define PARAMDOC_T_MAX_PHYS "maximum integration time (Gyr)"
/**
* @brief maximum integration time (Gyr)
*/
	int T_MAX_PHYS;
#define PARAMDOC_T_MAX_COUNT "maximum number of time steps"
/**
* @brief maximum number of time steps
*/
	int T_MAX_COUNT;
#define PARAMDOC_MAX_WCLOCK_TIME "maximum wall clock time (seconds)"
/**
* @brief maximum wall clock time (seconds)
*/
	int MAX_WCLOCK_TIME;
#define PARAMDOC_WIND_FACTOR "stellar evolution wind mass loss factor (0.5-2)"
/**
* @brief stellar evolution wind mass loss factor (0.5-2)
*/
	int WIND_FACTOR;
#define PARAMDOC_GAMMA "gamma in Coulomb logarithm"
/**
* @brief gamma in Coulomb logarithm
*/
	int GAMMA;
#define PARAMDOC_SEARCH_GRID "search grid (0=off, 1=on)"
/**
* @brief search grid (0=off, 1=on)
*/
        int SEARCH_GRID;
#define PARAMDOC_SG_STARSPERBIN "number of stars that should ideally be in each search bin"
/**
* @brief number of stars that should ideally be in each search bin
*/
        int SG_STARSPERBIN;
#define PARAMDOC_SG_MAXLENGTH "maximum length of the search grid"
/**
* @brief maximum length of the search grid
*/
        int SG_MAXLENGTH;
#define PARAMDOC_SG_MINLENGTH "minimum length of the search grid"
/**
* @brief minimum length of the search grid
*/
        int SG_MINLENGTH;
#define PARAMDOC_SG_POWER_LAW_EXPONENT "slope of the assumed power-law for r(N), where N is the number of stars within r. (0.5 hard coded)"
/**
* @brief slope of the assumed power-law for r(N), where N is the number of stars within r. (0.5 hard coded)
*/
        int SG_POWER_LAW_EXPONENT;
#define PARAMDOC_SG_MATCH_AT_FRACTION "fraction frac that adjusts the constant factor in the power-law for r(N) such that r_pl(frac*N_tot)=r(frac*N_tot) (0.5)"
/**
* @brief fraction frac that adjusts the constant factor in the power-law for r(N) such that r_pl(frac*N_tot)=r(frac*N_tot) (0.5)
*/
        int SG_MATCH_AT_FRACTION;
#define PARAMDOC_SG_PARTICLE_FRACTION "frac_p that defines the maximum Np= frac_p*N_tot for which r(N<Np) can be reasonably approximated as a power-law (0.95)"
/**
* @brief frac_p that defines the maximum Np= frac_p*N_tot for which r(N<Np) can be reasonably approximated as a power-law (0.95)
*/
        int SG_PARTICLE_FRACTION;
#define PARAMDOC_BH_LOSS_CONE "perform loss-cone physics for central black hole (0=off, 1=on)"
/**
* @brief perform loss-cone physics for central black hole (0=off, 1=on)
*/
        int BH_LOSS_CONE;
#define PARAMDOC_BH_R_DISRUPT_NB "central black hole disruption radius (N-body units)"
/**
* @brief central black hole disruption radius (N-body units)
*/
        int BH_R_DISRUPT_NB;
#define PARAMDOC_FORCE_RLX_STEP "force a relaxation step (useful when RELAXATION=0) (0=off, 1=on)"
/**
* @brief force a relaxation step (useful when RELAXATION=0) (0=off, 1=on)
*/
        int FORCE_RLX_STEP; 
#define PARAMDOC_DT_HARD_BINARIES "calculate the binary interaction time steps by only considering hard binaries (0=off, 1=on)"
/**
* @brief calculate the binary interaction time steps by only considering hard binaries (0=off, 1=on)
*/
        int DT_HARD_BINARIES; 
#define PARAMDOC_HARD_BINARY_KT "The minimum binary binding energy (in units of kT) for a binary to be considered 'hard' for the time step calculation."
/**
* @brief The minimum binary binding energy (in units of kT) for a binary to be considered 'hard' for the time step calculation.
*/
        int HARD_BINARY_KT; 
#ifdef EXPERIMENTAL
#define PARAMDOC_BH_LC_FDT "sub time step size on which the code tries to approximately advance particles that have a MC time step larger than BH_LC_FDT times the local relaxation time. In some versions of the code the particles are literally advanced while in other a simple scaling is used. None of them really work. (0)"
/**
* @brief sub time step size on which the code tries to approximately advance particles that have a MC time step larger than BH_LC_FDT times the local relaxation time. In some versions of the code the particles are literally advanced while in other a simple scaling is used. None of them really work. (0)
*/
        int BH_LC_FDT;
#define PARAMDOC_AVEKERNEL "one half the number of stars over which to average certain quantities"
/**
* @brief one half the number of stars over which to average certain quantities
*/
        int AVEKERNEL;
#endif
#define PARAMDOC_MIN_CHUNK_SIZE "minimum size of chunks that get partitioned across processors in the parallel code"
/**
* @brief minimum size of chunks that get partitioned across processors in the parallel code
*/
        int MIN_CHUNK_SIZE;
#define PARAMDOC_BH_AVEKERNEL "one half the number of stars over which to average 3bb related quantities"
        int BH_AVEKERNEL;
#define PARAMDOC_APSIDES_PRECISION "absolute precision of the roots of vr^2 for the numerical root finding algorithm."
/**
* @brief absolute precision of the roots of vr^2 for the numerical root finding algorithm.
*/
        int APSIDES_PRECISION;
#define PARAMDOC_APSIDES_MAX_ITER "maximum number of iterations to find the roots of vr^2 numerically"
/**
* @brief maximum number of iterations to find the roots of vr^2 numerically
*/
        int APSIDES_MAX_ITER;
#define PARAMDOC_APSIDES_CONVERGENCE "difference of the roots between two consecutive iterations of the numerical root finding algorithm below which the result is considered to be converged."
/**
* @brief difference of the roots between two consecutive iterations of the numerical root finding algorithm below which the result is considered to be converged.
*/
        int APSIDES_CONVERGENCE;
#define PARAMDOC_CIRC_PERIOD_THRESHOLD "A threshold value for the difference between apastron and periastron below which an orbit is considered to be circular. Currently this is only used for the period calculation in the loss-cone routine."
/**
* @brief A threshold value for the difference between apastron and periastron below which an orbit is considered to be circular. Currently this is only used for the period calculation in the loss-cone routine.
*/
        int CIRC_PERIOD_THRESHOLD;
#define PARAMDOC_WRITE_STELLAR_INFO "Write out information about stellar evolution for each single and binary star, (0=off, 1=on)"
/**
* @brief Write out information about stellar evolution for each single and binary star, (0=off, 1=on)
*/
        int WRITE_STELLAR_INFO;
#define PARAMDOC_WRITE_BH_INFO "Write out information about BHs each timestep, (0=off, 1=on)"
/**
* @brief Write out information about BHs each timestep, (0=off, 1=on)
*/
        int WRITE_BH_INFO;
#define PARAMDOC_WRITE_RWALK_INFO "Write out information about the random walk in J-space around the central black hole, (0=off, 1=on)"
/**
* @brief Write out information about the random walk in J-space around the central black hole, (0=off, 1=on)
*/
        int WRITE_RWALK_INFO;
#define PARAMDOC_WRITE_EXTRA_CORE_INFO "Write out information about cores that are defined differently from the standard (0=off, 1=on)"
/**
* @brief Write out information about cores that are defined differently from the standard (0=off, 1=on)
*/
        int WRITE_EXTRA_CORE_INFO;
#define PARAMDOC_WRITE_PULSAR_INFO "Write out information about pulsars (0=off, 1=on)"
/**
* @brief Write out information about pulsars (0=off, 1=on)
*/
        int WRITE_PULSAR_INFO;
#define PARAMDOC_CALCULATE10 "Write out information about 10\% lagrange radius (0=off, 1=on)"
/**
* @brief Write out information about 10\% lagrange radius (0=off, 1=on)
*/
        int CALCULATE10;
#define PARAMDOC_OVERWRITE_RVIR "Instead of reading the virial radius from the fits file use this value [pc]"
/**
* @brief Instead of reading the virial radius from the fits file use this value [pc]
*/
        int OVERWRITE_RVIR;
#define PARAMDOC_OVERWRITE_Z "Instead of reading the metallicity from the fits file use this value"
/**
* @brief Instead of reading the metallicity from the fits file use this value
*/
        int OVERWRITE_Z;
#define PARAMDOC_OVERWRITE_RTID "Instead of reading the tidal radius from the fits file use this value [pc]"
        int OVERWRITE_RTID;
#define PARAMDOC_OVERWRITE_MCLUS "Instead of reading the cluster mass from the fits file use this value [Msun]"
/**
* @brief Instead of reading the cluster mass from the fits file use this value [Msun]
*/
        int OVERWRITE_MCLUS;
#define PARAMDOC_BSE_NETA "neta > 0 turns wind mass-loss on, is also the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally)."
/**
* @brief neta > 0 turns wind mass-loss on, is also the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally).
*/
	int BSE_NETA;
#define PARAMDOC_BSE_BWIND "bwind is the binary enhanced mass-loss parameter (inactive for single, and normally 0.0 anyway)."
/**
* @brief bwind is the binary enhanced mass-loss parameter (inactive for single, and normally 0.0 anyway).
*/
	int BSE_BWIND;
#define PARAMDOC_BSE_HEWIND "hewind is a helium star mass-loss factor (0.5 normally)."
/**
* @brief hewind is a helium star mass-loss factor (0.5 normally).
*/
	int BSE_HEWIND;
#define PARAMDOC_BSE_WINDFLAG "windflag sets which wind prescription to use (0=BSE, 1=StarTrack, 2=Vink, 3=Vink+LBV for all stars)."
/**
* @brief windflag sets which wind prescription to use (0=BSE, 1=StarTrack, 2=Vink, 3=Vink+LBV for all stars).
*/
	int BSE_WINDFLAG;
#define PARAMDOC_BSE_ALPHA1 "alpha1 is the common-envelope efficiency parameter (1.0 or 3.0 depending upon what you like and if lambda is variable)"
/**
* @brief alpha1 is the common-envelope efficiency parameter (1.0 or 3.0 depending upon what you like and if lambda is variable)
*/
	int BSE_ALPHA1;
#define PARAMDOC_BSE_LAMBDA "labmda is the stellar binding energy factor for common-envelope evolution (0.5; +'ve allows it to vary, -'ve holds it constant at that value always)."
/**
* @brief labmda is the stellar binding energy factor for common-envelope evolution (0.5; +'ve allows it to vary, -'ve holds it constant at that value always).
*/
	int BSE_LAMBDA;
#define PARAMDOC_BSE_CEFLAG "ceflag sets CE prescription used. = 0 sets Tout et al. method, = 3 activates de Kool common-envelope models (normally 0)."
/**
* @brief ceflag sets CE prescription used. = 0 sets Tout et al. method, = 3 activates de Kool common-envelope models (normally 0).
*/
	int BSE_CEFLAG;
#define PARAMDOC_BSE_TFLAG "tflag > 0 activates tidal circularisation (1)."
/**
* @brief tflag > 0 activates tidal circularisation (1).
*/
	int BSE_TFLAG;
#define PARAMDOC_BSE_IFFLAG "ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0)."
/**
* @brief ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0).
*/
	int BSE_IFFLAG;
#define PARAMDOC_BSE_WDFLAG "wdflag > 0 uses modified-Mestel cooling for WDs (1)."
/**
* @brief wdflag > 0 uses modified-Mestel cooling for WDs (1).
*/
	int BSE_WDFLAG;
#define PARAMDOC_BSE_BHFLAG "bhflag > 0 allows velocity kick at BH formation (1)."
/**
* @brief bhflag > 0 allows velocity kick at BH formation (1).
*/
	int BSE_BHFLAG;
#define PARAMDOC_BSE_NSFLAG "0 is default BSE, 1 is Belczynski 2002 Model, 2 is Belczynski 2008, 3 is Fryer 2012 'Rapid' SN, 4 is Fryer 2012 'Delayed' SN"
/**
* @brief nsflag = 0 gives the NS/BH mass distribution from default BSE; 1 uses the distribution of Belczynski et al. 2002, ApJ, 572, 407 (1), while 2, 3, and 4 give you the standard, 'rapid' and 'delayed' models of Fryer et al. 2012, ApJ, 749, 91
*/
	int BSE_NSFLAG;
#define PARAMDOC_BSE_MXNS "mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1)."
/**
* @brief mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1).
*/
	int BSE_MXNS;
#define PARAMDOC_BSE_BCONST "bconst is the magnetic field decay timescale (-3000, although value and decay rate not really established...)."
/**
* @brief bconst is the magnetic field decay timescale (-3000, although value and decay rate not really established...).
*/
	int BSE_BCONST;
#define PARAMDOC_BSE_CK "CK is an accretion induced field decay constant (-1000, although again this isn't well established...)."
/**
* @brief CK is an accretion induced field decay constant (-1000, although again this isn't well established...).
*/
	int BSE_CK;
#define PARAMDOC_BSE_IDUM "idum in the random number seed used by kick.f and setting initial pulsar spin period and magnetic field."
/**
* @brief idum in the random number seed used by kick.f and setting initial pulsar spin period and magnetic field.
*/
	int BSE_IDUM;
#define PARAMDOC_BSE_SIGMA "sigma is the Maxwellian dispersion for SN kick speeds (265 km/s)."
/**
* @brief sigma is the Maxwellian dispersion for SN kick speeds (265 km/s).
*/
	int BSE_SIGMA;
#define PARAMDOC_BSE_BHSIGMAFRAC "bhsigmafrac is the factor to scale Maxwellian dispersion for BH SN kick speeds (1) compared to NS (265 km/s)."
/**
* @brief bhsigmafrac is the factor to scale Maxwellian dispersion for BH SN kick speeds (1) compared to NS (265 km/s).
*/
	int BSE_BHSIGMAFRAC;
#define PARAMDOC_BSE_OPENING_ANGLE "Opening angle (in degrees) for SN kicks around the north pole of the star.  Default (180) means isotropic kicks"
/**
* @brief Opening angle for SN kicks about the north pole of the star, in
* radians.  Default (180) means isotropic over whole sphere 
*/
	int BSE_OPENING_ANGLE;
#define PARAMDOC_BSE_BETA "beta is the wind velocity factor: proprotinal to vwind^2 (1/8)."
/**
* @brief beta is the wind velocity factor: proprotinal to vwind^2 (1/8).
*/
	int BSE_BETA;
#define PARAMDOC_BSE_EDDFAC "eddfac is Eddington limit factor for mass transfer (1.0, 10 turns Eddlimitation off)."
/**
* @brief eddfac is Eddington limit factor for mass transfer (1.0, 10 turns Eddlimitation off).
*/
	int BSE_EDDFAC;
#define PARAMDOC_BSE_GAMMA "gamma is the angular momentum factor for mass lost during Roche (-1.0, see evolv2.f for more details)."
/**
* @brief gamma is the angular momentum factor for mass lost during Roche (-1.0, see evolv2.f for more details).
*/
	int BSE_GAMMA;
#define PARAMDOC_TIMER "enable or disable timers. This would return a detailed profiling of the code, but uses barriers, so might slow down the code a bit."
/**
* @brief enable or disable timers. This would return a detailed profiling of the code, but uses barriers, so might slow down the code a bit.
*/
	int TIMER;
} parsed_t;


/**
* @brief  struct containing the units used
*/
typedef struct{
/**
* @brief  time
*/
	double t;
/**
* @brief  mass
*/
	double m;
/**
* @brief  length
*/
	double l;
/**
* @brief  energy
*/
	double E;
/**
* @brief  stars' masses are kept in different units
*/
	double mstar;
} units_t;

/**
* @brief  a struct for the potential, must be malloc'ed
*/
typedef struct{
/**
* @brief  number of elements
*/
	long n;
/**
* @brief  radius
*/
	double *r;
/**
* @brief  potential
*/
	double *phi;
} potential_t;

/**
* @brief a struct for the force, must be malloc'ed
*/
typedef struct{
/**
* @brief  number of elements
*/
	long n;
/**
* @brief  radius
*/
	double *r;
/**
* @brief  force
*/
	double *force;
} force_t;

/**
* @brief useful structure for central quantities
*/
typedef struct{
/**
* @brief  central mass density
*/
	double rho;
/**
* @brief  rms object velocity
*/
	double v_rms;
/**
* @brief  core radius
*/
	double rc;
/**
* @brief  average object mass
*/
	double m_ave;
/**
* @brief  number density of objects
*/
	double n;
/**
* @brief  Spitzer definition of core radius
*/
	double rc_spitzer;
/**
* @brief  number of objects that are single
*/
	long N_sin;
/**
* @brief  number of objects that are binary
*/
	long N_bin;
/**
* @brief  single star number density
*/
	double n_sin;
/**
* @brief  binary star number density
*/
	double n_bin;
/**
* @brief  central single star mass density
*/
	double rho_sin;
/**
* @brief  central binary star mass density
*/
	double rho_bin;
/**
* @brief  average single star mass
*/
	double m_sin_ave;
/**
* @brief  average binary star mass
*/
	double m_bin_ave;
/**
* @brief  rms single star velocity
*/
	double v_sin_rms;
/**
* @brief  rms binary star velocity
*/
	double v_bin_rms;
/**
* @brief average of 2*m*v^2 per average mass for all objects
*/
	double w2_ave;
/**
* @brief  average of R^2 for single stars
*/
	double R2_ave;
/**
* @brief  average of m*R for single stars
*/
	double mR_ave;
/**
* @brief  average of a for binaries
*/
	double a_ave;
/**
* @brief  average of a^2 for binaries
*/
	double a2_ave;
/**
* @brief  average of m*a for binaries
*/
	double ma_ave;
} central_t;

/* useful structure for core quantities */
/*typedef struct{
	double N;
	double kT;
} core_t;*/

/**
* @brief  to store the velocity dispersion profile
*/
typedef struct{
/**
* @brief ?
*/
	long n;
/**
* @brief ?
*/
	double *r;
/**
* @brief ?
*/
	double *sigma;
} sigma_t;

/**
* @brief parameters for orbit
*/
typedef struct{
/**
* @brief ?
*/
	double rp;
/**
* @brief ?
*/
	double ra;
/**
* @brief ?
*/
	double dQdrp;
/**
* @brief ?
*/
	double dQdra;
/**
* @brief ?
*/
	long kmax;
/**
* @brief ?
*/
	long kmin;
/**
* @brief ?
*/
	int circular_flag;
} orbit_rs_t;

/**
* @brief  parameters for calc_p_orb function
*/
typedef struct{
/**
* @brief ?
*/
	double E;
/**
* @brief ?
*/
	double J;
/**
* @brief ?
*/
	long index;
/**
* @brief ?
*/
	long kmin;
/**
* @brief ?
*/
	long kmax;
/**
* @brief ?
*/
	double rp;
/**
* @brief ?
*/
	double ra;
} calc_p_orb_params_t;

/* other useful structs */
/**
* @brief ?
*/
typedef struct{
/**
* @brief ?
*/
	long N_MAX;
/**
* @brief ?
*/
	long N_MAX_NEW;
/**
* @brief ?
*/
	long N_STAR;
/**
* @brief ?
*/
	long N_STAR_NEW;
/**
* @brief ?
*/
	long N_BINARY;
} clus_struct_t;

/**
* @brief ?
*/
typedef struct{
/**
* @brief  total = kinetic + potential + internal + binary binding + central mass energy
*/
	double tot;
/**
* @brief  ???
*/
	double New;
/**
* @brief  initial total energy
*/
	double ini;
/**
* @brief  total kinetic
*/
	double K;
/**
* @brief  total potential
*/
	double P;
/**
* @brief  total internal
*/
	double Eint;
/**
* @brief  total binary binding energy
*/
	double Eb;
} Etotal_struct_t;

/**
* @brief ?
*/
typedef struct{
/**
* @brief  half-mass radius
*/
	double rh;
} clusdyn_struct_t;

/********************** Function Declarations ************************/
double sqr(double x);
double cub(double x);
void tidally_strip_stars(void);

void remove_star_center(long j);
void print_results(void);
void print_conversion_script(void);
double potential(double r);	       /* get potential using star.phi */
double potential_serial(double r);
double fastpotential(double r, long kmin, long kmax);
long potential_calculate(void);	/* calculate potential at star locations in star.phi */
long potential_calculate_mimic(void); //to mimic the parallel version
#ifdef USE_MPI
long mpi_potential_calculate(void);
long mpi_potential_calculate2(void);
MPI_Comm inv_comm_create();
#endif

/* Bharath: Timing Functions */ 
double timeStartSimple();
void timeEndSimple(double timeStart, double *timeAccum);
void timeStart();
void timeEnd(char* fileName, char *funcName, double *tTime);
void timeStart2(double *st);
void timeEnd2(char* fileName, char *funcName, double *st, double *end, double *tot);
void create_timing_files();
/* End */

/* Bharath: Other refactored functions */
void mpiBcastGlobArrays();
void mpiInitGlobArrays();
void set_global_vars1();
void set_global_vars2();
void get_star_data(int argc, char *argv[], gsl_rng *rng);
void calc_sigma_new();
void calc_central_new();
void bin_vars_calculate();
void calc_potential_new();
void calc_potential_new2();
void compute_energy_new();
void set_energy_vars();
void reset_interaction_flags();
void calc_clusdyn_new();
void calc_timestep(gsl_rng *rng);
void energy_conservation1();
void energy_conservation2();
void new_orbits_calculate();
void toy_rejuvenation();
void pre_sort_comm();
void post_sort_comm();
void findIndices( long N, int blkSize, int i, int* begin, int* end );
void pulsar_write(long k, double kick);
void findLimits( long N, int blkSize );
int findProcForIndex( int j );
void set_rng_states();
int get_global_idx(int i);
int get_local_idx(int i);
/* End */

/* Bharath: Functions for handling of binaries for parallel version */
void distr_bin_data();
void collect_bin_data();
void alloc_bin_buf();
/* End */

/* Function to handle I/O in parallel version */
//void cat_and_rm_files(char* file_ext)
//void mpi_merge_files();
//void save_root_files();
//void save_root_files_helper(char* file_ext);
//void rm_files();
//void rm_files_helper(char* file_ext);
//void print_small_output();
void print_denprof_snapshot(char* infile);
/* End */

void comp_mass_percent(void);
void comp_multi_mass_percent(void);
orbit_rs_t calc_orbit_rs(long si, double E, double J);
double get_positions(void);	/* get positions and velocities */
void perturb_stars(double Dt);	/* take a time step (perturb E,J) */
long FindZero_r_serial(long kmin, long kmax, double r);
long FindZero_r(long x1, long x2, double r);
long FindZero_Q(long j, long x1, long x2, double E, double J);
double potentialDifference(int particleIndex);
void ComputeEnergy(void);

#ifdef USE_MPI
void mpi_ComputeEnergy(void);
void mpi_para_file_write(char* wrbuf, long long *len, long long *prev_cum_offset, MPI_File* fh);
void PrintParaFileOutput(void);
void mpi_close_node_buffers(void);
#endif

void PrintLogOutput(void);
double GetTimeStep(gsl_rng *rng);
long CheckStop(struct tms tmsbuf);
void ComputeIntermediateEnergy(void);
void energy_conservation3(void);
void set_velocities3(void);
void mpi_set_velocities3(void);
double rtbis(double (*func) (double), double x1, double x2, double xacc);
double trapzd(double (*func) (double), double a, double b, long n);
double qsimp(double (*func) (double), double a, double b);

void splint(double xa[], double ya[], double y2a[], long n, double x, double *y);
void spline(double x[], double y[], long n, double yp1, double ypn, double y2[]);
void print_2Dsnapshot(void);
void print_bh_snapshot(void);
void get_physical_units(void);
void update_vars(void);

void print_version(FILE *stream);
void cmc_print_usage(FILE *stream, char *argv[]);
void parser(int argc, char *argv[], gsl_rng *r);
void PrintFileOutput(void);
void find_nstars_within_r(double r, long *ns, long *nb);
char *sprint_star_dyn(long k, char string[MAX_STRING_LENGTH]);
char *sprint_bin_dyn(long k, char string[MAX_STRING_LENGTH]);

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
void do_stellar_evolution(gsl_rng *rng);
void write_stellar_data(void);
void handle_bse_outcome(long k, long kb, double *vs, double tphysf, int kprev0, int kprev1);
void cp_binmemb_to_star(long k, int kbi, long knew);
void cp_SEvars_to_newstar(long oldk, int kbi, long knew);
void cp_m_to_newstar(long oldk, int kbi, long knew);
void cp_SEvars_to_star(long oldk, int kbi, star_t *target_star);
void cp_m_to_star(long oldk, int kbi, star_t *target_star);
void cp_SEvars_to_newbinary(long oldk, int oldkbi, long knew, int kbinew);
void cp_starSEvars_to_binmember(star_t instar, long binindex, int bid);
void cp_starmass_to_binmember(star_t instar, long binindex, int bid);
double r_of_m(double M);
void cmc_bse_comenv(binary_t *tempbinary, double cmc_l_unit, double RbloodySUN, double *zpars, double *vs, int *fb, double *ecsnp, double *ecsn_mlow, int *ST_tide);

/* Fewbody stuff */
void destroy_obj(long i);
void destroy_obj_new(long i);
void destroy_binary(long i);
long create_star(int idx, int dyn_0_se_1);
long create_binary(int idx, int dyn_0_se_1);
void dynamics_apply(double dt, gsl_rng *rng);
void perturb_stars_fewbody(double dt, gsl_rng *rng);
void qsorts_new(void);
void qsorts(star_t *s, long n);
void units_set(void);
void central_calculate(void);

#ifdef USE_MPI
typedef star_t type;
typedef double keyType;
void remove_stripped_stars(type* buf, int* local_N);
int sample_sort( 	type			*buf,
						int			*local_N,
						MPI_Datatype dataType,
						binary_t		*b_buf,
						MPI_Datatype bDataType,
						MPI_Comm    commgroup,
						int			n_samples );
void load_balance( 	type 				*inbuf,
							type 				*outbuf,
							binary_t			*b_inbuf,
							binary_t			*b_outbuf,
							int 				*expected_count, 
							int				*actual_count, 
							int 				myid, 
							int 				procs, 
							MPI_Datatype 	dataType, 	
							MPI_Datatype 	b_dataType, 	
							MPI_Comm			commgroup	);
#endif

central_t central_hard_binary(double ktmin, central_t old_cent);
void clusdyn_calculate(void);
#ifdef USE_MPI
void mpi_clusdyn_calculate(void);
void copy_globals_to_locals(long k);
#endif

void print_interaction_status(char status_text[]);
void print_interaction_error(void);
long star_get_id_new(void);
long star_get_merger_id_new(long id1, long id2);
double calc_n_local(long k, long p, long N_LIMIT);
double calc_Ai_local(long k, long kp, long p, double W, long N_LIMIT);
void calc_encounter_dyns(long k, long kp, double v[4], double vp[4], double w[4], double *W, double *rcm, double vcm[4], gsl_rng *rng, int setY);
void set_star_EJ(long k);
void set_star_news(long k);
void set_star_olds(long k);
void zero_star(long j);
void zero_binary(long j);

void sscollision_do(long k, long kp, double rperi, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng);
void merge_two_stars(star_t *star1, star_t *star2, star_t *merged_star, double *vs, struct rng_t113_state* s);
double coll_CE(double Mrg, double Mint, double Mwd, double Rrg, double vinf);

void print_initial_binaries(void);

void bs_calcunits(fb_obj_t *obj[2], fb_units_t *bs_units);
fb_ret_t binsingle(double *t, long ksin, long kbin, double W, double bmax, fb_hier_t *hier, gsl_rng *rng);

void bb_calcunits(fb_obj_t *obj[2], fb_units_t *bb_units);
fb_ret_t binbin(double *t, long k, long kp, double W, double bmax, fb_hier_t *hier, gsl_rng *rng);

double binint_get_mass(long k, long kp, long id);
long binint_get_startype(long k, long kp, long id);
long binint_get_indices(long k, long kp, long id, int *bi);
void binint_log_obj(fb_obj_t *obj, fb_units_t units);
void binint_log_status(fb_ret_t retval);
void binint_log_collision(const char interaction_type[], long id, double mass, double r, fb_obj_t obj, long k, long kp, long startype);
void binint_do(long k, long kp, double rperi, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng);

double simul_relax(gsl_rng *rng);
double simul_relax_new(void);
#ifdef USE_MPI
double mpi_simul_relax_new(void);
void mpi_calc_sigma_r(long p, long N_LIMIT, double *sig_r, double *sig_sigma, long* sig_n, int r_0_mave_1);
void mpi_break_wide_binaries(void);
#endif

void break_wide_binaries(void);
void calc_sigma_r(long p, long N_LIMIT, double *sig_r, double *sig_sigma, long* sig_n, int r_0_mave_1);

double sigma_r(double r);

// Meagan
/* three-body binary formation */
void sort_three_masses(long sq, long *k1, long *k2, long *k3);
#ifdef USE_MPI
void mpi_sort_three_masses(long sq, long *k1, long *k2, long *k3);
#endif
double get_eta(double eta_min, long k1, long k2, long k3, double vrel12[4], double vrel3[4]);
void calc_3bb_encounter_dyns(long k1, long k2, long k3, double v1[4], double v2[4], double v3[4], double (*vrel12)[4], double (*vrel3)[4], gsl_rng *rng);
void make_threebodybinary(double P_3bb, long k1, long k2, long k3, long form_binary, double eta_min, double ave_local_mass, double n_local, double sigma_local, double v1[4], double v2[4], double v3[4], double vrel12[4], double vrel3[4], double delta_E_running, gsl_rng *rng);
void calc_sigma_local(long k1, long p, long N_LIMIT, double *ave_local_mass, double *sigma_local);
int remove_old_star(double time, long k);

// Meagan
/* extra output for bhs */
void bh_count(long k);
void print_bh_summary(void);
void count_esc_bhs(long j);
void print_esc_bh_summary(void);

/* signal/GSL error handling stuff */
void toggle_debugging(int signal);
void exit_cleanly_old(int signal);
void exit_cleanly(int signal, const char* fn);
void sf_gsl_errhandler(const char *reason, const char *file, int line, int gsl_errno);
void close_buffers(void);
void close_root_buffers(void); //files that are opened only by the root node in MPI version.
void close_node_buffers(void); //files that need to be opened by all nodes using MPI-IO.
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

/**
* @brief ?
*/
struct Search_Grid {
/* user modifiable parameters*/
/**
* @brief
* starsPerBin is the estimated (anticipated) number of stars per search
* bin. It is only used to calculate interpol_coeff and grid length.
*/
   long starsPerBin;
/**
* @brief ?
*/
   long max_length;
/**
* @brief ?
*/
   long min_length;
/**
* @brief
* fraction determines the interpol_coeff such that
* grid.radius[fraction*grid->length] is approx.
* star[clus.N_MAX*fraction].r. Depending how well the power_law_exponent
* fits the distribution of star[i::starsPerBin].r the better this is
* fulfilled
*/
   double fraction;
/**
* @brief
* Only particle_fraction * clus.N_MAX stars should be used to construct the
* grid (i.e. we ignore a possible outer halo).
*/
   double particle_fraction;
/**
* @brief ?
*/
   double power_law_exponent;
/**
* @brief
* You should not change any of the following variables yourself
* unless, of course, you know what you are doing.
*/
   long *radius;
/**
* @brief ?
*/
   long length;
/**
* @brief ?
*/
   double interpol_coeff;
};

/**
* @brief
*/
struct Interval {
/**
* @brief ?
*/
  long min;
/**
* @brief ?
*/
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
void search_grid_print_binsizes(struct Search_Grid *grid);

#ifdef DEBUGGING
#include <glib.h>
void load_id_table(GHashTable* ids, char *filename);
#endif
struct star_coords get_position(long index, double E, double J, double old_r, orbit_rs_t orbit);
orbit_rs_t calc_orbit_new_J(long index, double J, struct star_coords old_pos, orbit_rs_t orbit_old);
void set_a_b(long index, long k, double *a, double *b);
long find_zero_Q_slope(long index, long k, double E, double J, int positive);
//long get_positive_Q_index(long index, double E, double J);
long find_zero_Q(long j, long kmin, long kmax, long double E, long double J);
//extern inline long double function_q(long j, long double r, long double pot, long double E, long double J);
orbit_rs_t calc_orbit_new(long index, double E, double J);
double calc_average_mass_sqr(long index, long N_LIMIT);
double floateq(double a, double b);
double sigma_tc_nd(double n, double m1, double r1, double m2, double vinf);
double sigma_tc_nn(double na, double ma, double ra, double nb, double mb, double rb, double vinf);
double Tl(int order, double polytropicindex, double eta);
double Etide(double rperi, double Mosc, double Rosc, double nosc, double Mpert);
struct Interval get_r_interval(double r);
double calc_vr_within_interval(double r, void *p);
double calc_vr_in_interval(double r, long index, long k, double E, double J);
double calc_vr(double r, long index, double E, double J);
double find_root_vr(long index, long k, double E, double J);
double calc_pot_in_interval(double r, long k);
double local_kT(long si, int p);
void remove_star(long j, double phi_rtidal, double phi_zero);
double function_q(long j, long double r, long double pot, long double E, long double J);
void vt_add_kick(double *vt, double vs1, double vs2, struct rng_t113_state* rng_st);

void parse_snapshot_windows(char *option_string);
void print_snapshot_windows(void);
int valid_snapshot_window_units(void);
void write_snapshot(char *filename, int bh_only);

#include "cmc_bse_utils.h"

/* macros */ 
#define Q_function(j, r, pot, E, J) (2.0 * ((E) - ((pot) + PHI_S((r), j))) - SQR((J) / (r)) )
/* correction to potential due to subtracting star's contribution, and adding self-gravity */
/* #define PHI_S(rad, j) ( ((rad)>=star[(j)].r ? star[(j)].m/(rad) : star[(j)].m/star[(j)].r) * (1.0/clus.N_STAR) - 0.5*star[(j)].m/clus.N_STAR/(rad) ) */
/* correction to potential due to subtracting star's contribution (this is what Marc uses) */
#define PHI_S(rad, j) ( ((rad)>=star[(j)].r ? star[(j)].m/(rad) : star[(j)].m/star[(j)].r) * (1.0/clus.N_STAR) )
#define MPI_PHI_S(rad, j) ( ((rad)>=star_r[(j)] ? star_m[(j)]/(rad) : star_m[(j)]/star_r[(j)]) * (1.0/clus.N_STAR) )

#ifdef USE_MPI
#define function_Q(j, k, E, J) (2.0 * ((E) - (star_phi[(k)] + MPI_PHI_S(star_r[(k)], j))) - SQR((J) / star_r[(k)]))
#else
#define function_Q(j, k, E, J) (2.0 * ((E) - (star[(k)].phi + PHI_S(star[(k)].r, j))) - SQR((J) / star[(k)].r))
#endif

#define gprintf(args...) if (!quiet) {fprintf(stdout, args);}

#ifdef USE_MPI
#define diaprintf(args...) if (!quiet) { if(myid == 0) { fprintf(stdout, "DIAGNOSTIC: %s(): ", __FUNCTION__); fprintf(stdout, args); }}
#else
#define diaprintf(args...) if (!quiet) {fprintf(stdout, "DIAGNOSTIC: %s(): ", __FUNCTION__); fprintf(stdout, args);}
#endif

#ifdef USE_MPI
#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: in proc %d, %s(): ", myid, __FUNCTION__); fprintf(stderr, args);}
#define rootdprintf(args...) if (debug && myid==0) {fprintf(stderr, "DEBUG: in proc %d, %s(): ", myid, __FUNCTION__); fprintf(stderr, args);}
#else
#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#define rootdprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#endif

#define dmpiprintf(args...) if (mpi_debug) { fprintf(stderr, "DEBUG: in proc %d, %s(): ", myid, __FUNCTION__); fprintf(stderr, args); }

#ifdef USE_MPI
#define wprintf(args...) { if(myid==0) { fprintf(stderr, "WARNING: %s(): ", __FUNCTION__); fprintf(stderr, args);}}
#else
#define wprintf(args...) { fprintf(stderr, "WARNING: %s(): ", __FUNCTION__); fprintf(stderr, args);}
#endif

#ifdef USE_MPI
#define eprintf(args...) {fprintf(stderr, "ERROR: in proc %d: %s:%d in %s(): ", myid, __FILE__, __LINE__, __FUNCTION__); fprintf(stderr, args);}
#else
#define eprintf(args...) {fprintf(stderr, "ERROR: %s:%d in %s(): ", __FILE__, __LINE__, __FUNCTION__); fprintf(stderr, args);}
#endif

//MPI3-IO: This was the easiest way to convert hundreds of fprintf statements to do parallel IO without manually changing each one of them.
#ifdef USE_MPI

/**
* @brief Macro that prints given args into char buffer of the corresponding file.
*
* @param file File to be written to. This macro does not actually write to the file, but instead stores the data into the corresponding char buffer.
* @param args... arguments for standard printf
*/
#define parafprintf(file, args...)                       \
do { sprintf(mpi_ ## file ## _buf, args);               \
strcat(mpi_ ## file ## _wrbuf, mpi_ ## file ## _buf);   \
mpi_ ## file ## _len += strlen(mpi_ ## file ## _buf); }         \
while(0)

/**
* @brief Prints out given arguments into char buffer corresponding to given file, done only by the root node.
*
* @param file File to be written to
* @param args... arguments for standard printf
*/
#define pararootfprintf(file, args...) { if(myid==0) parafprintf(file, args); }

/**
* @brief Prints given arguments to given file only by the root node. Translates to regular printf in the serial version.
*
* @param file File to be written to
* @param args... arguments for standard printf
*/
#define rootfprintf(file, args...) { if(myid==0) fprintf(file, args); }

/**
* @brief Calls gprintf only for root node with given arguments.
*
* @param file File to be written to
* @param args... arguments for standard printf
*/
#define rootgprintf(args...) if (!quiet) {if(myid==0) fprintf(stdout, args);}

/**
* @brief Prints given arguments to stdout only by the root node. Translates to regular printf in the serial version.
*
* @param args... arguments for standard printf
*/
#define rootprintf(args...) {if(myid == 0) { fprintf(stdout, args); }}

#else

#define parafprintf(file, args...) { fprintf(file, args); }

#define pararootfprintf(file, args...) { parafprintf(file, args); }

#define rootfprintf(file, args...) { fprintf(file, args); }

#define rootgprintf(args...) if (!quiet) {fprintf(stdout, args);}

#define rootprintf(args...) { fprintf(stdout, args); }

#endif


#ifdef DEBUGGING
#undef MAX
#undef MIN
#endif
#define MAX(a, b) ((a)>=(b)?(a):(b))
#define MIN(a, b) ((a)<=(b)?(a):(b))
