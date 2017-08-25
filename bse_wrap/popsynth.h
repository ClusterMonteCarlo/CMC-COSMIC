#define NBIN 50
/* NBIN must be divisible by BATCHSIZE */
#define BATCHSIZE 5

typedef struct{
	long id1; /* unique id of star 1 */
	long id2; /* unique id of star 2 */
	double a; /* semimajor axis */
	double e; /* eccentricity */
	int bse_kw[2]; /* star types */
	double bse_mass0[2]; /* initial masses */
	double bse_mass[2]; /* masses */
	double bse_radius[2]; /* radii */
	double bse_lum[2]; /* luminosity */
	double bse_massc[2];
	double bse_radc[2];
	double bse_menv[2];
	double bse_renv[2];
	double bse_ospin[2]; /* original spin */
        double bse_B_0[2]; /* Pulsar magnetic field */
        double bse_bacc[2]; /* Amount of mass pulsar has accreted */
        double bse_tacc[2]; /* Amount of time pulsar has spent accreting */
	double bse_epoch[2];
	double bse_tms[2];
	double bse_bhspin[2];
	double bse_tphys; /* physical time */
	double bse_tb; /* binary orbital period */
	double bse_bcm_dmdt[2]; /* mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)] */
	double bse_bcm_radrol[2]; /* radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)] */
        double bse_bcm_B[2]; /* PK: Pulsar surface magnetic field strength */
/* parameter providing formation scenario 4-TypeII, 5-ECSN(when implemented), 6-AIC, 7-MIC*/
        double bse_bcm_formation[2];
	double m1init;
	double m2init;
	double ainit;
	double einit;
} binary_t;

#define RSUN 6.9599e10
#define AU 1.496e+13
