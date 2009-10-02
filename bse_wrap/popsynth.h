#define NBIN 500000
/* NBIN must be divisible by BATCHSIZE */
#define BATCHSIZE 500

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
	double bse_epoch[2];
	double bse_tms[2];
	double bse_tphys; /* physical time */
	double bse_tb; /* binary orbital period */
	double bse_bcm_dmdt[2]; /* mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)] */
	double bse_bcm_radrol[2]; /* radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)] */
	double m1init;
	double m2init;
	double ainit;
	double einit;
} binary_t;

#define RSUN 6.9599e10
#define AU 1.496e+13
