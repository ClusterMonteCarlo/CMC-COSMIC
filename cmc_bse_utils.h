/* vi: set filetype=c.doxygen: */

/**
* @brief sse structure?
*/
typedef struct {
/**
* @brief ?
*/
	double se_mass;
/**
* @brief ?
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
* @brief Pulsar initial magentif field
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
* @brief Pulsar surface magnetic field
*/
    double se_scm_B;
/**
* @brief formation pathway of NS
*/
    double se_scm_formation;
} sse_t;



/**
* @brief bse structure?
*/
typedef struct {
/**
* @brief star types
*/
	int bse_kw[2];
/**
* @brief initial masses
*/
	double bse_mass0[2];
/**
* @brief masses
*/
	double bse_mass[2];
/**
* @brief radii
*/
	double bse_radius[2];
/**
* @brief luminosity
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
* @brief original spin
*/
	double bse_ospin[2];
/**
* @brief Pulsar magnetic field
*/
    double bse_B_0[2];
/**
* @brief Amount of mass pulsar has accreted
*/
    double bse_bacc[2];
/**
* @brief Amount of time pulsar has spent accreting
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
* @brief physical time
*/
	double bse_tphys;
/**
* @brief binary orbital period
*/
	double bse_tb;
/**
* @brief mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)]
*/
	double bse_bcm_dmdt[2];
/**
* @brief radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)]
*/
	double bse_bcm_radrol[2];
/**
* @brief Pulsar magnetic field strength at surface
*/
    double bse_bcm_B[2];
/**
* @brief provides formation pathway of NS
*/
    double bse_bcm_formation[2];
} bse_t;

void update_bse_from_sse(bse_t *bvars, sse_t *svars, int bmember);
void update_sse_from_bse(bse_t *bvars, sse_t *svars, int bmember);
void get_sse_from_star(sse_t *svars, star_t *star);
void update_star_from_sse(star_t *star, sse_t svars);
void get_bse_from_binary(bse_t *bvars, binary_t *star);
void update_binary_from_bse(binary_t *star, bse_t *bvars);
void compress_binary(star_t *bincom, binary_t *bin);

