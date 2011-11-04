typedef struct {
	double se_mass;
	int se_k;
	double se_mt;
	double se_ospin;
        double se_B_0; /* Pulsar initial magentif field */
        double se_bacc;
        double se_tacc;
	double se_epoch;
	double se_tphys;
	double se_radius;
	double se_lum;
	double se_mc;
	double se_rc;
	double se_menv;
	double se_renv;
	double se_tms;
        double se_scm_B; /* Pulsar surface magnetic field */
        double se_scm_formation; /* formation pathway of NS */
} sse_t;

typedef struct {
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
	double bse_tphys; /* physical time */
	double bse_tb; /* binary orbital period */
	double bse_bcm_dmdt[2]; /* mass transfer rate for each star [bse_get_bcm(i,14), bse_get_bcm(i,28)] */
	double bse_bcm_radrol[2]; /* radius/roche_lobe_radius for each star [bse_get_bcm(i,15), bse_get_bcm(i,29)] */
        double bse_bcm_B[2]; /* Pulsar magnetic field strength at surface */
        double bse_bcm_formation[2]; /* provides formation pathway of NS */
} bse_t;

void update_bse_from_sse(bse_t *bvars, sse_t *svars, int bmember);
void update_sse_from_bse(bse_t *bvars, sse_t *svars, int bmember);
void get_sse_from_star(sse_t *svars, star_t *star);
void update_star_from_sse(star_t *star, sse_t svars);
void get_bse_from_binary(bse_t *bvars, binary_t *star);
void update_binary_from_bse(binary_t *star, bse_t *bvars);
void compress_binary(star_t *bincom, binary_t *bin);

