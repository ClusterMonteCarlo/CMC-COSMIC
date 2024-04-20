#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "cmc.h"
#include "fitslib.h"
#include "taus113-v2.h"
#include "bse_wrap.h"

#define INFILE "in.fits"
#define OUTFILE "debug.fits"

/* global variables needed by Star Track */
double *zpars, METALLICITY, WIND_FACTOR=1.0;

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  setstellar [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -i --infile <infile>   : input file [%s]\n", INFILE);
	fprintf(stream, "  -o --outfile <outfile> : output file [%s]\n", OUTFILE);
	fprintf(stream, "  -h --help              : display this help text\n");
}

/* set stellar type and radius */
void stellar_evolve(cmc_fits_data_t *cfd)
{
	long int k;
	double frac, dtp, tphysf;
	star_t star;
	double vs[20];

	/* initialize stellar evolution */
	/* SSE */
	/* bse_set_hewind(0.5); */
	
	/* BSE */
	/* These are the defaults from cmc_io.c
	 * Since we really only need the radii in this step, it shouldn't matter 
	 * that much, but best to be consistant*/
	double BSE_FPRIMC_ARRAY[16] = {0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095,0.095238095};
	double BSE_QCRIT_ARRAY[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double BSE_NATAL_KICK_ARRAY[10] = {-100.0,-100.0,-100.0,-100.0,0,-100.0,-100.0,-100.0,-100.0,0.0};
        int NO_BSE_NATAL_KICK_ARRAY = 10;
	int NO_BSE_FPRIMC_ARRAY = 16;
	int NO_BSE_QCRIT_ARRAY = 16;
        double BSE_NETA= 0.5;
        double BSE_BWIND= 0.0;
        double BSE_HEWIND= 1.0;
        double BSE_REMBAR_MASSLOSS = 0.5;
        int BSE_KICKFLAG = 0;
        int BSE_GRFLAG = 0;
        double BSE_ZSUN = 0.017;
        double BSE_BETA= -1.0;
        double BSE_XI= 0.5;
        double BSE_ACC2= 1.5;
        double BSE_ALPHA1= 1.0;
        double BSE_LAMBDAF= 0.5;
        int BSE_CEFLAG= 0;
        int BSE_CEKICKFLAG= 2;
        int BSE_CEMERGEFLAG= 0;
        int BSE_CEHESTARFLAG= 0;
        int BSE_QCFLAG= 2;
        double BSE_DON_LIM= 0;
        double BSE_ACC_LIM= 0;
        double BSE_SIGMA= 265.0;
        int BSE_BHFLAG= 1;
        int BSE_BHMS_COLL_FLAG= 0;
        double BSE_ECSN= 2.5;
        double BSE_ECSN_MLOW= 1.4;
        double BSE_SIGMADIV= -20.0;
        double BSE_PTS1= 0.05;
        double BSE_PTS2= 0.01;
        double BSE_PTS3= 0.02;
        int BSE_AIC= 1;
        int BSE_BDECAYFAC= 1;
        int BSE_HTPMB= 1;
        int BSE_ST_CR= 1;
        int BSE_ST_TIDE= 0;
        int BSE_REJUVFLAG= 0;
        int BSE_USSN= 0;
        double BSE_PISN= 45.00;
        double BSE_BHSIGMAFRAC= 1.00;
        double BSE_POLAR_KICK_ANGLE= 90.00;
        int BSE_REMNANTFLAG= 4;
        double BSE_MXNS= 3.00;
        int BSE_WD_MASS_LIM= 1;
        int BSE_BHSPINFLAG= 0;
        double BSE_WINDFLAG= 3;
        double BSE_EDDLIMFLAG= 0;
        double BSE_BHSPINMAG= 0.0;
        double BSE_EDDFAC= 1.0;
        double BSE_GAMMA= -2.0;
        int BSE_TFLAG= 1;
        int BSE_IFFLAG= 0;
        int BSE_WDFLAG= 1;
        double BSE_EPSNOV= 0.001;
        double BSE_BCONST= -3000.00;
        double BSE_CK= -1000.00;
        double BSE_REJUV_FAC= 0.1;
	int BSE_IDUM= -999;

	bse_set_using_cmc();
	bse_set_neta(BSE_NETA);
	bse_set_bwind(BSE_BWIND);
	bse_set_hewind(BSE_HEWIND);
	bse_set_windflag(BSE_WINDFLAG);
	bse_set_eddlimflag(BSE_EDDLIMFLAG);
	bse_set_pisn(BSE_PISN);
	bse_set_ecsn(BSE_ECSN);
	bse_set_ecsn_mlow(BSE_ECSN_MLOW);
	bse_set_aic(BSE_AIC);
	bse_set_bdecayfac(BSE_BDECAYFAC);
	bse_set_st_cr(BSE_ST_CR);
	bse_set_st_tide(BSE_ST_TIDE);
	bse_set_htpmb(BSE_HTPMB);
	bse_set_rejuvflag(BSE_REJUVFLAG);
	bse_set_ussn(BSE_USSN);

	bse_set_qcrit_array(BSE_QCRIT_ARRAY, NO_BSE_QCRIT_ARRAY); 
	bse_set_fprimc_array(BSE_FPRIMC_ARRAY, NO_BSE_FPRIMC_ARRAY);
	bse_set_natal_kick_array(BSE_NATAL_KICK_ARRAY, NO_BSE_NATAL_KICK_ARRAY); 
	bse_set_sigmadiv(BSE_SIGMADIV);
	bse_set_alpha1(BSE_ALPHA1);
	bse_set_lambdaf(BSE_LAMBDAF);
	bse_set_ceflag(BSE_CEFLAG);
	bse_set_cehestarflag(BSE_CEHESTARFLAG);
	bse_set_cemergeflag(BSE_CEMERGEFLAG);
	bse_set_cekickflag(BSE_CEKICKFLAG);
	bse_set_tflag(BSE_TFLAG);
	bse_set_qcflag(BSE_QCFLAG);
        bse_set_don_lim(BSE_DON_LIM);
        bse_set_acc_lim(BSE_ACC_LIM);
	bse_set_ifflag(BSE_IFFLAG);
	bse_set_wdflag(BSE_WDFLAG);
	bse_set_bhflag(BSE_BHFLAG);
        bse_set_bhms_coll_flag(BSE_BHMS_COLL_FLAG);
        bse_set_grflag(BSE_GRFLAG);
        bse_set_kickflag(BSE_KICKFLAG);
        bse_set_zsun(BSE_ZSUN);
        bse_set_rembar_massloss(BSE_REMBAR_MASSLOSS);
	bse_set_remnantflag(BSE_REMNANTFLAG);
	bse_set_bhspinflag(BSE_BHSPINFLAG);
	bse_set_bhspinmag(BSE_BHSPINMAG);
	bse_set_mxns(BSE_MXNS);
        bse_set_wd_mass_lim(BSE_WD_MASS_LIM);
	bse_set_bconst(BSE_BCONST);
	bse_set_CK(BSE_CK);
	bse_set_rejuv_fac(BSE_REJUV_FAC);
	bse_set_idum(BSE_IDUM);
	bse_set_pts1(BSE_PTS1);
	bse_set_pts2(BSE_PTS2);
	bse_set_pts3(BSE_PTS3);
	bse_set_sigma(BSE_SIGMA);
	bse_set_bhsigmafrac(BSE_BHSIGMAFRAC);
	bse_set_polar_kick_angle(BSE_POLAR_KICK_ANGLE);
	bse_set_beta(BSE_BETA);
	bse_set_xi(BSE_XI);
	bse_set_acc2(BSE_ACC2);
	bse_set_epsnov(BSE_EPSNOV);
	bse_set_eddfac(BSE_EDDFAC);
	bse_set_gamma(BSE_GAMMA);
	bse_set_merger(-1.0);
	
	/* set parameters relating to metallicity */
	zpars = (double *) malloc(20 * sizeof(double));
	bse_zcnsts(&METALLICITY, zpars);
	
	/* set collisions matrix */
	bse_instar();
	
	for(k=1; k<=cfd->NOBJ; k++){
		if (cfd->obj_binind[k] != 0) {
			cfd->obj_k[k] = NOT_A_STAR;
			/* set binary member properties here... */
		} else {
			/* need mass in solar masses */
			star.se_mass = cfd->obj_m[k] * cfd->Mclus;
			if(star.se_mass <= 0.7){
			  star.se_k = 0;
			} else {
			  star.se_k = 1;
			}
			/* setting the rest of the variables */
			star.se_mt = star.se_mass;
			star.se_ospin = 0.0;
			star.se_epoch = 0.0;
			star.se_tphys = 0.0;

			/* evolving stars a little bit to set luminosity and radius */
			tphysf = 1.0e-6;
			dtp = tphysf;
			bse_evolve_single(&(star.se_k), &(star.se_mass), &(star.se_mt), &(star.se_radius),
				   &(star.se_lum), &(star.se_mc), &(star.se_rc), &(star.se_menv), 
				   &(star.se_renv), &(star.se_ospin), &(star.se_epoch), &(star.se_tms), 
				   &(star.se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs, &(star.se_bhspin));
			
			/* setting star properties in FITS file, being careful with units */
			cfd->obj_k[k] = star.se_k;
			cfd->obj_m[k] = star.se_mass / cfd->Mclus;
			cfd->obj_Reff[k] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);
		}
	}
}

int main(int argc, char *argv[]){
	cmc_fits_data_t cfd;
	char infilename[1024], outfilename[1024];
	int i;
	const char *short_opts = "i:o:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	sprintf(infilename, "%s", INFILE);
	sprintf(outfilename, "%s", OUTFILE);

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'i':
			sprintf(infilename, "%s", optarg);
			break;
		case 'o':
			sprintf(outfilename, "%s", optarg);
			break;
		case 'h':
			print_usage(stdout);
			return(0);
		default:
			break;
		}
	}
	
	/* check to make sure there was nothing crazy on the command line */
	if (optind < argc) {
		print_usage(stdout);
		return(1);
	}

	cmc_read_fits_file(infilename, &cfd, 0);

	METALLICITY = cfd.Z;

	stellar_evolve(&cfd);

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

