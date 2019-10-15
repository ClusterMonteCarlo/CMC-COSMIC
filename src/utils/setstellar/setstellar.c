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
	double vs[12];

	/* initialize stellar evolution */
	/* SSE */
	/* bse_set_hewind(0.5); */
	
	/* BSE */
        bse_set_using_cmc();
        bse_set_neta(0.5);
        bse_set_bwind(0.0);
        bse_set_hewind(1.0);
        bse_set_windflag(3);
        bse_set_eddlimflag(0);
        bse_set_pisn(45.0);
        bse_set_ecsn(2.5);
        bse_set_ecsn_mlow(1.4);
        bse_set_aic(1);
        bse_set_bdecayfac(1);
        bse_set_st_cr(0);
        bse_set_st_tide(1);
        bse_set_htpmb(1);
        bse_set_rejuvflag(0);
        bse_set_ussn(0);
        double bse_qcrit_array[16];
        double bse_fprimc_array[16];
        double bse_natal_kick_array[6];
        int i;
        for(i=0;i<16;i++) bse_qcrit_array[i] = 0.0;
        for(i=0;i<16;i++) bse_fprimc_array[i] = 2.0/21.0;
        for(i=0;i<6;i++) bse_natal_kick_array[i] = -100.0;
        bse_set_qcrit_array(bse_qcrit_array, 16);
        bse_set_fprimc_array(bse_fprimc_array, 16);
        bse_set_natal_kick_array(bse_natal_kick_array, 6);
        bse_set_sigmadiv(-20.0);
        bse_set_alpha1(1.0); /* FIXME: is 3 too high? (normally 1.0) */
        bse_set_lambdaf(0.5);
        bse_set_ceflag(0);
        bse_set_cehestarflag(0);
        bse_set_cemergeflag(0);
        bse_set_cekickflag(2);
        bse_set_tflag(1);
        bse_set_qcflag(2);
        bse_set_ifflag(0);
        bse_set_wdflag(0);
        bse_set_bhflag(1);
        bse_set_nsflag(3);
        bse_set_bhspinflag(0);
        bse_set_bhspinmag(0.0);
        bse_set_mxns(2.5); //3 if nsflag=1 or 2, 1.8 if nsflag=0 (see evolv2.f)
        bse_set_bconst(-3000.0);
        bse_set_CK(-1000.0);
        bse_set_rejuv_fac(0.1);

        bse_set_idum(22);
        bse_set_pts1(0.001);
        bse_set_pts2(0.01);
        bse_set_pts3(0.02);
        bse_set_sigma(265.0);
        bse_set_bhsigmafrac(1.0);
        bse_set_polar_kick_angle(90.0);
        bse_set_beta(-0.125); //set -0.125 if variable beta (following startrack), otherwise 0.125 for bse.
        bse_set_xi(0.5);
        bse_set_acc2(1.5);
        bse_set_epsnov(0.001);
        bse_set_eddfac(1.0); /* (normally 1.0) */
        bse_set_gamma(-1.0);
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

