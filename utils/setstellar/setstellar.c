#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "../../cmc.h"
#include "../../common/fitslib.h"
#include "../../common/taus113-v2.h"
#include "../../bse_wrap/bse_wrap.h"

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
	double vs[3];

	/* initialize stellar evolution */
	/* SSE */
	/* bse_set_hewind(0.5); */
	
	/* BSE */
	bse_set_neta(0.5);
	bse_set_bwind(0.0);
	bse_set_hewind(1.0);
	bse_set_alpha1(3.0);
	bse_set_lambda(0.5);
	bse_set_ceflag(0);
	bse_set_tflag(1);
	bse_set_ifflag(0);
	bse_set_wdflag(1);
	bse_set_bhflag(0);
	bse_set_nsflag(1);
	bse_set_mxns(3.0);
	bse_set_idum(29769);
	bse_set_pts1(0.05);
	bse_set_pts2(0.01);
	bse_set_pts3(0.02);
	bse_set_sigma(190.0);
	bse_set_beta(0.125);
	bse_set_xi(1.0);
	bse_set_acc2(1.5);
	bse_set_epsnov(0.001);
	bse_set_eddfac(10.0);
	bse_set_gamma(-1.0);
	
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
			bse_evolv1(&(star.se_k), &(star.se_mass), &(star.se_mt), &(star.se_radius), 
				   &(star.se_lum), &(star.se_mc), &(star.se_rc), &(star.se_menv), 
				   &(star.se_renv), &(star.se_ospin), &(star.se_epoch), &(star.se_tms), 
				   &(star.se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);
			
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

	cmc_read_fits_file(infilename, &cfd);

	METALLICITY = cfd.Z;

	stellar_evolve(&cfd);

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

