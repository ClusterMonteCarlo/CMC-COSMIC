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
#include "../../startrack/sinbin.h"

#define INFILE "in.fits"
#define OUTFILE "debug.fits"

/* global variables needed by Star Track */
double METALLICITY, WIND_FACTOR=1.0;

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
	double frac;
	star_t star;

	/* the following function calls need to be done for
	 * Chris' stellar evolution code */
	/* -- begin black_magic -- */
	M_hook=M_hookf();    
	M_HeF=M_HeFf();
	M_FGB=M_FGBf();
	coef_aa();
	coef_bb(); 
	/* --  end  black_magic -- */
	
	for(k=1; k<=cfd->NOBJ; k++){
		if (cfd->obj_binind[k] != 0) {
			cfd->obj_k[k] = NOT_A_STAR;
			/* set binary member properties here... */
		} else {
			/* need mass in solar masses */
			star.mass = cfd->obj_m[k] * cfd->Mclus;
			if(star.mass <= 0.7){
				star.k = S_MASS_MS_STAR;
			} else {
				star.k = L_MASS_MS_STAR;
			}
			/* setting the rest of the variables */
			star.mzams = star.m0 = star.mass;
			star.tbeg = star.tvir = 0.0;
			star.tend = star.tbeg + 1e-11;
			star.mc = star.mcHe = star.mcCO = 0.0;
			star.dt = star.mpre = star.tstart = frac = 0.0;
			star.kpre = star.k;
			star.flag = 0;
			star.lum = star.rad = 1.0;

			/* evolving stars a little bit to set luminosity and radius */
			singl(&(star.mzams), &(star.m0), &(star.mass), 
			      &(star.k), &(star.tbeg), &(star.tvir),
			      &(star.tend), &(star.lum), &(star.rad),
			      &(star.mc), &(star.mcHe), &(star.mcCO),
			      &(star.flag), &(star.dt), &(star.mpre),
			      &(star.kpre), &(star.tstart), &frac);
			
			/* setting star properties in FITS file, being careful with units */
			cfd->obj_k[k] = star.k;
			cfd->obj_m[k] = star.mass / cfd->Mclus;
			cfd->obj_Reff[k] = star.rad / (cfd->Rvir * PARSEC / RSUN);
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

