#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "../../common/fitslib.h"
#include "../../common/taus113-v2.h"

#define INFILE "in.fits"
#define OUTFILE "debug.fits"

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  setunits [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -i --infile <infile>   : input file [%s]\n", INFILE);
	fprintf(stream, "  -o --outfile <outfile> : output file [%s]\n", OUTFILE);
	fprintf(stream, "  -M --Mclus <M_clus>    : set cluster mass in M_sun [from infile]\n");
	fprintf(stream, "  -R --Rvir <R_vir>      : set virial radius in parsecs [from infile]\n");
	fprintf(stream, "  -T --Rtid <R_tid>      : set tidal radius in N-body units [from infile]\n");
	fprintf(stream, "  -Z --Z <Z>             : set metallicity [from infile]\n");
	fprintf(stream, "  -h --help              : display this help text\n");
}

int main(int argc, char *argv[]){
	double Mclus=0.0, Rvir=0.0, Rtid=0.0, Z=0.0;
	int setMclus=0, setRvir=0, setRtid=0, setZ=0;
	cmc_fits_data_t cfd;
	char infilename[1024], outfilename[1024];
	int i;
	const char *short_opts = "i:o:M:R:T:Z:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"Mclus", required_argument, NULL, 'M'},
		{"Rvir", required_argument, NULL, 'R'},
		{"Rtid", required_argument, NULL, 'T'},
		{"Z", required_argument, NULL, 'Z'},
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
		case 'M':
			Mclus = strtod(optarg, NULL);
			setMclus = 1;
			break;
		case 'R':
			Rvir = strtod(optarg, NULL);
			setRvir = 1;
			break;
		case 'T':
			Rtid = strtod(optarg, NULL);
			setRtid = 1;
			break;
		case 'Z':
			Z = strtod(optarg, NULL);
			setZ = 1;
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

	if (setMclus) cfd.Mclus = Mclus;
	if (setRvir) cfd.Rvir = Rvir;
	if (setRtid) cfd.Rtid = Rtid;
	if (setZ) cfd.Z = Z;

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

