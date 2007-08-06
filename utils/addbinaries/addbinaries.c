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
#define NBIN 0
#define LIMITS 0
#define EBMIN 3.0
#define EBMAX 133.0

/* global variables needed by Star Track */
double METALLICITY, WIND_FACTOR=1.0;

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  addbinaries [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -i --infile <infile>           : input file [%s]\n", INFILE);
	fprintf(stream, "  -o --outfile <outfile>         : output file [%s]\n", OUTFILE);
	fprintf(stream, "  -N --Nbin <N_bin>              : number of binaries [%ld]\n", NBIN);
	fprintf(stream, "  -l --limits <limits algorithm> : algorithm for setting limits on binary semimajor axes (0=physical, 1=kT prescription) [%d]\n", LIMITS);
	fprintf(stream, "  -m --Ebmin <E_b,min>           : minimum binding energy, in kT [%g]\n", EBMIN);
	fprintf(stream, "  -M --Ebmax <E_b,max>           : maximum binding energy, in kT [%g]\n", EBMAX);
	fprintf(stream, "  -h --help              : display this help text\n");
}

void assign_binaries(cmc_fits_data_t *cfd, long Nbin, int limits, double Ebmin, double Ebmax)
{
	long i, j;
	double mass, Mmin, Mmax, amin, amax, W, vorb, emax, Mtotnew;
	double kTcore, Eb, Ebmin, Ebmax;

	/* pass through some header values */
	cfd2->NOBJ = cfd->NOBJ;
	if (cfd->NBINARY != 0) {
		fprintf(stderr, "Warning: NBINARY!=0 in input FITS file.  Be sure you know what you're doing!\n");
	}
	cfd2->NBINARY = Nbin;
	cfd2->Mclus = cfd->Mclus;
	cfd2->Rvir = cfd->Rvir;
	cfd2->Rtid = cfd->Rtid;
	cfd2->Z = cfd->Z;

	cmc_malloc_fits_data_t(cfd2);

	/* calculate minimum mass of mass function, plus other statistics */
	Mmin = GSL_POSINF;
	Mmax = GSL_NEGINF;
	for (i=1; i<=cfd->NOBJ; i++) {
		mass = cfd->obj_m[i];
		if (mass < Mmin) {
			Mmin = mass;
		}
		if (mass > Mmax) {
			Mmax = mass;
		}
	}

	/* first assign all secondary masses, since this will change the units */
	j = 0;
	while (j < Nbin) {
		i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
		
		/* make it a binary if it's not already */
		if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0) {
			j++;
			
			/* make this star a binary */
			star[i].binind = j;
			binary[j].inuse = 1;
			binary[j].id1 = star[i].id;
			binary[j].id2 = star_get_id_new();
			
			/* set secondary mass from dP/dq=1 distribution */
			binary[j].m1 = star[i].m;
			binary[j].m2 = Mmin + rng_t113_dbl() * (binary[j].m1 - Mmin);
			
			star[i].m = binary[j].m1 + binary[j].m2;
		}
	}
	
	/* rescale to new unit of mass (this keeps the units N-body, and also preserves the virial ratio
	   at 2T/|W|=1) */
        /* FIXME: I added the mass of a possible central object given in the fits file to Mtotnew 
         * by changing the start of the following two loops from i=1 to i=0. This is not entirely
         * consistent with how the central mass is treated in the other parts of the code, but preserves 
         * the mass units in the code if there are no binaries. When we include BH-stellar-binary 
         * interactions we need to rethink this. --- stefan: 05/23/07 */
	Mtotnew = 0.0;
	for (i=0; i<=clus.N_STAR; i++) {
		Mtotnew += star[i].m * madhoc;
	}
	
        dprintf("m[0]= %lg, Mtotnew= %lf\n", star[0].m* madhoc, Mtotnew);
	for (i=0; i<=clus.N_STAR; i++) {
		star[i].m /= Mtotnew;
		if (star[i].binind) {
			binary[star[i].binind].m1 /= Mtotnew;
			binary[star[i].binind].m2 /= Mtotnew;
		}
	}
	
	gprintf("Rescaling SOLAR_MASS_DYN from %g to", SOLAR_MASS_DYN);
	SOLAR_MASS_DYN /= Mtotnew;
	gprintf(" %g since creating binaries increases total cluster mass.\n", SOLAR_MASS_DYN);

	/* reset units */
	units_set();

	/* let's print out the mass statistics to make sure we did everything right */
	Mmin = GSL_POSINF;
	Mmax = GSL_NEGINF;
	for (i=1; i<=clus.N_STAR; i++) {
		if (star[i].m < Mmin) {
			Mmin = star[i].m;
		}
		if (star[i].m > Mmax) {
			Mmax = star[i].m;
		}
	}
	
	/* have to rescale some variables, due to new units */
	for (i=1; i<=clus.N_STAR; i++) {
		/* just single star radii, since binary stars' radii will be set below */
		if (star[i].binind == 0) {
			star[i].rad = r_of_m(star[i].m);
		}
	}

	/* calculate and store velocity dispersion profile, which we need to determine
	   the upper limit on binary semimajor axis */
	calc_sigma_r();

	/* calculate kT in cluster core */
	/* calculate core temperature to get scale for binary semimajor axis distribution */
	central_calculate();
	kTcore = (1.0/3.0) * central.m_ave * sqr(central.v_rms);
	/* kTcore = 0.0; */
	/* for (i=1; i<=NUM_CENTRAL_STARS; i++) { */
/* 		kTcore += (1.0/3.0) * (star[i].m/clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt)); */
/* 	} */
/* 	kTcore = kTcore/NUM_CENTRAL_STARS; */

	/* assign binary parameters */
	if (!BININITKT) {
		/* assign binaries physically */
		for (i=1; i<=clus.N_STAR; i++) {
			j = star[i].binind;
			if (j != 0) {
				/* set radii now that the masses are known */
				binary[j].rad1 = r_of_m(binary[j].m1);
				binary[j].rad2 = r_of_m(binary[j].m2);
				
				/* zero internal energies */
				binary[j].Eint1 = 0.0;
				binary[j].Eint2 = 0.0;
				
				/* choose a from a distribution uniform in 1/a from near contact to hard/soft boundary */
				amin = 5.0 * (binary[j].rad1 + binary[j].rad2);
				/* W = 4.0 * sigma_r(star[i].r) / sqrt(3.0 * PI); */
				/* use core velocity dispersion instead */
				W = 4.0 * central.v_rms / sqrt(3.0 * PI);
				vorb = XHS * W;
				amax = star[i].m * madhoc / sqr(vorb);
				binary[j].a = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
				
				/* get eccentricity from thermal distribution, truncated near contact */
				emax = 1.0 - amin/binary[j].a;
				binary[j].e = emax * sqrt(rng_t113_dbl());
			}
		}
	} else {
		/* assign binaries via kT description */
		/* set max and min binding energies */
		Ebmin = BININITEBMIN * kTcore;
		Ebmax = BININITEBMAX * kTcore;

		for (i=1; i<=clus.N_STAR; i++) {
			j = star[i].binind;
			if (j != 0) {
				/* set radii now that the masses are known */
				binary[j].rad1 = r_of_m(binary[j].m1);
				binary[j].rad2 = r_of_m(binary[j].m2);
				
				/* zero internal energies */
				binary[j].Eint1 = 0.0;
				binary[j].Eint2 = 0.0;
				
				/* choose binding energy uniformly in the log */
				Eb = pow(10.0, rng_t113_dbl()*(log10(Ebmax)-log10(Ebmin))+log10(Ebmin));
				binary[j].a = (binary[j].m1/clus.N_STAR) * (binary[j].m2/clus.N_STAR) / (2.0 * Eb);
				
				/* get eccentricity from thermal distribution */
				binary[j].e = sqrt(rng_t113_dbl());
			}
		}
	}

	/* update binary bookkeeping */
	M_b = 0.0;
	E_b = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		j = star[i].binind;
		if (j && binary[j].inuse) {
			M_b += star[i].m;
			E_b += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);
		}
	}

	diaprintf("done assigning binaries\n");
}

int main(int argc, char *argv[]){
	int limits;
	long Nbin;
	double Ebmin, Ebmax;
	cmc_fits_data_t cfd, cfd2;
	char infilename[1024], outfilename[1024];
	int i;
	const char *short_opts = "i:o:N:l:m:M:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"Nbin", required_argument, NULL, 'N'},
		{"limits", required_argument, NULL, 'l'},
		{"Ebmin", required_argument, NULL, 'm'},
		{"Ebmax", required_argument, NULL, 'M'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* argument defaults */
	sprintf(infilename, "%s", INFILE);
	sprintf(outfilename, "%s", OUTFILE);
	Nbin = NBIN;
	limits = LIMITS;
	Ebmin = EBMIN;
	Ebmax = EBMAX;

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'i':
			sprintf(infilename, "%s", optarg);
			break;
		case 'o':
			sprintf(outfilename, "%s", optarg);
			break;
		case 'N':
			Nbin = strtol(optarg, NULL, 10);
			if (Nbin < 0) {
				fprintf(stderr, "Nbin must be >=0!\n");
				exit(1);
			}
			break;
		case 'l':
			limits = strtol(optarg, NULL, 10);
			if (limits != 0 && limits != 1) {
				fprintf(stderr, "limits must be either 0 or 1!\n");
				exit(1);
			}
			break;
		case 'm':
			Ebmin = strtod(optarg, NULL);
			if (Ebmin <= 0.0) {
				fprintf(stderr, "Ebmin must be positive!\n");
				exit(1);
			}
			break;
		case 'M':
			Ebmax = strtod(optarg, NULL);
			if (Ebmax <= 0.0) {
				fprintf(stderr, "Ebmax must be positive!\n");
				exit(1);
			}
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

	assign_binaries(&cfd, &cfd2, Nbin, limits, Ebmin, Ebmax)

	cmc_write_fits_file(&cfd2, outfilename);

	cmc_free_fits_data_t(&cfd);
	cmc_free_fits_data_t(&cfd2);

	return 0;
}

