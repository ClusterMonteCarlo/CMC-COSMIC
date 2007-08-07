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
#define NBIN 0L
#define LIMITS 0
#define EBMIN 3.0
#define EBMAX 133.0
#define SEED 8732UL

/* global variables needed by Star Track */
double METALLICITY, WIND_FACTOR=1.0;
/* global variable for new star id */
long newstarid;

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
	fprintf(stream, "  -s --seed <seed>               : random seed [%ld]\n", SEED);
	fprintf(stream, "  -h --help                      : display this help text\n");
}

long star_get_id_new(void)
{
	newstarid++;
	return(newstarid);
}

/* calculate and store the velocity dispersion profile */
void calc_sigma_r(cmc_fits_data_t *cfd, double *r, double *sigma, double *mave)
{
	long si, k, p=AVEKERNEL, simin, simax;
	double Mv2ave, Mave, Sigma;

	/* calculate sliding average */
	for (si=1; si<=cfd->NOBJ; si++) {
		simin = si - p;
		simax = simin + (2 * p - 1);
		if (simin < 1) {
			simin = 1;
			simax = simin + (2 * p - 1);
		} else if (simax > cfd->NOBJ) {
			simax = cfd->NOBJ;
			simin = simax - (2 * p - 1);
		}

		Mv2ave = 0.0;
		Mave = 0.0;
		for (k=simin; k<=simax; k++) {
			Mv2ave += cfd->obj_m[k] * (cfd->obj_vr[k] * cfd->obj_vr[k] + cfd->obj_vt[k] * cfd->obj_vt[k]);
			Mave += cfd->obj_m[k];
		}
		Mv2ave /= (double) (2 * p);
		Mave /= (double) (2 * p);
		
		/* Sigma is the 3D velocity dispersion */
		Sigma = sqrt(Mv2ave/Mave);
		
		/* store sigma */
		r[si] = cfd->obj_r[si];
		sigma[si] = Sigma;
		mave[si] = Mave;
	}
}

void assign_binaries(cmc_fits_data_t *cfd, long Nbin, int limits, double EbminkT, double EbmaxkT)
{
	long i, j;
	double mass, Mmin, Mmax, amin, amax, W, vorb, emax, Mtotnew;
	double kTcore, vcore, Eb, Ebmin, Ebmax;
	double frac, *r, *sigma, *mave;
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

	/* set newstarid to emulate CMC function */
	newstarid = cfd->NOBJ;

	/* test for strange things */
	if (Nbin > cfd->NOBJ) {
		fprintf(stderr, "Error: You've requested more binaries than there are objects in the input file!\n");
		exit(1);
	}

	if (cfd->NBINARY != 0) {
		fprintf(stderr, "Warning: NBINARY!=0 in input FITS file.  Be sure you know what you're doing!\n");
	}

	cfd->NBINARY = Nbin;
	
	/* reallocate memory for binaries */
	cfd->bs_index = (long *) realloc(cfd->bs_index, (Nbin+1)*sizeof(long));
	cfd->bs_id1 = (long *) realloc(cfd->bs_id1, (Nbin+1)*sizeof(long));
	cfd->bs_k1 = (int *) realloc(cfd->bs_k1, (Nbin+1)*sizeof(int));
	cfd->bs_m1 = (double *) realloc(cfd->bs_m1, (Nbin+1)*sizeof(double));
	cfd->bs_Reff1 = (double *) realloc(cfd->bs_Reff1, (Nbin+1)*sizeof(double));
	cfd->bs_id2 = (long *) realloc(cfd->bs_id2, (Nbin+1)*sizeof(long));
	cfd->bs_k2 = (int *) realloc(cfd->bs_k2, (Nbin+1)*sizeof(int));
	cfd->bs_m2 = (double *) realloc(cfd->bs_m2, (Nbin+1)*sizeof(double));
	cfd->bs_Reff2 = (double *) realloc(cfd->bs_Reff2, (Nbin+1)*sizeof(double));
	cfd->bs_a = (double *) realloc(cfd->bs_a, (Nbin+1)*sizeof(double));
	cfd->bs_e = (double *) realloc(cfd->bs_e, (Nbin+1)*sizeof(double));

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

			/* make this object a binary */
			cfd->obj_binind[i] = j;
			cfd->bs_index[j] = j;

			/* copy properties for star 1 */
			cfd->bs_id1[j] = cfd->obj_id[i];
			cfd->bs_k1[j] = cfd->obj_k[i];
			cfd->bs_m1[j] = cfd->obj_m[i];
			cfd->bs_Reff1[j] = cfd->obj_Reff[i];

			/* assign properties for star 2;
			   this assumes all stars in input file are id'ed sequentially from 1 */
			cfd->bs_id2[j] = star_get_id_new();
			/* set secondary mass from dP/dq=1 distribution */
			cfd->bs_m2[j] = Mmin + rng_t113_dbl() * (cfd->bs_m1[j] - Mmin);
			
			/* stellar evolution stuff */
			star.mass = cfd->bs_m2[j] * cfd->Mclus;
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
			cfd->bs_k2[j] = star.k;
			cfd->bs_m2[j] = star.mass / cfd->Mclus;
			cfd->bs_Reff2[j] = star.rad / (cfd->Rvir * PARSEC / RSUN);

			/* set/unset stellar properties for obj */
			cfd->obj_id[i] = NOT_A_STAR;
			cfd->obj_k[i] = NOT_A_STAR;
			cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
			cfd->obj_Reff[i] = 0.0;
		}
	}
	
	/* rescale to new unit of mass (this keeps the units N-body, and also preserves the virial ratio
	   at 2T/|W|=1) */
        /* FIXME: I added the mass of a possible central object given in the fits file to Mtotnew 
         * by changing the start of the following two loops from i=1 to i=0. This is not entirely
         * consistent with how the central mass is treated in the other parts of the code, but preserves 
         * the mass units in the code if there are no binaries. When we include BH-stellar-binary 
         * interactions we need to rethink this. --- stefan: 05/23/07 */
	/* FIXME: Not sure how this should work with the new FITS format, but I'm leaving in Stefan's
	   change. --- John: 8/7/07
	 */
	/* Mtot is in N-body units */
	Mtotnew = 0.0;
	for (i=0; i<=cfd->NOBJ; i++) {
		Mtotnew += cfd->obj_m[i];
	}
	
	cfd->Mclus *= Mtotnew;
	
	for (i=0; i<=cfd->NOBJ; i++) {
		cfd->obj_m[i] /= Mtotnew;
	}
	
	for (i=1; i<=cfd->NBINARY; i++) {
	  cfd->bs_m1[i] /= Mtotnew;
	  cfd->bs_m2[i] /= Mtotnew;
	}
	
	/* calculate and store velocity dispersion profile, which we need to determine
	   the upper limit on binary semimajor axis */
	r = (double *) malloc((cfd->NOBJ+1)*sizeof(double));
	sigma = (double *) malloc((cfd->NOBJ+1)*sizeof(double));
	mave = (double *) malloc((cfd->NOBJ+1)*sizeof(double));
	calc_sigma_r(cfd, r, sigma, mave);

	/* calculate kT in cluster core to get scale for binary semimajor axis distribution */
	kTcore = 0.0;
	vcore = 0.0;
	for (i=1; i<=AVEKERNEL; i++) {
		vcore += sigma[i];
		kTcore += (1.0/3.0) * mave[i] * sigma[i] * sigma[i];
	}
	vcore /= AVEKERNEL;
	kTcore /= AVEKERNEL;

	/* assign binary parameters */
	if (!limits) {
		/* assign binaries physically */
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose a from a distribution uniform in 1/a from near contact to hard/soft boundary */
			amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
			W = 4.0 * vcore / sqrt(3.0 * PI);
			vorb = XHS * W;
			amax = cfd->obj_m[i] / (vorb * vorb);
			cfd->bs_a[i] = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
			
			/* get eccentricity from thermal distribution, truncated near contact */
			emax = 1.0 - amin / cfd->bs_a[i];
			cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
		}
	} else {
		/* assign binaries via kT description */
		/* set max and min binding energies */
		Ebmin = EbminkT * kTcore;
		Ebmax = EbminkT * kTcore;

		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose binding energy uniformly in the log */
			Eb = pow(10.0, rng_t113_dbl()*(log10(Ebmax)-log10(Ebmin))+log10(Ebmin));
			cfd->bs_a[i] = cfd->bs_m1[i] * cfd->bs_m2[i] / (2.0 * Eb);
			
			/* get eccentricity from thermal distribution */
			cfd->bs_e[i] = sqrt(rng_t113_dbl());
		}
	}

	free(sigma);
	free(r);
	free(mave);
}

int main(int argc, char *argv[]){
	int i, limits;
	long Nbin;
	unsigned long seed;
	double Ebmin, Ebmax;
	cmc_fits_data_t cfd;
	char infilename[1024], outfilename[1024];
	const char *short_opts = "i:o:N:l:m:M:s:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"Nbin", required_argument, NULL, 'N'},
		{"limits", required_argument, NULL, 'l'},
		{"Ebmin", required_argument, NULL, 'm'},
		{"Ebmax", required_argument, NULL, 'M'},
		{"seed", required_argument, NULL, 's'},
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
	seed = SEED;

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
		case 's':
			seed = strtol(optarg, NULL, 10);
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

	reset_rng_t113(seed);

	cmc_read_fits_file(infilename, &cfd);

	METALLICITY = cfd.Z;

	assign_binaries(&cfd, Nbin, limits, Ebmin, Ebmax);

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}
