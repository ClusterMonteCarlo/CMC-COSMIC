/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#include "cmc_vars.h"

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** %s %s (%s) [%s] **\n", PRETTYNAME, VERSION, NICK, DATE);
}

/* print the usage */
void print_usage(FILE *stream, char *argv[])
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  %s [options...] <input_file> <output_file_prefix>\n", argv[0]);
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -d --debug   : turn on debugging\n");
	fprintf(stream, "  -V --version : print version info\n");
	fprintf(stream, "  -h --help    : display this help text\n");
	fprintf(stream, "\n");
	fprintf(stream, "OUTPUT:\n");
	fprintf(stream, "\toutput_file_prefix_0 = log file\n");
	fprintf(stream, "\toutput_file_prefix_1 = lagrange radii\n");
	fprintf(stream, "\toutput_file_prefix_3 = total energy, virial ratio\n");
	fprintf(stream, "\toutput_file_prefix_4 = component-wise fraction of total mass\n");
}
void print_results(void){
	PrintLogOutput();
	PrintFileOutput();
	fflush(NULL);
}

/*********** Output 2D/3D snapshots **************/
void print_2Dsnapshot(void)
{
	long i;
	char outfile[100];

	/* open file for 2D snapshot */
	sprintf(outfile, "%ssnap%04ld.dat.gz", outprefix, snap_num);
	if ((snapfile = gzopen(outfile, "wb")) == NULL) {
		eprintf("cannot create 2D snapshot file %s\n", outfile);
		exit_cleanly(1);
	}

	/* print useful header */
	gzprintf(snapfile, "# t=%.8g [code units]\n", TotalTime);
	gzprintf(snapfile, "# id  m [MSUN]  r [code units]  vr [code units]  vt [code units]  E [code units]  J [code units]\n");

	/* then print data */
	for (i=1; i<=clus.N_MAX; i++) {
		gzprintf(snapfile, "%ld %.8g %.8g %.8g %.8g %.8g %.8g\n", 
			 star[i].id, star[i].m * initial_total_mass, star[i].r, star[i].vr, star[i].vt, star[i].E, star[i].J);
	}

	gzclose(snapfile);

	/* global counter for snapshot output file */
	snap_num++;
}

void PrintLogOutput(void)
{
	double m, rh, trh, conc_param, m_single, m_binary;
	long ih, k;
	
	/* Computing half-mass radii, and relaxation time */
	m = rh = trh = 0.0;
	rh_binary = rh_single = m_binary = m_single = 0.0;
	for (ih=1; ih<=clus.N_MAX; ih++) {
		k = ih;
		m += star[k].m / clus.N_STAR;
		
		if (star[k].binind > 0) {
			m_binary += star[k].m / clus.N_STAR;
		} else {
			m_single += star[k].m / clus.N_STAR;
		}
		
		if (m/Mtotal <= 0.5) {
			rh = star[k].r;
		}
		if (m_single / (Mtotal - (M_b / clus.N_STAR)) <= 0.5) {
			rh_single = star[k].r;
		}
		/* avoid dividing by zero if there are no binaries */
		if (M_b > 0) {
			if (m_binary / M_b * clus.N_STAR <= 0.5) {
				rh_binary = star[k].r;
			}
		}
	}
	
	/* t_rh calculated using r_h */
	trh = ((0.138 * clus.N_MAX) / log((double) clus.N_MAX)) * sqrt((rh * rh * rh) / Mtotal) * log((double) clus.N_MAX) / clus.N_MAX;

	/* Concentration parameter --- note that max_r is the max radius of all bound
	   stars, returned by get_positions(). When a finite R_MAX (i.e. initial Rtidal)
	   is specified, max_r should be approximately equal to Rtidal. But when 
	   Rtidal is very large (isolated cluster), r_max still gives the "actual"
	   size of the cluster, indicating the maximum radius.
	   The conc parameter is actually defined as LOG_10(rtidal/rcore). But for
	   easier reading, it is calculated here only as rtidal/rcore. */
	
	/* max_r & core_radius are output in out_3 file */
	if(max_r == 0.0){
		conc_param = 0.0;
	} else {
		conc_param = (max_r / core_radius);
	}

	fprintf(stdout, "******************************************************************************\n");
	fprintf(logfile, "******************************************************************************\n");

	fprintf(stdout, "tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);
	fprintf(logfile, "tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);

	fprintf(stdout, "sub.count=%ld sub.FACTOR=%ld sub.N_MAX=%ld sub.rmax=%g\n", sub.count, sub.FACTOR, sub.N_MAX, sub.rmax);
	fprintf(logfile, "sub.count=%ld sub.FACTOR=%ld sub.N_MAX=%ld sub.rmax=%g\n", sub.count, sub.FACTOR, sub.N_MAX, sub.rmax);
	
	fprintf(stdout, "Etotal=%g max_r=%g N(bound)=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);
	fprintf(logfile, "Etotal=%g max_r=%g N(bound)=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);
	
	fprintf(stdout, "Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);
	fprintf(logfile, "Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);
	
	fprintf(stdout, "TidalMassLoss=%g\n", TidalMassLoss);
	fprintf(logfile, "TidalMassLoss=%g\n", TidalMassLoss);
	
	fprintf(stdout, "core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	fprintf(logfile, "core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	
	fprintf(stdout, "Sin2Beta=%g\n", Sin2Beta);
	fprintf(logfile, "Sin2Beta=%g\n", Sin2Beta);
	
	fprintf(stdout, "trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	fprintf(logfile, "trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	
	fprintf(stdout, "N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);
	fprintf(logfile, "N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);

	fprintf(stdout, "******************************************************************************\n");
	fprintf(logfile, "******************************************************************************\n");
}

void PrintFileOutput(void) {
	long i, n_single, n_binary;
	double fb, fb_core;

	/* print lagrangian radii */
	/* Also note that there are only MASS_PC_COUNT-1 radii 
	 * from 0...MASS_PC_COUNT-1 
	 * So output only needs to go over i = 0 to MASS_PC_COUNT-2. */
	fprintf(out[0], "%.9e  ", TotalTime);
	fprintf(ave_mass_file, "%.9e  ",TotalTime);
	fprintf(no_star_file, "%.9e  ",TotalTime);
	fprintf(densities_file, "%.9e  ",TotalTime);
	fprintf(centmass_file, "%.9e  %.9e %.9e %.9e %.9e %.9e %.9e\n",
					TotalTime,
					cenma.m/clus.N_STAR/Mtotal, 
					Dt, rho_core,
					Etotal.tot, Etotal.K, Etotal.P);

	for (i = 0; i < MASS_PC_COUNT - 1; i++) {
		fprintf(out[0], "%e   ", mass_r[i]);
		fprintf(ave_mass_file,"%e ", ave_mass_r[i]);
		fprintf(no_star_file,"%g ", no_star_r[i]);
		fprintf(densities_file,"%e ", densities_r[i]);
	}
	fprintf(out[0], "\n");
	fprintf(ave_mass_file,"\n");
	fprintf(no_star_file,"\n");
	fprintf(densities_file,"\n");
	
	/* output Time,N_MAX,TotalE,TotalKE,TotalPE,Mtotal */
	fprintf(out[1], "%.8G  %8ld  %.8G  %.8G  %.8G  %.8G  %.8G  %.8G  %.8G  %8ld  %.8G  %.8G %.8G %.8G %.8G %.8G\n",
		TotalTime, clus.N_MAX, Etotal.tot, Etotal.K, Etotal.P, Mtotal, Etotal.New, Eescaped, Jescaped,
		tcount, max_r, core_radius, N_core, star[1].gravity, star[1].m, Ebescaped);
	
	/* Output binary data Note: N_BINARY counts ALL binaries (including escaped/destroyed ones)
	   whereas N_b only counts EXISTING BOUND binaries. */
	if (clus.N_BINARY > 0) {
		/* calculate core binary fraction */
		n_single = 0;
		n_binary = 0;
		for (i=1; star[i].r<=core_radius; i++) {
			if (star[i].binind > 0) {
				n_binary++;
			} else {
				n_single++;
			}
		}
		/* this is such a kludge: core_radius is not initialized on the first timestep */
		if (n_single + n_binary == 0) {
			fb_core = 0.0;
		} else {
			fb_core = ((double) n_binary)/((double) (n_single + n_binary));
		}
		
		/* calculate overall binary fraction */
		n_single = 0;
		n_binary = 0;
		for (i=1; i<=clus.N_MAX; i++) {
			if (star[i].binind > 0) {
				n_binary++;
			} else {
				n_single++;
			}
		}
		/* this is such a kludge: core_radius is not initialized on the first timestep */
		if (n_single + n_binary == 0) {
			fb = 0.0;
		} else {
			fb = ((double) n_binary)/((double) (n_single + n_binary));
		}
		
		/* print to file */
		fprintf(binaryfile,
			"%.6g %ld %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %ld %ld %.6g %.6g\n",
			TotalTime, N_b, M_b, E_b, rh_single, 
			rh_binary, rho_core_single, rho_core_bin, Delta_BE_bb, Delta_BE_bs, 
			DE_bb, DE_bs, N_bb, N_bs, fb_core, 
			fb);
	}

	/* also saves INITIAL snapshot (StepCount=0) */
	if (TotalTime >= T_PRINT_STEP * StepCount) {
		StepCount++;
		if (DUMPS == 1) {
			print_2Dsnapshot();
		}
	}
}

/*** Parsing of Input Parameters / Memory allocation / File I/O ***/
int parser(int argc, char *argv[], gsl_rng *r)
{
	char inputfile[1024], outfile[1024], outfilemode[5];
	char parameter_name[1024], values[1024], dummy[1024], line[2048];
	char *curr_mass;
	parsed_t parsed;
	parsed_t *spp;
	int i, allparsed=1;
	/* int *ip; */
	FILE *in;
	const char *short_opts = "dVh";
	const struct option long_opts[] = {
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* set parameters to default values */
	debug = 0;
	
	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'd':
			debug = 1;
			break;
		case 'V':
			print_version(stdout);
			exit(0);
		case 'h':
			print_version(stdout);
			fprintf(stdout, "\n");
			print_usage(stdout, argv);
			exit(0);
		default:
			break;
		}
	}

	/* check to make sure there was nothing crazy on the command line */
	if (argc - optind != 2) {
		print_usage(stdout, argv);
		exit(0);
	}
	
	/* set inputfile and outprefix now that the options have been parsed */
	sprintf(inputfile, "%s", argv[optind]);
	sprintf(outprefix, "%s_", argv[optind+1]);
	dprintf("inputfile=%s outprefix=%s\n", inputfile, outprefix);

	/*======= Opening of input & output files =======*/
	if ((in = fopen(inputfile, "r")) == NULL) {
		eprintf("Cannot open input file \"%s\".\n", inputfile);
		exit(1);
	}

	/* nothing is set yet (you rock Ato - nice code) */
	/* in C it is a complicated for loop that depends on the fact that
	 * all members of the struct is int, this is the fun way */
	//ip = (int *) &parsed;
	//for(i=0; i<(sizeof(parsed)/sizeof(int)); i++){
	//	*(ip++) = 0;
	//}
	/* faster and portable way to do this is using calloc, this does
	 * not depend on the structure of the struct */
	spp = calloc(1, sizeof(parsed_t));
	parsed = *spp;
	free(spp);
	
	while (fgets(line, 2048, in)) {
		/* the maximum size of an input line is limited by 
		 * the size of the line array */
		if (strlen(line) == 2047) {
			eprintf("input line too long: \"%s\".\n", line);
			exit(1);
		}
		/* strip comments and newline character */
		line[strcspn(line, "#\n")] = '\0';
		
		/* replace equal sign with a space */
		while (strchr(line, '=') != NULL) {
			line[strcspn(line, "=")] = ' ';
		}
		
		/* convert lowercase letters to uppercase in the keyword only */
		i = 0;
		while (!(isspace(line[i]))) {
			line[i] = toupper(line[i]);
			i++;
		}

		/* see if there are too many values for parameter */
		if (sscanf(line, "%s %s %s", parameter_name, values, dummy) == 3) {
			eprintf("too many values for parameter: \"%s\".\n", line);
			exit(1);
		} else if (sscanf(line, "%s %s", parameter_name, values) == 2) {
			if (strcmp(parameter_name, "MMIN") == 0) {
				sscanf(values, "%lf", &MMIN);
				parsed.MMIN = 1;
			} else if (strcmp(parameter_name, "BINBIN") == 0) {
				sscanf(values, "%d", &BINBIN);
				parsed.BINBIN = 1;
			} else if (strcmp(parameter_name, "BINBIN_FEWBODY") == 0) {
				sscanf(values, "%d", &BINBIN_FEWBODY);
				parsed.BINBIN_FEWBODY = 1;
			} else if (strcmp(parameter_name, "BINSINGLE") == 0) {
				wprintf("parameter BINSINGLE not yet implemented...\n");
				sscanf(values, "%d", &BINSINGLE);
				parsed.BINSINGLE = 1;
			} else if (strcmp(parameter_name, "BINSINGLE_FEWBODY") == 0) {
				sscanf(values, "%d", &BINSINGLE_FEWBODY);
				parsed.BINSINGLE_FEWBODY = 1;
			} else if (strcmp(parameter_name, "CENTRAL_MASS") == 0) {
				sscanf(values, "%lf", &cenma.m);
				cenma.E = 0.0;
				parsed.CENTRAL_MASS = 1;
			} else if (strcmp(parameter_name, "DT_FACTOR") == 0) {
				sscanf(values, "%lf", &DT_FACTOR);
				parsed.DT_FACTOR = 1;
			} else if (strcmp(parameter_name, "DUMPS") == 0) {
				sscanf(values, "%ld", &DUMPS);
				parsed.DUMPS = 1;
			} else if (strcmp(parameter_name, "E_CONS") == 0) {
				sscanf(values, "%ld", &E_CONS);
				parsed.E_CONS = 1;
			} else if (strcmp(parameter_name, "IDUM") == 0) {
				sscanf(values, "%ld", &IDUM);
				parsed.IDUM = 1;
			} else if (strcmp(parameter_name, "INDEX_UNIT") == 0) {
				sscanf(values, "%ld", &INDEX_UNIT);
				parsed.INDEX_UNIT = 1;
			} else if (strcmp(parameter_name, "INPUT_FILE") == 0) {
				sscanf(values, "%s", INPUT_FILE);
				parsed.INPUT_FILE = 1;
			} else if (strcmp(parameter_name, "MASS_PC") == 0) {
				strcpy(MASS_PC, values);
				curr_mass = (char *) strtok(values, ",; ");
				sscanf(curr_mass, "%ld", &NUM_MASS_RADII_BINS);
				for (MASS_PC_COUNT = 1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; MASS_PC_COUNT++);
				parsed.MASS_PC = 1;
			} else if (strcmp(parameter_name, "MAX_INDEX") == 0) {
				sscanf(values, "%ld", &MAX_INDEX);
				parsed.MAX_INDEX = 1;
			} else if (strcmp(parameter_name, "MEGA_YEAR") == 0) {
				sscanf(values, "%lf", &MEGA_YEAR);
				parsed.MEGA_YEAR = 1;
			} else if (strcmp(parameter_name, "METALLICITY") == 0) {
				sscanf(values, "%lf", &METALLICITY);
				parsed.METALLICITY = 1;
			} else if (strcmp(parameter_name, "MINIMUM_R") == 0) {
				sscanf(values, "%lf", &MINIMUM_R);
				parsed.MINIMUM_R = 1;
			} else if (strcmp(parameter_name, "MIN_LAGRANGIAN_RADIUS") == 0) {
				sscanf(values, "%lf", &MIN_LAGRANGIAN_RADIUS);
				parsed.MIN_LAGRANGIAN_RADIUS = 1;
			} else if (strcmp(parameter_name, "N_BINARY") == 0) {
				sscanf(values, "%ld", &clus.N_BINARY);
				parsed.N_BINARY = 1;
			} else if (strcmp(parameter_name, "NUM_CORE_STARS") == 0) {
				sscanf(values, "%ld", &NUM_CORE_STARS);
				parsed.NUM_CORE_STARS = 1;
			} else if (strcmp(parameter_name, "ORIGINAL_PERTURB_STARS") == 0) {
				sscanf(values, "%d", &ORIGINAL_PERTURB_STARS);
				parsed.ORIGINAL_PERTURB_STARS = 1;
			} else if (strcmp(parameter_name, "PERTURB") == 0) {
				sscanf(values, "%ld", &PERTURB);
				parsed.PERTURB = 1;
			} else if (strcmp(parameter_name, "R_MAX") == 0) {
				sscanf(values, "%lf", &R_MAX);
				parsed.R_MAX = 1;
			} else if (strcmp(parameter_name, "SIN2BETA_MAX") == 0) {
				sscanf(values, "%lf", &SIN2BETA_MAX);
				parsed.SIN2BETA_MAX = 1;
			} else if (strcmp(parameter_name, "SOLAR_MASS_DYN") == 0) {
				sscanf(values, "%lf", &SOLAR_MASS_DYN);
				parsed.SOLAR_MASS_DYN = 1;
			} else if (strcmp(parameter_name, "STELLAR_EVOLUTION") == 0) {
				sscanf(values, "%ld", &STELLAR_EVOLUTION);
				parsed.STELLAR_EVOLUTION = 1;
			} else if (strcmp(parameter_name, "TERMINAL_ENERGY_DISPLACEMENT") == 0) {
				sscanf(values, "%lf", &TERMINAL_ENERGY_DISPLACEMENT);
				parsed.TERMINAL_ENERGY_DISPLACEMENT = 1;
			} else if (strcmp(parameter_name, "T_MAX") == 0) {
				sscanf(values, "%lf", &T_MAX);
				parsed.T_MAX = 1;
			} else if (strcmp(parameter_name, "T_MAX_COUNT") == 0) {
				sscanf(values, "%ld", &T_MAX_COUNT);
				parsed.T_MAX_COUNT = 1;
			} else if (strcmp(parameter_name, "T_PRINT_STEP") == 0) {
				sscanf(values, "%lf", &T_PRINT_STEP);
				parsed.T_PRINT_STEP = 1;
			} else if (strcmp(parameter_name, "WIND_FACTOR") == 0) {
				sscanf(values, "%lf", &WIND_FACTOR);
				parsed.WIND_FACTOR = 1;
			} else if (strcmp(parameter_name, "GAMMA") == 0) {
				sscanf(values, "%lf", &GAMMA);
				parsed.GAMMA = 1;
			} else {
				wprintf("unknown parameter: \"%s\".\n", line);
			}
		} else if (sscanf(line, "%s", parameter_name) == 1) {
			eprintf("too few values for parameter: \"%s\".\n", line);
			exit(1);
		}
	}
	fclose(in);
	
	/* quit if some parameters are unset (Ato is a master programmer, 
	   but this would still be a *little* bit easier in python :)) */
#define CHECK_PARSED(A) \
	if (parsed.A == 0) { \
		eprintf("parameter unset: \"%s\".\n", #A); \
		allparsed = 0; \
	}
	
	CHECK_PARSED(MMIN);
	CHECK_PARSED(BINBIN);
	CHECK_PARSED(BINBIN_FEWBODY);
	CHECK_PARSED(BINSINGLE);
	CHECK_PARSED(BINSINGLE_FEWBODY);
	CHECK_PARSED(CENTRAL_MASS);
	CHECK_PARSED(DT_FACTOR);
	CHECK_PARSED(DUMPS);
	CHECK_PARSED(E_CONS);
	CHECK_PARSED(IDUM);
	CHECK_PARSED(INDEX_UNIT);
	CHECK_PARSED(INPUT_FILE);
	CHECK_PARSED(MASS_PC);
	CHECK_PARSED(MAX_INDEX);
	CHECK_PARSED(MEGA_YEAR);
	CHECK_PARSED(METALLICITY);
	CHECK_PARSED(MINIMUM_R);
	CHECK_PARSED(MIN_LAGRANGIAN_RADIUS);
	CHECK_PARSED(N_BINARY);
	CHECK_PARSED(NUM_CORE_STARS);
	CHECK_PARSED(ORIGINAL_PERTURB_STARS);
	CHECK_PARSED(PERTURB);
	CHECK_PARSED(R_MAX);
	CHECK_PARSED(SIN2BETA_MAX);
	CHECK_PARSED(SOLAR_MASS_DYN);
	CHECK_PARSED(STELLAR_EVOLUTION);
	CHECK_PARSED(TERMINAL_ENERGY_DISPLACEMENT);
	CHECK_PARSED(T_MAX);
	CHECK_PARSED(T_MAX_COUNT);
	CHECK_PARSED(T_PRINT_STEP);
	CHECK_PARSED(WIND_FACTOR);
	CHECK_PARSED(GAMMA);
	
#undef CHECK_PARSED

	/* exit if something is not set */
	if (!allparsed) {
		exit(1);
	}
	
	/* read the number of stars and possibly other parameters */
	read_fits_file_parameters(INPUT_FILE, r);
	
	if(!ReadSnapshot){
		clus.N_STAR_NEW = clus.N_STAR;
		/* add clus.N_STAR for collisions */
		N_STAR_DIM = clus.N_STAR + 2 + clus.N_BINARY + clus.N_STAR;

		/*********************************************/
		/* allocation of memory for global variables */
		/*********************************************/

		/* the main star array containing all star parameters */
		star = malloc(N_STAR_DIM * sizeof(star_t));
	
		/* make all stars initially uninteracted - this fixes an 
		   uninitialized memory read found by Valgrind */
		for (i=0; i<N_STAR_DIM; i++) {
			star[i].interacted = 0;
		}
		
		/* the main binary array containing all binary parameters */
		binary = malloc(N_STAR_DIM * sizeof(binary_t));

		/* This looks like the lookup table index for calculating pot */
		IndexTable = (long *) malloc((MAX_INDEX + 5) * sizeof(long));

		/* quantities calculated for various lagrange radii */
		mass_r = malloc((NUM_MASS_RADII_BINS + 1) 
					* MASS_PC_COUNT * sizeof(double));
		ave_mass_r = malloc((NUM_MASS_RADII_BINS + 1) 
					* MASS_PC_COUNT * sizeof(double));
		no_star_r = malloc((NUM_MASS_RADII_BINS + 1) 
					* MASS_PC_COUNT * sizeof(double));
		densities_r = malloc((NUM_MASS_RADII_BINS + 1) 
					* MASS_PC_COUNT * sizeof(double));
		mass_pc = malloc(MASS_PC_COUNT * sizeof(double));
		for(i=0;i<MASS_PC_COUNT; i++){
			mass_pc[i] = 0.0;
		}

		/*======= Reading of values for the Lagrange radii =======*/
		curr_mass = (char *) strtok(MASS_PC, ",; ");

		for (i = 0; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; i++)
			sscanf(curr_mass, "%lf", &mass_pc[i]);
	}

	/*======= Opening of output files =======*/
	sscanf("w", "%s", outfilemode);
	
	sprintf(outfile, "%s1", outprefix);
	if ((out[0] = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s3", outprefix);
	if ((out[1] = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s4", outprefix);
	if ((out[2] = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s6", outprefix);
	if ((ave_mass_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s7", outprefix);
	if ((no_star_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s8", outprefix);
	if ((densities_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s9", outprefix);
	if ((centmass_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}

	sprintf(outfile, "%s0", outprefix);
	if ((logfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create log output file \"%s\".\n", outfile);
		exit(1);
	}

	/* output files for binaries */
	if (clus.N_BINARY > 0) {
		/* general binary information */
		sprintf(outfile, "%sbinary", outprefix);
		if ((binaryfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create binary file \"%s\".\n", outfile);
			exit(1);
		}
		/* file for binary-single information */
		sprintf(outfile, "%sbinsinglelog.gz", outprefix);
		if ((binsinglefile = gzopen(outfile, "wb")) == NULL) {
			eprintf("cannot create binsinglelog file \"%s\".\n", outfile);
			exit(1);
		}
		/* file for fewbody information */
		sprintf(outfile, "%sbinbinlog.gz", outprefix);
		if ((binbinfile = gzopen(outfile, "wb")) == NULL) {
			eprintf("cannot create binbinlog file \"%s\".\n", outfile);
			exit(1);
		}
	}

	/* File for parameters of escaping stars */
	sprintf(outfile, "%sesc", outprefix);
	if ((escfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create escapers file \"%s\".\n", outfile);
		exit(1);
	}

	return(1);
}

/* close buffers */
void close_buffers(void)
{
	fclose(out[0]);
	fclose(out[1]);
	fclose(out[2]);
	fclose(logfile);
	fclose(ave_mass_file);
	fclose(no_star_file);
	fclose(densities_file);
	fclose(centmass_file);
	fclose(escfile);
	
	if (clus.N_BINARY > 0) {
		fclose(binaryfile);
		gzclose(binsinglefile);
		gzclose(binbinfile);
	}
}
