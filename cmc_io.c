/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include "cmc.h"
#include "cmc_vars.h"

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** %s %s (%s) [%s] **\n", CMCPRETTYNAME, CMCVERSION, CMCNICK, CMCDATE);
}

/* print the usage */
void cmc_print_usage(FILE *stream, char *argv[])
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  %s [options...] <input_file> <output_file_prefix>\n", argv[0]);
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -q --quiet   : do not print diagnostic info to stdout\n");
	fprintf(stream, "  -d --debug   : turn on debugging\n");
	fprintf(stream, "  -V --version : print version info\n");
	fprintf(stream, "  -h --help    : display this help text\n");
}

void print_results(void){
	PrintLogOutput();
	PrintFileOutput();
	fflush(NULL);
}

/*********** Output 2D/3D snapshots **************/
void print_2Dsnapshot(void)
{
	long i, j;
	char outfile[100];
	
	if (SNAPSHOTTING) {
		/* open file for 2D snapshot */
		sprintf(outfile, "%s.snap%04ld.dat.gz", outprefix, snap_num);
		if ((snapfile = (FILE *) gzopen(outfile, "wb")) == NULL) {
			eprintf("cannot create 2D snapshot file %s\n", outfile);
			exit_cleanly(1);
		}
		
		/* print useful header */
		gzprintf(snapfile, "# t=%.8g [code units]; All quantities below are in code units unless otherwise specified.\n", TotalTime);
		gzprintf(snapfile, "#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]\n");
		
		/* then print data */
		for (i=1; i<=clus.N_MAX; i++) {
			gzprintf(snapfile, "%ld %.8g %.8g %.8g %.8g %.8g %.8g ", 
				 star[i].id, star[i].m * (units.m / clus.N_STAR) / MSUN, 
				 star[i].r, star[i].vr, star[i].vt, 
				 star[i].E, star[i].J);
			if (star[i].binind) {
				j = star[i].binind;
				gzprintf(snapfile, "1 %.8g %.8g %ld %ld %.8g %.8g ", 
					 binary[j].m1 * (units.m / clus.N_STAR) / MSUN, 
					 binary[j].m2 * (units.m / clus.N_STAR) / MSUN, 
					 binary[j].id1, binary[j].id2,
					 binary[j].a * units.l / AU, binary[j].e);
			} else {
				gzprintf(snapfile, "0 0 0 0 0 0 0 ");	
			}
			
			if (star[i].binind == 0) {
				gzprintf(snapfile, "%d %.8g %.8g ", 
					 star[i].se_k, star[i].se_lum, star[i].rad * units.l / RSUN);
			} else {
				gzprintf(snapfile, "0 0 0 ");
			}
			gzprintf(snapfile, "\n");
		}
		
		gzclose(snapfile);
		
		/* global counter for snapshot output file */
		snap_num++;
	}
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

	gprintf("******************************************************************************\n");
	fprintf(logfile, "******************************************************************************\n");

	gprintf("tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);
	fprintf(logfile, "tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);

	gprintf("Etotal=%g max_r=%g N_bound=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);
	fprintf(logfile, "Etotal=%g max_r=%g N_bound=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);
	
	gprintf("Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);
	fprintf(logfile, "Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);
	
	gprintf("TidalMassLoss=%g\n", TidalMassLoss);
	fprintf(logfile, "TidalMassLoss=%g\n", TidalMassLoss);
	
	gprintf("core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	fprintf(logfile, "core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	
	gprintf("trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	fprintf(logfile, "trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	
	gprintf("N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);
	fprintf(logfile, "N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);

	gprintf("******************************************************************************\n");
	fprintf(logfile, "******************************************************************************\n");
}

void PrintFileOutput(void) {
	long i, j, n_single, n_binary;
	double fb, fb_core;
	int *multimassr_empty  = (int *) malloc((NO_MASS_BINS-1)*sizeof(int));

	/* print useful headers */
	if (tcount == 1) {
		fprintf(lagradfile, "# Lagrange radii [code units]\n");
		fprintf(ave_mass_file, "# Average mass within Lagrange radii [M_sun]\n");
		fprintf(no_star_file, "# Number of stars within Lagrange radii [dimensionless]\n");
		fprintf(densities_file, "# Density within Lagrange radii [code units]\n");
		fprintf(ke_rad_file, "# Total radial kinetic energy within Lagrange radii [code units]\n");
		fprintf(ke_tan_file, "# Total tangential kinetic energy within Lagrange radii [code units]\n");
		fprintf(v2_rad_file, "# Sum of v_r within Lagrange radii [code units]\n");
		fprintf(v2_tan_file, "# Sum of v_t within Lagrange radii [code units]\n");
		for(i=0; i<NO_MASS_BINS-1; i++){
			fprintf(mlagradfile[i], "# Lagrange radii for %g < m < %g range [code units]\n", mass_bins[i], mass_bins[i+1]);
		}
		
		fprintf(lagradfile, "# 1:t");
		fprintf(ave_mass_file, "# 1:t");
		fprintf(no_star_file, "# 1:t");
		fprintf(densities_file, "# 1:t");
		fprintf(ke_rad_file, "# 1:t");
		fprintf(ke_tan_file, "# 1:t");
		fprintf(v2_rad_file, "# 1:t");
		fprintf(v2_tan_file, "# 1:t");
		for(i=0; i<NO_MASS_BINS-1; i++){
			fprintf(mlagradfile[i], "# 1:t");
		}

		for (i=0; i<MASS_PC_COUNT; i++) {
			fprintf(lagradfile, " %ld:r(%g)", i+2, mass_pc[i]);
			fprintf(ave_mass_file, " %ld:<m>(%g)", i+2, mass_pc[i]);
			fprintf(no_star_file, " %ld:N(%g)", i+2, mass_pc[i]);
			fprintf(densities_file, " %ld:rho(%g)", i+2, mass_pc[i]);
			fprintf(ke_rad_file, " %ld:T_r(%g)", i+2, mass_pc[i]);
			fprintf(ke_tan_file, " %ld:T_t(%g)", i+2, mass_pc[i]);
			fprintf(v2_rad_file, " %ld:V2_r(%g)", i+2, mass_pc[i]);
			fprintf(v2_tan_file, " %ld:V2_t(%g)", i+2, mass_pc[i]);
			for(j=0; j<NO_MASS_BINS-1; j++){
				fprintf(mlagradfile[j], " %ld:r(%g)", i+2, mass_pc[i]);
			}
		}

		fprintf(lagradfile, "\n");
		fprintf(ave_mass_file, "\n");
		fprintf(no_star_file, "\n");
		fprintf(densities_file, "\n");
		fprintf(ke_rad_file, "\n");
		fprintf(ke_tan_file, "\n");
		fprintf(v2_rad_file, "\n");
		fprintf(v2_tan_file, "\n");
		for(i=0; i<NO_MASS_BINS-1; i++){
			fprintf(mlagradfile[i], "\n");
		}
	}

	/* print data */
	fprintf(lagradfile, "%.9e ", TotalTime);
	fprintf(ave_mass_file, "%.9e ",TotalTime);
	fprintf(no_star_file, "%.9e ",TotalTime);
	fprintf(densities_file, "%.9e ",TotalTime);
	fprintf(ke_rad_file, "%.9e ",TotalTime);
	fprintf(ke_tan_file, "%.9e ",TotalTime);
	fprintf(v2_rad_file, "%.9e ",TotalTime);
	fprintf(v2_tan_file, "%.9e ",TotalTime);
	for(i=0; i<NO_MASS_BINS-1; i++){
		multimassr_empty[i] = 1;
		for(j=0; j<NO_MASS_BINS-1; j++){
			if (multi_mass_r[j][i] > 0.0) multimassr_empty[i] = 0;
		}
	}
	for(i=0; i<NO_MASS_BINS-1; i++){
		if ( !multimassr_empty[i] ){
			fprintf(mlagradfile[i], "%.9e ", TotalTime);
		}
	}

	for (i = 0; i < MASS_PC_COUNT ; i++) {
		fprintf(lagradfile, "%e ", mass_r[i]);
		fprintf(ave_mass_file,"%e ", ave_mass_r[i] * units.m / MSUN);
		fprintf(no_star_file,"%g ", no_star_r[i]);
		fprintf(densities_file,"%e ", densities_r[i]);
		fprintf(ke_rad_file,"%e ", ke_rad_r[i]);
		fprintf(ke_tan_file,"%e ", ke_tan_r[i]);
		fprintf(v2_rad_file,"%e ", v2_rad_r[i]);
		fprintf(v2_tan_file,"%e ", v2_tan_r[i]);
		for(j=0; j<NO_MASS_BINS-1; j++){
			if ( !multimassr_empty[j] ){
				fprintf(mlagradfile[j], "%e ", multi_mass_r[j][i]);
			}
		}
	}
	fprintf(lagradfile, "\n");
	fprintf(ave_mass_file,"\n");
	fprintf(no_star_file,"\n");
	fprintf(densities_file,"\n");
	fprintf(ke_rad_file,"\n");
	fprintf(ke_tan_file,"\n");
	fprintf(v2_rad_file,"\n");
	fprintf(v2_tan_file,"\n");
	for(i=0; i<NO_MASS_BINS-1; i++){
		if ( !multimassr_empty[i] ){
			fprintf(mlagradfile[i], "\n");
		}
	}

	/* information on the central BH */
	/* print useful header */
	if (tcount == 1) {
		fprintf(centmass_file, "# Information on central black hole [code units unless otherwise noted]\n");
		fprintf(centmass_file, "#1:t #2:cenma.m #3:Dt #4:rho_core #5:Etotal.tot #6:Etotal.K #7:Etotal.P\n");
	}
	fprintf(centmass_file, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
		TotalTime, cenma.m * madhoc, Dt, rho_core, Etotal.tot, Etotal.K, Etotal.P);
	
	/* output Time,N_MAX,TotalE,TotalKE,TotalPE,Mtotal */
	/* print useful header */
	if (tcount == 1) {
		fprintf(dynfile, "# Dynamical information [code units]\n");
		fprintf(dynfile, "#1:t #2:Dt #3:tcount #4:N #5:M #6:VR #7:N_c #8:r_c #9:r_max #10:Etot #11:KE #12:PE #13:Etot_int #14:Etot_bin #15:E_cenma #16:Eesc #17:Ebesc #18:Eintesc #19:Eoops #20:Etot+Eoops #21:r_h #22:rho_0 #23:rc_spitzer #24:v0_rms\n");
	}
	fprintf(dynfile, "%.8g %.8g %ld %ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
		TotalTime, Dt, tcount, clus.N_MAX, Mtotal, -2.0*Etotal.K/Etotal.P, N_core, core_radius, max_r, 
		Etotal.tot, Etotal.K, Etotal.P, Etotal.Eint, Etotal.Eb, cenma.E, Eescaped, Ebescaped, Eintescaped, 
		Eoops, Etotal.tot+Eoops, clusdyn.rh, central.rho, central.rc_spitzer, central.v_rms);
	
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
		
		/* print useful header */
		if (tcount == 1) {
			fprintf(binaryfile, "# Binary information [code units]\n");
			fprintf(binaryfile, "# 1:t 2:N_b 3:M_b 4:E_b 5:r_h,s 6:r_h,b 7:rho_c,s 8:rho_c,b 9:N_bb 10:N_bs 11:f_b,c 12:f_b 13:E_bb 14:E_bs 15:DE_bb 16:DE_bs\n");
		}
		/* print data */
		fprintf(binaryfile,
			"%.6g %ld %.6g %.6g %.6g %.6g %.6g %.6g %ld %ld %.6g %.6g %.6g %.6g %.6g %.6g\n",
			TotalTime, N_b, M_b, E_b, rh_single, 
			rh_binary, rho_core_single, rho_core_bin, 
			N_bb, N_bs, fb_core, 
			fb, E_bb, E_bs, DE_bb, DE_bs);
	}

	/* also saves INITIAL snapshot (StepCount=0) */
	if (TotalTime >= SNAPSHOT_DELTAT * StepCount) {
		StepCount++;
		print_2Dsnapshot();
	}
}

/*** Parsing of Input Parameters / Memory allocation / File I/O ***/
void parser(int argc, char *argv[], gsl_rng *r)
{
	char inputfile[1024], outfile[1024], outfilemode[5];
	char parameter_name[1024], values[1024], dummy[1024], line[2048];
	char *curr_mass;
	parsed_t parsed;
	parsed_t *spp;
	long i;
	int allparsed=1;
	/* int *ip; */
	FILE *in, *parsedfp;
	const char *short_opts = "qdVh";
	const struct option long_opts[] = {
		{"quiet", no_argument, NULL, 'q'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* DEFAULT PARAMETER VALUES */
	quiet = 0;
	debug = 0;
	NO_MASS_BINS = 0;
	/* DEFAULT PARAMETER VALUES */
	
	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'q':
			quiet = 1;
			break;
		case 'd':
			debug = 1;
			break;
		case 'V':
			print_version(stdout);
			exit(0);
		case 'h':
			print_version(stdout);
			fprintf(stdout, "\n");
			cmc_print_usage(stdout, argv);
			exit(0);
		default:
			break;
		}
	}

	/* check to make sure there was nothing crazy on the command line */
	if (argc - optind != 2) {
		cmc_print_usage(stdout, argv);
		exit(0);
	}
	
	/* set inputfile and outprefix now that the options have been parsed */
	sprintf(inputfile, "%s", argv[optind]);
	sprintf(outprefix, "%s", argv[optind+1]);
	dprintf("inputfile=%s outprefix=%s\n", inputfile, outprefix);

	/*======= Opening of input & output files =======*/
	if ((in = fopen(inputfile, "r")) == NULL) {
		eprintf("Cannot open input file \"%s\".\n", inputfile);
		exit(1);
	}
	
	sprintf(outfile, "%s.cmc.parsed", outprefix);
	if ((parsedfp = fopen(outfile, "w")) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}

	/* nothing is set yet, so assigning all zeros to variable parsed */
	/* one way to do it is a complicated for loop (?which depends on the 
	 * fact that all members of the struct is int?), this is the fun way */
	//ip = (int *) &parsed;
	//for(i=0; i<(sizeof(parsed)/sizeof(int)); i++){
	//	*(ip++) = 0;
	//}
	/* faster and portable way to do this is using calloc, this does
	 * not depend on the structure of the struct */
	spp = (parsed_t *) calloc(1, sizeof(parsed_t));
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

#define PRINT_PARSED(DOC) fprintf(parsedfp, "# %s\n%s\n", DOC, line)

		/* see if there are too many values for parameter */
		if (sscanf(line, "%s %s %s", parameter_name, values, dummy) == 3) {
			eprintf("too many values for parameter: \"%s\".\n", line);
			exit(1);
		} else if (sscanf(line, "%s %s", parameter_name, values) == 2) {
			if (strcmp(parameter_name, "BINBIN") == 0) {
				PRINT_PARSED(PARAMDOC_BINBIN);
				sscanf(values, "%d", &BINBIN);
				parsed.BINBIN = 1;
			} else if (strcmp(parameter_name, "BINSINGLE") == 0) {
				PRINT_PARSED(PARAMDOC_BINSINGLE);
				sscanf(values, "%d", &BINSINGLE);
				parsed.BINSINGLE = 1;
			} else if (strcmp(parameter_name, "SNAPSHOTTING") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOTTING);
				sscanf(values, "%d", &SNAPSHOTTING);
				parsed.SNAPSHOTTING = 1;
                        } else if (strcmp(parameter_name, "SNAPSHOT_DELTACOUNT") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_DELTACOUNT);
				sscanf(values, "%ld", &SNAPSHOT_DELTACOUNT);
				parsed.SNAPSHOT_DELTACOUNT = 1;
                        } else if (strcmp(parameter_name, "SNAPSHOT_DELTAT") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_DELTAT);
				sscanf(values, "%lf", &SNAPSHOT_DELTAT);
				parsed.SNAPSHOT_DELTAT = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_CORE_BOUNCE") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_CORE_BOUNCE);
				sscanf(values, "%d", &SNAPSHOT_CORE_BOUNCE);
				parsed.SNAPSHOT_CORE_BOUNCE = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_CORE_COLLAPSE") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_CORE_COLLAPSE);
				sscanf(values, "%d", &SNAPSHOT_CORE_COLLAPSE);
				parsed.SNAPSHOT_CORE_COLLAPSE = 1;
			} else if (strcmp(parameter_name, "IDUM") == 0) {
				PRINT_PARSED(PARAMDOC_IDUM);
				sscanf(values, "%ld", &IDUM);
				parsed.IDUM = 1;
			} else if (strcmp(parameter_name, "INPUT_FILE") == 0) {
				PRINT_PARSED(PARAMDOC_INPUT_FILE);
				sscanf(values, "%s", INPUT_FILE);
				parsed.INPUT_FILE = 1;
			} else if (strcmp(parameter_name, "MASS_PC") == 0) {
				PRINT_PARSED(PARAMDOC_MASS_PC);
				strcpy(MASS_PC, values);
				curr_mass = (char *) strtok(values, ",; ");
				for (MASS_PC_COUNT = 1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; MASS_PC_COUNT++);
				parsed.MASS_PC = 1;
			} else if (strcmp(parameter_name, "MASS_BINS") == 0) {
				PRINT_PARSED(PARAMDOC_MASS_BINS);
				/* we recycle variable "curr_mass" for mass bins */
				strcpy(MASS_BINS, values);
				curr_mass = (char *) strtok(values, ",; ");
				for (NO_MASS_BINS = 1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; NO_MASS_BINS++);
				parsed.MASS_BINS = 1;
			} else if (strcmp(parameter_name, "MINIMUM_R") == 0) {
				PRINT_PARSED(PARAMDOC_MINIMUM_R);
				sscanf(values, "%lf", &MINIMUM_R);
				parsed.MINIMUM_R = 1;
			} else if (strcmp(parameter_name, "STOPATCORECOLLAPSE") == 0) {
				PRINT_PARSED(PARAMDOC_STOPATCORECOLLAPSE);
				sscanf(values, "%d", &STOPATCORECOLLAPSE);
				parsed.STOPATCORECOLLAPSE = 1;
			} else if (strcmp(parameter_name, "NUM_CENTRAL_STARS") == 0) {
				PRINT_PARSED(PARAMDOC_NUM_CENTRAL_STARS);
				sscanf(values, "%ld", &NUM_CENTRAL_STARS);
				parsed.NUM_CENTRAL_STARS = 1;
			} else if (strcmp(parameter_name, "PERTURB") == 0) {
				PRINT_PARSED(PARAMDOC_PERTURB);
				sscanf(values, "%ld", &PERTURB);
				parsed.PERTURB = 1;
			} else if (strcmp(parameter_name, "RELAXATION") == 0) {
				PRINT_PARSED(PARAMDOC_RELAXATION);
				sscanf(values, "%ld", &RELAXATION);
				parsed.RELAXATION = 1;
			} else if (strcmp(parameter_name, "THETASEMAX") == 0) {
				PRINT_PARSED(PARAMDOC_THETASEMAX);
				sscanf(values, "%lf", &THETASEMAX);
				parsed.THETASEMAX = 1;
			} else if (strcmp(parameter_name, "STELLAR_EVOLUTION") == 0) {
				PRINT_PARSED(PARAMDOC_STELLAR_EVOLUTION);
				sscanf(values, "%ld", &STELLAR_EVOLUTION);
				parsed.STELLAR_EVOLUTION = 1;
			} else if (strcmp(parameter_name, "SS_COLLISION") == 0) {
				PRINT_PARSED(PARAMDOC_SS_COLLISION);
				sscanf(values, "%ld", &SS_COLLISION);
				parsed.SS_COLLISION = 1;
			} else if (strcmp(parameter_name, "TERMINAL_ENERGY_DISPLACEMENT") == 0) {
				PRINT_PARSED(PARAMDOC_TERMINAL_ENERGY_DISPLACEMENT);
				sscanf(values, "%lf", &TERMINAL_ENERGY_DISPLACEMENT);
				parsed.TERMINAL_ENERGY_DISPLACEMENT = 1;
			} else if (strcmp(parameter_name, "T_MAX") == 0) {
				PRINT_PARSED(PARAMDOC_T_MAX);
				sscanf(values, "%lf", &T_MAX);
				parsed.T_MAX = 1;
			} else if (strcmp(parameter_name, "T_MAX_COUNT") == 0) {
				PRINT_PARSED(PARAMDOC_T_MAX_COUNT);
				sscanf(values, "%ld", &T_MAX_COUNT);
				parsed.T_MAX_COUNT = 1;
			} else if (strcmp(parameter_name, "MAX_WCLOCK_TIME") == 0) {
				PRINT_PARSED(PARAMDOC_MAX_WCLOCK_TIME);
				sscanf(values, "%ld", &MAX_WCLOCK_TIME);
				parsed.MAX_WCLOCK_TIME = 1;
			} else if (strcmp(parameter_name, "WIND_FACTOR") == 0) {
				PRINT_PARSED(PARAMDOC_WIND_FACTOR);
				sscanf(values, "%lf", &WIND_FACTOR);
				parsed.WIND_FACTOR = 1;
			} else if (strcmp(parameter_name, "GAMMA") == 0) {
				PRINT_PARSED(PARAMDOC_GAMMA);
				sscanf(values, "%lf", &GAMMA);
				parsed.GAMMA = 1;
			} else if (strcmp(parameter_name, "SEARCH_GRID")== 0) {
				PRINT_PARSED(PARAMDOC_SEARCH_GRID);
				sscanf(values, "%ld", &SEARCH_GRID);
				parsed.SEARCH_GRID = 1;
                        } else if (strcmp(parameter_name, "SG_STARSPERBIN")== 0) {
				PRINT_PARSED(PARAMDOC_SG_STARSPERBIN);
				sscanf(values, "%ld", &SG_STARSPERBIN);
				parsed.SG_STARSPERBIN = 1;
                        } else if (strcmp(parameter_name, "SG_MAXLENGTH")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MAXLENGTH);
				sscanf(values, "%ld", &SG_MAXLENGTH);
				parsed.SG_MAXLENGTH = 1;
                        } else if (strcmp(parameter_name, "SG_MINLENGTH")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MINLENGTH);
				sscanf(values, "%ld", &SG_MINLENGTH);
				parsed.SG_MINLENGTH = 1;
                        } else if (strcmp(parameter_name, "SG_POWER_LAW_EXPONENT")== 0) {
				PRINT_PARSED(PARAMDOC_SG_POWER_LAW_EXPONENT);
				sscanf(values, "%lf", &SG_POWER_LAW_EXPONENT);
				parsed.SG_POWER_LAW_EXPONENT = 1;
                        } else if (strcmp(parameter_name, "SG_MATCH_AT_FRACTION")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MATCH_AT_FRACTION);
				sscanf(values, "%lf", &SG_MATCH_AT_FRACTION);
				parsed.SG_MATCH_AT_FRACTION = 1;
                        } else if (strcmp(parameter_name, "SG_PARTICLE_FRACTION")== 0) {
				PRINT_PARSED(PARAMDOC_SG_PARTICLE_FRACTION);
				sscanf(values, "%lf", &SG_PARTICLE_FRACTION);
				parsed.SG_PARTICLE_FRACTION = 1;
			} else if (strcmp(parameter_name, "BH_LOSS_CONE")== 0) {
				PRINT_PARSED(PARAMDOC_BH_LOSS_CONE);
				sscanf(values, "%li", &BH_LOSS_CONE);
				parsed.BH_LOSS_CONE = 1;
                        } else if (strcmp(parameter_name, "BH_R_DISRUPT_NB")== 0) {
				PRINT_PARSED(PARAMDOC_BH_R_DISRUPT_NB);
				sscanf(values, "%lf", &BH_R_DISRUPT_NB);
				parsed.BH_R_DISRUPT_NB = 1;
        		} else if (strcmp(parameter_name, "FORCE_RLX_STEP")== 0) {
				PRINT_PARSED(PARAMDOC_FORCE_RLX_STEP);
				sscanf(values, "%i", &FORCE_RLX_STEP);
				parsed.FORCE_RLX_STEP = 1;
#ifdef EXPERIMENTAL
                        } else if (strcmp(parameter_name, "BH_LC_FDT")== 0) {
				PRINT_PARSED(PARAMDOC_BH_LC_FDT);
				sscanf(values, "%lf", &BH_LC_FDT);
				parsed.BH_LC_FDT = 1;
                        } else if (strcmp(parameter_name, "AVEKERNEL")== 0) {
				PRINT_PARSED(PARAMDOC_AVEKERNEL);
				sscanf(values, "%li", &AVEKERNEL);
				parsed.AVEKERNEL = 1;
#endif
			} else if (strcmp(parameter_name, "APSIDES_PRECISION")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_PRECISION);
				sscanf(values, "%lf", &APSIDES_PRECISION);
				parsed.APSIDES_PRECISION = 1;
        		} else if (strcmp(parameter_name, "APSIDES_MAX_ITER")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_MAX_ITER);
				sscanf(values, "%li", &APSIDES_MAX_ITER);
				parsed.APSIDES_MAX_ITER = 1;
                        } else if (strcmp(parameter_name, "APSIDES_CONVERGENCE")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_CONVERGENCE);
				sscanf(values, "%lf", &APSIDES_CONVERGENCE);
				parsed.APSIDES_CONVERGENCE = 1;
                        } else if (strcmp(parameter_name, "CIRC_PERIOD_THRESHOLD")== 0) {
				PRINT_PARSED(PARAMDOC_CIRC_PERIOD_THRESHOLD);
				sscanf(values, "%lf", &CIRC_PERIOD_THRESHOLD);
				parsed.CIRC_PERIOD_THRESHOLD = 1;
                        } else if (strcmp(parameter_name, "WRITE_STELLAR_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_STELLAR_INFO);
				sscanf(values, "%i", &WRITE_STELLAR_INFO);
				parsed.WRITE_STELLAR_INFO = 1;
			} else {
				wprintf("unknown parameter: \"%s\".\n", line);
			}
		} else if (sscanf(line, "%s", parameter_name) == 1) {
			eprintf("too few values for parameter: \"%s\".\n", line);
			exit(1);
		}
	}
	fclose(in);
	
	/* quit if some parameters are unset */
#define CHECK_PARSED(A) \
	if (parsed.A == 0) { \
		eprintf("parameter \"%s\" unset.\n", #A); \
		allparsed = 0; \
	}
	
	CHECK_PARSED(GAMMA);
	CHECK_PARSED(INPUT_FILE);
	CHECK_PARSED(MASS_PC);
	CHECK_PARSED(MASS_BINS);
	
#undef CHECK_PARSED

/* but only warn if some other parameters are unset and default values are used */
#define CHECK_PARSED(A,DEFAULT,DOC) \
	if (parsed.A == 0) { \
		wprintf("parameter \"%s\" unset: using default value \"%s\").\n", #A, #DEFAULT); \
                A=DEFAULT; \
                fprintf(parsedfp, "# %s\n%s %s     # default value\n", DOC, #A, #DEFAULT); \
	}
	
	CHECK_PARSED(PERTURB, 1, PARAMDOC_PERTURB);
	CHECK_PARSED(RELAXATION, 1, PARAMDOC_RELAXATION);
	CHECK_PARSED(THETASEMAX, 1.0, PARAMDOC_THETASEMAX);
	CHECK_PARSED(STELLAR_EVOLUTION, 0, PARAMDOC_STELLAR_EVOLUTION);
        CHECK_PARSED(WRITE_STELLAR_INFO, 1, PARAMDOC_WRITE_STELLAR_INFO);
	CHECK_PARSED(WIND_FACTOR, 1.0, PARAMDOC_WIND_FACTOR);
	CHECK_PARSED(SS_COLLISION, 0, PARAMDOC_SS_COLLISION);
	CHECK_PARSED(BINBIN, 1, PARAMDOC_BINBIN);
	CHECK_PARSED(BINSINGLE, 1, PARAMDOC_BINSINGLE);
	CHECK_PARSED(BH_LOSS_CONE, 0, PARAMDOC_BH_LOSS_CONE);
	CHECK_PARSED(MINIMUM_R, 0.0, PARAMDOC_MINIMUM_R);
	CHECK_PARSED(BH_R_DISRUPT_NB, 0., PARAMDOC_BH_R_DISRUPT_NB);
        CHECK_PARSED(CIRC_PERIOD_THRESHOLD, 1e-18, PARAMDOC_CIRC_PERIOD_THRESHOLD);
	
        CHECK_PARSED(T_MAX, 20.0, PARAMDOC_T_MAX);
	CHECK_PARSED(T_MAX_COUNT, 1000000, PARAMDOC_T_MAX_COUNT);
	CHECK_PARSED(MAX_WCLOCK_TIME, 2592000, PARAMDOC_MAX_WCLOCK_TIME);
	CHECK_PARSED(STOPATCORECOLLAPSE, 1, PARAMDOC_STOPATCORECOLLAPSE);
	CHECK_PARSED(TERMINAL_ENERGY_DISPLACEMENT, 0.5, PARAMDOC_TERMINAL_ENERGY_DISPLACEMENT);

	CHECK_PARSED(SNAPSHOTTING, 0, PARAMDOC_SNAPSHOTTING);
	CHECK_PARSED(SNAPSHOT_DELTACOUNT, 250, PARAMDOC_SNAPSHOT_DELTACOUNT);
	CHECK_PARSED(SNAPSHOT_DELTAT, 0.25, PARAMDOC_SNAPSHOT_DELTAT);
	CHECK_PARSED(SNAPSHOT_CORE_COLLAPSE, 0, PARAMDOC_SNAPSHOT_CORE_COLLAPSE);
        CHECK_PARSED(SNAPSHOT_CORE_BOUNCE, 0, PARAMDOC_SNAPSHOT_CORE_BOUNCE);

	CHECK_PARSED(NUM_CENTRAL_STARS, 300, PARAMDOC_NUM_CENTRAL_STARS);
	CHECK_PARSED(IDUM, 0, PARAMDOC_IDUM);

	CHECK_PARSED(SEARCH_GRID, 0, PARAMDOC_SEARCH_GRID);
        CHECK_PARSED(SG_STARSPERBIN, 100, PARAMDOC_SG_STARSPERBIN);
        CHECK_PARSED(SG_MAXLENGTH, 1000000, PARAMDOC_SG_MAXLENGTH);
        CHECK_PARSED(SG_MINLENGTH, 1000, PARAMDOC_SG_MINLENGTH);
        CHECK_PARSED(SG_POWER_LAW_EXPONENT, 0.5, PARAMDOC_SG_POWER_LAW_EXPONENT);
        CHECK_PARSED(SG_MATCH_AT_FRACTION, 0.5, PARAMDOC_SG_MATCH_AT_FRACTION);
        CHECK_PARSED(SG_PARTICLE_FRACTION, 0.95, PARAMDOC_SG_PARTICLE_FRACTION);
        CHECK_PARSED(FORCE_RLX_STEP, 0, PARAMDOC_FORCE_RLX_STEP);
#ifdef EXPERIMENTAL
        CHECK_PARSED(BH_LC_FDT, 0.0, PARAMDOC_BH_LC_FDT);
        CHECK_PARSED(AVEKERNEL, 20, PARAMDOC_AVEKERNEL);
#endif
        CHECK_PARSED(APSIDES_PRECISION, 1.0e-11, PARAMDOC_APSIDES_PRECISION);
        CHECK_PARSED(APSIDES_MAX_ITER, 100, PARAMDOC_APSIDES_MAX_ITER);
        CHECK_PARSED(APSIDES_CONVERGENCE, 5.e-13, PARAMDOC_APSIDES_CONVERGENCE);
#undef CHECK_PARSED

	/* exit if something is not set */
	if (!allparsed) {
		exit(1);
	}
	
	fclose(parsedfp);
	
	/* read the number of stars and possibly other parameters */
	cmc_read_fits_file(INPUT_FILE, &cfd);
	clus.N_STAR = cfd.NOBJ;
	clus.N_BINARY = cfd.NBINARY;
	R_MAX = cfd.Rtid;
	METALLICITY = cfd.Z;

	clus.N_STAR_NEW = clus.N_STAR;
	/* add 2 * clus.N_BINARY for binary disruptions */
	N_STAR_DIM = 2 + clus.N_STAR + 2 * clus.N_BINARY;
	N_BIN_DIM = 2 + clus.N_BINARY;
	
	/* safety factors, so we don't have to worry about memory management/garbage collection */
	N_STAR_DIM = (long) floor(1.5 * ((double) N_STAR_DIM));
	N_BIN_DIM = (long) floor(1.5 * ((double) N_BIN_DIM));
	
	/*********************************************/
	/* allocation of memory for global variables */
	/*********************************************/
	
	/* the main star array containing all star parameters */
	star = (star_t *) calloc(N_STAR_DIM, sizeof(star_t));
	
	/* allocate memory for velocity dispersion array */
	sigma_array.n = 0;
	sigma_array.r = (double *) calloc(N_STAR_DIM, sizeof(double));
	sigma_array.sigma = (double *) calloc(N_STAR_DIM, sizeof(double));
	
	/* the main binary array containing all binary parameters */
	binary = (binary_t *) calloc(N_BIN_DIM, sizeof(binary_t));
	
	/* quantities calculated for various lagrange radii */
	mass_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ave_mass_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	no_star_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	densities_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ke_rad_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ke_tan_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	v2_rad_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	v2_tan_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	mass_pc = (double *) calloc(MASS_PC_COUNT, sizeof(double));
	mass_bins = (double *) calloc(NO_MASS_BINS, sizeof(double));
	multi_mass_r = (double **) malloc(NO_MASS_BINS * sizeof(double *));
	for(i=0; i<NO_MASS_BINS; i++){
		multi_mass_r[i] = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	}
	
	/*======= Reading values for the Lagrange radii =======*/
	curr_mass = (char *) strtok(MASS_PC, ",; ");
	sscanf(curr_mass, "%lf", &mass_pc[0]);
	
	for (i=1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; i++){
		sscanf(curr_mass, "%lf", &mass_pc[i]);
	}
	
	/*======= Reading values for the mass bins =======*/
	curr_mass = (char *) strtok(MASS_BINS, ",; ");
	sscanf(curr_mass, "%lf", &mass_bins[0]);

	for (i=1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; i++){
		sscanf(curr_mass, "%lf", &mass_bins[i]);
	}

	/*======= Opening of output files =======*/
	sscanf("w", "%s", outfilemode);
	
	sprintf(outfile, "%s.lagrad.dat", outprefix);
	if ((lagradfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.dyn.dat", outprefix);
	if ((dynfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.avemass_lagrad.dat", outprefix);
	if ((ave_mass_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.nostar_lagrad.dat", outprefix);
	if ((no_star_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.rho_lagrad.dat", outprefix);
	if ((densities_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.ke_rad_lagrad.dat", outprefix);
	if ((ke_rad_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.ke_tan_lagrad.dat", outprefix);
	if ((ke_tan_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.v2_rad_lagrad.dat", outprefix);
	if ((v2_rad_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.v2_tan_lagrad.dat", outprefix);
	if ((v2_tan_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
	sprintf(outfile, "%s.centmass.dat", outprefix);
	if ((centmass_file = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}

	sprintf(outfile, "%s.log", outprefix);
	if ((logfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create log output file \"%s\".\n", outfile);
		exit(1);
	}

	/* output files for binaries */
	if (clus.N_BINARY > 0) {
		/* general binary information */
		sprintf(outfile, "%s.bin.dat", outprefix);
		if ((binaryfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create binary file \"%s\".\n", outfile);
			exit(1);
		}
		/* file for binary interaction information */
		sprintf(outfile, "%s.binint.log", outprefix);
		if ((binintfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create binintlog file \"%s\".\n", outfile);
			exit(1);
		}
	}

	/* File for parameters of escaping stars */
	sprintf(outfile, "%s.esc.dat", outprefix);
	if ((escfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create escapers file \"%s\".\n", outfile);
		exit(1);
	}
	/* print header */
	fprintf(escfile, "#1:tcount #2:t #3:m #4:r #5:vr #6:vt #7:r_peri #8:r_apo #9:Rtidal #10:phi_rtidal #11:phi_zero #12:E #13:J\n");

	/* Collision log file */
	sprintf(outfile, "%s.collision.log", outprefix);
	if ((collisionfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create collision log file \"%s\".\n", outfile);
		exit(1);
	}
	/* print header */
	fprintf(collisionfile, "# time interaction_type id_merger(mass_merger) id1(m1):id2(m2):id3(m3):... (r)\n");

	/* Relaxation data file */
	sprintf(outfile, "%s.relaxation.dat", outprefix);
	if ((relaxationfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create relaxation data file \"%s\".\n", outfile);
		exit(1);
	}

	/* lagrange radii for multiple mass bins */
	mlagradfile = (FILE **) malloc((NO_MASS_BINS-1)*sizeof(FILE *));
	for(i=0; i<NO_MASS_BINS-1; i++){ /* NO_MASS_BINS need to be >=2 to be
							meaningful */
		sprintf(outfile, "%s.lagrad%ld-%g-%g.dat", outprefix, i,
							mass_bins[i], mass_bins[i+1]);
		if ((mlagradfile[i] = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create output file \"%s\".\n", outfile);
			exit(1);
		}
	}
}

/* close buffers */
void close_buffers(void)
{
	int i;

	fclose(lagradfile);
	fclose(dynfile);
	fclose(logfile);
	fclose(ave_mass_file);
	fclose(no_star_file);
	fclose(densities_file);
	fclose(ke_rad_file);
	fclose(ke_tan_file);
	fclose(v2_rad_file);
	fclose(v2_tan_file);
	fclose(centmass_file);
	fclose(escfile);
	fclose(collisionfile);
	fclose(relaxationfile);

	if (clus.N_BINARY > 0) {
		fclose(binaryfile);
		fclose(binintfile);
	}

	for(i=0; i<NO_MASS_BINS-1; i++){
		fclose(mlagradfile[i]);
	}
}

/* trap signals */
void trap_sigs(void)
{
	/* Catch some signals */
	signal(SIGINT, exit_cleanly);
	signal(SIGTERM, exit_cleanly);
	signal(SIGQUIT, exit_cleanly);
	signal(SIGUSR1, toggle_debugging);
	
	/* override GSL error handler */
	//gsl_set_error_handler(&sf_gsl_errhandler);
}
