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
#define NBIN 0L
#define LIMITS 0
#define EBMIN 3.0
#define EBMAX 133.0
#define SEED 8732UL
#define AVEKERNEL 20
#define BINMF 0
#define PEAK_A 30.0  //in AU
#define MIN_A 0.01  //in AU
#define MAX_A 100.0  //in AU
#define MAX_A_HS 1.0 //max SMA in units of HS boundary (1 is default) 
#define NPLANETS 0
#define HIGHMASSFB 0.
#define MCRIT 15. //in Msun
#define HIGHMASSAINDEX -1.75  //default taken from Sana et al. 2012
#define HIGHMASSMRATIOLOW -1. //default is the lowest stellar mass in the IMF

/* global variables needed by Star Track */
double *zpars, METALLICITY, WIND_FACTOR=1.0;
/* global variable for new star id */
long newstarid;

/* return roche lobe radius (as fraction of semimajor axis) for object 1 */
double roche(double m1, double m2)
{
  /* Paczynski (1971) [1971ARA&A...9..183P] */
  if (m1/m2 > 0.523) {
    return(0.38 + 0.2 * log10(m1/m2));
  } else {
    return(2.0/pow(3.0, 4.0/3.0) * pow(m1/(m1+m2), 1.0/3.0));
  }
}

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
	fprintf(stream, "  -l --limits <limits algorithm> : algorithm for setting limits on binary semimajor axes (0=physical, 1=kT prescription, 2=M67 model of Hurley, et al. (2005), 3=Ivanova, et al. (2005)), 4=Eq. 16 Egleton, Fitchett & Tout 1989 with peak at 30AU [%d]\n", LIMITS);
	fprintf(stream, "  -p --peak_a <peak in a dist.>           : choosing where the peak in lognormal a dist. is; need to choose for -l 4: [%g] AU\n", PEAK_A);
	fprintf(stream, "  -a --min_a <minimum a in a dist.>           : choosing amin in lognormal dist.; need to choose for -l 4: [%g] AU\n", MIN_A);
	fprintf(stream, "  -A --max_a <max a in a dist.>           : choosing amax in lognormal dist.; need to choose for -l 4: [%g] AU\n", MAX_A);
	fprintf(stream, "  -H --max_a_hs <select max_a as multiple (H) of hard-soft boundary>           : amax is [%g] times HS boundary; use with -l 0 or -l 5\n", MAX_A_HS);
	fprintf(stream, "  -m --Ebmin <E_b,min>           : minimum binding energy, in kT [%g]\n", EBMIN);
	fprintf(stream, "  -M --Ebmax <E_b,max>           : maximum binding energy, in kT [%g]\n", EBMAX);
	fprintf(stream, "  -I --ignoreradii               : ignore radii when setting binary properties\n");
	fprintf(stream, "  -s --seed <seed>               : random seed [%ld]\n", SEED);
	fprintf(stream, "  -b --binary_mf <choosing bin MF> : set a different MF (KTG91 for now) for binaries [%d]\n", BINMF);
	fprintf(stream, "  -f --highmass_fb <choosing binary fraction for high-mass stars separately> : set a different binary fraction for M>mcrit [%g]\nNeeds this for binmf=98\n", HIGHMASSFB);
	fprintf(stream, "  -G --mcrit <choosing mcrit for binary fraction for high-mass stars separately> [%g]\nNeeds this for binmf=98\n", MCRIT);
	fprintf(stream, "  -q --highmass_mratio_low <choosing the lower limit in mass ratios for high-mass stars>: \nallows to change the lower limit in mass ratio dist. for high-mass stars ; if >1 stars with mass>mcrit are coupled among themselves [%g]. \nAt the moment required only for binmf 98\n", HIGHMASSMRATIOLOW);
	fprintf(stream, "  -k --highmass_aind <choosing power-law a-dist for high-mass stars> : set a power-law a-dist. in dn/da [%g].\nMust have l=5\n", HIGHMASSAINDEX);	
	fprintf(stream, "  -h --help                      : display this help text\n");
}

long star_get_id_new(void)
{
	newstarid++;
	return(newstarid);
}

/* calculate and store the velocity dispersion profile */
void addbin_calc_sigma_r(cmc_fits_data_t *cfd, double *r, double *sigma, double *mave)
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
        fprintf(stderr,"%g %g\n",Sigma,Sigma*(PARSEC*cfd->Rvir/1e5)/(pow(cfd->Rvir * PARSEC, 1.5) / sqrt(G * cfd->Mclus * MSUN)));
		
		/* store sigma */
		r[si] = cfd->obj_r[si];
		sigma[si] = Sigma;
		mave[si] = Mave;
	}
}

void assign_binaries(cmc_fits_data_t *cfd, long Nbin, int limits, double peak_a, double min_a, double max_a, double max_a_hs, double EbminkT, double EbmaxkT, int ignoreradii, int binmf, long nplanets, double highmassfb, double masscrit, double highmass_aindex, double highmass_mratio_low)
{
	int success, success_a;
	long i, j, planetntry;
	double mass, Mmin, Mmax, amin, amax, W, vorb, emax, Mtotnew, aminroche, mlow, mhigh, norm;
	double eta, temp, temp1, temp_a, temp_W, temp_k;
	double tempm0, tempm1, tempa, tempe;
	double dummy, dummy1, dummy2;
	double limval1, limval2, limfunc, tempo_a, tempo_fa, tempo_lim, tempamin_au, tempamax_au;
	double kTcore, vcore, Eb, Ebmin, Ebmax, timeunitcgs, mtotal, m1, m2, X, qbin;
	double *r, *sigma, *mave, dtp, tphysf;
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
		int BSE_RTMSFLAG= 0;
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
	bse_set_rtmsflag(BSE_RTMSFLAG);
	bse_set_bhflag(BSE_BHFLAG);
        bse_set_bhms_coll_flag(BSE_BHMS_COLL_FLAG);
	bse_set_remnantflag(BSE_REMNANTFLAG);
        bse_set_grflag(BSE_GRFLAG);
        bse_set_kickflag(BSE_KICKFLAG);
        bse_set_zsun(BSE_ZSUN);
        bse_set_rembar_massloss(BSE_REMBAR_MASSLOSS);
	bse_set_bhspinflag(BSE_BHSPINFLAG);
	bse_set_bhspinmag(BSE_BHSPINMAG);
	bse_set_mxns(BSE_MXNS);
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

	/*if reading binary properties from a file then open than file to read*/
	FILE *rfile;
	if (binmf==99){
		rfile = fopen("readbinary.dat", "r");
	}
	long nhighmass = 0;
	long nhighmassbinary = 0;
	if (binmf==98){
		for (i=0; i<=cfd->NOBJ; i++){
			//fprintf (stdout, "m=%g\n", cfd->obj_m[i]);
			if (cfd->obj_m[i]>masscrit/cfd->Mclus) {
				nhighmass++;
			}
		}
		nhighmassbinary = (long) floor(nhighmass*highmassfb);
		fprintf (stdout, "nhighmass=%ld nhighmassbinary=%ld clustermass=%g\n", nhighmass, nhighmassbinary, cfd->Mclus);
	}

	
	/* set parameters relating to metallicity */
	zpars = (double *) malloc(20 * sizeof(double));
	bse_zcnsts(&METALLICITY, zpars);
	
	/* set collisions matrix */
	bse_instar();

	/* set newstarid to emulate CMC function */
	newstarid = cfd->NOBJ;

	/* test for strange things */
	if (Nbin+nplanets+nhighmassbinary > cfd->NOBJ) {
		fprintf(stderr, "Error: You've requested more binaries than there are objects in the input file!\n");
		exit(1);
	}

	if (cfd->NBINARY != 0) {
		fprintf(stderr, "Warning: NBINARY!=0 in input FITS file.  Be sure you know what you're doing!\n");
	}

	cfd->NBINARY = Nbin+nplanets+nhighmassbinary;
	
	/* reallocate memory for binaries */
	cfd->bs_index = (long *) realloc(cfd->bs_index, (nhighmassbinary+nplanets+Nbin+1)*sizeof(long));
	cfd->bs_id1 = (long *) realloc(cfd->bs_id1, (nhighmassbinary+nplanets+Nbin+1)*sizeof(long));
	cfd->bs_k1 = (int *) realloc(cfd->bs_k1, (nhighmassbinary+nplanets+Nbin+1)*sizeof(int));
	cfd->bs_m1 = (double *) realloc(cfd->bs_m1, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));
	cfd->bs_Reff1 = (double *) realloc(cfd->bs_Reff1, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));
	cfd->bs_id2 = (long *) realloc(cfd->bs_id2, (nhighmassbinary+nplanets+Nbin+1)*sizeof(long));
	cfd->bs_k2 = (int *) realloc(cfd->bs_k2, (nhighmassbinary+nplanets+Nbin+1)*sizeof(int));
	cfd->bs_m2 = (double *) realloc(cfd->bs_m2, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));
	cfd->bs_Reff2 = (double *) realloc(cfd->bs_Reff2, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));
	cfd->bs_a = (double *) realloc(cfd->bs_a, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));
	cfd->bs_e = (double *) realloc(cfd->bs_e, (nhighmassbinary+nplanets+Nbin+1)*sizeof(double));

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
	long counthighmassbinary = 0;
	/*First assign high mass stars as binaries based on a separate binary fraction even before going to the random assignment. */
	if (binmf==98){
		//for (i=0; i<=cfd->NOBJ; i++){
		while (counthighmassbinary<nhighmassbinary){
			i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
			if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0 && cfd->obj_m[i]>masscrit/cfd->Mclus){
				counthighmassbinary++;
				j++;
				//dummy = rng_t113_dbl();
				//fprintf(stdout, "count=%ld; assigning binary for high mass=%g and i=%ld\n", counthighmassbinary, cfd->obj_m[i]*cfd->Mclus, i);
				//if (dummy<=highmassfb) {
				//fprintf(stdout, "count=%f; assigning binary for mass=%g", dummy, cfd->obj_m[i]*cfd->Mclus);
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
				if (highmass_mratio_low<=0){
					cfd->bs_m2[j] = Mmin + rng_t113_dbl() * (cfd->bs_m1[j] - Mmin);
					//fprintf (stderr, "highmass comp masses: %g %g %g\n", cfd->bs_m2[j]*cfd->Mclus, cfd->bs_m1[j]*cfd->Mclus, masscrit);
				} else if (highmass_mratio_low>1){
					dummy = masscrit;
					cfd->bs_m2[j] = dummy + rng_t113_dbl() * (cfd->bs_m1[j] - dummy);
					//fprintf (stderr, "highmass comp masses: %g %g %g\n", cfd->bs_m2[j]*cfd->Mclus, cfd->bs_m1[j]*cfd->Mclus, masscrit);
				} else {
					dummy = cfd->bs_m1[j]*highmass_mratio_low;
					cfd->bs_m2[j] = dummy + rng_t113_dbl() * (cfd->bs_m1[j] - dummy);
					//fprintf (stderr, "highmass comp masses: %g %g %g\n", cfd->bs_m2[j]*cfd->Mclus, cfd->bs_m1[j]*cfd->Mclus, masscrit);
				}
		
				/* stellar evolution stuff */
				star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			}
		}
		while (j < Nbin+nhighmassbinary) {
			i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
			if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0 && cfd->obj_m[i]<=masscrit/cfd->Mclus){
				j++;
				//dummy = rng_t113_dbl();
				//fprintf(stdout, "count=%ld; assigning binary for low mass=%g and i=%ld\n", counthighmassbinary, cfd->obj_m[i]*cfd->Mclus, i);
				//if (dummy<=highmassfb) {
				//fprintf(stdout, "count=%f; assigning binary for mass=%g", dummy, cfd->obj_m[i]*cfd->Mclus);
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
				//fprintf (stderr, "comp masses: %g %g %g\n", cfd->bs_m2[j]*cfd->Mclus, cfd->bs_m1[j]*cfd->Mclus,masscrit);	
					
				/* stellar evolution stuff */
				star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			}
		}
	}
	if (binmf==97){/*No BH progenotors are initially in a binary. */
		while (j < Nbin+nhighmassbinary){
			i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
			/* make it a binary if it's not already */
			if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0 && cfd->obj_m[i]<masscrit/cfd->Mclus) {
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
				star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;

			}
		}
	}
	if (binmf!=97 || binmf!=98){
	while (j < Nbin+nhighmassbinary) {
		i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
		
		/* make it a binary if it's not already */
		if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0) {
			j++;

			/* make this object a binary */
			cfd->obj_binind[i] = j;
			cfd->bs_index[j] = j;


			if (binmf==99){
				if(rfile==NULL){
					fprintf (stdout, "could not open file\n");
				} else{
					fscanf (rfile, "%lf %lf %lf %lf\n", &tempm0, &tempm1, &tempa, &tempe);
					fprintf (stdout, "binmf: %.6f %.6f %.6f %.6f\n", tempm0, tempm1, tempa, tempe);
					cfd->bs_id1[j] = cfd->obj_id[i];
					cfd->bs_id2[j] = star_get_id_new();

					cfd->obj_m[i] = tempm0 / cfd->Mclus;
					cfd->bs_m1[j] = cfd->obj_m[i];
					cfd->bs_m2[j] = tempm1 / cfd->Mclus;;
					
					/* stellar evolution stuff */
					/*set them one at a time: first for m2 then m1*/
					star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
					cfd->bs_k2[j] = star.se_k;
					cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
					cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);
					
					/*m2 stellar properties are set, now set m1 stellar properties*/
					star.se_mass = cfd->bs_m1[j] * cfd->Mclus;
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
					cfd->bs_k1[j] = star.se_k;
					cfd->bs_m1[j] = star.se_mass / cfd->Mclus;
					cfd->bs_Reff1[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);	

					/* set/unset stellar properties for obj */
					cfd->obj_id[i] = NOT_A_STAR;
					cfd->obj_k[i] = NOT_A_STAR;
					cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
					cfd->obj_Reff[i] = 0.0;

					/*assign orbital properties*/	
					cfd->bs_a[j] = tempa * AU / (cfd->Rvir * PARSEC);
					cfd->bs_e[j] = tempe;
				
				}
			}

			else if (binmf==0){   /*All single masses are chosen from the same MF.  Binary companion masses are chosen later.*/
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
				star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			} 
			/*new option*/ 
			else if (binmf == 1){  /*single masses are chosen from one MF. Binary total m is chosen from another MF*/
				/* copy properties for star 1 */
				success = 0;
				while (!success){
					X = rng_t113_dbl();
					mtotal = 0.33*( 1./( pow(1.-X,0.75) + 0.04 * pow(1-X,0.25) ) - pow(1-X,2.)/1.04 );
					mtotal = mtotal/cfd->Mclus;
					qbin = rng_t113_dbl();
					m1 = qbin*mtotal;
					m2 = mtotal - m1; 
					if (m1>=Mmin && m1 <=Mmax && m2>=Mmin && m2<=Mmax){
						success = 1;
					}
				}
				//fprintf(stderr, "m1=%g m2=%g\n", m1*cfd->Mclus, m2*cfd->Mclus);
				cfd->bs_id1[j] = cfd->obj_id[i];
				cfd->bs_id2[j] = star_get_id_new();
				if (m1>m2){
					cfd->obj_m[i] = m1;
					cfd->bs_m1[j] = cfd->obj_m[i];
					cfd->bs_m2[j] = m2;
				} else {
					cfd->obj_m[i] = m2;
					cfd->bs_m1[j] = cfd->obj_m[i];
					cfd->bs_m2[j] = m1;
				}
				//cfd->bs_id1[j] = cfd->obj_id[i];
				//cfd->bs_k1[j] = cfd->obj_k[i];
				//cfd->bs_m1[j] = cfd->obj_m[i];
				//cfd->bs_Reff1[j] = cfd->obj_Reff[i];

				/* assign properties for star 2;
			   	this assumes all stars in input file are id'ed sequentially from 1 */
				//cfd->bs_id2[j] = star_get_id_new();
				/* set secondary mass from dP/dq=1 distribution */
				//cfd->bs_m2[j] = Mmin + rng_t113_dbl() * (cfd->bs_m1[j] - Mmin);
			
				/* stellar evolution stuff */
				/*set them one at a time: first for m2 then m1*/
				star.se_mass = cfd->bs_m2[j] * cfd->Mclus;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/*m2 stellar properties are set, now set m1 stellar properties*/
				star.se_mass = cfd->bs_m1[j] * cfd->Mclus;
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
				cfd->bs_k1[j] = star.se_k;
				cfd->bs_m1[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff1[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			}
		}
	}
	}

	if (nplanets>0){
		planetntry = 0;
		while (j<(nplanets+Nbin) && planetntry<cfd->NOBJ){
			i = (long) floor(rng_t113_dbl() * cfd->NOBJ + 1.0);
			planetntry++;
			fprintf (stderr, "try=%ld j=%ld\n", planetntry, j);
			/* make it a binary if it's not already */
			if (i <= cfd->NOBJ && cfd->obj_binind[i] == 0 && cfd->obj_m[i]<(0.8/cfd->Mclus)) {
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
				/* set all secondary masses equal to 1 Jupiter*/
				//cfd->bs_m2[j] = 0.001 / cfd->Mclus;
				/*Assign masses so that df/dlogm = -0.48 equivalently df/dm = -1.48*/
				X = rng_t113_dbl();
				mlow = 3e-6;
				mhigh = 3e-3;
				norm = pow( (mhigh/mlow), -0.48 ) - 1.;
				m2 = mlow * pow( (norm*X + 1.), (1./-0.48));
				fprintf (stdout, "%.6g\n", m2);
				cfd->bs_m2[j] = m2 / cfd->Mclus ;
				/* stellar evolution stuff */
				star.se_mass = m2;
				star.se_k = 0;
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
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				//cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);
				cfd->bs_Reff2[j] = 0.16 / (cfd->Rvir * PARSEC / RSUN);  //forcing it to Jupiter radius

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			}
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
	fprintf (stderr, "Mtotal after addbinaries=%g\n", cfd->Mclus);
	
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
	addbin_calc_sigma_r(cfd, r, sigma, mave);

	/* calculate kT in cluster core to get scale for binary semimajor axis distribution */
	kTcore = 0.0;
	vcore = 0.0;
	for (i=1; i<=AVEKERNEL; i++) {
		vcore += sigma[i];
        fprintf(stderr,"sigma[%d]=%g\n",i,sigma[i]*(PARSEC*cfd->Rvir/1e5)/(pow(cfd->Rvir * PARSEC, 1.5) / sqrt(G * cfd->Mclus * MSUN)));
		kTcore += (1.0/3.0) * mave[i] * sigma[i] * sigma[i];
	}
	vcore /= AVEKERNEL;
	kTcore /= AVEKERNEL;
        printf("kTcore is %g\n", kTcore);

	/* assign binary parameters */
	if (limits == 0 && binmf!=99) {
		/* assign binaries physically */
		for (i=1; i<=cfd->NBINARY; i++) {
			if (cfd->bs_m2[i]*cfd->Mclus < 0.002){
				cfd->bs_a[i] = 5. * AU / (cfd->Rvir * PARSEC);
				cfd->bs_e[i] = 0.;
			} else {
				/* choose a from a distribution uniform in 1/a from near contact to hard/soft boundary */
				amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
				W = 4.0 * vcore / sqrt(3.0 * PI);
				vorb = XHS * W;
				amax = max_a_hs * (cfd->bs_m1[i]+cfd->bs_m2[i]) / (vorb * vorb);
                timeunitcgs = pow(cfd->Rvir * PARSEC, 1.5) / sqrt(G * cfd->Mclus * MSUN);
				if (amax <= amin && ignoreradii == 0) {
				  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g vorb = %g M=%g m1=%g m2=%g\n", amax*(PARSEC*cfd->Rvir)/AU, amin*(PARSEC*cfd->Rvir)/AU, vorb*(PARSEC*cfd->Rvir/1e5)/timeunitcgs, cfd->obj_m[i]*cfd->Mclus, cfd->bs_m1[i]*cfd->Mclus,cfd->bs_m2[i]*cfd->Mclus);
				  amax = amin;
				}
				cfd->bs_a[i] = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
	
				/* get eccentricity from thermal distribution, truncated near contact */
				if (ignoreradii == 0) {
				  emax = 1.0 - amin / cfd->bs_a[i];
				} else {
				  emax = 1.0;
				}
				cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
				//fprintf (stdout, "assigning flat in log a primary mass=%g a=%g e=%g\n", cfd->bs_m1[i]*cfd->Mclus, cfd->bs_a[i], cfd->bs_e[i]);
			}		
		}
	} else if (limits == 1 && binmf!=99) {
		/* assign binaries via kT description */
		/* set max and min binding energies */
		Ebmin = EbminkT * kTcore;
		Ebmax = EbmaxkT * kTcore;

		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose binding energy uniformly in the log */
			Eb = pow(10.0, rng_t113_dbl()*(log10(Ebmax)-log10(Ebmin))+log10(Ebmin));
			cfd->bs_a[i] = cfd->bs_m1[i] * cfd->bs_m2[i] / (2.0 * Eb);
			
			amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);

			if (cfd->bs_a[i] <= amin && ignoreradii == 0) {
			  fprintf(stderr, "WARNING: a <= amin! amax=%g amin=%g\n", cfd->bs_a[i], amin);
			  fprintf(stderr, "WARNING: setting a = amin\n");
			  cfd->bs_a[i] = amin;
			}
			
			if (ignoreradii == 0) {
			  emax = 1.0 - amin / cfd->bs_a[i];
			} else {
			  emax = 1.0;
			}
			cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
			/* get eccentricity from thermal distribution */
			/* cfd->bs_e[i] = sqrt(rng_t113_dbl()); */
		}
	} else if (limits == 2 && binmf!=99) {
		/* assign binaries sort of physically, similarly to the M67 model of Hurley, et al. (2005) */
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose a from a distribution uniform in 1/a from near contact to 50 AU */
			//amin = 2.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
			//amax = 50.0 * AU / (cfd->Rvir * PARSEC);
			amin = 5.0 * (cfd->bs_Reff1[i]);
			amax = 100.0 * AU / (cfd->Rvir * PARSEC);

			if (amax <= amin && ignoreradii == 0) {
			  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g\n", amax, amin);
			  fprintf(stderr, "WARNING: setting amax = amin\n");
			  amax = amin;
			}

			cfd->bs_a[i] = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
			
			/* get eccentricity from thermal distribution, truncated near contact */
			if (ignoreradii == 0) {
			  emax = 1.0 - amin / cfd->bs_a[i];
			} else {
			  emax = 1.0;
			}
			cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
		}
	} else if (limits == 3 && binmf!=99) {
		/* assign binaries in the same way as in Ivanova, et al. (2005), with many soft binaries */
	        timeunitcgs = pow(cfd->Rvir * PARSEC, 1.5) / sqrt(G * cfd->Mclus * MSUN);
		fprintf(stderr, "rvir=%g PC mclus=%g MSUN nbtime=%g YR\n", cfd->Rvir, cfd->Mclus, timeunitcgs/YEAR);
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose period from 0.1 d to 10^7 d, subject to roche-lobe limits */
		        aminroche = MAX(cfd->bs_Reff1[i] / roche(cfd->bs_m1[i], cfd->bs_m2[i]), 
					cfd->bs_Reff2[i] / roche(cfd->bs_m2[i], cfd->bs_m1[i]));
		        amin = MAX(pow((0.1/365.25*YEAR/timeunitcgs)*(0.1/365.25*YEAR/timeunitcgs)*(cfd->bs_m1[i]+cfd->bs_m2[i])/(4.0*PI*PI), 1.0/3.0), 
				   aminroche);
			amax = pow((1.0e7/365.25*YEAR/timeunitcgs)*(1.0e7/365.25*YEAR/timeunitcgs)*(cfd->bs_m1[i]+cfd->bs_m2[i])/(4.0*PI*PI), 1.0/3.0);

			if (amax <= amin && ignoreradii == 0) {
			  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g\n", amax, amin);
			  fprintf(stderr, "WARNING: setting amax = amin\n");
			  amax = amin;
			}

			cfd->bs_a[i] = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
			
			/* get eccentricity from thermal distribution, truncated at high end */
			if (ignoreradii == 0) {
			  emax = 1.0 - amin / cfd->bs_a[i];
			} else {
			  emax = 1.0;
			}
			cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
		}
	/* assign binary parameters */
	} else if (limits == 4 && binmf!=99) {
		/* assign binaries according to Hurley's prescription, lognormal */
		for (i=1; i<=cfd->NBINARY; i++) {
			if (cfd->bs_m2[i]*cfd->Mclus < 0.01){
				//fprintf (stdout, "came here\n");
				/*set a according to Eq. 16 Egleton, Fitchett & Tout 1989; use eta=0.33, mode of a=0.11, in the range 0.01-1.5 AU
				 * these choices worked reasonably well with the Kepler observed a dist.*/
				eta = 0.33;
				temp = 0.5*( pow(eta, 0.33) + pow(eta, -0.33) );
				temp_k = acos (1./temp);
	
				success_a = 0;
				while(!success_a){
					temp_W = (-1. + rng_t113_dbl()*(2.));
					temp1 = 1./cos(temp_k*temp_W) + tan(temp_k*temp_W) ;
					temp_a = 0.11 * pow(temp1, (1./0.33));
					temp_a = temp_a * AU / (cfd->Rvir * PARSEC);
					if (temp_a <= 1.5*AU/(cfd->Rvir*PARSEC) && temp_a >= 0.01*AU/(cfd->Rvir*PARSEC)){
						success_a = 1;
						fprintf (stdout, "a= %.6g\n", temp_a*cfd->Rvir*PARSEC/AU);
						cfd->bs_a[i] = temp_a;
					}
				}
				/*setting some fixed a value*/
				//cfd->bs_a[i] = 5. * AU / (cfd->Rvir * PARSEC);
				/*All planets have zero e; reasonable since there is a high abundance of zero e for close in planets. actually we don't know.*/
				cfd->bs_e[i] = 0.;
			} else {
				/* choose a from from near contact (5*(R1)) to 100 AU with 30 AU mode 
				 * this is for Hurley 07 papers for runs K100-5 and K100-10*/
				amin = 5.0 * (cfd->bs_Reff1[i]);
				amax = max_a * AU / (cfd->Rvir * PARSEC);
	
				if (amax <= amin && ignoreradii == 0) {
				  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g\n", amax, amin);
				  fprintf(stderr, "WARNING: setting amax = amin\n");
				  amax = amin;
				}
				
				/*set a according to Eq. 16 Egleton, Fitchett & Tout 1989*/
				eta = 0.001;
				temp = 0.5*( pow(eta, 0.33) + pow(eta, -0.33) );
				temp_k = acos (1./temp);
	
				success_a = 0;
				while(!success_a){
					temp_W = (-1. + rng_t113_dbl()*(2.));
					temp1 = 1./cos(temp_k*temp_W) + tan(temp_k*temp_W) ;
					temp_a = peak_a * pow(temp1, (1./0.33));
					temp_a = temp_a * AU / (cfd->Rvir * PARSEC);
					if (temp_a <= amax && temp_a >= amin){
						success_a = 1;
						cfd->bs_a[i] = temp_a;
					}
				}

				/* get eccentricity from thermal distribution, truncated near contact */
				if (ignoreradii == 0) {
				  emax = 1.0 - amin / cfd->bs_a[i];
				} else {
				  emax = 1.0;
				}
				/*Thermal eccentricity*/
				//cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
				/*Aaron insists not to put thermal e; Rather putting in flat e dist.*/
				cfd->bs_e[i] = emax * rng_t113_dbl();
			}
		}
	} else if (limits==5 && binmf!=99) { //Assign binary semimajor axis as a power law based on Sana et al. 2012, Science, 337, 444
		/* assign binaries physically */
		for (i=1; i<=cfd->NBINARY; i++) {
			if (cfd->bs_m2[i]*cfd->Mclus < 0.002){  //If planetary, then put them in circular orbits at 5 AU
				cfd->bs_a[i] = 5. * AU / (cfd->Rvir * PARSEC);
				cfd->bs_e[i] = 0.;
			} else if (cfd->bs_m1[i]*cfd->Mclus > masscrit) {  //If high-mass primary, > masscrit, use a user-defined power-law a-dist
				/* choose a from a powerlaw distribution from near contact to hard/soft boundary */
				
				amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
				W = 4.0 * vcore / sqrt(3.0 * PI);
				vorb = XHS * W;
				amax = max_a_hs * (cfd->bs_m1[i]+cfd->bs_m2[i]) / (vorb * vorb);
				if (amax <= amin && ignoreradii == 0) {
				  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g M=%g\n", amax, amin, cfd->obj_m[i]*cfd->Mclus);
				  fprintf(stderr, "WARNING: setting amax = amin\n");
				  amax = amin;
				  cfd->bs_a[i] = amin;
				} else {
					
					if (highmass_aindex==0) {
						dummy1 = rng_t113_dbl();
						cfd->bs_a[i] = amin + (amax-amin)*dummy1;  //If uniform in a
					} else if (highmass_aindex==-1) {
						fprintf (stderr, "Change your generating file and use limits 0 instead of 5.[%g]\n", highmass_aindex); //If uniform in log a
					} else { //use rejection method
						tempamin_au = amin * (cfd->Rvir * PARSEC) / AU;
						tempamax_au = amax * (cfd->Rvir * PARSEC) / AU;
						limval1 = pow(tempamin_au, highmass_aindex);
						limval2 = pow(tempamax_au, highmass_aindex);
						if (limval1>limval2){
							limfunc = limval1;
						} else {
							limfunc = limval2;
						}
						success_a = 0;
						while (success_a==0) {
							dummy1 = rng_t113_dbl();
							tempo_a = tempamin_au + (tempamax_au-tempamin_au)*dummy1;
							tempo_fa = pow(tempo_a, highmass_aindex);
							dummy2 = rng_t113_dbl();
							tempo_lim = limfunc*dummy2;
							if (tempo_lim<=tempo_fa) {
								success_a = 1;
								cfd->bs_a[i] = tempo_a * AU / (cfd->Rvir * PARSEC);
							}
							//norm1 = pow( (amax/amin), (highmass_aindex) ) - 1.;
							//cfd->bs_a[i] = amin * pow(1./(highmass_aindex-1), (norm1 * dummy1 + 1.)); 
							//fprintf (stdout, "dummy1=%g dummy2=%g tempa=%g tempfa=%g limfunc=%g amin=%g amax=%g aminpow=%g amaxpow=%g\n", dummy1, dummy2, tempo_a, tempo_fa, limfunc, amin, amax, limval1, limval2);
						}
					}
				}
				/* get eccentricity from thermal distribution, truncated near contact */
				if (ignoreradii == 0) {
				  emax = 1.0 - amin / cfd->bs_a[i];
				} else {
				  emax = 1.0;
				}
				cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
				fprintf (stdout, "assigning power-law a for high primary mass=%g a=%g e=%g\n", cfd->bs_m1[i]*cfd->Mclus, cfd->bs_a[i], cfd->bs_e[i]);
			} else {  //If not high-mass star, use the regular uniform in log a distribution. 
				/* choose a from a distribution uniform in 1/a from near contact to hard/soft boundary */
				amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
				W = 4.0 * vcore / sqrt(3.0 * PI);
				vorb = XHS * W;
				amax = max_a_hs * (cfd->bs_m1[i]+cfd->bs_m2[i]) / (vorb * vorb);
				if (amax <= amin && ignoreradii == 0) {
				  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g\n", amax, amin);
				  fprintf(stderr, "WARNING: setting amax = amin\n");
				  amax = amin;
				}
				cfd->bs_a[i] = pow(10.0, rng_t113_dbl()*(log10(amax)-log10(amin))+log10(amin));
	
				/* get eccentricity from thermal distribution, truncated near contact */
				if (ignoreradii == 0) {
				  emax = 1.0 - amin / cfd->bs_a[i];
				} else {
				  emax = 1.0;
				}
				cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
				fprintf (stdout, "assigning power-law a for low primary mass=%g a=%g e=%g\n", cfd->bs_m1[i]*cfd->Mclus, cfd->bs_a[i], cfd->bs_e[i]);
			}

		
		}
		

	} else if (binmf!=99) {
		fprintf(stderr, "limits=%d unknown!\n", limits);
		exit(1);
	}


	/*if binary properties were read from a file then close that file*/
	if (binmf==99){
		fclose(rfile);
	}

	free(sigma);
	free(r);
	free(mave);
}

int main(int argc, char *argv[]){
        int i, limits, ignoreradii, binary_mf;
	long Nbin, nplanets;
	unsigned long seed;
	double Ebmin, Ebmax, peak_a, min_a, max_a, max_a_hs, highmass_fb, mcrit, highmass_aind, highmass_mratio_low;
	cmc_fits_data_t cfd;
	char infilename[1024], outfilename[1024];
	const char *short_opts = "i:o:N:l:p:a:A:H:m:M:I:s:b:n:f:G:k:q:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"Nbin", required_argument, NULL, 'N'},
		{"limits", required_argument, NULL, 'l'},
		{"peak_a", required_argument, NULL, 'p'},
		{"min_a", required_argument, NULL, 'a'},
		{"max_a", required_argument, NULL, 'A'},
		{"max_a_hs", required_argument, NULL, 'H'},
		{"Ebmin", required_argument, NULL, 'm'},
		{"Ebmax", required_argument, NULL, 'M'},
		{"ignoreradii", no_argument, NULL, 'I'},
		{"seed", required_argument, NULL, 's'},
		{"binary_mf", required_argument, NULL, 'b'},
		{"Nplanets", required_argument, NULL, 'n'},
		{"highmass_fb", required_argument, NULL, 'f'},
		{"mcrit", required_argument, NULL, 'G'},
		{"highmass_aind", required_argument, NULL, 'k'},
		{"highmass_mratio_low", required_argument, NULL, 'q'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* argument defaults */
	sprintf(infilename, "%s", INFILE);
	sprintf(outfilename, "%s", OUTFILE);
	Nbin = NBIN;
	limits = LIMITS;
	peak_a = PEAK_A;
	min_a = MIN_A;
	max_a = MAX_A;
	max_a_hs = MAX_A_HS;
	Ebmin = EBMIN;
	Ebmax = EBMAX;
	ignoreradii = 0;
	seed = SEED;
	nplanets = NPLANETS;
	binary_mf = BINMF;
	highmass_fb = HIGHMASSFB;
	mcrit = MCRIT;
	highmass_aind = HIGHMASSAINDEX;
	highmass_mratio_low = HIGHMASSMRATIOLOW;
	

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
			break;
		case 'p':
			peak_a = strtod(optarg, NULL);
			fprintf (stderr, "peak of a = %g\n", peak_a);
			break;
		case 'a':
			min_a = strtod(optarg, NULL);
			fprintf (stderr, "minimum of a = %g\n", min_a);
			break;
		case 'A':
			max_a = strtod(optarg, NULL);
			fprintf (stderr, "maximum of a = %g\n", max_a);
			break;
		case 'H':
                        max_a_hs = strtod(optarg, NULL);
                        //fprintf (stderr, "maximum a is %g times HS boundary\n", max_a_hs);
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
		case 'I':
		        ignoreradii = 1;
			break;
		case 's':
			seed = strtol(optarg, NULL, 10);
			break;
		case 'b':
			binary_mf = strtol(optarg, NULL, 10);
			if (binary_mf!=0 && binary_mf!=1 && binary_mf!=99 && binary_mf!=98 && binary_mf!=97){
				fprintf(stderr, "WARNING: binary_mf=%d must be either 0, 1, 97, 98, 99!\n", binary_mf);
			}
			break;
		case 'n':
			nplanets = strtol(optarg, NULL, 10);
			if (nplanets < 0) {
				fprintf(stderr, "nplanets must be >=0!\n");
				exit(1);
			}
			break;
		case 'f':
			highmass_fb = strtod(optarg, NULL);
			if (highmass_fb<0 || highmass_fb>1){
				fprintf (stderr, "highmass binary fraction must bee >=0 and <=1![%g]\n", highmass_fb);
			} else if (binary_mf!=98) {
				fprintf (stderr, "binary_mf is expected to be 98. Do you know what your are doing?\n");
			}
			break;
		case 'G':
			mcrit = strtod(optarg, NULL);
			if (mcrit<1 || mcrit>150){
				fprintf (stderr, "mcrit for highmass binary fraction must bee >=1 and <=150![%g]\n", mcrit);
			} else if (binary_mf!=98) {
				fprintf (stderr, "binary_mf is expected to be 98. Do you know what your are doing?\n");
			}
			break;
		case 'k':
			highmass_aind = strtod(optarg, NULL);
			if (highmass_aind==-1){
				fprintf (stderr, "should directly use l=0\n");
			} else if (binary_mf!=98) {
				fprintf (stderr, "looks like power-law dist. is desired. Make sure to choose b=98.\n");
			}
			break;
		case 'q':
			highmass_mratio_low = strtod(optarg, NULL);
			if (highmass_mratio_low<=0){
				fprintf (stderr, "using default value of minimum stellar mass as secondary mass lower limit\n");
			}
			else{
				fprintf(stderr, "using the user defined lower limit in secondary mass based on the limiting mass ratio. Only for high-mass stars. %g\n", highmass_mratio_low);
			}
			break;
		case 'h':
			fprintf (stdout, "Came here: 2\n");
			print_usage(stdout);
			return(0);
		default:
			break;
		}
	}

	/* check to make sure there was nothing crazy on the command line */
	if (optind < argc) {
		fprintf (stdout, "Came here: 1\n");
		print_usage(stdout);
		return(1);
	}

	reset_rng_t113(seed);

	cmc_read_fits_file(infilename, &cfd, 0);

	METALLICITY = cfd.Z;
	//fprintf(stderr, "l=%d, p=%g AU, a=%g AU, A=%g AU, m=%g M=%g b=%d\n", limits, peak_a, min_a, max_a, Ebmin, Ebmax, binary_mf);
	assign_binaries(&cfd, Nbin, limits, peak_a, min_a, max_a, max_a_hs, Ebmin, Ebmax, ignoreradii, binary_mf, nplanets, highmass_fb, mcrit, highmass_aind, highmass_mratio_low);

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

