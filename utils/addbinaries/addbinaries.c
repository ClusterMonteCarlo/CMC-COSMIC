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
#define NBIN 0L
#define LIMITS 0
#define EBMIN 3.0
#define EBMAX 133.0
#define SEED 8732UL
#define AVEKERNEL 20
#define BINMF 0

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
	fprintf(stream, "  -l --limits <limits algorithm> : algorithm for setting limits on binary semimajor axes (0=physical, 1=kT prescription, 2=M67 model of Hurley, et al. (2005), 3=Ivanova, et al. (2005)) [%d]\n", LIMITS);
	fprintf(stream, "  -m --Ebmin <E_b,min>           : minimum binding energy, in kT [%g]\n", EBMIN);
	fprintf(stream, "  -M --Ebmax <E_b,max>           : maximum binding energy, in kT [%g]\n", EBMAX);
	fprintf(stream, "  -I --ignoreradii               : ignore radii when setting binary properties\n");
	fprintf(stream, "  -s --seed <seed>               : random seed [%ld]\n", SEED);
	fprintf(stream, "  -b --binary_mf <choosing bin MF> : set a different MF (KTG91 for now) for binaries [%d]\n", BINMF);	
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
		
		/* store sigma */
		r[si] = cfd->obj_r[si];
		sigma[si] = Sigma;
		mave[si] = Mave;
	}
}

void assign_binaries(cmc_fits_data_t *cfd, long Nbin, int limits, double EbminkT, double EbmaxkT, int ignoreradii, int binmf)
{
	int success, success_a;
	long i, j, it;
	double mass, Mmin, Mmax, amin, amax, W, vorb, emax, Mtotnew, aminroche;
	double norm_fa, binwidth, temp_a, temp_fa, temp_n;
	double kTcore, vcore, Eb, Ebmin, Ebmax, timeunitcgs, mtotal, m1, m2, X, qbin;
	double *r, *sigma, *mave, dtp, tphysf;
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

			if (binmf==0){   /*All single masses are chosen from the same MF.  Binary companion masses are chosen later.*/
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
				bse_evolv1(&(star.se_k), &(star.se_mass), &(star.se_mt), &(star.se_radius), 
					   &(star.se_lum), &(star.se_mc), &(star.se_rc), &(star.se_menv), 
				   	&(star.se_renv), &(star.se_ospin), &(star.se_epoch), &(star.se_tms), 
				   	&(star.se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);

				/* setting star properties in FITS file, being careful with units */
				cfd->bs_k2[j] = star.se_k;
				cfd->bs_m2[j] = star.se_mass / cfd->Mclus;
				cfd->bs_Reff2[j] = star.se_radius / (cfd->Rvir * PARSEC / RSUN);

				/* set/unset stellar properties for obj */
				cfd->obj_id[i] = NOT_A_STAR;
				cfd->obj_k[i] = NOT_A_STAR;
				cfd->obj_m[i] = cfd->bs_m1[j] + cfd->bs_m2[j];
				cfd->obj_Reff[i] = 0.0;
			} /*new option*/ else if (binmf == 1){  /*single masses are chosen from one MF. Binary total m is chosen from another MF*/
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
				bse_evolv1(&(star.se_k), &(star.se_mass), &(star.se_mt), &(star.se_radius), 
					   &(star.se_lum), &(star.se_mc), &(star.se_rc), &(star.se_menv), 
				   	&(star.se_renv), &(star.se_ospin), &(star.se_epoch), &(star.se_tms), 
				   	&(star.se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);

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
				bse_evolv1(&(star.se_k), &(star.se_mass), &(star.se_mt), &(star.se_radius), 
					   &(star.se_lum), &(star.se_mc), &(star.se_rc), &(star.se_menv), 
				   	&(star.se_renv), &(star.se_ospin), &(star.se_epoch), &(star.se_tms), 
				   	&(star.se_tphys), &tphysf, &dtp, &METALLICITY, zpars, vs);

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
	addbin_calc_sigma_r(cfd, r, sigma, mave);

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
	if (limits == 0) {
		/* assign binaries physically */
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose a from a distribution uniform in 1/a from near contact to hard/soft boundary */
			amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
			W = 4.0 * vcore / sqrt(3.0 * PI);
			vorb = XHS * W;
			amax = cfd->obj_m[i] / (vorb * vorb);
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
	} else if (limits == 1) {
		/* assign binaries via kT description */
		/* set max and min binding energies */
		Ebmin = EbminkT * kTcore;
		Ebmax = EbminkT * kTcore;

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
	} else if (limits == 2) {
		/* assign binaries sort of physically, similarly to the M67 model of Hurley, et al. (2005) */
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose a from a distribution uniform in 1/a from near contact to 50 AU */
			amin = 2.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
			amax = 50.0 * AU / (cfd->Rvir * PARSEC);

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
	} else if (limits == 3) {
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
	} else if (limits == 4) {
		/* assign binaries physically */
		for (i=1; i<=cfd->NBINARY; i++) {
			/* choose a from from near contact (5*(R1+R2)) to 100 AU with 30 AU mode 
			 * this is for Hurley 07 papers for runs K100-5 and K100-10*/
			//amin = 5.0 * (cfd->bs_Reff1[i] + cfd->bs_Reff2[i]);
			amin = 5.0 * (cfd->bs_Reff1[i]);
			amax = 100.0 * AU / (cfd->Rvir * PARSEC);
			if (amax <= amin && ignoreradii == 0) {
			  fprintf(stderr, "WARNING: amax <= amin! amax=%g amin=%g\n", amax, amin);
			  fprintf(stderr, "WARNING: setting amax = amin\n");
			  amax = amin;
			}
			
			/*set a according to Eq. 16 Egleton, Fitchett & Tout 1989*/
			
			/*Find normalization in the range amin-amax for Eq. 16 EFT 89*/
			//norm_fa = 0.;
			//binwidth = (amin-amax)/100.;
			//for (it=0; it<100; it++){
			//	temp_a = amin + i*binwidth;
			//	temp_fa = 0.33/( pow(temp_a/30., 0.33) + pow(30./temp_a, 0.33) );
			//	norm_fa += temp_fa*binwidth;
			//}
			/*now find a distribution from the normalized formula*/
			success_a = 0;
			while(!success_a){
				temp_a = (amin + rng_t113_dbl()*(amax-amin)) * (cfd->Rvir * PARSEC)/AU;
				temp_fa = 0.33/( pow(temp_a/30., 0.33) + pow(30./temp_a, 0.33) );
				temp_n = rng_t113_dbl() * 0.2;
				if (temp_n <= temp_fa){
					success_a = 1;
					cfd->bs_a[i] = temp_a * AU / (cfd->Rvir * PARSEC);
				}
			}

			/* get eccentricity from thermal distribution, truncated near contact */
			if (ignoreradii == 0) {
			  emax = 1.0 - amin / cfd->bs_a[i];
			} else {
			  emax = 1.0;
			}
			cfd->bs_e[i] = emax * sqrt(rng_t113_dbl());
		}
	} else {
		fprintf(stderr, "limits=%d unknown!\n", limits);
		exit(1);
	}

	free(sigma);
	free(r);
	free(mave);
}

int main(int argc, char *argv[]){
        int i, limits, ignoreradii, binary_mf;
	long Nbin;
	unsigned long seed;
	double Ebmin, Ebmax;
	cmc_fits_data_t cfd;
	char infilename[1024], outfilename[1024];
	const char *short_opts = "i:o:N:l:m:M:Is:b:h";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"Nbin", required_argument, NULL, 'N'},
		{"limits", required_argument, NULL, 'l'},
		{"Ebmin", required_argument, NULL, 'm'},
		{"Ebmax", required_argument, NULL, 'M'},
		{"ignoreradii", no_argument, NULL, 'I'},
		{"seed", required_argument, NULL, 's'},
		{"binary_mf", required_argument, NULL, 'b'},
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
	ignoreradii = 0;
	seed = SEED;
	binary_mf = BINMF;

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
			fprintf(stderr, "binary_mf=%d must be either 0 or 1!\n", binary_mf);
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
	fprintf(stderr, "l=%d, m=%g M=%g b=%d\n", limits, Ebmin, Ebmax, binary_mf);
	assign_binaries(&cfd, Nbin, limits, Ebmin, Ebmax, ignoreradii, binary_mf);

	cmc_write_fits_file(&cfd, outfilename);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

