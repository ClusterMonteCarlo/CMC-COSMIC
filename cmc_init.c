/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"

void assign_binaries(void)
{
	long i, j;
	double Mmin, Mave, Mmax, amin, amax, W, vorb, emax, Mtotnew;
	double kTcore, Eb, Ebmin, Ebmax;
	
	/* initialize a variable that should probably be initialized in main() */
	N_b = clus.N_BINARY;
	
	if (clus.N_BINARY > 0) {
		diaprintf("assigning %ld binaries...\n", clus.N_BINARY);
	}

	/* calculate minimum mass of mass function, plus other statistics */
	Mmin = GSL_POSINF;
	Mmax = GSL_NEGINF;
	Mave = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		if (star[i].m < Mmin) {
			Mmin = star[i].m;
		}
		if (star[i].m > Mmax) {
			Mmax = star[i].m;
		}
		Mave += star[i].m;
	}
	Mave /= clus.N_STAR;

	/* report on statistics */
	diaprintf("Mmin=%g MSUN (%g) Mmax=%g MSUN (%g) Mave=%g MSUN (%g)\n", 
		Mmin*units.mstar/MSUN, Mmin,
		Mmax*units.mstar/MSUN, Mmax,
		Mave*units.mstar/MSUN, Mave);

	/* first assign all secondary masses, since this will change the units */
	j = 0;
	while (j < clus.N_BINARY) {
		i = (long) floor(rng_t113_dbl()*clus.N_STAR+1.0);
		
		/* make it a binary if it's not already */
		if (i <= clus.N_STAR && star[i].binind == 0) {
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
	Mtotnew = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		Mtotnew += star[i].m * madhoc;
	}
	
	for (i=1; i<=clus.N_STAR; i++) {
		star[i].m /= Mtotnew;
		if (star[i].binind) {
			binary[star[i].binind].m1 /= Mtotnew;
			binary[star[i].binind].m2 /= Mtotnew;
		}
	}
	
	fprintf(stdout, "Rescaling SOLAR_MASS_DYN from %g to", SOLAR_MASS_DYN);
	SOLAR_MASS_DYN /= Mtotnew;
	fprintf(stdout, " %g since creating binaries increases total cluster mass.\n", SOLAR_MASS_DYN);

	/* reset units */
	units_set();

	/* let's print out the mass statistics to make sure we did everything right */
	Mmin = GSL_POSINF;
	Mmax = GSL_NEGINF;
	Mave = 0.0;
	for (i=1; i<=clus.N_STAR; i++) {
		if (star[i].m < Mmin) {
			Mmin = star[i].m;
		}
		if (star[i].m > Mmax) {
			Mmax = star[i].m;
		}
		Mave += star[i].m;
	}
	Mave /= clus.N_STAR;
	
	/* report on statistics */
	diaprintf("Mmin=%g MSUN (%g) Mmax=%g MSUN (%g) Mave=%g MSUN (%g)\n", 
		Mmin*units.mstar/MSUN, Mmin,
		Mmax*units.mstar/MSUN, Mmax,
		Mave*units.mstar/MSUN, Mave);

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
	if (1) {
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
		Ebmin = 10.0 * kTcore;
		Ebmax = 100.0 * kTcore;

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

/* print out initial binary paramaeters to data file */
void print_initial_binaries(void)
{
	long i, j;
	char outfile[1024];
	FILE *initbinfile;

	/* open file for writing */
	sprintf(outfile, "%s.initbin.dat", outprefix);
	if ((initbinfile = fopen(outfile, "w")) == NULL) {
		eprintf("cannot create initbin file \"%s\".\n", outfile);
		exit(1);
	}
	
	/* and write data */
	fprintf(initbinfile, "# m0 [MSUN]  m1 [MSUN]  R0 [RSUN]  R1 [RSUN]  id0  id1  a [AU]  e\n");
	for (i=1; i<=clus.N_STAR; i++) {
		j = star[i].binind;
		if (j != 0) {
			fprintf(initbinfile, "%g %g %g %g %ld %ld %g %g\n",
				binary[j].m1 * units.mstar / MSUN, binary[j].m2 * units.mstar / MSUN, 
				binary[j].rad1 * units.l / RSUN, binary[j].rad2 * units.l / RSUN, 
				binary[j].id1, binary[j].id2,
				binary[j].a * units.l / AU, 
				binary[j].e);
		}
	}

	fclose(initbinfile);
}
