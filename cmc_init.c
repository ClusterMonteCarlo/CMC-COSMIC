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
	double kTcore, Msec, Eb, Ebmin, Ebmax;
	
	dprintf("assigning %ld binaries...\n", clus.N_BINARY);

	/* initialize a few variables that should probably be initialized in main() */
	N_b = clus.N_BINARY;
	
	/* calculate core temperature to get scale for binary semimajor axis distribution */
	kTcore = 0.0;
	for (i=1; i<=NUM_CENTRAL_STARS; i++) {
		kTcore += (1.0/3.0) * (star[i].m/clus.N_STAR) * (sqr(star[i].vr) + sqr(star[i].vt));
	}
	kTcore = kTcore/NUM_CENTRAL_STARS;
	
	/* set max and min binding energies */
	Ebmin = 1.0 * kTcore;
	Ebmax = 133.0 * kTcore;
	
	/* assign binaries, uniformly distributed */
	M_b = 0.0;
	Msec = 0.0;
	j = 0;
	while (j < clus.N_BINARY) {
		i = (long) floor(rng_t113_dbl()*(clus.N_STAR-1)+1.0);
		
		/* make it a binary if it's not already */
		/* the extra paranoia i!=clus.N_STAR is just in case rng_t113_dbl() returns 
		   [0,1] and not [0,1) */
		if (i != clus.N_STAR && star[i].binind == 0) {
			j++;
			
			/* make this star a binary */
			star[i].binind = j;
			binary[j].inuse = 1;
			binary[j].id1 = star[i].id;
			binary[j].id2 = star_get_id_new();

			/* set masses from mass ratio distribution */
			binary[j].m1 = star[i].m;
			binary[j].m2 = MMIN / (initial_total_mass / SOLAR_MASS_DYN) 
			             + rng_t113_dbl() * (binary[j].m1 - MMIN / (initial_total_mass / SOLAR_MASS_DYN));
			
			star[i].m = binary[j].m1 + binary[j].m2;
			M_b += star[i].m;
			
			/* Msec is in dynamical units */
			Msec += binary[j].m2 * initial_total_mass / clus.N_STAR;
			
			/* set radii now that the masses are known */
			binary[j].rad1 = r_of_m(binary[j].m1);
			binary[j].rad2 = r_of_m(binary[j].m2);

			binary[j].Eint1 = 0.0;
			binary[j].Eint2 = 0.0;
			
			/* choose binding energy uniformly in the log */
			/* FIXME: should truncate at contact */
			Eb = pow(10.0, rng_t113_dbl()*(log10(Ebmax)-log10(Ebmin))+log10(Ebmin));
			binary[j].a = (binary[j].m1/clus.N_STAR) * (binary[j].m2/clus.N_STAR) / (2.0 * Eb);
			
			/* get eccentricity from thermal distribution */
			/* FIXME: should truncate at contact */
			binary[j].e = sqrt(rng_t113_dbl());
		}
	}

	/* Adding mass due to secondary stars causes the TOTAL cluster mass to 
	   increase. Hence the average mass also changes. Hence ALL masses must be
	   changed to the new unit of average mass. */
	for (i=0; i<=clus.N_STAR; i++) {
		star[i].m *= initial_total_mass / (initial_total_mass + Msec);
		if (star[i].binind != 0) {
			binary[star[i].binind].m1 *= initial_total_mass / (initial_total_mass + Msec);
			binary[star[i].binind].m2 *= initial_total_mass / (initial_total_mass + Msec);
		}
	}
	M_b *= initial_total_mass / (initial_total_mass + Msec);
	
	initial_total_mass += Msec;
	
	/* update binding energy */
	E_b = 0.0;
	for (i=0; i<=clus.N_STAR; i++) {
		if (star[i].binind != 0 && binary[star[i].binind].inuse) {
			E_b += (binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
				(2.0 * binary[star[i].binind].a);
		}
	}

	dprintf("done assigning binaries\n");
}

/* set binary parameters to specific values */
void assign_binaries_test(void)
{
	long i, j;
	
	dprintf("assigning %ld binaries (test)...\n", clus.N_BINARY);

	/* initialize a few variables that should probably be initialized in main() */
	N_b = clus.N_BINARY;
	
	/* assign binaries, uniformly distributed */
	M_b = 0.0;
	j = 0;
	while (j < clus.N_BINARY) {
		i = (long) floor(rng_t113_dbl()*(clus.N_STAR-1)+1.0);
		
		/* make it a binary if it's not already */
		/* the extra paranoia i!=clus.N_STAR is just in case rng_t113_dbl() returns 
		   [0,1] and not [0,1) */
		if (i != clus.N_STAR && star[i].binind == 0) {
			j++;
			
			/* make this star a binary */
			star[i].binind = j;
			binary[j].inuse = 1;
			binary[j].id1 = star[i].id;
			binary[j].id2 = star_get_id_new();

			/* set masses */
			binary[j].m1 = star[i].m/2.0;
			binary[j].m2 = binary[j].m1;

			star[i].m = binary[j].m1 + binary[j].m2;
			M_b += star[i].m;
			
			/* set radii now that the masses are known */
			binary[j].rad1 = r_of_m(binary[j].m1);
			binary[j].rad2 = r_of_m(binary[j].m2);
			
			binary[j].Eint1 = 0.0;
			binary[j].Eint2 = 0.0;
			
			binary[j].a = 1.0 * AU / units.l;
			binary[j].e = 0.0;
		}
	}

	/* update binding energy */
	E_b = 0.0;
	for (i=0; i<=clus.N_STAR; i++) {
		if (star[i].binind != 0 && binary[star[i].binind].inuse) {
			E_b += (binary[star[i].binind].m1/clus.N_STAR) * (binary[star[i].binind].m2/clus.N_STAR) / 
				(2.0 * binary[star[i].binind].a);
		}
	}

	dprintf("done assigning binaries\n");
}
