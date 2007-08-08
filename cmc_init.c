/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"


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
