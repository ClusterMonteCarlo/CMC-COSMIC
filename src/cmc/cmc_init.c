/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"


/**
* @brief print out initial binary paramaeters to data file
*/
void print_initial_binaries(void)
{
	long i, j;
	char outfile[1024];

    sprintf(outfile, "%s.initbin.dat", outprefix);

#ifdef USE_MPI
	 //MPI: Open corresponding MPI files, and declare buffers reqd for parallel write.
    MPI_File mpi_initbinfile;
    char mpi_initbinfile_buf[10000], mpi_initbinfile_wrbuf[10000000];
    long long mpi_initbinfile_len=0, mpi_initbinfile_ofst_total=0;
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_initbinfile);
    MPI_File_set_size(mpi_initbinfile, 0);
#else
    /* open file for writing */
	FILE *initbinfile;
    if ((initbinfile = fopen(outfile, "w")) == NULL) {
        eprintf("cannot create initbin file \"%s\".\n", outfile);
        exit(1);
    }
#endif

	/* and write data */
	//MPI: Header printed only by the root node.
	pararootfprintf(initbinfile, "# m0 [MSUN]  m1 [MSUN]  R0 [RSUN]  R1 [RSUN]  id0  id1  a [AU]  e\n");

	for (i=1; i<=clus.N_MAX_NEW; i++) {
		j = star[i].binind;
		if (j != 0) {
			parafprintf(initbinfile, "%g %g %g %g %ld %ld %g %g\n",
				binary[j].m1 * units.mstar / MSUN, binary[j].m2 * units.mstar / MSUN, 
				binary[j].rad1 * units.l / RSUN, binary[j].rad2 * units.l / RSUN, 
				binary[j].id1, binary[j].id2,
				binary[j].a * units.l / AU, 
				binary[j].e);
		}
	}

#ifdef USE_MPI
	 //MPI: Write in parallel
    mpi_para_file_write(mpi_initbinfile_wrbuf, &mpi_initbinfile_len, &mpi_initbinfile_ofst_total, &mpi_initbinfile);
    MPI_File_close(&mpi_initbinfile);
#else
	fclose(initbinfile);
#endif
}
