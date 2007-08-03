/* -*- linux-c -*- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <fitsio.h>
#include "../common/fitslib.h"

int main(int argc, char *argv[]){
	cmc_fits_data_t cfd;
	int i;

	cfd.NOBJ=100;
	cfd.NBINARY=10;

	cmc_malloc_fits_data_t(&cfd);

	cfd.Mclus = 100.0;
	cfd.Rvir = 1.0;
	cfd.Rtid = 100.0;
	cfd.Z = 0.01;

	/* assign objs */
	cfd.obj_id[0] = 0;
	cfd.obj_k[0] = 0;
	cfd.obj_m[0] = 0.0;
	cfd.obj_Reff[0] = 0.0;
	cfd.obj_r[0] = 0.0;
	cfd.obj_vr[0] = 0.0;
	cfd.obj_vt[0] = 0.0;
	cfd.obj_binind[0] = 0;

	for (i=1; i<=cfd.NOBJ; i++) {
		cfd.obj_id[i] = i;
		cfd.obj_k[i] = 0;
		cfd.obj_m[i] = 0.01;
		cfd.obj_Reff[i] = 1.0e-6;
		cfd.obj_r[i] = 1.0;
		cfd.obj_vr[i] = -1.3;
		cfd.obj_vt[i] = 3.4;
		cfd.obj_binind[i] = 3;
	}

	cfd.obj_id[cfd.NOBJ+1] = 0;
	cfd.obj_k[cfd.NOBJ+1] = 0;
	cfd.obj_m[cfd.NOBJ+1] = 0.0;
	cfd.obj_Reff[cfd.NOBJ+1] = 0.0;
	cfd.obj_r[cfd.NOBJ+1] = 0.0;
	cfd.obj_vr[cfd.NOBJ+1] = 0.0;
	cfd.obj_vt[cfd.NOBJ+1] = 0.0;
	cfd.obj_binind[cfd.NOBJ+1] = 0;

	/* assign binary stars */
	cfd.bs_index[0] = 0;
	cfd.bs_id1[0] = 0;
	cfd.bs_k1[0] = 0;
	cfd.bs_m1[0] = 0.0;
	cfd.bs_Reff1[0] = 0.0;
	cfd.bs_id2[0] = 0;
	cfd.bs_k2[0] = 0;
	cfd.bs_m2[0] = 0.0;
	cfd.bs_Reff2[0] = 0.0;
	cfd.bs_a[0] = 0.0;
	cfd.bs_e[0] = 0.0;
	
	for (i=1; i<=cfd.NBINARY; i++) {
		cfd.bs_index[i] = i;
		cfd.bs_id1[i] = cfd.obj_id[i];
		cfd.bs_k1[i] = 0;
		cfd.bs_m1[i] = 0.005;
		cfd.bs_Reff1[i] = 0.0;
		cfd.bs_id2[i] = cfd.NOBJ+i;
		cfd.bs_k2[i] = 0;
		cfd.bs_m2[i] = 0.005;
		cfd.bs_Reff2[i] = 0.0;
		cfd.bs_a[i] = 1.0e-2;
		cfd.bs_e[i] = 0.67;
	}

	cmc_write_fits_file(&cfd, "test.fits");

	cmc_free_fits_data_t(&cfd);

	return 0;
}
