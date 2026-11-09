/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#ifdef USE_FITS
#include <fitsio.h>
#endif
#include <gsl/gsl_rng.h>
#include "cmc.h"
#include "cmc_vars.h"

/* Rank-local copy of the binary input rows this rank needs, indexed by LOCAL binary row:
   position 0 is global binary row 0 (the sentinel); positions 1..n-1 are the consecutive
   global rows first..first+n-2. Only the fields cmc uses are kept; types match cfd's bs_*
   arrays so the double -> long conversions stay two-step as before. For FITS input the
   slice aliases cfd's full-length arrays instead (owned == 0; see bin_slice_alias_cfd). */
typedef struct {
	long    n;
	long    first;
	int     owned;
	long    *id1;
	long    *id2;
	double  *m1;
	double  *m2;
	double  *Reff1;
	double  *Reff2;
	double  *a;
	double  *e;
} bin_slice_t;

static void bin_slice_alloc(bin_slice_t *s, long n, long first)
{
	s->n     = n;
	s->first = first;
	s->owned = 1;
	s->id1   = (long *)   calloc(n, sizeof(long));
	s->id2   = (long *)   calloc(n, sizeof(long));
	s->m1    = (double *) calloc(n, sizeof(double));
	s->m2    = (double *) calloc(n, sizeof(double));
	s->Reff1 = (double *) calloc(n, sizeof(double));
	s->Reff2 = (double *) calloc(n, sizeof(double));
	s->a     = (double *) calloc(n, sizeof(double));
	s->e     = (double *) calloc(n, sizeof(double));
	if (!s->id1 || !s->id2 || !s->m1 || !s->m2 || !s->Reff1 || !s->Reff2 || !s->a || !s->e) {
		eprintf("cannot allocate binary slice of %ld rows\n", n);
		exit_cleanly(-1, __FUNCTION__);
	}
}

/* Points the slice at cfd's full-length binary arrays (FITS input) without copying, so
   local row g is global row g. cfd stays allocated until after get_star_data() returns. */
static void bin_slice_alias_cfd(bin_slice_t *s)
{
	s->n     = cfd.NBINARY + 1;
	s->first = 1;
	s->owned = 0;
	s->id1   = cfd.bs_id1;
	s->id2   = cfd.bs_id2;
	s->m1    = cfd.bs_m1;
	s->m2    = cfd.bs_m2;
	s->Reff1 = cfd.bs_Reff1;
	s->Reff2 = cfd.bs_Reff2;
	s->a     = cfd.bs_a;
	s->e     = cfd.bs_e;
}

static void bin_slice_free(bin_slice_t *s)
{
	if (!s->owned)
		return;
	free(s->id1); free(s->id2); free(s->m1); free(s->m2);
	free(s->Reff1); free(s->Reff2); free(s->a); free(s->e);
}

/* Local slice row of global binary row g, or -1 if outside this rank's range. */
static long bin_local_idx(const bin_slice_t *s, long g)
{
	if (g == 0)
		return 0;
	if (s->n <= 1 || g < s->first || g > s->first + s->n - 2)
		return -1;
	return g - s->first + 1;
}

/**
* @brief Copies the properties of global binary input row s_k to the index d_k of the binary array.
*
* @param bins this rank's binary input slice
* @param d_k destination index in the binary array where the data needs to be copied to
* @param s_k global binary row (as stored in the object table's binind) to copy from
*/
static void load_binary_data(const bin_slice_t *bins, int d_k, int s_k) {
		long l_k = bin_local_idx(bins, s_k);

		if (l_k < 0) {
			eprintf("binary row %d not in rank %d binary slice [%ld, %ld]\n", s_k, myid, bins->first, bins->first + bins->n - 2);
			exit_cleanly(-1, __FUNCTION__);
		}
		binary[d_k].id1 = bins->id1[l_k];
		binary[d_k].id2 = bins->id2[l_k];
		if (binary[d_k].id1 > newstarid) {
			newstarid = binary[d_k].id1;
		}
		if (binary[d_k].id2 > newstarid) {
			newstarid = binary[d_k].id2;
		}
		binary[d_k].rad1 = bins->Reff1[l_k];
		binary[d_k].rad2 = bins->Reff2[l_k];
		binary[d_k].m1 = bins->m1[l_k] * ((double) clus.N_STAR);
		binary[d_k].m2 = bins->m2[l_k] * ((double) clus.N_STAR);
		binary[d_k].a = bins->a[l_k];
		binary[d_k].e = bins->e[l_k];
		binary[d_k].inuse = 1;
		/*Sourav: assign lifetimes to the binary components*/
		if (STAR_AGING_SCHEME==1 ||STAR_AGING_SCHEME==3){
			if (PREAGING){
				binary[d_k].createtime_m1 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
				binary[d_k].createtime_m2 = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			}else {
				binary[d_k].createtime_m1 = 0.0;
				binary[d_k].createtime_m2 = 0.0;
			}
			binary[d_k].lifetime_m1 = pow(10.0,9.921)*pow((binary[d_k].m1 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			binary[d_k].lifetime_m2 = pow(10.0,9.921)*pow((binary[d_k].m2 * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
		}
		else {
			binary[d_k].createtime_m1 = 0.0;
			binary[d_k].createtime_m2 = 0.0;
			binary[d_k].lifetime_m1 = GSL_POSINF;
			binary[d_k].lifetime_m2 = GSL_POSINF;
		}
}

/* Rank-local copy of the per-object input fields load_fits_file_data needs, indexed by
   LOCAL star index: position 0 is global row 0 (the sentinel), position i >= 1 is the i-th
   row this rank owns, i.e. global row Start[myid] + i - 1 (see get_global_idx). Field types
   match cfd's obj_* arrays so the double -> long/int conversions stay two-step as before. */
typedef struct {
	long    n;
	long    *id;
	int     *k;
	double  *Reff;
	double  *vr;
	double  *vt;
	long    *binind;
} obj_slice_t;

static void obj_slice_alloc(obj_slice_t *s, long n)
{
	s->n      = n;
	s->id     = (long *)   calloc(n, sizeof(long));
	s->k      = (int *)    calloc(n, sizeof(int));
	s->Reff   = (double *) calloc(n, sizeof(double));
	s->vr     = (double *) calloc(n, sizeof(double));
	s->vt     = (double *) calloc(n, sizeof(double));
	s->binind = (long *)   calloc(n, sizeof(long));
	if (!s->id || !s->k || !s->Reff || !s->vr || !s->vt || !s->binind) {
		eprintf("cannot allocate object slice of %ld rows\n", n);
		exit_cleanly(-1, __FUNCTION__);
	}
}

static void obj_slice_free(obj_slice_t *s)
{
	free(s->id); free(s->k); free(s->Reff); free(s->vr); free(s->vt); free(s->binind);
}

/* Local star index of global object row g on this rank, or -1 if another rank owns it.
   Row 0 maps to local 0 on every rank. Inverse of get_global_idx() for owned rows. */
static long obj_local_idx(hsize_t g)
{
	if (g == 0)
		return 0;
	if ((long) g < Start[myid] || (long) g > End[myid])
		return -1;
	return (long) g - Start[myid] + 1;
}

/**
* @brief Block-reads the initial object table once: m and r for every row straight into
* star_m/star_r, and the other fields only for row 0 and the rows this rank owns.
*/
static void load_hdf5_objects(obj_slice_t *objs)
{
	hsize_t nss = (hsize_t) cfd.NOBJ + 2;
	hsize_t dims[2], nrows, ncols, nblock, start, r, offset[2], count[2];
	hid_t   file, dset, space, memspace;
	double  *blk;

	/* Same defaults cmc_malloc_fits_data_t gives obj_m/obj_r; the file overwrites them. */
	star_m[0] = 0.0;
	star_r[0] = DBL_MIN;
	star_m[nss-1] = 0.0;
	star_r[nss-1] = LARGE_DISTANCE;

	file = H5Fopen(INPUT_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
	dset = (file < 0) ? -1 : H5Dopen(file, "CLUS_OBJ_DATA/block0_values", H5P_DEFAULT);
	if (dset < 0) {
		eprintf("cannot open CLUS_OBJ_DATA/block0_values in %s\n", INPUT_FILE);
		exit_cleanly(-1, __FUNCTION__);
	}
	space = H5Dget_space(dset);
	H5Sget_simple_extent_dims(space, dims, NULL);
	ncols = dims[1];
	if (ncols < 8) {
		eprintf("object table has %llu columns, need 8\n", (unsigned long long) ncols);
		exit_cleanly(-1, __FUNCTION__);
	}
	nrows  = (dims[0] < nss) ? dims[0] : nss;
	nblock = (nrows < (hsize_t) CMC_READ_BLOCK_ROWS) ? nrows : (hsize_t) CMC_READ_BLOCK_ROWS;
	if (nblock == 0)
		nblock = 1;

	blk = (double *) malloc((size_t) nblock * (size_t) ncols * sizeof(double));
	if (blk == NULL) {
		eprintf("cannot allocate object read block\n");
		exit_cleanly(-1, __FUNCTION__);
	}

	count[0] = nblock;
	count[1] = ncols;
	memspace = H5Screate_simple(2, count, NULL);

	for (start = 0; start < nrows; start += nblock) {
		hsize_t nthis = (nrows - start < nblock) ? nrows - start : nblock;

		offset[0] = start;
		offset[1] = 0;
		count[0]  = nthis;
		count[1]  = ncols;

		H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
		H5Sset_extent_simple(memspace, 2, count, NULL);
		H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, blk);

		for (r = 0; r < nthis; r++) {
			double  *row = blk + r * ncols;
			hsize_t  g   = start + r;           /* global row */
			long     li  = obj_local_idx(g);    /* local star index, -1 if not ours */

			star_m[g] = row[2] * ((double) clus.N_STAR);
			star_r[g] = row[4];

			if (li >= 0) {
				objs->id[li]     = row[0];
				objs->k[li]      = row[1];
				objs->Reff[li]   = row[3];
				objs->vr[li]     = row[5];
				objs->vt[li]     = row[6];
				objs->binind[li] = row[7];
			}
		}
	}

	free(blk);
	H5Sclose(memspace);
	H5Sclose(space);
	H5Dclose(dset);
	H5Fclose(file);
}

/**
* @brief Reads binary input row 0 and the consecutive rows bins->first.. this rank needs from
* the HDF5 binary table into bins, in bounded blocks.
*/
static void load_hdf5_binaries(bin_slice_t *bins)
{
	hsize_t dims[2], ncols, nblock, maxrows, start, r, offset[2], count[2];
	hsize_t rng_start[2], rng_rows[2];
	hid_t   file, dset, space, memspace;
	double  *blk;
	int     k;

	rng_start[0] = 0;
	rng_rows[0]  = 1;
	rng_start[1] = (hsize_t) bins->first;
	rng_rows[1]  = (hsize_t) (bins->n - 1);

	file = H5Fopen(INPUT_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
	dset = (file < 0) ? -1 : H5Dopen(file, "CLUS_BINARY_DATA/block0_values", H5P_DEFAULT);
	if (dset < 0) {
		eprintf("cannot open CLUS_BINARY_DATA/block0_values in %s\n", INPUT_FILE);
		exit_cleanly(-1, __FUNCTION__);
	}
	space = H5Dget_space(dset);
	H5Sget_simple_extent_dims(space, dims, NULL);
	ncols = dims[1];
	if (ncols < 11 || dims[0] < 1 || (bins->n > 1 && rng_start[1] + rng_rows[1] > dims[0])) {
		eprintf("binary table (%llu x %llu) too small for rows %ld..%ld\n",
			(unsigned long long) dims[0], (unsigned long long) ncols, bins->first, bins->first + bins->n - 2);
		exit_cleanly(-1, __FUNCTION__);
	}

	maxrows = (rng_rows[1] > rng_rows[0]) ? rng_rows[1] : rng_rows[0];
	nblock  = (maxrows < (hsize_t) CMC_READ_BLOCK_ROWS) ? maxrows : (hsize_t) CMC_READ_BLOCK_ROWS;
	blk = (double *) malloc((size_t) nblock * (size_t) ncols * sizeof(double));
	if (blk == NULL) {
		eprintf("cannot allocate binary read block\n");
		exit_cleanly(-1, __FUNCTION__);
	}
	count[0] = nblock;
	count[1] = ncols;
	memspace = H5Screate_simple(2, count, NULL);

	for (k = 0; k < 2; k++) {
		for (start = 0; start < rng_rows[k]; start += nblock) {
			hsize_t nthis = (rng_rows[k] - start < nblock) ? rng_rows[k] - start : nblock;

			offset[0] = rng_start[k] + start;
			offset[1] = 0;
			count[0]  = nthis;
			count[1]  = ncols;

			H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
			H5Sset_extent_simple(memspace, 2, count, NULL);
			H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, blk);

			for (r = 0; r < nthis; r++) {
				double *row = blk + r * ncols;
				long    g   = (long) (offset[0] + r);      /* global binary row */
				long    li  = bin_local_idx(bins, g);      /* local slice row   */

				bins->id1[li]   = row[1];
				bins->m1[li]    = row[3];
				bins->Reff1[li] = row[4];
				bins->id2[li]   = row[5];
				bins->m2[li]    = row[7];
				bins->Reff2[li] = row[8];
				bins->a[li]     = row[9];
				bins->e[li]     = row[10];
			}
		}
	}

	free(blk);
	H5Sclose(memspace);
	H5Sclose(space);
	H5Dclose(dset);
	H5Fclose(file);
}

/**
* @brief Populates the star and binary arrays from the data obtained from the cfd structure. Each processor just reads and stores the corresponding data from the cfd structure it is responsible for, and ignores the rest.
*/
void load_fits_file_data(void)
{
	long i, g_i=0, nown;
	//0th element of the binary holds nothing
	int b_i=1;
	obj_slice_t objs;
	bin_slice_t bins;
	long bmin=0, bmax=0, nbin=0, prev=0, consecutive=1;

	newstarid = 0;

	/* set units */
	units_set();

	/* MPI: this rank owns global object rows Start[myid]..End[myid] (possibly none). */
	nown = End[myid] - Start[myid] + 1;
	if (nown < 0 || (nown > 0 && (get_global_idx(1) != Start[myid] || get_global_idx(nown) != End[myid]))) {
		eprintf("inconsistent object ownership on rank %d: Start=%d End=%d\n", myid, Start[myid], End[myid]);
		exit_cleanly(-1, __FUNCTION__);
	}
	obj_slice_alloc(&objs, nown + 1);

	//MPI: All procs read data for global arrays.
	/* star_m/star_r are filled before star[], which is sourced from them below; the
	   sentinel and central-mass overrides must still run after both. The HDF5 reader
	   leaves cfd's obj_* arrays NULL and the table is read here; the FITS reader still
	   stores full-length cfd arrays. */
	if (cfd.obj_m != NULL) {
		for (i=0; i<=cfd.NOBJ+1; i++) {
			star_m[i] = cfd.obj_m[i] * ((double) clus.N_STAR);
			star_r[i] = cfd.obj_r[i];
		}
		for (i=0; i<=nown; i++) {
			g_i = (i==0) ? 0 : get_global_idx(i);
			objs.id[i]     = cfd.obj_id[g_i];
			objs.k[i]      = cfd.obj_k[g_i];
			objs.Reff[i]   = cfd.obj_Reff[g_i];
			objs.vr[i]     = cfd.obj_vr[g_i];
			objs.vt[i]     = cfd.obj_vt[g_i];
			objs.binind[i] = cfd.obj_binind[g_i];
		}
	} else {
		load_hdf5_objects(&objs);
	}

	/* MPI: the binary rows this rank needs are the nonzero binind values of its objects.
	   For generated inputs they form one consecutive run, so a direct-mapped slice
	   [bmin, bmax] holds exactly them. The slice is correct for any binind set; if the
	   run is not consecutive it is only larger than necessary, and we say so. FITS input
	   already holds the full binary table in cfd, so it is used in place instead. */
	for (i=0; i<=nown; i++) {
		long b = objs.binind[i];
		if (!b)
			continue;
		if (b < 1 || b > cfd.NBINARY) {
			eprintf("object binind %ld out of range 1..%ld on rank %d\n", b, cfd.NBINARY, myid);
			exit_cleanly(-1, __FUNCTION__);
		}
		if (nbin == 0) {
			bmin = bmax = b;
		} else {
			if (b != prev + 1)
				consecutive = 0;
			if (b < bmin) bmin = b;
			if (b > bmax) bmax = b;
		}
		prev = b;
		nbin++;
	}
	if (cfd.bs_id1 != NULL) {	/* FITS reader: full-length cfd binary arrays */
		bin_slice_alias_cfd(&bins);
	} else {
		if (nbin > 0 && !consecutive)
			fprintf(stderr, "WARNING: rank %d: owned binind values are not one consecutive run; reading binary rows %ld..%ld (%ld rows) for %ld binaries\n",
				myid, bmin, bmax, bmax - bmin + 1, nbin);
		bin_slice_alloc(&bins, (nbin > 0) ? bmax - bmin + 2 : 1, (nbin > 0) ? bmin : 0);
		load_hdf5_binaries(&bins);
	}

	//MPI: Just to load the 0th element (never used) with whatever is in the input.
	load_binary_data(&bins, 0, 0);

	//MPI: Copying only the data each process needs.
	for (i=0; i<=nown; i++) {
		//MPI: Getting global index to read only stars that belong to this processor.
		if(i==0) g_i = 0;
		else g_i = get_global_idx(i);
		star[i].id = objs.id[i];
		if (star[i].id > newstarid) {
			newstarid = star[i].id;
		}
		star[i].se_k = objs.k[i];
		star[i].rad = objs.Reff[i];
		star[i].vr = objs.vr[i];
		star[i].vt = objs.vt[i];

		//Copying binary info
      if(objs.binind[i])
      {
         load_binary_data(&bins, b_i, objs.binind[i]);
         star[i].binind = b_i;
         b_i++;
      }

		star[i].m = star_m[g_i];
		star[i].r = star_r[g_i];
		/*Sourav: putting creation time and lifetime as a variable*/
		/*Sourav: this is ongoing changes*/
		if (!star[i].binind){
			if (STAR_AGING_SCHEME==1 ||STAR_AGING_SCHEME==3){
				if(PREAGING){
					star[i].createtime = - pow(10.0,9.921)*pow(PREAGING_MASS,-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
				} else {
					star[i].createtime = 0.0;
				}	
				star[i].lifetime = pow(10.0,9.921)*pow((star[i].m * units.mstar / MSUN),-3.6648)*YEAR*log(GAMMA*clus.N_STAR)/units.t/clus.N_STAR;
			}
			else {
				star[i].createtime = 0.0;
				star[i].lifetime = GSL_POSINF;
			}
		}
	}

	obj_slice_free(&objs);
	bin_slice_free(&bins);

	N_b_local = b_i-1;

	/* some assignments so the code won't break */
	star_r[0] = ZERO; 
	star_r[clus.N_STAR + 1] = SF_INFINITY;
	Mtotal = 1.0;

	// central mass business, read in from file
	// I believe the normalization should be correct, from above
	cenma.m = star_m[0];
	cenma.m_new = 0.0;
	star_m[0] = 0.0;
	//MPI: The way newstarid is assigned above works in the serial version, but it wont be the same on all procs in the parallel version. Assuming it will be the sum of N_STAR and N_BINARY.
	newstarid = clus.N_STAR+clus.N_BINARY;
    //Maia: Initialize new star ID variables
    idblocksize = newstarid/procs; // This way, idblocksize*procs is guaranteed to be strictly equal to or less than newstarid

	if (BH_R_DISRUPT_NB> 0) {
		diaprintf("R_disrupt in NB-units for all stars, Rdisr=%lg\n", BH_R_DISRUPT_NB);
	} else {
		diaprintf("R_disrupt for a solar mass star in NB-units, Rdisr=%lg\n", 
				pow(2.*cenma.m/MSUN*units.mstar, 1./3.)*RSUN/units.l);
	};
}
