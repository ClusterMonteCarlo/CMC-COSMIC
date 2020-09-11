/* -*- linux-c -*- */

#ifndef _FITSLIB_H
#define _FITSLIB_H 1

typedef struct{
	long NOBJ;
	long NBINARY;
	double Mclus;
	double Rvir;
	double Rtid;
	double Z;
	/* object info */
	long *obj_id;
	int *obj_k;
	double *obj_m;
	double *obj_Reff;
	double *obj_r;
	double *obj_vr;
	double *obj_vt;
	long *obj_binind;
	/* binary info */
	long *bs_index;
	long *bs_id1;
	int *bs_k1;
	double *bs_m1;
	double *bs_Reff1;
	long *bs_id2;
	int *bs_k2;
	double *bs_m2;
	double *bs_Reff2;
	double *bs_a;
	double *bs_e;
} cmc_fits_data_t;

/* error handling */
void cmc_fits_printerror(int status);
/* malloc big data structure; NOBJ and NBINARY must be set first */
void cmc_malloc_fits_data_t(cmc_fits_data_t *cfd);
/* free big data structure */
void cmc_free_fits_data_t(cmc_fits_data_t *cfd);
/* read in FITS file and assign to data structure */
void cmc_read_fits_file(char *filename, cmc_fits_data_t *cfd, long RESTART_TCOUNT);
/* read in HDF5 file and assign to data structure */
void cmc_read_hdf5_file(char *filename, cmc_fits_data_t *cfd, long RESTART_TCOUNT);
/* write FITS file using data structure */
void cmc_write_fits_file(cmc_fits_data_t *cfd, char *filename);

#endif /* _FITSLIB_H */
