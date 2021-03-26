/* -*- linux-c -*- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#ifdef USE_FITS
	#include <fitsio.h>
#endif
#include "fitslib.h"
#include "hdf5.h"

#define LARGE_DISTANCE 1.0e40


#ifdef USE_FITS
void cmc_fits_printerror(int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
	
	if (status) {
		fits_report_error(stderr, status); /* print error report */
		exit(status);    /* terminate the program, returning error status */
	}
	return;
}
#endif

/**
* @brief malloc the input cmc_fits_data_t structure based on the NOBJ and NBINARY parameters of the struct.
*
* @param cfd Structure of type cmc_fits_data_t which needs to be malloc'ed. Please note that the NOBJ and NBINARY parameters must of the input struct must be set first.
*/
void cmc_malloc_fits_data_t(cmc_fits_data_t *cfd){
	long NSS = cfd->NOBJ+2;
	long NBS = cfd->NBINARY+1;
	
	/* single star stuff */
	cfd->obj_id = (long *) malloc(NSS * sizeof(long));
	cfd->obj_k = (int *) malloc(NSS * sizeof(int));
	cfd->obj_m = (double *) malloc(NSS * sizeof(double));
	cfd->obj_Reff = (double *) malloc(NSS * sizeof(double));
	cfd->obj_r = (double *) malloc(NSS * sizeof(double));
	cfd->obj_vr = (double *) malloc(NSS * sizeof(double));
	cfd->obj_vt = (double *) malloc(NSS * sizeof(double));
	cfd->obj_binind = (long *) malloc(NSS * sizeof(long));

	/* zero out sentinel objects */
	cfd->obj_id[0] = 0;
	cfd->obj_k[0] = 0;
	cfd->obj_m[0] = 0.0;
	cfd->obj_Reff[0] = 0.0;
	cfd->obj_r[0] = DBL_MIN;
	cfd->obj_vr[0] = 0.0;
	cfd->obj_vt[0] = 0.0;
	cfd->obj_binind[0] = 0;

	cfd->obj_id[NSS-1] = 0;
	cfd->obj_k[NSS-1] = 0;
	cfd->obj_m[NSS-1] = 0.0;
	cfd->obj_Reff[NSS-1] = 0.0;
	cfd->obj_r[NSS-1] = LARGE_DISTANCE;
	cfd->obj_vr[NSS-1] = 0.0;
	cfd->obj_vt[NSS-1] = 0.0;
	cfd->obj_binind[NSS-1] = 0;

	/* binary star stuff */
	cfd->bs_index = (long *) malloc(NBS * sizeof(long));
	cfd->bs_id1 = (long *) malloc(NBS * sizeof(long));
	cfd->bs_k1 = (int *) malloc(NBS * sizeof(int));
	cfd->bs_m1 = (double *) malloc(NBS * sizeof(double));
	cfd->bs_Reff1 = (double *) malloc(NBS * sizeof(double));
	cfd->bs_id2 = (long *) malloc(NBS * sizeof(long));
	cfd->bs_k2 = (int *) malloc(NBS * sizeof(int));
	cfd->bs_m2 = (double *) malloc(NBS * sizeof(double));
	cfd->bs_Reff2 = (double *) malloc(NBS * sizeof(double));
	cfd->bs_a = (double *) malloc(NBS * sizeof(double));
	cfd->bs_e = (double *) malloc(NBS * sizeof(double));

	/* zero out sentinel binary */
	cfd->bs_index[0] = 0;
	cfd->bs_id1[0] = 0;
	cfd->bs_k1[0] = 0;
	cfd->bs_m1[0] = 0.0;
	cfd->bs_Reff1[0] = 0.0;
	cfd->bs_id2[0] = 0;
	cfd->bs_k2[0] = 0;
	cfd->bs_m2[0] = 0.0;
	cfd->bs_Reff2[0] = 0.0;
	cfd->bs_a[0] = 0.0;
	cfd->bs_e[0] = 0.0;
}

/**
* @brief  frees the given cmc_fits_data_t data structure
*
* @param cfd Structure of type cmc_fits_data_t which needs to be freed.
*/
void cmc_free_fits_data_t(cmc_fits_data_t *cfd){
	free(cfd->obj_id);
	free(cfd->obj_k);
	free(cfd->obj_m);
	free(cfd->obj_Reff);
	free(cfd->obj_r);
	free(cfd->obj_vr);
	free(cfd->obj_vt);
	free(cfd->obj_binind);
	free(cfd->bs_index);
	free(cfd->bs_id1);
	free(cfd->bs_k1);
	free(cfd->bs_m1);
	free(cfd->bs_Reff1);
	free(cfd->bs_id2);
	free(cfd->bs_k2);
	free(cfd->bs_m2);
	free(cfd->bs_Reff2);
	free(cfd->bs_a);
	free(cfd->bs_e);
}

#ifdef USE_FITS
/**
* @brief Writes the given cmc_fits_data_t data into the given fits file
*
* @param cfd Struct of type cmc_fits_data_t containing the input data to be writted
* @param filename fits file to which the data needs to be written to
*/
void cmc_write_fits_file(cmc_fits_data_t *cfd, char *filename){
	fitsfile *fptr;
	int status=0;
	int firstrow=1, firstelem=1;
	int tfields;
	int nrows;
	char extname1[] = "CLUS_OBJ_DATA";
	char *ttype1[] = { "id",   "k",    "m",  "Reff", "r",  "vr", "vt", "binind" };
	char *tform1[] = { "1K",   "1K",   "1D", "1D",   "1D", "1D", "1D", "1K" };
	char *tunit1[] = { "none", "none", "NB", "NB",   "NB", "NB", "NB", "none" };
	char extname2[] = "CLUS_BS_DATA";
	char *ttype2[] = { "index", "id1",  "k1",   "m1", "Reff1", "id2", "k2",    "m2", "Reff2", "a",  "e" };
	char *tform2[] = { "1K",    "1K",   "1K",   "1D", "1D",    "1K",  "1K",    "1D", "1D",    "1D", "1D" };
	char *tunit2[] = { "none",  "none", "none", "NB", "NB",    "none", "none", "NB", "NB",    "NB", "NB" };

	fits_create_file(&fptr, filename, &status);
	cmc_fits_printerror(status);
	fits_open_file(&fptr, filename, READWRITE, &status);
	cmc_fits_printerror(status);

	/* create first real HDU */
	tfields = 8;
	nrows = cfd->NOBJ+2;

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype1, tform1, tunit1, extname1, &status);
	cmc_fits_printerror(status);

	fits_update_key(fptr, TLONG, "NOBJ", &(cfd->NOBJ), "number of objects", &status);
	fits_update_key(fptr, TLONG, "NBINARY", &(cfd->NBINARY), "number of binaries", &status);
	fits_update_key(fptr, TDOUBLE, "MCLUS", &(cfd->Mclus), "cluster mass in Msun", &status);
	fits_update_key(fptr, TDOUBLE, "RVIR", &(cfd->Rvir), "virial radius in parsecs", &status);
	fits_update_key(fptr, TDOUBLE, "RTID", &(cfd->Rtid), "tidal radius in NB units", &status);
	fits_update_key(fptr, TDOUBLE, "Z", &(cfd->Z), "metallicity", &status);
	cmc_fits_printerror(status);

	fits_write_col(fptr, TLONG, 1, firstrow, firstelem, nrows, cfd->obj_id, &status);
	fits_write_col(fptr, TINT, 2, firstrow, firstelem, nrows, cfd->obj_k,  &status);
	fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nrows, cfd->obj_m, &status);
	fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows, cfd->obj_Reff, &status);
	fits_write_col(fptr, TDOUBLE, 5, firstrow, firstelem, nrows, cfd->obj_r, &status);
	fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, nrows, cfd->obj_vr, &status);
	fits_write_col(fptr, TDOUBLE, 7, firstrow, firstelem, nrows, cfd->obj_vt, &status);
	fits_write_col(fptr, TLONG, 8, firstrow, firstelem, nrows, cfd->obj_binind, &status);
	cmc_fits_printerror(status);
	
	/* create second real HDU */
	tfields = 11;
	nrows = cfd->NBINARY+1;

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype2, tform2, tunit2, extname2, &status);
	cmc_fits_printerror(status);

	fits_write_col(fptr, TLONG, 1, firstrow, firstelem, nrows, cfd->bs_index, &status);
	fits_write_col(fptr, TLONG, 2, firstrow, firstelem, nrows, cfd->bs_id1,  &status);
	fits_write_col(fptr, TINT, 3, firstrow, firstelem, nrows, cfd->bs_k1, &status);
	fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows, cfd->bs_m1, &status);
	fits_write_col(fptr, TDOUBLE, 5, firstrow, firstelem, nrows, cfd->bs_Reff1, &status);
	fits_write_col(fptr, TLONG, 6, firstrow, firstelem, nrows, cfd->bs_id2, &status);
	fits_write_col(fptr, TINT, 7, firstrow, firstelem, nrows, cfd->bs_k2, &status);
	fits_write_col(fptr, TDOUBLE, 8, firstrow, firstelem, nrows, cfd->bs_m2, &status);
	fits_write_col(fptr, TDOUBLE, 9, firstrow, firstelem, nrows, cfd->bs_Reff2, &status);
	fits_write_col(fptr, TDOUBLE, 10, firstrow, firstelem, nrows, cfd->bs_a, &status);
	fits_write_col(fptr, TDOUBLE, 11, firstrow, firstelem, nrows, cfd->bs_e, &status);
	cmc_fits_printerror(status);

	fits_close_file(fptr, &status);
	cmc_fits_printerror(status);
}
#endif

/**
* @brief Reads from the given hdf5 file and stores into the corrsponding members of the given cmc_hdf5_data_t data structure
*
* @param cfd Struct of type cmc_hdf5_data_t where the data is stored after reading from the file
* @param filename input fits file that needs to be read
*/
void cmc_read_hdf5_file(char *filename, cmc_fits_data_t *cfd, long RESTART_TCOUNT){
    hid_t       file, filetype, memtype, space, dset, attr1, dset1;
    hsize_t     dims[2] = {1, 1};

    herr_t      status;
    double      **read_stellar_objects;                       /* Read buffer */
    double      **read_binary_objects;                       /* Read buffer */
    int         ndims,
                i, j, k;

    char **read_ext_name = (char **) malloc (dims[0] * sizeof (char *));

    double **read_z = (double **) malloc (dims[0] * sizeof (double *));
    read_z[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_z[i] = read_z[0] + i * dims[1];

    double **read_rtid = (double **) malloc (dims[0] * sizeof (double *));
    read_rtid[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_rtid[i] = read_rtid[0] + i * dims[1];

    double **read_rvir = (double **) malloc (dims[0] * sizeof (double *));
    read_rvir[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_rvir[i] = read_rvir[0] + i * dims[1];


    double **read_mclus = (double **) malloc (dims[0] * sizeof (double *));
    read_mclus[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_mclus[i] = read_mclus[0] + i * dims[1];

    long *read_nobj = (long *) malloc (dims[0] * sizeof (long));

    long *read_nbinary = (long *) malloc (dims[0] * sizeof (long));

    /*
     * Now we simply read back the data and output it to the screen.
     */
    hsize_t dims_attr1[1] = {1};

    /*
     * Open file and dataset using the default properties.
     */
    if(access( filename, F_OK ) != 0 ){
        fprintf(stderr,"ERROR: input file %s does not exist!\n",filename);
        exit(1);
    }

    file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, "CLUS_OBJ_DATA/block0_values", H5P_DEFAULT);
    dset1 = H5Dopen (file, "CLUS_BINARY_DATA/block0_values", H5P_DEFAULT);
    attr1 = H5Aopen (dset, "EXTNAME", H5P_DEFAULT);

    /*
     * Get the datatype and its size.
     */
    filetype = H5Aget_type (attr1);

    /*
     * Get dataspace and allocate memory for read buffer.
     */
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);

    /*
     * Create the memory datatype.
     */
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, H5T_VARIABLE);
    status = H5Tset_cset(memtype, H5T_CSET_UTF8); // Specify UTF8 here

    /*
     * Read the data.
     */
    status = H5Aread (attr1, memtype, read_ext_name);
    status = H5Aclose (attr1);


    /*
     * Read next attribute
     */
    attr1 = H5Aopen (dset, "Z", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_NATIVE_DOUBLE, read_z[0]);
    status = H5Aclose (attr1);

    attr1 = H5Aopen (dset, "NOBJ", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_STD_I64LE, read_nobj);
    status = H5Aclose (attr1);

    attr1 = H5Aopen (dset, "NBINARY", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_STD_I64LE, read_nbinary);
    status = H5Aclose (attr1);

    attr1 = H5Aopen (dset, "RVIR", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_NATIVE_DOUBLE, read_rvir[0]);
    status = H5Aclose (attr1);

    attr1 = H5Aopen (dset, "RTID", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_NATIVE_DOUBLE, read_rtid[0]);
    status = H5Aclose (attr1);

    attr1 = H5Aopen (dset, "MCLUS", H5P_DEFAULT);
    filetype = H5Aget_type (attr1);
    space = H5Aget_space (attr1);
    ndims = H5Sget_simple_extent_dims(space, dims_attr1, NULL);
    status = H5Aread (attr1, H5T_IEEE_F64LE, read_mclus[0]);
    status = H5Aclose (attr1);

    cfd->NOBJ = read_nobj[0];
    cfd->NBINARY = read_nbinary[0];
    cfd->Mclus = read_mclus[0][0];
    cfd->Rvir = read_rvir[0][0];
    cfd->Rtid = read_rtid[0][0];
    cfd->Z = read_z[0][0];

    if(RESTART_TCOUNT != 0){
        free (read_z[0]);
        free (read_z);
        free (read_mclus[0]);
        free (read_mclus);
        free (read_rvir[0]);
        free (read_rvir);
        free (read_rtid[0]);
        free (read_rtid);


        free (read_ext_name);
        free (read_nobj);
        free (read_nbinary);

        status = H5Dclose (dset);
        status = H5Dclose (dset1);
        status = H5Sclose (space);
        status = H5Tclose (memtype);
        status = H5Tclose (filetype);
        status = H5Fclose (file);
        return;
    }

    /* Allocate cfd */
    cmc_malloc_fits_data_t(cfd);

    /*
     * Get dataspace and allocate memory for the read buffer as before.
     */
    space = H5Dget_space (dset);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    read_stellar_objects = (double **) malloc (dims[0] * sizeof (double *));
    read_stellar_objects[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_stellar_objects[i] = read_stellar_objects[0] + i * dims[1];

    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                read_stellar_objects[0]);

    for (i=0; i<dims[0]; i++) {
        cfd->obj_id[i] = read_stellar_objects[i][0];
        cfd->obj_k[i] = read_stellar_objects[i][1];
        cfd->obj_m[i] = read_stellar_objects[i][2];
        cfd->obj_Reff[i] = read_stellar_objects[i][3];
        cfd->obj_r[i] = read_stellar_objects[i][4];
        cfd->obj_vr[i] = read_stellar_objects[i][5];
        cfd->obj_vt[i] = read_stellar_objects[i][6];
        cfd->obj_binind[i] = read_stellar_objects[i][7];
    }

    space = H5Dget_space (dset1);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    read_binary_objects = (double **) malloc (dims[0] * sizeof (double *));
    read_binary_objects[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    for (i=1; i<dims[0]; i++)
        read_binary_objects[i] = read_binary_objects[0] + i * dims[1];

    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      read_binary_objects[0]);

    for (i=0; i<dims[0]; i++) {
        cfd->bs_index[i] = read_binary_objects[i][0];
        cfd->bs_id1[i] = read_binary_objects[i][1];
        cfd->bs_k1[i] = read_binary_objects[i][2];
        cfd->bs_m1[i] = read_binary_objects[i][3];
        cfd->bs_Reff1[i] = read_binary_objects[i][4];
        cfd->bs_id2[i] = read_binary_objects[i][5];
        cfd->bs_k2[i] = read_binary_objects[i][6];
        cfd->bs_m2[i] = read_binary_objects[i][7];
        cfd->bs_Reff2[i] = read_binary_objects[i][8];
        cfd->bs_a[i] = read_binary_objects[i][9];
        cfd->bs_e[i] = read_binary_objects[i][10]; 
    }

    /*
     * Close and release resources.
     */
    free (read_stellar_objects[0]);
    free(read_stellar_objects);
    free (read_binary_objects[0]);
    free(read_binary_objects);
    free (read_z[0]);
    free (read_z);
    free (read_mclus[0]);
    free (read_mclus);
    free (read_rvir[0]);
    free (read_rvir);
    free (read_rtid[0]);
    free (read_rtid);


    free (read_ext_name);
    free (read_nobj);
    free (read_nbinary);

    status = H5Dclose (dset);
    status = H5Dclose (dset1);
    status = H5Sclose (space);
    status = H5Tclose (memtype);
    status = H5Tclose (filetype);
    status = H5Fclose (file);


}

#ifdef USE_FITS
/**
* @brief Reads from the given fits file and stores into the corrsponding members of the given cmc_fits_data_t data structure
*
* @param cfd Struct of type cmc_fits_data_t where the data is stored after reading from the file
* @param filename input fits file that needs to be read
*/
void cmc_read_fits_file(char *filename, cmc_fits_data_t *cfd, long RESTART_TCOUNT){
	int status=0, hdunum, hdutype, anynull;
	fitsfile *fptr;
	long frow=1, felem=1, nelem;
	float floatnull=0.0;

	/* open file for reading */
	fits_open_file(&fptr, filename, READONLY, &status);
	cmc_fits_printerror(status);
	hdunum = 2;
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
	cmc_fits_printerror(status);

	/* read keys and then malloc data structure */
	fits_read_key(fptr, TLONG, "NOBJ", &(cfd->NOBJ), NULL, &status);
	fits_read_key(fptr, TLONG, "NBINARY", &(cfd->NBINARY), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "MCLUS", &(cfd->Mclus), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RVIR", &(cfd->Rvir), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RTID", &(cfd->Rtid), NULL, &status);
	fits_read_key(fptr, TDOUBLE, "Z", &(cfd->Z), NULL, &status);
	cmc_fits_printerror(status);

	/*if we're restarting from a checkpoint, we don't need to allocate or read
	 * in the original star files from the fits file...just the header info*/
	if(RESTART_TCOUNT != 0){
		fits_close_file(fptr, &status);
		return;
	}

	cmc_malloc_fits_data_t(cfd);

	/* read in data columns */
	nelem = cfd->NOBJ+2;
	fits_read_col(fptr, TLONG, 1, frow, felem, nelem, &floatnull, cfd->obj_id, &anynull, &status);
	fits_read_col(fptr, TINT, 2, frow, felem, nelem, &floatnull, cfd->obj_k, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 3, frow, felem, nelem, &floatnull, cfd->obj_m, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, &floatnull, cfd->obj_Reff, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 5, frow, felem, nelem, &floatnull, cfd->obj_r, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 6, frow, felem, nelem, &floatnull, cfd->obj_vr, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 7, frow, felem, nelem, &floatnull, cfd->obj_vt, &anynull, &status);
	fits_read_col(fptr, TLONG, 8, frow, felem, nelem, &floatnull, cfd->obj_binind, &anynull, &status);
	cmc_fits_printerror(status);
	
	hdunum = 3;
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
	cmc_fits_printerror(status);

	/* read in data columns */
	nelem = cfd->NBINARY+1;
	fits_read_col(fptr, TLONG, 1, frow, felem, nelem, &floatnull, cfd->bs_index, &anynull, &status);
	fits_read_col(fptr, TLONG, 2, frow, felem, nelem, &floatnull, cfd->bs_id1, &anynull, &status);
	fits_read_col(fptr, TINT, 3, frow, felem, nelem, &floatnull, cfd->bs_k1, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, &floatnull, cfd->bs_m1, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 5, frow, felem, nelem, &floatnull, cfd->bs_Reff1, &anynull, &status);
	fits_read_col(fptr, TLONG, 6, frow, felem, nelem, &floatnull, cfd->bs_id2, &anynull, &status);
	fits_read_col(fptr, TINT, 7, frow, felem, nelem, &floatnull, cfd->bs_k2, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 8, frow, felem, nelem, &floatnull, cfd->bs_m2, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 9, frow, felem, nelem, &floatnull, cfd->bs_Reff2, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 10, frow, felem, nelem, &floatnull, cfd->bs_a, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 11, frow, felem, nelem, &floatnull, cfd->bs_e, &anynull, &status);
	cmc_fits_printerror(status);
	
	fits_close_file(fptr, &status);
	cmc_fits_printerror(status);
}
#endif
