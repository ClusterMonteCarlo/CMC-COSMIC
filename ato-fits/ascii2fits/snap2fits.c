#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <fitsio.h>
#include "taus113-v2.h"

#define LARGE_DISTANCE 1.0e40
#define PI 3.14159265358979323

/* default parameters */
#define RMAX 300.0
#define NSTAR 100000UL
#define OUTFILE_FORMAT "%s.fits"
#define SEED 0UL

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Thu Feb 10 16:58:35 CST 2005"

int debug=0, mtot_given=0, scale_units=1;

double Mrtotal, Mtot;

#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}

/* print the usage */
void print_usage(FILE *stream)
{
	char outfile[1024];
	
	sprintf(outfile, OUTFILE_FORMAT, "<snapshot file>");
	
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  snap2fits [options...] <snapshot file>\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -o --outfile <outfile> : set name of outfile [<snapshot file>.fits]\n");
        fprintf(stream, "  -m --total-mass        : the total mass of the cluster [calculated from snapshot]\n");
        fprintf(stream, "  -s --no-scaling        : turn off scaling of positions and velocities\n");
	fprintf(stream, "  -d --debug             : turn on debugging\n");
	fprintf(stream, "  -V --version           : print version info\n");
	fprintf(stream, "  -h --help              : display this help text\n");
}

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** plumcik %s (%s) [%s] **\n", VERSION, NICK, DATE);
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

void write_output_file(double *m, double *r, double *vr, double *vt, 
		long int N, char *filename){

	struct rng_t113_state rng_st;

	fitsfile *fptr;
	int status;
	long firstrow, firstelem;
	int tfields;       /* table will have n columns */
	long nrows;	

	char extname[] = "CLUSTER_STARS";          /* extension name */
	char *ttype[] = { "Mass",  "Position", "vr",    "vt" };
	char *tform[] = { "1D",    "1D",       "1D",    "1D" };
	char *tunit[] = { "Nbody", "Nbody",    "Nbody", "Nbody" };
	
	/* these go to header */
	int tstep = 0;
	double time = 0.0;
	
	status = 0;
	tfields = 4;
	nrows = N+2;

	get_rng_t113(&rng_st);

	fits_create_file(&fptr, filename, &status);
	fits_open_file(&fptr, filename, READWRITE, &status);

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status);
	fits_update_key(fptr, TLONG, "NSTAR", &N, 
			"No of Stars", &status);
	fits_update_key(fptr, TDOUBLE, "Time", &time, 
			"Age of cluster", &status);
	fits_update_key(fptr, TLONG, "Step", &tstep, 
			"Iteration Step", &status);
	fits_update_key(fptr, TULONG, "RNG_Z1", &(rng_st.z1), 
			"RNG STATE Z1", &status);
	fits_update_key(fptr, TULONG, "RNG_Z2", &(rng_st.z2), 
			"RNG STATE Z2", &status);
	fits_update_key(fptr, TULONG, "RNG_Z3", &(rng_st.z3), 
			"RNG STATE Z3", &status);
	fits_update_key(fptr, TULONG, "RNG_Z4", &(rng_st.z4), 
			"RNG STATE Z4", &status);

	firstrow  = 1;  /* first row in table to write   */
	firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

	fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, m,
                   &status);
	fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, r,
                   &status);
	fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nrows, vr,
                   &status);
	fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nrows, vt,
                   &status);

	fits_close_file(fptr, &status);
	printerror(status);
}

void check_for_file(char *filename){
	if(access(filename, F_OK)==0){
		if( strncmp(filename, "debug", 6)==0 ){
			system("rm debug");
		} else if (strncmp(filename, "debug.fit", 9)==0 ){
			system("rm debug.fit");
		} else {
			printf("the given filename, %s, exists!\n", filename);
			exit(EXIT_FAILURE);
		}
	}
}

#define GENSORT_NAME                 ind_sort
#define GENSORT_ARGS                 ,r
#define GENSORT_ARGSPROTO            ,double *r
#define GENSORT_TYPE                 long
#define GENSORT_KEYTYPE              double 
#define GENSORT_GETKEY(a)            r[a]
#define GENSORT_COMPAREKEYS(k1,k2)   k1<k2

#define GENSORT_NISINT

#include "gensort.h"

void scale_pos_and_vel(double *m, double *r, double *vr, 
		double *vt, long int N, double rcr){
	long int i, *ind;
	double PEtot, KEtot, U, T;
	double MM, rfac, vfac;

        ind= (long int *) calloc(N+2, sizeof(long));
        for (i=0; i<N+2; i++) {
          ind[i]= i;
        };

        ind_sort(ind, N+1, r);

	PEtot = KEtot = 0.0;
	U = 0.0;
	MM = 1.0; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
                dprintf("At index %li %li\n", ind[i], ind[i+1])
		U -= MM*(1.0/r[ind[i]] - 1.0/r[ind[i+1]]);
		T = 0.5 * (vr[ind[i]]*vr[ind[i]] + vt[ind[i]]*vt[ind[i]]);
		MM -= m[ind[i]]/Mtot;
		PEtot += 0.5*U*m[ind[i]]/Mtot;
		KEtot += T*m[ind[i]]/Mtot;
	}
	printf("Before scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
	/* scaling position and velocity */
	rfac = -PEtot*2.0;
	vfac = 1.0/sqrt(4.0*KEtot);
	for(i=1; i<=N; i++){
		r[i] *= rfac;
		vr[i] *= vfac;
		vt[i] *= vfac;
	}

	PEtot = KEtot = 0.0;
	U = 0.0;
	MM = 1.0; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
		U -= MM*(1.0/r[ind[i]] - 1.0/r[ind[i+1]]);
		T = 0.5 * (vr[ind[i]]*vr[ind[i]] + vt[ind[i]]*vt[ind[i]]);
		MM -= m[ind[i]]/Mtot;
		PEtot += 0.5*U*m[ind[i]]/Mtot;
		KEtot += T*m[ind[i]]/Mtot;
	}
	printf("After  scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
	printf("value of rcr is: %e\n", rcr*rfac);
}

void readSnapshot(const char *filename, double **r, double **vr, double **vt, double **m, unsigned long *N) {
  FILE *inputfile;
  int items;
  long i, block, lines;
  char *line;
  size_t n;

  line= NULL;
  block= 1000;
  lines= 0;
  n= 1;
  items= 1;

  *r= (double *) malloc(sizeof(double));
  *vr= (double *) malloc(sizeof(double));
  *vt= (double *) malloc(sizeof(double));
  *m= (double *) malloc(sizeof(double));
  **r= DBL_MIN;
  **vr= **vt= 0;
  **m= 0;

  inputfile= fopen(filename, "r");
  /* discard the first two lines */
  for (i=0; i<2; i++) {
    getline(&line, &n, inputfile);
    dprintf("The line skipped is %s.\n", line);
  };
  
  i=1;
  do {
    lines+= block;
    *r= (double *) realloc(*r, sizeof(double)*lines);
    *vr= (double *) realloc(*vr, sizeof(double)*lines);
    *vt= (double *) realloc(*vt, sizeof(double)*lines);
    *m= (double *) realloc(*m, sizeof(double)*lines);
    for (; items!=EOF && i<lines-1; i++) {
      items=fscanf(inputfile, 
          "%*i %lf %lf %lf %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f\n",
          &(*m)[i], &(*r)[i], &(*vr)[i], &(*vt)[i]);
    };
  } while (items!=EOF);

  *N= i-2;
  (*r)[*N+1]= LARGE_DISTANCE;
  (*vr)[*N+1]= (*vt)[*N+1]= 0.;
  (*m)[*N+1]= 0.;
}

int main(int argc, char *argv[]){
	double *r, *vr, *vt, *m;
	double rmax=RMAX;
	double rcr = 1.0/(6.0*(PI/32.0));
	unsigned long int N=NSTAR, seed=SEED;
	char filename[1024]; 
        char *input;
	int i;
	const char *short_opts = "so:m:dVh";
	const struct option long_opts[] = {
		{"outfile", required_argument, NULL, 'o'},
                {"total-mass", required_argument, NULL, 'm'},
                {"no-scaling", required_argument, NULL, 's'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};

        r= vr= vt= m= NULL;

	filename[0] = 0;	/* setting a meaningless default name */

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'o':
			sprintf(filename, "%s", optarg);
			break;
                case 'm':
                        mtot_given=1;
                        Mtot= strtod(optarg, NULL);
                        break;
                case 's':
                        scale_units=0;
                        break;
		case 'd':
			debug = 1;
			break;
		case 'V':
			print_version(stdout);
			return(0);
		case 'h':
			print_version(stdout);
			fprintf(stdout, "\n");
			print_usage(stdout);
			return(0);
		default:
			break;
		}
	}
	
	/* check to make sure there was nothing crazy on the command line */
	if (optind >= argc) {
		print_usage(stdout);
		return(EXIT_FAILURE);
	}

        input= argv[optind];
	
	if (filename[0]==0){
		sprintf(filename, OUTFILE_FORMAT, input);
	}

	dprintf("rmax=%g N=%ld filename=%s seed=%ld\n", rmax, N, filename, seed);

        if (!mtot_given) {
          for (i=1; i<N+1; i++) {
            Mtot+= m[i];
          };
        };
        dprintf("The total mass is %lf\n", Mtot);

	check_for_file(filename);
        readSnapshot(input, &r, &vr, &vt, &m, &N);
        if (scale_units) {
          scale_pos_and_vel(m, r, vr, vt, N, rcr);
        };

        /* scale the masses such that Mtot= 1 (N-body unit)*/
        for (i=1; i<N+1; i++) {
          m[i]/= Mtot;
        };

        write_output_file(m, r, vr, vt, N, filename);

	/*for(i=1; i<=N; i++){
		printf("%.30e %.30e %.30e\n", 
				r[i], X[i], calc_M(r[i], rcr, rmax, cms));
	}*/
	return 0;
}
