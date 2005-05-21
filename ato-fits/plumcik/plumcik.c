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
#define OUTFILE_FORMAT "plummer_n%.3g.fits"
#define SEED 0UL

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Thu Feb 10 16:58:35 CST 2005"

int debug=0;

double Mrtotal;

#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}

/* print the usage */
void print_usage(FILE *stream)
{
	char outfile[1024];
	
	sprintf(outfile, OUTFILE_FORMAT, (double) NSTAR);
	
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  plumcik [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -r --rmax <r_max>      : set maximum radius [%.6g]\n", RMAX);
	fprintf(stream, "  -N --N <N>             : set number of stars [%ld]\n", NSTAR);
	fprintf(stream, "  -o --outfile <outfile> : set name of outfile [plummer_n<N>.fits]\n");
	fprintf(stream, "                           where <N> is replaced by its value\n");
	fprintf(stream, "  -s --seed <seed>       : set random seed [%ld]\n", SEED);
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

/* smaller cms means less number of stars near the center  */
/* rcr is the radius within which there is a different IMF 
 * it is given in units of plummer radius and converted
 * to Nboody units by dividing out with (3.0*PI/16.0)      */

double calc_M(double r, double rcr, double rmax, double cms){
	double M;
	double rcr2, rcr3, r2, r3, rmax2, rmax3;

	r2 = r*r;
	r3 = r2*r;
	rcr2 = rcr*rcr;
	rcr3 = rcr2*rcr;
	rmax2 = rmax*rmax;
	rmax3 = rmax2*rmax;
	if(r<rcr) {
		M = cms*r3*pow(1.0+r2,-3.0/2.0);
	} else {
		M = (cms-1.0)*rcr3*pow(1.0+rcr2,-3.0/2.0);
		M += r3*pow(1.0+r2,-3.0/2.0);
	}
	M /= Mrtotal;
	return M;
}

double find_r(double X, double rcr, double cms, double r_abs_max, double tol){
	double rmin, rmax, rtry;
	double M;
	int k;
	
	k = 0;
	rmin = 0.0;
	rmax = r_abs_max;
	while(rmax-rmin>tol){
		rtry = (rmax+rmin)/2.0;
		M = calc_M(rtry, rcr, r_abs_max, cms);
		if(M>X){
			rmax = rtry;
		} else {
			rmin = rtry;
		}
		if(k++>80){
			fprintf(stderr,"too many iterations...\n");
			fprintf(stderr,"X = %.10e\n", X);
			break;
		}
	}
	return (rmax+rmin)/2.0;
}

#define GENSORT_NAME                 sort_X
#define GENSORT_TYPE                 double
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            a
#define GENSORT_COMPAREKEYS(k1,k2)   k1 < k2

#include "gensort.h"

void create_random_array(double *X, long int N){
	long int i, j;
	int duplicate;
	double tempX;

	X[0]=0.0; X[N+1]=1.0;
	for(i=1; i<=N; i++){
		X[i] = rng_t113_dbl();
	}
	sort_X(X+1, N);
	/* removing the duplicates from the random number set 
	 * since we use binary search this does NOT guarantee
	 * unique positions 					*/
	do{
		duplicate = 0;
		for(i=1; i<=N-1; i++){
			if(X[i]==X[i+1]) {
				duplicate = 1;
				break;
			}
		}
		if (duplicate){
			for(j=i; j<=N-1; j++){
				X[j] = X[j+1];
			}
			X[N] = rng_t113_dbl();
			for(j=N; j>1; j--){
				if (X[j]<X[j-1]){
					tempX = X[j];
					X[j] = X[j-1];
					X[j-1] = tempX;
				}
			}
		} 
	} while (duplicate);
}

void find_positions(double *X, double *r, long int N, 
		double rcr, double cms, double rmax, double tol){
	long int i;

	r[0] = DBL_MIN;
	for(i=1; i<=N; i++){
		r[i] = find_r(X[i], rcr, cms, rmax, tol);
	}
	r[N+1] = LARGE_DISTANCE;
}

void find_velocities(double *r, double *vr, double *vt, long int N){
	long int i;
	double X1, X2, X3, vesc, v, g, q;
	
	vr[0] = vt[0] = 0.0;
	for(i=1; i<=N; i++){
		vesc = sqrt(2.0/sqrt(1.0+r[i]*r[i])); 
		do {
			X1 = rng_t113_dbl();
			X2 = rng_t113_dbl();
			q = X1;
			g = q*q*pow(1-q*q, 7.0/2.0);
		} while (0.1*X2 > g);
		v = q * vesc;
		X3 = rng_t113_dbl();
		vr[i] = (1.0 - 2.0*X3) * v;
		vt[i] = sqrt(v*v - vr[i]*vr[i]);
	}
	vr[N+1] = vt[N+1] = 0.0;
}

void set_masses(double *m, double *r, long int N, double rcr, double cms){
	long int i;
	double mtot;

	mtot = 0.0;
	m[0] = 0.0;
	for(i=1; i<=N; i++){
		if(r[i]>rcr){
			m[i] = cms;
		} else {
			m[i] = 1.0;
		}
		mtot += m[i];
	}
	m[N+1] = 0.0; 
	for(i=1; i<=N; i++){
		m[i] /= mtot;
	}
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

void scale_pos_and_vel(double *m, double *r, double *vr, 
		double *vt, long int N, double rcr){
	long int i;
	double PEtot, KEtot, U, T;
	double MM, rfac, vfac;
	
	PEtot = KEtot = 0.0;
	U = 0.0;
	MM = 1.0; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
		U -= MM*(1.0/r[i] - 1.0/r[i+1]);
		T = 0.5 * (vr[i]*vr[i] + vt[i]*vt[i]);
		MM -= m[i];
		PEtot += 0.5*U*m[i];
		KEtot += T*m[i];
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
		U -= MM*(1.0/r[i] - 1.0/r[i+1]);
		T = 0.5 * (vr[i]*vr[i] + vt[i]*vt[i]);
		MM -= m[i];
		PEtot += 0.5*U*m[i];
		KEtot += T*m[i];
	}
	printf("After  scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
	printf("value of rcr is: %e\n", rcr*rfac);
}
	
int main(int argc, char *argv[]){
	double *X, *r, *vr, *vt, *m;
	double rmax=RMAX;
	double rcr = 1.0/(6.0*(PI/32.0));
	double cms=1.0;
	unsigned long int N=NSTAR, seed=SEED;
	char filename[1024];
	int i;
	const char *short_opts = "r:N:o:s:dVh";
	const struct option long_opts[] = {
		{"rmax", required_argument, NULL, 'r'},
		{"N", required_argument, NULL, 'N'},
		{"outfile", required_argument, NULL, 'o'},
		{"seed", required_argument, NULL, 's'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};

	filename[0] = 0;	/* setting a meaningless default name */

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'r':
			rmax = atof(optarg);
			break;
		case 'N':
			N = atol(optarg);
			break;
		case 'o':
			sprintf(filename, "%s", optarg);
			break;
		case 's':
			seed = atol(optarg);
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
	if (optind < argc) {
		print_usage(stdout);
		return(1);
	}
	
	if (filename[0]==0){
		sprintf(filename, OUTFILE_FORMAT, (double) N);
	}

	dprintf("rmax=%g N=%ld filename=%s seed=%ld\n", rmax, N, filename, seed);

	X = malloc((N+2)*sizeof(double));
	r = malloc((N+2)*sizeof(double));
	vr = malloc((N+2)*sizeof(double));
	vt = malloc((N+2)*sizeof(double));
	m = malloc((N+2)*sizeof(double));
	reset_rng_t113(seed);

	Mrtotal = (cms-1.0)*rcr*rcr*rcr*pow(1.0+rcr*rcr,-3.0/2.0)
		+ rmax*rmax*rmax*pow(1.0+rmax*rmax,-3.0/2.0);
	check_for_file(filename);
	create_random_array(X, N);
	find_positions(X, r, N, rcr, cms, rmax, 1e-12);
	find_velocities(r, vr, vt, N);
	set_masses(m, r, N, rcr, cms);
	scale_pos_and_vel(m, r, vr, vt, N, rcr);
	write_output_file(m, r, vr, vt, N, filename);

	/*for(i=1; i<=N; i++){
		printf("%.30e %.30e %.30e\n", 
				r[i], X[i], calc_M(r[i], rcr, rmax, cms));
	}*/
	return 0;
}
