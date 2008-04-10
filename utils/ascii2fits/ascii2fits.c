/* -*- linux-c -*- */
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

/* default parameters */
#define COLMASS 1
#define COLRADIUS 2
#define COLVR 3
#define COLVT 4

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Mon Feb 14 17:28:33 CST 2005"

int debug=0;

#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  ascii2fits [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -i --infile <infile>   : ASCII input file [no default]\n");
	fprintf(stream, "  -o --outfile <outfile> : FITS output file [no default]\n");
	fprintf(stream, "  -m --mass <column #>   : column number of input file containing mass [%d]\n", COLMASS);
	fprintf(stream, "  -R --radius <column #> : column number of input file containing radial position [%d]\n", COLRADIUS);
	fprintf(stream, "  -r --vr <column #>     : column number of input file containing radial velocity [%d]\n", COLVR);
	fprintf(stream, "  -t --vt <column #>     : column number of input file containing tangential velocity [%d]\n", COLVT);
	fprintf(stream, "  -d --debug             : turn on debugging\n");
	fprintf(stream, "  -V --version           : print version info\n");
	fprintf(stream, "  -h --help              : display this help text\n");
}

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** ascii2fits %s (%s) [%s] **\n", VERSION, NICK, DATE);
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

double f0(double t, double x, double y, double w0){
	return y;
}

double f1(double t, double x, double y, double w0){

	return result;
}

double calc_rho_rho0(double t, double x, double y, double w0){
	double rho, rho0;

	rho  = exp(x)*gsl_sf_erf(sqrt(x))-sqrt(4.0*x/M_PI)*(1.0+2.0*x/3.0);
	rho0 = exp(w0)*gsl_sf_erf(sqrt(w0))-sqrt(4.0*w0/M_PI)*(1.0+2.0*w0/3.0);

	return (rho/rho0);
}

double incr(double k1, double k2, double k3, double k4, double k5, double k6){
       	return (25.0/216.0*k1 + 1408.0/2565.0*k3
		+ 2197.0/4104.0*k4 - 1.0/5.0*k5);
}

void set_masses(double *m, long int N){
	long int i;
	double mtot;

	mtot = 0.0;
	m[0] = 0.0;
	for(i=1; i<=N; i++){
		m[i] = 1.0/N;
	}
	m[N+1] = 0.0; 
}

void scale_pos_and_vel(double *m, double *r, double *vr, 
			double *vt, long int N, struct strpar *struct_par){
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
	struct_par->rk = rfac;
	struct_par->rt *= rfac;
	printf("value of r_0 is: %e\n", struct_par->rk);
	printf("value of r_t is: %e\n", struct_par->rt);
}

void write_output_file(double *m, double *r, double *vr, double *vt, 
		long int N, char *filename, struct strpar struct_par){

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
	fits_update_key(fptr, TDOUBLE, "rking", &(struct_par.rk), 
			"King Radius", &status);
	fits_update_key(fptr, TDOUBLE, "rtidal", &(struct_par.rt), 
			"Tidal Radius", &status);

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

void check_for_file(char *filename){
	if(access(filename, F_OK)==0){
		if( strncmp(filename, "debug", 6)==0 ){
			system("rm debug");
		} else if (strncmp(filename, "debug.fit", 9)==0 ){
			system("rm debug.fit");
		} else if (strncmp(filename, "debug.fits", 9)==0 ){
			system("rm debug.fits");
		} else {
			printf("the given filename, %s, exists!\n", filename);
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char *argv[]){
	int colmass=COLMASS, colrad=COLRAD, colvr=COLVR, colvt=COLVT;
	unsigned long int N=NSTAR, seed=SEED;
	char infilename[1024], outfilename[1024], line[16384]="";
	double tol = 1e-10, hmax = 0.01, hmin = 1e-11, error;
	double t, x, y;
	double t2, t3, t4, t5, t6, x2, x3, x4, x5, x6, y2, y3, y4, y5, y6;
	double h;
	double k1, k2, k3, k4, k5, k6, l1, l2, l3, l4, l5, l6;
	double rho_rho0;
	int i, array_incr = 10000, array_size=array_incr;
	int num_point;
	double *pos, *den, *pot, *mR;
	double X1, X2, X3, *X;
	double *r, *m, *vt, *vr, *psi;
	double v_0, f_0, f, F;
	int jmin, jmax, jtry;
	FILE *ifp;
	struct strpar struct_par;
	const char *short_opts = "i:o:m:R:r:t:dVh";
	const struct option long_opts[] = {
		{"infile", required_argument, NULL, 'i'},
		{"outfile", required_argument, NULL, 'o'},
		{"mass", required_argument, NULL, 'm'},
		{"radius", required_argument, NULL, 'R'},
		{"vr", required_argument, NULL, 'r'},
		{"vt", required_argument, NULL, 't'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* setting a meaningless default name */
	infilename[0] = 0;
	outfilename[0] = 0;

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'i':
			sprintf(infilename, "%s", optarg);
			break;
		case 'o':
			sprintf(outfilename, "%s", optarg);
			break;
		case 'm':
			colmass = atol(optarg);
			break;
		case 'R':
			colrad = atol(optarg);
			break;
		case 'r':
			colvr = atol(optarg);
			break;
		case 't':
			colvt = atol(optarg);
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
	if (optind < argc || infilename[0]==0 || outfilename[0] == 0) {
		print_usage(stdout);
		return(1);
	}
	
	dprintf("infilename=%s outfilename=%s colmass=%d colrad=%d colvr=%d colvt=%d\n", 
		infilename, outfilename, colmass, colrad, colvr, colvt);

	check_for_file(outfilename);
	
	/* reset RNG to something */
	reset_rng_t113(2718234592UL);

	ifp = fopen(infilename, "r");
	
	while ((ret = fgets(line, 16384, ifp)) != NULL) {
		/* the maximum size of an input line is limited by 
		 * the size of the line array */
		if (strlen(line) == 16383) {
			fprintf(stderr, "input line too long: \"%s\".\n", line);
			exit(1);
		}
		
		/* strip comments and newline character */
		line[strcspn(line, "#\n")] = '\0';
		
		
		
		if (line[0] != '#' && ) {
			
		}
	}

	num_point = i;
	struct_par.rt = pos[num_point-1];
	mR = malloc(num_point*sizeof(double));
	mR[0] = 0.0;
	for(i=1; i<num_point; i++){
		mR[i] = mR[i-1] + (pos[i]-pos[i-1])
			*(den[i]*pos[i]*pos[i]+den[i-1]*pos[i-1]*pos[i-1])/2.0;
	}
	for(i=0; i<num_point; i++){
		mR[i] /= mR[num_point-1];
	}

	r = malloc((N+2)*sizeof(double));
	m = malloc((N+2)*sizeof(double));
	psi = malloc((N+2)*sizeof(double));
	vr = malloc((N+2)*sizeof(double));
	vt = malloc((N+2)*sizeof(double));
	X = malloc((N+2)*sizeof(double));
	create_random_array(X, N);
	r[0] = DBL_MIN;
	for(i=1; i<=N; i++){
		/* XXX below is uniformly spaced points, use RNG for
		 * a more random distribution distribution in r XXX */
		/* changed */
		jmin = 0; jmax = num_point;
		while(jmin != jmax){ 
			/* loop invariant is: mR[jmin]<X[i]<mR[jmax+1] */
			jtry = (jmin+jmax+1)/2;
			if(mR[jtry]>X[i]){
				jmax = jtry-1;
			} else {
				jmin = jtry;
			}
		}
		if(mR[jmin]>X[i] || mR[jmin+1]<=X[i]){ /* binary search failed! */
			printf("binary search failed! ");
			printf("i = %d, jmin = %d, ",i, jmin);
			printf("mR[jmin]=%e, mR[jmin+1]=%e, X[i] = %e\n",
			 		mR[jmin], mR[jmin+1], X[i]);
			exit(EXIT_FAILURE);
		}
		r[i] = pos[jmin] + (pos[jmin+1]-pos[jmin])*(X[i]-mR[jmin])
			/(mR[jmin+1]-mR[jmin]);
		psi[i] = pot[jmin] + (pot[jmin+1]-pot[jmin])*(X[i]-mR[jmin])
			/(mR[jmin+1]-mR[jmin]);
		/* rho[i] =  den[jmin] + (den[jmin+1]-den[jmin])*(X[i]-mR[jmin]) */
		/*	/(mR[jmin+1]-mR[jmin]); */
	}
	r[N+1] = LARGE_DISTANCE;

	vr[0]=vt[0]=0.0;
	for(i=1; i<=N; i++){
		F = (exp(psi[i])-1.0)*2.0*psi[i];
		do {
			X1 = rng_t113_dbl();
			X2 = rng_t113_dbl();
			v_0 = X1*sqrt(2.0*psi[i]);
			f_0 = X2*F;
			f = (exp(psi[i]-v_0*v_0/2.0)-1.0)*v_0*v_0;
		} while (f_0>f);
		X3 = rng_t113_dbl();
		vr[i] = (1.0 - 2.0*X3) * v_0;
		vt[i] = sqrt(v_0*v_0 - vr[i]*vr[i]);
		/* v[i] = v_0; */
		/* printing data relevant to fig 4-9a in Benn-Trim */
		/* printf("%e %e\n", r[i], rho[i]); */
		/* printing data relevant to fig 4-11 in Benn-Trim */
		/* printf("%e %e\n", r[i], v[i]/sqrt(3.0)); */
		/* printf("%e %e\n", r[i], v_0/sqrt(3.0)); */
	}
	vr[N+1]=vt[N+1]=0.0;
	
	set_masses(m, N);
	scale_pos_and_vel(m, r, vr, vt, N, &struct_par);
	write_output_file(m, r, vr, vt, N, filename, struct_par);
	
	free(r); free(m); free(psi); free(vr); free(vt); free(X);

	return 0;
}
