#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_sf_erf.h>
#include <fitsio.h>
#include "taus113-v2.h"

#define LARGE_DISTANCE 1.0e40

/* default parameters */
#define W0 6.0
#define NSTAR 100000UL
#define OUTFILE_FORMAT "king_w%.3g_n%.3g.fits"
#define SEED 0UL

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Wed Jul 28 15:43:12 CDT 2004"

int debug=0;

struct strpar {
	double rt;
	double rk;
};

/* print the usage */
void print_usage(FILE *stream)
{
	char outfile[1024];
	
	sprintf(outfile, OUTFILE_FORMAT, W0, (double) NSTAR);
	
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  kingfit [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -w --W0 <W_0 parameter> : set King model W_0 parameter [%.6g]\n", W0);
	fprintf(stream, "  -N --N <N>              : set number of stars [%ld]\n", NSTAR);
	fprintf(stream, "  -o --outfile <outfile>  : set name of outfile [king_w<w0>_n<N>.fits]\n");
	fprintf(stream, "                            where <w0> and <N> are replaced by their respective values\n");
	fprintf(stream, "  -s --seed <seed>        : set random seed [%ld]\n", SEED);
	fprintf(stream, "  -d --debug              : turn on debugging\n");
	fprintf(stream, "  -V --version            : print version info\n");
	fprintf(stream, "  -h --help               : display this help text\n");
}

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** kingfit %s (%s) [%s] **\n", VERSION, NICK, DATE);
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
	double rho0, result;

	rho0 = exp(w0)*gsl_sf_erf(sqrt(w0))-sqrt(4.0*w0/M_PI)*(1.0+2.0*w0/3.0);
	result = -9.0*(exp(x)*gsl_sf_erf(sqrt(x))-sqrt(4.0*x/M_PI)*(1.0+2.0*x/3.0))/rho0;
	if(fabs(y)<DBL_EPSILON && fabs(t)<DBL_EPSILON){
		result /= 3.0;
	} else {
		result -= 2.0*y/t;
	}

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
	double w0=W0;
	unsigned long int N=NSTAR, seed=SEED;
	char filename[1024];
	double tol = 1e-10, hmax = 0.01, hmin = 1e-11, error;
	double t, x, y;
	double t2, t3, t4, t5, t6, x2, x3, x4, x5, x6, y2, y3, y4, y5, y6;
	double h;
	double k1, k2, k3, k4, k5, k6, l1, l2, l3, l4, l5, l6;
	double rho_rho0;
	int i=0, array_incr = 10000, array_size=array_incr;
	int num_point;
	double *pos, *den, *pot, *mR;
	double X1, X2, X3, *X;
	double *r, *m, *vt, *vr, *psi;
	double v_0, f_0, f, F;
	int jmin, jmax, jtry;
	struct strpar struct_par;
	const char *short_opts = "w:N:o:s:dVh";
	const struct option long_opts[] = {
		{"W0", required_argument, NULL, 'w'},
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
		case 'w':
			w0 = atof(optarg);
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
		sprintf(filename, OUTFILE_FORMAT, w0, (double) N);
	}

	check_for_file(filename);
	/* reset_rng_t113(2718234592UL); */
	reset_rng_t113(seed);
	/* t is r/r0, x is Psi/sigma^2, y is dx/dt */
	h = tol;
	t = 0.0; 
	x = w0; y = 0.0; 
	pos = malloc(array_size*sizeof(double));
	den = malloc(array_size*sizeof(double));
	pot = malloc(array_size*sizeof(double));
	while(x>0 && (rho_rho0=calc_rho_rho0(t,x,y,w0))>0){
		pos[i] = t; 
		den[i] = rho_rho0; 
		pot[i] = x;
		i++;
		if(i>=array_size){
			array_size += array_incr;
			pos = realloc(pos, array_size*sizeof(double));
			den = realloc(den, array_size*sizeof(double));
			pot = realloc(pot, array_size*sizeof(double));
		}

		k1 = h*f0(t,x,y,w0);
		if (isnan(k1) && h>hmin){ h*=0.5; continue; }
		l1 = h*f1(t,x,y,w0);
		
		t2 = t + h/4.0; 
		x2 = x + k1/4.0; 
		y2 = y + l1/4.0;
		k2 = h*f0(t2,x2,y2,w0);
		if (isnan(k2) && h>hmin){ h*=0.5; continue; }
		l2 = h*f1(t2,x2,y2,w0);
		
		t3 = t + 3.0/8.0*h; 
		x3 = x + 3.0/32.0*k1 + 9.0/32.0*k2; 
		y3 = y + 3.0/32.0*l1 + 9.0/32.0*l2; 
		k3 = h*f0(t3,x3,y3,w0);
		if (isnan(k3) && h>hmin){ h*=0.5; continue; }
		l3 = h*f1(t3,x3,y3,w0);
		
		t4 = t + 12.0/13.0*h;
		x4 = x + 1932.0/2197*k1 - 7200.0/2197*k2 + 7296.0/2197*k3;
		y4 = y + 1932.0/2197*l1 - 7200.0/2197*l2 + 7296.0/2197*l3;
		k4 = h*f0(t4,x4,y4,w0);
		if (isnan(k4) && h>hmin){ h*=0.5; continue; }
		l4 = h*f1(t4,x4,y4,w0);

		t5 = t + h;
		x5 = x + 439.0/216.0*k1 - 8.0*k2 
		       + 3680.0/513.0*k3 - 845.0/4104*k4;
		y5 = y + 439.0/216.0*l1 - 8.0*l2 
		       + 3680.0/513.0*l3 - 845.0/4104*l4;
		k5 = h*f0(t5,x5,y5,w0);
		if (isnan(k5) && h>hmin){ h*=0.5; continue; }
		l5 = h*f1(t5,x5,y5,w0);


		t6 = t + h/2.0;
		x6 = x - 8.0/27.0*k1 + 2.0*k2 - 3544.0/2565.0*k3
		       + 1859.0/4104.0*k4 - 11.0/40.0*k5;
		y6 = y - 8.0/27.0*l1 + 2.0*l2 - 3544.0/2565.0*l3
		       + 1859.0/4104.0*l4 - 11.0/40.0*l5;
		k6 = h*f0(t6,x6,y6,w0);
		if (isnan(k6) && h>hmin){ h*=0.5; continue; }
		l6 = h*f1(t6,x6,y6,w0);
		
		error = 1.0/360.0*k1 - 128.0/4275.0*k3 - 2197.0/75240.0*k4
			+ 1.0/50.0*k5 + 2.0/55.0*k6;
		if (error>tol){
			h *= 0.5;
			continue;
		}
		
		x += incr(k1, k2, k3, k4, k5, k6);
		y += incr(l1, l2, l3, l4, l5, l6);
		t += h;
		if (error<tol/2.0 && h<hmax/2.0){
			h *= 2.0;
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
