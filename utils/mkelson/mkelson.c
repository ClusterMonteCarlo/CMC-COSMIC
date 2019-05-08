#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include "../../common/fitslib.h"
#include "../../common/taus113-v2.h"
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define LARGE_DISTANCE 1.0e40
#define PI 3.14159265358979323

/* default parameters */
#define RMAX 300.0
#define NSTAR 100000UL
#define GAMMA 3.
#define OUTFILE_FORMAT "elson_n%.3g.fits"
#define SEED 0UL

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Thu Feb 10 16:58:35 CST 2005"

int debug=0;

#define dprintf(args...) if (debug) {fprintf(stderr, "DEBUG: %s(): ", __FUNCTION__); fprintf(stderr, args);}

/* the gsl version of the 2F1 hypergeometric function doesn't include the 
 * analytic continuation to X<-1 for the real numbers.  Implement it here (adapted
 * from the scipy special functions) */
double extended_hyperg_2F1(double a, double b, double c, double x){

	//fprintf(stderr,"a=%g b=%g c=%g x=%g\n",a,b,c,x);

	double p,q,t1,s,y;
	double EPS=1.0e-13;

	s = 1. - x;	
	t1 = fabs(b - a);

	if (x < -1.0) {
		if (fabs(a) < fabs(b)) {
			//fprintf(stderr,"%g %g\n",pow(s, -a), extended_hyperg_2F1(a, c - b, c, x / (x - 1)));
			return pow(s, -a) * extended_hyperg_2F1(a, c - b, c, x / (x - 1));
		}
		else {
			//fprintf(stderr,"%g %g\n",pow(s, -b), extended_hyperg_2F1(b, c - a, c, x / (x - 1)));
			return pow(s, -b) * extended_hyperg_2F1(b, c - a, c, x / (x - 1));
		}
	} else {
		return gsl_sf_hyperg_2F1(a,b,c,x);
	}
}

/* print the usage */
void print_usage(FILE *stream)
{
	char outfile[1024];
	
	sprintf(outfile, OUTFILE_FORMAT, (double) NSTAR);
	
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  mkelson [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -r --rmax <r_max>      : set maximum radius [%.6g]\n", RMAX);
	fprintf(stream, "  -N --N <N>             : set number of stars [%ld]\n", NSTAR);
	fprintf(stream, "  -g --gamma <g>         : set gamma for Elson profile [%.6g]\n", GAMMA);
	fprintf(stream, "  -o --outfile <outfile> : set name of outfile [elson_n<N>.fits]\n");
	fprintf(stream, "                           where <N> is replaced by its value\n");
	fprintf(stream, "  -s --seed <seed>       : set random seed [%ld]\n", SEED);
	fprintf(stream, "  -d --debug             : turn on debugging\n");
	fprintf(stream, "  -V --version           : print version info\n");
	fprintf(stream, "  -h --help              : display this help text\n");
}

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** mkelson %s (%s) [%s] **\n", VERSION, NICK, DATE);
}

/* Compute the mass enclosed in an Elson profile at radius r with slope 
 * gamma, central concentration rho_0, and assumed scale factor a = 1
 *
 * In practice, this is only used to sample the positions of stars, so rho_0 is 
 * just picked to normalize the distribution (i.e. rho_0 s.t. M_enclosed(rmax) = 1) */
double M_enclosed(double r, double gamma, double rho_0){

	double prefactor, hypergeometric;

	prefactor = 4*PI*rho_0*r*r*r / 3.;

	hypergeometric = extended_hyperg_2F1(1.5, (gamma+1.)/2., 2.5, -r*r);

	return prefactor*hypergeometric;
}

/* Compute the density of the Elson profile at radius r 
 * Best to use the same normalized rho_0 from M_enclosed */
double rho_r(double r, double gamma, double rho_0){
	return rho_0*pow(1+r*r,-(gamma+1.)/2.);
}


/* Setup for integrating the Jeans equations to get sigma */
struct jeans_params {double gamma; double rho_0;};

double jeans_equation(double r, void *params){
	struct jeans_params *p = (struct jeans_params *)params;
	double gamma = p->gamma;
	double rho_0 = p->rho_0;

	return rho_r(r,gamma,rho_0)*M_enclosed(r,gamma,rho_0)/r/r; 
}


/* Find the 1D velocity dispersion at a given radius r using one of the 
 * spherial Jeans equations (and assuming velocity is isotropic) */
double find_sigma2(double r, double r_max_cluster, double gamma){
	double integral, error;
	double rho_0;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;

	rho_0 = 1. / M_enclosed(r_max_cluster,gamma,1);
	struct jeans_params parameters = {gamma,rho_0};

	F.function = &jeans_equation;
	F.params = &parameters;

	gsl_integration_qags(&F,r,r_max_cluster,1e-7,1e-8,1000, w,&integral, &error);

	gsl_integration_workspace_free(w);

	return integral / rho_r(r,gamma,rho_0);

}

/* Sample from the cumulative mass profile of an Elson profile; pick a random 
 * number X, and convert that to a radius using the cumulative mass distribution */
double find_r(double X, double r_max_cluster, double gamma, double tol){
	double rmin, rmax, rtry;
	double M;
	double rho_0;
	int k;

	rho_0 = 1. / M_enclosed(r_max_cluster,gamma,1);
	
	k = 0;
	rmin = 0.0;
	rmax = r_max_cluster;
	while(rmax-rmin>tol){
		rtry = (rmax+rmin)/2.0;
		M = M_enclosed(rtry, gamma, rho_0);
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

#include "../../common/gensort.h"

void create_random_array(double *X, long int N){
	long int i;
	int duplicate;

	X[0]=0.0; X[N+1]=1.0;
	for(i=1; i<=N; i++){
		X[i] = rng_t113_dbl();
	}
	/* removing the duplicates from the random number set 
	 * since we use binary search this does NOT guarantee
	 * unique positions 					*/
	do{
		sort_X(X+1, N);
		duplicate = 0;
		for(i=1; i<=N-1; i++){
			if(X[i]==X[i+1]){
				duplicate++;
				X[i] = rng_t113_dbl();
			}
		}
	} while (duplicate);
}

void find_positions(double *X, double *r, long int N, 
		double gamma, double rmax, double tol){
	long int i;

	r[0] = DBL_MIN;
	for(i=1; i<=N; i++){
		r[i] = find_r(X[i], rmax, gamma, tol);
	}
	r[N+1] = LARGE_DISTANCE;
}

void find_velocities(double *r, double *vr, double *vt, long int N, double gamma, double r_max_cluster){
	long int i;
	double X1, X2, X3, X4, sigma2, vx, vy, vz;
	
	vr[0] = vt[0] = 0.0;
	for(i=1; i<=N; i++){
		X1 = rng_t113_dbl();
		X2 = rng_t113_dbl();
		X3 = rng_t113_dbl();
		X4 = rng_t113_dbl();

		// Get the 1D velocity dispersion from the Elson profile
		sigma2 = find_sigma2(r[i],r_max_cluster,gamma);	

		// Use a Box-Muller transform to get Gaussian random numbers
		vx = sqrt(sigma2)*sqrt(-2*log(X1))*cos(2*PI*X2);
		vy = sqrt(sigma2)*sqrt(-2*log(X1))*sin(2*PI*X2);
		vz = sqrt(sigma2)*sqrt(-2*log(X3))*cos(2*PI*X4);

		vr[i] = vx; 
		vt[i] = sqrt(vy*vy +vz*vz);
	}
	vr[N+1] = vt[N+1] = 0.0;
}

void set_masses(double *m, double *r, long int N ){
	long int i;
	double mtot;

	mtot = 0.0;
	m[0] = 0.0;
	for(i=1; i<=N; i++){
		m[i] = 1.0;
		mtot += m[i];
	}
	m[N+1] = 0.0; 
	for(i=1; i<=N; i++){
		m[i] /= mtot;
	}
}

void write_output_file(double *m, double *r, double *vr, double *vt, long int N, char *filename){
	long i;
	cmc_fits_data_t cfd;
	
	cfd.NOBJ = N;
	cfd.NBINARY = 0;
	cfd.Mclus = N;
	cfd.Rvir = 1.0;
	cfd.Rtid = 1.0e6;
	cfd.Z = 0.02;

	cmc_malloc_fits_data_t(&cfd);

	for (i=1; i<=cfd.NOBJ; i++) {
		cfd.obj_id[i] = i;
		cfd.obj_k[i] = 0;
		cfd.obj_m[i] = m[i];
		cfd.obj_Reff[i] = 0.0;
		cfd.obj_r[i] = r[i];
		cfd.obj_vr[i] = vr[i];
		cfd.obj_vt[i] = vt[i];
		cfd.obj_binind[i] = 0;
	}
	
	cmc_write_fits_file(&cfd, filename);

	cmc_free_fits_data_t(&cfd);
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
		double *vt, long int N){
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
	fprintf(stdout,"Before scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
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
	fprintf(stdout,"After scaling:  PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
}
	
int main(int argc, char *argv[]){
	double *X, *r, *vr, *vt, *m;
	double rmax=RMAX, gamma=GAMMA;
	unsigned long int N=NSTAR, seed=SEED;
	char filename[1024];
	int i;
	const char *short_opts = "r:g:N:o:s:dVh";
	const struct option long_opts[] = {
		{"rmax", required_argument, NULL, 'r'},
		{"gamma", required_argument, NULL, 'g'},
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
		case 'g':
			gamma = atof(optarg);
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

	if (gamma == 4.){
		fprintf(stdout,"\nERROR: gamma=4 produces a special case of the hypergeometric function\n");		     fprintf(stdout,"\t\twhich lazy Carl has not implemented yet.\n");
		fprintf(stdout,"\nBut gamma=4 is just a Plummer profile, so just use cmc_mkplummer.\n\n");
	}	

	if (filename[0]==0){
		sprintf(filename, OUTFILE_FORMAT, (double) N);
	}

	dprintf("rmax=%g N=%ld filename=%s seed=%ld\n", rmax, N, filename, seed);

	X = (double *) malloc((N+2)*sizeof(double));
	r = (double *) malloc((N+2)*sizeof(double));
	vr = (double *) malloc((N+2)*sizeof(double));
	vt = (double *) malloc((N+2)*sizeof(double));
	m = (double *) malloc((N+2)*sizeof(double));
	reset_rng_t113(seed);

	check_for_file(filename);
	create_random_array(X, N);
	find_positions(X, r, N, gamma, rmax, 1e-12);
	find_velocities(r, vr, vt, N, gamma, rmax);
	set_masses(m, r, N);
	scale_pos_and_vel(m, r, vr, vt, N);
	write_output_file(m, r, vr, vt, N, filename);

	return 0;
}
