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
#include "../../common/fitslib.h"

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
	fprintf(stream, "  nbody2fits [options...] <snapshot file>\n");
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

void write_output_file(double *m, double *r, double *vr, double *vt, long int N, char *filename){
	long i;
	cmc_fits_data_t cfd;
	
	cfd.NOBJ = N;
	cfd.NBINARY = 0;
	cfd.Mclus = Mtot;
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

#define GENSORT_NAME                 ind_sort
#define GENSORT_ARGS                 ,r
#define GENSORT_ARGSPROTO            ,double *r
#define GENSORT_TYPE                 size_t
#define GENSORT_KEYTYPE              double 
#define GENSORT_GETKEY(a)            r[a]
#define GENSORT_COMPAREKEYS(k1,k2)   k1<k2

/*
 *#define GENSORT_NISINT
 */

#include "gensort.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_double.h>

void scale_pos_and_vel(double *m, double *r, double *vr, 
		double *vt, long int N, double rcr){
	long int i;
        size_t *ind;
	double PEtot, KEtot, U, T;
	double MM, rfac, vfac;

        ind= (size_t *) calloc(N+2, sizeof(size_t));
        for (i=0; i<N+2; i++) {
          ind[i]= i;
        };

        ind_sort(ind, N+1, r);

	PEtot = KEtot = 0.0;
	U = 0.0;
	MM = 1.0; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
                dprintf("At index %zi %zi\n", ind[i], ind[i+1])
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

void calculate_mc_vars(double *v3d, double *pos, double *r, double *vr, double *vt) {
  /* Calculates the tangential and radial velocity components from the star's 
   * 3d velocity and position*/

  double n[3], vr3d[3];
  int i;

  *r= sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  for (i=0; i<3; i++) {
    n[i]= pos[i]/(*r);
  }

  *vr= n[0]*v3d[0]+n[1]*v3d[1]+n[2]*v3d[2];
  for (i=0; i<3; i++) {
    vr3d[i]= n[i]*(*vr);
  }

  *vt=0.;
  for (i=0; i<3; i++) {
    *vt+= (v3d[i]-vr3d[i])*(v3d[i]-vr3d[i]);
  }

  *vt= sqrt(*vt);
};

void sort_radially(double *r, double *vr, double *vt, double *m, size_t N) {
  int res;
  gsl_permutation *ind;

  ind= gsl_permutation_calloc(N+2);
  gsl_permutation_init(ind);

/*
 *  ind= (long int *) calloc(N+2, sizeof(long));
 *  for (i=0; i<N+2; i++) {
 *    ind[i]= i;
 *  };
 *
 */
  ind_sort(ind->data, N+1, r);
  
  res= gsl_permute(ind->data, r, 1, N+1);
  res= gsl_permute(ind->data, vr, 1, N+1);
  res= gsl_permute(ind->data, vt, 1, N+1);
  res= gsl_permute(ind->data, m, 1, N+1);

};


void readNbodySnapshot(const char *filename, double **r, double **vr, double **vt, double **m, unsigned long *N) {
  FILE *inputfile;
  int items;
  long i, block, lines;
  char *line;
  size_t n;
  double pos[3], vel[3];

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

  i=1;
  do {
    lines+= block;
    *r= (double *) realloc(*r, sizeof(double)*lines);
    *vr= (double *) realloc(*vr, sizeof(double)*lines);
    *vt= (double *) realloc(*vt, sizeof(double)*lines);
    *m= (double *) realloc(*m, sizeof(double)*lines);
    for (; items!=EOF && i<lines-1; i++) {
      items=fscanf(inputfile, 
          "%lf %lf %lf %lf %lf %lf %lf\n",
          &(*m)[i], &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
      calculate_mc_vars(vel, pos, &(*r)[i], &(*vr)[i], &(*vt)[i]);
    };
  } while (items!=EOF);

  *N= i-2;
  (*r)[*N+1]= LARGE_DISTANCE;
  (*vr)[*N+1]= (*vt)[*N+1]= 0.;
  (*m)[*N+1]= 0.;
  sort_radially(*r, *vr, *vt, *m, *N);
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

	check_for_file(filename);
        readNbodySnapshot(input, &r, &vr, &vt, &m, &N);
        
        if (!mtot_given) {
          for (i=1; i<N+1; i++) {
            Mtot+= m[i];
          };
        };
        dprintf("The total mass is %lf\n", Mtot);

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
