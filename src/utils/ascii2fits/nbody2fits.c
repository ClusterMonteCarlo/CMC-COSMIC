#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <fitsio.h>
#include "../../common/taus113-v2.h"
#include "../../common/fitslib.h"

#define LARGE_DISTANCE 1.0e40
#define PI 3.14159265358979323

/* default parameters */
#define RMAX 300.0
#define NSTAR 10000000UL
#define OUTFILE_FORMAT "%s.fits"
#define SEED 0UL

/* define version */
#define VERSION "0.0.0"
#define NICK "Bad Horsie"
#define DATE "Thu Feb 10 16:58:35 CST 2005"

int debug=0, mtot_given=0, scale_units=1;
int BINFILE=0;

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
	fprintf(stream, "  -i --infile <singlefile> : set name of infile (singles) [<nbody file>.dat]\n");
	fprintf(stream, "  -b --binaries <binfile> : set name of infile (binaries) [<binary file>.dat]\n");
        fprintf(stream, "  -m --total-mass        : the total mass of the cluster [calculated from snapshot]\n");
        fprintf(stream, "  -s --no-scaling        : turn off scaling of positions and velocities\n");
	fprintf(stream, "  -d --debug             : turn on debugging\n");
	fprintf(stream, "  -V --version           : print version info\n");
	fprintf(stream, "  -h --help              : display this help text\n");
	fprintf(stream, " nbody infile should be of the form 'm x y z vx vy vz'\n"); 
	fprintf(stream, " binary infile should be of the form 'e a m1 m2 x y z vx vy vz'\n"); 
	fprintf(stream, " \ta must be in code units (pc/rVir), since scale will not adjust it later\n");
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

void write_output_file(double *m, double *r, double *vr, double *vt, double *m1, double *m2, double *a, double *e, double *binind, long int N, long int NBIN, char *filename){
	long i,j=1;
	cmc_fits_data_t cfd;
	
	cfd.NOBJ = N;
	cfd.NBINARY = NBIN;
	cfd.Mclus = Mtot;
	cfd.Rvir = 1.0;
	cfd.Rtid = 1.0e6;
	cfd.Z = 0.02;
    fprintf(stdout, "%ld %ld\n",N,NBIN);

	cmc_malloc_fits_data_t(&cfd);

	for (i=1; i<=cfd.NOBJ; i++) {
		cfd.obj_id[i] = i;
		cfd.obj_k[i] = 0;
		cfd.obj_m[i] = m[i];
		cfd.obj_Reff[i] = 0.0;
		cfd.obj_r[i] = r[i];
		cfd.obj_vr[i] = vr[i];
		cfd.obj_vt[i] = vt[i];
        //fprintf(stderr, "%ld %lg %ld\n",cfd.obj_id[i],cfd.obj_m[i],cfd.obj_binind[i]);
        if(binind[i] == -1)
            cfd.obj_binind[i] = 0;
        else{
            cfd.obj_id[i] = -1;
            cfd.obj_binind[i] = j;
            cfd.bs_index[j] = i;
            cfd.bs_m1[j] = m1[(long)binind[i]];
            cfd.bs_m2[j] = m2[(long)binind[i]];
            cfd.bs_id1[j] = i;
            cfd.bs_id2[j] = N+j;
            cfd.bs_a[j] = a[(long)binind[i]];
            cfd.bs_e[j] = e[(long)binind[i]];
            //fprintf(stderr, "%lg %lg %lg\n",cfd.obj_m[i],cfd.bs_m1[j],cfd.bs_m2[j]);
            j++;
        }
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

#include "../../common/gensort.h"
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

void sort_radially(double *r, double *vr, double *vt, double *m, double *binind, size_t N) {
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
  if(BINFILE)
      res= gsl_permute(ind->data, binind, 1, N+1);

};


void readNbodySnapshot(const char *filename, const char *binaries, double **r, double **vr, double **vt, double **m, unsigned long *N, unsigned long *NBIN, double **m1, double **m2, double **a, double **e, double **binind) {
  FILE *inputfile;
  FILE *binariesFile;
  int itemsSingle, itemsBinary;
  long i, block, linesSingle, linesBin;
  char *line;
  size_t n;
  double pos[3], vel[3];

  line= NULL;
  block= 1000;
  linesSingle = 0;
  linesBin = 0;
  n= 1;
  itemsSingle = itemsBinary = 1;
  if(!BINFILE)
      itemsBinary = EOF;

  *r= (double *) malloc(sizeof(double));
  *vr= (double *) malloc(sizeof(double));
  *vt= (double *) malloc(sizeof(double));
  *m= (double *) malloc(sizeof(double));
  *m1= (double *) malloc(sizeof(double));
  *m2= (double *) malloc(sizeof(double));
  *a= (double *) malloc(sizeof(double));
  *e= (double *) malloc(sizeof(double));
  *binind= (double *) malloc(sizeof(double));
  **r= DBL_MIN;
  **vr= **vt= 0;
  **m= 0;
  **m1=0;
  **m2=0;
  **a=0;
  **e=0;
  **binind=0;

  inputfile= fopen(filename, "r");
  if(BINFILE)
      binariesFile = fopen(binaries, "r");

  i=1;
  int bin=1;
  while ((itemsSingle!=EOF) || (itemsBinary!=EOF)){
    linesSingle+= block;
    linesBin+= block;
    *r= (double *) realloc(*r, sizeof(double)*linesSingle);
    *vr= (double *) realloc(*vr, sizeof(double)*linesSingle);
    *vt= (double *) realloc(*vt, sizeof(double)*linesSingle);
    *m= (double *) realloc(*m, sizeof(double)*linesSingle);
    *m1= (double *) realloc(*m1, sizeof(double)*linesBin);
    *m2= (double *) realloc(*m2, sizeof(double)*linesBin);
    *a= (double *) realloc(*a, sizeof(double)*linesBin);
    *e= (double *) realloc(*e, sizeof(double)*linesBin);
    *binind= (double *) realloc(*binind, sizeof(double)*linesSingle);
    for (; itemsSingle!=EOF && i<linesSingle; i++) {
      itemsSingle=fscanf(inputfile, 
          "%lf %lf %lf %lf %lf %lf %lf\n",
          &(*m)[i], &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
      if(itemsSingle==EOF){
          i--;
          continue;
        }
      calculate_mc_vars(vel, pos, &(*r)[i], &(*vr)[i], &(*vt)[i]);
      (*binind)[i]=-1;
    }
    for (; itemsBinary!=EOF && i<linesBin; i++) {
      itemsBinary=fscanf(binariesFile, 
          "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
          &(*e)[bin], &(*a)[bin], &(*m1)[bin], &(*m2)[bin], &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
      (*m)[i] = (*m1)[bin]+(*m2)[bin];
      calculate_mc_vars(vel, pos, &(*r)[i], &(*vr)[i], &(*vt)[i]);
      if(itemsBinary==EOF){
          i--;
          continue;
        }
      (*binind)[i]=bin;
        bin++;
    }
  } 

  *N= i-1;
  *NBIN = bin-1;
  (*r)[*N+1]= LARGE_DISTANCE;
  (*vr)[*N+1]= (*vt)[*N+1]= 0.;
  (*m)[*N+1]= 0.;
  (*binind)[*N+1]=0;
  sort_radially(*r, *vr, *vt, *m, *binind, *N);
}

int main(int argc, char *argv[]){
	double *r, *vr, *vt, *m, *m1, *m2, *a, *e, *binind;
	double rmax=RMAX;
	double rcr = 1.0/(6.0*(PI/32.0));
	unsigned long int N=NSTAR, NBIN=NSTAR, seed=SEED;
	char filename[1024]; 
        char input[1024];
        char binaries[1024];
	int i;
	const char *short_opts = "si:b:o:m:dVh";
	const struct option long_opts[] = {
		{"input", required_argument, NULL, 'i'},
		{"binaries", required_argument, NULL, 'b'},
		{"outfile", required_argument, NULL, 'o'},
                {"total-mass", required_argument, NULL, 'm'},
                {"no-scaling", required_argument, NULL, 's'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};

        r= vr= vt= m= m1= m2= a= e= binind= NULL;

	filename[0] = 0;	/* setting a meaningless default name */

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
        case 'i':
			sprintf(input, "%s", optarg);
            break;
        case 'b':
            BINFILE = 1;
			sprintf(binaries, "%s", optarg);
            break;
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
	/*if (optind >= argc) {
		print_usage(stdout);
		return(EXIT_FAILURE);
	}*/

        //input= argv[optind];
	
	if (filename[0]==0){
		sprintf(filename, OUTFILE_FORMAT, input);
	}

	dprintf("rmax=%g N=%ld filename=%s seed=%ld\n", rmax, N, filename, seed);

	check_for_file(filename);
        readNbodySnapshot(input, binaries, &r, &vr, &vt, &m, &N,&NBIN,&m1,&m2,&a,&e,&binind);

        
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
        for(i=1; i<NBIN+1; i++){
          m1[i]/= Mtot;
          m2[i]/= Mtot;
        }

        write_output_file(m, r, vr, vt, m1,m2,a,e, binind, N,NBIN, filename);

	/*for(i=1; i<=N; i++){
		printf("%.30e %.30e %.30e\n", 
				r[i], X[i], calc_M(r[i], rcr, rmax, cms));
	}*/
	return 0;
}
