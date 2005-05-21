#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <fitsio.h>
#include "taus113-v2.h"

struct imf_param {
	int imf;
	int Ntrace;
	int oflag;
	double mmin;
	double mmax;
	double pl_index;
	double rcr;
	double Cms;
	char infile[1024];
	char outfile[1024];
};

struct star {
	double		r;		/* position			*/
	double		vr;		/* rad. comp. of velocity	*/
	double		vt;		/* tang. comp. of velocity	*/
	double		m;		/* mass				*/
};

struct cluster {
	unsigned long int	NSTAR; 	/* number of stars		*/
};

struct code_par {
	double		time;		/* cluster's age		*/
	int		step;		/* time step we are in		*/
	struct rng_t113_state rng_st;	/* RNG state			*/
};

void printerror( int status) {
	/*****************************************************/
	/* Print out cfitsio error messages and exit program */
	/*****************************************************/
	if (status) {
		fits_report_error(stderr, status); /* print error report */
		exit( status );    /* terminate the program,
				      returning error status */
	}
	return;
}

void write_usage(void){
	printf("Valid options are:\n");
	printf("-i <file> : input file name\n");
	printf("-o <file> : output file name\n");
	printf("-I <int>  : IMF model\n");
	printf("            1: Power-law (default)\n");
	printf("            2: Miller & Scalo (a la starlab)\n");
	printf("            3: Scalo (a la starlab)\n");
	printf("            4: Kroupa (a la starlab)\n");
	printf("            5: Scalo (from Kroupa etal 1993 eq.15)\n");
	printf("            6: Kroupa (old form, for comparison reasons)\n");
	printf("            7: Tracers (1 species, the rest are 1 Msun)\n");
	printf("-w        : overwrite flag\n");
	printf("-m <dbl>  : minimum mass, or tracer mass\n");
	printf("-M <dbl>  : maximum mass\n");
	printf("-p <dbl>  : power-law index (-2.35 is Salpeter)\n");
	printf("-N <int>  : number of tracers\n");
	printf("-r <dbl>  : rcr, the value of r within which average mass");
	printf(" is different\n");
	printf("-C <dbl>  : Cms, how much will the masses be different\n");
	printf("-h        : prints this message and exits\n");
}

void parse_options(struct imf_param *param, int argc, char *argv[]){
	int c;

	opterr = 0;
	
	param->imf = 1;
	param->oflag = 0;
	param->mmin = 0.2;
	param->mmax = 120.0;
	param->pl_index = -2.35;
	param->Ntrace = 0;
	param->rcr = 1.00;
	param->Cms = 1.00;
	(param->infile)[0] = '\0';
	(param->outfile)[0] = '\0';

	while((c = getopt(argc, argv, "i:o:I:wm:M:p:N:r:C:h")) != -1)
		switch(c) {
		case 'i':
			strncpy(param->infile, optarg, 1024);
			break;
		case 'o':
			strncpy(param->outfile, optarg, 1024);
			break;
		case 'I':
			param->imf = strtol(optarg, NULL, 10);
			break;
		case 'N':
			param->Ntrace = strtol(optarg, NULL, 10);
			break;
		case 'w':
			param->oflag = 1;
			break;
		case 'm':
			param->mmin = strtod(optarg, NULL);
			break;
		case 'M':
			param->mmax = strtod(optarg, NULL);
			break;
		case 'p':
			param->pl_index = strtod(optarg, NULL);
			break;
		case 'r':
			param->rcr = strtod(optarg, NULL);
			break;
		case 'C':
			param->Cms = strtod(optarg, NULL);
			break;
		case 'h':
			write_usage();
			exit(EXIT_SUCCESS);
			break;
		case '?':
			if (isprint (optopt)) 
				printf( "Unknown option `-%c'.\n", optopt);
			else 
				printf("Unknown option character `\\x%x'.\n",
					       optopt);
			write_usage();
			exit(EXIT_FAILURE);
      		default:
			printf("This can't happen!\n");
        		exit(EXIT_FAILURE);
		}
	if ((param->infile)[0] == '\0'){
		printf("input file not given or invalid\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if ((param->outfile)[0] == '\0'){
		printf("output file not given or invalid\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if (param->mmin<=0){
		printf("Invalid mmin value\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if (param->mmax<=param->mmin){
		printf("Invalid mmax value, less than mmin\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if (param->imf<1 || param->imf>7){
		printf("Invalid IMF model value\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if (param->imf==7 && param->Ntrace<1 ){
		printf("Number of tracers has to be positive\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	if (param->Cms != 1.00 && param->imf != 1){
		printf("mass segregation is implemented only for");
		printf(" power-law IMF\n");
		exit(EXIT_FAILURE);
	}
	if (param->Cms > 1.00){
		printf("parameter Cms has to be less than 1.0\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	
	fprintf(stderr, "infile=%s outfile=%s imf=%d Mmin=%g Mmax=%g pl_index=%g rcr=%g Cms=%g\n", 
		param->infile, param->outfile, param->imf, param->mmin, param->mmax, param->pl_index, 
		param->rcr, param->Cms);
}

void read_input(struct imf_param param, struct star *s[],
		struct cluster *c, struct code_par *cp){
	fitsfile *fptr;
	int status, hdunum, hdutype, anynull;
	long frow, felem, nelem;
	float floatnull;
	unsigned long int i, NSTAR;
	double *mass, *r, *vr, *vt;

	status = 0;
	
	/* open file */
	fits_open_file(&fptr, param.infile, READONLY, &status);
	/* move to proper part of file */
	hdunum = 2; 		/* data table is in second HDU */
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);

	/* read headers
	 * set cluster parameters and 
	 * allocate memory */
	fits_read_key(fptr, TULONG, "NSTAR", &NSTAR, NULL, &status);
	c->NSTAR = NSTAR;
	fits_read_key(fptr, TDOUBLE, "Time", &(cp->time), NULL, &status);
	fits_read_key(fptr, TINT, "Step", &(cp->step), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z1", &(cp->rng_st.z1), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z2", &(cp->rng_st.z2), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z3", &(cp->rng_st.z3), NULL, &status);
	fits_read_key(fptr, TULONG, "RNG_Z4", &(cp->rng_st.z4), NULL, &status);
	printerror(status);
	*s = malloc((NSTAR+2)*sizeof(struct star));
	mass = malloc((NSTAR+2)*sizeof(double));
	vr = malloc((NSTAR+2)*sizeof(double));
	vt = malloc((NSTAR+2)*sizeof(double));
	r = malloc((NSTAR+2)*sizeof(double));

	/* read the columns into variables and copy variables to
	 * struct star s[]                                       */
	frow = 1; felem = 1; floatnull = 0.;  
	nelem = NSTAR+2;
	fits_read_col(fptr, TDOUBLE, 1, frow, felem, nelem, &floatnull, mass, 
			&anynull, &status);
	fits_read_col(fptr, TDOUBLE, 2, frow, felem, nelem, &floatnull, r, 
			&anynull, &status);
	fits_read_col(fptr, TDOUBLE, 3, frow, felem, nelem, &floatnull, vr, 
			&anynull, &status);
	fits_read_col(fptr, TDOUBLE, 4, frow, felem, nelem, &floatnull, vt, 
			&anynull, &status);
	for(i=0; i<=NSTAR+1; i++){
		(*s)[i].m = mass[i];
		(*s)[i].r = r[i];
		(*s)[i].vr = vr[i];
		(*s)[i].vt = vt[i];
	}
			
	free(mass); free(vr); free(vt); free(r);
	fits_close_file(fptr, &status);
	printerror(status);
}

double calc_f_mminp(double mminp,
	       	double mmin, double mmax, double alpha, double cms){
	/* XXX alpha<0 XXX */
	/* this buddy here returns <m>_mminp - (2/cms-1)*<m>_mmin */
	/* supposedly the function f=<m>_mminp - (2/cms-1)*<m>_mmin
	 * is <0 for mminp=mmin and >0 mminp=mmax */
	double avem, avemp;
	double f;

	avem = (alpha+2)*(pow(mmax, alpha+2)-pow(mmin, alpha+2))
		/ ((alpha+1)*(pow(mmax, alpha+1)-pow(mmin, alpha+1)));
	if (mminp==mmax){
		f = mmax - (2.0/cms-1.0)*avem;
		return f;
	}
	avemp = (alpha+2)*(pow(mmax, alpha+2)-pow(mminp, alpha+2))
		/ ((alpha+1)*(pow(mmax, alpha+1)-pow(mminp, alpha+1)));
	f = avemp - (2.0/cms-1.0)*avem;
	return f;
}

double find_mminp(double cms, double m_abs_min, double m_abs_max,
	       double alpha){
	/* XXX alpha<0 XXX */
	/* this function returns a mminp value such that
	 * <m>_mminp = (2/cms-1)*<m>mmin
	 * that 2 is due to selection in the actual code 
	 * uses binary search (bisection) and external function */
	/* supposedly the function f=<m>_mminp - (2/cms-1)*<m>_mmin
	 * is <0 for mminp=mmin and >0 mminp=mmax 
	 * this brings restrictions on cms of course!!! */
	double mmin, mmax, mtry, f;
	double tol = 1e-8;
	int k;

	k = 0;
	mmin = m_abs_min;
	mmax = m_abs_max;
	while(mmax-mmin > tol){
		mtry = (mmax+mmin)/2.0;
		f = calc_f_mminp(mtry,  m_abs_min, m_abs_max, alpha, cms);
		if(f>0){
			mmax=mtry;
		} else {
			mmin=mtry;
		}
		if(k++>80){
			fprintf(stderr,"too many iterations...\n");
			fprintf(stderr,"in find_mminp bisection.\n");
			break;
		}
	}
	return (mmax+mmin)/2.0;
}

void set_masses(struct imf_param param, struct star *s[],
			struct cluster c, struct code_par cp){
	unsigned long int i;
	double m, X, X2, norm, tmp;
	double total_mass, mmin, mmin_ms;
	double Xcrit, C1, C2, C3, C4;

	set_rng_t113(cp.rng_st);
	/* to change the RND seed we make call(s) to RNG */
	X = rng_t113_dbl();
	X = rng_t113_dbl();
	X = rng_t113_dbl();
	X = rng_t113_dbl();
	total_mass = 0.0;
	if(param.imf==1 && param.Cms==1.0){  /* Power-law w/o ms     */
		for(i=1; i<=c.NSTAR; i++){
			tmp = param.pl_index+1.0;
			X = rng_t113_dbl();
			if(param.pl_index==-1.0){
				norm = log(param.mmax/param.mmin);
				m = param.mmin*exp(norm*X);
			} else {
				norm = pow(param.mmax/param.mmin, tmp) - 1.0;
				m = param.mmin*pow(norm*X+1, 1.0/tmp);
			}
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==1){ /* Power-law w/ mass segregation    */
		/* XXX maybe there should be a warning for Cms~=1.0 XXX */
		mmin_ms = find_mminp(param.Cms, param.mmin, param.mmax,
						param.pl_index);
		for(i=1; i<=c.NSTAR; i++){
			tmp = param.pl_index+1.0;
			X = rng_t113_dbl();
			X2 = rng_t113_dbl();
			if (X2<0.5 && (*s)[i].r<param.rcr){
				mmin = mmin_ms;
			} else {
				mmin = param.mmin;
			}
			if(param.pl_index==-1.0){
				norm = log(param.mmax/mmin);
				m = mmin*exp(norm*X);
			} else {
				norm = pow(param.mmax/mmin, tmp) - 1.0;
				m = mmin*pow(norm*X+1, 1.0/tmp);
			}
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==2){ /* Miller & Scalo */
		for(i=1; i<=c.NSTAR; i++){
			do {
				X = rng_t113_dbl();
				m = 0.19*X
	    			/ (pow(1-X, 0.75) + 0.032*pow(1-X, 0.25));
			} while (m<param.mmin || m>param.mmax) ;
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==3){ /* Scalo          */
		for(i=1; i<=c.NSTAR; i++){
			do {
				X = rng_t113_dbl();
				m = 0.3*X / pow(1-X, 0.55);
			} while (m<param.mmin || m>param.mmax) ;
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==4){ /* Kroupa         */
		for(i=1; i<=c.NSTAR; i++){
			do {
				X = rng_t113_dbl();
				m = 0.08 + (0.19*pow(X,1.55) + 0.05*pow(X,0.6))
	    			/  pow(1-X,0.58);
			} while (m<param.mmin || m>param.mmax) ;
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==5){ /* Scalo (from Kroupa etal 1993 eq.15) */
		for(i=1; i<=c.NSTAR; i++){
			do {
				X = rng_t113_dbl();
				m = 0.284*pow(X,0.337)
				/(pow(1-X,0.5)-0.015*pow(1.0-X,0.085));
			} while (m<param.mmin || m>param.mmax);
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==6){ /* Kroupa (old)         */
		C1 = 1.903988317884160571490719612429355684076;
		C2 = 0.9542546378990059541285953683909362471882;
		C3 = 1.000276574531978107760263208804720707748;
		C4 = 0.1101063043729622254763763886604926439063;
		Xcrit = 0.7291630515263233686411262877173302148501;
		for(i=1; i<=c.NSTAR; i++){
			do {
				X = rng_t113_dbl();
				if (X<Xcrit){
					m = pow(C2/(C1-X), 0.3);
				} else {
					m = pow(C4/(C3-X), 1.3);
				}
			} while (m<param.mmin || m>param.mmax);
			(*s)[i].m = m;
			total_mass += m;
		}
	} else if (param.imf==7){ /* Tracers */
		if(c.NSTAR <= param.Ntrace){
			printf("Number of tracers = %d\n", param.Ntrace);
			printf("Number of stars = %ld\n", c.NSTAR);
			printf("something wrong\n");
			exit(EXIT_FAILURE);
		}
		for(i=1; i<=c.NSTAR; i++){
			m = 1.0;
			(*s)[i].m = m;
			total_mass += m;
		}
		/* note that the following algorithm puts trace stars 
		 * uniformly. for a random selection algorithm, see:
		 *  Knuth, vol2, section 3.4.3, Algorithm S */
		for(i=0; i<param.Ntrace; i++){
			m = param.mmin;
			(*s)[c.NSTAR / (2*param.Ntrace) 
				+ i*c.NSTAR/param.Ntrace + 1].m = m;
			total_mass += m - 1.0;
		}
	} else {
		printf("This can't happen!\n");
		exit(EXIT_FAILURE);
	}

	for(i=0; i<=c.NSTAR; i++){
		(*s)[i].m /= total_mass;
	}
}

void write_output(struct star *s, struct cluster c, 
		struct code_par cp, struct imf_param param){
	unsigned long int i;
	
	fitsfile *fptr;
	int status;
	long firstrow, firstelem;
	int tfields;       /* table will have n columns */
	long nrows;	

	char extname[] = "CLUSTER_STARS";          /* extension name */
	char *filename;
	char *ttype[] = { "Mass",  "Position", "vr",    "vt" };
	char *tform[] = { "1D",    "1D",       "1D",    "1D" };
	char *tunit[] = { "Nbody", "Nbody",    "Nbody", "Nbody" };
	
	double *mass, *r, *vr, *vt;

	/* these go to header */
	int tstep = 0;
	double time = 0.0;
	
	get_rng_t113(&(cp.rng_st));
	mass = malloc((c.NSTAR+2)*sizeof(double));
	r = malloc((c.NSTAR+2)*sizeof(double));
	vr = malloc((c.NSTAR+2)*sizeof(double));
	vt = malloc((c.NSTAR+2)*sizeof(double));
	for(i=0; i<=c.NSTAR+1; i++){
		mass[i] = s[i].m;
		r[i] = s[i].r;
		vr[i] = s[i].vr;
		vt[i] = s[i].vt;
	}
	
	status = 0;
	tfields = 4;
	filename = param.outfile;
	nrows = c.NSTAR+2;
	fits_create_file(&fptr, filename, &status);
	fits_open_file(&fptr, filename, READWRITE, &status);

	fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status);
	fits_update_key(fptr, TLONG, "NSTAR", &(c.NSTAR), 
			"No of Stars", &status);
	fits_update_key(fptr, TDOUBLE, "Time", &time, 
			"Age of cluster", &status);
	fits_update_key(fptr, TLONG, "Step", &tstep, 
			"Iteration Step", &status);
	fits_update_key(fptr, TULONG, "RNG_Z1", &(cp.rng_st.z1), 
			"RNG STATE Z1", &status);
	fits_update_key(fptr, TULONG, "RNG_Z2", &(cp.rng_st.z2), 
			"RNG STATE Z2", &status);
	fits_update_key(fptr, TULONG, "RNG_Z3", &(cp.rng_st.z3), 
			"RNG STATE Z3", &status);
	fits_update_key(fptr, TULONG, "RNG_Z4", &(cp.rng_st.z4), 
			"RNG STATE Z4", &status);

	firstrow  = 1;  /* first row in table to write   */
	firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

	fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, mass,
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


int main(int argc, char *argv[]){
	struct imf_param param;
	struct star *s;
	struct cluster clus;
	struct code_par cp;
	/* int i; */
	
	parse_options(&param, argc, argv);
	read_input(param, &s, &clus, &cp);
	set_masses(param, &s, clus, cp);
	/* I might add scaling to E0 = -1/4 here */
	write_output(s, clus, cp, param);

	/*for(i=1; i<=clus.NSTAR; i++){
		printf("%e\n", s[i].m);
	}*/

	return 0;
}

