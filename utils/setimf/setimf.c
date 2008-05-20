#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "../../common/fitslib.h"
#include "../../common/taus113-v2.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_double.h>


#define SEED 768364873UL

/* a fast square function */
double sqr(double x)
{
        return(x*x);
}

struct imf_param {
	int imf;
	int Ntrace;
	int oflag;
        int n_neutron;
	double mmin;
	double mmax;
	double pl_index;
	double rcr;
	double Cms;
        double bhmass;
        int scale;
        long N_mass_seg;
	char infile[1024];
	char outfile[1024];
};

void write_usage(void){
	printf("Valid options are:\n");
	printf("-i <file> : input file name\n");
	printf("-o <file> : output file name\n");
	printf("-I <int>  : IMF model\n");
	printf("            0: Kroupa, 2001. MNRAS 322, 231-246 (Mmin=0.01)\n");
	printf("            1: Power-law (default)\n");
	printf("            2: Miller & Scalo (a la starlab)\n");
	printf("            3: Scalo (a la starlab)\n");
	printf("            4: Kroupa (a la starlab; Kroupa, Tout & Gilmore 1993, MNRAS 262, 545; Mmin=0.08)\n");
	printf("            5: Scalo (from Kroupa etal 1993 eq.15)\n");
	printf("            6: Kroupa (old form, for comparison reasons)\n");
	printf("            7: Tracers (1 species, the rest are 1 Msun)\n");
	printf("            8: Power-law with neutron stars (objects at 1.4 Msun)\n");
        printf("            9: single mass (set by -m)\n");
	printf("-w        : overwrite flag\n");
	printf("-m <dbl>  : minimum mass, or tracer mass\n");
	printf("-M <dbl>  : maximum mass\n");
	printf("-p <dbl>  : power-law index (-2.35 is Salpeter)\n");
	printf("-u <int>  : number of neutron stars\n");
	printf("-N <int>  : number of tracers\n");
        printf("-b <dbl>  : add a central black hole with mass <dbl>\n");
        printf("-s        : scale positions and velocities such that the cluster is in exact ");
        printf("virial equilibrium (T/2W = 1)\n");
        printf("-S <int>  : Create a maximally mass segregated cluster according to Baumgardt et al. (2008).\n");
        printf("          : <int> is the number of stars which must be smaller or equal to NSTAR * m_low/<m>,\n");
        printf("          : where NSTAR is the number of stars in the input file, m_low the lowest mass star,\n"); 
        printf("          : and <m> the mean stellar mass. For <int> = -1, the largest value is assumed.\n");
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
	param->n_neutron = 0;
	param->Ntrace = 0;
	param->rcr = 1.00;
	param->Cms = 1.00;
        param->bhmass = 0.0;
        param->scale = 0;
        param->N_mass_seg= 0;
	(param->infile)[0] = '\0';
	(param->outfile)[0] = '\0';

	while((c = getopt(argc, argv, "i:o:I:wm:M:p:u:N:r:C:hsb:S:")) != -1)
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
		case 'u':
		        param->n_neutron = strtol(optarg, NULL, 10);
			break;
		case 'r':
			param->rcr = strtod(optarg, NULL);
			break;
		case 'C':
			param->Cms = strtod(optarg, NULL);
			break;
		case 'b':
			param->bhmass = strtod(optarg, NULL);
			break;
                case 's':
                        param->scale = 1;
                        break;
                case 'S':
                        param->N_mass_seg = strtol(optarg, NULL, 10);
                        if (param->N_mass_seg == 0) {
                          printf("WARNING: Setting -S to 0 disables the generation of ");
                          printf("a mass segregated cluster\n");
                          printf("       : Set it to -1 if you want a mass segregated cluster with the\n");
                          printf("       : maximum number of stars (see setimf --help)\n");
                        };
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
	if (param->imf<0 || param->imf>9){
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
	if (param->imf==8 && param->n_neutron < 0){
	        printf("Number of neutron stars must be positive");
		write_usage();
		exit(EXIT_FAILURE);
	}
        if (param->bhmass<0){
		printf("ANTIGRAVITY! Wooaa... cool! (black hole mass is negative *grin*)\n");
		write_usage();
		exit(EXIT_FAILURE);
	}
	
	fprintf(stderr, "infile=%s outfile=%s imf=%d Mmin=%g Mmax=%g pl_index=%g rcr=%g Cms=%g\n", 
		param->infile, param->outfile, param->imf, param->mmin, param->mmax, param->pl_index, 
		param->rcr, param->Cms);
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

double set_masses(struct imf_param param, cmc_fits_data_t *cfd){
	long i, j;
	double Mass[4], alpha[4], Cons[4], Xlim[5];
	double m, X, X2, norm, tmp, n_rat, ncheck;
	double total_mass, mmin, mmin_ms;
	double Xcrit, C1, C2, C3, C4;
        
	reset_rng_t113(SEED);
	total_mass = 0.0;
	n_rat = ((double) param.n_neutron / (double) cfd->NOBJ);
	if(param.imf==1 && param.Cms==1.0){  /* Power-law w/o ms     */
		for(i=1; i<=cfd->NOBJ; i++){
		        
			tmp = param.pl_index+1.0;
			X = rng_t113_dbl();
			if(param.pl_index==-1.0){
				norm = log(param.mmax/param.mmin);
				m = param.mmin*exp(norm*X);
			} else {
				norm = pow(param.mmax/param.mmin, tmp) - 1.0;
				m = param.mmin*pow(norm*X+1, 1.0/tmp);
			}
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==0){ /* Kroupa, 2001. MNRAS 322, 231-246 */
	        Cons[0] = 1.986846095810826;
		Cons[1] = 0.15894768766486606;
		Cons[2] = 0.07947384383243306;
		Cons[3] = 0.07947384383243306;
		alpha[0] = 0.3;
		alpha[1] = 1.3;
		alpha[2] = 2.3;
		alpha[3] = 2.3;
		Xlim[0] = 0.0;
		Xlim[1] = 0.371431122772297;
		Xlim[2] = 0.8494711094748518;
		Xlim[3] = 0.9388662739750515;
		Xlim[4] = 1.0;
		Mass[0] = 0.01;
		Mass[1] = 0.08;
		Mass[2] = 0.5;
		Mass[3] = 1.0;
		
		/* implement broken power-law, and use rejection for new limits */
		for(i=1; i<=cfd->NOBJ; i++){
		        do {
				X = rng_t113_dbl();
				j = 0;
				while (X > Xlim[j+1]) {
				  j++;
				}
				m = pow((1.0-alpha[j])/Cons[j]*(X-Xlim[j])+pow(Mass[j],1.0-alpha[j]), 1.0/(1.0-alpha[j]));
				if (isnan(m)) {
				  fprintf(stderr, "Oops!  m=NaN.  Please make coefficients more precise.\n");
				  exit(-127);
				}
			} while (m<param.mmin || m>param.mmax) ;
			/* fprintf(stderr, "X=%g Xlim[j]=%g j=%d m=%g\n", X, Xlim[j], j, m); */
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==1){ /* Power-law w/ mass segregation    */
		/* XXX maybe there should be a warning for Cms~=1.0 XXX */
		mmin_ms = find_mminp(param.Cms, param.mmin, param.mmax,
						param.pl_index);
		for(i=1; i<=cfd->NOBJ; i++){
			tmp = param.pl_index+1.0;
			X = rng_t113_dbl();
			X2 = rng_t113_dbl();
			if (X2<0.5 && cfd->obj_r[i]<param.rcr){
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
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==2){ /* Miller & Scalo */
		for(i=1; i<=cfd->NOBJ; i++){
			do {
				X = rng_t113_dbl();
				m = 0.19*X
	    			/ (pow(1-X, 0.75) + 0.032*pow(1-X, 0.25));
			} while (m<param.mmin || m>param.mmax) ;
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==3){ /* Scalo          */
		for(i=1; i<=cfd->NOBJ; i++){
			do {
				X = rng_t113_dbl();
				m = 0.3*X / pow(1-X, 0.55);
			} while (m<param.mmin || m>param.mmax) ;
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==4){ /* Kroupa         */
		for(i=1; i<=cfd->NOBJ; i++){
			do {
				X = rng_t113_dbl();
				m = 0.08 + (0.19*pow(X,1.55) + 0.05*pow(X,0.6))
	    			/  pow(1-X,0.58);
			} while (m<param.mmin || m>param.mmax) ;
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==5){ /* Scalo (from Kroupa etal 1993 eq.15) */
		for(i=1; i<=cfd->NOBJ; i++){
			do {
				X = rng_t113_dbl();
				m = 0.284*pow(X,0.337)
				/(pow(1-X,0.5)-0.015*pow(1.0-X,0.085));
			} while (m<param.mmin || m>param.mmax);
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==6){ /* Kroupa (old)         */
		C1 = 1.903988317884160571490719612429355684076;
		C2 = 0.9542546378990059541285953683909362471882;
		C3 = 1.000276574531978107760263208804720707748;
		C4 = 0.1101063043729622254763763886604926439063;
		Xcrit = 0.7291630515263233686411262877173302148501;
		for(i=1; i<=cfd->NOBJ; i++){
			do {
				X = rng_t113_dbl();
				if (X<Xcrit){
					m = pow(C2/(C1-X), 0.3);
				} else {
					m = pow(C4/(C3-X), 1.3);
				}
			} while (m<param.mmin || m>param.mmax);
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	} else if (param.imf==7){ /* Tracers */
		if(cfd->NOBJ <= param.Ntrace){
			printf("Number of tracers = %d\n", param.Ntrace);
			printf("Number of stars = %ld\n", cfd->NOBJ);
			printf("something wrong\n");
			exit(EXIT_FAILURE);
		}
		for(i=1; i<=cfd->NOBJ; i++){
			m = 1.0;
			cfd->obj_m[i] = m;
			total_mass += m;
		}
		/* note that the following algorithm puts trace stars 
		 * uniformly. for a random selection algorithm, see:
		 *  Knuth, vol2, section 3.4.3, Algorithm S */
		for(i=0; i<param.Ntrace; i++){
			m = param.mmin;
			cfd->obj_m[cfd->NOBJ / (2*param.Ntrace) 
				+ i*cfd->NOBJ/param.Ntrace + 1] = m;
			total_mass += m - 1.0;
		}
	  } else if(param.imf==8 && param.Cms==1.0){  /* Power-law w neutron stars */
	    for(i=1; i<=cfd->NOBJ; i++){
	    ncheck = rng_t113_dbl();
	    tmp = param.pl_index+1.0;
	    X = rng_t113_dbl();
	    /* Is this a regular star? */
	    if (ncheck > n_rat){
	    if(param.pl_index==-1.0){
	      norm = log(param.mmax/param.mmin);
	      m = param.mmin*exp(norm*X);
	    } else {
	      norm = pow(param.mmax/param.mmin, tmp) - 1.0;
	      m = param.mmin*pow(norm*X+1, 1.0/tmp);
	    }
	    /* If neutron star, assign to 1.4 solar masses */
	      } else
		{
		  m = 1.4;
		} 
	    cfd->obj_m[i] = m;
	    total_mass += m;
	    }
	  }
          else if (param.imf==9){ 
		for(i=1; i<=cfd->NOBJ; i++){
			m = param.mmin;
			cfd->obj_m[i] = m;
			total_mass += m;
		}
	  }
	else {
		printf("This can't happen!\n");
		exit(EXIT_FAILURE);
	}

	cfd->Mclus = total_mass;
	
	// code commented out by Stefan
        //total_mass+= param.bhmass;
	// rescale masses to N-body units
	for(i=0; i<=cfd->NOBJ; i++){
		cfd->obj_m[i] /= total_mass;
	}

	return(total_mass);
}

void scale_pos_and_vel(struct imf_param param, cmc_fits_data_t *cfd, double total_mass){
	long int i, N;
	double PEtot, KEtot, U, T;
	double MM, rfac, vfac;
	
	PEtot = KEtot = 0.0;
        N= cfd->NOBJ;
	U = 0.0;
	MM = 1.0+ param.bhmass/total_mass; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
		U -= MM*(1.0/cfd->obj_r[i] - 1.0/cfd->obj_r[i+1]);
		T = 0.5 * (cfd->obj_vr[i] * cfd->obj_vr[i] + cfd->obj_vt[i] * cfd->obj_vt[i]);
		MM -= cfd->obj_m[i];
		PEtot += 0.5*U* cfd->obj_m[i];
		KEtot += T* cfd->obj_m[i];
	}
        if (param.bhmass>0.) {
          U-= -MM/cfd->obj_r[1];
          PEtot+= 0.5*U * cfd->obj_m[0];
        };

	printf("Before scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
	/* scaling position and velocity */
	rfac = -PEtot*2.0;
	vfac = 1.0/sqrt(4.0*KEtot);
        printf("rfac= %lf, vfac= %lf\n", rfac, vfac);
	for(i=1; i<=N; i++){
		cfd->obj_r[i] *= rfac;
		cfd->obj_vr[i] *= vfac;
		cfd->obj_vt[i] *= vfac;
	}

	PEtot = KEtot = 0.0;
	U = 0.0;
	MM = 1.0+ param.bhmass/total_mass; /* because of units, the total mass has to be 1 initially */
	for(i=N; i>=1; i--){
		U -= MM*(1.0/cfd->obj_r[i] - 1.0/cfd->obj_r[i+1]);
		T = 0.5 * (cfd->obj_vr[i] * cfd->obj_vr[i] + cfd->obj_vt[i] * cfd->obj_vt[i]);
		MM -= cfd->obj_m[i];
		PEtot += 0.5*U* cfd->obj_m[i];
		KEtot += T* cfd->obj_m[i];
	}
        if (param.bhmass>0.) {
          U-= -MM/cfd->obj_r[1];
          PEtot+= 0.5*U*cfd->obj_m[0];
        };

	printf("After  scaling: PEtot = %f, KEtot = %f, vir rat = %f\n", 
			PEtot, KEtot, KEtot/PEtot);
}

/** 
 * @brief Calculates the maximum number of stars for a mass segregated cluster 
 * for the recipe  of Baumgardt et al. (2008)
 *
 * The maximum number of stars the generated mass segregated cluster can have depends on
 * the number N' of stars of the 'configuration' cluster (the cluster that determines the 
 * density profile), as well as the minimum and average mass of the stars, m_min and 
 * m_ave respectively. N_max= N' * m_min/m_ave .
 * 
 * @param s Array containing the cluster members of the 'configuration' cluster.
 * @param clus Cluster mass and number of stars.
 * 
 * @return Maximium number of stars for mass segregated cluster.
 */
long calc_n_mass_seg_max(cmc_fits_data_t *cfd) {
  long i;
  double m_ave, m_min;

  m_ave= cfd->Mclus/cfd->NOBJ;
  printf("Mclus= %g, NOBJ=%li\n", cfd->Mclus, cfd->NOBJ);

  m_min= cfd->obj_m[1];

  for (i=2; i< cfd->NOBJ+1; i++) {
    if (cfd->obj_m[i] < m_min) m_min= cfd->obj_m[i];
  }

  m_min*= cfd->Mclus;
  printf("m_ave= %g, m_min=%g\n", m_ave, m_min);
  /* The "-2" is for safety - keep fingers crossed that this is enough */
  return ((long)(cfd->NOBJ*m_min/m_ave)-2);
}

/** 
 * @brief Calculates the total mechanical energy of each star.
 * 
 * This routine assumes N-body units for all quantities, i.e.
 * G=1, M_total=1, E_total=-1/4.
 *
 * @param s Pointer to star array.
 * @param clus Cluster paramters. Only NSTARS is used.
 * 
 * @return Returns an array of total energies for each star in (*s)[1:NSTARS+1]
 * (python slice notation!)
 */
double * calc_energy(cmc_fits_data_t *cfd) {
  double *E, phi;
  long double m_total, mprev;
  long i;

  E= (double *) calloc(cfd->NOBJ, sizeof(double));

  /* This is just for safety, the masses should already be normalized to 1 */
  m_total= 0.;
  for (i=1; i<cfd->NOBJ+1; i++) {
    m_total+= cfd->obj_m[i];
  }

  /* calculate energy for outermost star */
  mprev= 1.; phi= 0.; i= cfd->NOBJ;

  phi= - mprev/cfd->obj_r[i];
  E[i-1]= phi+ 0.5* (sqr(cfd->obj_vr[i]) + sqr(cfd->obj_vt[i]));
  mprev-= cfd->obj_m[i]/m_total;

  for (i=cfd->NOBJ-1; i>0; i--) {
    phi= phi - mprev * (1. / cfd->obj_r[i] - 1. / cfd->obj_r[i+1] );
    E[i-1]= phi+ 0.5  * (sqr(cfd->obj_vr[i]) + sqr(cfd->obj_vt[i]));
    mprev-= cfd->obj_m[i]/m_total;
  }

  return (E);
}

/* routine to sort the mass array in decreasing order */
#define GENSORT_NAME                 mass_sort
#define GENSORT_TYPE                 double
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            a
#define GENSORT_COMPAREKEYS(k1,k2)   k1 > k2

#include "../../common/gensort.h"

/* routine to sort an array E in increasing order indirectly*/
#define GENSORT_NAME                 table_sort
#define GENSORT_ARGS                 ,E
#define GENSORT_ARGSPROTO            ,double *E
#define GENSORT_TYPE                 size_t
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            E[a]
#define GENSORT_COMPAREKEYS(k1,k2)   k1<k2

#include "../../common/gensort.h"

/** 
 * @brief Create a fully mass segregated cluster from the 'configuration' cluster s.
 *
 * This is based on the recipe given in Baumgardt et al. (2008). The 'configuration' 
 * cluster determines the density distribution and is generated the standard way, e.g.,
 * mkplummer, mkking. 
 *  
 * @param param Provides N_mass_seg, the number of stars for the mass segregated cluster.
 * @param cfd Contains the 'configuration' cluster.
 * 
 * @return The fully mass segregated cluster.
 */
cmc_fits_data_t * create_mass_seg(struct imf_param param, cmc_fits_data_t *cfd)
{
  size_t N_config= cfd->NOBJ, N_mass= param.N_mass_seg, *ind, i;
  size_t lower_i, higher_i, index;
  double *E, *Mcum, X;
  cmc_fits_data_t *ms;
  gsl_permutation *perm;
  int res;

  /* Allocate fits structure */
  ms= (cmc_fits_data_t *) calloc(1, sizeof(cmc_fits_data_t));

  /* Initialize the new mass segregated cluster */
  ms->NOBJ= N_mass;
  ms->NBINARY= cfd->NBINARY;
  cmc_malloc_fits_data_t(ms);

  E= calc_energy(cfd);

  /* Initialize the indirection table ... */
  ind= (size_t *) calloc(N_config, sizeof(size_t));
  for (i=0; i< N_config; i++) ind[i]= i;

  /* ... and get the ordered indices for the energy */
  table_sort(ind, N_config, E);

  /* Set the masses for the new particle number N_mass */
  /* Note: set_masses already normalizes them wrt. the total cluster
   * mass.
   */
  set_masses(param, ms);

  /* Read the masses into the masses array and sort them in descending 
   * order. (This is a bit ugly but spares me to rewrite 
   * substantial parts of the code)
   */
  Mcum = (double *) calloc(N_mass, sizeof(double)); 

  mass_sort(&(ms->obj_m[1]), N_mass);

  Mcum[0]= ms->obj_m[1];
  for (i=1; i< N_mass; i++) {
    Mcum[i]= Mcum[i-1]+ ms->obj_m[i+1];
  }

  /* For each mass in the mass array choose a random position and velocity 
   * between N_config*Mcum[i-1] and N_config*Mcum[i]. 
   */
  for (i=0; i< N_mass; i++) {
    long j;

    if (i==0) {
      lower_i= 0;
    } else {
      lower_i= N_config* Mcum[i-1];
    }
    higher_i= N_config*Mcum[i];

    /* choose random index between [lower_i, higher_i) */
    /* FIXME: I have chosen the inclusive/exclusive boundaries according 
     * to what I thought might make sense but can't tell if this choice 
     * was intended by Baumgardt et al. (2008).
     */
    X = rng_t113_dbl();
    index= (long) floor(lower_i + (higher_i-lower_i)*X);

    /* Transfer positions and velocities from s[index] to ms[i] */
    ms->obj_r [i+1] = cfd->obj_r [ind[index]+1];
    ms->obj_vr[i+1] = cfd->obj_vr[ind[index]+1];
    ms->obj_vt[i+1] = cfd->obj_vt[ind[index]+1];

    /* copy the rest of the structure */
    ms->obj_id[i+1] = i+1;
    ms->obj_k [i+1] = cfd->obj_k [ind[index]+1];
    ms->obj_Reff[i+1] = cfd->obj_Reff[ind[index]+1];
    ms->obj_binind[i+1] = cfd->obj_binind[ind[index]+1];
    for (j=0; j< cfd->NBINARY+1; j++) {
      ms->bs_index[j] = cfd->bs_index[j];
      ms->bs_id1[j] = cfd->bs_id1[j]; 
      ms->bs_k1[j] = cfd->bs_k1[j];
      ms->bs_m1[j] = cfd->bs_m1[j];
      ms->bs_Reff1[j] = cfd->bs_Reff1[j];
      ms->bs_id2[j] = cfd->bs_id2[j];
      ms->bs_k2[j] = cfd->bs_k2[j];
      ms->bs_m2[j] = cfd->bs_m2[j];
      ms->bs_Reff2[j] = cfd->bs_Reff2[j];
      ms->bs_a[j] = cfd->bs_a[j];
      ms->bs_e[j] = cfd->bs_e[j];
    }
    ms->Rvir = cfd->Rvir;
    ms->Rtid = cfd->Rtid;
    ms->Z    = cfd->Z;
  }

  /* need to sort radially for cmc*/
  perm= gsl_permutation_calloc(N_mass);
  gsl_permutation_init(perm);
  
  table_sort(perm->data, N_mass, &(ms->obj_r[1]));
  
  res= gsl_permute(perm->data, &(ms->obj_r[1]), 1, N_mass);
  res= gsl_permute(perm->data, &(ms->obj_vr[1]), 1, N_mass);
  res= gsl_permute(perm->data, &(ms->obj_vt[1]), 1, N_mass);
  res= gsl_permute(perm->data, &(ms->obj_m[1]), 1, N_mass);

  free(Mcum); free(ind); free(E);

  return (ms);
}

int main(int argc, char *argv[]){
	cmc_fits_data_t *cfd, *ms;
	struct imf_param param;
	double total_mass;
	/* int i; */
	
	parse_options(&param, argc, argv);

        /* Allocate fits structure */
        cfd= (cmc_fits_data_t *) calloc(1, sizeof(cmc_fits_data_t));

	cmc_read_fits_file(param.infile, cfd);
	
	/* add central BH, hidden in 0th star */
        if (param.bhmass > 0.0) {
		cfd->obj_m[0] = param.bhmass;
        }

	total_mass = set_masses(param, cfd);

        if (param.N_mass_seg!=0) {

          if (param.N_mass_seg<0) {
            param.N_mass_seg= calc_n_mass_seg_max(cfd);
            printf("Using N=%li for mass segregated cluster\n", param.N_mass_seg);
          } else {
            /* Check if NSTAR is large enough to produce a segr. cluster 
             * with N_mass_seg stars*/
            long nmax;
            nmax= calc_n_mass_seg_max(cfd);
            if (param.N_mass_seg> nmax) {
              printf("The requested number of stars (%li) is larger than", param.N_mass_seg);
              printf("the maximum (%li).\n", nmax);
              write_usage();
              exit(EXIT_FAILURE);
            }
          }

          ms= create_mass_seg(param, cfd);
	  cmc_free_fits_data_t(cfd);
          cfd= ms;
        }

	/* I might add scaling to E0 = -1/4 here */
        if (param.scale) {
		scale_pos_and_vel(param, cfd, total_mass);
        };
	
	cmc_write_fits_file(cfd, param.outfile);
	cmc_free_fits_data_t(cfd);

	return 0;
}

