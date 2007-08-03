#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "../../common/fitslib.h"
#include "../../common/taus113-v2.h"

#define SEED 768364873UL

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
	(param->infile)[0] = '\0';
	(param->outfile)[0] = '\0';

	while((c = getopt(argc, argv, "i:o:I:wm:M:p:u:N:r:C:hsb:")) != -1)
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
        //total_mass+= param.bhmass;
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
int main(int argc, char *argv[]){
	cmc_fits_data_t cfd;
	struct imf_param param;
	double total_mass;
	
	parse_options(&param, argc, argv);

	cmc_read_fits_file(param.infile, &cfd);
	
        if (param.bhmass > 0.0) {
		cfd.obj_m[0] = param.bhmass;
        }

	total_mass = set_masses(param, &cfd);

	/* I might add scaling to E0 = -1/4 here */
        if (param.scale) {
		scale_pos_and_vel(param, &cfd, total_mass);
        };

	cmc_write_fits_file(&cfd, param.outfile);

	cmc_free_fits_data_t(&cfd);

	return 0;
}

