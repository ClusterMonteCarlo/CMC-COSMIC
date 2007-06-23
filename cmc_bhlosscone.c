/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "cmc.h"
#include "cmc_vars.h"

void bh_rand_walk(long index, double v[4], double vcm[4], double beta, double dt)
{ 
	/*This is the random walk procedure as outlined by Freitag & Benz (2002).
	  Change of notation: beta here is theta in the paper. */
	double w[3], n_orb, P_orb, deltabeta_orb, L2, Rdisr, Jlc, vlc;
	double deltamax, deltasafe, delta, dbeta;
	double vt, vr, w_mag;
	int i;
#ifdef DEBUGGING
        FILE *rwalk_file;
        char fname[80];
        long is_in_ids;
        
        is_in_ids= 0;
        if (g_hash_table_lookup(star_ids, &star[index].id)!=NULL) {
          sprintf(fname, "%s.rwalk_steps_%.4i", outprefix, star[index].id);
          rwalk_file= fopen(fname, "a");
          fprintf(rwalk_file, "# Time deltabeta_orb deltasafe sqrt(L2) delta\n");
          is_in_ids=1;
        };
#endif
 	/* simulate loss cone physics for central mass */
	P_orb = calc_P_orb(index);
	n_orb = dt * ((double) clus.N_STAR)/log(GAMMA * ((double) clus.N_STAR)) / P_orb; 
	deltabeta_orb = 1.0/sqrt(n_orb) * beta; 
	L2 = fb_sqr(beta); 
        if (BH_R_DISRUPT_NB>0.) {
          Rdisr= BH_R_DISRUPT_NB;
	} else {
          Rdisr= pow(2.*cenma.m/star[index].m, 1./3.)*star[index].rad;
        };
	Jlc= sqrt(2.*cenma.m*madhoc*Rdisr);
	vlc= Jlc/star[index].r;
  	for (i=0; i<3; i++) {
    		w[i]= v[i+1]- vcm[i+1];
  	}
        w_mag= sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
	while (L2 > 0.0) { 
		if (sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(w[1]+vcm[2])) <= vlc) { 
			dprintf("index=%ld, id=%ld: star eaten by BH\n", index, star[index].id);
			cenma.m += star[index].m; 
			destroy_obj(index);
			L2 = 0.0; 
		} else { 
			deltamax= 0.1*FB_CONST_PI;
			deltasafe= CSAFE*(sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(vcm[2]+w[1]))-vlc)/w_mag;
			delta = MAX(deltabeta_orb, MIN(deltamax, MIN(deltasafe, sqrt(L2)))); 
			dbeta = 2.0 * PI * rng_t113_dbl(); 
			do_random_step(w, dbeta, delta); 
			L2 -= fb_sqr(delta); 
#ifdef DEBUGGING
                        if (is_in_ids) {
                          fprintf(rwalk_file, "%f %g %g %g %g\n", deltabeta_orb, deltasafe, sqrt(L2), delta);
                        }
#endif
		} 
	}; 
#ifdef DEBUGGING
        if (is_in_ids) {
          fclose(rwalk_file);
        }
#endif
}; 


void get_3d_velocities(double *w, double vr, double vt) {
   double phi;

   phi= rng_t113_dbl()*2.*FB_CONST_PI;
   w[0]= vt* cos(phi);
   w[1]= vt* sin(phi);
   w[2]= vr;
};

/* here the notation of Freitag & Benz (2002) is used */
void do_random_step(double *w, double beta, double delta) {
   double theta, phi, w_mag, new_w_dir[3];
   double a; /* this is a variable to store an intermediate result*/ 

   /* calculate direction of w*/
   w_mag= sqrt(fb_sqr(w[0])+ fb_sqr(w[1])+ fb_sqr(w[2]));
   theta= acos(w[2]/w_mag);
   phi= atan2(w[1], w[0]);
   
   /* rotate new vector (w_mag, beta, delta) into the direction of w */
   a= cos(theta)* sin(beta)* sin(delta)+ sin(theta)* cos(delta);
   new_w_dir[0]= sin(phi)* cos(beta) * sin(delta) + cos(phi)* a;
   new_w_dir[1]= -cos(phi)* cos(beta)* sin(delta) + sin(phi)* a;
   new_w_dir[2]= -sin(theta)* sin(beta)* sin(delta) + cos(theta)* cos(delta);

   w[0]= w_mag* new_w_dir[0];
   w[1]= w_mag* new_w_dir[1];
   w[2]= w_mag* new_w_dir[2];
};

double check_angle_w_w_new(double *w, double *w_new, double delta) {
   /* calculate the angle between w and w_new and compare to deltat */
   double angle;

   angle=acos(w[0]*w_new[0]+w[1]*w_new[1]+w[2]*w_new[2]);

   return(angle-delta);
};

/* calculate star's radial orbital period */
double calc_P_orb(long index)
{
	double E, J, Porb, error, Porbapproxmin, Porbapproxmax, Porbapprox;
	//double Porbtmp;
	orbit_rs_t orbit_rs;
	calc_p_orb_params_t params;
	gsl_integration_workspace *w;
	gsl_integration_qaws_table *tab;
	gsl_function F;
        struct Interval star_interval;
        
        /* default values for star_interval */
	star_interval.min= 1;
        star_interval.max= clus.N_MAX+1;

	E = star[index].E + PHI_S(star[index].r, index);
	J = star[index].J;
	
	//dprintf("index=%ld ", index);

	orbit_rs = calc_orbit_rs(index, E, J);
	
	//dprintf("rp=%g ra=%g ", orbit_rs.rp, orbit_rs.ra);

	Porbapproxmin = 2.0 * PI * orbit_rs.rp / (J / orbit_rs.rp);
	Porbapproxmax = 2.0 * PI * orbit_rs.ra / (J / orbit_rs.ra);
	Porbapprox = sqrt(Porbapproxmin * Porbapproxmax);

	/* Return the approximate value of the radial period.  If this is commented out
	   the radial period will be calculated properly. */
	//return(Porbapprox);

	if (orbit_rs.circular_flag == 1) {
		/* We're returning the azimuthal period here, which is not the same as the
		   radial period for the general cluster potential.  This shouldn't make
		   any difference for the BH loss cone stuff, since the orbit is circular. */
		return(2.0 * PI * star[index].r / star[index].vt);		
	} else {
		w = gsl_integration_workspace_alloc(1000);
		tab = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0.0, 0.0);

		params.E = E;
		params.J = J;
		params.index = index;
                if (SEARCH_GRID)
                  star_interval= search_grid_get_interval(r_grid, orbit_rs.rp);
		params.kmin = FindZero_r(star_interval.min, star_interval.max, orbit_rs.rp);
                if (SEARCH_GRID)
                  star_interval= search_grid_get_interval(r_grid, orbit_rs.ra);
		params.kmax = FindZero_r(star_interval.min, star_interval.max, orbit_rs.ra) + 1;
		params.rp = orbit_rs.rp;
		params.ra = orbit_rs.ra;
		F.params = &params;

		if (0) { /* use standard potential function with Stefan's speedup trick here */
			F.function = &calc_p_orb_f;
			gsl_integration_qags(&F, orbit_rs.rp, orbit_rs.ra, 0, 1.0e-3, 1000, w, &Porb, &error);
			//dprintf("Porb=%g Porb/Porbapprox=%g intervals=%d\n", Porb, Porb/Porbapprox, w->size);
		}

		if (0) { /* use fast potential function here (not much of a speedup over Stefan's technique in practice) */
			F.function = &calc_p_orb_f2;
			gsl_integration_qags(&F, orbit_rs.rp, orbit_rs.ra, 0, 1.0e-3, 1000, w, &Porb, &error);
			//dprintf("\tFast: Porb=%g Porb/Porbapprox=%g intervals=%d\n", Porb, Porb/Porbapprox, w->size);
		}

		if (1) { /* use Gauss-Chebyshev for factor of ~few speedup over standard method */
			//Porbtmp = Porb;
			F.function = &calc_p_orb_gc;
			gsl_integration_qaws(&F, orbit_rs.rp, orbit_rs.ra, tab, 1.0e-3, 1.0e-3, 1000, w, &Porb, &error);
			//dprintf("Porb=%g Porb/Porbtmp=%g Porb/Porbapprox=%g intervals=%d\n", 
			//	Porb, Porb/Porbtmp, Porb/Porbapprox, w->size);
		}
		
		gsl_integration_qaws_table_free(tab);
		gsl_integration_workspace_free(w);
		return(Porb);
	}
}

/* integrand for calc_P_orb */
double calc_p_orb_f(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (potential(x) + PHI_S(x, myparams.index));
	
	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/* integrand for calc_P_orb */
double calc_p_orb_f2(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (fastpotential(x, myparams.kmin, myparams.kmax) + PHI_S(x, myparams.index));
	
	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/* integrand for calc_P_orb, using Gauss-Chebyshev for regularizing the integrand near the endpoints */
double calc_p_orb_gc(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;
	double phik, phik1, phi0, phi1, rk, rk1, rminus, rplus;
	double E, J, rp, ra;
	long index, kmin, kmax;
	
	E = myparams.E;
	J = myparams.J;
	index = myparams.index;
	kmin = myparams.kmin;
	kmax = myparams.kmax;
	rp = myparams.rp;
	ra = myparams.ra;

	if (x < star[kmin+1].r) { /* return integrand regularized at r=rp */
		//dprintf("regularizing near rp...\n");
		phik = star[kmin].phi + PHI_S(star[kmin].r, index);
		phik1 = star[kmin+1].phi + PHI_S(star[kmin+1].r, index);
		rk = star[kmin].r;
		rk1 = star[kmin+1].r;
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1);
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1);
		/* rplus = rperi, rminus = rapo-primed */
		rminus = (phi1 - sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//rplus = (phi1 + sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//dprintf("rplus/rp=%g rminus/ra=%g\n", rplus/rp, rminus/ra);
		if (kmax == kmin + 1) {
			/* then rminus = ra, so must cancel (ra-x)/(rminus-x) term analytically */
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
			return(2.0*x*sqrt((ra-x)/((2.0*phi0-2.0*E)*(rminus-x))));
		}
	} else if (x > star[kmax-1].r) { /* return integrand regularized at r=ra*/
		//dprintf("regularizing near ra...\n");
		phik = star[kmax-1].phi + PHI_S(star[kmax-1].r, index);
		phik1 = star[kmax].phi + PHI_S(star[kmax].r, index);
		rk = star[kmax-1].r;
		rk1 = star[kmax].r;
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1);
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1);
		/* rplus = rperi-primed, rminus = rapo */
		//rminus = (phi1 - sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		rplus = (phi1 + sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//dprintf("rplus/rp=%g rminus/ra=%g\n", rplus/rp, rminus/ra);
		if (kmax == kmin + 1) {
			/* then rplus = rp, so must cancel (x-rp)/(x-rplus) term analytically */
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
			return(2.0*x*sqrt((x-rp)/((2.0*phi0-2.0*E)*(x-rplus))));
		}
	} else {
		radicand = 2.0 * E - fb_sqr(J/x) - 2.0 * (potential(x) + PHI_S(x, index));
		if (radicand < 0.0) {
			dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, index);
			radicand = 0.0;
		}
		return(2.0 * sqrt((x-rp)*(ra-x)/radicand));
	}
}
