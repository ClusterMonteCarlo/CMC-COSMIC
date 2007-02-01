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

void bh_rand_walk(long index, double beta, double dt)
{ 
	/*This is the random walk procedure as outlined by Freitag & Benz (2002).
	  Change of notation: beta here is theta in the paper. */
	double w[3], n_orb, P_orb, deltabeta_orb, L2, Rdisr, Jlc, vlc;
	double deltamax, deltasafe, delta, dbeta;
	double vt, vr;
	
 	/* simulate loss cone physics for central mass */
	vt = star[index].vt;
	vr = star[index].vr;
	P_orb = calc_P_orb(index);
	n_orb = dt * ((double) clus.N_STAR)/log(GAMMA * ((double) clus.N_STAR)) / P_orb; 
	deltabeta_orb = 1.0/sqrt(n_orb) * beta; 
	L2 = fb_sqr(beta); 
	Rdisr= pow(2.*cenma.m/star[index].m, 1./3.)*star[index].rad;
	Jlc= sqrt(2.*cenma.m*madhoc*Rdisr);
	vlc= Jlc/star[index].r;
	get_3d_velocities(w, vr, vt);
	while (L2 > 0.0) { 
		if (sqrt(fb_sqr(w[0])+fb_sqr(w[1])) <= vlc) { 
			dprintf("index=%ld: star eaten by BH\n", index);
			cenma.m += star[index].m; 
			destroy_obj(index);
			L2 = 0.0; 
		} else { 
			deltamax= 0.1*FB_CONST_PI;
			deltasafe= CSAFE*(sqrt(fb_sqr(w[0])+fb_sqr(w[1]))-vlc)/sqrt(vt*vt+vr*vr);
			delta = MAX(deltabeta_orb, MIN(deltamax, MIN(deltasafe, sqrt(L2)))); 
			dbeta = 2.0 * PI * rng_t113_dbl(); 
			do_random_step(w, dbeta, delta); 
			L2 -= fb_sqr(delta); 
		} 
	}; 
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
	orbit_rs_t orbit_rs;
	calc_p_orb_params_t params;
	gsl_integration_workspace *w;
	gsl_function F;
	
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
	return(Porbapprox);

	if (orbit_rs.circular_flag == 1) {
		/* We're returning the azimuthal period here, which is not the same as the
		   radial period for the general cluster potential.  This shouldn't make
		   any difference for the BH loss cone stuff, since the orbit is circular. */
		return(2.0 * PI * star[index].r / star[index].vt);		
	} else {
		w = gsl_integration_workspace_alloc(1000);

		params.E = E;
		params.J = J;
		params.index = index;

		F.function = &calc_p_orb_f;
		F.params = &params;

		gsl_integration_qags(&F, orbit_rs.rp, orbit_rs.ra, 0, 1e-3, 1000, w, &Porb, &error);

		//dprintf("Porb=%g Porb/Porbapprox=%g\n", Porb, Porb/Porbapprox);

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
