/* -*- linux-c -*- */
/* fewbody_coll.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "fewbody.h"

/* the main collision criterion */
int fb_is_collision(double r, double R1, double R2)
{
	if (r < R1 + R2) {
		return(1);
	} else {
		return(0);
	}
}

// PAU int fb_collide(fb_hier_t *hier, double f_exp)
int fb_collide(fb_hier_t *hier, double f_exp, fb_units_t units, gsl_rng *rng, struct rng_t113_state *curr_st, double bh_reff)
{
	int i, j=-1, k, retval=0, cont=1;
	double R[3], peinit;

	/* this is a non-recursive way to perform a recursive operation: keep going until there are no more
	   mergers */
	while (cont) {
		cont = 0;
		for (i=0; i<hier->nstar-1; i++) {
			for (j=i+1; j<hier->nstar; j++) {
				/* calculate relative separation */
				for (k=0; k<3; k++) {
					R[k] = hier->hier[hier->hi[1]+i].x[k] - hier->hier[hier->hi[1]+j].x[k];
				}

				/* test collision criterion */
				if (fb_is_collision(fb_mod(R), hier->hier[hier->hi[1]+i].R, hier->hier[hier->hi[1]+j].R)) {
					cont = 1;
					/* break out of the double loop if there is a collision, so we can merge
					   the stars immediately */
					break;
				}
			}
			/* break out of the double loop if there is a collision, so we can merge
			   the stars immediately */
			if (cont) {
				break;
			}
		}
		
		/* merge two stars if necessary */
		if (cont) {
			/* return 1 if there is a collision */
			retval = 1;

			/* calculate the potential energy before the collision, since we're going to need
			   to account for the change in potential energy, and put it in Eint */
			peinit = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);			
			
			/* do the actual merger */
			fb_dprintf("fewbody: collide(): merging stars: i=%d j=%d\n", i, j);
			// PAU fb_merge(&(hier->hier[hier->hi[1]+i]), &(hier->hier[hier->hi[1]+j]), hier->nstarinit, f_exp);
			fb_merge(&(hier->hier[hier->hi[1]+i]), &(hier->hier[hier->hi[1]+j]), hier->nstarinit, f_exp, units, rng, curr_st, bh_reff);
			fb_objcpy(&(hier->hier[hier->hi[1]+j]), &(hier->hier[hier->hi[1]+hier->nstar-1]));
			hier->nstar--;

			/* calculate the difference in potential energy before and after the collision, and put it
			   in Eint for accounting */
			hier->hier[hier->hi[1]+i].Eint += peinit - fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);
		}
	}

	return(retval);
}

// PAU void fb_merge(fb_obj_t *obj1, fb_obj_t *obj2, int nstarinit, double f_exp)
void fb_merge(fb_obj_t *obj1, fb_obj_t *obj2, int nstarinit, double f_exp, fb_units_t units, gsl_rng *rng, struct rng_t113_state *curr_st, double bh_reff)
{
	int i;
	double x1[3], x2[3], v1[3], v2[3], l1[3], l2[3], A[3], L[3], l[3];
    double E;
	double clight, xrel[3], vrel[3], x[3], y[3], z[3], theta;
	double v_perp, v_para, afinal, mass_frac;
	fb_obj_t tmpobj;
	clight = FB_CONST_C / units.v;

	/* sanity check */
	if (obj1->n != 1 || obj2->n != 1) {
		fprintf(stderr, "fb_merge: trying to merge an object that isn't single?!\n");
		exit(1);
	}

	/* need some temporary storage here */
	tmpobj.id = (long *) malloc(nstarinit * sizeof(long));
	tmpobj.vkick = (double *) malloc(nstarinit * sizeof(double));
	tmpobj.a_merger = (double *) malloc(nstarinit * sizeof(double));
	tmpobj.e_merger = (double *) malloc(nstarinit * sizeof(double));

	/* merge id's */
	tmpobj.ncoll = obj1->ncoll + obj2->ncoll;
	for (i=0; i<obj1->ncoll; i++) {
		tmpobj.id[i] = obj1->id[i];
		tmpobj.vkick[i] = obj1->vkick[i];
		tmpobj.a_merger[i] = obj1->a_merger[i];
		tmpobj.e_merger[i] = obj1->e_merger[i];
	}
	for (i=0; i<obj2->ncoll; i++) {
		tmpobj.id[obj1->ncoll + i] = obj2->id[i];
		tmpobj.vkick[obj1->ncoll + i] = obj2->vkick[i];
		tmpobj.a_merger[obj1->ncoll + i] = obj2->a_merger[i];
		tmpobj.e_merger[obj1->ncoll + i] = obj2->e_merger[i];
	}
    tmpobj.vkick[tmpobj.ncoll-1] = 0; 

	/* create idstring */
	snprintf(tmpobj.idstring, FB_MAX_STRING_LENGTH, "%s:%s", obj1->idstring, obj2->idstring);

	/* assume no mass loss */
	/* If a BH merger, we'll apply g-wave mass loss below*/
	tmpobj.m = obj1->m + obj2->m;

	/* this is just a simple prescription */
	tmpobj.R = f_exp * (obj1->R + obj2->R);

    /* this will either get reset below or by CMC */
    tmpobj.k_type = -1;

	/* set new position and velocity, calculate relative positions and velocities */
	for (i=0; i<3; i++) {
		tmpobj.x[i] = (obj1->m * obj1->x[i] + obj2->m * obj2->x[i]) / tmpobj.m;
		tmpobj.v[i] = (obj1->m * obj1->v[i] + obj2->m * obj2->v[i]) / tmpobj.m;
		x1[i] = obj1->x[i] - tmpobj.x[i];
		x2[i] = obj2->x[i] - tmpobj.x[i];
		v1[i] = obj1->v[i] - tmpobj.v[i];
		v2[i] = obj2->v[i] - tmpobj.v[i];
        xrel[i] = obj1->x[i] - obj2->x[i];
        vrel[i] = obj1->v[i] - obj2->v[i];
	}

    /* compute the orbital energy and semi-major axis at merger*/
    E = 0.5 * (obj1->m * fb_dot(obj1->v, obj1->v) + obj2->m * fb_dot(obj2->v, obj2->v)) - 
        obj1->m * obj2->m/fb_mod(xrel);
	tmpobj.a_merger[tmpobj.ncoll-1] = -obj1->m * obj2->m / (2.0 * E);

	/* set internal energy, using the difference in kinetic energy; the difference in potential energy
	   depends on the positions of the other stars, and will be calculated later and added to Eint */
	tmpobj.Eint = obj1->Eint + obj2->Eint +
		0.5 * (obj1->m * fb_dot(obj1->v, obj1->v) + obj2->m * fb_dot(obj2->v, obj2->v)) -
		0.5 * tmpobj.m * fb_dot(tmpobj.v, tmpobj.v);

	/* set internal angular momentum */
	fb_cross(x1, v1, l1);
	fb_cross(x2, v2, l2);
	for (i=0; i<3; i++) {
		tmpobj.Lint[i] = obj1->Lint[i] + obj2->Lint[i] + obj1->m * l1[i] + obj2->m * l2[i];
	}

    /* compute angular momenta for LRL vector*/
	for (i=0; i<3; i++) {
		L[i] = obj1->m * l1[i] + obj2->m * l2[i];
		l[i] = L[i] * (obj1->m + obj2->m)/(obj1->m * obj2->m);
	}
	
	/* -A = l x v + G M \hat r */
	fb_cross(vrel, l, A);
	for (i=0; i<3; i++) {
		A[i] -= (obj1->m + obj2->m) * xrel[i]/fb_mod(xrel);
	}
	
    /* magnitude of A gives the eccentricity at merger*/
	tmpobj.e_merger[tmpobj.ncoll-1] = fb_mod(A)/(obj1->m + obj2->m);
	

	/* Apply a change in mass/spin/recoil speed for compact-object mergers */
	if (obj1->k_type == 14 && obj2->k_type == 14){

		/* first calculate relative velocity, which we'll use for the y-vector */
		for (i=0; i<3; i++) {
			vrel[i] = obj1->v[i] - obj2->v[i];
		}
		/* the z-vector is the unit vector pointing in the direction of the angular momentum */
		for (i=0; i<3; i++) {
			y[i] = vrel[i]/fb_mod(vrel);
			z[i] = tmpobj.Lint[i]/fb_mod(tmpobj.Lint);
		}
		/* x equals y cross z */
		fb_cross(y, z, x);

		/* we randomize the angle in the orbital plane */
		theta = 2.0 * FB_CONST_PI * rng_t113_dbl_new(curr_st);

		/*Compute the actual merger; returns the final spin, change in mass, and the
		 * in-plane and z-axis components of the kick*/
		fb_bh_merger(obj1->m, obj2->m, obj1->chi, obj2->chi, &mass_frac, &afinal, &v_para, &v_perp, curr_st);

		/* Add the kick to our coordinate system*/
		for (i=0; i<3; i++) {
            // CGS units for the velocity divided by the critical velocity
			tmpobj.v[i] += 1.e5*v_perp/units.v * (cos(theta) * x[i] + sin(theta) * y[i]);
			tmpobj.v[i] += 1.e5*v_para/units.v * z[i];
		}

		/* Update the spin and mass of the merger product */
		tmpobj.chi = afinal;
        tmpobj.k_type = 14;
		tmpobj.m *= mass_frac;
        tmpobj.vkick[tmpobj.ncoll-1] = sqrt(v_perp*v_perp + v_para*v_para);

        /* Also update the radius with the new ISCO radius */ 
        /* We could use the Kerr ISCO here instead of the Schwarzschild,
         * but any case where that's relavant will almost always merge...*/
        tmpobj.R = bh_reff*2*(tmpobj.m * units.m)*FB_CONST_G / FB_CONST_C / FB_CONST_C / units.l;

    /* Note: BSE will be called for collision, but only AFTER the fewbody integration is complete
     * On the off chance we have a repeated merger, we need to know what to do now.  For now,
     * any BBH merger becomes a new BH, and the BH/non-BH also becomes BH, with the same spin as 
     * the progenitor BH  
	 *
	 * Note that it does NOT allow for any accretion; the star is simply
	 * destroyed
     *
     * TODO: Josh, Kyle, or possibly Future Carl: if you ever want to consider BH spin-up during
     * mergers, you'll have to change this as well as the BSE merger matrix*/
	} else if (obj1->k_type == 14){
        tmpobj.chi = obj1->chi;
		tmpobj.m = obj1->m;
        tmpobj.k_type = 14;
        tmpobj.R = bh_reff*2*(tmpobj.m * units.m)*FB_CONST_G / FB_CONST_C / FB_CONST_C / units.l;
    } else if (obj2->k_type == 14){
        tmpobj.chi = obj2->chi;
		tmpobj.m = obj2->m;
        tmpobj.k_type = 14;
        tmpobj.R = bh_reff*2*(tmpobj.m * units.m)*FB_CONST_G / FB_CONST_C / FB_CONST_C / units.l;
    }

	/* and better set these, too... */
	tmpobj.n = 1;
	
	tmpobj.obj[0] = NULL;
	tmpobj.obj[1] = NULL;

	/* finally, copy over the merger from temporary storage */
	fb_objcpy(obj1, &tmpobj);

	free(tmpobj.id);
	free(tmpobj.vkick);
	free(tmpobj.a_merger);
	free(tmpobj.e_merger);
}

/* How to deal with a merger: takes as input m1,m2 (in any units), a1, a2
 * (dimmensionless Kerr parameter).
 *
 * Returns the mass fraction (i.e. M_final/M_initial), the final spin
 * (dimmensionless), and the recoil kick (in km/s)*/ 
void fb_bh_merger(double m1, double m2, double a1, double a2, 
				  double *mass_frac, double *afinal, 
				  double *v_para, double *v_perp, struct rng_t113_state *curr_st)
{

	double delta_par, delta_perp, chi_par, chi_perp;
	double z1,z2,rISCO,eISCO;
	/*Keeping this seperate in case we want to implement non-uniform spins
	 * at some point*/
	double vm, vs_perp, vs_par, vk;
	double lil_L;
	double eta = m2*m1 / pow(m2+m1,2);
	double theta1,theta2,phi1,phi2,Theta;
    double delta_x,delta_y,delta_z,chi_x,chi_y,chi_z;
	double X;
    double q,chi1,chi2;

    if (m1>m2){
        q = m2/m1;
        chi1 = a1;
        chi2 = a2;
    } else{
        q = m1/m2;
        chi1 = a2;
        chi2 = a1;
    }

	X = rng_t113_dbl_new(curr_st);
	theta1 = acos(2*X - 1.);
	X = rng_t113_dbl_new(curr_st);
	theta2 = acos(2*X - 1.);
	X = rng_t113_dbl_new(curr_st);
    phi1 = X * 2 * FB_CONST_PI; 
	X = rng_t113_dbl_new(curr_st);
    phi2 = X * 2 * FB_CONST_PI; 
	X = rng_t113_dbl_new(curr_st);
	Theta = X * FB_CONST_PI; 

    /*Compute the approprite mass-weighted spin combinations and their
	 * projections parallel and perpendicular to L*/
    delta_x = (q*chi2*sin(theta2)*cos(phi2) - chi1*sin(theta1)*cos(phi1)) /(1+q);
    delta_y = (q*chi2*sin(theta2)*sin(phi2) - chi1*sin(theta1)*sin(phi1)) /(1+q);
    delta_z = (q*chi2*cos(theta2) - chi1*cos(theta1)) /(1+q);
               
    chi_x = (q*q*chi2*sin(theta2)*cos(phi2) + chi1*sin(theta1)*cos(phi1)) /pow(1+q,2);
    chi_y = (q*q*chi2*sin(theta2)*sin(phi2) + chi1*sin(theta1)*sin(phi1)) /pow(1+q,2);
    chi_z = (q*q*chi2*cos(theta2) + chi1*cos(theta1)) /pow(1+q,2);

    delta_par = delta_z;
    chi_par = chi_z;
    delta_perp = sqrt(delta_x*delta_x + delta_y*delta_y);
    chi_perp = sqrt(chi_x*chi_x + chi_y*chi_y);

    /*The compute the energy-per-mass at the Kerr ISCO of an effective particle
	 * with that spin*/
	z1 = 1 + pow(1-chi_par*chi_par,0.3333333333)*(pow(1+chi_par,0.3333333333)+pow(1-chi_par,0.3333333333));
	z2 = sqrt(3*chi_par*chi_par + z1*z1);
	rISCO = 3 + z2 - copysignf(sqrt((3-z1)*(3+z1+2*z2)),chi_par);
	eISCO = sqrt(1 - 2. / 3. / rISCO);

	/*The final remnant mass is, based on the extrapolation between equal-mass
	 * and test mass limits (Barausse et al, Apj, 758, 63 (2012))*/
	*mass_frac = 1 - eta*(1-4*eta)*(1-eISCO) - 16*eta*eta*(0.04827 + 4*0.01707*chi_par*(chi_par + 1));
	/*The final spin is taken from Barausse et al. Apj, 704, L40 (2009)
	 * Note that I've assumed L-hat lies along the Z-axis*/
	lil_L = 3.464102 - 3.51712*eta + 2.5763*eta*eta - 0.1229*pow(1+q,4.)*(chi_par*chi_par + chi_perp*chi_perp)
	          / pow(1+q*q,2.) + (0.4537*eta - 2.8904 + 2)*pow(1+q,2)*chi_par / (1+q*q);
	*afinal = FB_MIN(1.,fabs(((q*q*chi2*cos(theta2) + chi1*cos(theta1))/pow(1+q,2)) + q * lil_L / pow(1+q,2.)));

	/*Then compute velocity of the recoil kick, both from the assymetric mass
	 * ratio and the misalignment of the spins.  These fits to NR simulations
	 * are taken from:
	 *
	 *	Campanelli et al, APJ 659, L5 (2007)
	 *	Gonzalez et al, PRL, 98, 091101 (2007)
	 *	Lousto et al, PRD, 77, 044028 (2008)
	 *	Lousto et al, PRD, 85, 084015 (2012)
	 *	Lousto and Zlochower, PRD, 87, 084027 (2013)
	 *
	 *	Though I just copied them all from the Gerosa paper...
	 * */
	vm = 1.2e4 * eta*eta * ((1-q) / (1+q)) * (1 - 0.93*eta);
	vs_perp = 6.9e3 * eta*eta * delta_par;
	vs_par = 16.*eta*eta *(delta_perp *(3677.76 + 2*2481.21*chi_par + 4*1792.45*chi_par*chi_par + 
	         8*1506.52*pow(chi_par,3.)) + 2.*chi_perp*delta_par*(1140. + 2*2481.*chi_par))*cos(Theta);
	*v_para = sqrt(vs_par*vs_par);
	*v_perp = sqrt(vm*vm - 2*vm*vs_perp*0.81915 + vs_perp*vs_perp);
}
