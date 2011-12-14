/* -*- linux-c -*- */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cmc.h"
#include "cmc_vars.h"

void zero_star(long j)
{
	star[j].r = 0.0;
	star[j].vr = 0.0;
	star[j].vt = 0.0;
	star[j].m = 0.0;
	star[j].E = 0.0;
	star[j].J = 0.0;
	star[j].EI = 0.0;
	star[j].Eint = 0.0;
	star[j].rnew = 0.0;
	star[j].vrnew = 0.0;
	star[j].vtnew = 0.0;
	star[j].rOld = 0.0;
	star[j].X = 0.0;
	star[j].Y = 0.0;
	star[j].r_peri = 0.0;
	star[j].r_apo = 0.0;
	star[j].phi = 0.0;
	star[j].interacted = 0;
	star[j].binind = 0;
	star[j].id = 0;
	star[j].rad = 0.0;
	star[j].Uoldrold = 0.0;
	star[j].Uoldrnew = 0.0;
	//Sourav: toy rejuvenation- some more variables
	star[j].createtime = 0.0;
	star[j].lifetime = GSL_POSINF;
	star[j].vtold = 0.0;
	star[j].vrold = 0.0;
	star[j].se_mass = 0.0;
	star[j].se_k = 0;
	star[j].se_mt = 0.0;
	star[j].se_ospin = 0.0;
	star[j].se_epoch = 0.0;
	star[j].se_tphys = 0.0;
	star[j].se_radius = 0.0;
	star[j].se_lum = 0.0;
	star[j].se_mc = 0.0;
	star[j].se_rc = 0.0;
	star[j].se_menv = 0.0;
	star[j].se_renv = 0.0;
	star[j].se_tms = 0.0;
}

void zero_binary(long j)
{
	binary[j].id1 = 0;
	binary[j].id2 = 0;
	binary[j].rad1 = 0.0;
	binary[j].rad2 = 0.0;
	binary[j].m1 = 0.0;
	binary[j].m2 = 0.0;
	binary[j].Eint1 = 0.0;
	binary[j].Eint2 = 0.0;
	binary[j].a = 0.0;
	binary[j].e = 0.0;
	binary[j].inuse = 0;
	//Sourav: toy rejuvenation- some new variables
	binary[j].createtime_m1 = 0.0;
	binary[j].createtime_m2 = 0.0;
	binary[j].lifetime_m1 = GSL_POSINF;
	binary[j].lifetime_m2 = GSL_POSINF;
	// Meagan: to keep track of three-body binaries
	//binary[j].threebodybinary = 0;
}

void print_interaction_status(char status_text[])
{
	gprintf(" %s", status_text);
	if (!quiet) {
		fflush(stdout);
	}
	fprintf(logfile, " %s", status_text);
}

void print_interaction_error(void)
{
	gprintf("!");
	if (!quiet) {
		fflush(stdout);
	}
	fprintf(logfile, "!");
}

/* destroy an object (can be star or binary) */
void destroy_obj(long i)
{
	double r, phi;

	if (star[i].binind) {
		destroy_binary(star[i].binind);
	}
	
	/* need to zero out E's, J, but can't zero out potential---this is the easiest way */
	r = star[i].r;
	phi = star[i].phi;
	zero_star(i);
	star[i].r = r;
	star[i].phi = phi;

	remove_star_center(i);
}

/* destroy a binary */
void destroy_binary(long i)
{
	/* set inuse flag to zero, and zero out all other properties for safety */
	zero_binary(i);
}

/* create a new star, returning its index */
long create_star(void)
{
	long i;
	
	/* account for new star */
	clus.N_STAR_NEW++;
	clus.N_MAX_NEW++;
	
	/* put new star at end; the +1 is to not overwrite the boundary star */
	i = clus.N_MAX_NEW + 1;

	/* initialize to zero for safety */
	zero_star(i);
	
	return(i);
}

/* create a new binary, returning its index */
long create_binary(void)
{
	long i, j;
	
	/* find first free binary */
	i = 1;
	while (i <= N_BIN_DIM && binary[i].inuse) {
		i++;
	}
	
	/* problem! */
	if (i > N_BIN_DIM) {
		eprintf("cannot find unused binary.\n");
		exit_cleanly(1);
	}
	
	/* initialize to zero for safety */
	zero_binary(i);

	/* mark binary as being in use */
	binary[i].inuse = 1;

	/* create the star that points to the binary */
	j = create_star();
	
	star[j].binind = i;
	
	return(j);
}


void sort_three_masses(long sq, long *k1, long *k2, long *k3) {
/* Sort by mass: k1, k2, k3 in order of most massive to least massive. */

	//fprintf(threebbfile, "beginning of sort_three_masses function" );
	if ((star[sq].m >= star[sq+1].m) && (star[sq].m >= star[sq + 2].m)) {
		*k1 = sq;
		if (star[sq + 1].m >= star[sq + 2].m) {
			*k2 = sq + 1;
			*k3 = sq + 2;
		}
		else {
			*k2 = sq + 2;
			*k3 = sq + 1;
		}
	}
	else if ((star[sq + 1].m >= star[sq].m) && (star[sq + 1].m >= star[sq + 2].m)) {
		*k1 = sq + 1;
		if (star[sq].m >= star[sq + 2].m) {
			*k2 = sq;
			*k3 = sq + 2;
		}
		else {
			*k2 = sq + 2;
			*k3 = sq;
		}
	}
	else if ((star[sq + 2].m >= star[sq].m) && (star[sq + 2].m >= star[sq + 1].m)) {
		*k1 = sq + 2;
		if (star[sq].m > star[sq + 1].m) {
			*k2 = sq;
			*k3 = sq + 1;
		}
		else {
			*k2 = sq + 1;
			*k3 = sq;

		}
	}
	//fprintf(threebbfile, "end of sort_three_masses function" );
}


void calc_3bb_encounter_dyns(long k1, long k2, long k3, double angle1, double angle2, double v1[4], double v2[4], double v3[4], double (*vrel12)[4], double (*vrel3)[4], gsl_rng *rng) {
	/* Function for calculating quantities relevant for three-body binary formation
		Generates random direction for vt for all stars (so then have 3D veloc. for each)
		Then, for the two possible binaries that could form (1,2 or 1,3) we calculate:
		(1) COM veloc of candidate binary pair  (2) v3: veloc of 3rd star relative to binary
		In COM frame of the system of three stars, calculate total momentum and energy

	*/
	long j;
	double vcm12[4];

	/* set random angle between vt's */
	/* first store random variable */
//	star[k].Y = rng_t113_dbl();
//	star[kp].Y = star[k].Y;
	angle1 = rng_t113_dbl() * 2.0 * PI;
	angle2 = rng_t113_dbl() * 2.0 * PI;

//	phi = star[k].Y * 2.0 * PI;
	
	// Set velocities of three stars
	// Set velocities of three stars
	v1[1] = star[k1].vt;
	v1[2] = 0.0;
	v1[3] = star[k1].vr;
	v2[1] = star[k2].vt * cos(angle1);
	v2[2] = star[k2].vt * sin(angle1);
	v2[3] = star[k2].vr;	
	v3[1] = star[k3].vt * cos(angle2);
	v3[2] = star[k3].vt * sin(angle2);
	v3[3] = star[k3].vr;

	v1[0] = sqrt(sqr(v1[1]) + sqr(v1[2]) + sqr(v1[3]));
	v2[0] = sqrt(sqr(v2[1]) + sqr(v2[2]) + sqr(v2[3]));
	v3[0] = sqrt(sqr(v3[1]) + sqr(v3[2]) + sqr(v3[3]));

	for (j=1; j<=3; j++) {

		// Quantities needed for calculating 3bb formation rate:
		// relative veloc between the two candidates stars for binary
		// COM velocity of candidate binary pair
		// velocity of third star wrt com of binary pair

		vcm12[j] = (star[k1].m * v1[j] + star[k2].m * v2[j])/(star[k1].m + star[k2].m);
		(*vrel12)[j] = v2[j] - v1[j];
		(*vrel3)[j] = v3[j] - vcm12[j];
		//vrel12[j] = v2[j] - v1[j];
		//vrel3[j] = v3[j] - vcm12[j];
	}
		// Calculate magnitudes
		(*vrel12)[0] = sqrt(sqr((*vrel12)[1]) + sqr((*vrel12)[2]) + sqr((*vrel12)[3]));
		(*vrel3)[0] = sqrt(sqr((*vrel3)[1]) + sqr((*vrel3)[2]) + sqr((*vrel3)[3]));
		//vrel12[0] = sqrt(sqr(vrel12[1]) + sqr(vrel12[2]) + sqr(vrel12[3]));
		//vcm12[0] = sqrt(sqr(vcm12[1]) + sqr(vcm12[2]) + sqr(vcm12[3]));
		//vrel3[0] = sqrt(sqr(vrel3[1]) + sqr(vrel3[2]) + sqr(vrel3[3]));

}


double get_eta(double eta_min, long k1, long k2, long k3, double vrel12[4], double vrel3[4])
{
	double kk, eta, comp_value, norm, area_max, area, eta_test, comp_ymax, comp_y, true_y;
	double eta_max=50.0;
	long found_eta;
	// choose eccentricity of new binary from thermal distribution: f(e)=2e; 
	// e can range from zero to e_max, and for now set emax to 1; 
	// note when we set up initial binaries, have the option to truncate the 
	// e distribution near contact, which sets emax to something less than 1.
	// Since we're most likely to form BH binaries, their radii are so small 
	// that large eccentricities are possible, and maybe we don't have to worry about this?

	// Tested: sampled from distribution and graphed - the distribution is properly sampled

	// start by setting comp_value to something that will be rejected (above distribution for all eta)
	kk = 2.0 * ( (star[k1].m + star[k2].m + star[k3].m) / (star[k1].m + star[k2].m) ) * sqr(vrel12[0] / vrel3[0]);
	comp_y=1e6;
	found_eta = 0;

	/* normalize distribution, d(Rate)/d(eta) (differential rate of 3bb formation); take abs. value*/
	norm = fabs(1.0/(-1*(pow(eta_max, -5.5)) + (pow(eta_min, -5.5)) + (kk + 2) * ( -1*(pow(eta_max, -4.5)) + (pow(eta_min, -4.5)) ) + 2 * kk * ( -1*(pow(eta_max, -3.5)) + (pow(eta_min, -3.5)))));
	/* now start a while loop that will keep restarting as long as the last computed comp_value was rejected */
	while (found_eta == 0) {

		/* total area under comparison curve in range of eta from eta_min to eta_max */
	//	area_max = fabs(3.25*norm*((pow(eta_max,-4.0)) - (pow(eta_min, -4.0))));
		area_max = fabs((2.0/7.0) * norm * kk * ((pow(eta_max,-3.5)) - (pow(eta_min, -3.5))));
		/* choose area in range (0, area_max); this area corresponds to a particular value of eta */
		area = rng_t113_dbl()*area_max;
		/* find the eta that corresponds to chosen area */
		/* test coordinate 1 (x) */
	//	eta_test = pow(((-2.0*area/(13.0*norm)) + pow(eta_min, -4.0)), -0.25);

	// TODO: figure out this negative sign in front of 2 below - where did it come from??
		eta_test = pow(((-7.0*area/(2.0*norm*kk)) + pow(eta_min, -3.5)), (-2.0/7.0));
		/* now choose second test coordinate (y) by choosing a number between zero and the value of the comparison function at eta_test, comp_ymax */
		comp_ymax = kk * norm * pow(eta_test, -4.5);
		comp_y = rng_t113_dbl()*comp_ymax;
		/* now accept or reject the point if it lies below/above true distribution: 
		d(rate)/d(eta) ~ norm*(5*eta^-6 + 8*eta^-5) */
	//	true_y = norm*(5.0*pow(eta_test, -6.0) + 8.0*pow(eta_test, -5.0));
		true_y = norm*(5.5*pow(eta_test, -6.5) + (4.5 * kk +9.0) * (pow(eta_test, -5.5)) +7.0 * kk * (pow(eta_test, -4.5)));
		if (comp_y < true_y) { /* fall within true distribution */
			found_eta = 1;
			eta = eta_test;
		}
	}
	return(eta);
}


void make_threebodybinary(double P_3bb, long k1, long k2, long k3, long form_binary, double eta_min, double ave_local_mass, double n_local, double sigma_local, double v1[4], double v2[4], double v3[4], double vrel12[4], double vrel3[4], double delta_E_3bb, gsl_rng *rng)
{
	double PE_i, PE_f, KE_i, KE_f, KE_cmf_i, KE_cmf_f, delta_PE, delta_KE, delta_E, system_cm, binary_cm;
	double m1, m2, m3, ms, mb;
	double cm_vel[4], v1_cmf[4], v2_cmf[4], v3_cmf[4];
	double eta, Eb, ecc_max, ecc, r_p, semi_major;
	double vs_cmf[4], vb_cmf[4], vs[4], vb[4], angle3, angle4;
	long j, knew;
	// Form new binary, set new binary/stellar properties, destroy old stars	
	// Calculate Total potential energy of three stars (enforce cons. of potential energy at end to get new positions)
	// Find energy of new binary
	// Move to COM frame: for initial system of 3 stars, compute com veloc, linear momentum and Energy
	// Choose random direction for kick to single in COM frame - given direction of vs, know that vb points oppositely
	// Solve cons. of momentum and cons. of energy to obtain scalar velocities of binary and single, vb and vs

	
	// INITIAL VALUES
	m1=star[k1].m;
	m2=star[k2].m;
	m3=star[k3].m;

	// Initial energy
	PE_i = madhoc*(star[k1].m*star[k1].phi + star[k2].m*star[k2].phi + star[k3].m*star[k3].phi);
	KE_i = 0.5*madhoc*(star[k1].m * sqr(v1[0]) + star[k2].m * sqr(v2[0]) + star[k3].m * sqr(v3[0]));

	// COM position for three star system
	//system_cm = (star[k1].m * star[k1].r + star[k2].m * star[k2].r + star[k3].m * star[k3].r)/(star[k1].m + star[k2].m + star[k3].m);

	// COM of candidate binary pair
	binary_cm = (star[k1].m * star[k1].r + star[k2].m * star[k2].r)/(star[k1].m + star[k2].m);

	for (j=1; j<=3; j++) {
		cm_vel[j] = (star[k1].m * v1[j] + star[k2].m * v2[j] + star[k3].m * v3[j])/(star[k1].m + star[k2].m + star[k3].m);
	}

	cm_vel[0] = sqrt(sqr(cm_vel[1]) + sqr(cm_vel[2]) + sqr(cm_vel[3]));

	// v1, v2, v3 were calculated in function calc_3bb_encounter_dyns()

	// Move to COM frame
	for (j=1; j<=3; j++) {
		v1_cmf[j] = v1[j] - cm_vel[j];
		v2_cmf[j] = v2[j] - cm_vel[j];
		v3_cmf[j] = v3[j] - cm_vel[j];
	} 

	// scalar velocities, in COM frame
	v1_cmf[0] = sqrt(sqr(v1_cmf[1]) + sqr(v1_cmf[2]) + sqr(v1_cmf[3]));
	v2_cmf[0] = sqrt(sqr(v2_cmf[1]) + sqr(v2_cmf[2]) + sqr(v2_cmf[3]));
	v3_cmf[0] = sqrt(sqr(v3_cmf[1]) + sqr(v3_cmf[2]) + sqr(v3_cmf[3]));

	// Initial kinetic energy of three stars in COM frame (relative to COM motion)
	KE_cmf_i = 0.5*madhoc*(star[k1].m * sqr(v1_cmf[0]) + star[k2].m * sqr(v2_cmf[0]) + star[k3].m * sqr(v3_cmf[0]));

	// COMPUTE NEW QUANTITIES
		// Binary orbital properties: 
		//	binding energy: choose eta, then Eb = eta * <m> / sigma^2
		//	eccentricity: chosen from thermal distribution
	eta = get_eta(eta_min, k1, k2, k3, vrel12, vrel3);
	Eb = -0.5*eta * ave_local_mass * sqr(sigma_local); //Note: ave_local_mass is already multiplied by madhoc
	// Note binding energy is NEGATIVE
//	Eb = -1.0*eta * madhoc * sqr(sigma_local); // Note binding energy is NEGATIVE
	ecc_max = 1.0;
	ecc = ecc_max * sqrt(rng_t113_dbl()); 					
	r_p = star[k1].m * star[k2].m * madhoc / (eta * sqr(sigma_local)); //note, for units, have madhoc^2 in numerator, and madhoc in the denominator ====> single power madhoc in numerator
//	semi_major = r_p/(1.0-ecc);
	semi_major = star[k1].m * star[k2].m * sqr(madhoc) / (-2*Eb);
//	fprintf(threebbdebugfile, "eta=%g madhoc=%g sigma_local=%g\n", eta, madhoc, sigma_local);
		// Using cons. of momentum and energy, can calculate scalar velocities of binary and single
	vs_cmf[0] = sqrt((2.0*(KE_cmf_i-Eb)*((star[k1].m + star[k2].m)/star[k3].m))/((madhoc)*(star[k1].m + star[k2].m + star[k3].m)));
	vb_cmf[0] = 1.0*vs_cmf[0]*star[k3].m/(star[k1].m + star[k2].m);

		// Choose random direction for motion of single star
	angle3 = rng_t113_dbl() * PI; // polar angle - [0, PI)
	angle4 = rng_t113_dbl() * 2.0 * PI; // azimuthal angle - [0, 2*PI)
		// Compute vector velocities of single and binary
	vs_cmf[1] = vs_cmf[0] * sin(angle3) * cos(angle4);
	vs_cmf[2] = vs_cmf[0] * sin(angle3) * sin(angle4);
	vs_cmf[3] = vs_cmf[0] * cos(angle3);

	// Fix direction of velocity of binary to be opposite to that of single - add Pi to each of the angles (polar angle and azimuthal) that were used to orient the velocity of the single
	vb_cmf[1] = vb_cmf[0] * sin(PI - angle3) * cos(angle4 + PI);
	vb_cmf[2] = vb_cmf[0] * sin(PI - angle3) * sin(angle4 + PI);
	vb_cmf[3] = vb_cmf[0] * cos(PI - angle3);

	vs_cmf[0] = sqrt(sqr(vs_cmf[1]) + sqr(vs_cmf[2]) + sqr(vs_cmf[3]));
	vb_cmf[0] = sqrt(sqr(vb_cmf[1]) + sqr(vb_cmf[2]) + sqr(vb_cmf[3]));

	// Final COM energy - Binding energy of new binary plus kinetic energy of single and binary, calculated using the magnitudes of their velocities (vs_cmf[0] and vb_cmf[0])
	KE_cmf_f = Eb + 0.5*madhoc*(star[k3].m*sqr(vs_cmf[0]) + (star[k1].m + star[k2].m)*sqr(vb_cmf[0]));


	// DEBUG
	//fprintf(threebbdebugfile, "r1=%g  r2=%g  r3=%g\n", star[k1].r, star[k2].r, star[k3].r);
//	fprintf(threebbdebugfile, "m1=%g  m2=%g  m3=%g  madhoc=%g\n", star[k1].m, star[k2].m, star[k3].m, madhoc);
	//fprintf(threebbdebugfile, "cm_vel[1]=%g  cm_vel[2]=%g  cm_vel[3]=%g\n", cm_vel[1], cm_vel[2], cm_vel[3]);
	//fprintf(threebbdebugfile, "cm_vel=(%g, %g, %g, %g)  v1_cmf=(%g, %g, %g, %g)  v2_cmf=(%g, %g, %g, %g)  v3_cmf=(%g, %g, %g, %g)\n", cm_vel[0], cm_vel[1], cm_vel[2], cm_vel[3], v1_cmf[0], v1_cmf[1], v1_cmf[2], v1_cmf[3], v2_cmf[0], v2_cmf[1], v2_cmf[2], v2_cmf[3], v3_cmf[0], v3_cmf[1], v3_cmf[2], v3_cmf[3]);
	// convert to physical units
//	fprintf(threebbdebugfile, "Double check binding energy: Eb= %g = which should = %g = m1*madhoc*m2*madhoc / 2*a\n", Eb, star[k1].m * star[k2].m * sqr(madhoc) / (2*semi_major));
//	fprintf(threebbdebugfile, "code units:  eta=%g  sigma=%g  Eb=%g  ecc=%g  rp=%g  a=%g  KE_cmf_i=%g\n", eta, sigma_local, Eb, ecc, r_p, semi_major, KE_cmf_i);
//	fprintf(threebbdebugfile, "Physical units:  sigma=%g  Eb=%g  rp=%g  a=%g  KE_cmf_i=%g\n", sigma_local*(units.l/units.t), Eb*(units.m)*sqr(units.l / units.t), r_p * 6.674e-8 * units.m / sqr(units.l / units.t), r_p * 6.674e-8 * units.m / (sqr(units.l / units.t) * (1.0-ecc)), KE_cmf_i*(units.m)*sqr(units.l/units.t));
	//fprintf(threebbdebugfile, "CHECK CONSERVATION OF MOMENTUM:  BEFORE\n");
	//fprintf(threebbdebugfile, "p1=%g  p2=%g  p3=%g\n", (star[k1].m*madhoc*v1_cmf[1]) + (star[k2].m*madhoc*v2_cmf[1]) + (star[k3].m*madhoc*v3_cmf[1]), (star[k1].m*madhoc*v1_cmf[2]) + (star[k2].m*madhoc*v2_cmf[2]) + (star[k3].m*madhoc*v3_cmf[2]), (star[k1].m*madhoc*v1_cmf[3]) + (star[k2].m*madhoc*v2_cmf[3]) + (star[k3].m*madhoc*v3_cmf[3]));
	//fprintf(threebbdebugfile, "CODE:  v1=(%g, %g, %g, %g)  v2=(%g, %g, %g, %g)  v3=(%g, %g, %g, %g)\n", v1[0], v1[1], v1[2], v1[3], v2[0], v2[1], v2[2], v2[3], v3[0], v3[1], v3[2], v3[3]);
	//fprintf(threebbdebugfile, "CODE:  vs_cmf[0]=%g  vb_cmf[0]=%g  cm_vel[0]=%g\n", vs_cmf[0], vb_cmf[0], cm_vel[0]);
	//fprintf(threebbdebugfile, "PHYSICAL:  v1=(%g, %g, %g, %g)  v2=(%g, %g, %g, %g)  v3=(%g, %g, %g, %g)\n", v1[0]*units.l/units.t, v1[1]*units.l/units.t, v1[2]*units.l/units.t, v1[3]*units.l/units.t, v2[0]*units.l/units.t, v2[1]*units.l/units.t, v2[2]*units.l/units.t, v2[3]*units.l/units.t, v3[0]*units.l/units.t, v3[1]*units.l/units.t, v3[2]*units.l/units.t, v3[3]*units.l/units.t);
	//fprintf(threebbdebugfile, "PHYSICAL:  vs_cmf[0]=%g  vb_cmf[0]=%g  cm_vel[0]=%g\n", vs_cmf[0]*units.l/units.t, vb_cmf[0]*units.l/units.t, cm_vel[0]*units.l/units.t);
	//fprintf(threebbdebugfile, "sin(a3)*cos(a4)=%g  sin(a3)*sin(a4)=%g  cos(a3)=%g\n",sin(angle3) * cos(angle4), sin(angle3) * sin(angle4), cos(angle3));
	//fprintf(threebbdebugfile, "sin(a3 + pi)*cos(a4 + pi)=%g  sin(a3 + pi)*sin(a4 + pi)=%g  cos(a3 + pi)=%g\n", sin(PI - angle3) * cos(angle4 - PI), sin(PI - angle3) * sin(angle4 + PI), cos(PI -angle3));
	//fprintf(threebbdebugfile, "CODE:  vs_cmf[0]=%g  vb_cmf[0]=%g  cm_vel[0]=%g\n", vs_cmf[0], vb_cmf[0], cm_vel[0]);
	//fprintf(threebbdebugfile, "CHECK CONSERVATION OF MOMENTUM:  AFTER\n");
	//fprintf(threebbdebugfile, "(star[k1].m + star[k2].m)=%g\n  vb_cmf[1]=%g  vb_cmf[2]=%g  vb_cmf[3]=%g  vs_cmf[1]=%g  vs_cmf[2]=%g  vs_cmf[3]=%g\n", star[k1].m + star[k2].m,  vb_cmf[1],   vb_cmf[2], vb_cmf[3], vs_cmf[1], vs_cmf[2], vs_cmf[3]);
	//fprintf(threebbdebugfile, "p1=%g  p2=%g  p3=%g\n", ((star[k1].m+star[k2].m)*madhoc*vb_cmf[1]) + (star[k3].m*madhoc*vs_cmf[1]), ((star[k1].m+star[k2].m)*madhoc*vb_cmf[2]) + (star[k3].m*madhoc*vs_cmf[2]), ((star[k1].m+star[k2].m)*madhoc*vb_cmf[3]) + (star[k3].m*madhoc*vs_cmf[3]));
	//fprintf(threebbdebugfile, "CODE:  KE_cmf_i=%g  KE_cmf_f+Eb=%g  Delta_E_cmf=%g\n", KE_cmf_i, E_b + 0.5*madhoc*( (star[k1].m + star[k2].m)*sqr(vb_cmf[0]) + star[k3].m*sqr(vs_cmf[0]) ),  0.5*madhoc*( (star[k1].m + star[k2].m)*sqr(vb_cmf[0]) + star[k3].m*sqr(vs_cmf[0]) )-KE_cmf_i );
	//fprintf(threebbdebugfile, "CODE:  KE_cmf_i=%g  KE_cmf_f=%g  Delta_E_cmf=%g\n", KE_cmf_i, KE_cmf_f, KE_cmf_f - KE_cmf_i );
	//fprintf(threebbdebugfile, "PHYSICAL:  KE_cmf_i=%g  KE_cmf_f+Eb=%g  Delta_E_cmf=%g\n", KE_cmf_i*(units.m)*sqr(units.l/units.t), (E_b + 0.5*madhoc*( (star[k1].m + star[k2].m)*sqr(vb_cmf[0]) + star[k3].m*sqr(vs_cmf[0])))*units.m*sqr(units.l/units.t),  (0.5*madhoc*( (star[k1].m + star[k2].m)*sqr(vb_cmf[0]) + star[k3].m*sqr(vs_cmf[0]) )-KE_cmf_i)*units.m*sqr(units.l/units.t) );
	//fprintf(threebbfile, "vs_cmf[1]=%g  vs_cmf[2]=%g  vs_cmf[3]=%g\n", vs_cmf[1], vs_cmf[2], vs_cmf[3]);
//	fprintf(threebbfile, "vb_cmf[1]=%g  vb_cmf[2]=%g  vb_cmf[3]=%g\n", vb_cmf[1], vb_cmf[2], vb_cmf[3]);
//	fprintf(threebbfile, "cm_vel[1]=%g  cm_vel[2]=%g  cm_vel[3]=%g\n", cm_vel[1], cm_vel[2], cm_vel[3]);


		// Transform back to lab frame: final velocities of single and binary are vs and vb
	for (j=1; j<=3; j++) {
		vs[j] = vs_cmf[j] + cm_vel[j];
		vb[j] = vb_cmf[j] + cm_vel[j];	
	}
	vs[0] = sqrt(sqr(vs[1]) + sqr(vs[2]) + sqr(vs[3]));
	vb[0] = sqrt(sqr(vb[1]) + sqr(vb[2]) + sqr(vb[3]));
	// Only if binary will be formed, set star and binary properties for new binary, and change properties of single star
	if (form_binary == 0) { // leave three stars alone; allow them to interact in two-body loop
		star[k1].threebb_interacted = 0;
		star[k2].threebb_interacted = 0;
		star[k3].threebb_interacted = 0;

		fprintf(lightcollisionfile, "%.16g %ld %ld %ld %ld %ld %ld %g %g %g %d %d %d %g %g %g %g %g %g %g\n", TotalTime, k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, star[k1].m*(units.m / clus.N_STAR / MSUN), star[k2].m * (units.m / clus.N_STAR / MSUN), star[k3].m *(units.m / clus.N_STAR / MSUN), star[k1].se_k, star[k2].se_k, star[k3].se_k, star[k1].rad * units.l / AU, star[k2].rad * units.l / AU, star[k3].rad * units.l / AU, Eb, ecc, semi_major * units.l / AU, r_p * units.l / AU);
	}
	else if (form_binary == 1) {

		// DEBUG extra output
		//fprintf(threebbdebugfile, "angle3=%g  angle4=%g\n", angle3, angle4);
		//fprintf(threebbdebugfile, "vs_cmf=(%g, %g, %g, %g)  vb_cmf=(%g, %g, %g, %g)\n", vs_cmf[0], vs_cmf[1], vs_cmf[2], vs_cmf[3], vb_cmf[0], vb_cmf[1], vb_cmf[2], vb_cmf[3]);
		//fprintf(threebbdebugfile, "vs=(%g, %g, %g, %g)  vb=(%g, %g, %g, %g)\n", vs[0], vs[1], vs[2], vs[3], vb[0], vb[1], vb[2], vb[3]);
		//fprintf(threebbdebugfile, "PE_o=%g  E_o-PE_o=%g  KE_cmf_i=%g  eta=%g  Eb=%g\n", PE_i, 0.5*madhoc*((star[k1].m*(sqr(star[k1].vr)+sqr(star[k1].vt))) + (star[k2].m*(sqr(star[k2].vr)+sqr(star[k2].vt))) + (star[k3].m*(sqr(star[k3].vr)+sqr(star[k3].vt)))), KE_cmf_i, eta, Eb);
		//fprintf(threebbdebugfile, "Eint-PE=%g\n", 0.5*madhoc*((star[k1].m + star[k2].m)*sqr(vb[0]) + star[k3].m*sqr(vs[0])));

	// Only if binary will be formed, set star and binary properties for new binary, and change properties of single star
		knew = create_binary(); /* returns an unused binary id */

		star[knew].threebb_interacted = 1;
		star[k1].threebb_interacted = 1;
		star[k2].threebb_interacted = 1;
		star[k3].threebb_interacted = 1;


		star[knew].interacted = 1;
		star[k1].interacted = 1;
		star[k2].interacted = 1;
		star[k3].interacted = 1;

		binary[star[knew].binind].id1 = star[k1].id;
		binary[star[knew].binind].id2 = star[k2].id;
		binary[star[knew].binind].rad1 = star[k1].rad;
		binary[star[knew].binind].rad2 = star[k2].rad;
		binary[star[knew].binind].m1 = star[k1].m;
		binary[star[knew].binind].m2 = star[k2].m;
		binary[star[knew].binind].Eint1 = star[k1].Eint;
		binary[star[knew].binind].Eint2 = star[k2].Eint;
		binary[star[knew].binind].a = semi_major;
		binary[star[knew].binind].e = ecc;
		//Sourav: toy rejuvenation- some new variables
		//binary[knew].createtime_m1 = star[k1].createtime; // These are all set in cp_SEvars_to_newbinary()
		//binary[knew].createtime_m2 = star[k2].createtime;
		//binary[knew].lifetime_m1 = star[k1].lifetime;;
		//binary[knew].lifetime_m2 = star[k2].lifetime;
			// ????? For star properties that do not apply for binary (such as rad...since we don't 
			// have a single star with a single radius), what do I set to? These are quantities
			// that are printed as "-100" in the output files. 


		// star properties for the new binary
//		star[knew].binind = 1;  OOPS! I thought that binind=1 for binary. Actually, this is the binary index (when you run "knew = create_binary()", it finds the first unused binary index, and assigns star[knew].binind the value of this available binary index. So for a single star k, star[k].binind = 0, and for a binary, star[k].binind >0.
		star[knew].r = binary_cm;
		star[knew].vr = vb[3];
		star[knew].vt = sqrt(sqr(vb[1]) + sqr(vb[2]));
		star[knew].m = star[k1].m + star[k2].m;
		star[knew].Eint = star[k1].Eint + star[k2].Eint;  // Set new internal energy for the binary, sum of Eint1 and Eint2 of the two components


	//	star[knew].E = 0.0;
	//	star[knew].J = 0.0;
	//	star[knew].EI = 0.0; 		// Intermediate energy: seems that I can leave this set to zero?
	//	star[knew].Eint = 0.0;
	//	star[knew].rnew = 0.0;		// set in set_star_news() function
	//	star[knew].vrnew = 0.0;
	//	star[knew].vtnew = 0.0;
	//	star[knew].rOld = 0.0;		// do I have to set this? What is it?
	//	star[knew].X = 0.0;
	//	star[knew].Y = 0.0;
	//	star[knew].r_peri = 0.0; 	// these are set later
	//	star[knew].r_apo = 0.0;
	//	star[knew].phi = 0.0;
	//	star[knew].interacted = 0;
	//	star[knew].id = 0;		// what do I set ID to in binary
	//	star[knew].rad = 0.0;		// what about Uoldrold and Uoldrnew?
	//	star[knew].Uoldrold = 0.0;
	//	star[knew].Uoldrnew = 0.0;

	//Sourav: toy rejuvenation- some more variables
	//	star[knew].vtold = 0.0;		// What about these?
	//	star[knew].vrold = 0.0;

	//	fprintf(threebbfile, "new binary properties for object 'knew':  r=%g  vr=%g  vt=%g\n", star[knew].r, star[knew].vr, star[knew].vt);

		// Set new properties of single star
	//		star[k3].interacted = 1;
		star[k3].r = star[k3].r;  // later I can adjust the positions/potentials
		star[k3].vr = vs[3];
		star[k3].vt = sqrt(sqr(vs[1]) + sqr(vs[2]));

		// copy SE variables over to new binary from old single stars
		cp_SEvars_to_newbinary(k1, -1, knew, 0);
		cp_SEvars_to_newbinary(k2, -1, knew, 1);

		// DEBUG
		//fprintf(threebbdebugfile, "k1=%ld  k2=%ld  knew=%ld  star[k1].binind=%ld  star[k2].binind=%ld  star[knew].binind=%ld\n", k1, k2, knew, star[k1].binind, star[k2].binind, star[knew].binind);

		// radii, and tb : rhs was set in calls to cp_SEvars_to_newbinary()
		binary[star[knew].binind].rad1 = binary[star[knew].binind].bse_radius[0] * RSUN / units.l;
		binary[star[knew].binind].rad2 = binary[star[knew].binind].bse_radius[1] * RSUN / units.l;
		binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
		//binary[star[knew].binind].threebodybinary = 1;
		// DEBUG
		//fprintf(threebbdebugfile, "binary[  ].rad1=%g  binary[  ].rad2=%g  binary[   ].bse_tb=%g\n", binary[star[knew].binind].rad1, binary[star[knew].binind].rad2, binary[star[knew].binind].bse_tb);

		/* track binding energy */
		//BEf += binary[star[knew].binind].m1 * binary[star[knew].binind].m2 * sqr(madhoc)
		// Do I need to worry about this binding energy??

		star[knew].phi = potential(star[knew].r); // set
		star[k3].phi = potential(star[k3].r); // set potential for single to be same as original...should I adjust?
		PE_f = madhoc*(star[knew].m*star[knew].phi + star[k3].m*star[k3].phi);
		KE_f = Eb + 0.5*madhoc*(star[knew].m * sqr(vb[0]) + star[k3].m * sqr(vs[0]));


		// In the right units for comparing to the Energies computed in cmc_utils.c ComputeEnergy(), and plotted in quick_cluster_plot
		delta_PE = (PE_f - PE_i);
		delta_KE = (KE_f - KE_i);
		delta_E = delta_KE + delta_PE;

		fprintf(threebbdebugfile, "%ld %ld %ld %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, star[k1].r, star[k2].r, star[k3].r, star[k1].m, star[k2].m, star[k3].m, v1[0], v1[1], v1[2], v1[3], v2[0], v2[1], v2[2], v2[3], v3[0], v3[1], v3[2], v3[3], v1_cmf[0], v1_cmf[1], v1_cmf[2], v1_cmf[3], v2_cmf[0], v2_cmf[1], v2_cmf[2], v2_cmf[3], v3_cmf[0], v3_cmf[1], v3_cmf[2], v3_cmf[3]);

		fprintf(threebbdebugfile, "%ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", knew, star[knew].id, star[knew].r, star[knew].m, vs_cmf[0], vs_cmf[1], vs_cmf[2], vs_cmf[3], vb_cmf[0], vb_cmf[1], vb_cmf[2], vb_cmf[3], vs[0], vs[1], vs[2], vs[3], vb[0], vb[1], vb[2], vb[3], ave_local_mass, n_local,sigma_local, eta, Eb, ecc, r_p, semi_major, PE_i, PE_f, KE_cmf_i, KE_cmf_f, KE_i, KE_f, delta_PE, delta_KE, delta_E);
		//fprintf(threebbdebugfile, "Internal Energies:  Eint1=%g  Eint2=%g  Eint_binary=%g\n", star[k1].Eint, star[k2].Eint, star[knew].Eint);



		// running total for simulation - change in energy occuring during 3bb formation procedure
		delta_E_3bb += delta_E;

		set_star_EJ(knew);	// binary
		set_star_EJ(k3);	// single

		binary[star[knew].binind].a = semi_major;
		binary[star[knew].binind].e = ecc;

		// destroy the two former single stars (which have now formed a binary)
		// leave the remaining single star (properties have already been updated)
		fprintf(threebbfile, "%.16g %ld %ld %ld %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %ld\n", TotalTime, k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, m1*(units.m / clus.N_STAR / MSUN), m2*(units.m / clus.N_STAR / MSUN), m3*(units.m / clus.N_STAR / MSUN), ave_local_mass, n_local, sigma_local, eta, Eb, binary[star[knew].binind].e, binary[star[knew].binind].a * units.l / AU, r_p * units.l / AU, star[knew].r, star[k3].r, star[knew].vr, star[knew].vt, star[k3].vr, star[k3].vt, star[knew].phi, star[k3].phi, delta_PE, delta_KE, delta_E, delta_E_3bb, N3bbformed);

		dprintf("Energies after binary is formed");
		ComputeEnergy();
		destroy_obj(k1);
		destroy_obj(k2);
	//	fprintf(threebbfile, "KE_cmf_i=%g  eta=%g  Eb=%g\n", KE_cmf_i, eta, Eb);
	//	fprintf(threebbfile, "vs_cmf[0]=%g  vb_cmf[0]=%g  cm_vel[0]=%g\n", vs_cmf[0], vb_cmf[0], cm_vel[0]);
	//	fprintf(threebbfile, "new binary properties for object 'knew':  r=%g  vr=%g  vt=%g\n", star[knew].r, star[knew].vr, star[knew].vt);

		// DEBUG extra output
		//fprintf(threebbdebugfile, "Phi_f=%g  E_f-phi_f=%g\n", star[knew].phi + star[k3].phi, Eb + 0.5*madhoc*(star[k3].m*(sqr(star[k3].vr) + sqr(star[k3].vt)) + star[knew].m*(sqr(star[knew].vr) + sqr(star[knew].vt))));

		//fprintf(threebbdebugfile, "ms=%g  mb=%g  vs_r=%g  vs_t=%g  vb_r=%g  vb_t=%g\n\n\n", star[k3].m, star[knew].m, star[k3].vr, star[k3].vt, star[knew].vr, star[knew].vt);
		//fprintf(threebbdebugfile, "binary[star[knew].binind].inuse=%ld\n", binary[star[knew].binind].inuse);
		//fprintf(threebbfile, "%g %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", TotalTime, k1, k2, k3, m1, m2, m3, eta, Eb, binary[star[knew].binind].e, binary[star[knew].binind].a * units.l / AU, r_p * units.l / AU, star[knew].r, star[knew].vr, star[knew].vt, star[knew].phi, star[k3].r, star[k3].vr, star[k3].vt, star[k3].phi, delta_PE, delta_E, delta_E_3bb);


		}

}


void calc_sigma_local(long k, long p, long N_LIMIT, double *ave_local_mass, double *sigma_local)
{
        long si, simin, simax;
        double Mv2ave, Mave, M2ave, sigma;

//        N_LIMIT = clus.N_MAX;
	simin = k - p;
        simax = simin + (2 * p + 1);
        if (simin < 1) {
		simin = 1;
                simax = simin + (2 * p + 1);
        } else if (simax > N_LIMIT) {
                simax = N_LIMIT;
                simin = simax - (2 * p + 1);
        }

        Mv2ave = 0.0;
        Mave = 0.0;
        M2ave = 0.0;
        for (si=simin; si<=simax; si++) {
		Mv2ave += star[si].m * madhoc * (sqr(star[si].vr) + sqr(star[si].vt));
                Mave += star[si].m * madhoc;
                M2ave += sqr(star[si].m * madhoc);
	}
        Mv2ave /= (double) (2 * p);
        Mave /= (double) (2 * p);
        M2ave /= (double) (2 * p);

        /* sigma is the 3D velocity dispersion */
        sigma = sqrt(Mv2ave/Mave);

	*sigma_local = sigma;
	*ave_local_mass = Mave;
        /* average relative speed for a Maxwellian, from Binney & Tremaine */
  //      vrel_ave = 4.0 * sigma / sqrt(3.0 * PI);
//	return(sigma);
}


// end functions used in three-body binary formation


/* generate unique star id's */
long star_get_id_new(void)
{
	newstarid++;
	return(newstarid);
}

/* calculate local density by averaging */
double calc_n_local(long k, long p, long N_LIMIT)
{
	long kmin, kmax;

	kmin = k - p;
	kmax = k + p + 1;

	/* shouldn't the boundary here be kmin=1, not 0? */
	if (kmin < 0) {
		kmin = 0;
		kmax = 2 * p + 1;
	} else if (kmax > N_LIMIT) {
		kmax = N_LIMIT;
		kmin = N_LIMIT - 2 * p - 1;
	}
	
	return((2.0 * ((double) p)) * 3.0 / (4.0 * PI * (cub(star[kmax].r) - cub(star[kmin].r))));
}

double calc_Ai_local(long k, long kp, long p, double W, long N_LIMIT)
{
	long kmin, kmax;

	kmin = k - p;
	kmax = k + p + 1;

	if (kmin < 0) {
		kmin = 0;
		kmax = 2 * p + 1;
	} else if (kmax > N_LIMIT) {
		kmax = N_LIMIT;
		kmin = N_LIMIT - 2 * p - 1;
	}

	return(3.0 * ((double) p) * sqr(star[k].m + star[kp].m) / (cub(W) * (cub(star[kmax].r) - cub(star[kmin].r))));
}

void calc_encounter_dyns(long k, long kp, double v[4], double vp[4], double w[4], double *W, double *rcm, double vcm[4], gsl_rng *rng, int setY)
{
	int j;
	double phi;

	/* set random angle between vt's */
	/* first store random variable */
	if (setY) {
		star[k].Y = rng_t113_dbl();
		star[kp].Y = star[k].Y;
	}
	phi = star[k].Y * 2.0 * PI;
	v[1] = star[k].vt;
	v[2] = 0.0;
	v[3] = star[k].vr;
	vp[1] = star[kp].vt * cos(phi);
	vp[2] = star[kp].vt * sin(phi);
	vp[3] = star[kp].vr;
		
	for (j=1; j<=3; j++) {
		w[j] = vp[j] - v[j];
	}
	
	*W = sqrt(sqr(w[1]) + sqr(w[2]) + sqr(w[3]));

	if (*W == 0.0) {
		eprintf("W = 0!\n");
		exit_cleanly(1);
	}
		
	/* compute CM quantities */
	*rcm = (star[k].m * star[k].r + star[kp].m * star[kp].r) / (star[k].m + star[kp].m);
	for (j=1; j<=3; j++) {
		vcm[j] = (star[k].m * v[j] + star[kp].m * vp[j]) / (star[k].m + star[kp].m);
	}
}

void set_star_EJ(long k)
{
	star[k].E = star[k].phi + 0.5 * (sqr(star[k].vr) + sqr(star[k].vt));
	star[k].J = star[k].r * star[k].vt;
}

void set_star_news(long k)
{
	star[k].rnew = star[k].r;
	star[k].vrnew = star[k].vr;
	star[k].vtnew = star[k].vt;
}

void set_star_olds(long k)
{
	star[k].rOld = star[k].r;
	star[k].r_peri = star[k].r;
	star[k].r_apo = star[k].r;
}

/* find masses of merging stars from binary interaction components */
double binint_get_mass(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].m);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].m1);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].m2);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].m);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].m1);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].m2);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1);
	/* this is just for the compiler */
	exit(1);
}

/*Sourav: find stellar types of merging stars from binary interaction components */
long binint_get_startype(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_k);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_kw[0]);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_kw[1]);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_k);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_kw[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_kw[1]);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1);
	/* this is just for the compiler */
	exit(1);
}


//Sourav: toy rejuvenation- finding creation times of binary interaction components
double binint_get_createtime(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].createtime);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].createtime_m1);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].createtime_m2);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].createtime);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].createtime_m1);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].createtime_m2);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1);
	/* this is just for the compiler */
	exit(1);
}

//Sourav: toy rejuvenation- finding MS lifetimes of binary interaction components
double binint_get_lifetime(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].lifetime);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].lifetime_m1);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].lifetime_m2);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].lifetime);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].lifetime_m1);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].lifetime_m2);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1);
	/* this is just for the compiler */
	exit(1);
}


/* return star index of star with id "id" from binary interaction components,
   along with which member "bi=0,1" if a binary */
long binint_get_indices(long k, long kp, long id, int *bi)
{
	*bi = -1;

	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(k);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			*bi = 0;
			return(k);
		} else if (binary[star[k].binind].id2 == id) {
			*bi = 1;
			return(k);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(kp);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			*bi = 0;
			return(kp);
		} else if (binary[star[kp].binind].id2 == id) {
			*bi = 1;
			return(kp);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1);
	/* this is just for the compiler */
	exit(1);
}

void binint_log_obj(fb_obj_t *obj, fb_units_t units)
{
	int bid, sid, i;
	char dumstring[FB_MAX_STRING_LENGTH], idstring1[FB_MAX_STRING_LENGTH], idstring2[FB_MAX_STRING_LENGTH], idstring3[FB_MAX_STRING_LENGTH];
	
	if (fb_n_hier(obj) == 1) {
		/* first write id string */
		snprintf(idstring1, FB_MAX_STRING_LENGTH, "%ld", obj->id[0]);
		for (i=1; i<obj->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->id[i]);
			strncat(idstring1, dumstring, FB_MAX_STRING_LENGTH);
		}
		/* then print to log */
		fprintf(binintfile, "type=single m=%g R=%g Eint=%g id=%s\n", obj->m*units.m/MSUN, obj->R*units.l/RSUN, obj->Eint*units.E, idstring1);
	} else if (fb_n_hier(obj) == 2) {
		/* first write id strings */
		snprintf(idstring1, FB_MAX_STRING_LENGTH, "%ld", obj->obj[0]->id[0]);
		for (i=1; i<obj->obj[0]->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->obj[0]->id[i]);
			strncat(idstring1, dumstring, FB_MAX_STRING_LENGTH);
		}
		snprintf(idstring2, FB_MAX_STRING_LENGTH, "%ld", obj->obj[1]->id[0]);
		for (i=1; i<obj->obj[1]->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->obj[1]->id[i]);
			strncat(idstring2, dumstring, FB_MAX_STRING_LENGTH);
		}
		/* then print to log */
		fprintf(binintfile, "type=binary m0=%g m1=%g R0=%g R1=%g Eint1=%g Eint2=%g id0=%s id1=%s a=%g e=%g\n", 
			obj->obj[0]->m*units.m/MSUN, obj->obj[1]->m*units.m/MSUN, 
			obj->obj[0]->R*units.l/RSUN, obj->obj[1]->R*units.l/RSUN, 
			obj->obj[0]->Eint*units.E, obj->obj[1]->Eint*units.E, 
			idstring1, idstring2, 
			obj->a*units.l/AU, obj->e);
	} else if (fb_n_hier(obj) == 3) {
		/* identify inner binary */
		if (obj->obj[0]->n==2) {
			bid = 0;
			sid = 1;
		} else {
			bid = 1;
			sid = 0;
		}
		/* first write id strings */
		snprintf(idstring1, FB_MAX_STRING_LENGTH, "%ld", obj->obj[bid]->obj[0]->id[0]);
		for (i=1; i<obj->obj[bid]->obj[0]->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->obj[bid]->obj[0]->id[i]);
			strncat(idstring1, dumstring, FB_MAX_STRING_LENGTH);
		}
		snprintf(idstring2, FB_MAX_STRING_LENGTH, "%ld", obj->obj[bid]->obj[1]->id[0]);
		for (i=1; i<obj->obj[bid]->obj[1]->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->obj[bid]->obj[1]->id[i]);
			strncat(idstring2, dumstring, FB_MAX_STRING_LENGTH);
		}
		snprintf(idstring3, FB_MAX_STRING_LENGTH, "%ld", obj->obj[sid]->id[0]);
		for (i=1; i<obj->obj[sid]->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->obj[sid]->id[i]);
			strncat(idstring3, dumstring, FB_MAX_STRING_LENGTH);
		}
		/* then print to log */
		fprintf(binintfile, "type=triple min0=%g min1=%g mout=%g Rin0=%g Rin1=%g Rout=%g Eintin0=%g Eintin1=%g Eintout=%g idin1=%s idin2=%s idout=%s ain=%g aout=%g ein=%g eout=%g\n",
			obj->obj[bid]->obj[0]->m*units.m/MSUN, obj->obj[bid]->obj[1]->m*units.m/MSUN, obj->obj[sid]->m*units.m/MSUN,
			obj->obj[bid]->obj[0]->R*units.l/RSUN, obj->obj[bid]->obj[1]->R*units.l/RSUN, obj->obj[sid]->R*units.l/RSUN,
			obj->obj[bid]->obj[0]->Eint*units.E, obj->obj[bid]->obj[1]->Eint*units.E, obj->obj[sid]->Eint*units.E,
			idstring1, idstring2, idstring3, 
			obj->obj[bid]->a*units.l/AU, obj->a*units.l/AU,
			obj->obj[bid]->e, obj->e);
	} else {
		/* thankfully won't need to print out quads */
		eprintf("Don't know how to print out object with >3 stars!\n");
		exit_cleanly(1);
	}
}

void binint_log_status(fb_ret_t retval)
{
	/* must print out Nosc when upgraded to latest Fewbody */
	fprintf(binintfile, "status: DE/E=%g DE=%g DL/L=%g DL=%g tcpu=%g\n", 
		retval.DeltaEfrac, retval.DeltaE, retval.DeltaLfrac, retval.DeltaL, retval.tcpu);
}

void binint_log_collision(const char interaction_type[], long id,
			  double mass, double r, fb_obj_t obj, long k, long kp, long startype)
{
	int j;
	
	fprintf(collisionfile, "t=%g %s idm=%ld(mm=%g) id1=%ld(m1=%g)",
		TotalTime, interaction_type, id, 
		mass * units.mstar / FB_CONST_MSUN, obj.id[0], 
		binint_get_mass(k, kp, obj.id[0]) * units.mstar 
						  / FB_CONST_MSUN);
	for (j=1; j<obj.ncoll; j++) {
		fprintf(collisionfile, ":id%d=%ld(m%d=%g)", 
			j+1, obj.id[j], j+1,
			binint_get_mass(k, kp, obj.id[j]) * units.mstar 
							  / FB_CONST_MSUN);
	}
	fprintf(collisionfile," (r=%g) ", r);
//Sourav
	fprintf(collisionfile, "typem=%ld ", startype);
	for (j=0; j<obj.ncoll; j++) {
		fprintf(collisionfile, "type%d=%ld ", j+1, 
				binint_get_startype(k, kp, obj.id[j]));
	}
	fprintf(collisionfile, "\n");
}

/* do binary interaction (bin-bin or bin-single) */
void binint_do(long k, long kp, double rperi, double w[4], double W, double rcm, double vcm[4], gsl_rng *rng)
{
	int i, j, isbinsingle=0, isbinbin=0, sid=-1, bid=-1, istriple, bi, nmerged;
	long ksin=-1, kbin=-1, jbin, jbinp, knew, knewp=-1, oldk;
	double t, bmax, wp, wx[4], wy[4], wz[4], vnew[4], alpha, BEi, BEf=0.0;
	fb_hier_t hier;
	fb_units_t cmc_units, printing_units;
	fb_ret_t retval;
	fb_obj_t threeobjs[3];
	char string1[1024], string2[1024];
	star_t tempstar, tempstar2;
	double vs[12], VK0;

	/* perform actions that are specific to the type of binary interaction */
	if (star[k].binind != 0 && star[kp].binind != 0) {
		/* binary-binary */
		isbinbin = 1;
		N_bb++;
		hier.nstarinit = 4;

		jbin = star[k].binind;
		jbinp = star[kp].binind;

		BEi = binary[jbin].m1 * binary[jbin].m2 * sqr(madhoc) / (2.0 * binary[jbin].a)
			+ binary[jbinp].m1 * binary[jbinp].m2 * sqr(madhoc) / (2.0 * binary[jbinp].a)
			- binary[jbin].Eint1 - binary[jbin].Eint2
			- binary[jbinp].Eint1 - binary[jbinp].Eint2;
		
		cmc_units.v = sqrt((star[k].m+star[kp].m)/(star[k].m*star[kp].m) * 
				   (binary[jbin].m1*binary[jbin].m2/binary[jbin].a + 
				    binary[jbinp].m1*binary[jbinp].m2/binary[jbinp].a) * madhoc);
		cmc_units.l = binary[jbin].a + binary[jbinp].a;
		cmc_units.t = cmc_units.l / cmc_units.v;
		cmc_units.m = cmc_units.l * sqr(cmc_units.v);
		cmc_units.E = cmc_units.m * sqr(cmc_units.v);
	} else if ((star[k].binind == 0 && star[kp].binind != 0) || 
		   (star[k].binind != 0 && star[kp].binind == 0)) {
		/* binary-single */
		isbinsingle = 1;
		N_bs++;
		hier.nstarinit = 3;

		if (star[k].binind == 0) {
			ksin = k;
			kbin = kp;
			jbin = star[kp].binind;
		} else {
			ksin = kp;
			kbin = k;
			jbin = star[k].binind;
		}

		BEi = binary[jbin].m1 * binary[jbin].m2 * sqr(madhoc) / (2.0 * binary[jbin].a)
			- binary[jbin].Eint1 - binary[jbin].Eint2 - star[ksin].Eint;

		cmc_units.v = sqrt((star[ksin].m+star[kbin].m)/(star[ksin].m*star[kbin].m) * 
				   (binary[jbin].m1 * binary[jbin].m2 / binary[jbin].a) * madhoc);
		cmc_units.l = binary[jbin].a;
		cmc_units.t = cmc_units.l / cmc_units.v;
		cmc_units.m = cmc_units.l * sqr(cmc_units.v);
		cmc_units.E = cmc_units.m * sqr(cmc_units.v);
	} else {
		eprintf("no binaries!");
		exit_cleanly(1);
		exit(1);
	}
	
	/* malloc hier (based on value of hier.nstarinit) */
	fb_malloc_hier(&hier);

	bmax = rperi * sqrt(1.0 + 2.0 * ((star[k].m + star[kp].m) * madhoc) / (rperi * sqr(W)));
	
	/* call fewbody! */
	if (isbinbin) {
		retval = binbin(&t, k, kp, W, bmax, &hier, rng);
	} else {
		retval = binsingle(&t, ksin, kbin, W, bmax, &hier, rng);
	}
	
	/* set up axes */
	wp = sqrt(sqr(w[1]) + sqr(w[2]));
	if (wp == 0.0) {
		eprintf("wp = 0!\n");
		exit_cleanly(1);
	}

	/* wx, wy, and wz are the x, y, and z axes (unit vectors) of the fewbody
	   coordinate system represented in the cluster frame */
	wx[0] = 1.0;
	wx[1] = w[1]/W;
	wx[2] = w[2]/W;
	wx[3] = w[3]/W;

	wy[0] = 1.0;
	wy[1] = -w[2] / wp;
	wy[2] = w[1] / wp;
	wy[3] = 0.0;
	
	wz[0] = 1.0;
	wz[1] = -w[1] * w[3] / (wp * W);
	wz[2] = -w[2] * w[3] / (wp * W);
	wz[3] = wp / W;

	/* alpha is the factor by which to scale the velocities, as an artificial
	   way of bringing the objects to infinity while conserving energy and angular momentum */
	if (hier.nobj == 1) {
		alpha = 1.0;
	} else {
		alpha = sqrt(1.0 + fb_outerpetot(hier.obj, hier.nobj)/fb_outerketot(hier.obj, hier.nobj));
	}

	/* logging */
	binint_log_status(retval);
	printing_units.v = cmc_units.v * units.l / units.t;
	printing_units.l = cmc_units.l * units.l;
	printing_units.t = cmc_units.t * units.t;
	printing_units.m = cmc_units.m * units.m;
	printing_units.E = cmc_units.E * units.E;

	/* now do something with the Fewbody result */
	if ( !( (fabs(retval.DeltaEfrac) < 1.0e-3 || fabs(retval.DeltaE) < 1.0e-3) && 
		 (fabs(retval.DeltaLfrac) < 1.0e-3 || fabs(retval.DeltaL) < 1.0e-3) ) ) {
		/* energy error; ignore for now */
		fprintf(binintfile, "outcome: energy and/or angular momentum error\n");
		print_interaction_error();
	} else if (retval.retval == 0) {
		/* bad outcome; ignore for now */
		fprintf(binintfile, "outcome: stopped\n");
		print_interaction_error();
	} else if (hier.obj[0]->n == 4) {
		/* outcome is a quadruple */
		fprintf(binintfile, "outcome: error\n");
		print_interaction_error();
	} else {
		fprintf(binintfile, "outcome: %s (%s)\n", fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2));
		
		for (i=0; i<hier.nobj; i++) {
			/* logging */
			fprintf(binintfile, "output: ");
			binint_log_obj(hier.obj[i], printing_units);

			/* single/binary/triple stars */
			istriple = 0;
			if (hier.obj[i]->n == 1) {
				knew = create_star();
			} else if (hier.obj[i]->n == 2) {
				knew = create_binary();
			} else if (hier.obj[i]->n == 3) {
				istriple = 1;
				/* break triple for now */
				knew = create_binary();
				knewp = create_star();
				
				if (hier.obj[i]->obj[0]->n == 1) {
					sid = 0;
					bid = 1;
				} else {
					sid = 1;
					bid = 0;
				}
			} else {
				eprintf("object with n=%d!\n", hier.obj[i]->n);
				exit_cleanly(1);
				/* this is just for the compiler */
				exit(1);
			}
			
			/* generic properties */
			/* set radial position */
			star[knew].r = rcm;
			if (istriple) {
				star[knewp].r = rcm;
			}

			/* figure out new velocities */
			for (j=1; j<=3; j++) {
				vnew[j] = vcm[j] + cmc_units.v * alpha * 
					(hier.obj[i]->v[0] * wx[j] + hier.obj[i]->v[1] * wy[j]+ hier.obj[i]->v[2] * wz[j]);
			}
			
			/* set new velocities */
			star[knew].vr = vnew[3];
			star[knew].vt = sqrt(sqr(vnew[1]) + sqr(vnew[2]));
			if (istriple) {
				star[knewp].vr = vnew[3];
				star[knewp].vt = sqrt(sqr(vnew[1]) + sqr(vnew[2]));
			}
			
			/* set mass; this gets overwritten later for collisions */
			if (istriple) {
				star[knew].m = hier.obj[i]->obj[bid]->m * cmc_units.m / madhoc;
				star[knewp].m = hier.obj[i]->obj[sid]->m * cmc_units.m / madhoc;
			} else {
				star[knew].m = hier.obj[i]->m * cmc_units.m / madhoc;
			}

			/* set potential */
			star[knew].phi = potential(star[knew].r);
			if (istriple) {
				star[knewp].phi = potential(star[knewp].r);
			}
			
			/* Calculate new energies by recomputing E = PE + KE using new velocity */
			set_star_EJ(knew);
			if (istriple) {
				set_star_EJ(knewp);
			}
			
			/* set rnew, vrnew, vtnew */
			set_star_news(knew);
			if (istriple) {
				set_star_news(knewp);
			}
			
			/* I don't even know if this is necessary */
			set_star_olds(knew);
			if (istriple) {
				set_star_olds(knewp);
			}
			
			/* mark stars as interacted */
			star[knew].interacted = 1;
			if (istriple) {
				star[knewp].interacted = 1;
			}
			
			/* properties specific to single/binary/triple stars */
			if (hier.obj[i]->n == 1) {
				/* single star */
				/* internal energy */
				star[knew].Eint = hier.obj[i]->Eint * cmc_units.E;

				/* id */
				if (hier.obj[i]->ncoll == 1) {
					star[knew].id = hier.obj[i]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, star[knew].id, &bi);
					cp_SEvars_to_newstar(oldk, bi, knew);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->id[0], &bi);
					cp_SEvars_to_newstar(oldk, bi, knew);
					cp_m_to_newstar(oldk, bi, knew);
					nmerged = 1;
					while (nmerged < hier.obj[i]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar);
						cp_m_to_star(oldk, bi, &tempstar);
						merge_two_stars(&(star[knew]), &tempstar, &(star[knew]), vs);
                                                /* Owing to merger only useful vs's are v[1-3] */
						star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
						vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
						//star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knew);
					
					star[knew].id = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge1: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,star[knew].id);
					}
					
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						star[knew].id, star[knew].m, 
						star[knew].r,
						*(hier.obj[i]), k, kp, star[knew].se_k);
				}
				
				star[knew].rad = star[knew].se_radius * RSUN / units.l;

				/* track binding energy */
				BEf += -star[knew].Eint;
			} else if (hier.obj[i]->n == 2) {
				/* binary star */
				/* semimajor axis and eccentricity */
				binary[star[knew].binind].a = hier.obj[i]->a * cmc_units.l;
				binary[star[knew].binind].e = hier.obj[i]->e;
				
				/* masses */
				binary[star[knew].binind].m1 = hier.obj[i]->obj[0]->m * cmc_units.m / madhoc;
				binary[star[knew].binind].m2 = hier.obj[i]->obj[1]->m * cmc_units.m / madhoc;
				
				/* internal energies */
				binary[star[knew].binind].Eint1 = hier.obj[i]->obj[0]->Eint * cmc_units.E;
				binary[star[knew].binind].Eint2 = hier.obj[i]->obj[1]->Eint * cmc_units.E;
				
				/* id's */
				if (hier.obj[i]->obj[0]->ncoll == 1) {
					binary[star[knew].binind].id1 = hier.obj[i]->obj[0]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, binary[star[knew].binind].id1, &bi);
					cp_SEvars_to_newbinary(oldk, bi, knew, 0);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->obj[0]->id[0], &bi);
					cp_SEvars_to_star(oldk, bi, &tempstar);
					cp_m_to_star(oldk, bi, &tempstar);
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[0]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[0]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
						merge_two_stars(&tempstar, &tempstar2, &tempstar, vs);
						/* FIXME: really we're supposed to add the kick to each binary
						   member separately, then calculate the systemic kick to the binary,
						   but hopefully this doesn't happen too much. */
                                                /* The kick routine within /bse_wrap/bse/ correctly updates COM velocity... */
						if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
							wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
								sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
						}
						star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
						vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
						//star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knew);

					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 0);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 0);

					binary[star[knew].binind].id1 = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge2: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id1);
					}
					
                                        /* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id1,
						binary[star[knew].binind].m1, 
						star[knew].r,
						*(hier.obj[i]->obj[0]), k, kp, binary[star[knew].binind].bse_kw[0]);

                                        if (binary[star[knew].binind].m1==0.) {
                                          dprintf("Zero mass remnant! Parameters: knew=%li, binind=%li, kw[0]=%i, kw[1]=%i\n",
                                              knew, star[knew].binind, binary[star[knew].binind].bse_kw[0], 
                                              binary[star[knew].binind].bse_kw[1]);
                                        }
                                        
				}
				if (hier.obj[i]->obj[1]->ncoll == 1) {
					binary[star[knew].binind].id2 = hier.obj[i]->obj[1]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, binary[star[knew].binind].id2, &bi);
					cp_SEvars_to_newbinary(oldk, bi, knew, 1);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->obj[1]->id[0], &bi);
					cp_SEvars_to_star(oldk, bi, &tempstar);
					cp_m_to_star(oldk, bi, &tempstar);
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[1]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[1]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
						merge_two_stars(&tempstar, &tempstar2, &tempstar, vs);
						/* FIXME: really we're supposed to add the kick to each binary
						   member separately, then calculate the systemic kick to the binary,
						   but hopefully this doesn't happen too much. */
						if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
							wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
								sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
						}
						star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
						vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
						//star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 1);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 1);

					binary[star[knew].binind].id2 = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge3: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id2);
					}
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id2,
						binary[star[knew].binind].m2, 
						star[knew].r,
						*(hier.obj[i]->obj[1]), k, kp, binary[star[knew].binind].bse_kw[1]);
                                        if (binary[star[knew].binind].m2==0.) 
                                          dprintf("Zero mass remnant! Parameters: knew=%li, binind=%li, kw[0]=%i, kw[1]=%i\n",
                                              knew, star[knew].binind, binary[star[knew].binind].bse_kw[0], 
                                              binary[star[knew].binind].bse_kw[1]);
				}
				
				star[knew].m = binary[star[knew].binind].m1 + binary[star[knew].binind].m2;

				/* radii, and tb */
				binary[star[knew].binind].rad1 = binary[star[knew].binind].bse_radius[0] * RSUN / units.l;
				binary[star[knew].binind].rad2 = binary[star[knew].binind].bse_radius[1] * RSUN / units.l;
				binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
				/* track binding energy */
				BEf += binary[star[knew].binind].m1 * binary[star[knew].binind].m2 * sqr(madhoc) 
					/ (2.0 * binary[star[knew].binind].a) 
					- binary[star[knew].binind].Eint1 - binary[star[knew].binind].Eint2;
                                compress_binary(&star[knew], &binary[star[knew].binind]);
			} else if (hier.obj[i]->n == 3) {
				/******************************************/
				/* break triple by shrinking inner binary */
				/******************************************/
				/* put triple at origin */
				hier.obj[i]->x[0] = 0.0;
				hier.obj[i]->x[1] = 0.0;
				hier.obj[i]->x[2] = 0.0;
				
				/* and set to zero velocity */
				hier.obj[i]->v[0] = 0.0;
				hier.obj[i]->v[1] = 0.0;
				hier.obj[i]->v[2] = 0.0;
				
				/* trickle down properties */
				fb_downsync(hier.obj[i], t);
				fb_downsync(hier.obj[i]->obj[bid], t);
				
				/* temporarily store triple's stars' information so that 
				   we can calculate the triple's energy */
				threeobjs[0] = *(hier.obj[i]->obj[sid]);
				threeobjs[1] = *(hier.obj[i]->obj[bid]->obj[0]);
				threeobjs[2] = *(hier.obj[i]->obj[bid]->obj[1]);
				
				/* bring outer member of triple to zero energy, decreasing inner binary's semimajor axis
				   in the process, but preserving its eccentricity */
				hier.obj[i]->obj[bid]->a = -(hier.obj[i]->obj[bid]->obj[0]->m)*(hier.obj[i]->obj[bid]->obj[1]->m)/
					(2.0 * (fb_ketot(threeobjs, 3) + fb_petot(threeobjs, 3)));
				
				/********************************/
				/* set single star's properties */
				/********************************/
				/* internal energy */
				star[knewp].Eint = hier.obj[i]->obj[sid]->Eint * cmc_units.E;

				/* id */
				if (hier.obj[i]->obj[sid]->ncoll == 1) {
					star[knewp].id = hier.obj[i]->obj[sid]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, star[knewp].id, &bi);
					cp_SEvars_to_newstar(oldk, bi, knewp);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->obj[sid]->id[0], &bi);
					cp_SEvars_to_newstar(oldk, bi, knewp);
					cp_m_to_newstar(oldk, bi, knewp);
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[sid]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[sid]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar);
						cp_m_to_star(oldk, bi, &tempstar);
						merge_two_stars(&(star[knewp]), &tempstar, &(star[knewp]), vs);
						star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
						vt_add_kick(&(star[knewp].vt),vs[1],vs[2]);
						//star[knewp].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knewp);
					
					star[knewp].id = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge4: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,star[knewp].id);
					}
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						star[knewp].id, star[knewp].m, 
						star[knewp].r,
						*(hier.obj[i]->obj[sid]), k, kp, star[knewp].se_k);
				}

				/* radius */
				star[knewp].rad = star[knewp].se_radius * RSUN / units.l;

				/***************************/
				/* set binary's properties */
				/***************************/
				/* semimajor axis and eccentricity */
				binary[star[knew].binind].a = hier.obj[i]->obj[bid]->a * cmc_units.l;
				binary[star[knew].binind].e = hier.obj[i]->obj[bid]->e;
				
				/* masses */
				binary[star[knew].binind].m1 = hier.obj[i]->obj[bid]->obj[0]->m * cmc_units.m / madhoc;
				binary[star[knew].binind].m2 = hier.obj[i]->obj[bid]->obj[1]->m * cmc_units.m / madhoc;
				
				/* internal energies */
				binary[star[knew].binind].Eint1 = hier.obj[i]->obj[bid]->obj[0]->Eint * cmc_units.E;
				binary[star[knew].binind].Eint2 = hier.obj[i]->obj[bid]->obj[1]->Eint * cmc_units.E;
				
				/* id's */
				if (hier.obj[i]->obj[bid]->obj[0]->ncoll == 1) {
					binary[star[knew].binind].id1 = hier.obj[i]->obj[bid]->obj[0]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, binary[star[knew].binind].id1, &bi);
					cp_SEvars_to_newbinary(oldk, bi, knew, 0);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[0], &bi);
					cp_SEvars_to_star(oldk, bi, &tempstar);
					cp_m_to_star(oldk, bi, &tempstar);
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[bid]->obj[0]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
						merge_two_stars(&tempstar, &tempstar2, &tempstar, vs);
						/* FIXME: really we're supposed to add the kick to each binary
						   member separately, then calculate the systemic kick to the binary,
						   but hopefully this doesn't happen too much. */
						if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
							wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
								sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
						}
						star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
						vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
						//star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 0);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 0);

					binary[star[knew].binind].id1 = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge5: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id1);
					}
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id1,
						binary[star[knew].binind].m1, 
						star[knew].r,
						*(hier.obj[i]->obj[bid]->obj[0]), k, kp, binary[star[knew].binind].bse_kw[0]);
				}
				if (hier.obj[i]->obj[bid]->obj[1]->ncoll == 1) {
					binary[star[knew].binind].id2 = hier.obj[i]->obj[bid]->obj[1]->id[0];
					/* copy SE variables over to new star from old single star or binary member */
					oldk = binint_get_indices(k, kp, binary[star[knew].binind].id2, &bi);
					cp_SEvars_to_newbinary(oldk, bi, knew, 1);
				} else {
					/* merge progenitor stars one by one */
					oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[0], &bi);
					cp_SEvars_to_star(oldk, bi, &tempstar);
					cp_m_to_star(oldk, bi, &tempstar);
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[bid]->obj[1]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[nmerged], &bi);
						nmerged++;
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
						merge_two_stars(&tempstar, &tempstar2, &tempstar, vs);
						/* FIXME: really we're supposed to add the kick to each binary
						   member separately, then calculate the systemic kick to the binary,
						   but hopefully this doesn't happen too much. */
						if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
							wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
								sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
						}
						star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
						vt_add_kick(&(star[knew].vt),vs[1],vs[2]);
						//star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 1);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 1);

					binary[star[knew].binind].id2 = star_get_id_new();
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id2,
						binary[star[knew].binind].m2, 
						star[knew].r,
						*(hier.obj[i]->obj[bid]->obj[1]), k, kp, binary[star[knew].binind].bse_kw[1]);
				}
				
				star[knew].m = binary[star[knew].binind].m1 + binary[star[knew].binind].m2;

				/* radii and tb */
				binary[star[knew].binind].rad1 = binary[star[knew].binind].bse_radius[0] * RSUN / units.l;
				binary[star[knew].binind].rad2 = binary[star[knew].binind].bse_radius[1] * RSUN / units.l;
				binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;
				
				/* track binding energy */
				BEf += binary[star[knew].binind].m1 * binary[star[knew].binind].m2 * sqr(madhoc) / 
					(2.0 * binary[star[knew].binind].a) 
					- binary[star[knew].binind].Eint1 - binary[star[knew].binind].Eint2
					- star[knewp].Eint;
                                compress_binary(&star[knew], &binary[star[knew].binind]);
			}
		}
		
		/* destroy two progenitors */
		destroy_obj(k);
		destroy_obj(kp);

		/* update binding energy */
		if (isbinbin) {
			E_bb += BEf - BEi;
			DE_bb += BEf - BEi;
		} else {
			E_bs += BEf - BEi;
			DE_bs += BEf - BEi;
		}
	}
	
	/* free Fewbody memory */
	fb_free_hier(hier);
}

/* simulate relaxation to set timestep */
double simul_relax(gsl_rng *rng)
{
	long si, k, p=AVEKERNEL, N_LIMIT, simin, simax;
	double dt, dtmin=GSL_POSINF, W, n_local;
	double Mv2ave, Mave, M2ave, sigma;
	
	N_LIMIT = clus.N_MAX;

	/* calculate sliding average timesteps */
	p = MAX((long) (1.0e-4 * ((double) clus.N_STAR) / 2.0), AVEKERNEL);
	for (si=1; si<=N_LIMIT; si++) {
		simin = si - p;
		simax = simin + (2 * p - 1);
		if (simin < 1) {
			simin = 1;
			simax = simin + (2 * p - 1);
		} else if (simax > N_LIMIT) {
			simax = N_LIMIT;
			simin = simax - (2 * p - 1);
		}

		Mv2ave = 0.0;
		Mave = 0.0;
		M2ave = 0.0;
		for (k=simin; k<=simax; k++) {
			Mv2ave += star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt));
			Mave += star[k].m * madhoc;
			M2ave += sqr(star[k].m * madhoc);
		}
		Mv2ave /= (double) (2 * p);
		Mave /= (double) (2 * p);
		M2ave /= (double) (2 * p);
		
		/* sigma is the 3D velocity dispersion */
		sigma = sqrt(Mv2ave/Mave);
		/* average relative speed for a Maxwellian, from Binney & Tremaine */
		W = 4.0 * sigma / sqrt(3.0 * PI);

		/* Compute local density */
		n_local = calc_n_local(si, AVEKERNEL, clus.N_MAX);
		
		/* remember that code time units are t_cross * N/log(GAMMA*N) */
		/* this expression is from Freitag & Benz (2001), eqs. (8) and (9), we're just
		   inputting locally-averaged quantities */
		dt = sqr(2.0*THETASEMAX/PI) * (PI/32.0) * 
			cub(W) / ( ((double) clus.N_STAR) * n_local * (4.0 * M2ave) );

		// Meagan: shorten timestep while in core-collapsed state
		//if (rho_core > 50.0) {
		//if (rho_core > 35.0) {
		//	dt = sqr(2.0*0.707106781*THETASEMAX/PI) * (PI/32.0) *
		//	dt = sqr(2.0*0.5*THETASEMAX/PI) * (PI/32.0) *
		//		cub(W) / ( ((double) clus.N_STAR) * n_local * (4.0 * M2ave) );
		//}
		//if (N_core < 100) {
		//	dt = sqr(2.0*0.707106781*THETASEMAX/PI) * (PI/32.0) *
		//		cub(W) / ( ((double) clus.N_STAR) * n_local * (4.0 * M2ave) );
		//}


		dtmin = MIN(dtmin, dt);
	}

	return(dtmin);
}

/* Since the binary interactions are done in a vaccuum, it is possible for them
   to produce pathologically wide binaries, which must be broken by hand, else 
   they shorten the timestep to a crawl. */
void break_wide_binaries(void)
{
	long j, k, knew, knewp;
	double W, vorb, Eexcess=0.0, exc_ratio, nlocal, llocal;
	
	for (k=1; k<=clus.N_MAX_NEW; k++) {
		if (star[k].binind) {
			/* binary index */
			j = star[k].binind;
			
			/* get relative velocity from velocity dispersion at binary's radial position */
			W = 4.0 * sigma_r(star[k].r) / sqrt(3.0 * PI);
			
			/* this is an order of magnitude estimate for the orbital speed */
			vorb = sqrt(star[k].m * madhoc / binary[j].a);

			nlocal = calc_n_local(k, AVEKERNEL, clus.N_MAX);
			llocal = 0.1 * pow(nlocal, -1.0/3.0);

			/* Destroy binary if its orbital speed is less than some fraction of the 
			   local relative velocity. */
			/* if (vorb <= XHS*W) {*/
			/* break if apocenter is larger than interparticle separation */
			if (binary[j].a*(1.0+binary[j].e) >= llocal) {
				dprintf("breaking wide binary: vorb=%g W=%g\n", vorb, W);
				
				Eexcess += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);

				/* create two stars for the binary components */
				knew = create_star();
				knewp = create_star();
				
				cp_binmemb_to_star(k, 0, knew);
				cp_binmemb_to_star(k, 1, knewp);
				
				/* destroy this binary */
				destroy_obj(k);
			} else {
				/* take excess energy from nearby field star (single or binary) */
				if(Eexcess > 0 && star[k].interacted == 0 && Eexcess < 0.5*(sqr(star[k].vt)+sqr(star[k].vr))*star[k].m*madhoc) {
					exc_ratio = 
						sqrt( (sqr(star[k].vt)+sqr(star[k].vr)-2.0*Eexcess/(star[k].m*madhoc))/
						      (sqr(star[k].vt)+sqr(star[k].vr)) );
					star[k].vr *= exc_ratio;
					star[k].vt *= exc_ratio;
					set_star_EJ(k);
					Eexcess = 0.0;
					star[k].interacted = 1;
				}
			}
		}
	}
	
	/* keep track of the energy that's vanishing due to our negligence */
	Eoops += -Eexcess;
}

/* calculate and store the velocity dispersion profile */
void calc_sigma_r(void)
{
	long si, k, p=AVEKERNEL, N_LIMIT, simin, simax, siminlast, simaxlast;
	double Mv2ave, Mave;
	
	N_LIMIT = clus.N_MAX;
	sigma_array.n = N_LIMIT;

	/* p = MAX((long) (1.0e-4 * ((double) clus.N_STAR) / 2.0), AVEKERNEL); */
	siminlast = 1;
	simaxlast = 0;
	Mv2ave = 0.0;
	Mave = 0.0;
	for (si=1; si<=N_LIMIT; si++) {
		// determine appropriate bounds for summing
		simin = si - p;
		simax = simin + (2 * p - 1);
		if (simin < 1) {
			simin = 1;
			simax = simin + (2 * p - 1);
		} else if (simax > N_LIMIT) {
			simax = N_LIMIT;
			simin = simax - (2 * p - 1);
		}

		// do sliding sum
		for (k=siminlast; k<simin; k++) {
			Mv2ave -= star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt));
			Mave -= star[k].m * madhoc;
		}

		for (k=simaxlast+1; k<=simax; k++) {
			Mv2ave += star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt));
			Mave += star[k].m * madhoc;
		}
		// don't need to average since one gets divided by the other
		//Mv2ave /= (double) (2 * p);
		//Mave /= (double) (2 * p);
		
		/* store sigma (sigma is the 3D velocity dispersion) */
		sigma_array.r[si] = star[si].r;
		sigma_array.sigma[si] = sqrt(Mv2ave/Mave);
		
		siminlast = simin;
		simaxlast = simax;
	}
}

double calc_average_mass_sqr(long index, long N_LIMIT) {
  long simin, simax, si, p, k;
  double M2ave;

  si= index;

  /* calculate sliding average timesteps */
  p = MAX((long) (1.0e-4 * ((double) clus.N_STAR) / 2.0), AVEKERNEL);
  simin = si - p;
  simax = simin + (2 * p - 1);
  if (simin < 1) {
    simin = 1;
    simax = simin + (2 * p - 1);
  } else if (simax > N_LIMIT) {
    simax = N_LIMIT;
    simin = simax - (2 * p - 1);
  }

  M2ave = 0.0;
  for (k=simin; k<=simax; k++) {
    M2ave += sqr(star[k].m * madhoc);
  }
  M2ave /= (double) (2 * p);

  return(M2ave);
};

/* generic routine for testing for approximate equality of floating point numbers */
double floateq(double a, double b) {
	double diff=a-b, mean=0.5*(a+b);

	if (a == 0.0 && b == 0.0) {
		return(1);
	} else if (fabs(diff)/mean < 1.0e-6) {
		return(1);
	} else {
		return(0);
	}
}

/* tidal capture (including merger) cross section; inputs are assumed to be in (self-consistent) code units */
double sigma_tc_nd(double n, double m1, double r1, double m2, double vinf) {
	double a, beta=2.2, vstar1=sqrt(2.0*m1/r1);

	if (floateq(n, 1.5)) {
		a = 6.60 * pow(m2/m1, 0.242) + 5.06 * pow(m2/m1, 1.33);
	} else if (floateq(n, 3.0)) {
		a = 3.66 * pow(m2/m1, 0.200) + 2.94 * pow(m2/m1, 1.32);
	} else {
		eprintf("unknown polytropic index n=%g!\n", n);
		exit_cleanly(-1);
		exit(1);
	}
	
	return(a*pow(vinf/vstar1,-beta)*r1*r1);
}

/* tidal capture (including merger) cross section; inputs are assumed to be in (self-consistent) code units */
double sigma_tc_nn(double na, double ma, double ra, double nb, double mb, double rb, double vinf) {
	double n1, m1, r1, n2, m2, r2;
	double a, beta=2.2, gamma, vstar1;

	/* make sure m2 >= m1 */
	if (mb >= ma) {
		n1 = na;
		m1 = ma;
		r1 = ra;
		n2 = nb;
		m2 = mb;
		r2 = rb;
	} else {
		n1 = nb;
		m1 = mb;
		r1 = rb;
		n2 = na;
		m2 = ma;
		r2 = ra;
	}
	
	gamma=log(r2/r1)/log(m2/m1);
	vstar1=sqrt(2.0*m1/r1);

	if (floateq(n1, 1.5) && floateq(n2, 1.5)) {
		a = 6.05 * pow(m2/m1, 0.835*log(gamma)+0.468) + 6.50 * pow(m2/m1, 0.563*log(gamma)+1.75);
	} else if (floateq(n1, 3.0) && floateq(n2, 3.0)) {
		a = 3.50 * pow(m2/m1, 0.814*log(gamma)+0.551) + 3.53 * pow(m2/m1, 0.598*log(gamma)+1.80);
	} else if (floateq(n1, 1.5) && floateq(n2, 3.0)) {
		a = 7.98 * pow(m2/m1, -1.23*log(gamma)-0.232) + 3.57 * pow(m2/m1, 0.625*log(gamma)+1.81);
	} else if (floateq(n1, 3.0) && floateq(n2, 1.5)) {
		/* this combo is not presented in Kim & Lee (1999), so we'll use the largest a to make sure we
		   capture everything */
		a = 6.05 * pow(m2/m1, 0.835*log(gamma)+0.468) + 6.50 * pow(m2/m1, 0.563*log(gamma)+1.75);
	} else {
		eprintf("unknown polytropic indexes n1=%g n2=%g!\n", n1, n2);
		exit_cleanly(-1);
		exit(1);
	}
	
	return(a*pow(vinf/vstar1,-beta)*r1*r1);
}

/* T_l function for use in tidal capture calculations */
double Tl(int order, double polytropicindex, double eta)
{
	int l=order;
	double n=polytropicindex, x=log10(eta), x2=x*x, x3=x*x2, x4=x2*x2, x5=x2*x3;

	if (l != 2 && l != 3) {
		eprintf("unknown order l=%d\n", l);
		exit_cleanly(-1);
		exit(1);
	}

	if (floateq(n, 1.5)) {
		if (l == 2) {
			return(pow(10.0, -0.397 + 1.678*x + 1.277*x2 - 12.42*x3 + 9.446*x4 - 5.550*x5));
		} else {
			return(pow(10.0, -0.909 + 1.574*x + 12.37*x2 - 57.40*x3 + 80.10*x4 - 46.43*x5));	
		}
	} else if (floateq(n, 2.0)) {
		if (l == 2) {
			return(pow(10.0, -0.517 - 0.906*x + 23.88*x2 - 93.49*x3 + 112.3*x4 - 44.15*x5));
		} else {
			return(pow(10.0, -1.040 - 1.354*x + 37.64*x2 - 139.9*x3 + 168.2*x4 - 66.53*x5));	
		}
	} else if (floateq(n, 3.0)) {
		if (l == 2) {
			return(pow(10.0, -1.124 + 0.877*x - 13.37*x2 + 21.55*x3 - 16.48*x4 + 4.124*x5));
		} else {
			return(pow(10.0, -1.703 + 2.653*x - 14.34*x2 + 12.85*x3 - 0.492*x4 - 3.600*x5));	
		}
	} else {
		eprintf("unknown polytropic index n=%g\n", n);
		exit_cleanly(-1);
		exit(1);
	}
}

/* tidal energy of Mpert acting on Mosc, in the polytropic approximation */
double Etide(double rperi, double Mosc, double Rosc, double nosc, double Mpert)
{
	double eta=sqrt(Mosc/(Mosc+Mpert))*pow(rperi/Rosc, 1.5);

	return(sqr(Mpert)/Rosc * (pow(Rosc/rperi, 6.0) * Tl(2, nosc, eta) + pow(Rosc/rperi, 8.0) * Tl(3, nosc, eta)));
}

/* Add kicks to vt */
 void vt_add_kick(double *vt, double vs1, double vs2)
{
	double X, theta, vtx, vty;
	X = rng_t113_dbl();
	theta = 2.0 * PI * X;
	vtx = *vt*sin(theta);
	vty = *vt*cos(theta);
	vtx = vtx + vs1*1.0e5/(units.l/units.t);
	vty = vty + vs2*1.0e5/(units.l/units.t);
	*vt = sqrt(vtx*vtx+vty*vty);
}


