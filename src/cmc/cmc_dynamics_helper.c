/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cmc.h"
#include "cmc_vars.h"


/**
* @brief Zero'es out all star variables corresponding to given index
*
* @param j star index
*/
void zero_star(long j)
{
	int g_j = get_global_idx(j);
	star_r[g_j] = 0.0;
	star_m[g_j] = 0.0;
	star_phi[g_j] = 0.0;
	//OPT: Try using memset or bzero for better readable code.
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
        star[j].se_B_0 = 0.0; /* PK */
        star[j].se_bacc = 0.0;
        star[j].se_tacc = 0.0;
        star[j].se_scm_B = 0.0;
        star[j].se_scm_formation = 0.0;
	star[j].se_epoch = 0.0;
	star[j].se_tphys = 0.0;
	star[j].se_radius = 0.0;
	star[j].se_lum = 0.0;
	star[j].se_mc = 0.0;
	star[j].se_rc = 0.0;
	star[j].se_menv = 0.0;
	star[j].se_renv = 0.0;
	star[j].se_tms = 0.0;
	star[j].se_bhspin = 0.0;
}

/**
* @brief Zero'es out all variables of the binary corresponding to the given index
*
* @param j index of binary
*/
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
}

/**
* @brief prints the specified interaction status text
*
* @param status_text[] status text
*/
void print_interaction_status(char status_text[])
{
	gprintf(" %s", status_text);
	if (!quiet) {
		fflush(stdout);
	}
	parafprintf(logfile, " %s", status_text);
}

/**
* @brief prints error message during interaction
*/
void print_interaction_error(void)
{
	gprintf("!");
	if (!quiet) {
		fflush(stdout);
	}
	parafprintf(logfile, "!");
}

/**
* @brief destroy an object (can be star or binary)
*
* @param i Index of star to be destroyed
*/
void destroy_obj(long i)
{
	double r, phi;

	if (star[i].binind) {
		dprintf("Star %ld destroyed with id %ld, binary binind=%ld id1=%ld id2=%ld!\n", i, star[i].id, star[i].binind, binary[star[i].binind].id1, binary[star[i].binind].id2);
		destroy_binary(star[i].binind);
	} else {
		dprintf("Star %ld destroyed with id %ld!\n", i, star[i].id);
	}

	/* need to zero out E's, J, but can't zero out potential---this is the easiest way */
	int g_i = get_global_idx(i);
	r = star_r[g_i];
	phi = star_phi[g_i];

	zero_star(i);

	star_r[g_i] = r;
	star_phi[g_i] = phi;

	remove_star_center(i);
}

/**
* @brief destroy an object (can be star or binary) - pure serial version
*
* @param i index of object/star to be destroyed
*/
void destroy_obj_new(long i)
{
	if (star[i].binind) {
		dprintf("Star %ld destroyed, binary %ld!\n", i, star[i].binind);
		destroy_binary(star[i].binind);
	}

	zero_star(i);

	star[i].r = SF_INFINITY;	/* send star to infinity */
	star[i].m = DBL_MIN;		/* set mass to very small number */
	star[i].vrnew = 0.0;		/* setup vr and vt for           */
	star[i].vtnew = 0.0;		/*		future calculations  */
}

/**
* @brief destroy a binary
*
* @param i index of binary (in the binary array) to be destroyed
*/
void destroy_binary(long i)
{
	/* set inuse flag to zero, and zero out all other properties for safety */
	zero_binary(i);
	dprintf("Binary %ld with id1=%ld id2=%ld destroyed!\n", i, binary[i].id1, binary[i].id2);
	N_b_local--;
}

/**
* @brief create a new star, returning its index
*
* @param idx index of star that creates new star
* @param dyn_0_se_1 if 0, created by dynamics, if 1 created by stellar evolution
*
* @return index of new star
*/
long create_star(int idx, int dyn_0_se_1)
{
	long i;

	/* account for new star */
	clus.N_STAR_NEW++;
	clus.N_MAX_NEW++;
	

	i = clus.N_MAX_NEW;

	dprintf("star created!!!, idx=%ld by star idx=%d\ton node=%d,\tdyn_or_se=%d\t\n", i, idx, myid, dyn_0_se_1);

	/* initialize to zero for safety */
	zero_star(i);
	
	return(i);
}

/**
* @brief create a new binary, returning its index. A new binary creation is done as follows. Since binaries get destroyed, so, first we scan the existing binary array and look for holes i.e. which were left behind by destroyed stars. If one is found, we insert a new binary there, if not we insert it at the end.
*
* @param idx index of the star that is creating the binary
* @param dyn_0_se_1 0 if created by dynamics, 1 if created by stellar evolution
*
* @return index of new binary
*/
long create_binary(int idx, int dyn_0_se_1)
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
		exit_cleanly(1, __FUNCTION__);
	}

	dprintf("HOLE FOUND at %ld\t INSERTING STAR %d\n", i, idx);

	/* account for new binary */
	++N_b_local;
	
	/* initialize to zero for safety */
	zero_binary(i);

	/* mark binary as being in use */
	binary[i].inuse = 1;

	/* create the star that points to the binary */
	j = create_star(idx, dyn_0_se_1);
	
	star[j].binind = i;
	
	dprintf("Binary Created on node %d!! single star idx = %ld binind = %ld\n", myid, j, i);

	return(j);
}



/**
* @brief Sort by mass: k1, k2, k3 in order of most massive to least massive. (parallel version of sort_three_masses)
*
* @param sq index of the first of the three consecutive stars
* @param k1 variable to store the most massive star's index
* @param k2 variable to store the 2nd most massive star's index
* @param k3 variable to store the least massive star's index
*/
void sort_three_masses(long sq, long *k1, long *k2, long *k3) {
/* Sort by mass: k1, k2, k3 in order of most massive to least massive. */

    long g_sq0 = get_global_idx(sq);
    long g_sq1 = get_global_idx(sq+1);
    long g_sq2 = get_global_idx(sq+2);

	//parafprintf(threebbfile, "beginning of sort_three_masses function" );
	if ((star_m[g_sq0] >= star_m[g_sq1]) && (star_m[g_sq0] >= star_m[g_sq2])) {
		*k1 = sq;
		if (star_m[g_sq1] >= star_m[g_sq2]) {
			*k2 = sq + 1;
			*k3 = sq + 2;
		}
		else {
			*k2 = sq + 2;
			*k3 = sq + 1;
		}
	}
	else if ((star_m[g_sq1] >= star_m[g_sq0]) && (star_m[g_sq1] >= star_m[g_sq2])) {
		*k1 = sq + 1;
		if (star_m[g_sq0] >= star_m[g_sq2]) {
			*k2 = sq;
			*k3 = sq + 2;
		}
		else {
			*k2 = sq + 2;
			*k3 = sq;
		}
	}
	else if ((star_m[g_sq2] >= star_m[g_sq0]) && (star_m[g_sq2] >= star_m[g_sq1])) {
		*k1 = sq + 2;
		if (star_m[g_sq0] > star_m[g_sq1]) {
			*k2 = sq;
			*k3 = sq + 1;
		}
		else {
			*k2 = sq + 1;
			*k3 = sq;

		}
	}
	//parafprintf(threebbfile, "end of sort_three_masses function" );
}

/**
* @brief
	   Function for calculating quantities relevant for three-body binary formation
		Generates random direction for vt for all stars (so then have 3D veloc. for each)
		Then, for the two possible binaries that could form (1,2 or 1,3) we calculate:
		(1) COM veloc of candidate binary pair  (2) v3: veloc of 3rd star relative to binary
		In COM frame of the system of three stars, calculate total momentum and energy.
*
* @param k1 index of 1st star
* @param k2 index of 2nd star
* @param k3 index of 3rd star
* @param v1[4] ?
* @param v2[4] ?
* @param v3[4] ?
* @param vrel12 ?
* @param vrel3 ?
* @param rng gsl rng
*/
void calc_3bb_encounter_dyns(long k1, long k2, long k3, double v1[4], double v2[4], double v3[4], double (*vrel12)[4], double (*vrel3)[4], gsl_rng *rng) {
	long j;
	double vcm12[4];

	/* set random angle between vt's */
	/* first store random variable */
//	star[k].Y = rng_t113_dbl_new(curr_st);
//	star[kp].Y = star[k].Y;
	double angle1 = rng_t113_dbl_new(curr_st) * 2.0 * PI;
	double angle2 = rng_t113_dbl_new(curr_st) * 2.0 * PI;

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

        vcm12[j] = (star_m[get_global_idx(k1)] * v1[j] + star_m[get_global_idx(k2)] * v2[j])/(star_m[get_global_idx(k1)] + star_m[get_global_idx(k2)]);
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


/**
* @brief ?
*
* @param eta_min ?
* @param k1 index of star 1
* @param k2 index of star 2
* @param k3 index of star 3
* @param vrel12[4] ?
* @param vrel3[4] ?
*
* @return ?
*/
double get_eta(double eta_min, long k1, long k2, long k3, double vrel12[4], double vrel3[4])
{
	double kk, eta, norm, area_max, area, eta_test, comp_ymax, comp_y, true_y;
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

	kk = 2.0 * ( (star_m[get_global_idx(k1)] + star_m[get_global_idx(k2)] + star_m[get_global_idx(k3)]) / (star_m[get_global_idx(k1)] + star_m[get_global_idx(k2)]) ) * sqr(vrel12[0] / vrel3[0]);
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
		area = rng_t113_dbl_new(curr_st)*area_max;
		/* find the eta that corresponds to chosen area */
		/* test coordinate 1 (x) */
	//	eta_test = pow(((-2.0*area/(13.0*norm)) + pow(eta_min, -4.0)), -0.25);

	// TODO: figure out this negative sign in front of 2 below - where did it come from??
		eta_test = pow(((-7.0*area/(2.0*norm*kk)) + pow(eta_min, -3.5)), (-2.0/7.0));
		/* now choose second test coordinate (y) by choosing a number between zero and the value of the comparison function at eta_test, comp_ymax */
		comp_ymax = kk * norm * pow(eta_test, -4.5);
		comp_y = rng_t113_dbl_new(curr_st)*comp_ymax;
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


/**
* @brief ?
*
* @param P_3bb ?
* @param k1 index of star 1
* @param k2 index of star 2
* @param k3 index of star 3
* @param form_binary ?
* @param eta_min ?
* @param ave_local_mass ?
* @param n_local ?
* @param sigma_local ?
* @param v1[4] ?
* @param v2[4] ?
* @param v3[4] ?
* @param vrel12[4] ?
* @param vrel3[4] ?
* @param delta_E_3bb ?
* @param rng gsl rng
*/
void make_threebodybinary(double P_3bb, long k1, long k2, long k3, long form_binary, double eta_min, double ave_local_mass, double n_local, double sigma_local, double v1[4], double v2[4], double v3[4], double vrel12[4], double vrel3[4], double delta_E_3bb, gsl_rng *rng)
{
	double PE_i, PE_f, KE_i, KE_f, KE_cmf_i, KE_cmf_f, delta_PE, delta_KE, delta_E, binary_cm;
	double m1, m2, m3;
//  double ms, mb;
	double r1, r2, r3;
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
    m1=star_m[get_global_idx(k1)];
    m2=star_m[get_global_idx(k2)];
    m3=star_m[get_global_idx(k3)];
    r1=star_r[get_global_idx(k1)];
    r2=star_r[get_global_idx(k2)];
    r3=star_r[get_global_idx(k3)];

	// Initial energy
    PE_i = madhoc*(m1*star_phi[get_global_idx(k1)] + m2*star_phi[get_global_idx(k2)] + m3*star_phi[get_global_idx(k3)]);
	KE_i = 0.5*madhoc*(m1 * sqr(v1[0]) + m2 * sqr(v2[0]) + m3 * sqr(v3[0]));

	// COM of candidate binary pair
	binary_cm = (m1 * r1 + m2 * r2)/(m1 + m2);

	for (j=1; j<=3; j++) {
		cm_vel[j] = (m1 * v1[j] + m2 * v2[j] + m3 * v3[j])/(m1 + m2 + m3);
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
	KE_cmf_i = 0.5*madhoc*(m1 * sqr(v1_cmf[0]) + m2 * sqr(v2_cmf[0]) + m3 * sqr(v3_cmf[0]));

	// COMPUTE NEW QUANTITIES
		// Binary orbital properties: 
		//	binding energy: choose eta, then Eb = eta * <m> / sigma^2
		//	eccentricity: chosen from thermal distribution
	eta = get_eta(eta_min, k1, k2, k3, vrel12, vrel3);
	Eb = -0.5*eta * ave_local_mass * sqr(sigma_local); //Note: ave_local_mass is already multiplied by madhoc
	// **Note binding energy is NEGATIVE
	ecc_max = 1.0;
	ecc = ecc_max * sqrt(rng_t113_dbl_new(curr_st)); 					
	r_p = m1 * m2 * madhoc / (eta * sqr(sigma_local)); //note, for units, have madhoc^2 in numerator, and madhoc in the denominator ====> single power madhoc in numerator
	semi_major = m1 * m2 * sqr(madhoc) / (-2*Eb);
	// Using cons. of momentum and energy, can calculate scalar velocities of binary and single
	vs_cmf[0] = sqrt((2.0*(KE_cmf_i-Eb)*((m1 + m2)/m3))/((madhoc)*(m1 + m2 + m3)));
	vb_cmf[0] = 1.0*vs_cmf[0]*m3/(m1 + m2);

	// Choose random direction for motion of single star
	angle3 = rng_t113_dbl_new(curr_st) * PI; // polar angle - [0, PI)
	angle4 = rng_t113_dbl_new(curr_st) * 2.0 * PI; // azimuthal angle - [0, 2*PI)
	// Compute vector velocities of single and binary
	vs_cmf[1] = vs_cmf[0] * sin(angle3) * cos(angle4);
	vs_cmf[2] = vs_cmf[0] * sin(angle3) * sin(angle4);
	vs_cmf[3] = vs_cmf[0] * cos(angle3);

	/* Fix direction of velocity of binary to be opposite to that of single - add Pi to each of the angles 
		(polar angle and azimuthal) that were used to orient the velocity of the single  */
	vb_cmf[1] = vb_cmf[0] * sin(PI - angle3) * cos(angle4 + PI);
	vb_cmf[2] = vb_cmf[0] * sin(PI - angle3) * sin(angle4 + PI);
	vb_cmf[3] = vb_cmf[0] * cos(PI - angle3);

	vs_cmf[0] = sqrt(sqr(vs_cmf[1]) + sqr(vs_cmf[2]) + sqr(vs_cmf[3]));
	vb_cmf[0] = sqrt(sqr(vb_cmf[1]) + sqr(vb_cmf[2]) + sqr(vb_cmf[3]));

	/* Final COM energy - Binding energy of new binary plus kinetic energy of single and binary, 
		calculated using the magnitudes of their velocities (vs_cmf[0] and vb_cmf[0]) */
	KE_cmf_f = Eb + 0.5*madhoc*(m3*sqr(vs_cmf[0]) + (m1 + m2)*sqr(vb_cmf[0]));

	// Transform back to lab frame: final velocities of single and binary are vs and vb
	for (j=1; j<=3; j++) {
		vs[j] = vs_cmf[j] + cm_vel[j];
		vb[j] = vb_cmf[j] + cm_vel[j];	
	}
	// magnitudes
	vs[0] = sqrt(sqr(vs[1]) + sqr(vs[2]) + sqr(vs[3]));
	vb[0] = sqrt(sqr(vb[1]) + sqr(vb[2]) + sqr(vb[3]));

	// IF binary should NOT actually be formed, make sure stars are not marked as interacted,
	// so that they will interact in two-body loop
	if (form_binary == 0) { 
		star[k1].threebb_interacted = 0;
		star[k2].threebb_interacted = 0;
		star[k3].threebb_interacted = 0;

		parafprintf(lightcollisionfile, "%.16g %ld %ld %ld %ld %ld %ld %g %g %g %d %d %d %g %g %g %g %g %g %g\n", TotalTime, k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, m1*(units.m / clus.N_STAR / MSUN), m2 * (units.m / clus.N_STAR / MSUN), m3 *(units.m / clus.N_STAR / MSUN), star[k1].se_k, star[k2].se_k, star[k3].se_k, star[k1].rad * units.l / AU, star[k2].rad * units.l / AU, star[k3].rad * units.l / AU, Eb, ecc, semi_major * units.l / AU, r_p * units.l / AU);
                //append_to_table(LIGHTCOLLISION_ID, NFIELDS_LIGHT_COLLISION, TotalTime, k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, m1*(units.m / clus.N_STAR / MSUN), m2 * (units.m / clus.N_STAR / MSUN), m3 *(units.m / clus.N_STAR / MSUN),  star[k1].se_k,  star[k2].se_k,  star[k3].se_k,  star[k1].rad * units.l / AU,  star[k2].rad * units.l / AU,  star[k3].rad * units.l / AU,  Eb,  ecc,  semi_major * units.l / AU,  r_p * units.l / AU);
	}
	// IF binary IS TO BE FORMED, set star and binary properties for new binary, as well as properties of single star
	else if (form_binary == 1) {

		knew = create_binary(k1, 0); /* returns an unused binary id */

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
		binary[star[knew].binind].m1 = m1;
		binary[star[knew].binind].m2 = m2;
		binary[star[knew].binind].Eint1 = star[k1].Eint;
		binary[star[knew].binind].Eint2 = star[k2].Eint;
		binary[star[knew].binind].a = semi_major;
		binary[star[knew].binind].e = ecc;
	
		// star properties for the new binary
        star_r[get_global_idx(knew)] = binary_cm;
        star_m[get_global_idx(knew)] = m1 + m2;

		star[knew].vr = vb[3];
		star[knew].vt = sqrt(sqr(vb[1]) + sqr(vb[2]));
		star[knew].Eint = star[k1].Eint + star[k2].Eint;  // Set new internal energy for the binary, sum of Eint1 and Eint2 of the two components

	//	parafprintf(threebbfile, "new binary properties for object 'knew':  r=%g  vr=%g  vt=%g\n", star[knew].r, star[knew].vr, star[knew].vt);

		// Set new properties of single star
	//		star[k3].interacted = 1;
	//	star[k3].r = star[k3].r;  // later I can adjust the positions/potentials
		star[k3].vr = vs[3];
		star[k3].vt = sqrt(sqr(vs[1]) + sqr(vs[2]));

		// copy SE variables over to new binary from old single stars
		cp_SEvars_to_newbinary(k1, -1, knew, 0);
		cp_SEvars_to_newbinary(k2, -1, knew, 1);

		// radii, and tb : rhs was set in calls to cp_SEvars_to_newbinary()
		binary[star[knew].binind].rad1 = binary[star[knew].binind].bse_radius[0] * RSUN / units.l;
		binary[star[knew].binind].rad2 = binary[star[knew].binind].bse_radius[1] * RSUN / units.l;
		binary[star[knew].binind].bse_tb = sqrt(cub(binary[star[knew].binind].a * units.l / AU)/(binary[star[knew].binind].bse_mass[0]+binary[star[knew].binind].bse_mass[1]))*365.25;

		/* track binding energy */
		//BEf += binary[star[knew].binind].m1 * binary[star[knew].binind].m2 * sqr(madhoc)
		// Do I need to worry about this binding energy??

		star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]); // set potential for binary to be potential at its new cm location (r)
//		star_phi[get_global_idx(k3)] = potential(r3); // set potential for single to be same as original...should I adjust?
		PE_f = madhoc*(star_m[get_global_idx(knew)]*star_phi[get_global_idx(knew)] + m3*star_phi[get_global_idx(k3)]);
		KE_f = Eb + 0.5*madhoc*(star_m[get_global_idx(knew)] * sqr(vb[0]) + m3 * sqr(vs[0]));


		// In the right units for comparing to the Energies computed in cmc_utils.c ComputeEnergy(), and plotted in quick_cluster_plot
		delta_PE = (PE_f - PE_i);
		delta_KE = (KE_f - KE_i);
		delta_E = delta_KE + delta_PE;

		parafprintf(threebbdebugfile, "%ld %ld %ld %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, r1, r2, r3, m1, m2, m3, v1[0], v1[1], v1[2], v1[3], v2[0], v2[1], v2[2], v2[3], v3[0], v3[1], v3[2], v3[3], v1_cmf[0], v1_cmf[1], v1_cmf[2], v1_cmf[3], v2_cmf[0], v2_cmf[1], v2_cmf[2], v2_cmf[3], v3_cmf[0], v3_cmf[1], v3_cmf[2], v3_cmf[3]);

		parafprintf(threebbdebugfile, "%ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", knew, star[knew].id, star_r[get_global_idx(knew)], star_m[get_global_idx(knew)], vs_cmf[0], vs_cmf[1], vs_cmf[2], vs_cmf[3], vb_cmf[0], vb_cmf[1], vb_cmf[2], vb_cmf[3], vs[0], vs[1], vs[2], vs[3], vb[0], vb[1], vb[2], vb[3], ave_local_mass, n_local,sigma_local, eta, Eb, ecc, r_p, semi_major, PE_i, PE_f, KE_cmf_i, KE_cmf_f, KE_i, KE_f, delta_PE, delta_KE, delta_E);
		//parafprintf(threebbdebugfile, "Internal Energies:  Eint1=%g  Eint2=%g  Eint_binary=%g\n", star[k1].Eint, star[k2].Eint, star[knew].Eint);

		// running total for simulation - change in energy occuring during 3bb formation procedure
		delta_E_3bb += delta_E;

		// Set E and J for new objects
		set_star_EJ(knew);	// binary
		set_star_EJ(k3);	// single

		// Set binary properties
		binary[star[knew].binind].a = semi_major;
		binary[star[knew].binind].e = ecc;

		// destroy the two former single stars (which have now formed a binary)
		// leave the remaining single star (properties have already been updated)
		parafprintf(threebbfile, "%.16g %ld %ld %ld %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %ld\n", TotalTime, k1, k2, k3, star[k1].id, star[k2].id, star[k3].id, m1*(units.m / clus.N_STAR / MSUN), m2*(units.m / clus.N_STAR / MSUN), m3*(units.m / clus.N_STAR / MSUN), ave_local_mass, n_local, sigma_local, eta, Eb, binary[star[knew].binind].e, binary[star[knew].binind].a * units.l / AU, r_p * units.l / AU, star_r[get_global_idx(knew)], r3, star[knew].vr, star[knew].vt, star[k3].vr, star[k3].vt, star_phi[get_global_idx(knew)], star_phi[get_global_idx(k3)], delta_PE, delta_KE, delta_E, delta_E_3bb, N3bbformed);
		destroy_obj(k1);
		destroy_obj(k2);
	//	parafprintf(threebbfile, "KE_cmf_i=%g  eta=%g  Eb=%g\n", KE_cmf_i, eta, Eb);
	//	parafprintf(threebbfile, "vs_cmf[0]=%g  vb_cmf[0]=%g  cm_vel[0]=%g\n", vs_cmf[0], vb_cmf[0], cm_vel[0]);
	//	parafprintf(threebbfile, "new binary properties for object 'knew':  r=%g  vr=%g  vt=%g\n", star[knew].r, star[knew].vr, star[knew].vt);

		// DEBUG extra output
		//parafprintf(threebbdebugfile, "Phi_f=%g  E_f-phi_f=%g\n", star[knew].phi + star[k3].phi, Eb + 0.5*madhoc*(m3*(sqr(star[k3].vr) + sqr(star[k3].vt)) + star[knew].m*(sqr(star[knew].vr) + sqr(star[knew].vt))));

		//parafprintf(threebbdebugfile, "ms=%g  mb=%g  vs_r=%g  vs_t=%g  vb_r=%g  vb_t=%g\n\n\n", m3, star[knew].m, star[k3].vr, star[k3].vt, star[knew].vr, star[knew].vt);
		//parafprintf(threebbdebugfile, "binary[star[knew].binind].inuse=%ld\n", binary[star[knew].binind].inuse);
		//parafprintf(threebbfile, "%g %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", TotalTime, k1, k2, k3, m1, m2, m3, eta, Eb, binary[star[knew].binind].e, binary[star[knew].binind].a * units.l / AU, r_p * units.l / AU, star[knew].r, star[knew].vr, star[knew].vt, star[knew].phi, r3, star[k3].vr, star[k3].vt, star[k3].phi, delta_PE, delta_E, delta_E_3bb);


		}

}


/**
* @brief ?
*
* @param k ?
* @param p ?
* @param N_LIMIT ?
* @param ave_local_mass ?
* @param sigma_local ?
*/
void calc_sigma_local(long k, long p, long N_LIMIT, double *ave_local_mass, double *sigma_local)
// Modified by Meagan 2/10/12. Needed to calculate average local mass for 3bb formation, as well as sigma
// Now function passes both values back as global variables: ave_local_mass, sigma_local
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
}


// end functions created for three-body binary formation


/**
* @brief generate unique star id's (original version). Returns successive ids after N.
*
* @return unique id for newly created star
*/
long star_get_id_new(void)
{
	newstarid++;
	return(newstarid);
}

/**
* @brief generate unique star id's. The serial version of this function star_get_id_new generated new ids incrementally starting with N. But, we cant do this in the parallel version since this would need additional communication to keep the latest assigned id synchronized across processors, and also will result in a race condition. So, we use the following formula new_id = N + min(id1%N, id2%N) + id1/N + id2/N. Here / implies integer division, and % is the modulo.
*
* @param id1 id of first star
* @param id2 id of second star
*
* @return new id
*/
long star_get_merger_id_new(long id1, long id2)
{
	return (MIN(id1%newstarid, id2%newstarid) + id1/newstarid + id2/newstarid + newstarid);
}

/**
* @brief calculate local density by averaging
*
* @param k star index
* @param p number of stars to average over
* @param N_LIMIT typically the total number of stars
*
* @return local density for the given star
*/
double calc_n_local(long k, long p, long N_LIMIT)
{
	long kmin, kmax;

	kmin = k - p;
	//Bharath: Uneven no.of stars on either sides? Setting kmax = k+p rather than k+p+1.
	kmax = k + p;

	/* shouldn't the boundary here be kmin=1, not 0? */
	if (kmin < 0) {
		kmin = 0;
		kmax = 2 * p + 1;
	} else if (kmax > N_LIMIT) {
		kmax = N_LIMIT;
		kmin = N_LIMIT - 2 * p - 1;
	}

	return((2.0 * ((double) p)) * 3.0 / (4.0 * PI * (cub(star_r[kmax]) - cub(star_r[kmin]))));
}

/**
* @brief ?
*
* @param k index of star 1
* @param kp index of star 2
* @param p ?
* @param W ?
* @param N_LIMIT ?
*
* @return ?
*/
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

	return(3.0 * ((double) p) * sqr(star_m[k] + star_m[kp]) / (cub(W) * (cub(star_r[kmax]) - cub(star_r[kmin]))));
}

/**
* @brief ?
*
* @param k index of star 1
* @param kp index of star 2
* @param v[4] ?
* @param vp[4] ?
* @param w[4] ?
* @param W ?
* @param rcm ?
* @param vcm[4] ?
* @param rng ?
* @param setY ?
*/
void calc_encounter_dyns(long k, long kp, double v[4], double vp[4], double w[4], double *W, double *rcm, double vcm[4], gsl_rng *rng, int setY)
{
	int j;
	double phi;

	/* set random angle between vt's */
	/* first store random variable */
	if (setY) {
		star[k].Y = rng_t113_dbl_new(curr_st);
		//star[k].Y = rng_t113_dbl();

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
		// If you're seeing this and it's during the first timestep, most likely the initial conditions are wrong 
		// (e.g. if you've created a binary that is destroyed during stellar evolution initialization)
		eprintf("W = 0! for star index = %ld\tv1k = %g\tv2k = %g\tv1kp = %g\tv2kp = %g\nIf you're seeing this and it's during the first timestep, most likely the initial conditions are wrong\n (e.g. if you've created a binary that is destroyed during stellar evolution initialization)", star[k].id, star[k].vr, star[k].vt, star[kp].vr, star[kp].vt);
		exit_cleanly(1, __FUNCTION__);
	}

	int g_k = get_global_idx(k);
	int g_kp = get_global_idx(kp);
	*rcm = (star_m[g_k] * star_r[g_k] + star_m[g_kp] * star_r[g_kp]) / (star_m[g_k] + star_m[g_kp]);
	for (j=1; j<=3; j++) {
		vcm[j] = (star_m[g_k] * v[j] + star_m[g_kp] * vp[j]) / (star_m[g_k] + star_m[g_kp]);
	}
}

/**
* @brief Sets E and J values of the star based on vr and vt values, current position and potential.
*
* @param k index of star
*/
void set_star_EJ(long k)
{
	star[k].E = star_phi[get_global_idx(k)] + 0.5 * (sqr(star[k].vr) + sqr(star[k].vt));
	star[k].J = star_r[get_global_idx(k)] * star[k].vt;
}

/**
* @brief copies some variables of the star to the new variables
*
* @param k index of star
*/
void set_star_news(long k)
{
	star[k].rnew = star_r[get_global_idx(k)];
	star[k].vrnew = star[k].vr;
	star[k].vtnew = star[k].vt;
}

/**
* @brief copies some variables of the star to the old variables
*
* @param k index of star
*/
void set_star_olds(long k)
{
	int g_k = get_global_idx(k);
	star[k].rOld = star_r[g_k];
	star[k].r_peri = star_r[g_k];
	star[k].r_apo = star_r[g_k];
}

/**
* @brief copies the global/duplicated array values to the local star structure
*
* @param k index of star
*/
void copy_globals_to_locals(long k)
{
    int g_k = get_global_idx(k);
    star[k].m = star_m[g_k];
    star[k].r = star_r[g_k];
    star[k].phi = star_phi[g_k];
}

/**
* @brief find masses of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return masses of merging stars
*/
double binint_get_mass(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star_m[get_global_idx(k)]);
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
			return(star_m[get_global_idx(kp)]);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].m1);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].m2);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**Elena
* @brief find Radii of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return masses of merging stars
*/
double binint_get_radii(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].rad);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].rad1);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].rad2);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].rad);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].rad1);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].rad2);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**Elena
* @brief find core mass of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return core  masses of merging stars
*/
double binint_get_core_mass(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_mc);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_massc[0]);
			
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_massc[1]);
			
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_mc);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_massc[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_massc[1]);
			
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**Elena
* @brief find environment mass of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return env masses of merging stars
*/
double binint_get_env_mass(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_menv);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_menv[0]);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_menv[1]);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_menv);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_menv[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_menv[1]);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**Elena
* @brief find core Radii of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return core radii of merging stars
*/
double binint_get_core_radii(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_rc);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_radc[0]);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_radc[1]);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_rc);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_radc[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_radc[1]);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**Elena
* @brief find core Radii of merging stars from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return core radii of merging stars
*/
double binint_get_env_radii(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_renv);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_renv[0]);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_renv[1]);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_renv);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_renv[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_renv[1]);
		}
	}
	
	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**
* @brief find spins of merging black holes from binary interaction components
*
* @param k index of first star
* @param kp index of second star
* @param id id of star
*
* @return spin of merging black holes
*/
double binint_get_spins(long k, long kp, long id)
{
	/* first look at k */
	if (star[k].binind == 0) {
		if (star[k].id == id) {
			return(star[k].se_bhspin);
		}
	} else {
		if (binary[star[k].binind].id1 == id) {
			return(binary[star[k].binind].bse_bhspin[0]);
		} else if (binary[star[k].binind].id2 == id) {
			return(binary[star[k].binind].bse_bhspin[1]);
		}
	}
	
	/* then at kp */
	if (star[kp].binind == 0) {
		if (star[kp].id == id) {
			return(star[kp].se_bhspin);
		}
	} else {
		if (binary[star[kp].binind].id1 == id) {
			return(binary[star[kp].binind].bse_bhspin[0]);
		} else if (binary[star[kp].binind].id2 == id) {
			return(binary[star[kp].binind].bse_bhspin[1]);
		}
	}

	eprintf("cannot find matching id %ld!\n", id);
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**
* @brief Sourav: find stellar types of merging stars from binary interaction components
*
* @param k index of star 1
* @param kp index of star 2
* @param id star id to be matched
*
* @return ?
*/
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
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}


/**
* @brief Sourav: Sourav: toy rejuvenation- finding creation times of binary interaction components
*
* @param k index of star 1
* @param kp index of star 2
* @param id star id to be matched
*
* @return ?
*/
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
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**
* @brief Sourav: Sourav: toy rejuvenation- finding MS lifetimes of binary interaction components
*
* @param k index of star 1
* @param kp index of star 2
* @param id star id to be matched
*
* @return ?
*/
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
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}


/**
* @brief return star index of star with id "id" from binary interaction components, along with which member
*
* @param k index of star 1
* @param kp index of star 2
* @param id id of star to be matched
* @param bi "bi=0,1" if a binary, else single
*
* @return return star index of star with id "id" from binary interaction components, along with which member
*/
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
	exit_cleanly(1, __FUNCTION__);
	/* this is just for the compiler */
	exit(1);
}

/**
* @brief ?
*
* @param obj ?
* @param units ?
*/
void binint_log_obj(fb_obj_t *obj, fb_units_t units)
{
	int bid, sid, i;
	char dumstring[FB_MAX_STRING_LENGTH], idstring1[FB_MAX_STRING_LENGTH], idstring2[FB_MAX_STRING_LENGTH], idstring3[FB_MAX_STRING_LENGTH];
	double mtriplein, mtriple, pouttriple2, pintriple, tlkquad, epsoct, tlkoct, tgr, epsgr;
	double l0[3], l1[3], LL1[3], LL2[3], i12;
	
	if (fb_n_hier(obj) == 1) {
		/* first write id string */
		snprintf(idstring1, FB_MAX_STRING_LENGTH, "%ld", obj->id[0]);
		for (i=1; i<obj->ncoll; i++) {
			snprintf(dumstring, FB_MAX_STRING_LENGTH, ":%ld", obj->id[i]);
			strncat(idstring1, dumstring, FB_MAX_STRING_LENGTH);
		}
		/* then print to log */
		parafprintf(binintfile, "type=single m=%g R=%g Eint=%g id=%s ktype=%d\n", obj->m*units.m/MSUN, obj->R*units.l/RSUN, obj->Eint*units.E, idstring1, obj->k_type);
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
		parafprintf(binintfile, "type=binary m0=%g m1=%g R0=%g R1=%g Eint1=%g Eint2=%g id0=%s id1=%s a=%g e=%g ktype1=%d ktype2=%d\n", 
			obj->obj[0]->m*units.m/MSUN, obj->obj[1]->m*units.m/MSUN, 
			obj->obj[0]->R*units.l/RSUN, obj->obj[1]->R*units.l/RSUN, 
			obj->obj[0]->Eint*units.E, obj->obj[1]->Eint*units.E, 
			idstring1, idstring2, 
			obj->a*units.l/AU, obj->e,
			obj->obj[0]->k_type, obj->obj[1]->k_type);
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
		parafprintf(binintfile, "type=triple min0=%g min1=%g mout=%g Rin0=%g Rin1=%g Rout=%g Eintin0=%g Eintin1=%g Eintout=%g idin1=%s idin2=%s idout=%s ain=%g aout=%g ein=%g eout=%g ktypein1=%d ktypein1=%d ktypeout=%d\n",
			obj->obj[bid]->obj[0]->m*units.m/MSUN, obj->obj[bid]->obj[1]->m*units.m/MSUN, obj->obj[sid]->m*units.m/MSUN,
			obj->obj[bid]->obj[0]->R*units.l/RSUN, obj->obj[bid]->obj[1]->R*units.l/RSUN, obj->obj[sid]->R*units.l/RSUN,
			obj->obj[bid]->obj[0]->Eint*units.E, obj->obj[bid]->obj[1]->Eint*units.E, obj->obj[sid]->Eint*units.E,
			idstring1, idstring2, idstring3, 
			obj->obj[bid]->a*units.l/AU, obj->a*units.l/AU,
			obj->obj[bid]->e, obj->e,
			obj->obj[bid]->obj[0]->k_type, obj->obj[bid]->obj[1]->k_type, obj->obj[sid]->k_type);

		/* In addition, for triples we may want the full configuration of the system*/

		/* Fist compute the inclincation of the two binaries*/
		fb_cross(obj->obj[bid]->obj[0]->x, obj->obj[bid]->obj[0]->v, l0);
                fb_cross(obj->obj[bid]->obj[1]->x, obj->obj[bid]->obj[1]->v, l1);
                for (i=0; i<3; i++) {
                	LL1[i] = obj->obj[bid]->obj[0]->m * l0[i] + obj->obj[bid]->obj[1]->m * l1[i];
                }
                fb_cross(obj->obj[0]->x, obj->obj[0]->v, l0);
                fb_cross(obj->obj[1]->x, obj->obj[1]->v, l1);
                for (i=0; i<3; i++) {
                	LL2[i] = obj->obj[0]->m * l0[i] + obj->obj[1]->m * l1[i];
                }
                i12 = 180.0/3.14159*acos(fb_dot(LL1, LL2)/(fb_mod(LL1)*fb_mod(LL2)));

		/* And some of the relevant timescales */
		mtriplein=obj->obj[bid]->obj[0]->m*units.m/MSUN+obj->obj[bid]->obj[1]->m*units.m/MSUN;
		mtriple=mtriplein+obj->obj[sid]->m*units.m/MSUN;
		pouttriple2=pow(obj->a*units.l/AU,3)/mtriple;
		pintriple=pow(obj->obj[bid]->a*units.l/AU,3)/mtriplein;
                tlkquad=(8.0/3.1415/15.0)*mtriple/(obj->obj[sid]->m*units.m/MSUN)*pouttriple2/pintriple*pow(1.0-pow(obj->e,2),1.5);
		epsoct=abs(obj->obj[bid]->obj[0]->m*units.m/MSUN-obj->obj[bid]->obj[1]->m*units.m/MSUN)/mtriplein*(obj->obj[bid]->a*units.l/AU)/(obj->a*units.l/AU)*obj->e/(1.0-pow(obj->e,2));
		tlkoct=tlkquad/epsoct;
		tgr=pow(10.0,8)*pow(obj->obj[bid]->a*units.l/AU,2.5)/(3.0*pow(mtriplein,1.5))*(1.0-pow(obj->obj[bid]->e,2));
		epsgr=tlkquad/tgr;

                parafprintf(triplefile, "%.18g %g %g %g %g %g %g %g %g %g %g %g %d %d %d %g %g %g %g %g\n",
                        TotalTime, obj->obj[bid]->obj[0]->m*units.m/MSUN, obj->obj[bid]->obj[1]->m*units.m/MSUN, obj->obj[sid]->m*units.m/MSUN,
                        obj->obj[bid]->obj[0]->R*units.l/RSUN, obj->obj[bid]->obj[1]->R*units.l/RSUN, obj->obj[sid]->R*units.l/RSUN,
                        obj->obj[bid]->a*units.l/AU, obj->a*units.l/AU,
                        obj->obj[bid]->e, obj->e, i12,
                        obj->obj[bid]->obj[0]->k_type, obj->obj[bid]->obj[1]->k_type, obj->obj[sid]->k_type,
			tlkquad,tlkoct,epsoct,tgr,epsgr);
	} else {
		/* thankfully won't need to print out quads */
		eprintf("Don't know how to print out object with >3 stars!\n");
		exit_cleanly(1, __FUNCTION__);
	}
}

/**
* @brief ?
*
* @param retval ?
*/
void binint_log_status(fb_ret_t retval, double vesc)
{
	/* must print out Nosc when upgraded to latest Fewbody */
	parafprintf(binintfile, "status: DE/E=%g DE=%g DL/L=%g DL=%g DE_GW/E=%g DE_GW=%g v_esc_cluster[km/s]=%g tcpu=%g\n", 
		retval.DeltaEfrac, retval.DeltaE, retval.DeltaLfrac, retval.DeltaL, retval.DeltaE_GWfrac, retval.DeltaE_GW, vesc, retval.tcpu);
}

/**
* @brief ?
*
* @param interaction_type[] ?
* @param id ?
* @param mass ?
* @param r ?
* @param obj ?
* @param k index of star 1
* @param kp index of star 2
* @param startype star type
*/
void binint_log_collision(const char interaction_type[], long id,
			  double mass, double r, fb_obj_t obj, long k, long kp, long startype)
{
	int j;

	parafprintf(collisionfile, "t=%g %s idm=%ld(mm=%g) id1=%ld(m1=%g)",
		TotalTime, interaction_type, id, 
		mass * units.mstar / FB_CONST_MSUN, obj.id[0], 
		binint_get_mass(k, kp, obj.id[0]) * units.mstar 
						  / FB_CONST_MSUN);
	for (j=1; j<obj.ncoll; j++) {
		parafprintf(collisionfile, ":id%d=%ld(m%d=%g)", 
			j+1, obj.id[j], j+1,
			binint_get_mass(k, kp, obj.id[j]) * units.mstar 
							  / FB_CONST_MSUN);
	}
	parafprintf(collisionfile," (r=%g) ", r);
//Sourav
	parafprintf(collisionfile, "typem=%ld ", startype);
	for (j=0; j<obj.ncoll; j++) {
		parafprintf(collisionfile, "type%d=%ld ", j+1, 
				binint_get_startype(k, kp, obj.id[j]));// Use this, not the Fewbody type, since this is changed by BSE after mergers
	}
//Elena: extra output for bs and bb interactions
	
	for (j=0; j<obj.ncoll; j++) {
	        parafprintf(collisionfile, "rad%d[RSUN]=%g ", j+1, binint_get_radii(k, kp, obj.id[j])*units.l/RSUN);
	}

	parafprintf(collisionfile, "\n");
}

/**
* @brief Store additional information in morecollisions file
*
* @param interaction_type[] ?
* @param id ?
* @param mass ?
* @param r ?
* @param obj ?
* @param k index of star 1
* @param kp index of star 2
* @param startype star type
*/
void binint_log_morecollision(const char interaction_type[], long remnant_id,
			  double remnant_mass, double remnant_radius, long remnant_type, double remnant_mc, double remnant_menv, double remnant_rc, 				double remnant_renv, fb_obj_t obj, long k, long kp, double W, double rperi)
{

	/*remnant radii are already in the right units*/
	
	double m0_core = binint_get_core_mass(k, kp, obj.id[0]);
	double m1_core = binint_get_core_mass(k, kp, obj.id[1]);
	double m0_env = binint_get_env_mass(k, kp, obj.id[0]);
	double m1_env = binint_get_env_mass(k, kp, obj.id[1]);
	
	double r0_core = binint_get_core_radii(k, kp, obj.id[0]);
	double r1_core = binint_get_core_radii(k, kp, obj.id[1]);
	double r0_env = binint_get_env_radii(k, kp, obj.id[0]);
	double r1_env = binint_get_env_radii(k, kp, obj.id[1]);
	
	
	double rho0_c   = (m0_core) / ((4/3)* PI * pow((r0_core),3));
	double rho1_c   = (m1_core) / ((4/3)* PI * pow((r1_core),3));
	double rho0_env = (m0_env)  / ((4/3)* PI * pow((r0_env),3));
	double rho1_env = (m1_env)  / ((4/3)* PI * pow((r1_env),3));
	double rhor_c   = (remnant_mc)   / ((4/3)* PI * pow((remnant_rc ),3));
	double rhor_env = (remnant_menv) / ((4/3)* PI * pow((remnant_renv),3));

	if(isnan(rho0_c)){rho0_c = -100;}
	if(isnan(rho1_c)){rho1_c = -100;}
	if(isnan(rhor_c)){rhor_c = -100;}
	if(isnan(rho0_env)){rho0_env = -100;}
	if(isnan(rho1_env)){rho1_env = -100;}
	if(isnan(rhor_env)){rhor_env = -100;}
	
	// Elena: For some stars, COSMIC assigns default renv and menv values of of e-10, which makes my densities exactly 3.1831e-19. I will change these vales to output a -100 intead, since it is not physical //
	

	if(rho0_env >= 1.0e19){rho0_env = -100;}
	if(rho1_env >= 1.0e19){rho1_env = -100;}
	if(rhor_env >= 1.0e19){rhor_env = -100;}

	parafprintf(morecollfile, "%g %s %ld %ld %g %g %g %g %g %g %g %g %d %d %ld %g %g %g %g %d %g %g\n",
				    TotalTime, interaction_type, obj.id[0], obj.id[1], 
				    binint_get_mass(k, kp, obj.id[0]) * units.mstar / FB_CONST_MSUN, 
				    binint_get_mass(k, kp, obj.id[1]) * units.mstar / FB_CONST_MSUN,
				    binint_get_radii(k, kp, obj.id[0]) * units.l/RSUN, 
				    binint_get_radii(k, kp, obj.id[1]) * units.l/RSUN,
				    rho0_c,rho1_c,rho0_env, rho1_env, 
				    binint_get_startype(k,kp, obj.id[0]), binint_get_startype(k,kp, obj.id[1]), 
				    remnant_id, remnant_mass * units.mstar / FB_CONST_MSUN, remnant_radius*units.l/RSUN, 
				    rhor_c, rhor_env,remnant_type, W*units.l/units.t/1.e5, rperi*units.l/RSUN);
				    
}

/**
* @brief do binary interaction (bin-bin or bin-single)
*
* @param k index of 1st star
* @param kp index of 2nd star
* @param rperi ?
* @param w[4] ?
* @param W ?
* @param rcm ?
* @param vcm[4] ?
* @param rng gsl rng
*/
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
	double vs[20], VK0;
	double energy_from_outer=0.;

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
		
		cmc_units.v = sqrt((star_m[get_global_idx(k)]+star_m[get_global_idx(kp)])/(star_m[get_global_idx(k)]*star_m[get_global_idx(kp)]) * 
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

		cmc_units.v = sqrt((star_m[get_global_idx(ksin)]+star_m[get_global_idx(kbin)])/(star_m[get_global_idx(ksin)]*star_m[get_global_idx(kbin)]) * 
				   (binary[jbin].m1 * binary[jbin].m2 / binary[jbin].a) * madhoc);
		cmc_units.l = binary[jbin].a;
		cmc_units.t = cmc_units.l / cmc_units.v;
		cmc_units.m = cmc_units.l * sqr(cmc_units.v);
		cmc_units.E = cmc_units.m * sqr(cmc_units.v);
	} else {
		eprintf("no binaries!");
		exit_cleanly(1, __FUNCTION__);
		exit(1);
	}
	
	/* malloc hier (based on value of hier.nstarinit) */
	fb_malloc_hier(&hier);

	t=0;
	bmax = rperi * sqrt(1.0 + 2.0 * ((star_m[get_global_idx(k)] + star_m[get_global_idx(kp)]) * madhoc) / (rperi * sqr(W)));

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
		exit_cleanly(1, __FUNCTION__);
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

        /* record the escape speed of the cluser where the encounter occured*/
        int g_k = get_global_idx(k);
        double vesc;
        vesc = sqrt(-2*star_phi[g_k]) * (units.l/units.t) / 1.0e5;

	/* logging */
	binint_log_status(retval,vesc);
	printing_units.v = cmc_units.v * units.l / units.t;
	printing_units.l = cmc_units.l * units.l;
	printing_units.t = cmc_units.t * units.t;
	printing_units.m = cmc_units.m * units.m;
	printing_units.E = cmc_units.E * units.E;

	/* now do something with the Fewbody result */
	if ( !( (fabs(retval.DeltaEfrac) < 1.0e-3 || fabs(retval.DeltaE) < 1.0e-3) && 
		 (fabs(retval.DeltaLfrac) < 1.0e-3 || fabs(retval.DeltaL) < 1.0e-3) ) && 
         (!((fabs(retval.DeltaE_GWfrac) > 1.0e-3 || fabs(retval.DeltaE_GW > 1.0e-3)
            ) && retval.PN_ON == 1))) /* did we have a significant energy error that wasn't from gravitational waves? */
    {
		parafprintf(binintfile, "outcome: energy and/or angular momentum error\n");
		print_interaction_error();
	} else if ( isnan(retval.DeltaE) || isnan(retval.DeltaL) ) {
		parafprintf(binintfile, "outcome: NaN returned by fewbody\n");
		print_interaction_error();
	} else if (retval.retval == 0) {
		/* bad outcome; ignore for now */
		parafprintf(binintfile, "outcome: stopped\n");
		print_interaction_error();
	} else if (hier.obj[0]->n == 4) {
		/* outcome is a quadruple */
		parafprintf(binintfile, "outcome: error\n");
		print_interaction_error();
	} else {
		parafprintf(binintfile, "outcome: %s (%s)\n", fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2));
		
		for (i=0; i<hier.nobj; i++) {
			/* logging */
			parafprintf(binintfile, "output: ");
			binint_log_obj(hier.obj[i], printing_units);

			/* single/binary/triple stars */
			istriple = 0;
			if (hier.obj[i]->n == 1) {
				knew = create_star(k, 0);
			} else if (hier.obj[i]->n == 2) {
				knew = create_binary(k, 0);
			} else if (hier.obj[i]->n == 3) {
				istriple = 1;
				/* break triple for now */
				knew = create_binary(k, 0);
				knewp = create_star(k, 0);
				
				if (hier.obj[i]->obj[0]->n == 1) {
					sid = 0;
					bid = 1;
				} else {
					sid = 1;
					bid = 0;
				}
			} else {
				eprintf("object with n=%d!\n", hier.obj[i]->n);
				exit_cleanly(1, __FUNCTION__);
				/* this is just for the compiler */
				exit(1);
			}
			
			/* generic properties */
			/* set radial position */


			star_r[get_global_idx(knew)] = rcm;
			if (istriple) {
				star_r[get_global_idx(knewp)] = rcm;
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
				star_m[get_global_idx(knew)] = hier.obj[i]->obj[bid]->m * cmc_units.m / madhoc;
				star_m[get_global_idx(knewp)] = hier.obj[i]->obj[sid]->m * cmc_units.m / madhoc;
			} else {
				star_m[get_global_idx(knew)] = hier.obj[i]->m * cmc_units.m / madhoc;
			}

			/* set potential */
			star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);
			if (istriple) {
				star_phi[get_global_idx(knewp)] = potential(star_r[get_global_idx(knewp)]);
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
					star[knew].id = hier.obj[i]->id[0];
					nmerged = 1;

					while (nmerged < hier.obj[i]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->id[nmerged], &bi);
						star[knew].id = star_get_merger_id_new(star[knew].id, hier.obj[i]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar);
						cp_m_to_star(oldk, bi, &tempstar);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(tempstar.se_k == 14 || star[knew].se_k == 14 ){
                            star[knew].m = hier.obj[i]->m * cmc_units.m / madhoc;
                            star[knew].se_mass = star[knew].m * units.mstar / MSUN;
                            star[knew].se_mt  = star[knew].m * units.mstar / MSUN;
                            star[knew].se_mc = star[knew].m * units.mstar / MSUN;
                            star[knew].se_bhspin = hier.obj[i]->chi;
                            star[knew].se_radius = hier.obj[i]->R * cmc_units.l / BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            star[knew].Eint = 0;
                            if(WRITE_BH_INFO && tempstar.se_k == 14 && star[knew].se_k == 14)
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->id[0],hier.obj[i]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
					    binint_get_spins(k, kp, hier.obj[i]->id[0]), binint_get_spins(k, kp, hier.obj[i]->id[nmerged]), star[knew].id, 
                                                          star[knew].m*units.mstar/MSUN,hier.obj[i]->chi,hier.obj[i]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->e_merger[nmerged],
														  hier.obj[i]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->e_50M[nmerged],
														  hier.obj[i]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->e_100M[nmerged],
														  hier.obj[i]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->e_500M[nmerged]);
                            star[knew].se_k = 14;
                        } else{
                            merge_two_stars(&(star[knew]), &tempstar, &(star[knew]), vs, curr_st);
                                                    /* Owing to merger only useful vs's are v[1-3] */
                            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);

                            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                            //star[knew].vt += sqrt(vs[1]*vs[1]+vs[2]*vs[2]) * 1.0e5 / (units.l/units.t);
                        }
						nmerged++;
					}
					set_star_EJ(knew);
					
					//star[knew].id = star_get_id_new();

					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge1: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,star[knew].id);
					}
					
					/* log collision */
                    star_m[get_global_idx(knew)] = star[knew].m;
					binint_log_collision(isbinbin?"binary-binary":"binary-single",
						star[knew].id, star_m[get_global_idx(knew)],
						star_r[get_global_idx(knew)],
						*(hier.obj[i]), k, kp, star[knew].se_k);

					if (WRITE_MORECOLL_INFO  && hier.obj[i]-> ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", star[knew].id,
			  			star_m[get_global_idx(knew)], star[knew].rad, star[knew].se_k, star[knew].se_mc, 
			  			star[knew].se_menv, star[knew].se_rc,star[knew].se_renv,*(hier.obj[i]),k,kp, W, rperi);
			  			}

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
					tempstar.id = hier.obj[i]->obj[0]->id[0];
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[0]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[0]->id[nmerged], &bi);
						tempstar.id = star_get_merger_id_new(tempstar.id, hier.obj[i]->obj[0]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(tempstar2.se_k == 14 || tempstar.se_k == 14){
                            tempstar.m = hier.obj[i]->obj[0]->m * cmc_units.m / madhoc;
                            tempstar.se_mass = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mt  = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mc = tempstar.m * units.mstar / MSUN;
                            tempstar.se_bhspin = hier.obj[i]->obj[0]->chi;
                            tempstar.se_radius = hier.obj[i]->obj[0]->R * cmc_units.l/ BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            tempstar.Eint = 0;
                            if(WRITE_BH_INFO && tempstar2.se_k == 14 && tempstar.se_k == 14)
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->obj[0]->id[0],hier.obj[i]->obj[0]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[0]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[0]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
                                                          binint_get_spins(k, kp, hier.obj[i]->obj[0]->id[0]), binint_get_spins(k, kp, hier.obj[i]->obj[0]->id[nmerged]), 
					    tempstar.id, tempstar.m*units.mstar/MSUN,hier.obj[i]->obj[0]->chi,hier.obj[i]->obj[0]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->obj[0]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[0]->e_merger[nmerged],
														  hier.obj[i]->obj[0]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[0]->e_50M[nmerged],
														  hier.obj[i]->obj[0]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[0]->e_100M[nmerged],
														  hier.obj[i]->obj[0]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[0]->e_500M[nmerged]);
                            tempstar.se_k = 14;
                        } else{
                            merge_two_stars(&tempstar, &tempstar2, &tempstar, vs, curr_st);
                            /* FIXME: really we're supposed to add the kick to each binary
                               member separately, then calculate the systemic kick to the binary,
                               but hopefully this doesn't happen too much. */
                                                    /* The kick routine within /bse_wrap/bse/ correctly updates COM velocity... */
                            if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
                                wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
                                    sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
                            }
                            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
                            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        }
						nmerged++;
					}
					set_star_EJ(knew);

					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 0);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 0);

					//binary[star[knew].binind].id1 = star_get_id_new();
					binary[star[knew].binind].id1 = tempstar.id;

					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge2: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id1);
					}
					
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single",
						binary[star[knew].binind].id1,
						binary[star[knew].binind].m1,
						star_r[get_global_idx(knew)],
						*(hier.obj[i]->obj[0]), k, kp, binary[star[knew].binind].bse_kw[0]);
					if (binary[star[knew].binind].m1==0.) {
						dprintf("Zero mass remnant! Parameters: knew=%li, binind=%li, kw[0]=%i, kw[1]=%i\n",
								knew, star[knew].binind, binary[star[knew].binind].bse_kw[0], 
								binary[star[knew].binind].bse_kw[1]);
					}
					

					if (WRITE_MORECOLL_INFO  && hier.obj[i]->obj[0]->ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", binary[star[knew].binind].id1,
			  			binary[star[knew].binind].m1, binary[star[knew].binind].rad1, binary[star[knew].binind].bse_kw[0], 							binary[star[knew].binind].bse_massc[0], binary[star[knew].binind].bse_menv[0], 							binary[star[knew].binind].bse_radc[0],binary[star[knew].binind].bse_renv[0],
			  			*(hier.obj[i]->obj[0]),k,kp, W, rperi);}

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
					tempstar.id = hier.obj[i]->obj[1]->id[0];
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[1]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[1]->id[nmerged], &bi);
						tempstar.id = star_get_merger_id_new(tempstar.id, hier.obj[i]->obj[1]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(tempstar2.se_k == 14 || tempstar.se_k == 14){
                            tempstar.m = hier.obj[i]->obj[1]->m * cmc_units.m / madhoc;
                            tempstar.se_mass = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mt  = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mc = tempstar.m * units.mstar / MSUN;
                            tempstar.se_bhspin = hier.obj[i]->obj[1]->chi;
                            tempstar.se_radius = hier.obj[i]->obj[1]->R * cmc_units.l/ BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            tempstar.Eint = 0;
                            if(WRITE_BH_INFO && tempstar2.se_k == 14 && tempstar.se_k == 14)
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->obj[1]->id[0],hier.obj[i]->obj[1]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[1]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[1]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
                                                          binint_get_spins(k, kp, hier.obj[i]->obj[1]->id[0]), binint_get_spins(k, kp, hier.obj[i]->obj[1]->id[nmerged]), 
					    tempstar.id, tempstar.m*units.mstar/MSUN,hier.obj[i]->obj[1]->chi,hier.obj[i]->obj[1]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->obj[1]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[1]->e_merger[nmerged],
														  hier.obj[i]->obj[1]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[1]->e_50M[nmerged],
														  hier.obj[i]->obj[1]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[1]->e_100M[nmerged],
														  hier.obj[i]->obj[1]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[1]->e_500M[nmerged]);
                            tempstar.se_k = 14;
                        } else{
                            merge_two_stars(&tempstar, &tempstar2, &tempstar, vs, curr_st);
                            /* FIXME: really we're supposed to add the kick to each binary
                               member separately, then calculate the systemic kick to the binary,
                               but hopefully this doesn't happen too much. */
                            if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
                                wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
                                    sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
                            }
                            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
                            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        }
						nmerged++;
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 1);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 1);

					//binary[star[knew].binind].id2 = star_get_id_new();
					binary[star[knew].binind].id2 = tempstar.id;

					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge3: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id2);
					}

					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single",
						binary[star[knew].binind].id2,
						binary[star[knew].binind].m2,
						star_r[get_global_idx(knew)],
						*(hier.obj[i]->obj[1]), k, kp, binary[star[knew].binind].bse_kw[1]);

					if (binary[star[knew].binind].m2==0.)
						dprintf("Zero mass remnant! Parameters: knew=%li binind=%li kw[0]=%i kw[1]=%i\n",
								knew, star[knew].binind, binary[star[knew].binind].bse_kw[0],
								binary[star[knew].binind].bse_kw[1]);

					if (WRITE_MORECOLL_INFO  && hier.obj[i]->obj[1]->ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", binary[star[knew].binind].id2,
			  			binary[star[knew].binind].m2, binary[star[knew].binind].rad2, binary[star[knew].binind].bse_kw[1], 							binary[star[knew].binind].bse_massc[1], binary[star[knew].binind].bse_menv[1], 
			  			binary[star[knew].binind].bse_radc[1], binary[star[knew].binind].bse_renv[1],
			  			*(hier.obj[i]->obj[1]),k,kp, W, rperi);}
				}
				
				star_m[get_global_idx(knew)] = binary[star[knew].binind].m1 + binary[star[knew].binind].m2;

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
				
				/* determine the difference in energy between the inner binary and the 
 				   triple as a whole; this will be added to the cluster in break_wide_binaries */
				energy_from_outer = (-(hier.obj[i]->obj[bid]->obj[0]->m)*(hier.obj[i]->obj[bid]->obj[1]->m)/
							(2.0 * hier.obj[i]->obj[bid]->a)) - (fb_ketot(threeobjs, 3) + fb_petot(threeobjs, 3));

				/* Unless fewbody has screwed up the classification, this can't happen, but better safe than sorry */
				if (energy_from_outer < 0) {
					eprintf("energy_from_outer is negative; fewbody classification has clearly screwed up\n");
					energy_from_outer = 0.;
				}
				
				/********************************/
				/* set single star's properties */
				/********************************/
				/* internal energy */
				star[knewp].Eint = hier.obj[i]->obj[sid]->Eint * cmc_units.E;
				star[knewp].E_excess = energy_from_outer * cmc_units.E;

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
					star[knew].id = hier.obj[i]->obj[sid]->id[0];
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[sid]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[sid]->id[nmerged], &bi);
						star[knew].id = star_get_merger_id_new(star[knew].id, hier.obj[i]->obj[sid]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar);
						cp_m_to_star(oldk, bi, &tempstar);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(star[knewp].se_k == 14 || tempstar.se_k == 14){
                            star[knewp].m = hier.obj[i]->obj[sid]->m * cmc_units.m / madhoc;
                            star[knewp].se_mass = star[knewp].m * units.mstar / MSUN;
                            star[knewp].se_mt  = star[knewp].m * units.mstar / MSUN;
                            star[knewp].se_mc = star[knewp].m * units.mstar / MSUN;
                            star[knewp].se_bhspin = hier.obj[i]->obj[sid]->chi;
                            star[knewp].se_radius = hier.obj[i]->obj[sid]->R * cmc_units.l/ BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            star[knewp].Eint = 0;
                            if(WRITE_BH_INFO && star[knewp].se_k == 14 && tempstar.se_k == 14 )
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->obj[sid]->id[0],hier.obj[i]->obj[sid]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[sid]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[sid]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
                                                          binint_get_spins(k, kp, hier.obj[i]->obj[sid]->id[0]), binint_get_spins(k, kp, hier.obj[i]->obj[sid]->id[nmerged]), 
					    star[knewp].id, star[knewp].m*units.mstar/MSUN,hier.obj[i]->obj[sid]->chi,hier.obj[i]->obj[sid]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->obj[sid]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[sid]->e_merger[nmerged],
														  hier.obj[i]->obj[sid]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[sid]->e_50M[nmerged],
														  hier.obj[i]->obj[sid]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[sid]->e_100M[nmerged],
														  hier.obj[i]->obj[sid]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[sid]->e_500M[nmerged]);
                            star[knewp].se_k = 14;
                        } else{
                            merge_two_stars(&(star[knewp]), &tempstar, &(star[knewp]), vs, curr_st);
                            star[knewp].vr += vs[3] * 1.0e5 / (units.l/units.t);					       

                            //MPI2: parallel rng mimicking removed due to parent function which takes k as parameter?
                            vt_add_kick(&(star[knewp].vt),vs[1],vs[2], curr_st);
                        }
						nmerged++;
					}
					set_star_EJ(knewp);

					//star[knewp].id = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge4: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,star[knewp].id);
					}
					/* log collision */
                    star_m[get_global_idx(knewp)] = star[knewp].m;
					binint_log_collision(isbinbin?"binary-binary":"binary-single",
						star[knewp].id, star_m[get_global_idx(knewp)],
						star_r[get_global_idx(knewp)],
						*(hier.obj[i]->obj[sid]), k, kp, star[knewp].se_k);

					if (WRITE_MORECOLL_INFO  && hier.obj[i]->obj[sid]->ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", star[knewp].id,
			  			star_m[get_global_idx(knewp)], star[knewp].rad, star[knewp].se_k, star[knewp].se_mc, 
			  			star[knewp].se_menv, star[knewp].se_rc,star[knewp].se_renv,*(hier.obj[i]->obj[sid]),k,kp, W, rperi);}	
		
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
					tempstar.id = hier.obj[i]->obj[bid]->obj[0]->id[0];
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[bid]->obj[0]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[nmerged], &bi);
                        tempstar.id = star_get_merger_id_new(tempstar.id, hier.obj[i]->obj[bid]->obj[0]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(tempstar2.se_k == 14 || tempstar.se_k == 14){
                            tempstar.m = hier.obj[i]->obj[bid]->obj[0]->m * cmc_units.m / madhoc;
                            tempstar.se_mass = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mt  = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mc = tempstar.m * units.mstar / MSUN;
                            tempstar.se_bhspin = hier.obj[i]->obj[bid]->obj[0]->chi;
                            tempstar.se_radius = hier.obj[i]->obj[bid]->obj[0]->R * cmc_units.l/ BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            tempstar.Eint = 0;
                            if(WRITE_BH_INFO && tempstar2.se_k == 14 && tempstar.se_k == 14)
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->obj[bid]->obj[0]->id[0],hier.obj[i]->obj[bid]->obj[0]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
                                                          binint_get_spins(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[0]), binint_get_spins(k, kp, hier.obj[i]->obj[bid]->obj[0]->id[nmerged]), 
					    tempstar.id, tempstar.m*units.mstar/MSUN,hier.obj[i]->obj[bid]->obj[0]->chi,hier.obj[i]->obj[bid]->obj[0]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->obj[bid]->obj[0]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[0]->e_merger[nmerged],
														  hier.obj[i]->obj[bid]->obj[0]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[0]->e_50M[nmerged],
														  hier.obj[i]->obj[bid]->obj[0]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[0]->e_100M[nmerged],
														  hier.obj[i]->obj[bid]->obj[0]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[0]->e_500M[nmerged]);
                            tempstar.se_k = 14;
                        } else{
                            merge_two_stars(&tempstar, &tempstar2, &tempstar, vs, curr_st);
                            /* FIXME: really we're supposed to add the kick to each binary
                               member separately, then calculate the systemic kick to the binary,
                               but hopefully this doesn't happen too much. */
                            if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
                                wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
                                    sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
                            }
                            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
                            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        }
						nmerged++;
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 0);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 0);

					binary[star[knew].binind].id1 = tempstar.id;
					//binary[star[knew].binind].id1 = star_get_id_new();
					if(vs[1]!=0.0){
						VK0 = sqrt(sqr(vs[1])+sqr(vs[2])+sqr(vs[3]));
						dprintf("dynhelp_merge5: TT=%.18g vs[0]=%.18g vs[1]=%.18g vs[2]=%.18g vs[3]=%.18g vs[4]=%.18g vs[5]=%.18g vs[6]=%.18g VK0=%.18g star_id=%ld\n",TotalTime,vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],VK0,binary[star[knew].binind].id1);
					}
					/* log collision */
					
					binint_log_collision(isbinbin?"binary-binary":"binary-single",
						binary[star[knew].binind].id1,
						binary[star[knew].binind].m1,
						star_r[get_global_idx(knew)],
						*(hier.obj[i]->obj[bid]->obj[0]), k, kp, binary[star[knew].binind].bse_kw[0]);

					if (WRITE_MORECOLL_INFO  && hier.obj[i]->obj[bid]->obj[0]->ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", binary[star[knew].binind].id1,
			  			binary[star[knew].binind].m1, binary[star[knew].binind].rad1, 
			  			binary[star[knew].binind].bse_kw[0], 
			  			binary[star[knew].binind].bse_massc[0], binary[star[knew].binind].bse_menv[0], 
			  			binary[star[knew].binind].bse_radc[0], binary[star[knew].binind].bse_renv[0],
			  			*(hier.obj[i]->obj[bid]->obj[0]),k,kp, W, rperi);}
			  				
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
					tempstar.id = hier.obj[i]->obj[bid]->obj[1]->id[0];
					nmerged = 1;
					while (nmerged < hier.obj[i]->obj[bid]->obj[1]->ncoll) {
						oldk = binint_get_indices(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[nmerged], &bi);
                        tempstar.id = star_get_merger_id_new(tempstar.id, hier.obj[i]->obj[bid]->obj[1]->id[nmerged]);
						cp_SEvars_to_star(oldk, bi, &tempstar2);
						cp_m_to_star(oldk, bi, &tempstar2);
                        /* NOTE: if I have a BH/star or BBH merger, this will overwrite the 
                         * BSE results with the fewbody dynamical results */
                        if(tempstar2.se_k == 14 || tempstar.se_k == 14){
                            tempstar.m = hier.obj[i]->obj[bid]->obj[1]->m * cmc_units.m / madhoc;
                            tempstar.se_mass = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mt  = tempstar.m * units.mstar / MSUN;
                            tempstar.se_mc = tempstar.m * units.mstar / MSUN;
                            tempstar.se_bhspin = hier.obj[i]->obj[bid]->obj[1]->chi;
                            tempstar.se_radius = hier.obj[i]->obj[bid]->obj[1]->R * cmc_units.l/ BH_RADIUS_MULTIPLYER * units.l / RSUN;
                            tempstar.Eint = 0;
                            if(WRITE_BH_INFO && tempstar2.se_k == 14 && tempstar.se_k == 14)
                                parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g %g %g %g %g %g %g\n",
                                                          TotalTime, (isbinbin?"binary-binary":"binary-single"),
                                                          star_r[get_global_idx(knew)], hier.obj[i]->obj[bid]->obj[1]->id[0],hier.obj[i]->obj[bid]->obj[1]->id[nmerged], 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[0]) * units.mstar / FB_CONST_MSUN, 
                                                          binint_get_mass(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[nmerged]) * units.mstar / FB_CONST_MSUN,
                                                          binint_get_spins(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[0]), binint_get_spins(k, kp, hier.obj[i]->obj[bid]->obj[1]->id[nmerged]), 
					    tempstar.id, tempstar.m*units.mstar/MSUN,hier.obj[i]->obj[bid]->obj[1]->chi,hier.obj[i]->obj[bid]->obj[1]->vkick[nmerged],
                                                          sqrt(-2*star_phi[get_global_idx(knew)]) * (units.l/units.t) / 1.0e5,
														  hier.obj[i]->obj[bid]->obj[1]->a_merger[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[1]->e_merger[nmerged],
														  hier.obj[i]->obj[bid]->obj[1]->a_50M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[1]->e_50M[nmerged],
														  hier.obj[i]->obj[bid]->obj[1]->a_100M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[1]->e_100M[nmerged],
														  hier.obj[i]->obj[bid]->obj[1]->a_500M[nmerged]*cmc_units.l*units.l / FB_CONST_AU ,hier.obj[i]->obj[bid]->obj[1]->e_500M[nmerged]);
                            tempstar.se_k = 14;
                        } else{
                            merge_two_stars(&tempstar, &tempstar2, &tempstar, vs, curr_st);
                            /* FIXME: really we're supposed to add the kick to each binary
                               member separately, then calculate the systemic kick to the binary,
                               but hopefully this doesn't happen too much. */
                            if (sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]) != 0.0) {
                                wprintf("Adding merger-induced kick of %g km/s to binary CoM instead of binary member!\n",
                                    sqrt(vs[1]*vs[1]+vs[2]*vs[2]+vs[3]*vs[3]));
                            }
                            star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);					       
                            vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);
                        }
						nmerged++;
					}
					set_star_EJ(knew);
					
					cp_starSEvars_to_binmember(tempstar, star[knew].binind, 1);
					cp_starmass_to_binmember(tempstar, star[knew].binind, 1);

					binary[star[knew].binind].id2 = tempstar.id;
					//binary[star[knew].binind].id2 = star_get_id_new();
					/* log collision */
					binint_log_collision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id2,
						binary[star[knew].binind].m2, 
						star_r[get_global_idx(knew)],
						*(hier.obj[i]->obj[bid]->obj[1]), k, kp, binary[star[knew].binind].bse_kw[1]);

					if (WRITE_MORECOLL_INFO  && hier.obj[i]->obj[bid]->obj[1]->ncoll == 2){
						/*Elena: Creating a file with additional collision information */
						binint_log_morecollision(isbinbin?"binary-binary":"binary-single", 
						binary[star[knew].binind].id2, binary[star[knew].binind].m2,
			  			binary[star[knew].binind].rad2, binary[star[knew].binind].bse_kw[1],
			  			binary[star[knew].binind].bse_massc[1], binary[star[knew].binind].bse_menv[1], 							binary[star[knew].binind].bse_radc[1], binary[star[knew].binind].bse_renv[1],
			  			*(hier.obj[i]->obj[bid]->obj[1]),k,kp, W, rperi);}	
						
				}

				star_m[get_global_idx(knew)] = binary[star[knew].binind].m1 + binary[star[knew].binind].m2;

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

/**
* @brief Parallel version of simul_relax_new
*
* @return relaxation timestep
*/
double simul_relax_new(void)
{
	long si, k, j, p, N_LIMIT, simin, simax;
	double dt, dtmin=GSL_POSINF, DTrel=0.0, W, n_local;
	double Mv2ave, Mave, M2ave, sigma;

	N_LIMIT = clus.N_MAX;
	p = AVEKERNEL; //For this value, the results are very close to the original simul_relax() function.

	//MPI: Earlier the divion of stars among processors for this part was different, (and was similar to the original simul_relax function), and the one for the main code was different. But, in that case this function required communication with neighbors. So, it was changed such that both the main code and this function use the same kind of division of stars among processors. Now, stars are divided in sets of AVEKERNEL which is typically set to MIN_CHUNK_SIZE to avoid communication caused due to this function.
	//for (si=mpiBegin+p; si<mpiEnd-p; si+=2*p) {
	for (si=1+p; si<mpiEnd-mpiBegin+1-p; si+=2*p) {
		simin = si - p;
		simax = simin + (2 * p - 1);

		Mv2ave = 0.0;
		Mave = 0.0;
		M2ave = 0.0;
		for (k=simin; k<=simax; k++) {
			j = get_global_idx(k);
			double tmp = star_m[j] * madhoc;
			Mv2ave += tmp * (sqr(star[k].vr) + sqr(star[k].vt));
			Mave += tmp;
			M2ave += sqr(tmp);
		}
		//OPT: Remove double?
		Mv2ave /= (double) (2 * p);
		Mave /= (double) (2 * p);
		M2ave /= (double) (2 * p);
		
		/* sigma is the 3D velocity dispersion */
		sigma = sqrt(Mv2ave/Mave);
		/* average relative speed for a Maxwellian, from Binney & Tremaine */
		W = 4.0 * sigma / sqrt(3.0 * PI);

		/* Compute local density */
		n_local = calc_n_local(get_global_idx(si), p, clus.N_MAX);
		
		/* remember that code time units are t_cross * N/log(GAMMA*N) */
		/* this expression is from Freitag & Benz (2001), eqs. (8) and (9), we're just
		   inputting locally-averaged quantities */

		dt = sqr(2.0*THETASEMAX/PI) * (PI/32.0) * 
			cub(W) / ( ((double) clus.N_STAR) * n_local * (4.0 * M2ave) );

		dtmin = MIN(dtmin, dt);
	}

	double tmpTimeStart = timeStartSimple();
	MPI_Allreduce(&dtmin, &DTrel, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);     
	timeEndSimple(tmpTimeStart, &t_comm);

	return(DTrel);
}

/**
* @brief simulate relaxation to get timestep (original serial version). Calculates timestep for each star using some average quantities taken around the star, and then returns the minimum of these.
*
* @param rng gsl rng
*
* @return relaxation timestep
*/
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

/**
* @brief
   Here we need to compare the inspiral time for binary black holes to the typical
   timescale for another encounter.  Note that I'm not being careful with units, since we
   only need to compare timescales
*/

int destroy_bbh(double m1, double m2,double a,double e,double nlocal,double sigma,struct rng_t113_state* rng_st)
{
    /* if the binary isn't eccentric, don't bother*/
    if(e < 0.9)
        return 1;

    double T_bs, lambda_bs, T_gw, X, beta;
	double clight = 2.9979e10 / (units.l/units.t); 

    /* first compute the typical rate between encounters 
     * 4*sqrt(pi) = 7.089815
     * Note that sigma is the 3D velocity dispersion; we want 1D */
    lambda_bs = 7.089815 * nlocal * sqr(a) * (sigma/sqrt(3)) * (1 + madhoc*(m1+m2)/(2*a*sqr(sigma)/3));

    /* then the gravitational-wave timescale to merger (in the high e limit)
     * from Peters 1964 
     * 64/5 = 12.8
     * 768/425 = 1.807058 */
    beta = 12.8 * m1*m2*(m1+m2)*madhoc*madhoc*madhoc / pow(clight,5);
    T_gw = 1.807058 * (pow(a,4) / (4*beta)) * pow(1-sqr(e),3.5);

    /* Then compute a random time for the next encounter, taken from the 
     * exponential distribution */
	X = rng_t113_dbl_new(rng_st);
    T_bs = -log(1-X) / lambda_bs;

    /* now just compare the two; if T_bs < T_gw, break the binary*/
    if(T_bs < T_gw)
        return 1;
    else
        return 0;
    
}

/**
* @brief
   Since the binary interactions are done in a vaccuum, it is possible for them
   to produce pathologically wide binaries, which must be broken by hand, else 
   they shorten the timestep to a crawl.
*/
void break_wide_binaries(struct rng_t113_state* rng_st)
{
	long i,g_i,j, k, g_k, knew, knewp;
    int breakBinary = 0;
	double W, vorb, m, v2, Eexcess=0.0, exc_ratio, nlocal, llocal;
    double E_dump, E_dump_capacity, E_dump_factor=0.8;
    double length_factor = BINARY_DISTANCE_BREAKING;
    double Eexcess_prev, Eexcess_check;
    double hardness, mAveLocal, sigma2;
	MPI_Status stat;

	for (k=1; k<=clus.N_MAX_NEW; k++)
	{
		g_k = get_global_idx(k);

		/* Add in any excess energy from breaking triples in fewbody here  */
		Eexcess += star[k].E_excess;
		star[k].E_excess = 0.;

		if (star[k].binind) {

			/* binary index */
			j = star[k].binind;
			
			nlocal = calc_n_local(g_k, AVEKERNEL, clus.N_MAX);
            mAveLocal = sqrt(calc_average_mass_sqr(k,clus.N_MAX));
            sigma2 = sqr(sigma_array.sigma[k]);
			llocal = length_factor * pow(nlocal, -1.0/3.0);
            hardness = (binary[j].m1 * binary[j].m2 * sqr(madhoc)) /
                         (binary[j].a * mAveLocal * sigma2);

            /*if Binary_breaking_min is set, then use the hardness and MIN_BINARY_HARDNESS as a breaking criterion
            otherwise, use the length of the apoastron compared to the interparticule seperation as the criterion*/
            if(BINARY_BREAKING_MIN)
                breakBinary = (hardness <= MIN_BINARY_HARDNESS);
            else
                breakBinary = (binary[j].a*(1.0+binary[j].e) >= llocal);

            /* Special care must be taken for binary black holes, which can have very wide orbits
             * but very small merger timescales if the eccentricity is high */
            if(breakBinary && binary[j].bse_kw[0] == 14 && binary[j].bse_kw[1] == 14){
                breakBinary = destroy_bbh(binary[j].m1, binary[j].m2, binary[j].a, binary[j].e,nlocal,sqrt(sigma2),rng_st);
            }

			if (breakBinary){
				Eexcess += binary[j].m1 * binary[j].m2 * sqr(madhoc) / (2.0 * binary[j].a);

				/* create two stars for the binary components */
				knew = create_star(k, 0);
				knewp = create_star(k, 0);
				cp_binmemb_to_star(k, 0, knew);
				cp_binmemb_to_star(k, 1, knewp);
				
				/* destroy this binary */
				destroy_obj(k);
                breakBinary = 0;
            }
        } else if(Eexcess > 0. && star[k].interacted == 0 ){
            m = star_m[g_k];
            /* take excess energy from nearby field star (single or binary) */
            v2 = sqr(star[k].vt)+sqr(star[k].vr);
            E_dump_capacity = E_dump_factor * 0.5 * m * madhoc * v2; 
            E_dump = MIN(E_dump_capacity,Eexcess);
            exc_ratio = sqrt( (v2-2.0*E_dump/(m*madhoc) ) / v2 );
            star[k].vr *= exc_ratio;
            star[k].vt *= exc_ratio;
            set_star_EJ(k);
            Eexcess -= E_dump;
            star[k].interacted = 1;
		   }
	}
	double tmpTimeStart;
    MPI_Allreduce(&Eexcess, &Eexcess_check, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
    Eexcess_prev = 0.0;
    while(Eexcess_check > 0.)
    {
        tmpTimeStart = timeStartSimple();
        if(myid != procs-1)
            MPI_Send(&Eexcess, 1, MPI_DOUBLE, ( myid + 1 ), 0, MPI_COMM_WORLD);
        if(myid != 0)
            MPI_Recv(&Eexcess_prev, 1, MPI_DOUBLE, ( myid - 1), 0, MPI_COMM_WORLD, &stat);
        timeEndSimple(tmpTimeStart, &t_comm);

        Eexcess = Eexcess_prev;
        for (i = 1; i <= clus.N_MAX_NEW; i++) {
            g_i = get_global_idx(i);
            m = star_m[g_i];

            v2 = sqr(star[i].vr)+sqr(star[i].vt);
            if (star[i].interacted == 0) {
                    E_dump_capacity = E_dump_factor * 0.5 * m * madhoc * v2;
                    E_dump = MIN( E_dump_capacity, Eexcess );
                    if(Eexcess > 0) {
                        exc_ratio = sqrt( (v2 - 2 * E_dump / (m*madhoc)) / v2 );
                        star[i].vr *= exc_ratio;
                        star[i].vt *= exc_ratio;
                        Eexcess -= E_dump;
                        set_star_EJ(k);
                        star[k].interacted = 1;
                    }
            }
        }
        tmpTimeStart = timeStartSimple();
        MPI_Allreduce(&Eexcess, &Eexcess_check, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
        timeEndSimple(tmpTimeStart, &t_comm);
	}
	/* keep track of the energy that's vanishing due to our negligence */
	Eoops += -Eexcess;
}

/**
* @brief Computes the local average velocity dispersion value for each star (parallel version of calc_sigma_r).
*
* @param p the window to be averaged over
* @param N_LIMIT total number of stars in the processor
* @param sig_r_or_mave gets filled up with either the radial positions or average masses
* @param sig_sigma sigma array
* @param sig_n n value of sigma structure
* @param r_0_mave_1 if 0, sig_r_or_mave is filled with r values, if not, average mass values
*/
void calc_sigma_r(long p, long N_LIMIT, double *sig_r_or_mave, double *sig_sigma, long* sig_n, int r_0_mave_1)
{
	long si, k, simin, simax, siminlast, simaxlast;
	double Mv2ave, Mave;
	
//	N_LIMIT = mpiEnd-mpiBegin+1; //Its an input now.
	*sig_n = N_LIMIT;

	//MPI: Structure for ghost particles
	struct ghost_pts {
		double* prev;
		double* next;
	};

	struct ghost_pts ghost_pts_vr;
	struct ghost_pts ghost_pts_vt;

	//MPI: Needs p ghost particles on either side, so allocating memory
	ghost_pts_vr.prev = (double *) malloc(p * sizeof(double));
	ghost_pts_vr.next = (double *) malloc(p * sizeof(double));
	ghost_pts_vt.prev = (double *) malloc(p * sizeof(double));
	ghost_pts_vt.next = (double *) malloc(p * sizeof(double));
	//MPI: Buffer for communication
	double* buf_v = (double *) malloc(2 * p * sizeof(double));

	MPI_Status stat;

	/* MPI: Communicating ghost particles */
	double tmpTimeStart = timeStartSimple();

	//MPI: The 0 to p-1 elements of the buffer holds vr values, and p to 2p-1 holds the vt values.
	for(k=0; k<2*p; k++)
	{
		if( k < p )
			buf_v[k] = star[1 + k].vr;
		else
			buf_v[k] = star[1 + k - p].vt;
	}

	//OPT: Replace with MPI_Sendrecv. Does not optimize, but just code becomes compact.
	//MPI: Send to previous processor
	MPI_Send(buf_v, 2 * p, MPI_DOUBLE, ( myid + procs - 1 ) % procs, 0, MPI_COMM_WORLD);
	MPI_Recv(buf_v, 2 * p, MPI_DOUBLE, ( myid + 1) % procs, 0, MPI_COMM_WORLD, &stat);

	//MPI: Copying the received ghost particles to the allocated structure
	for(k=0; k<2*p; k++)
	{
		if(myid != procs-1)
		{
			if( k < p )
				ghost_pts_vr.next[k] = buf_v[k];
			else
				ghost_pts_vt.next[k-p] = buf_v[k];
		}
	}

	/*****************/

	//MPI: Forward communication, similar to above
	for(k=0; k<2*p; k++)
	{
		if( k < p )
			buf_v[k] = star[N_LIMIT - p + k + 1].vr;
		else
			buf_v[k] = star[N_LIMIT - p + k + 1 - p].vt;
	}

	MPI_Send(buf_v, 2 * p, MPI_DOUBLE, ( myid + 1 ) % procs, 0, MPI_COMM_WORLD);
	MPI_Recv(buf_v, 2 * p, MPI_DOUBLE, ( myid + procs - 1) % procs, 0, MPI_COMM_WORLD, &stat);

	for(k=0; k<2*p; k++)
	{
		if( myid != 0 )
		{
			if( k < p )
				ghost_pts_vr.prev[k] = buf_v[k];
			else
				ghost_pts_vt.prev[k-p] = buf_v[k];
		}
	}

	free(buf_v);
	timeEndSimple(tmpTimeStart, &t_comm);
	/* End of communication */

	siminlast = 1;//set to min index
	if (myid!=0)
		simaxlast = - p;
	else
		simaxlast = 0;

	Mv2ave = 0.0;
	Mave = 0.0;
	for (si=1; si<=N_LIMIT; si++) {

		//Also find the global index to figure out special cases
		int g_si = get_global_idx(si);
		simin = si - p;
		int g_simin = g_si - p;
		simax = simin + (2 * p - 1);
		int g_simax = g_simin + (2 * p - 1); 

		if (g_simin < 1) {
			//Special case for the root node
			simin = 1;
			simax = simin + (2 * p - 1);
		} else if (g_simax > clus.N_MAX) {
			//Special case for the last node
			simax = N_LIMIT;
			simin = simax - (2 * p - 1);
		}

		double vr=0.0, vt=0.0;
		// do sliding sum
		for (k=siminlast; k<simin; k++) {
			if (k < 1) {
				vr = ghost_pts_vr.prev[k+p-1];
				vt = ghost_pts_vt.prev[k+p-1];
			} else {
				vr = star[k].vr;
				vt = star[k].vt;
			}

			//MPI: Using a direct expression instead of get_global_idx() since it was changed to return the global index for stars outside local subset.
			int g_k = Start[myid] + k - 1; //get_global_idx(k);
			/*MPI: Using the global mass array*/
			Mv2ave -= star_m[g_k] * madhoc * (sqr(vr) + sqr(vt));
			Mave -= star_m[g_k] * madhoc;
		}

		for (k=simaxlast+1; k<=simax; k++) {
			int g_k = Start[myid] + k - 1; //get_global_idx(k);

			if (k > N_LIMIT) {
				vr = ghost_pts_vr.next[k-N_LIMIT-1];
				vt = ghost_pts_vt.next[k-N_LIMIT-1];
			} else if (k < 1) {
				vr = ghost_pts_vr.prev[k+p-1];
				vt = ghost_pts_vt.prev[k+p-1];
			} else {
				vr = star[k].vr;
				vt = star[k].vt;
			}

			/*MPI: Using the global mass array*/
			Mv2ave += star_m[g_k] * madhoc * (sqr(vr) + sqr(vt));
			Mave += star_m[g_k] * madhoc;
		}
		
		/* Storing r or average mass based on input parameter */
		if(r_0_mave_1 == 0)
			sig_r_or_mave[si] = star_r[g_si];
		else
			sig_r_or_mave[si] = Mave/2./p;

		/* store sigma (sigma is the 3D velocity dispersion) */
		sig_sigma[si] = sqrt(Mv2ave/Mave);
		
		siminlast = simin;
		simaxlast = simax;
	}
	free(ghost_pts_vr.prev);
	free(ghost_pts_vt.prev);
	free(ghost_pts_vr.next);
	free(ghost_pts_vt.next);
}


/**
* @brief calculates sliding averages of mass^2 around given index
*
* @param index index of star around which average is needed
* @param N_LIMIT total number of stars
*
* @return average of mass^2
*/
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
    M2ave += sqr(star_m[get_global_idx(k)] * madhoc);
  }
  M2ave /= (double) (2 * p);

  return(M2ave);
};

/**
* @brief generic routine for testing for approximate equality of floating point numbers
*
* @param a first floating point number
* @param b second floating point number
*
* @return ?
*/
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

/**
* @brief tidal capture (including merger) cross section; inputs are assumed to be in (self-consistent) code units
*
* @param n ?
* @param m1 ?
* @param r1 ?
* @param m2 ?
* @param vinf ?
*
* @return ?
*/
double sigma_tc_nd(double n, double m1, double r1, double m2, double vinf) {
	double a, beta=2.2, vstar1=sqrt(2.0*m1/r1);

	if (floateq(n, 1.5)) {
		a = 6.60 * pow(m2/m1, 0.242) + 5.06 * pow(m2/m1, 1.33);
	} else if (floateq(n, 3.0)) {
		a = 3.66 * pow(m2/m1, 0.200) + 2.94 * pow(m2/m1, 1.32);
	} else {
		eprintf("unknown polytropic index n=%g!\n", n);
		exit_cleanly(-1, __FUNCTION__);
		exit(1);
	}
	
	return(a*pow(vinf/vstar1,-beta)*r1*r1);
}

/**
* @brief tidal capture (including merger) cross section; inputs are assumed to be in (self-consistent) code units
*
* @param na ?
* @param ma ?
* @param ra ?
* @param nb ?
* @param mb ?
* @param rb ?
* @param vinf ?
*
* @return ?
*/
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
		exit_cleanly(-1, __FUNCTION__);
		exit(1);
	}
	
	return(a*pow(vinf/vstar1,-beta)*r1*r1);
}

/**
* @brief T_l function for use in tidal capture calculations
*
* @param order ?
* @param polytropicindex ?
* @param eta ?
*
* @return ?
*/
double Tl(int order, double polytropicindex, double eta)
{
	int l=order;
	double n=polytropicindex, x=log10(eta), x2=x*x, x3=x*x2, x4=x2*x2, x5=x2*x3;

	if (l != 2 && l != 3) {
		eprintf("unknown order l=%d\n", l);
		exit_cleanly(-1, __FUNCTION__);
		exit(1);
	}

        /* From Portegies Zwart & Meinen 1993 */
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
		exit_cleanly(-1, __FUNCTION__);
		exit(1);
	}
}

/**
* @brief tidal energy of Mpert acting on Mosc, in the polytropic approximation
*
* @param rperi ?
* @param Mosc ?
* @param Rosc ?
* @param nosc ?
* @param Mpert ?
*
* @return ?
*/
double Etide(double rperi, double Mosc, double Rosc, double nosc, double Mpert)
{
	double eta=sqrt(Mosc/(Mosc+Mpert))*pow(rperi/Rosc, 1.5);

	return(sqr(Mpert)/Rosc * (pow(Rosc/rperi, 6.0) * Tl(2, nosc, eta) + pow(Rosc/rperi, 8.0) * Tl(3, nosc, eta)));
}

/**
* @brief Add kicks to vt
*
* @param vt transverse velocity
* @param vs1 ?
* @param vs2 ?
* @param rng_st random state
*/
void vt_add_kick(double *vt, double vs1, double vs2, struct rng_t113_state* rng_st)
{
	double X, theta, vtx, vty;
	//X = rng_t113_dbl();
	X = rng_t113_dbl_new(rng_st);
	theta = 2.0 * PI * X;
	vtx = *vt*sin(theta);
	vty = *vt*cos(theta);
	vtx = vtx + vs1*1.0e5/(units.l/units.t);
	vty = vty + vs2*1.0e5/(units.l/units.t);
	*vt = sqrt(vtx*vtx+vty*vty);
}


/**
* @brief add GW recoil kicks and mass loss for mergers of BBHs
*
* These are all based on fits to NR simulations, though the functional forms
* I've taken from Section V of Gerosa and Kesden, PRD, 93, 12, 124066 (2016)
*/
void binary_bh_merger(long k, long kb, long knew, int kprev0, int kprev1, struct rng_t113_state* rng_st){
	double m1 = binary[kb].m1*units.mstar / FB_CONST_MSUN;//Only need the mass ratio anyway
	double m2 = binary[kb].m2*units.mstar / FB_CONST_MSUN;
	double chi1 = binary[kb].bse_bhspin[0];
	double chi2 = binary[kb].bse_bhspin[1];
	double afinal, mass_frac, v_para, v_perp;
    double X, theta, phi, vk, vx, vy, vz;

	/* The actual kick routine is in fewbody_coll.c */
	fb_bh_merger(m1,m2,chi1,chi2,&mass_frac,&afinal,&v_para,&v_perp,rng_st);

	/*Finally, apply these to the newly formed (single!) BH*/
    /* Pick a random 3D vector for the kick (probably overkill...)*/
	X = rng_t113_dbl_new(rng_st);
	phi = 2.0 * PI * X;
	X = rng_t113_dbl_new(rng_st);
    theta = acos(2*X - 1.);
    vk = sqrt(v_para*v_para + v_perp*v_perp);
    vx = vk*sin(theta)*sin(phi);
	vy = vk*sin(theta)*cos(phi);
    vz = vk*cos(theta);

    /* Then add the 3D vector to the stars velocity
     * Note that the kick is in km/s; convert to CGS then code units*/
    star[knew].vr += vx * 1.0e5 / (units.l/units.t);					       
	vt_add_kick(&(star[knew].vt),vy,vz,rng_st);

	star[knew].m *= mass_frac;
    star[knew].se_mt *= mass_frac;
	star[knew].rad = (binary[kb].rad1 + binary[kb].rad2) * mass_frac; 
	star[knew].se_bhspin = afinal;

    if(WRITE_BH_INFO)
        parafprintf(bhmergerfile, "%.18g %s %g %ld %ld %g %g %g %g %ld %g %g %g %g %g %g -100 -100 -100 -100 -100 -100\n", TotalTime, "isolat-binary",
		    star_r[get_global_idx(knew)], binary[kb].id1,binary[kb].id2, m1,m2,chi1,chi2,star[knew].id,
                                        (m1+m2)*mass_frac, afinal,vk, 
                                        sqrt(-2*star_phi[get_global_idx(knew)])*(units.l/units.t) / 1.0e5, 
                                        binary[kb].a*units.l/AU,binary[kb].e);
}

