/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief ?
*
* @param fname file name
*/
void create_rwalk_file(char *fname) {

#ifdef USE_MPI
    MPI_File mpi_rwalk_file;
    char mpi_rwalk_file_buf[10000];
    char mpi_rwalk_file_wrbuf[10000000];
    int mpi_rwalk_file_len=0, mpi_rwalk_file_ofst_total=0;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_rwalk_file);
	 if(tcount==1)
		 MPI_File_set_size(mpi_rwalk_file, 0);
#else
  FILE *rwalk_file= NULL;
  rwalk_file= fopen(fname, "a");
#endif

  pararootfprintf(rwalk_file, "\n");
  pararootfprintf(rwalk_file, 
          "# 1:index, 2:Time, 3:r, 4:Trel, 5:dt, 6:l2_scale, 7:n_steps, 8:beta 9:n_local, 10:W, 11:P_orb, 12:n_orb\n");
#ifdef USE_MPI
  mpi_para_file_write(mpi_rwalk_file_wrbuf, &mpi_rwalk_file_len, &mpi_rwalk_file_ofst_total, &mpi_rwalk_file);
  MPI_File_close(&mpi_rwalk_file);
#else
  fclose(rwalk_file);
#endif
}

/**
* @brief ?
*
* @param fname ?
* @param index star index
* @param Trel ?
* @param dt ?
* @param l2_scale ?
* @param n_steps ?
* @param beta ?
* @param n_local ?
* @param W ?
* @param P_orb ?
* @param n_orb ?
*/
void write_rwalk_data(char *fname, long index, double Trel, double dt, 
    double l2_scale, double n_steps, double beta, double n_local, double W, 
    double P_orb, double n_orb) {

#ifdef USE_MPI
	double r = star_r[get_global_idx(index)];
    MPI_File mpi_rwalk_file;
    char mpi_rwalk_file_buf[10000];
    char mpi_rwalk_file_wrbuf[10000000];
    int mpi_rwalk_file_len=0, mpi_rwalk_file_ofst_total=0;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_rwalk_file);
	 if(tcount==1)
		 MPI_File_set_size(mpi_rwalk_file, 0);
#else
	 double r = star[index].r;
  FILE *rwalk_file= NULL;
  rwalk_file= fopen(fname, "a");
#endif

  parafprintf(rwalk_file, "%li %g %g %g %g %g %g %g %g %g %g %g\n", 
      index, TotalTime, r, Trel, dt, sqrt(l2_scale), n_steps, beta, n_local, W, P_orb, n_orb);

#ifdef USE_MPI
  mpi_para_file_write(mpi_rwalk_file_wrbuf, &mpi_rwalk_file_len, &mpi_rwalk_file_ofst_total, &mpi_rwalk_file);
  MPI_File_close(&mpi_rwalk_file);
#else
  fclose(rwalk_file);
#endif
}

/**
* @brief This is the random walk procedure as outlined by Freitag & Benz (2002). Change of notation: beta here is theta in the paper.
*
* @param index star index
* @param v[4] ?
* @param vcm[4] ?
* @param beta ?
* @param dt ?
*/
void bh_rand_walk(long index, double v[4], double vcm[4], double beta, double dt)
{ 
	double w[3], n_orb, P_orb, deltabeta_orb, L2, Rdisr, Jlc, vlc;
	double deltamax, deltasafe, delta, dbeta;
	double w_mag, l2_scale;
	int i;
    int g_index;
#ifdef USE_MPI
	g_index = get_global_idx(index);
#else
    g_index = index;
#endif

#ifdef EXPERIMENTAL
	char fname[80];
	long is_in_ids;
	double Trel, n_local, M2ave; 
	double W, n_steps= 1.;

	is_in_ids= 0;
	sprintf(fname, "%s.rwalk_steps.dat", outprefix);
	n_local= calc_n_local(g_index, AVEKERNEL, clus.N_MAX);
	W = 4.0 * sigma_array.sigma[index] / sqrt(3.0*PI);
	M2ave= calc_average_mass_sqr(g_index, clus.N_MAX);
	Trel= (PI/32.)*cub(W)/ ( ((double) clus.N_STAR) * n_local * (4.0* M2ave) );
	//if (g_hash_table_lookup(star_ids, &star[index].id)!=NULL) {
	if (index==1 && TotalTime>= SNAPSHOT_DELTAT*(StepCount) && SNAPSHOTTING && WRITE_RWALK_INFO) {
		is_in_ids=1;
		create_rwalk_file(fname);
	};
#endif
	/* simulate loss cone physics for central mass */
	//MPI: Parallelized, but might have mistakes since I am not clear as to what some functions are doing.
	P_orb = calc_P_orb(index);
	n_orb = dt * ((double) clus.N_STAR)/log(GAMMA * ((double) clus.N_STAR)) / P_orb; 
	l2_scale= 1.;
#ifdef EXPERIMENTAL
	/* scale down L2 if the time step is larger than BH_LC_FDT*Trel */
	/* This is inconsistent, as for stars with dt< BH_LC_FDT*Trel the probability
	 * of hitting the loss cone becomes smaller, compared to the case with 
	 * dt=BH_LC_FDT*Trel
	 */
	/* if (BH_LC_FDT>0. && dt> BH_LC_FDT*Trel) { */
	if (BH_LC_FDT>0.) {
		n_steps= dt/BH_LC_FDT/Trel;
		l2_scale= 1./n_steps;
	};
#endif
	deltabeta_orb = 1.0/sqrt(n_orb) * sqrt(l2_scale)*beta;
	L2 = l2_scale*fb_sqr(beta);
	if (BH_R_DISRUPT_NB>0.) {
		Rdisr= BH_R_DISRUPT_NB;
	} else if (STELLAR_EVOLUTION){
		double Rss;
		//dprintf("cenma.m= %g, star[%li].m= %g\n", cenma.m, index, star[index].m);
#ifdef USE_MPI
		Rdisr= pow(2.*cenma.m/star_m[g_index], 1./3.)*star[index].rad*RSUN/units.l;
#else
		Rdisr= pow(2.*cenma.m/star[index].m, 1./3.)*star[index].rad*RSUN/units.l;
#endif
		Rss= 4.24e-06*cenma.m/SOLAR_MASS_DYN*RSUN/units.l;
		Rdisr= MAX(Rdisr, Rss);
	} else {
#ifdef USE_MPI
		Rdisr= pow(2.*cenma.m/star_m[g_index], 1./3.)*star[index].rad;
#else
		Rdisr= pow(2.*cenma.m/star[index].m, 1./3.)*star[index].rad;
#endif
	};
	Jlc= sqrt(2.*cenma.m*madhoc*Rdisr);
#ifdef USE_MPI
	vlc= Jlc/star_r[g_index];
#else
	vlc= Jlc/star[index].r;
#endif
	for (i=0; i<3; i++) {
		w[i]= v[i+1]- vcm[i+1];
	}
	w_mag= sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
	delta= 0.0;
	while (L2 > 0.0) { 
		L2 -= fb_sqr(delta); 
		if (sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(w[1]+vcm[2])) <= vlc) { 
			dprintf("index=%ld, id=%ld: star eaten by BH\n", g_index, star[index].id);
#ifdef USE_MPI
			cenma.m_new += star_m[g_index]; 
#else
			cenma.m_new += star[index].m; 
#endif
			destroy_obj(index);
			L2 = 0.0; 
		} else { 
			deltamax= 0.1*FB_CONST_PI;
			deltasafe= CSAFE*(sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(vcm[2]+w[1]))-vlc)/w_mag;
			delta = MAX(deltabeta_orb, MIN(deltamax, MIN(deltasafe, sqrt(L2)))); 
			//if (delta>sqrt(L2)) delta=sqrt(L2);
			//delta = MAX(deltabeta_orb, MIN(deltamax, sqrt(L2)));

#ifndef USE_MPI
			curr_st = &st[findProcForIndex(index)];
#endif
			dbeta = 2.0 * PI * rng_t113_dbl_new(st); 
			//dbeta = 2.0 * PI * rng_t113_dbl(); 

			do_random_step(w, dbeta, delta); 
#ifdef EXPERIMENTAL
			if (is_in_ids) {
				//fprintf(rwalk_file, "%f %g %g %g %g %g %g %g %g %g\n",
				//TotalTime, deltabeta_orb, deltasafe, sqrt(L2), delta, n_orb, beta, star[index].r, dt/Trel, l2_scale);
			};
#endif
		} 
	}; 
#ifdef EXPERIMENTAL
	if (TotalTime>= SNAPSHOT_DELTAT*(StepCount) && SNAPSHOTTING && WRITE_RWALK_INFO) {
		write_rwalk_data(fname, g_index, Trel, dt, l2_scale, n_steps, beta,
				n_local, W, P_orb, n_orb);
	}
#endif
};


/**
* @brief ?
*
* @param w ?
* @param vr ?
* @param vt ?
*/
void get_3d_velocities(double *w, double vr, double vt) {
   double phi;

	//MPI: Dont know how what to do with this rng call, but this function itself is not used anywhere! So, we can safely ignore this as of now.
   phi= rng_t113_dbl()*2.*FB_CONST_PI;
   w[0]= vt* cos(phi);
   w[1]= vt* sin(phi);
   w[2]= vr;
};

/* here the notation of Freitag & Benz (2002) is used */
/**
* @brief ?
*
* @param w ?
* @param beta ?
* @param delta ?
*/
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

/**
* @brief calculate the angle between w and w_new and compare to deltat
*
* @param w ?
* @param w_new ?
* @param delta ?
*
* @return ?
*/
double check_angle_w_w_new(double *w, double *w_new, double delta) {
   double angle;

   angle=acos(w[0]*w_new[0]+w[1]*w_new[1]+w[2]*w_new[2]);

   return(angle-delta);
};

/**
* @brief calculate star's radial orbital period
*
* @param index star index
*
* @return star's radial orbital period
*/
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
    int g_index;
#ifdef USE_MPI
    g_index = get_global_idx(index);
#else
    g_index = index;
#endif

        /* default values for star_interval */
	star_interval.min= 1;
    star_interval.max= clus.N_MAX+1;

#ifdef USE_MPI
	E = star[index].E + MPI_PHI_S(star_r[g_index], g_index);
#else
	E = star[index].E + PHI_S(star[index].r, index);
#endif
	J = star[index].J;
	
	//dprintf("index=%ld ", index);

#ifdef EXPERIMENTAL
        //if (index>80000) dprintf("Aaaaahhh index= %li\n", index);
	orbit_rs = calc_orbit_new(index, E, J);
#else
	orbit_rs = calc_orbit_rs(index, E, J);
#endif
	
	//dprintf("rp=%g ra=%g ", orbit_rs.rp, orbit_rs.ra);

	Porbapproxmin = 2.0 * PI * orbit_rs.rp / (J / orbit_rs.rp);
	Porbapproxmax = 2.0 * PI * orbit_rs.ra / (J / orbit_rs.ra);
	Porbapprox = sqrt(Porbapproxmin * Porbapproxmax);

	/* Return the approximate value of the radial period.  If this is commented out
	   the radial period will be calculated properly. */
	//return(Porbapprox);

        if (orbit_rs.ra-orbit_rs.rp< CIRC_PERIOD_THRESHOLD) {
          dprintf("Orbit is considered circular for period calculation.\n");
          dprintf("ra-rp= %g and is less than the threshold %g\n", 
              orbit_rs.ra-orbit_rs.rp, CIRC_PERIOD_THRESHOLD);
          orbit_rs.circular_flag = 1;
        }

	if (orbit_rs.circular_flag == 1) {
		/* We're returning the azimuthal period here, which is not the same as the
		   radial period for the general cluster potential.  This shouldn't make
		   any difference for the BH loss cone stuff, since the orbit is circular. */
#ifdef USE_MPI
		return(2.0 * PI * star_r[g_index] / star[index].vt);
#else
		return(2.0 * PI * star[index].r / star[index].vt);		
#endif
	} else {
		w = gsl_integration_workspace_alloc(1000);
		tab = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0.0, 0.0);

		params.E = E;
		params.J = J;
		params.index = g_index;
		//if (SEARCH_GRID)
		//  star_interval= search_grid_get_interval(r_grid, orbit_rs.rp);
		//params.kmin = FindZero_r(star_interval.min, star_interval.max, orbit_rs.rp);
		//if (SEARCH_GRID)
		//  star_interval= search_grid_get_interval(r_grid, orbit_rs.ra);
		//params.kmax = FindZero_r(star_interval.min, star_interval.max, orbit_rs.ra) + 1;
		params.kmax= orbit_rs.kmax+1;
		params.kmin= orbit_rs.kmin;
		params.rp = orbit_rs.rp;
		params.ra = orbit_rs.ra;
		F.params = &params;

                /* test if the interval of rmax is the same as [kmax,kmax+1] */
		if (!orbit_rs.circular_flag) {
			if (function_Q(g_index, params.kmax, E, J)>0 || function_Q(g_index, params.kmax-1,E, J)<0) {
			  dprintf("r and phi interval do not match: id= %li, r_kmax= %li\n",
			    star[index].id, params.kmax-1);
			  dprintf("star index is %li\n", index);
			  dprintf("f_Q[r_kmax]= %g; f_Q[r_kmax+1]= %g\n", function_Q(g_index, params.kmax-1, E, J),
			    function_Q(g_index, params.kmax, E, J));
			  dprintf("phi_kmax= %li, phi_kmax+1= %li\n", orbit_rs.kmax, orbit_rs.kmax+1);
			  dprintf("f_Q[phi_kmax]= %g; f_Q[phi_kmax+1]= %g\n", 
			    function_Q(g_index, orbit_rs.kmax, E, J),
			    function_Q(g_index, orbit_rs.kmax+1, E, J));
			  dprintf("(r[r_kmax]-rmax=%g, r[r_kmax+1]-rmax)= %g\n", 
			    star[params.kmax-1].r-orbit_rs.ra,star[params.kmax].r-orbit_rs.ra);
			};
                        if (calc_vr(params.rp, index, E, J)< 0.) {
                          dprintf("Harrrg: vr(rmin)< 0.! Damn it! Index: %li, Id: %li\n", index, 
                            star[index].id);
                        };
                        if (calc_vr(params.ra, index, E, J)< 0.) {
                          dprintf("Harrrg: vr(rmax)< 0.! Damn it! Index: %li, Id: %li\n", index, 
                            star[index].id);
                        };
                        if (params.kmax!=orbit_rs.kmax+1) 
                          dprintf("kmax in orbit_rs and params differ! kmax_o= %li, kmax_p=%li, Index: %li, Id: %li\n", 
                            orbit_rs.kmax, params.kmax, index, star[index].id);
                        if ((params.kmin!=orbit_rs.kmin)&& (params.kmin>1)) 
                          dprintf("kmin in orbit_rs and params differ! kmin_o= %li, kmin_p=%li, Index: %li, Id: %li\n", 
                            orbit_rs.kmin, params.kmin, index, star[index].id);
		};

        //MPI: These seem to be never used, so they are not parallelized as of now. Also unclear how these functions will be used, so dont understand if index transformation is reqd or not. Might later have to double check if the parallelization is correct.
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

/**
* @brief integrand for calc_P_orb
*
* @param x ?
* @param params ?
*
* @return ?
*/
double calc_p_orb_f(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

#ifdef USE_MPI
	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (potential(x) + MPI_PHI_S(x, myparams.index));
#else
	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (potential(x) + PHI_S(x, myparams.index));
#endif

	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/**
* @brief integrand for calc_P_orb
*
* @param x ?
* @param params ?
*
* @return ?
*/
double calc_p_orb_f2(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

#ifdef USE_MPI
	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (fastpotential(x, myparams.kmin, myparams.kmax) + MPI_PHI_S(x, myparams.index));
#else
	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (fastpotential(x, myparams.kmin, myparams.kmax) + PHI_S(x, myparams.index));
#endif

	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/**
* @brief integrand for calc_P_orb, using Gauss-Chebyshev for regularizing the integrand near the endpoints
*
* @param x ?
* @param params ?
*
* @return ?
*/
double calc_p_orb_gc(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;
	double phik, phik1, phi0, phi1, rk, rk1, rminus, rplus;
	double E, J, rp, ra;
	double result;
	long index, kmin, kmax;

	E = myparams.E;
	J = myparams.J;
	index = myparams.index;
	kmin = myparams.kmin;
	kmax = myparams.kmax;
	rp = myparams.rp;
	ra = myparams.ra;
        
#ifdef USE_MPI
	if (x <= star_r[kmin+1]) { /* return integrand regularized at r=rp */
#else
	if (x <= star[kmin+1].r) { /* return integrand regularized at r=rp */
#endif
		//dprintf("regularizing near rp...\n");
#ifdef USE_MPI
		phik = star_phi[kmin] + MPI_PHI_S(star_r[kmin], index);
		phik1 = star_phi[kmin+1] + MPI_PHI_S(star_r[kmin+1], index);
		rk = star_r[kmin];
		rk1 = star_r[kmin+1];
#else
		phik = star[kmin].phi + PHI_S(star[kmin].r, index);
		phik1 = star[kmin+1].phi + PHI_S(star[kmin+1].r, index);
		rk = star[kmin].r;
		rk1 = star[kmin+1].r;
#endif 
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1);
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1);
		/* rplus = rperi, rminus = rapo-primed */
		if (E-phi0==0.) {
		  dprintf("E is phi0 near rp! Damn it!");
		};
		rminus = (phi1 - sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//rplus = (phi1 + sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//dprintf("rplus/rp=%g rminus/ra=%g\n", rplus/rp, rminus/ra);
		if (kmax == kmin + 1) {
			/* then rminus = ra, so must cancel (ra-x)/(rminus-x) term analytically */
			result= 2.0*x*sqrt(1.0/((2.0*phi0-2.0*E)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near rp! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near rp! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
		        result= 2.0*x*sqrt((ra-x)/((2.0*phi0-2.0*E)*(rminus-x)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near rp! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near rp! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			return(2.0*x*sqrt((ra-x)/((2.0*phi0-2.0*E)*(rminus-x))));
		}
#ifdef USE_MPI
	} else if (x >= star_r[kmax-1]) { /* return integrand regularized at r=ra*/
#else
	} else if (x >= star[kmax-1].r) { /* return integrand regularized at r=ra*/
#endif
		//dprintf("regularizing near ra...\n");
#ifdef USE_MPI
		phik = star_phi[kmax-1] + MPI_PHI_S(star_r[kmax-1], index);
		phik1 = star_phi[kmax] + MPI_PHI_S(star_r[kmax], index);
		rk = star_r[kmax-1];
		rk1 = star_r[kmax];
#else
		phik = star[kmax-1].phi + PHI_S(star[kmax-1].r, index);
		phik1 = star[kmax].phi + PHI_S(star[kmax].r, index);
		rk = star[kmax-1].r;
		rk1 = star[kmax].r;
#endif
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1);
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1);
		/* rplus = rperi-primed, rminus = rapo */
		//rminus = (phi1 - sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		rplus = (phi1 + sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		//dprintf("rplus/rp=%g rminus/ra=%g\n", rplus/rp, rminus/ra);
		if (kmax == kmin + 1) {
			result=2.0*x*sqrt(1.0/((2.0*phi0-2.0*E)));
			/* then rplus = rp, so must cancel (x-rp)/(x-rplus) term analytically */
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
			result= 2.0*x*sqrt((x-rp)/((2.0*phi0-2.0*E)*(x-rplus)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near ra! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			return(2.0*x*sqrt((x-rp)/((2.0*phi0-2.0*E)*(x-rplus))));
		}
	} else {
#ifdef USE_MPI
		radicand = 2.0 * (E - (potential(x) + MPI_PHI_S(x, index)))- fb_sqr(J/x);
#else
		radicand = 2.0 * (E - (potential(x) + PHI_S(x, index)))- fb_sqr(J/x);
#endif
		if (radicand < 0.0) {
			dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, index);
			dprintf("kmin= %li, kmax= %li, rp=%g, ra=%g, Id: %li\n",
			  kmin, kmax, rp, ra, star[index].id);
			radicand = 0.0;
		};
		result= 2.0 * sqrt((x-rp)*(ra-x)/radicand);
		if (gsl_isinf(result)) {
		  dprintf("result is infinite! Damn it!\n");
		  dprintf("kmax=%li, kmin=%li, index=%li\n",
		    kmax, kmin, index);
		  dprintf("E= %g, J=%g, id=%li\n", E, J, star[index].id);
  	          dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
#ifdef USE_MPI
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star_r[kmin+1]);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star_r[kmax-1]);
#else
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star[kmin+1].r);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star[kmax-1].r);
#endif
		};
		if (gsl_isnan(result)) {
		  dprintf("result is NaN! Damn it!\n");
		  dprintf("kmax=%li, kmin=%li, index=%li\n",
		    kmax, kmin, index);
		  dprintf("E= %g, J=%g, id=%li\n", E, J, star[index].id);
  	          dprintf("x=%g, ra= %g, rp= %g, radicand= %g\n", x, ra, rp, radicand);
#ifdef USE_MPI
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star_r[kmin+1]);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star_r[kmax-1]);
#else
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star[kmin+1].r);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star[kmax-1].r);
#endif
		};
		return(2.0 * sqrt((x-rp)*(ra-x)/radicand));
	}
}

/**
* @brief ?
*
* @param r ?
*
* @return ?
*/
struct Interval get_r_interval(double r) {
  long kmax, kmin, i;
  struct Interval star_interval;

  if (SEARCH_GRID) {
   star_interval= search_grid_get_interval(r_grid, r);
   kmax= star_interval.max;
   kmin= star_interval.min;
  } else {
   kmax= clus.N_MAX+1;
   kmin= 1;
  };
  if (kmin==kmax-1) {
   i= kmin;
  } else {
   i =  FindZero_r(kmin, kmax, r);
  };

  star_interval.max= i+1;
  star_interval.min= i;
  return (star_interval);
}


