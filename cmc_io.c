/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief print the version
*
* @param stream stream to be printed to
*/
void print_version(FILE *stream)
{
	fprintf(stream, "** %s %s (%s) [%s] **\n", CMCPRETTYNAME, CMCVERSION, CMCNICK, CMCDATE);
}

/**
* @brief print the usage
*
* @param stream stream to be printed to
* @param argv[] input arg list
*/
void cmc_print_usage(FILE *stream, char *argv[])
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  %s [options...] <input_file> <output_file_prefix> <old_output_file_prefix[if -R specified]>\n", argv[0]);
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -q --quiet   : do not print diagnostic info to stdout\n");
	fprintf(stream, "  -d --debug   : turn on debugging\n");
	fprintf(stream, "  -V --version : print version info\n");
	fprintf(stream, "  -R --hard-restart : start from a saved checkpoint (specify number) with prefix <old...prefix> and write to new file prefix\n");
	fprintf(stream, "  -r --soft-restart : start from a saved checkpoint (specify number) and write to the same files in the same place\n");
	fprintf(stream, "  -h --help    : display this help text\n");
}

/**
* @brief Writes output to stdout and all files required for post-simulation analysis.
*/
void print_results(void){
#ifdef USE_MPI
	//MPI: Writes out files that need contribution from all/many processors.
    PrintParaFileOutput();
#endif
	 //MPI: These two routines mostly write files that are need data only from the root node.
    PrintLogOutput();
    PrintFileOutput();
    fflush(NULL);
}

/*********** Output 2D/3D snapshots **************/
/**
* @brief prints a 2D snapshot
*/
void print_2Dsnapshot(void)
{
	long i, j;
	j=0;
	char outfile[100];

	if (SNAPSHOTTING) {
		// open file for 2D snapshot 
		sprintf(outfile, "%s.snap%04ld.dat.gz", outprefix, snap_num);
		write_snapshot(outfile, 0);
		// global counter for snapshot output file 
		snap_num++;
	}
}


/**
* @brief prints BH snapshots
*/
void print_bh_snapshot(void) {
	long i, j;
	char outfile[100];
	
	if (BH_SNAPSHOTTING) {
		/* open file for BH snapshot */
	
		sprintf(outfile, "%s.bhinfo%04ld.dat.gz", outprefix, bh_snap_num);
		write_snapshot(outfile, 1);

		/* global counter for snapshot output file */
		bh_snap_num++;
	}
}



/**
* @brief Calculates some quantities required for output and prints out log files. These are mostly files that need data only from the root node to be written out.
*/
void PrintLogOutput(void)
{
	double m, rh, trh, conc_param, m_single, m_binary;
	long ih, k;
	
	/* Computing half-mass radii, and relaxation time */
	m = rh = trh = 0.0;
	rh_binary = rh_single = m_binary = m_single = 0.0;

#ifdef USE_MPI
	//MPI: Do computations on duplicated arrays separately.
	for (ih=1; ih<=clus.N_MAX; ih++) {
		k = ih;

		m += star_m[k] / clus.N_STAR;
		
		if (m/Mtotal <= 0.5) {
			rh = star_r[k];
		}
    }

	//MPI: Arrays to store all intermediate values while cumulating these sums.
	double *m_binary_arr = (double*) calloc(clus.N_MAX_NEW+1, sizeof(double));
	double *m_single_arr = (double*) calloc(clus.N_MAX_NEW+1, sizeof(double));

	for (ih=1; ih<=clus.N_MAX_NEW; ih++) {
		k = ih;
		int g_k = get_global_idx(k);

		if (star[k].binind > 0) {
			m_binary_arr[k] = m_binary_arr[k-1] + star_m[g_k] / clus.N_STAR;
			m_single_arr[k] = m_single_arr[k-1];
		} else {
			m_single_arr[k] = m_single_arr[k-1] + star_m[g_k] / clus.N_STAR;
			m_binary_arr[k] = m_binary_arr[k-1];
		}
    }

    double buf_comm[2];
    double buf_comm_recv[2];
    buf_comm[0] = m_binary_arr[clus.N_MAX_NEW];
    buf_comm[1] = m_single_arr[clus.N_MAX_NEW];

	 //MPI: Reduce to find out the total sums.
    //MPI: Since only the root node will print out these stuff, Allreduce is not reqd, just Reduce will do.
    //MPI_Allreduce(buf_comm, buf_comm_recv, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	 double tmpTimeStart = timeStartSimple();
    MPI_Reduce(buf_comm, buf_comm_recv, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);
    if(myid==0)
    {
        m_binary = buf_comm_recv[0];
        m_single = buf_comm_recv[1];
    }

	 //MPI: Find out cumulative sums of values in all processors with ids less than current processor.
	 tmpTimeStart = timeStartSimple();
    MPI_Exscan(buf_comm, buf_comm_recv, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);

	 for (ih=1; ih<=clus.N_MAX_NEW; ih++) {
		k = ih;
		int g_k = get_global_idx(k);

		//MPI: Now, we can find out when this condition is satisfied in each processor with the intermediate values.
		if ((m_single_arr[k] + buf_comm_recv[1]) / (Mtotal - (M_b / clus.N_STAR)) <= 0.5) {
			rh_single = star_r[g_k];
		}

		// avoid dividing by zero if there are no binaries
		if (M_b > 0) {
			if ((m_binary_arr[k] + buf_comm_recv[0])/ M_b * clus.N_STAR <= 0.5) {
				rh_binary = star_r[g_k];
			}
		}
	 }

    buf_comm[0] = rh_binary;
    buf_comm[1] = rh_single;

    //MPI: Since r's are always monotonically increasing since they are sorted, I can just take the max.
    //MPI_Allreduce(buf_comm, buf_comm_recv, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	 tmpTimeStart = timeStartSimple();
    MPI_Reduce(buf_comm, buf_comm_recv, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);

    rh_binary = buf_comm_recv[0];
    rh_single = buf_comm_recv[1];

	 free(m_binary_arr);
	 free(m_single_arr);
#else
	for (ih=1; ih<=clus.N_MAX; ih++) {
		k = ih;
		m += star[k].m / clus.N_STAR;
		
		if (star[k].binind > 0) {
			m_binary += star[k].m / clus.N_STAR;
		} else {
			m_single += star[k].m / clus.N_STAR;
		}
		
		if (m/Mtotal <= 0.5) {
			rh = star[k].r;
		}
		if (m_single / (Mtotal - (M_b / clus.N_STAR)) <= 0.5) {
			rh_single = star[k].r;
		}
		// avoid dividing by zero if there are no binaries
		if (M_b > 0) {
			if (m_binary / M_b * clus.N_STAR <= 0.5) {
				rh_binary = star[k].r;
			}
		}
	}
#endif

	/* t_rh calculated using r_h */
	trh = ((0.138 * clus.N_MAX) / log((double) clus.N_MAX)) * sqrt((rh * rh * rh) / Mtotal) * log((double) clus.N_MAX) / clus.N_MAX;

	/* Concentration parameter --- note that max_r is the max radius of all bound
	   stars, returned by get_positions(). When a finite R_MAX (i.e. initial Rtidal)
	   is specified, max_r should be approximately equal to Rtidal. But when 
	   Rtidal is very large (isolated cluster), r_max still gives the "actual"
	   size of the cluster, indicating the maximum radius.
	   The conc parameter is actually defined as LOG_10(rtidal/rcore). But for
	   easier reading, it is calculated here only as rtidal/rcore. */
	
	/* max_r & core_radius are output in out_3 file */
	if(max_r == 0.0){
		conc_param = 0.0;
	} else {
		conc_param = (max_r / core_radius);
	}

	//MPI: Print out files.
	rootgprintf("******************************************************************************\n");
	pararootfprintf(logfile, "******************************************************************************\n");

	rootgprintf("tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);
	pararootfprintf(logfile, "tcount=%ld TotalTime=%.16e Dt=%.16e\n", tcount, TotalTime, Dt);

	rootgprintf("Etotal=%g max_r=%g N_bound=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);
	pararootfprintf(logfile, "Etotal=%g max_r=%g N_bound=%ld Rtidal=%g\n", Etotal.tot, max_r, clus.N_MAX, Rtidal);

	rootgprintf("Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);
	pararootfprintf(logfile, "Mtotal=%g Etotal.P=%g Etotal.K=%g VRatio=%g\n", Mtotal, Etotal.P, Etotal.K, -2.0 * Etotal.K / Etotal.P);

	rootgprintf("TidalMassLoss=%g\n", TidalMassLoss);
	pararootfprintf(logfile, "TidalMassLoss=%g\n", TidalMassLoss);
	
	rootgprintf("core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	pararootfprintf(logfile, "core_radius=%g rho_core=%g v_core=%g Trc=%g conc_param=%g N_core=%g\n",
		core_radius, rho_core, v_core, Trc, conc_param, N_core);
	
	rootgprintf("trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	pararootfprintf(logfile, "trh=%g rh=%g rh_single=%g rh_binary=%g\n", trh, rh, rh_single, rh_binary);
	
	rootgprintf("N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);
	pararootfprintf(logfile, "N_b=%ld M_b=%g E_b=%g\n", N_b, M_b/clus.N_STAR, E_b);

	rootgprintf("******************************************************************************\n");
	pararootfprintf(logfile, "******************************************************************************\n");

#ifdef USE_MPI
	//MPI: The log file is written both in parallel in PrintParaFileOutput before this where details of interactions between stars etc are printed out, as well as here by the root node where the summary of the timestep is printed out.
    mpi_para_file_write(mpi_logfile_wrbuf, &mpi_logfile_len, &mpi_logfile_ofst_total, &mpi_logfile);
#endif

}

/**
* @brief prints out some more files
*/
void PrintFileOutput(void) {
	long i, j, n_single, n_binary, n_single_c, n_binary_c, n_single_nb, n_binary_nb, N_core_binary, N_core_binary_nb, n_10=1, n_sing_10=0, n_bin_10=0;
	double fb, fb_core, fb_core_nb, m_sing_10=0.0, m_bin_10=0.0, m_10=0.0, r_10=0.0, rho_10=0.0;
	int *multimassr_empty = (int *) malloc((NO_MASS_BINS-1)*sizeof(int));

	/* print data */
	rootfprintf(lagradfile, "%.9e ", TotalTime);
	rootfprintf(ave_mass_file, "%.9e ",TotalTime);
	rootfprintf(no_star_file, "%.9e ",TotalTime);
	rootfprintf(densities_file, "%.9e ",TotalTime);
	rootfprintf(ke_rad_file, "%.9e ",TotalTime);
	rootfprintf(ke_tan_file, "%.9e ",TotalTime);
	rootfprintf(v2_rad_file, "%.9e ",TotalTime);
	rootfprintf(v2_tan_file, "%.9e ",TotalTime);
	for(i=0; i<NO_MASS_BINS-1; i++){
		multimassr_empty[i] = 1;
		for(j=0; j<NO_MASS_BINS-1; j++){
			if (multi_mass_r[j][i] > 0.0) multimassr_empty[i] = 0;
		}
	}
	for(i=0; i<NO_MASS_BINS-1; i++){
		if ( !multimassr_empty[i] ){
			rootfprintf(mlagradfile[i], "%.9e ", TotalTime);
		}
	}

	for (i = 0; i < MASS_PC_COUNT ; i++) {
		rootfprintf(lagradfile, "%e ", mass_r[i]);
		rootfprintf(ave_mass_file,"%e ", ave_mass_r[i] * units.m / MSUN);
		rootfprintf(no_star_file,"%g ", no_star_r[i]);
		rootfprintf(densities_file,"%e ", densities_r[i]);
		rootfprintf(ke_rad_file,"%e ", ke_rad_r[i]);
		rootfprintf(ke_tan_file,"%e ", ke_tan_r[i]);
		rootfprintf(v2_rad_file,"%e ", v2_rad_r[i]);
		rootfprintf(v2_tan_file,"%e ", v2_tan_r[i]);
		for(j=0; j<NO_MASS_BINS-1; j++){
			if ( !multimassr_empty[j] ){
				rootfprintf(mlagradfile[j], "%e ", multi_mass_r[j][i]);
			}
		}
	}
	rootfprintf(lagradfile, "\n");
	rootfprintf(ave_mass_file,"\n");
	rootfprintf(no_star_file,"\n");
	rootfprintf(densities_file,"\n");
	rootfprintf(ke_rad_file,"\n");
	rootfprintf(ke_tan_file,"\n");
	rootfprintf(v2_rad_file,"\n");
	rootfprintf(v2_tan_file,"\n");
	for(i=0; i<NO_MASS_BINS-1; i++){
		if ( !multimassr_empty[i] ){
			rootfprintf(mlagradfile[i], "\n");
		}
	}

	/* information on the central BH */
	rootfprintf(centmass_file, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
		TotalTime, cenma.m * madhoc, Dt, rho_core, Etotal.tot, Etotal.K, Etotal.P);
	
	/* output Time,N_MAX,TotalE,TotalKE,TotalPE,Mtotal */
	rootfprintf(dynfile, "%.8g %.8g %ld %ld %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
		TotalTime, Dt, tcount, clus.N_MAX, Mtotal, -2.0*Etotal.K/Etotal.P, N_core, core_radius, max_r, 
		Etotal.tot, Etotal.K, Etotal.P, Etotal.Eint, Etotal.Eb, cenma.E, Eescaped, Ebescaped, Eintescaped, 
		Eoops, Etotal.tot+Eoops, clusdyn.rh, central.rho, central.rc_spitzer, central.v_rms, rc_nb, DMse*units.m/MSUN, DMrejuv*units.m/MSUN, N_core_nb);


	//Sourav: printing properties at 10% lagrange radius
	if (CALCULATE10){
#ifdef USE_MPI
        int i;
        double m_10_prev=0;

        for(i=1; i<Start[myid]; i++)
            m_10_prev += star_m[i] / clus.N_STAR;

		n_10=1;
		m_bin_10 = 0.0;
		m_sing_10 = 0.0;
		m_10 = m_10_prev;

		while (m_10 < 0.1 * Mtotal && n_10 <= clus.N_MAX_NEW) {
            int g_n_10 = get_global_idx(n_10);
			m_10 += star_m[g_n_10] / clus.N_STAR;

			if (star[n_10].binind>0){
				n_bin_10++;
				m_bin_10 += star_m[g_n_10] / clus.N_STAR;
			}
			else{
				n_sing_10++;
				m_sing_10 += star_m[g_n_10] / clus.N_STAR;
			}
			n_10++;
		}

        n_10--; //MPI: since n_10 is initialized to 1 for all procs. So if summed, it'll give the wrong index.
        long buf_comm_long[3];
        long buf_comm_long_recv[3];
        double buf_comm_dbl[3];
        double buf_comm_dbl_recv[3];

        buf_comm_long[0] = n_sing_10;
        buf_comm_long[1] = n_bin_10;
        buf_comm_long[2] = n_10;

        buf_comm_dbl[0] = m_sing_10;
        buf_comm_dbl[1] = m_bin_10;
        buf_comm_dbl[2] = m_10 - m_10_prev; //MPI: bugfix: subrtaction is reqd, otherwise m_10 of proc 0 will be added proc times, of proc 1 added proc-1 times and so on.

        //MPI: Since only root node is printing Allreduce is not reqd, just Reduce will do.
		  double tmpTimeStart = timeStartSimple();
        MPI_Reduce(buf_comm_long, buf_comm_long_recv, 3, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(buf_comm_dbl, buf_comm_dbl_recv, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		  timeEndSimple(tmpTimeStart, &t_comm);

        n_sing_10 = buf_comm_long_recv[0];
        n_bin_10 = buf_comm_long_recv[1];
        n_10 = buf_comm_long_recv[2];

        m_sing_10 = buf_comm_dbl_recv[0];
        m_bin_10 = buf_comm_dbl_recv[1];
        m_10 = buf_comm_dbl_recv[2];

		/* exit if not enough stars */
        if(myid==0)
        {
            if (n_10 <= 6 || n_10 >= clus.N_STAR-6) {
                eprintf("clus.N_STAR <= 2*J || n_10 >= clus.N_STAR-6\n");
                exit_cleanly(-1, __FUNCTION__);
            }
            else{
                r_10=star_r[n_10];
                rho_10 = m_10/(4.0 * 3.0 * PI * fb_cub(r_10));	
            }
        }

        rootfprintf(lagrad10file, "%.8g %.8g %.ld %ld %.8g %ld %.8g %ld %.8g %.8g %.8g\n",
                TotalTime, Dt, tcount, n_10, m_10, n_sing_10, m_sing_10, n_bin_10, m_bin_10, r_10, rho_10);

#else
		n_10=1;
		m_bin_10 = 0.0;
		m_sing_10 = 0.0;
		m_10 = 0.0;
		while (m_10 < 0.1 * Mtotal) {
			m_10 += star[n_10].m / clus.N_STAR;
			if (star[n_10].binind>0){
				n_bin_10++;
				m_bin_10 += star[n_10].m / clus.N_STAR;
			}
			else{
				n_sing_10++;
				m_sing_10 += star[n_10].m / clus.N_STAR;
			}
			n_10++;
		}
		/* exit if not enough stars */
		if (n_10 <= 6 || n_10 >= clus.N_STAR-6) {
			eprintf("clus.N_STAR <= 2*J || n_10 >= clus.N_STAR-6\n");
			exit_cleanly(-1, __FUNCTION__);
		}
		else{
			r_10=star[n_10].r;
			rho_10 = m_10/(4.0 * 3.0 * PI * fb_cub(r_10));	
		}
		fprintf(lagrad10file, "%.8g %.8g %.ld %ld %.8g %ld %.8g %ld %.8g %.8g %.8g\n",
					TotalTime, Dt, tcount, n_10, m_10, n_sing_10, m_sing_10, n_bin_10, m_bin_10, r_10, rho_10);
#endif
	}

	/* Output binary data Note: N_BINARY counts ALL binaries (including escaped/destroyed ones)
	   whereas N_b only counts EXISTING BOUND binaries. */
	/* calculate core binary fraction */

	n_single = 0;
	n_binary = 0;
	n_single_c = 0;
	n_binary_c = 0;
	//Sourav:initialize nb core properties
	n_single_nb = 0;
	n_binary_nb = 0;

    find_nstars_within_r(core_radius, &n_single_c, &n_binary_c);
	//Sourav:calculate n_sing and n_bin for nb core
    find_nstars_within_r(rc_nb, &n_single_nb, &n_binary_nb);

	// calculate overall binary fraction
#ifdef USE_MPI
	for (i=1; i<=clus.N_MAX_NEW; i++) {
#else
	for (i=1; i<=clus.N_MAX; i++) {
#endif
		if (star[i].binind > 0) {
			n_binary++;
		} else {
			n_single++;
		}
	}

#ifdef USE_MPI
    long buf_comm[6];
    long buf_comm_recv[6];

    buf_comm[0] = n_single_c;
    buf_comm[1] = n_binary_c;
    buf_comm[2] = n_single_nb;
    buf_comm[3] = n_binary_nb;
    buf_comm[4] = n_single;
    buf_comm[5] = n_binary;

    //MPI: Since only root node is printing Allreduce is not reqd.
	 double tmpTimeStart = timeStartSimple();
    MPI_Reduce(buf_comm, buf_comm_recv, 6, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);

    if(myid==0)
    {
        n_single_c = buf_comm_recv[0];
        n_binary_c = buf_comm_recv[1];
        n_single_nb = buf_comm_recv[2];
        n_binary_nb = buf_comm_recv[3];
        n_single = buf_comm_recv[4];
        n_binary = buf_comm_recv[5];
    }
#endif

	N_core_binary = n_binary_c;
	N_core_binary_nb = n_binary_nb;
	
	// this is such a kludge: core_radius is not initialized on the first timestep
	if (n_single_c + n_binary_c == 0) {
		fb_core = 0.0;
	} else {
		fb_core = ((double) n_binary_c)/((double) (n_single_c + n_binary_c));
	}
	
	//calculate the same for nb core
	if (n_single_nb + n_binary_nb == 0) {
		fb_core_nb = 0.0;
	} else {
		fb_core_nb = ((double) n_binary_nb)/((double) (n_single_nb + n_binary_nb));
	}
	
	// this is such a kludge: core_radius is not initialized on the first timestep
	if (n_single + n_binary == 0) {
		fb = 0.0;
	} else {
		fb = ((double) n_binary)/((double) (n_single + n_binary));
	}
	
	// print data
	rootfprintf(binaryfile,
		"%.6g %ld %.6g %.6g %.6g %.6g %.6g %.6g %ld %ld %.6g %.6g %.6g %.6g %.6g %.6g %ld %.8g %ld\n",
		TotalTime, N_b, M_b, E_b, rh_single, 
		rh_binary, rho_core_single, rho_core_bin, 
		N_bb, N_bs, fb_core, 
		fb, E_bb, E_bs, DE_bb, DE_bs, 
		N_core_binary_nb, fb_core_nb, N_core_binary);

    if (WRITE_EXTRA_CORE_INFO) {
        write_core_data(corefile, no_remnants);
        rootfprintf(corefile, "\n");
    }

	/* also saves INITIAL snapshot (StepCount=0) */
	if (TotalTime >= SNAPSHOT_DELTAT * StepCount) {
		StepCount++;
		print_2Dsnapshot();
		if (WRITE_STELLAR_INFO){
			write_stellar_data();	
		}
	}
	print_snapshot_windows();

	free(multimassr_empty);

	/* Meagan - extra output for bhs */
	if (WRITE_BH_INFO) {
		print_bh_summary();
		print_esc_bh_summary();
	}
}

/**
* @brief Meagan: extra output for bhs
*/
void print_bh_summary() {

#ifdef USE_MPI
    long buf_comm[12];
    long buf_comm_recv[12];
    buf_comm[0] = bhsingle;
    buf_comm[1] = bhbinary;
    buf_comm[2] = bhbh;
    buf_comm[3] = bhnonbh;
    buf_comm[4] = bhwd;
    buf_comm[5] = bh01;
    buf_comm[6] = bh7;
    buf_comm[7] = bh10;
    buf_comm[8] = bh12;
    buf_comm[9] = bh13;
    buf_comm[10] = bh26;
    buf_comm[11] = bh89;

    //MPI:OPT: bhstar, and bhwd might be calculated after the reduce instead of in bh_count to save 2 communication calls.
	 double tmpTimeStart = timeStartSimple();
    MPI_Reduce(buf_comm, buf_comm_recv, 12, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);

    bhsingle = buf_comm_recv[0];
    bhbinary = buf_comm_recv[1];
    bhbh = buf_comm_recv[2];
    bhnonbh  = buf_comm_recv[3];
    bhwd = buf_comm_recv[4];
    bh01 = buf_comm_recv[5];
    bh7 = buf_comm_recv[6];
    bh10 = buf_comm_recv[7];
    bh12 = buf_comm_recv[8];
    bh13 = buf_comm_recv[9];
    bh26 = buf_comm_recv[10];
    bh89 = buf_comm_recv[11];
#endif

	double fb_bh;	
	if ((bhbinary + bhsingle) > 0) {
		fb_bh = ((double) (bhbinary))/((double) (bhsingle + bhbinary));
	} else {
		fb_bh = 0.0;
	}
	rootfprintf(bhsummaryfile, "%ld %.8g %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %.6g\n",
		tcount, TotalTime, (bhsingle + bhnonbh + 2*bhbh), bhsingle, bhbinary, 
		bhbh, bhnonbh, bh13, bhwd, (bhnonbh-bhwd-bh13), bh01, (bhnonbh-bhwd-bh13-bh01), fb_bh);

	// reset all counts for next timestep
	bhbinary=0;
	bhsingle=0;
	bhbh=0;
	bhnonbh=0;
	bh13=0;
	bh10=0;
	bh11=0;
	bh12=0;
	bhwd=0;
	bhstar=0;
	bh01=0;
	bh26=0;
	bh7=0;
	bh89=0;
}

/**
* @brief Meagan - extra output for bhs
*/
void print_esc_bh_summary() {
    // Meagan: log info about escaped bhs

#ifdef USE_MPI
    long buf_comm[12];
    long buf_comm_recv[12];
    buf_comm[0] = esc_bhsingle;
    buf_comm[1] = esc_bhbinary;
    buf_comm[2] = esc_bhbh;
    buf_comm[3] = esc_bhnonbh;
    buf_comm[4] = esc_bhwd;
    buf_comm[5] = esc_bh01;
    buf_comm[6] = esc_bh7;
    buf_comm[7] = esc_bh10;
    buf_comm[8] = esc_bh12;
    buf_comm[9] = esc_bh13;
    buf_comm[10] = esc_bh26;
    buf_comm[11] = esc_bh89;

    //MPI: esc_bhstar, and esc_bhwd might be calculated after the reduce instead of in bh_count to save 2 communication calls.
	 double tmpTimeStart = timeStartSimple();
    MPI_Reduce(buf_comm, buf_comm_recv, 12, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	 timeEndSimple(tmpTimeStart, &t_comm);

    esc_bhsingle = buf_comm_recv[0];
    esc_bhbinary = buf_comm_recv[1];
    esc_bhbh = buf_comm_recv[2];
    esc_bhnonbh  = buf_comm_recv[3];
    esc_bhwd = buf_comm_recv[4];
    esc_bh01 = buf_comm_recv[5];
    esc_bh7 = buf_comm_recv[6];
    esc_bh10 = buf_comm_recv[7];
    esc_bh12 = buf_comm_recv[8];
    esc_bh13 = buf_comm_recv[9];
    esc_bh26 = buf_comm_recv[10];
    esc_bh89 = buf_comm_recv[11];
#endif

    esc_bhsingle_tot += esc_bhsingle;
    esc_bhbinary_tot += esc_bhbinary;
    esc_bhbh_tot += esc_bhbh;
    esc_bhnonbh_tot += esc_bhnonbh;
    esc_bh13_tot += esc_bh13;
    esc_bhwd_tot += esc_bhwd;
    esc_bhstar_tot = esc_bhnonbh_tot - esc_bhwd_tot - esc_bh13_tot;
    esc_bh10_tot += esc_bh10;
    esc_bh11_tot += esc_bh11;
    esc_bh12_tot += esc_bh12;
    esc_bh01_tot += esc_bh01;
    esc_bh26_tot += esc_bh26;
    esc_bh7_tot += esc_bh7;
    esc_bh89_tot += esc_bh89;
    /*if ((esc_bhbinary + esc_bhsingle)>0) {
        esc_fb_bh = ((double) (esc_bhbinary))/((double) (esc_bhbinary + esc_bhsingle));
    } else {
        esc_fb_bh = 0.0;
    }*/
    if (esc_bhbinary_tot>0 || esc_bhsingle_tot>0) {
        esc_fb_bh_tot = ((double) (esc_bhbinary_tot))/((double) (esc_bhbinary_tot + esc_bhsingle_tot));
    } else {
        esc_fb_bh_tot = 0.0;
    }
    rootfprintf(escbhsummaryfile, "%ld %.8g %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %.6g\n",
            tcount, TotalTime, 
            // counts for this timestep
            (esc_bhsingle_tot + esc_bhnonbh_tot + 2*esc_bhbh_tot), esc_bhsingle_tot, 
            esc_bhbinary_tot, esc_bhbh_tot, esc_bhnonbh_tot, esc_bh13_tot, esc_bhwd_tot,
	    esc_bhstar_tot, esc_bh01_tot, (esc_bhstar_tot-esc_bh01_tot),esc_fb_bh_tot);
            // cumulative counts for escaping bhs
            /*(esc_bhsingle_tot + esc_bhbinary_tot + 2*esc_bhbh_tot), esc_bhsingle_tot, 
            esc_bhbinary_tot, esc_bhbh_tot, esc_bhnonbh_tot, esc_bh13_tot, esc_bhwd_tot, 
            (esc_bh10_tot+esc_bh11_tot+esc_bh12_tot), esc_bhstar_tot, esc_bh01_tot, 
            esc_bh26_tot, esc_bh7_tot, esc_bh89_tot, esc_fb_bh_tot);*/

    // reset counts for next timestep
    esc_bhsingle = 0;
    esc_bhbinary = 0;
    esc_bhbh = 0;
    esc_bhnonbh = 0;
    esc_bh13 = 0;
    esc_bhwd = 0;
    esc_bh10 = 0;
    esc_bh11 = 0;
    esc_bh12 = 0;
    esc_bhstar = 0;
    esc_bh01 = 0;
    esc_bh26 = 0;
    esc_bh7 = 0;
    esc_bh89 = 0;
    esc_fb_bh = 0.0;
}

/**
* @brief finds number of single and binary stars within a given radial position
*
* @param r radial position
* @param ns variable to hold/return number of single stars
* @param nb variable to hold/return number of binary stars
*/
void find_nstars_within_r(double r, long *ns, long *nb)
{
    long i;
#ifdef USE_MPI
	for (i=1; star_r[get_global_idx(i)]<=r && i<=clus.N_MAX_NEW; i++) {
#else
	for (i=1; star[i].r<=r; i++) {
#endif
		if (star[i].binind > 0) {
			(*nb)++;
		} else {
			(*ns)++;
		}
	}
}

#ifdef USE_MPI
/**
* @brief This is the function which actually writes the parallel file buffers into the file in parallel using MPI IO.
*/
void PrintParaFileOutput(void)
{
	//This macro writes out the corresponding buffer into the corresponding file in parallel using MPI-IO. Here we write out all the files that need contribution from more than one processor.
    mpi_para_file_write(mpi_logfile_wrbuf, &mpi_logfile_len, &mpi_logfile_ofst_total, &mpi_logfile);
    mpi_para_file_write(mpi_escfile_wrbuf, &mpi_escfile_len, &mpi_escfile_ofst_total, &mpi_escfile);
    mpi_para_file_write(mpi_binintfile_wrbuf, &mpi_binintfile_len, &mpi_binintfile_ofst_total, &mpi_binintfile);
    mpi_para_file_write(mpi_collisionfile_wrbuf, &mpi_collisionfile_len, &mpi_collisionfile_ofst_total, &mpi_collisionfile);
    mpi_para_file_write(mpi_tidalcapturefile_wrbuf, &mpi_tidalcapturefile_len, &mpi_tidalcapturefile_ofst_total, &mpi_tidalcapturefile);
    mpi_para_file_write(mpi_semergedisruptfile_wrbuf, &mpi_semergedisruptfile_len, &mpi_semergedisruptfile_ofst_total, &mpi_semergedisruptfile);
    mpi_para_file_write(mpi_removestarfile_wrbuf, &mpi_removestarfile_len, &mpi_removestarfile_ofst_total, &mpi_removestarfile);
    mpi_para_file_write(mpi_relaxationfile_wrbuf, &mpi_relaxationfile_len, &mpi_relaxationfile_ofst_total, &mpi_relaxationfile);

	 if(WRITE_PULSAR_INFO)
		 mpi_para_file_write(mpi_pulsarfile_wrbuf, &mpi_pulsarfile_len, &mpi_pulsarfile_ofst_total, &mpi_pulsarfile);

    /* Meagan's 3bb files */
    if (WRITE_BH_INFO){
        mpi_para_file_write(mpi_newbhfile_wrbuf, &mpi_newbhfile_len, &mpi_newbhfile_ofst_total, &mpi_newbhfile);
        mpi_para_file_write(mpi_bhmergerfile_wrbuf, &mpi_bhmergerfile_len, &mpi_bhmergerfile_ofst_total, &mpi_bhmergerfile);
    }

    if (THREEBODYBINARIES)
    {
        mpi_para_file_write(mpi_threebbfile_wrbuf, &mpi_threebbfile_len, &mpi_threebbfile_ofst_total, &mpi_threebbfile);
        mpi_para_file_write(mpi_threebbprobabilityfile_wrbuf, &mpi_threebbprobabilityfile_len, &mpi_threebbprobabilityfile_ofst_total, &mpi_threebbprobabilityfile);
        mpi_para_file_write(mpi_lightcollisionfile_wrbuf, &mpi_lightcollisionfile_len, &mpi_lightcollisionfile_ofst_total, &mpi_lightcollisionfile);
        mpi_para_file_write(mpi_threebbdebugfile_wrbuf, &mpi_threebbdebugfile_len, &mpi_threebbdebugfile_ofst_total, &mpi_threebbdebugfile);
    }
}

/**
* @brief Flushes out data in parallel present in the char buffer to the corresponding file using MPI-IO
*
* @param wrbuf write buffer containing the data to be flushed out
* @param len length/size of the data in the buffer
* @param prev_cum_offset offset of the file where the data needs to be written
* @param fh MPI-IO File handle
*/
void mpi_para_file_write(char* wrbuf, long long* len, long long* prev_cum_offset, MPI_File* fh)
{
    MPI_Offset mpi_offset=0;
	 long long offset=0;
	 long long tot_offset=0;
    MPI_Status mpistat;

	 //First find out the offset for this processor based on the buffer lengths of other procecessors
    MPI_Exscan(len, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	 //Add this offset to the previous cumulative file offset
    offset += *prev_cum_offset;
	 mpi_offset = offset;

	 //Write data to file in parallel
    MPI_File_write_at_all(*fh, mpi_offset, wrbuf, *len, MPI_CHAR, &mpistat);

	 //Update cumulative file offset for next flush
    MPI_Allreduce (len, &tot_offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    *prev_cum_offset += tot_offset;

	 //Reset buffer and length variables
    wrbuf[0] = '\0';
    *len=0;
}
#endif

/**
* @brief Parsing of Input Parameters / Memory allocation / File I/O
*
* @param argc Input argument count from main()
* @param argv[] Input argument list from main()
* @param r gsl rng
*/
void parser(int argc, char *argv[], gsl_rng *r)
{
	char inputfile[1024], outfile[1024], outfilemode[5];
	char parameter_name[1024], values[1024], dummy[1024], line[2048];
	char *curr_mass;
	parsed_t parsed;
	parsed_t *spp;
	long i, j;
	int allparsed=1;
	int hard_restart=0;
	/* int *ip; */
	FILE *in, *parsedfp;
	const char *short_opts = "qdVhs:R:r:";
	const struct option long_opts[] = {
		{"quiet", no_argument, NULL, 'q'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{"restart", required_argument, NULL, 'R'},
		{"streams", required_argument, NULL, 's'}, //Run with multiple random streams. To mimic the parallel version with the given number of processors
		{NULL, 0, NULL, 0}
	};
	
	/* DEFAULT PARAMETER VALUES */
	quiet = 0;
	debug = 0;
	NO_MASS_BINS = 0;
	snapshot_window_count=0;
	/* DEFAULT PARAMETER VALUES */
	
	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'q':
			quiet = 1;
			break;
		case 'd':
			debug = 1;
			break;
		case 'V':
			print_version(stdout);
			exit(0);
		case 's':
			procs = atoi(optarg);
			break;
		case 'R':
			RESTART_TCOUNT = atol(optarg);
			hard_restart = 1;
			break;
		case 'r':
			RESTART_TCOUNT = atol(optarg);
			hard_restart = 0;
			break;
		case 'h':
			print_version(stdout);
			fprintf(stdout, "\n");
			cmc_print_usage(stdout, argv);
			exit(0);
		default:
			break;
		}
	}

	/* check to make sure there was nothing crazy on the command line */
	if (argc - optind != 2 + hard_restart) {
		cmc_print_usage(stdout, argv);
		exit(0);
	}

	/* set inputfile and outprefix now that the options have been parsed */
	sprintf(inputfile, "%s", argv[optind]);
	sprintf(outprefix, "%s", argv[optind+1]);
	if(hard_restart)
		sprintf(oldoutprefix, "%s", argv[optind+2]);
	else
		sprintf(oldoutprefix, "%s", argv[optind+1]);
		
/*
#ifdef USE_MPI
	strcpy(outprefix_bak, outprefix);
	sprintf(outprefix, "%s%d", outprefix, myid);
#endif
*/
	dprintf("inputfile=%s outprefix=%s\n", inputfile, outprefix);

	/*======= Opening of input & output files =======*/
	if ((in = fopen(inputfile, "r")) == NULL) {
		eprintf("Cannot open input file \"%s\".\n", inputfile);
		exit(1);
	}
	
//MPI: File printed out only by the root node.
#ifdef USE_MPI
if(myid==0) {
#endif
	sprintf(outfile, "%s.cmc.parsed", outprefix);
	if ((parsedfp = fopen(outfile, "w")) == NULL) {
		eprintf("cannot create output file \"%s\".\n", outfile);
		exit(1);
	}
#ifdef USE_MPI
}
#endif

	/* nothing is set yet, so assigning all zeros to variable parsed */
	/* one way to do it is a complicated for loop (?which depends on the 
	 * fact that all members of the struct is int?), this is the fun way */
	//ip = (int *) &parsed;
	//for(i=0; i<(sizeof(parsed)/sizeof(int)); i++){
	//	*(ip++) = 0;
	//}
	/* faster and portable way to do this is using calloc, this does
	 * not depend on the structure of the struct */
	spp = (parsed_t *) calloc(1, sizeof(parsed_t));
	parsed = *spp;
	free(spp);
	
	while (fgets(line, 2048, in)) {
		/* the maximum size of an input line is limited by 
		 * the size of the line array */
		if (strlen(line) == 2047) {
			eprintf("input line too long: \"%s\".\n", line);
			exit(1);
		}
		/* strip comments and newline character */
		line[strcspn(line, "#\n")] = '\0';
		
		/* replace equal sign with a space */
		while (strchr(line, '=') != NULL) {
			line[strcspn(line, "=")] = ' ';
		}
		
		/* convert lowercase letters to uppercase in the keyword only */
		i = 0;
		while (!(isspace(line[i]))) {
			line[i] = toupper(line[i]);
			i++;
		}

#ifdef USE_MPI
#define PRINT_PARSED(DOC) if(myid==0) fprintf(parsedfp, "# %s\n%s\n", DOC, line)
#else
#define PRINT_PARSED(DOC) fprintf(parsedfp, "# %s\n%s\n", DOC, line)
#endif

		/* see if there are too many values for parameter */
		if (sscanf(line, "%s %s %s", parameter_name, values, dummy) == 3) {
			eprintf("too many values for parameter: \"%s\".\n", line);
			exit(1);
		} else if (sscanf(line, "%s %s", parameter_name, values) == 2) {
			if (strcmp(parameter_name, "SAMPLESIZE") == 0) {
				PRINT_PARSED(PARAMDOC_SAMPLESIZE);
				sscanf(values, "%d", &SAMPLESIZE);
				parsed.SAMPLESIZE = 1;
			} else if (strcmp(parameter_name, "BINBIN") == 0) {
				PRINT_PARSED(PARAMDOC_BINBIN);
				sscanf(values, "%d", &BINBIN);
				parsed.BINBIN = 1;
			} else if (strcmp(parameter_name, "BINSINGLE") == 0) {
				PRINT_PARSED(PARAMDOC_BINSINGLE);
				sscanf(values, "%d", &BINSINGLE);
				parsed.BINSINGLE = 1;
			} else if (strcmp(parameter_name, "STREAMS") == 0) {
				PRINT_PARSED(PARAMDOC_STREAMS);
				sscanf(values, "%d", &procs);
				parsed.STREAMS = 1;
	/* Meagan: new parameters for three-body binaries: THREEBODYBINARIES, MIN_BINARY_HARDNESS, ONLY_FORM_BH_THREEBODYBINARIES */
			} else if (strcmp(parameter_name, "THREEBODYBINARIES") == 0) {
				PRINT_PARSED(PARAMDOC_THREEBODYBINARIES);
				sscanf(values, "%d", &THREEBODYBINARIES);
				parsed.THREEBODYBINARIES = 1;
			} else if (strcmp(parameter_name, "MIN_BINARY_HARDNESS") == 0) {
				PRINT_PARSED(PARAMDOC_MIN_BINARY_HARDNESS);
				sscanf(values, "%lf", &MIN_BINARY_HARDNESS);
				parsed.MIN_BINARY_HARDNESS = 1;
			} else if (strcmp(parameter_name, "ONLY_FORM_BH_THREEBODYBINARIES") == 0) {
				PRINT_PARSED(PARAMDOC_ONLY_FORM_BH_THREEBODYBINARIES);
				sscanf(values, "%d", &ONLY_FORM_BH_THREEBODYBINARIES);
				parsed.ONLY_FORM_BH_THREEBODYBINARIES = 1;
			} else if (strcmp(parameter_name, "BH_SNAPSHOTTING") == 0) {
				PRINT_PARSED(PARAMDOC_BH_SNAPSHOTTING);
				sscanf(values, "%d", &BH_SNAPSHOTTING);
				parsed.BH_SNAPSHOTTING = 1;
            } else if (strcmp(parameter_name, "BH_SNAPSHOT_DELTACOUNT") == 0) {
				PRINT_PARSED(PARAMDOC_BH_SNAPSHOT_DELTACOUNT);
				sscanf(values, "%ld", &BH_SNAPSHOT_DELTACOUNT);
				parsed.BH_SNAPSHOT_DELTACOUNT = 1;
			} else if (strcmp(parameter_name, "SNAPSHOTTING") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOTTING);
				sscanf(values, "%d", &SNAPSHOTTING);
				parsed.SNAPSHOTTING = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_DELTACOUNT") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_DELTACOUNT);
				sscanf(values, "%ld", &SNAPSHOT_DELTACOUNT);
				parsed.SNAPSHOT_DELTACOUNT = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_DELTAT") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_DELTAT);
				sscanf(values, "%lf", &SNAPSHOT_DELTAT);
				parsed.SNAPSHOT_DELTAT = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_CORE_BOUNCE") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_CORE_BOUNCE);
				sscanf(values, "%d", &SNAPSHOT_CORE_BOUNCE);
				parsed.SNAPSHOT_CORE_BOUNCE = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_CORE_COLLAPSE") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_CORE_COLLAPSE);
				sscanf(values, "%d", &SNAPSHOT_CORE_COLLAPSE);
				parsed.SNAPSHOT_CORE_COLLAPSE = 1;
            } else if (strcmp(parameter_name, "SNAPSHOT_WINDOWS") == 0) {
                PRINT_PARSED(PARAMDOC_SNAPSHOT_WINDOWS);
                if (strncmp(values, "NULL", 4) == 0) {
                    SNAPSHOT_WINDOWS=NULL;
                } else {
                    SNAPSHOT_WINDOWS= (char *) malloc(sizeof(char)*300);
                    strncpy(SNAPSHOT_WINDOWS, values, 300);
                }
				parsed.SNAPSHOT_WINDOWS = 1;
			} else if (strcmp(parameter_name, "SNAPSHOT_WINDOW_UNITS") == 0) {
				PRINT_PARSED(PARAMDOC_SNAPSHOT_WINDOW_UNITS);
				SNAPSHOT_WINDOW_UNITS= (char *) malloc(sizeof(char)*10);
				strncpy(SNAPSHOT_WINDOW_UNITS,values, 10);
				if (!valid_snapshot_window_units()) {
					eprintf("Unrecognized snapshot window time unit %s.", values);
					free_arrays();
					exit(-1);
				}
				parsed.SNAPSHOT_WINDOW_UNITS = 1;
			} else if (strcmp(parameter_name, "IDUM") == 0) {
				PRINT_PARSED(PARAMDOC_IDUM);
				sscanf(values, "%ld", &IDUM);
				parsed.IDUM = 1;
			} else if (strcmp(parameter_name, "INPUT_FILE") == 0) {
				PRINT_PARSED(PARAMDOC_INPUT_FILE);
				sscanf(values, "%s", INPUT_FILE);
				parsed.INPUT_FILE = 1;
			} else if (strcmp(parameter_name, "BINARY_DISTANCE_BREAKING") == 0) {
				PRINT_PARSED(PARAMDOC_BINARY_DISTANCE_BREAKING);
				sscanf(values, "%lf", &BINARY_DISTANCE_BREAKING);
				parsed.BINARY_DISTANCE_BREAKING = 1;
			} else if (strcmp(parameter_name, "BINARY_BREAKING_MIN") == 0) {
				PRINT_PARSED(PARAMDOC_BINARY_BREAKING_MIN);
				sscanf(values, "%d", &BINARY_BREAKING_MIN);
				parsed.BINARY_BREAKING_MIN = 1;
			} else if (strcmp(parameter_name, "MASS_PC_BH_INCLUDE") == 0) {
				PRINT_PARSED(PARAMDOC_MASS_PC_BH_INCLUDE);
				sscanf(values, "%d", &MASS_PC_BH_INCLUDE);
				parsed.MASS_PC_BH_INCLUDE = 1;
			} else if (strcmp(parameter_name, "MASS_PC") == 0) {
				PRINT_PARSED(PARAMDOC_MASS_PC);
				strcpy(MASS_PC, values);
				curr_mass = (char *) strtok(values, ",; ");
				for (MASS_PC_COUNT = 1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; MASS_PC_COUNT++);
				parsed.MASS_PC = 1;
			} else if (strcmp(parameter_name, "MASS_BINS") == 0) {
				PRINT_PARSED(PARAMDOC_MASS_BINS);
				/* we recycle variable "curr_mass" for mass bins */
				strcpy(MASS_BINS, values);
				curr_mass = (char *) strtok(values, ",; ");
				for (NO_MASS_BINS = 1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; NO_MASS_BINS++);
				parsed.MASS_BINS = 1;
			} else if (strcmp(parameter_name, "MINIMUM_R") == 0) {
				PRINT_PARSED(PARAMDOC_MINIMUM_R);
				sscanf(values, "%lf", &MINIMUM_R);
				parsed.MINIMUM_R = 1;
			} else if (strcmp(parameter_name, "STOPATCORECOLLAPSE") == 0) {
				PRINT_PARSED(PARAMDOC_STOPATCORECOLLAPSE);
				sscanf(values, "%d", &STOPATCORECOLLAPSE);
				parsed.STOPATCORECOLLAPSE = 1;
			} else if (strcmp(parameter_name, "NUM_CENTRAL_STARS") == 0) {
				PRINT_PARSED(PARAMDOC_NUM_CENTRAL_STARS);
				sscanf(values, "%ld", &NUM_CENTRAL_STARS);
				parsed.NUM_CENTRAL_STARS = 1;
			} else if (strcmp(parameter_name, "PERTURB") == 0) {
				PRINT_PARSED(PARAMDOC_PERTURB);
				sscanf(values, "%ld", &PERTURB);
				parsed.PERTURB = 1;
			} else if (strcmp(parameter_name, "RELAXATION") == 0) {
				PRINT_PARSED(PARAMDOC_RELAXATION);
				sscanf(values, "%ld", &RELAXATION);
				parsed.RELAXATION = 1;
			} else if (strcmp(parameter_name, "TIDALLY_STRIP_STARS") == 0) {
				PRINT_PARSED(PARAMDOC_TIDALLY_STRIP_STARS);
				sscanf(values, "%ld", &TIDALLY_STRIP_STARS);
				parsed.TIDALLY_STRIP_STARS = 1;
			} else if (strcmp(parameter_name, "THETASEMAX") == 0) {
				PRINT_PARSED(PARAMDOC_THETASEMAX);
				sscanf(values, "%lf", &THETASEMAX);
				parsed.THETASEMAX = 1;
			} else if (strcmp(parameter_name, "STELLAR_EVOLUTION") == 0) {
				PRINT_PARSED(PARAMDOC_STELLAR_EVOLUTION);
				sscanf(values, "%ld", &STELLAR_EVOLUTION);
				parsed.STELLAR_EVOLUTION = 1;
			} else if (strcmp(parameter_name, "TIDAL_TREATMENT") == 0) {
				PRINT_PARSED(PARAMDOC_TIDAL_TREATMENT);
				sscanf(values, "%ld", &TIDAL_TREATMENT);
				parsed.TIDAL_TREATMENT = 1;
			} else if (strcmp(parameter_name, "SS_COLLISION") == 0) {
				PRINT_PARSED(PARAMDOC_SS_COLLISION);
				sscanf(values, "%ld", &SS_COLLISION);
				parsed.SS_COLLISION = 1;
			} else if (strcmp(parameter_name, "TIDAL_CAPTURE") == 0) {
				PRINT_PARSED(PARAMDOC_TIDAL_CAPTURE);
				sscanf(values, "%ld", &TIDAL_CAPTURE);
				parsed.TIDAL_CAPTURE = 1;
			} /*Sourav:new parameter*/
			else if (strcmp(parameter_name, "STAR_AGING_SCHEME") == 0) {
			 	PRINT_PARSED(PARAMDOC_STAR_AGING_SCHEME);
				sscanf(values, "%ld", &STAR_AGING_SCHEME);
				parsed.STAR_AGING_SCHEME = 1;
			} else if  (strcmp(parameter_name, "PREAGING") == 0) {
			 	PRINT_PARSED(PARAMDOC_PREAGING);
				sscanf(values, "%ld", &PREAGING);
				parsed.PREAGING = 1;
			} else if (strcmp(parameter_name, "TERMINAL_ENERGY_DISPLACEMENT") == 0) {
				PRINT_PARSED(PARAMDOC_TERMINAL_ENERGY_DISPLACEMENT);
				sscanf(values, "%lf", &TERMINAL_ENERGY_DISPLACEMENT);
				parsed.TERMINAL_ENERGY_DISPLACEMENT = 1;
			} else if (strcmp(parameter_name, "T_MAX") == 0) {
				PRINT_PARSED(PARAMDOC_T_MAX);
				sscanf(values, "%lf", &T_MAX);
				parsed.T_MAX = 1;
			} else if (strcmp(parameter_name, "T_MAX_PHYS") == 0) {
				PRINT_PARSED(PARAMDOC_T_MAX_PHYS);
				sscanf(values, "%lf", &T_MAX_PHYS);
				parsed.T_MAX_PHYS = 1;
			} else if (strcmp(parameter_name, "T_MAX_COUNT") == 0) {
				PRINT_PARSED(PARAMDOC_T_MAX_COUNT);
				sscanf(values, "%ld", &T_MAX_COUNT);
				parsed.T_MAX_COUNT = 1;
			} else if (strcmp(parameter_name, "MAX_WCLOCK_TIME") == 0) {
				PRINT_PARSED(PARAMDOC_MAX_WCLOCK_TIME);
				sscanf(values, "%ld", &MAX_WCLOCK_TIME);
				parsed.MAX_WCLOCK_TIME = 1;
			} else if (strcmp(parameter_name, "CHECKPOINT_INTERVAL") == 0) {
				PRINT_PARSED(PARAMDOC_CHECKPOINT_INTERVAL);
				sscanf(values, "%ld", &CHECKPOINT_INTERVAL);
				parsed.CHECKPOINT_INTERVAL = 1;
			} else if (strcmp(parameter_name, "CHECKPOINTS_TO_KEEP") == 0) {
				PRINT_PARSED(PARAMDOC_CHECKPOINTS_TO_KEEP);
				sscanf(values, "%ld", &CHECKPOINTS_TO_KEEP);
				parsed.CHECKPOINTS_TO_KEEP = 1;
			} else if (strcmp(parameter_name, "WIND_FACTOR") == 0) {
				PRINT_PARSED(PARAMDOC_WIND_FACTOR);
				sscanf(values, "%lf", &WIND_FACTOR);
				parsed.WIND_FACTOR = 1;
			} else if (strcmp(parameter_name, "GAMMA") == 0) {
				PRINT_PARSED(PARAMDOC_GAMMA);
				sscanf(values, "%lf", &GAMMA);
				parsed.GAMMA = 1;
			} else if (strcmp(parameter_name, "SEARCH_GRID")== 0) {
				PRINT_PARSED(PARAMDOC_SEARCH_GRID);
				sscanf(values, "%ld", &SEARCH_GRID);
				parsed.SEARCH_GRID = 1;
			} else if (strcmp(parameter_name, "SG_STARSPERBIN")== 0) {
				PRINT_PARSED(PARAMDOC_SG_STARSPERBIN);
				sscanf(values, "%ld", &SG_STARSPERBIN);
				parsed.SG_STARSPERBIN = 1;
			} else if (strcmp(parameter_name, "SG_MAXLENGTH")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MAXLENGTH);
				sscanf(values, "%ld", &SG_MAXLENGTH);
				parsed.SG_MAXLENGTH = 1;
			} else if (strcmp(parameter_name, "SG_MINLENGTH")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MINLENGTH);
				sscanf(values, "%ld", &SG_MINLENGTH);
				parsed.SG_MINLENGTH = 1;
			} else if (strcmp(parameter_name, "SG_POWER_LAW_EXPONENT")== 0) {
				PRINT_PARSED(PARAMDOC_SG_POWER_LAW_EXPONENT);
				sscanf(values, "%lf", &SG_POWER_LAW_EXPONENT);
				parsed.SG_POWER_LAW_EXPONENT = 1;
			} else if (strcmp(parameter_name, "SG_MATCH_AT_FRACTION")== 0) {
				PRINT_PARSED(PARAMDOC_SG_MATCH_AT_FRACTION);
				sscanf(values, "%lf", &SG_MATCH_AT_FRACTION);
				parsed.SG_MATCH_AT_FRACTION = 1;
			} else if (strcmp(parameter_name, "SG_PARTICLE_FRACTION")== 0) {
				PRINT_PARSED(PARAMDOC_SG_PARTICLE_FRACTION);
				sscanf(values, "%lf", &SG_PARTICLE_FRACTION);
				parsed.SG_PARTICLE_FRACTION = 1;
			} else if (strcmp(parameter_name, "BH_LOSS_CONE")== 0) {
				PRINT_PARSED(PARAMDOC_BH_LOSS_CONE);
				sscanf(values, "%li", &BH_LOSS_CONE);
				parsed.BH_LOSS_CONE = 1;
			} else if (strcmp(parameter_name, "BH_R_DISRUPT_NB")== 0) {
				PRINT_PARSED(PARAMDOC_BH_R_DISRUPT_NB);
				sscanf(values, "%lf", &BH_R_DISRUPT_NB);
				parsed.BH_R_DISRUPT_NB = 1;
			} else if (strcmp(parameter_name, "FORCE_RLX_STEP")== 0) {
				PRINT_PARSED(PARAMDOC_FORCE_RLX_STEP);
				sscanf(values, "%i", &FORCE_RLX_STEP);
				parsed.FORCE_RLX_STEP = 1;
            } else if (strcmp(parameter_name, "DT_HARD_BINARIES")== 0) {
				PRINT_PARSED(PARAMDOC_DT_HARD_BINARIES);
				sscanf(values, "%i", &DT_HARD_BINARIES);
				parsed.DT_HARD_BINARIES = 1;
            } else if (strcmp(parameter_name, "HARD_BINARY_KT")== 0) {
				PRINT_PARSED(PARAMDOC_HARD_BINARY_KT);
				sscanf(values, "%lf", &HARD_BINARY_KT);
				parsed.HARD_BINARY_KT = 1;
#ifdef EXPERIMENTAL
			} else if (strcmp(parameter_name, "BH_LC_FDT")== 0) {
				PRINT_PARSED(PARAMDOC_BH_LC_FDT);
				sscanf(values, "%lf", &BH_LC_FDT);
				parsed.BH_LC_FDT = 1;
			} else if (strcmp(parameter_name, "AVEKERNEL")== 0) {
				PRINT_PARSED(PARAMDOC_AVEKERNEL);
				sscanf(values, "%li", &AVEKERNEL);
				parsed.AVEKERNEL = 1;
#endif
			} else if (strcmp(parameter_name, "MIN_CHUNK_SIZE")== 0) {
				PRINT_PARSED(PARAMDOC_MIN_CHUNK_SIZE);
				sscanf(values, "%li", &MIN_CHUNK_SIZE);
				parsed.MIN_CHUNK_SIZE = 1;
                        } else if (strcmp(parameter_name, "BH_AVEKERNEL")== 0) {
				PRINT_PARSED(PARAMDOC_BH_AVEKERNEL);
				sscanf(values, "%li", &BH_AVEKERNEL);
				parsed.BH_AVEKERNEL = 1;
			} else if (strcmp(parameter_name, "APSIDES_PRECISION")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_PRECISION);
				sscanf(values, "%lf", &APSIDES_PRECISION);
				parsed.APSIDES_PRECISION = 1;
			} else if (strcmp(parameter_name, "APSIDES_MAX_ITER")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_MAX_ITER);
				sscanf(values, "%li", &APSIDES_MAX_ITER);
				parsed.APSIDES_MAX_ITER = 1;
			} else if (strcmp(parameter_name, "APSIDES_CONVERGENCE")== 0) {
				PRINT_PARSED(PARAMDOC_APSIDES_CONVERGENCE);
				sscanf(values, "%lf", &APSIDES_CONVERGENCE);
				parsed.APSIDES_CONVERGENCE = 1;
			} else if (strcmp(parameter_name, "CIRC_PERIOD_THRESHOLD")== 0) {
				PRINT_PARSED(PARAMDOC_CIRC_PERIOD_THRESHOLD);
				sscanf(values, "%lf", &CIRC_PERIOD_THRESHOLD);
				parsed.CIRC_PERIOD_THRESHOLD = 1;
			} else if (strcmp(parameter_name, "WRITE_STELLAR_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_STELLAR_INFO);
				sscanf(values, "%i", &WRITE_STELLAR_INFO);
				parsed.WRITE_STELLAR_INFO = 1;
            } else if (strcmp(parameter_name, "WRITE_BH_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_BH_INFO);
				sscanf(values, "%i", &WRITE_BH_INFO);
				parsed.WRITE_BH_INFO = 1;
            } else if (strcmp(parameter_name, "WRITE_RWALK_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_RWALK_INFO);
				sscanf(values, "%i", &WRITE_RWALK_INFO);
				parsed.WRITE_RWALK_INFO = 1;
			} else if (strcmp(parameter_name, "WRITE_EXTRA_CORE_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_EXTRA_CORE_INFO);
				sscanf(values, "%i", &WRITE_EXTRA_CORE_INFO);
				parsed.WRITE_EXTRA_CORE_INFO = 1;
			} else if (strcmp(parameter_name, "WRITE_PULSAR_INFO")== 0) {
				PRINT_PARSED(PARAMDOC_WRITE_PULSAR_INFO);
				sscanf(values, "%i", &WRITE_PULSAR_INFO);
				parsed.WRITE_PULSAR_INFO = 1;
			} else if (strcmp(parameter_name, "CALCULATE10")== 0) {
				PRINT_PARSED(PARAMDOC_CALCULATE10);
				sscanf(values, "%i", &CALCULATE10);
				parsed.CALCULATE10 = 1;
			} else if (strcmp(parameter_name, "OVERWRITE_RVIR")== 0) {
				PRINT_PARSED(PARAMDOC_OVERWRITE_RVIR);
				sscanf(values, "%lf", &OVERWRITE_RVIR);
				parsed.OVERWRITE_RVIR = 1;
			} else if (strcmp(parameter_name, "OVERWRITE_Z")== 0) {
				PRINT_PARSED(PARAMDOC_OVERWRITE_Z);
				sscanf(values, "%lf", &OVERWRITE_Z);
				parsed.OVERWRITE_Z = 1;
			} else if (strcmp(parameter_name, "OVERWRITE_RTID")== 0) {
				PRINT_PARSED(PARAMDOC_OVERWRITE_RTID);
				sscanf(values, "%lf", &OVERWRITE_RTID);
				parsed.OVERWRITE_RTID = 1;
			} else if (strcmp(parameter_name, "OVERWRITE_MCLUS")== 0) {
				PRINT_PARSED(PARAMDOC_OVERWRITE_MCLUS);
				sscanf(values, "%lf", &OVERWRITE_MCLUS);
				parsed.OVERWRITE_MCLUS = 1;
				// Begin reading in stellar and binary evolution assumptions.
			} else if (strcmp(parameter_name, "BSE_NETA")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_NETA);
				sscanf(values, "%lf", &BSE_NETA);
				parsed.BSE_NETA = 1;
			} else if (strcmp(parameter_name, "BSE_BWIND")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BWIND);
				sscanf(values, "%lf", &BSE_BWIND);
				parsed.BSE_BWIND = 1;
			} else if (strcmp(parameter_name, "BSE_HEWIND")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_HEWIND);
				sscanf(values, "%lf", &BSE_HEWIND);
				parsed.BSE_HEWIND = 1;
			} else if (strcmp(parameter_name, "BSE_WINDFLAG")==0) {
				PRINT_PARSED(PARAMDOC_BSE_WINDFLAG);
				sscanf(values, "%d", &BSE_WINDFLAG);
				parsed.BSE_WINDFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_PPSN")==0) {
				PRINT_PARSED(PARAMDOC_BSE_PPSN);
				sscanf(values, "%d", &BSE_PPSN);
				parsed.BSE_PPSN = 1;
			} else if (strcmp(parameter_name, "BSE_ALPHA1")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_ALPHA1);
				sscanf(values, "%lf", &BSE_ALPHA1);
				parsed.BSE_ALPHA1 = 1;
			} else if (strcmp(parameter_name, "BSE_LAMBDA")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_LAMBDA);
				sscanf(values, "%lf", &BSE_LAMBDA);
				parsed.BSE_LAMBDA = 1;
			} else if (strcmp(parameter_name, "BSE_CEFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_CEFLAG);
				sscanf(values, "%i", &BSE_CEFLAG);
				parsed.BSE_CEFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_TFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_TFLAG);
				sscanf(values, "%i", &BSE_TFLAG);
				parsed.BSE_TFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_IFFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_IFFLAG);
				sscanf(values, "%i", &BSE_IFFLAG);
				parsed.BSE_IFFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_WDFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_WDFLAG);
				sscanf(values, "%i", &BSE_WDFLAG);
				parsed.BSE_WDFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_BHFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BHFLAG);
				sscanf(values, "%i", &BSE_BHFLAG);
				parsed.BSE_BHFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_BHSPINFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BHSPINFLAG);
				sscanf(values, "%i", &BSE_BHSPINFLAG);
				parsed.BSE_BHSPINFLAG = 1;
			} else if (strcmp(parameter_name, "BH_RADIUS_MULTIPLYER")== 0) {
				PRINT_PARSED(PARAMDOC_BH_RADIUS_MULTIPLYER);
				sscanf(values, "%i", &BH_RADIUS_MULTIPLYER);
				parsed.BH_RADIUS_MULTIPLYER = 1;
			} else if (strcmp(parameter_name, "BSE_NSFLAG")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_NSFLAG);
				sscanf(values, "%i", &BSE_NSFLAG);
				parsed.BSE_NSFLAG = 1;
			} else if (strcmp(parameter_name, "BSE_MXNS")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_MXNS);
				sscanf(values, "%lf", &BSE_MXNS);
				parsed.BSE_MXNS = 1;
			} else if (strcmp(parameter_name, "BSE_BCONST")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BCONST);
				sscanf(values, "%lf", &BSE_BCONST);
				parsed.BSE_BCONST = 1;
			} else if (strcmp(parameter_name, "BSE_CK")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_CK);
				sscanf(values, "%lf", &BSE_CK);
				parsed.BSE_CK = 1;
			} else if (strcmp(parameter_name, "BSE_IDUM")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_IDUM);
				sscanf(values, "%i", &BSE_IDUM);
				parsed.BSE_IDUM = 1;
			} else if (strcmp(parameter_name, "BSE_BETA")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BETA);
				sscanf(values, "%lf", &BSE_BETA);
				parsed.BSE_BETA = 1;
			} else if (strcmp(parameter_name, "BSE_SIGMA")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_SIGMA);
				sscanf(values, "%lf", &BSE_SIGMA);
				parsed.BSE_SIGMA = 1;
			} else if (strcmp(parameter_name, "BSE_BHSIGMAFRAC")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_BHSIGMAFRAC);
				sscanf(values, "%lf", &BSE_BHSIGMAFRAC);
				parsed.BSE_BHSIGMAFRAC = 1;
			} else if (strcmp(parameter_name, "BSE_OPENING_ANGLE")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_OPENING_ANGLE);
				sscanf(values, "%d", &BSE_OPENING_ANGLE);
				parsed.BSE_OPENING_ANGLE = 1;
			} else if (strcmp(parameter_name, "BSE_EDDFAC")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_EDDFAC);
				sscanf(values, "%lf", &BSE_EDDFAC);
				parsed.BSE_EDDFAC = 1;
			} else if (strcmp(parameter_name, "BSE_GAMMA")== 0) {
				PRINT_PARSED(PARAMDOC_BSE_GAMMA);
				sscanf(values, "%lf", &BSE_GAMMA);
				parsed.BSE_GAMMA = 1;
			} else if (strcmp(parameter_name, "TIMER")== 0) {
				PRINT_PARSED(PARAMDOC_TIMER);
				sscanf(values, "%d", &TIMER);
				parsed.TIMER = 1;
			} else {
				wprintf("unknown parameter: \"%s\".\n", line);
			}
		} else if (sscanf(line, "%s", parameter_name) == 1) {
			eprintf("too few values for parameter: \"%s\".\n", line);
			exit(1);
		}
	}
	fclose(in);
	
	/* quit if some parameters are unset */
#define CHECK_PARSED(A) \
	if (parsed.A == 0) { \
		eprintf("parameter \"%s\" unset.\n", #A); \
		allparsed = 0; \
	}
	
	CHECK_PARSED(GAMMA);
	CHECK_PARSED(INPUT_FILE);
	CHECK_PARSED(MASS_PC);
	CHECK_PARSED(MASS_BINS);
	
#undef CHECK_PARSED

/* but only warn if some other parameters are unset and default values are used */
#define CHECK_PARSED(A,DEFAULT,DOC) \
	if (parsed.A == 0) { \
		wprintf("parameter \"%s\" unset: using default value \"%s\").\n", #A, #DEFAULT); \
                A=DEFAULT; \
                rootfprintf(parsedfp, "# %s\n%s %s     # default value\n", DOC, #A, #DEFAULT); \
	}
	
	CHECK_PARSED(MASS_PC_BH_INCLUDE, 1, PARAMDOC_MASS_PC_BH_INCLUDE);
	CHECK_PARSED(PERTURB, 1, PARAMDOC_PERTURB);
	CHECK_PARSED(RELAXATION, 1, PARAMDOC_RELAXATION);
	CHECK_PARSED(TIDALLY_STRIP_STARS, 1, PARAMDOC_TIDALLY_STRIP_STARS);
	CHECK_PARSED(THETASEMAX, 1.0, PARAMDOC_THETASEMAX);
	CHECK_PARSED(STELLAR_EVOLUTION, 0, PARAMDOC_STELLAR_EVOLUTION);
    CHECK_PARSED(WRITE_STELLAR_INFO, 0, PARAMDOC_WRITE_STELLAR_INFO);
    CHECK_PARSED(WRITE_BH_INFO, 0, PARAMDOC_WRITE_BH_INFO);
    CHECK_PARSED(WRITE_RWALK_INFO, 0, PARAMDOC_WRITE_RWALK_INFO);
    CHECK_PARSED(WRITE_EXTRA_CORE_INFO, 0, PARAMDOC_WRITE_EXTRA_CORE_INFO);
    CHECK_PARSED(WRITE_PULSAR_INFO, 0, PARAMDOC_WRITE_PULSAR_INFO);
	CHECK_PARSED(CALCULATE10, 0, PARAMDOC_CALCULATE10);
	CHECK_PARSED(WIND_FACTOR, 1.0, PARAMDOC_WIND_FACTOR);
	CHECK_PARSED(TIDAL_TREATMENT, 0, PARAMDOC_TIDAL_TREATMENT);
	CHECK_PARSED(SS_COLLISION, 0, PARAMDOC_SS_COLLISION);
	CHECK_PARSED(TIDAL_CAPTURE, 0, PARAMDOC_TIDAL_CAPTURE);
	/*Sourav: new parameter*/
	CHECK_PARSED(STAR_AGING_SCHEME, 0, PARAMDOC_STAR_AGING_SCHEME);
	CHECK_PARSED(SAMPLESIZE, 1024, PARAMDOC_SAMPLESIZE);
	CHECK_PARSED(PREAGING, 0, PARAMDOC_PREAGING);
	CHECK_PARSED(BINBIN, 1, PARAMDOC_BINBIN);
	CHECK_PARSED(BINSINGLE, 1, PARAMDOC_BINSINGLE);
	CHECK_PARSED(STREAMS, 1, PARAMDOC_STREAMS);
	/*Meagan: new parameters for 3-body binary formation*/
	CHECK_PARSED(THREEBODYBINARIES, 0, PARAMDOC_THREEBODYBINARIES);
	CHECK_PARSED(MIN_BINARY_HARDNESS, 5.0, PARAMDOC_MIN_BINARY_HARDNESS);
	CHECK_PARSED(ONLY_FORM_BH_THREEBODYBINARIES, 1, PARAMDOC_ONLY_FORM_BH_THREEBODYBINARIES);
	// default - 1: three-body binary formation only allowed for black holes
	CHECK_PARSED(BH_LOSS_CONE, 0, PARAMDOC_BH_LOSS_CONE);
	CHECK_PARSED(MINIMUM_R, 0.0, PARAMDOC_MINIMUM_R);
	CHECK_PARSED(BH_R_DISRUPT_NB, 0., PARAMDOC_BH_R_DISRUPT_NB);
	CHECK_PARSED(CIRC_PERIOD_THRESHOLD, 1e-18, PARAMDOC_CIRC_PERIOD_THRESHOLD);
	CHECK_PARSED(BINARY_DISTANCE_BREAKING,0.1, PARAMDOC_BINARY_DISTANCE_BREAKING);
	CHECK_PARSED(BINARY_BREAKING_MIN,0.0, PARAMDOC_BINARY_BREAKING_MIN);

	CHECK_PARSED(T_MAX, 20.0, PARAMDOC_T_MAX);
	CHECK_PARSED(T_MAX_PHYS, 12.0, PARAMDOC_T_MAX_PHYS);
	CHECK_PARSED(T_MAX_COUNT, 1000000, PARAMDOC_T_MAX_COUNT);
	CHECK_PARSED(MAX_WCLOCK_TIME, 2592000, PARAMDOC_MAX_WCLOCK_TIME);
	CHECK_PARSED(CHECKPOINT_INTERVAL, 43200, PARAMDOC_CHECKPOINT_INTERVAL);
	CHECK_PARSED(CHECKPOINTS_TO_KEEP, 1, PARAMDOC_CHECKPOINTS_TO_KEEP);
	CHECK_PARSED(STOPATCORECOLLAPSE, 1, PARAMDOC_STOPATCORECOLLAPSE);
	CHECK_PARSED(TERMINAL_ENERGY_DISPLACEMENT, 0.5, PARAMDOC_TERMINAL_ENERGY_DISPLACEMENT);

	CHECK_PARSED(BH_SNAPSHOTTING, 0, PARAMDOC_BH_SNAPSHOTTING);
	CHECK_PARSED(BH_SNAPSHOT_DELTACOUNT, 50, PARAMDOC_BH_SNAPSHOT_DELTACOUNT);
	CHECK_PARSED(SNAPSHOTTING, 0, PARAMDOC_SNAPSHOTTING);
	CHECK_PARSED(SNAPSHOT_DELTACOUNT, 250, PARAMDOC_SNAPSHOT_DELTACOUNT);
	CHECK_PARSED(SNAPSHOT_DELTAT, 0.25, PARAMDOC_SNAPSHOT_DELTAT);
	CHECK_PARSED(SNAPSHOT_CORE_COLLAPSE, 0, PARAMDOC_SNAPSHOT_CORE_COLLAPSE);
	CHECK_PARSED(SNAPSHOT_CORE_BOUNCE, 0, PARAMDOC_SNAPSHOT_CORE_BOUNCE);
	CHECK_PARSED(SNAPSHOT_WINDOWS, NULL, PARAMDOC_SNAPSHOT_WINDOWS);
	CHECK_PARSED(SNAPSHOT_WINDOW_UNITS, "Trel", PARAMDOC_SNAPSHOT_WINDOW_UNITS);

	CHECK_PARSED(NUM_CENTRAL_STARS, 300, PARAMDOC_NUM_CENTRAL_STARS);
	CHECK_PARSED(IDUM, 0, PARAMDOC_IDUM);

	CHECK_PARSED(SEARCH_GRID, 0, PARAMDOC_SEARCH_GRID);
	CHECK_PARSED(SG_STARSPERBIN, 100, PARAMDOC_SG_STARSPERBIN);
	CHECK_PARSED(SG_MAXLENGTH, 1000000, PARAMDOC_SG_MAXLENGTH);
	CHECK_PARSED(SG_MINLENGTH, 1000, PARAMDOC_SG_MINLENGTH);
	CHECK_PARSED(SG_POWER_LAW_EXPONENT, 0.5, PARAMDOC_SG_POWER_LAW_EXPONENT);
	CHECK_PARSED(SG_MATCH_AT_FRACTION, 0.5, PARAMDOC_SG_MATCH_AT_FRACTION);
	CHECK_PARSED(SG_PARTICLE_FRACTION, 0.95, PARAMDOC_SG_PARTICLE_FRACTION);
	CHECK_PARSED(FORCE_RLX_STEP, 0, PARAMDOC_FORCE_RLX_STEP);
    CHECK_PARSED(DT_HARD_BINARIES, 0, PARAMDOC_DT_HARD_BINARIES);
    CHECK_PARSED(HARD_BINARY_KT, 1, PARAMDOC_HARD_BINARY_KT);
#ifdef EXPERIMENTAL
	CHECK_PARSED(BH_LC_FDT, 0.0, PARAMDOC_BH_LC_FDT);
	CHECK_PARSED(AVEKERNEL, 20, PARAMDOC_AVEKERNEL);
#endif
	CHECK_PARSED(MIN_CHUNK_SIZE, 20, PARAMDOC_MIN_CHUNK_SIZE);
    CHECK_PARSED(BH_AVEKERNEL, 10, PARAMDOC_BH_AVEKERNEL);
	CHECK_PARSED(APSIDES_PRECISION, 1.0e-11, PARAMDOC_APSIDES_PRECISION);
	CHECK_PARSED(APSIDES_MAX_ITER, 100, PARAMDOC_APSIDES_MAX_ITER);
	CHECK_PARSED(APSIDES_CONVERGENCE, 5.e-13, PARAMDOC_APSIDES_CONVERGENCE);
	CHECK_PARSED(OVERWRITE_RVIR, 0., PARAMDOC_OVERWRITE_RVIR);
	CHECK_PARSED(OVERWRITE_Z, 0., PARAMDOC_OVERWRITE_Z);
	CHECK_PARSED(OVERWRITE_RTID, 0., PARAMDOC_OVERWRITE_RTID);
	CHECK_PARSED(OVERWRITE_MCLUS, 0., PARAMDOC_OVERWRITE_MCLUS);
	CHECK_PARSED(BSE_NETA, 0.50, PARAMDOC_BSE_NETA);
	CHECK_PARSED(BSE_BWIND, 0.00, PARAMDOC_BSE_BWIND);
	CHECK_PARSED(BSE_HEWIND, 0.50, PARAMDOC_BSE_HEWIND);
	CHECK_PARSED(BSE_WINDFLAG, 1, PARAMDOC_BSE_WINDFLAG);
	CHECK_PARSED(BSE_PPSN, 0, PARAMDOC_BSE_PPSN);
	CHECK_PARSED(BSE_ALPHA1, 3.00, PARAMDOC_BSE_ALPHA1);
	CHECK_PARSED(BSE_LAMBDA, 0.50, PARAMDOC_BSE_LAMBDA);
	CHECK_PARSED(BSE_CEFLAG, 0, PARAMDOC_BSE_CEFLAG);
	CHECK_PARSED(BSE_TFLAG, 1, PARAMDOC_BSE_TFLAG);
	CHECK_PARSED(BSE_IFFLAG, 0, PARAMDOC_BSE_IFFLAG);
	CHECK_PARSED(BSE_WDFLAG, 1, PARAMDOC_BSE_WDFLAG);
	CHECK_PARSED(BSE_BHFLAG, 1, PARAMDOC_BSE_BHFLAG);
	CHECK_PARSED(BSE_BHSPINFLAG, 4, PARAMDOC_BSE_BHSPINFLAG);
	CHECK_PARSED(BH_RADIUS_MULTIPLYER, 3, PARAMDOC_BH_RADIUS_MULTIPLYER);
	CHECK_PARSED(BSE_NSFLAG, 1, PARAMDOC_BSE_NSFLAG);
	CHECK_PARSED(BSE_MXNS, 3.00, PARAMDOC_BSE_MXNS);
	CHECK_PARSED(BSE_BCONST, -3000.00, PARAMDOC_BSE_BCONST);
	CHECK_PARSED(BSE_CK, -1000.00, PARAMDOC_BSE_CK);
	CHECK_PARSED(BSE_IDUM, -999, PARAMDOC_BSE_IDUM);
	CHECK_PARSED(BSE_SIGMA, 265.00, PARAMDOC_BSE_SIGMA);
	CHECK_PARSED(BSE_BHSIGMAFRAC, 1.00, PARAMDOC_BSE_BHSIGMAFRAC);
	CHECK_PARSED(BSE_OPENING_ANGLE, 180., PARAMDOC_BSE_OPENING_ANGLE);
	CHECK_PARSED(BSE_BETA, 0.12500, PARAMDOC_BSE_BETA);
	CHECK_PARSED(BSE_EDDFAC, 1.00, PARAMDOC_BSE_EDDFAC);
	CHECK_PARSED(BSE_GAMMA, -1.00, PARAMDOC_BSE_GAMMA);
	CHECK_PARSED(TIMER, 0, PARAMDOC_TIMER);
#undef CHECK_PARSED

	/* exit if something is not set */
	if (!allparsed) {
		exit(1);
	}
	
#ifdef USE_MPI
if(myid==0)
#endif
	fclose(parsedfp);
	
    /* set-up snapshot window variables */
    parse_snapshot_windows(SNAPSHOT_WINDOWS);

	/* read the number of stars and possibly other parameters */
	/* MPI: Currently, all processors read the entire data (entire list of stars and binaries) from the file. The data partitioning is done in load_fits_file_data(). This might limit scalability since each node requires enough memory to store the entire data set. */
	cmc_read_fits_file(INPUT_FILE, &cfd, RESTART_TCOUNT);
	clus.N_STAR = cfd.NOBJ;
	clus.N_BINARY = cfd.NBINARY;
	if (OVERWRITE_RVIR>0.) {
		cfd.Rvir= OVERWRITE_RVIR;
	}
	if (OVERWRITE_RTID>0.) {
		cfd.Rtid= OVERWRITE_RTID;
	}
	if (OVERWRITE_Z>0.) {
		cfd.Z= OVERWRITE_Z;
	}
	if (OVERWRITE_MCLUS>0.) {
		cfd.Mclus= OVERWRITE_MCLUS;
	}

	R_MAX = cfd.Rtid;
	METALLICITY = cfd.Z;

	clus.N_STAR_NEW = clus.N_STAR;
	/* add 2 * clus.N_BINARY for binary disruptions */
	N_STAR_DIM = 2 + clus.N_STAR + 2 * clus.N_BINARY;
	/* remember we can form binaries from tidal capture */
	N_BIN_DIM = 2 + clus.N_STAR / 2 + clus.N_BINARY;
	/* safety factors, so we don't have to worry about memory management/garbage collection */
	N_STAR_DIM = (long) floor(1.1 * ((double) N_STAR_DIM));
	N_BIN_DIM = (long) floor(1.1 * ((double) N_BIN_DIM));

	
	/*********************************************/
	/* allocation of memory for global variables */
	/*********************************************/
	
	/* allocate memory for velocity dispersion array */
	sigma_array.n = 0;

#ifdef USE_MPI
	//MPI: Allocating only enough memory per processor.
	N_STAR_DIM_OPT = 1 + clus.N_STAR / procs + 2 * clus.N_BINARY / procs;
	N_BIN_DIM_OPT = clus.N_STAR / (2 * procs) + clus.N_BINARY / procs;
	//MPI: Giving a larger safety factor for the parallel version.
	N_STAR_DIM_OPT = (long) floor(1.5 * ((double) N_STAR_DIM_OPT));
	N_BIN_DIM_OPT = (long) floor(1.5 * ((double) N_BIN_DIM_OPT));
	if(RESTART_TCOUNT == 0){
		/* the main star array containing all star parameters */
		star = (star_t *) calloc(N_STAR_DIM_OPT, sizeof(star_t));
		/* the binary array containing all binary parameters */
		binary = (binary_t *) calloc(N_BIN_DIM_OPT, sizeof(binary_t));
	}
	sigma_array.r = (double *) calloc(N_STAR_DIM_OPT, sizeof(double));
	sigma_array.sigma = (double *) calloc(N_STAR_DIM_OPT, sizeof(double));
#else
	/* the main star array containing all star parameters */
	star = (star_t *) calloc(N_STAR_DIM, sizeof(star_t));
	/* the binary array containing all binary parameters */
	binary = (binary_t *) calloc(N_BIN_DIM, sizeof(binary_t));
	sigma_array.r = (double *) calloc(N_STAR_DIM, sizeof(double));
	sigma_array.sigma = (double *) calloc(N_STAR_DIM, sizeof(double));
#endif


	/* quantities calculated for various lagrange radii */
	mass_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ave_mass_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	no_star_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	densities_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ke_rad_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	ke_tan_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	v2_rad_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	v2_tan_r = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	mass_pc = (double *) calloc(MASS_PC_COUNT, sizeof(double));
	mass_bins = (double *) calloc(NO_MASS_BINS, sizeof(double));
	multi_mass_r = (double **) malloc(NO_MASS_BINS * sizeof(double *));
	for(i=0; i<NO_MASS_BINS; i++){
		multi_mass_r[i] = (double *) malloc(MASS_PC_COUNT * sizeof(double));
	}
	
	/*======= Reading values for the Lagrange radii =======*/
	curr_mass = (char *) strtok(MASS_PC, ",; ");
	sscanf(curr_mass, "%lf", &mass_pc[0]);
	
	for (i=1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; i++){
		sscanf(curr_mass, "%lf", &mass_pc[i]);
	}
	
	/*======= Reading values for the mass bins =======*/
	curr_mass = (char *) strtok(MASS_BINS, ",; ");
	sscanf(curr_mass, "%lf", &mass_bins[0]);

	for (i=1; (curr_mass = (char *) strtok(NULL, " ,;")) != NULL; i++){
		sscanf(curr_mass, "%lf", &mass_bins[i]);
	}

	/*======= Opening of output files =======*/
	if(RESTART_TCOUNT == 0)
		sscanf("w", "%s", outfilemode);
	else 
		sscanf("a", "%s", outfilemode);
	
/*
MPI: In the parallel version, IO is done in the following way. Some files require data only from the root node, and others need data from all nodes. The former are opened and written to only by the root node using C IO APIs. However, for the latter, the files are opened by all processors using MPI-IO. At places when the files are suposed to be written to in the serial version, in the parallel version, each processor writes the data into a string/char buffer. At the end of the timestep, all processors flush the data from the buffers into the corresponding files in parallel using MPI-IO. The code uses 5 variables for this process - the MPI-IO file pointer, which follows the format mpi_<serial fptr name>, 2 char buffers, which hav the format mpi_<serial fptr name>_buf and mpi_<ser fptr name>_wrbuf, and an int/longlong variables to maintain the length of the buffer (format mpi_<ser fptr name>_len) and the offset in the file (format mpi_<ser fptr name>_ofst_total) where data has to be written.
*/

#ifdef USE_MPI
    //MPI-IO: Following are files that require data only from the root node, and are opened only by the root node using standard C IO APIs.
    if(myid==0) {
#endif

        sprintf(outfile, "%s.lagrad.dat", outprefix);
        if ((lagradfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.dyn.dat", outprefix);
        if ((dynfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        //Sourav:new file containing info at 10% lagrange radius
        sprintf(outfile, "%s.lagrad_10_info.dat", outprefix);
        if ((lagrad10file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.avemass_lagrad.dat", outprefix);
        if ((ave_mass_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.nostar_lagrad.dat", outprefix);
        if ((no_star_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.rho_lagrad.dat", outprefix);
        if ((densities_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.ke_rad_lagrad.dat", outprefix);
        if ((ke_rad_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.ke_tan_lagrad.dat", outprefix);
        if ((ke_tan_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.v2_rad_lagrad.dat", outprefix);
        if ((v2_rad_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.v2_tan_lagrad.dat", outprefix);
        if ((v2_tan_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }
        sprintf(outfile, "%s.centmass.dat", outprefix);
        if ((centmass_file = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create output file \"%s\".\n", outfile);
            exit(1);
        }

        /* general binary information */
        sprintf(outfile, "%s.bin.dat", outprefix);
        if ((binaryfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create binary file \"%s\".\n", outfile);
            exit(1);
        }

        /* lagrange radii for multiple mass bins */
        mlagradfile = (FILE **) malloc((NO_MASS_BINS-1)*sizeof(FILE *));
        for(i=0; i<NO_MASS_BINS-1; i++){ /* NO_MASS_BINS need to be >=2 to be
                                            meaningful */
            sprintf(outfile, "%s.lagrad%ld-%g-%g.dat", outprefix, i,
                    mass_bins[i], mass_bins[i+1]);
            if ((mlagradfile[i] = fopen(outfile, outfilemode)) == NULL) {
                eprintf("cannot create output file \"%s\".\n", outfile);
                exit(1);
            }
        }

        /* File that contains data for various definitions of core radii */
        if (WRITE_EXTRA_CORE_INFO) {
            sprintf(outfile, "%s.core.dat", outprefix);
            if ((corefile = fopen(outfile, outfilemode)) == NULL) {
                eprintf("cannot create output file \"%s\".\n", outfile);
                exit(1);
            }

            /* the core that contains no remnants */
            append_core_header(corefile, "norem", 0);
            rootfprintf(corefile, "\n");
        }

        /* Meagan: output file for three-body binary formation */
		// file for bh summary, at each timestep
		sprintf(outfile, "%s.bh.dat", outprefix);
		if ((bhsummaryfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create bh.dat file %s\n", outfile);
			exit(1);
		}
	
		// file for escaping bh summary, at each timestep
		sprintf(outfile, "%s.esc.bh.dat", outprefix);
		if ((escbhsummaryfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create esc.bh.dat file %s\n", outfile);
			exit(1);
		}
	
		if(TIMER)
		{
			sprintf(outfile, "%s.timer.dat", outprefix);
			if ((timerfile = fopen(outfile, outfilemode)) == NULL) {
				eprintf("cannot create output file \"%s\".\n", outfile);
				exit(1);
			}
		}

		if(RESTART_TCOUNT == 0){
			/* Printing our headers */
			fprintf(lagradfile, "# Lagrange radii [code units]\n");
			fprintf(ave_mass_file, "# Average mass within Lagrange radii [M_sun]\n");
			fprintf(no_star_file, "# Number of stars within Lagrange radii [dimensionless]\n");
			fprintf(densities_file, "# Density within Lagrange radii [code units]\n");
			fprintf(ke_rad_file, "# Total radial kinetic energy within Lagrange radii [code units]\n");
			fprintf(ke_tan_file, "# Total tangential kinetic energy within Lagrange radii [code units]\n");
			fprintf(v2_rad_file, "# Sum of v_r within Lagrange radii [code units]\n");
			fprintf(v2_tan_file, "# Sum of v_t within Lagrange radii [code units]\n");
			for(i=0; i<NO_MASS_BINS-1; i++){
				fprintf(mlagradfile[i], "# Lagrange radii for %g < m < %g range [code units]\n", mass_bins[i], mass_bins[i+1]);
			}

			fprintf(lagradfile, "# 1:t");
			fprintf(ave_mass_file, "# 1:t");
			fprintf(no_star_file, "# 1:t");
			fprintf(densities_file, "# 1:t");
			fprintf(ke_rad_file, "# 1:t");
			fprintf(ke_tan_file, "# 1:t");
			fprintf(v2_rad_file, "# 1:t");
			fprintf(v2_tan_file, "# 1:t");
			for(i=0; i<NO_MASS_BINS-1; i++){
				fprintf(mlagradfile[i], "# 1:t");
			}

			for (i=0; i<MASS_PC_COUNT; i++) {
				fprintf(lagradfile, " %ld:r(%g)", i+2, mass_pc[i]);
				fprintf(ave_mass_file, " %ld:<m>(%g)", i+2, mass_pc[i]);
				fprintf(no_star_file, " %ld:N(%g)", i+2, mass_pc[i]);
				fprintf(densities_file, " %ld:rho(%g)", i+2, mass_pc[i]);
				fprintf(ke_rad_file, " %ld:T_r(%g)", i+2, mass_pc[i]);
				fprintf(ke_tan_file, " %ld:T_t(%g)", i+2, mass_pc[i]);
				fprintf(v2_rad_file, " %ld:V2_r(%g)", i+2, mass_pc[i]);
				fprintf(v2_tan_file, " %ld:V2_t(%g)", i+2, mass_pc[i]);
				for(j=0; j<NO_MASS_BINS-1; j++){
					fprintf(mlagradfile[j], " %ld:r(%g)", i+2, mass_pc[i]);
				}
			}

			fprintf(lagradfile, "\n");
			fprintf(ave_mass_file, "\n");
			fprintf(no_star_file, "\n");
			fprintf(densities_file, "\n");
			fprintf(ke_rad_file, "\n");
			fprintf(ke_tan_file, "\n");
			fprintf(v2_rad_file, "\n");
			fprintf(v2_tan_file, "\n");
			for(i=0; i<NO_MASS_BINS-1; i++){
				fprintf(mlagradfile[i], "\n");
			}

			fprintf(centmass_file, "# Information on central black hole [code units unless otherwise noted]\n");
			fprintf(centmass_file, "#1:t #2:cenma.m #3:Dt #4:rho_core #5:Etotal.tot #6:Etotal.K #7:Etotal.P\n");
			fprintf(dynfile, "# Dynamical information [code units]\n");
			fprintf(dynfile, "#1:t #2:Dt #3:tcount #4:N #5:M #6:VR #7:N_c #8:r_c #9:r_max #10:Etot #11:KE #12:PE #13:Etot_int #14:Etot_bin #15:E_cenma #16:Eesc #17:Ebesc #18:Eintesc #19:Eoops #20:Etot+Eoops #21:r_h #22:rho_0 #23:rc_spitzer #24:v0_rms #25:rc_nb #26.DMse(MSUN) #27.DMrejuv(MSUN) #28.N_c_nb\n");
			//Sourav:printing properties at 10 lagrange radii
			if (CALCULATE10){
				fprintf(lagrad10file, "#Dynamical information at 0.1 lagrange radius\n");
				fprintf(lagrad10file, "#1.t #2.Dt #3.tcount #4.N_10 #5.M_10 #6.N_s,10 #7.M_s,10 #8.N_b,10 #9.M_b_10 #10.r_10 #11.rho_10\n");
			}

			fprintf(binaryfile, "# Binary information [code units]\n");
			fprintf(binaryfile, "# 1:t 2:N_b 3:M_b 4:E_b 5:r_h,s 6:r_h,b 7:rho_c,s 8:rho_c,b 9:N_bb 10:N_bs 11:f_b,c 12:f_b 13:E_bb 14:E_bs 15:DE_bb 16:DE_bs 17:N_bc,nb 18:f_b,c,nb 19:N_bc \n");

			// print header
			fprintf(bhsummaryfile, "#1:tcount  #2:TotalTime  #3:Nbh,tot  #4:Nbh,single  #5:Nbinarybh  #6:Nbh-bh  #7:Nbh-nonbh  #8:Nbh-ns  #9:N_bh-wd  #10:N_bh-star  #11:Nbh-ms  #12:Nbh-postms #13:fb_bh [(# binaries containing a bh)/(total # systems containing a bh)\n");

			// print header
			fprintf(escbhsummaryfile, "# Ejected BHs\n#1:tcount  #2:TotalTime  #3:Nbh,tot  #4:Nbh,single  #5:Nbinarybh  #6:Nbh-bh  #7:Nbh-nonbh  #8:Nbh-ns  #9:N_bh-wd  #10:N_bh-star  #11:Nbh-ms  #12:Nbh-postms #13:fb_bh [(# binaries containing a bh)/(total # systems containing a bh)]\n");

			if(TIMER)
				fprintf(timerfile, "#1:tcount\t#2:t_cen_calc\t#3:t_timestep\t#4:t_dyn\t#5:t_se\t#6:t_orb\t#7:t_tid_str\t#8:t_sort\t#9:t_postsort_comm\t#10:t_pot_cal\t#11:t_ener_con3\t#12:t_calc_io_vars1\t#13:t_calc_io_vars1\t#14:t_comp_ener\t#15:t_upd_vars\t#16:t_io\t#17:t_io_ignore\t#18:t_oth\t#19:t_sort_lsort1\t#20:t_sort_splitters\t#21:t_sort_a2a\t#22:t_sort_lsort2\t#23:t_sort_oth\t#24:t_sort_lb\t#25:t_sort_only\n");
		}/*if (RESTARTING_TCOUNT == 0)*/

#ifdef USE_MPI
    }
#endif


#ifdef USE_MPI
    //MPI3-IO: Files that might require data from all nodes, and are opened by all procs using MPI-IO. In the serial version, these are just opened as normal (see under #else below).
	
	int MPI_MODE_RESTART;
	if(RESTART_TCOUNT > 0)
		MPI_MODE_RESTART = (MPI_MODE_APPEND | MPI_MODE_WRONLY);
	else
		MPI_MODE_RESTART = (MPI_MODE_CREATE | MPI_MODE_WRONLY);

    sprintf(outfile, "%s.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_logfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_logfile, 0);

    // output files for binaries 
    sprintf(outfile, "%s.binint.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_binintfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_binintfile, 0);

    sprintf(outfile, "%s.esc.dat", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_escfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_escfile, 0);

    sprintf(outfile, "%s.collision.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_collisionfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_collisionfile, 0);

    sprintf(outfile, "%s.tidalcapture.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_tidalcapturefile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_tidalcapturefile, 0);

    sprintf(outfile, "%s.semergedisrupt.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_semergedisruptfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_semergedisruptfile, 0);

    sprintf(outfile, "%s.removestar.log", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_removestarfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_removestarfile, 0);

    sprintf(outfile, "%s.relaxation.dat", outprefix);
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_relaxationfile);
	if(RESTART_TCOUNT == 0)
		MPI_File_set_size(mpi_relaxationfile, 0);

    if (THREEBODYBINARIES)
    {
        sprintf(outfile, "%s.3bb.log", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_threebbfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_threebbfile, 0);

        sprintf(outfile, "%s.3bbprobability.log", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_threebbprobabilityfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_threebbprobabilityfile, 0);

        sprintf(outfile, "%s.lightcollision.log", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_lightcollisionfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_lightcollisionfile, 0);

        sprintf(outfile, "%s.3bbdebug.log", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_threebbdebugfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_threebbdebugfile, 0);
    }

    // Meagan: extra output for black holes
    if (WRITE_BH_INFO) {
        sprintf(outfile, "%s.bhformation.dat", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_newbhfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_newbhfile, 0);

        sprintf(outfile, "%s.bhmerger.dat", outprefix);
        MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_bhmergerfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_bhmergerfile, 0);
        
    }

	if(WRITE_PULSAR_INFO)
	{
		sprintf(outfile, "%s.pulsars.dat", outprefix);
		MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_RESTART, MPI_INFO_NULL, &mpi_pulsarfile);
		if(RESTART_TCOUNT == 0)
			MPI_File_set_size(mpi_pulsarfile, 0);
	}

#else

	/* output files for binaries */
	/* file for binary interaction information */
	sprintf(outfile, "%s.binint.log", outprefix);
	if ((binintfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create binintlog file \"%s\".\n", outfile);
		exit(1);
	}
	
	/* File for parameters of escaping stars */
	sprintf(outfile, "%s.esc.dat", outprefix);
	if ((escfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create escapers file \"%s\".\n", outfile);
		exit(1);
	}

	/* Collision log file */
	sprintf(outfile, "%s.collision.log", outprefix);
	if ((collisionfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create collision log file \"%s\".\n", outfile);
		exit(1);
	}

	/* Tidal capture log file */
	sprintf(outfile, "%s.tidalcapture.log", outprefix);
	if ((tidalcapturefile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create tidal capture log file \"%s\".\n", outfile);
		exit(1);
	}

	/* Stellar evolution merger log file */
	sprintf(outfile, "%s.semergedisrupt.log", outprefix);
	if ((semergedisruptfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create stellar evolution merger/disruption log file \"%s\".\n", outfile);
		exit(1);
	}

	/*Sourav: old star removal info file*/
	sprintf(outfile, "%s.removestar.log", outprefix);
	if ((removestarfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create old star removal log file \"%s\".\n", outfile);
		exit(1);
	}

	/* Relaxation data file */
	sprintf(outfile, "%s.relaxation.dat", outprefix);
	if ((relaxationfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create relaxation data file \"%s\".\n", outfile);
		exit(1);
	}

	sprintf(outfile, "%s.log", outprefix);
	if ((logfile = fopen(outfile, outfilemode)) == NULL) {
		eprintf("cannot create log output file \"%s\".\n", outfile);
		exit(1);
	}

    if (THREEBODYBINARIES)
    {
        /* Meagan: output file for three-body binary formation */
        sprintf(outfile, "%s.3bb.log", outprefix);
        if ((threebbfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create 3bb log output file \"%s\".\n", outfile);
            exit(1);
        }

        sprintf(outfile, "%s.3bbprobability.log", outprefix);
        if ((threebbprobabilityfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create three-body formation probability log output file \"%s\".\n", outfile);
            exit(1);
        }

        sprintf(outfile, "%s.lightcollision.log", outprefix);
        if ((lightcollisionfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create light collision log output file \"%s\".\n", outfile);
            exit(1);
        }

        sprintf(outfile, "%s.3bbdebug.log", outprefix);
        if ((threebbdebugfile = fopen(outfile, outfilemode)) == NULL) {
            eprintf("cannot create 3bb debug output file \"%s\".\n", outfile);
            exit(1);
        }
    }

	// Meagan: extra output for black holes
	if (WRITE_BH_INFO) {
		// file for info about newly formed BHs
		sprintf(outfile, "%s.bhformation.dat", outprefix);
		if ((newbhfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create bhformation.dat file %s\n", outfile);
			exit(1);
		}
		sprintf(outfile, "%s.bhmerger.dat", outprefix);
		if ((bhmergerfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create bhmerger.dat file %s\n", outfile);
			exit(1);
		}
    }

	/* Pulsar log file */
	if(WRITE_PULSAR_INFO)
	{
		sprintf(outfile, "%s.pulsars.dat", outprefix);
		//if ((pulsarfile = gzopen(outfile, outfilemode)) == NULL) {
		if ((pulsarfile = fopen(outfile, outfilemode)) == NULL) {
			eprintf("cannot create pulsar file \"%s\".\n", outfile);
			exit(1);
		}
	}
#endif

	//MPI: Headers are written out only by the root node.
   // print header
    if(RESTART_TCOUNT == 0){
		pararootfprintf(escfile, "#1:tcount #2:t #3:m #4:r #5:vr #6:vt #7:r_peri #8:r_apo #9:Rtidal #10:phi_rtidal #11:phi_zero #12:E #13:J #14:id #15:binflag #16:m0[MSUN] #17:m1[MSUN] #18:id0 #19:id1 #20:a #21:e #22:startype #23:bin_startype0 #24:bin_startype1 #25:rad0 #26:rad1 #27:tb #28:lum0 #29:lum1 #30:massc0 #31:massc1 #32:radc0 #33:radc1 #34:menv0 #35:menv1 #36:renv0 #37:renv1 #38:tms0 #39:tms1 #40:dmdt0 #41:dmdt1 #42:radrol0 #43:radrol1 #44:ospin0 #45:ospin1 #46:B0 #47:B1 #48:formation0 #49:formation1 #50:bacc0 #51:bacc1 #52:tacc0 $53:tacc1 #54:mass0_0 #55:mass0_1 #56:epoch0 #57:epoch1 \n");
	   // print header
		pararootfprintf(collisionfile, "# time interaction_type id_merger(mass_merger) id1(m1):id2(m2):id3(m3):... (r) type_merger type1 ...\n");
	   // print header
		pararootfprintf(tidalcapturefile, "# time interaction_type (id1,m1,k1)+(id2,m2,k2)->[(id1,m1,k1)-a,e-(id2,m2,k2)]\n");
	   // print header
		pararootfprintf(semergedisruptfile, "# time interaction_type id_rem(mass_rem) id1(m1):id2(m2) (r)\n");
	   //Sourav:  print header
		pararootfprintf(removestarfile, "#single destroyed: time star_id star_mass(MSun) star_age(Gyr) star_birth(Gyr) star_lifetime(Gyr)\n");
		pararootfprintf(removestarfile, "#binary destroyed: time obj_id bin_id removed_comp_id left_comp_id m1(MSun) m2(MSun) removed_m(MSun) left_m(MSun) left_m_sing(MSun) star_age(Gyr) star_birth(Gyr) star_lifetime(Gyr)\n");

		if (THREEBODYBINARIES)
			{
			// print header
			pararootfprintf(threebbfile, "#1:time #2:k1 #3:k2 #4:k3 #5:id1 #6:id2 #7:id3 #8:m1 #9:m2 #10:m3 #11:ave_local_mass #12:n_local #13:sigma_local #14:eta #15:Eb #16:ecc #17:a[AU] #18:r_peri[AU] #19:r(bin) #20:r(single) #21:vr(bin) #22:vt(bin) #23:vr(single) #24:vt(single) #25:phi(bin) #26:phi(single) #27:delta_PE #28:delta_KE #29:delta_E(interaction) #30:delta_E(cumulative) #31:N_3bb\n");
			// print header
			pararootfprintf(threebbprobabilityfile, "#1:time #2:dt #3:dt*N/log(gamma*N) #3:Rate_3bb #4:P_3bb #5:r\n### average rate and probability of three-body binary formation in the timestep; calculated from the innermost 300 triplets of single stars considered for three-body binary formation\n");
			// print header
			pararootfprintf(lightcollisionfile, "#1:time #2:k1 #3:k2 #4:k3 #5:id1 #6:id2 #7:id3 #8:m1 #9:m2 #10:m3 #11:type1 #12:type2 #13:type3 #14:rad1 #15:rad2 #16:rad3 #17:Eb #18:ecc #19:a(au) #20:rp(au)\n");
			// print header
			pararootfprintf(threebbdebugfile, "#1:k1 #2:k2 #3:k3 #4:id1 #5:id2 #6:id3 #7:r1 #8:r2 #9:r3 #10:m1 #11:m2 #12:m3 #13:v1 #14:v1[1] #15:v1[2] #16:v1[2] #17:v2 #18:v2[1] #19:v2[2] #20:v2[3] #21:v3 #22:v3[1] #23:v3[2] #24:v3[3] #25:v1_cmf #26:v1_cmf[1] #27:v1_cmf[2] #28:v1_cmf[3] #29:v2_cmf #30:v2_cmf[1] #31:v2_cmf[2] #32:v2_cmf[3] #33:v3_cmf #34:v3_cmf[1] #35:v3_cmf[2] #36:v3_cmf[3] #37:knew #38:bin_id #39:bin_r #40:bin_m #41:vs_cmf #42:vs_cmf[1] #43:vs_cmf[2] #44:vs_cmf[3] #45:vb_cmf #46:vb_cmf[1] #47:vb_cmf[2] #48:vb_cmf[3] #49:vs #50:vs[1] #51:vs[2] #52:vs[3] #53:vb #55:vb[1] #56:vb[2] #57:vb[3] #58:ave_local_m #59:sigma #60:eta #61:Eb #62:ecc #63:rp #64:a #65:PE_i #66:PE_f #67:KE_cmf_i #68:KE_cmf_f #69:KE_i #70:KE_f #71:delta_PE #72:delta_KE #73:delta_E\n");
		}

		// print header
		/*pararootfprintf(escbhsummaryfile, "# Ejected BHs\n#1:tcount  #2:TotalTime  #3:Nbh,tot  #4:Nbh,single  #5:Nbinarybh  #6:Nbh-bh  #7:Nbh-nonbh  #8:Nbh-ns  #9:N_bh-wd  #10:N_bh-star  #11:Nbh-ms  #12:Nbh-postms #13:fb_bh [(# binaries containing a bh)/(total # systems containing a bh)]\n");
		*/

		// print header
		if (WRITE_BH_INFO)
			pararootfprintf(newbhfile,"#1:time #2:r #3.binary? #4:ID #5:zams_m #6:m_progenitor #7:bh mass #8:bh_spin #9:birth-kick(km/s) #10-25:vsarray\n");
			pararootfprintf(bhmergerfile,"#1:time #2:r #3.type #4:id1 #5:id2 #6:m1[MSUN] #7:m2[MSUN] #8:a1 #9:a2 #10:m_final[MSUN] 11:a_final 12:vkick[km/s]\n");
			pararootfprintf(bhmergerfile,"#NOTE: if repeated mergers occur in fewbody (binary-single or binary-binary), the initial masses will be wrong; check collision.log\n");
	//"#1:tcount  #2:TotalTime  #3:bh  #4:bh_single  #5:bh_binary  #6:bh-bh  #7:bh-ns  #8:bh-wd  #9:bh-star  #10:bh-nonbh  #11:fb_bh  #12:bh_tot  #13:bh_single_tot  #14:bh_binary_tot  #15:bh-bh_tot  #16:bh-ns_tot  #17:bh-wd_tot  #18:bh-star_tot  #19:bh-nonbh_tot  #20:fb_bh_tot\n");

		/* print header */
		if(WRITE_PULSAR_INFO)
			pararootfprintf(pulsarfile, "tcount    TotalTime    Star_id      Rperi    Rapo    R     VR    VT    PHI    PHIr0    PHIrt    kick    Binary_id1    Binary_id2    kw2     P     B    formation     bacc    tacc    B0   TB     M2    M1     e     R2/RL2     dm1/dt   \n");
	} /*if(RESTART_TCOUNT == 0)*/
}


/**
* @brief close file buffers/pointers
*/
void close_buffers(void)
{
//MPI: These files are written to only by the root, and hence are closed only by root.
#ifdef USE_MPI
    if(myid==0)
#endif
        close_root_buffers();

	 //MPI: These include the rest that are written to by all procs, and hence depending on whether the serial or parallel version is being compiled, close the corresponding file pointers.
#ifdef USE_MPI
    mpi_close_node_buffers();
#else
    close_node_buffers();
#endif
}

/**
* @brief Closes some of the file pointers - of files which require writing only by the root node.
*/
void close_root_buffers(void)
{
	int i;

	fclose(lagradfile);
	fclose(dynfile);
	fclose(lagrad10file);
	fclose(ave_mass_file);
	fclose(no_star_file);
	fclose(densities_file);
	fclose(ke_rad_file);
	fclose(ke_tan_file);
	fclose(v2_rad_file);
	fclose(v2_tan_file);
	fclose(centmass_file);
	fclose(binaryfile);
	for(i=0; i<NO_MASS_BINS-1; i++)
		fclose(mlagradfile[i]);
	if (WRITE_EXTRA_CORE_INFO) 
		fclose(corefile);
    /* Meagan's 3bb stuff */
    fclose(bhsummaryfile);
    fclose(escbhsummaryfile);
	 if(TIMER)
		 fclose(timerfile);
}

/**
* @brief Closes some of the file pointers - of files which require writing by the all nodes.
*/
void close_node_buffers(void)
{
	fclose(logfile);
	fclose(binintfile);
	fclose(escfile);
	fclose(collisionfile);
	fclose(tidalcapturefile);
	fclose(semergedisruptfile);
	fclose(relaxationfile);
	/*Sourav: closing the file I opened*/
	fclose(removestarfile);
    /* Meagan: close 3bb log file */
    if (THREEBODYBINARIES)
    {
        fclose(threebbfile);
        fclose(threebbprobabilityfile);
        fclose(lightcollisionfile);
        fclose(threebbdebugfile);
    }
    if (WRITE_BH_INFO) {
        fclose(newbhfile);
        fclose(bhmergerfile);
    }
	 if(WRITE_PULSAR_INFO)
	 {
		 fclose(pulsarfile);
	 }
}

#ifdef USE_MPI
/**
* @brief Closes the MPI file pointers - of files which require writing only by the all nodes.
*/
void mpi_close_node_buffers(void)
{
	MPI_File_close(&mpi_logfile);
	MPI_File_close(&mpi_binintfile);
	MPI_File_close(&mpi_escfile);
	MPI_File_close(&mpi_collisionfile);
	MPI_File_close(&mpi_tidalcapturefile);
	MPI_File_close(&mpi_semergedisruptfile);
	MPI_File_close(&mpi_relaxationfile);
	/*Sourav: closing the file I opened*/
	MPI_File_close(&mpi_removestarfile);
    /* Meagan: close 3bb log file */
    if (THREEBODYBINARIES)
    {
        MPI_File_close(&mpi_threebbfile);
        MPI_File_close(&mpi_threebbprobabilityfile);
        MPI_File_close(&mpi_lightcollisionfile);
        MPI_File_close(&mpi_threebbdebugfile);
    }
    if (WRITE_BH_INFO) {
        MPI_File_close(&mpi_newbhfile);
        MPI_File_close(&mpi_bhmergerfile);
    }
	 if(WRITE_PULSAR_INFO)
	 {
		 MPI_File_close(&mpi_pulsarfile);
	 }

	//TEMPORARY
/*
#ifdef USE_MPI
	if(myid==0)
	{
		fclose(fp_lagrad);
		fclose(fp_log);
	}
#endif
*/
}
#endif

/**
* @brief traps signals
*/
void trap_sigs(void)
{
	/* Catch some signals */
	signal(SIGINT, exit_cleanly_old);
	signal(SIGTERM, exit_cleanly_old);
	signal(SIGQUIT, exit_cleanly_old);
	signal(SIGUSR1, toggle_debugging);
	
	/* override GSL error handler */
	//gsl_set_error_handler(&sf_gsl_errhandler);
}

/**
* @brief print handy script for converting output files to physical units
*/
void print_conversion_script(void)
{
//MPI: Just letting the root node handle this.
#ifdef USE_MPI
if(myid==0)
{
#endif
	char dummystring[1024];
	FILE *ofp;
	/* BOOKMARK */

	sprintf(dummystring, "%s.conv.sh", outprefix);
	if ((ofp = fopen(dummystring, "w")) == NULL) {
		eprintf("cannot create convsersion script file \"%s\".\n", dummystring);
		exit(1);
	}
	
	fprintf(ofp, "#!/bin/bash\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "# outfile prefix\n");
	fprintf(ofp, "outprefix=%s\n", outprefix);
	fprintf(ofp, "# code unit of mass (cgs)\n");
	fprintf(ofp, "massunitcgs=%g\n", units.m);
	fprintf(ofp, "# code unit of mass (M_sun)\n");
	fprintf(ofp, "massunitmsun=%g\n", units.m/MSUN);
	fprintf(ofp, "# code unit of stellar mass (cgs)\n");
	fprintf(ofp, "mstarunitcgs=%g\n", units.mstar);
	fprintf(ofp, "# code unit of stellar mass (M_Sun)\n");
	fprintf(ofp, "mstarunitmsun=%g\n", units.mstar/MSUN);
	fprintf(ofp, "# code unit of length (cgs)\n");
	fprintf(ofp, "lengthunitcgs=%g\n", units.l);
	fprintf(ofp, "# code unit of length (parsecs)\n");	
	fprintf(ofp, "lengthunitparsec=%g\n", units.l/PARSEC);
	fprintf(ofp, "# code unit of time (cgs)\n");
	fprintf(ofp, "timeunitcgs=%g\n", units.t * clus.N_STAR / log(GAMMA * clus.N_STAR));
	fprintf(ofp, "# code unit of time (Myr)\n");
	fprintf(ofp, "timeunitsmyr=%g\n", units.t * clus.N_STAR / log(GAMMA * clus.N_STAR) / (1.0e6 * YEAR));
	fprintf(ofp, "# N-body  unit of time (cgs)\n");
	fprintf(ofp, "nbtimeunitcgs=%g\n", units.t);
	fprintf(ofp, "# N-body unit of time (Myr)\n");
	fprintf(ofp, "nbtimeunitsmyr=%g\n", units.t / (1.0e6 * YEAR));
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.dyn.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $8/$21}' > $outprefix.tmyr_rcrh.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_rcrh.dat > $outprefix.tmyr_rcrh-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.dyn.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $25/$21}' > $outprefix.tmyr_rcnbrh.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_rcnbrh.dat > $outprefix.tmyr_rcnbrh-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.dyn.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $7/(4.0/3.0*3.14159265*$8*'$lengthunitparsec')^3}' > $outprefix.tmyr_nc.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_nc.dat > $outprefix.tmyr_nc-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.dyn.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $5}' > $outprefix.tmyr_m.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_m.dat > $outprefix.tmyr_m-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.bin.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $11}' > $outprefix.tmyr_fbc.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_fbc.dat > $outprefix.tmyr_fbc-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.bin.dat | grep -vE '^#' | awk '{print $1*'$timeunitsmyr', $12}' > $outprefix.tmyr_fb.dat\n");
	fprintf(ofp, "prunedata.pl -d 30 $outprefix.tmyr_fb.dat > $outprefix.tmyr_fb-pruned.dat\n");
	fprintf(ofp, "\n");
	fprintf(ofp, "cat $outprefix.tmyr_nc.dat | awk '{print NR, $2}' > $outprefix.NR_nc.dat\n");
	fprintf(ofp, "cat $outprefix.tmyr_fbc.dat | awk '{print NR, $2}' > $outprefix.NR_fbc.dat\n");
	fprintf(ofp, "join $outprefix.NR_fbc.dat $outprefix.NR_nc.dat | awk '{print $2, $3}' > $outprefix.fbc_nc.dat\n");
	fprintf(ofp, "prunedata.pl -n 300 $outprefix.fbc_nc.dat > $outprefix.fbc_nc-pruned.dat\n");
	fprintf(ofp, "head -n 1 $outprefix.fbc_nc.dat > $outprefix.fbc_nc-arrow.dat\n");
	fprintf(ofp, "tail -n 1 $outprefix.fbc_nc.dat >> $outprefix.fbc_nc-arrow.dat\n");
	fclose(ofp);

	chmod(dummystring, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#ifdef USE_MPI
}
#endif
}

/**
* @brief routines for printing star/binary info in a unified log format
*
* @param k index of star
* @param string[MAX_STRING_LENGTH] ?
*
* @return ?
*/
char *sprint_star_dyn(long k, char string[MAX_STRING_LENGTH])
{
#ifdef USE_MPI
	snprintf(string, MAX_STRING_LENGTH, "(%ld,%.3g,%d)",
		 star[k].id, star_m[get_global_idx(k)]*units.mstar/FB_CONST_MSUN, star[k].se_k);
#else
	snprintf(string, MAX_STRING_LENGTH, "(%ld,%.3g,%d)",
		 star[k].id, star[k].m*units.mstar/FB_CONST_MSUN, star[k].se_k);
#endif
	return(string);
}

/**
* @brief ?
*
* @param k index of star
* @param string[MAX_STRING_LENGTH] ?
*
* @return ?
*/
char *sprint_bin_dyn(long k, char string[MAX_STRING_LENGTH])
{
	long bi=star[k].binind;

	snprintf(string, MAX_STRING_LENGTH, "[(%ld,%.3g,%d)-%.3g,%.3g-(%ld,%.3g,%d)]", 
		 binary[bi].id1, binary[bi].m1*units.mstar/FB_CONST_MSUN, binary[bi].bse_kw[0],
		 binary[bi].a*units.l/FB_CONST_AU, binary[bi].e,
		 binary[bi].id2, binary[bi].m2*units.mstar/FB_CONST_MSUN, binary[bi].bse_kw[1]);

	return(string);
}

/**
* @brief ?
*
* @param param_string ?
*/
void parse_snapshot_windows(char *param_string) {
  char *cur_window, *intern_window=NULL, *intern_param=NULL;
  char *cur_wstring, *cur_pstring;
  int j;

  if (param_string==NULL) {
    return;
  }
  snapshot_window_count= 0;
  snapshot_windows= NULL;
  cur_wstring= param_string;
  while ((cur_window = strtok_r(cur_wstring,";", &intern_window))!= NULL) {
    snapshot_window_count++;
    snapshot_windows = (double *) realloc(snapshot_windows, 3*snapshot_window_count*sizeof(double));
    cur_wstring= NULL;
    for (j=0, cur_pstring=cur_window; j< 3; j++, cur_pstring=NULL) {
      char *wparam;
      int cur_wpidx;
      wparam= strtok_r(cur_pstring, ",", &intern_param);
      if (wparam==NULL) {
        eprintf("Error parsing snapshot window list.\n");
        eprintf("The current window is %s, and the config parameter is %s.\n", 
            cur_window, param_string);
        free_arrays();
        exit(-1);
      }
      cur_wpidx= (snapshot_window_count-1)*3+j;
      sscanf(wparam, "%lf", &(snapshot_windows[cur_wpidx]));
    }
  }

  dprintf("Finished parsing window list.\n");
  for (j=0; j< snapshot_window_count; j++) {
    dprintf("Window %i: ", j);
    dprintf("start %g, step %g, stop %g\n", snapshot_windows[3*j], snapshot_windows[3*j+1], snapshot_windows[3*j+2]);
  }
  snapshot_window_counters= (int *) calloc(snapshot_window_count, sizeof(int));
}

/**
* @brief ?
*/
void print_snapshot_windows(void) {
  int i, step_counter;
  double start, stop, step, total_time;
  total_time = 0;

  if (!snapshot_window_count || !SNAPSHOTTING) return;

  if (strncmp(SNAPSHOT_WINDOW_UNITS, "Trel", 5)==0) {
    total_time= TotalTime;
  } else if (strncmp(SNAPSHOT_WINDOW_UNITS, "Gyr", 4)==0) {
    total_time= TotalTime * units.t * clus.N_STAR / log(GAMMA * clus.N_STAR) / YEAR/1e9;
  } else if (strncmp(SNAPSHOT_WINDOW_UNITS, "Tcr", 4)==0) {
    total_time= TotalTime * log(GAMMA * clus.N_STAR) / clus.N_STAR;
  } else {
    eprintf("Unrecognized unit %s.", SNAPSHOT_WINDOW_UNITS);
    exit_cleanly(-1, __FUNCTION__);
  }

  for (i=0; i<snapshot_window_count; i++) {
    step_counter= snapshot_window_counters[i];
    start= snapshot_windows[i*3];
    step=  snapshot_windows[i*3+1];
    stop=  snapshot_windows[i*3+2];
    if (total_time>= start+step_counter*step && total_time<=stop) {
      char outfile[100];
      sprintf(outfile, "%s.w%02i_snap%04d.dat.gz", outprefix, i+1, step_counter+1);
      write_snapshot(outfile, 0);
//		print_denprof_snapshot(outfile);

      snapshot_window_counters[i]++;
      dprintf("Wrote snapshot #%i for time window %i (%s) at time %g %s.\n", step_counter+1, i+1, outfile, 
          total_time, SNAPSHOT_WINDOW_UNITS);
    }
  }
}

/**
* @brief ?
*
* @return ?
*/
int valid_snapshot_window_units(void) {
  int valid;

  valid=0;
  if (strncmp(SNAPSHOT_WINDOW_UNITS, "Trel", 5)==0) {
    valid=1;
  } else if (strncmp(SNAPSHOT_WINDOW_UNITS, "Gyr", 4)) {
    valid=1;
  } else if (strncmp(SNAPSHOT_WINDOW_UNITS, "Tcr", 4)) {
    valid=1;
  } 

  return (valid);
}

/**
* @brief writes out snapshot to the given file
*
* @param filename name of the file
* @param bh_only if bh_only>0 this'll print only BHs.
*/
void write_snapshot(char *filename, int bh_only) {
	long i, j;
	j=0;
	double m, r, phi;

#ifdef USE_MPI
	//Serializing the snapshot printing.
	int k;
	for(k=0; k<procs; k++)
	{
		if(myid==k)
		{
			//removing file if already exists.
			if(myid==0)
				remove(filename);
#endif

			if ((snapfile = (FILE *) gzopen(filename, "ab")) == NULL) {
				eprintf("cannot create 2D snapshot file %s\n", filename);
				exit_cleanly(1, __FUNCTION__);
			}

#ifdef USE_MPI
			//Header printed only by root node.
			if(myid==0)
			{
#endif
				// print useful header
				gzprintf(snapfile, "# t=%.8g [code units]; All quantities below are in code units unless otherwise specified.\n", TotalTime);
				gzprintf(snapfile, "#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.bin.Eb 25.eta 26.star.phi#27:rad0 #28:rad1 #29:tb #30:lum0 #31:lum1 #32:massc0 #33:massc1 #34:radc0 #35:radc1 #36:menv0 #37:menv1 #38:renv0 #39:renv1 #40:tms0 #41:tms1 #42:dmdt0 #43:dmdt1 #44:radrol0 #45:radrol1 #46:ospin0 #47:ospin1 #48:B0 #49:B1 #50:formation0 #51:formation1 #52:bacc0 #53:bacc1 #54:tacc0 $55:tacc1 #56:mass0_0 #57:mass0_1 #58:epoch0 #59:epoch1 \n");
#ifdef USE_MPI
			}
#endif

			// then print data
#ifdef USE_MPI
			for (i=1; i<=clus.N_MAX_NEW; i++) {
				long g_i = get_global_idx(i);
				m = star_m[g_i];
				r = star_r[g_i];
				phi = star_phi[g_i];
#else
			for (i=1; i<=clus.N_MAX; i++) {
				m = star[i].m;
				r = star[i].r;
				phi = star[i].phi;
#endif
				j=star[i].binind;
				//if bh_only>0, print only BHs
				if( (bh_only==0) || ( (bh_only!=0) && (star[i].se_k==14 || binary[j].bse_kw[0]==14 || binary[j].bse_kw[1]==14) ) )
				{
					gzprintf(snapfile, "%ld %.8g %.8g %.8g %.8g %.8g %.8g ",
							star[i].id, m * (units.m / clus.N_STAR) / MSUN,
							r, star[i].vr, star[i].vt,
							star[i].E, star[i].J);
					if (j) {
						gzprintf(snapfile, "1 %.8g %.8g %ld %ld %.8g %.8g ",
								binary[j].m1 * (units.m / clus.N_STAR) / MSUN,
								binary[j].m2 * (units.m / clus.N_STAR) / MSUN,
								binary[j].id1, binary[j].id2,
								binary[j].a * units.l / AU, binary[j].e);
					} else {
						gzprintf(snapfile, "0 0 0 0 0 0 0 ");
					}

					if (j == 0) {
						gzprintf(snapfile, "%d %.8g %.8g -100 -100 -100 -100 -100 -100 -100 -100 ",
								star[i].se_k, star[i].se_lum, star[i].rad * units.l / RSUN);
					} else {
						gzprintf(snapfile, "0 0 0 %d %d %.8g %.8g %.8g %.8g %.8g %.8g ",
								binary[j].bse_kw[0], binary[j].bse_kw[1],
								binary[j].bse_lum[0], binary[j].bse_lum[1],
								binary[j].rad1*units.l/RSUN, binary[j].rad2*units.l/RSUN,
								-(binary[j].m1/clus.N_STAR)*(binary[j].m2/clus.N_STAR)/(2*binary[j].a),
                                 (binary[j].m1 * binary[j].m2 * sqr(madhoc)) /
                                 (binary[j].a * sqrt(calc_average_mass_sqr(i,clus.N_MAX)) * sqr(sigma_array.sigma[i])));
					}
					gzprintf(snapfile, "%0.12g ", phi);
					if (j == 0) {
						gzprintf(snapfile, "na na na na na na na na na na na na na na na na na na na na na na na na na na na na na na na na na \n");
					} else {
						gzprintf(snapfile, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
								binary[j].bse_radius[0], binary[j].bse_radius[1], binary[j].bse_tb, binary[j].bse_lum[0], binary[j].bse_lum[1], binary[j].bse_massc[0], binary[j].bse_massc[1], binary[j].bse_radc[0], binary[j].bse_radc[1], binary[j].bse_menv[0], binary[j].bse_menv[1], binary[j].bse_renv[0], binary[j].bse_renv[1], binary[j].bse_tms[0], binary[j].bse_tms[1], binary[j].bse_bcm_dmdt[0], binary[j].bse_bcm_dmdt[1], binary[j].bse_bcm_radrol[0], binary[j].bse_bcm_radrol[1], binary[j].bse_ospin[0], binary[j].bse_ospin[1], binary[j].bse_bcm_B[0], binary[j].bse_bcm_B[1], binary[j].bse_bcm_formation[0], binary[j].bse_bcm_formation[1], binary[j].bse_bacc[0], binary[j].bse_bacc[1], binary[j].bse_tacc[0], binary[j].bse_tacc[1], binary[j].bse_mass0[0], binary[j].bse_mass0[1], binary[j].bse_epoch[0], binary[j].bse_epoch[1]);
					}
				}
			}
			gzclose(snapfile);

#ifdef USE_MPI
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

}

/**
* @brief smaller snapshot outputting limited data.
*
* @param infile file name
*/
void print_denprof_snapshot(char* infile)
{
#ifdef USE_MPI
	//char filenm[150];
	//sprintf(filenm, "%s.%s", outprefix, "small.denprof.dat");

	//if(tcount == 1)
	fp_denprof = fopen(infile, "w");

	int i, j;
	double m, den;
	//fprintf(fp_denprof, "%ld\t", tcount);
	for(i=1; i<clus.N_MAX-MIN_CHUNK_SIZE; i+=MIN_CHUNK_SIZE)
	{
		m=0;
		for(j=i; j<i+MIN_CHUNK_SIZE; j++)
			m += star_m[j];

		den = m * madhoc / (4 * PI * ( cub(star_r[i+19]) - cub(star_r[i]) ) / 3);
		fprintf(fp_denprof, "%.16e\t%.16e\n", star_r[i+10], den);
	}			
	//fprintf(fp_denprof, "\n");
	fclose(fp_denprof);
#endif
}

/**
* @brief Does a miscellaneous set of things, primarily copying the data read from the input file into the star and binary data structures. Also allocates duplicate array for the parallel version. Also, initializes a few global variables and the random number generator.
*
* @param argc input arg count
* @param argv[] input arg list
* @param rng random number generator
*/
void get_star_data(int argc, char *argv[], gsl_rng *rng)
{
	/* print version information to log file */
	pararootfprintf(logfile, "** %s %s (%s) [%s] **\n", CMCPRETTYNAME, CMCVERSION, CMCNICK, CMCDATE);
#ifdef USE_MPI
    mpi_para_file_write(mpi_logfile_wrbuf, &mpi_logfile_len, &mpi_logfile_ofst_total, &mpi_logfile);
#endif

	/* initialize the Search_Grid r_grid */
	//If we use the GPU code, we dont need the SEARCH_GRID. So commenting it out
	/*        if (SEARCH_GRID) { //Parallelize Search Grid
				 r_grid= search_grid_initialize(SG_POWER_LAW_EXPONENT, \
				 SG_MATCH_AT_FRACTION, SG_STARSPERBIN, SG_PARTICLE_FRACTION);
				 };
	 */
	/* initialize the root finder algorithm */
    q_root = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

#ifdef DEBUGGING
	// create a new hash table for the star ids 
	star_ids= g_hash_table_new(g_int_hash, g_int_equal);
	// load a trace list 
	//load_id_table(star_ids, "trace_list");
#endif

	/* MPI: Allocate global arrays */
	mpiInitGlobArrays();

	/* Set up initial conditions */
	//MPI: This function populates the star and binary arrays with the cfd struct data.
	load_fits_file_data(); 

	/* set some important global variables */
	set_global_vars2();

	//Resets the global rng based on input seed. Does not alter the parallel rng.
	reset_rng_t113(IDUM);

	/* binary remainders */
	clus.N_MAX = clus.N_STAR;
	N_b = clus.N_BINARY;

	//MPI3: Check with Stefan if the sentinel is reqd.
	//star[clus.N_MAX+1].E = star[clus.N_MAX+1].J = 0.0;
}

/**
* @brief Allocate global arrays based on total number of stars
*/
void mpiInitGlobArrays()
{
#ifdef USE_MPI
	/*MPI: Allocating global/duplicated arrays that will be needed by all processors.*/
	star_r = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_m = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_phi = (double *) malloc(N_STAR_DIM * sizeof(double));
#endif
}

typedef struct{
	double s_Eescaped;
	double s_Jescaped;
	double s_Eintescaped;
	double s_Ebescaped;
	double s_TidalMassLoss;
	double s_Etidal;
    double s_Prev_Dt;

    long long s_mpi_logfile_len;
    long long s_mpi_escfile_len;
    long long s_mpi_binintfile_len;
    long long s_mpi_collisionfile_len;
    long long s_mpi_tidalcapturefile_len;
    long long s_mpi_semergedisruptfile_len;
    long long s_mpi_removestarfile_len;
    long long s_mpi_relaxationfile_len;
    long long s_mpi_pulsarfile_len;
    long long s_mpi_bhmergerfile_len;
    long long s_mpi_logfile_ofst_total;
    long long s_mpi_escfile_ofst_total;
    long long s_mpi_binaryfile_ofst_total;
    long long s_mpi_binintfile_ofst_total;
    long long s_mpi_collisionfile_ofst_total;
    long long s_mpi_tidalcapturefile_ofst_total;
    long long s_mpi_semergedisruptfile_ofst_total;
    long long s_mpi_removestarfile_ofst_total;
    long long s_mpi_relaxationfile_ofst_total;
    long long s_mpi_pulsarfile_ofst_total;
    long long s_mpi_bhmergerfile_ofst_total;

	double s_OldTidalMassLoss;
    double s_TidalMassLoss_old;
	long   s_N_bb;		
	long   s_N_bs;
	double s_E_bb;
	double s_E_bs;
	long   s_Echeck; 		
	int    s_se_file_counter; 	
	long   s_snap_num; 		
	long   s_bh_snap_num;
	long   s_StepCount; 		
	long   s_tcount;
	double s_TotalTime;                           
	long   s_newstarid;
    double s_cenma_m;
    double s_cenma_m_new;
    double s_cenma_e;
} restart_struct_t;

void save_global_vars(restart_struct_t *rest){
	rest->s_Eescaped                           =Eescaped;
	rest->s_Jescaped                           =Jescaped;
	rest->s_Eintescaped                        =Eintescaped;
	rest->s_Ebescaped                          =Ebescaped;
	rest->s_TidalMassLoss                      =TidalMassLoss;
	rest->s_Etidal                             =Etidal;
    rest->s_Prev_Dt                            =Prev_Dt;

	rest->s_mpi_logfile_len                    =mpi_logfile_len;
	rest->s_mpi_escfile_len                    =mpi_escfile_len;
	rest->s_mpi_binintfile_len                 =mpi_binintfile_len;
	rest->s_mpi_collisionfile_len              =mpi_collisionfile_len;
	rest->s_mpi_tidalcapturefile_len           =mpi_tidalcapturefile_len;
	rest->s_mpi_semergedisruptfile_len         =mpi_semergedisruptfile_len;
	rest->s_mpi_removestarfile_len             =mpi_removestarfile_len;
	rest->s_mpi_relaxationfile_len             =mpi_relaxationfile_len;
	rest->s_mpi_pulsarfile_len                 =mpi_pulsarfile_len;
	rest->s_mpi_bhmergerfile_len               =mpi_bhmergerfile_len;
	rest->s_mpi_logfile_ofst_total             =mpi_logfile_ofst_total;
	rest->s_mpi_escfile_ofst_total             =mpi_escfile_ofst_total;
	rest->s_mpi_binaryfile_ofst_total          =mpi_binaryfile_ofst_total;
	rest->s_mpi_binintfile_ofst_total          =mpi_binintfile_ofst_total;
	rest->s_mpi_collisionfile_ofst_total       =mpi_collisionfile_ofst_total;
	rest->s_mpi_tidalcapturefile_ofst_total    =mpi_tidalcapturefile_ofst_total;
	rest->s_mpi_semergedisruptfile_ofst_total  =mpi_semergedisruptfile_ofst_total;
	rest->s_mpi_removestarfile_ofst_total      =mpi_removestarfile_ofst_total;
	rest->s_mpi_relaxationfile_ofst_total      =mpi_relaxationfile_ofst_total;
	rest->s_mpi_pulsarfile_ofst_total          =mpi_pulsarfile_ofst_total;
	rest->s_mpi_bhmergerfile_ofst_total        =mpi_bhmergerfile_ofst_total;

    rest->s_OldTidalMassLoss                   =OldTidalMassLoss;
    rest->s_TidalMassLoss_old                  =TidalMassLoss_old;
	rest->s_N_bb		                       =N_bb;		
	rest->s_N_bs                               =N_bs;
	rest->s_E_bb                               =E_bb;
	rest->s_E_bs                               =E_bs;
	rest->s_Echeck 		                       =Echeck; 		
	rest->s_se_file_counter 	               =se_file_counter; 	
	rest->s_snap_num 		                   =snap_num; 		
	rest->s_bh_snap_num                        =bh_snap_num;
	rest->s_StepCount 		                   =StepCount; 		
	rest->s_tcount                             =tcount;
	rest->s_TotalTime                          =TotalTime;
	rest->s_newstarid                          =newstarid;
	rest->s_cenma_m                            =cenma.m;                 
	rest->s_cenma_m_new                        =cenma.m_new;                
	rest->s_cenma_e                            =cenma.E;                       
}

void load_global_vars(restart_struct_t *rest){
	Eescaped                           =rest->s_Eescaped;
	Jescaped                           =rest->s_Jescaped;
	Eintescaped                        =rest->s_Eintescaped;
	Ebescaped                          =rest->s_Ebescaped;
	TidalMassLoss                      =rest->s_TidalMassLoss;
	Etidal                             =rest->s_Etidal;
    Prev_Dt                            =rest->s_Prev_Dt;

	mpi_logfile_len                    =rest->s_mpi_logfile_len;
	mpi_escfile_len                    =rest->s_mpi_escfile_len;
	mpi_binintfile_len                 =rest->s_mpi_binintfile_len;
	mpi_collisionfile_len              =rest->s_mpi_collisionfile_len;
	mpi_tidalcapturefile_len           =rest->s_mpi_tidalcapturefile_len;
	mpi_semergedisruptfile_len         =rest->s_mpi_semergedisruptfile_len;
	mpi_removestarfile_len             =rest->s_mpi_removestarfile_len;
	mpi_relaxationfile_len             =rest->s_mpi_relaxationfile_len;
	mpi_pulsarfile_len                 =rest->s_mpi_pulsarfile_len;
	mpi_bhmergerfile_len               =rest->s_mpi_bhmergerfile_len;
	mpi_logfile_ofst_total             =rest->s_mpi_logfile_ofst_total;
	mpi_escfile_ofst_total             =rest->s_mpi_escfile_ofst_total;
	mpi_binaryfile_ofst_total          =rest->s_mpi_binaryfile_ofst_total;
	mpi_binintfile_ofst_total          =rest->s_mpi_binintfile_ofst_total;
	mpi_collisionfile_ofst_total       =rest->s_mpi_collisionfile_ofst_total;
	mpi_tidalcapturefile_ofst_total    =rest->s_mpi_tidalcapturefile_ofst_total;
	mpi_semergedisruptfile_ofst_total  =rest->s_mpi_semergedisruptfile_ofst_total;
	mpi_removestarfile_ofst_total      =rest->s_mpi_removestarfile_ofst_total;
	mpi_relaxationfile_ofst_total      =rest->s_mpi_relaxationfile_ofst_total;
	mpi_pulsarfile_ofst_total          =rest->s_mpi_pulsarfile_ofst_total;
	mpi_bhmergerfile_ofst_total        =rest->s_mpi_bhmergerfile_ofst_total;

    OldTidalMassLoss                   =rest->s_OldTidalMassLoss;
    TidalMassLoss_old                  =rest->s_TidalMassLoss_old;
	N_bb		                       =rest->s_N_bb;		
	N_bs                               =rest->s_N_bs;
	E_bb                               =rest->s_E_bb;
	E_bs                               =rest->s_E_bs;
	Echeck 		                       =rest->s_Echeck; 		
	se_file_counter 	               =rest->s_se_file_counter; 	
	snap_num 		                   =rest->s_snap_num; 		
	bh_snap_num                        =rest->s_bh_snap_num;
	StepCount 		                   =rest->s_StepCount; 		
	tcount                             =rest->s_tcount;
	TotalTime                          =rest->s_TotalTime;
	newstarid                          =rest->s_newstarid;
	cenma.m                            =rest->s_cenma_m;
	cenma.m_new                        =rest->s_cenma_m_new;
	cenma.E                            =rest->s_cenma_e;
}

void save_restart_file(){
	FILE *my_restart_file;
	char restart_file[200];
	char delete_file[200];
	char restart_folder[200];
	int i,j;
	int num_bin=0;
	struct stat folder_thing = {0};
	restart_struct_t restart_struct;
	
	sprintf(restart_folder, "./%s-RESTART", outprefix);
	sprintf(restart_file, "%s/%s.restart.%ld-%d.bin",restart_folder,outprefix,NEXT_RESTART,myid);

	/*If it's our first restart, we need to create the folder; if not, we can
	 * just save it in the existing folder*/
	if (stat(restart_folder,&folder_thing) == -1) {
		mkdir(restart_folder, 0700);
	}

	my_restart_file = fopen(restart_file,"wb");
	if (!my_restart_file){
		eprintf("can't open restart file %s for writing!\n",restart_file);
		exit_cleanly(-1, __FUNCTION__);
	}

	
	/*Save the entire star and binary arrays, including the many empty stars at
	 * the end; easier this way, and it ensures a bit-by-bit restart (also, the
	 * N_*_DIM_OPT are never updated, and should be the same as when the
	 * original star and binary structures were allocated)*/
	clus.N_BINARY = N_b;
	save_global_vars(&restart_struct);

	fwrite(curr_st, sizeof(struct rng_t113_state), 1, my_restart_file);
	fwrite(&restart_struct, sizeof(restart_struct_t), 1, my_restart_file);
	fwrite(&clus, sizeof(clus_struct_t), 1, my_restart_file);
	fwrite(star, sizeof(star_t), N_STAR_DIM_OPT, my_restart_file);
	fwrite(binary, sizeof(binary_t), N_BIN_DIM_OPT, my_restart_file);

	fclose(my_restart_file);

	/*Delete the last restart (or the last after however many we want to keep)*/
	long restart_to_delete = NEXT_RESTART - CHECKPOINTS_TO_KEEP;
	if ((restart_to_delete > 0) && (CHECKPOINTS_TO_KEEP != 0)){
		sprintf(delete_file, "%s/%s.restart.%ld-%d.bin",restart_folder,outprefix,restart_to_delete,myid);
		remove(delete_file);
	}

	rootprintf("******************************************************************************\n");
	rootprintf("Saving checkpoint %d at time %g in folder %s\n",NEXT_RESTART,TotalTime,restart_folder);
	rootprintf("******************************************************************************\n");
		
	NEXT_RESTART += 1;

}

void load_restart_file(){
	FILE *my_restart_file;
	char restart_file[200];
	char restart_folder[200];
	int i,j;
	int num_bin=0;
	struct stat folder_thing = {0};
	restart_struct_t restart_struct;

	sprintf(restart_folder, "./%s-RESTART", oldoutprefix);
	sprintf(restart_file, "%s/%s.restart.%ld-%d.bin",restart_folder,oldoutprefix,RESTART_TCOUNT,myid);

	if (stat(restart_folder,&folder_thing) == -1) {
		eprintf("can't find the restart folder %s\n",restart_folder);
		exit_cleanly(-1, __FUNCTION__);
	}

	my_restart_file = fopen(restart_file,"rb");
	if (!my_restart_file){
		eprintf("can't open restart file %s for writing!\n",restart_file);
		exit_cleanly(-1, __FUNCTION__);
	}

	/*These must be allocated here for the binary files to load correctly*/
	star = (star_t *) malloc(N_STAR_DIM_OPT*sizeof(star_t));
	binary = (binary_t *) malloc(N_BIN_DIM_OPT*sizeof(binary_t));
	curr_st = (struct rng_t113_state*) malloc(sizeof(struct rng_t113_state));

	/*Set the units using the original data from the fits file*/
	units_set();
	
	/*Load the entire star and binaries arrays at once.  Because this is done in
	 * a single chunk of memory and with the same size of arrays as was
	 * generated from the FITS file, this should load the exact local state into
	 * each file*/
	fread(curr_st, sizeof(struct rng_t113_state), 1, my_restart_file);
	fread(&restart_struct, sizeof(restart_struct_t), 1, my_restart_file);
	fread(&clus, sizeof(clus_struct_t), 1, my_restart_file);
	fread(star, sizeof(star_t), N_STAR_DIM_OPT, my_restart_file);
	fread(binary, sizeof(binary_t), N_BIN_DIM_OPT, my_restart_file);

	fclose(my_restart_file);

	/*Set the random number generator back where it was*/
	set_rng_t113(*curr_st);

	/*Load the global variables from the structure (replaces set_global_vars 1
	 * and 2 in the code)*/
	load_global_vars(&restart_struct);

	/*Set this and allocate the global arrays*/
	N_b = clus.N_BINARY;
    q_root = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	mpiInitGlobArrays();

	/*Set the MPI files back to exactly where they were, using the saved offsets*/
	MPI_File_seek(mpi_logfile,mpi_logfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_binintfile,mpi_binintfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_escfile,mpi_escfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_collisionfile,mpi_collisionfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_tidalcapturefile,mpi_tidalcapturefile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_semergedisruptfile,mpi_semergedisruptfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_relaxationfile,mpi_relaxationfile_ofst_total,MPI_SEEK_SET);
	MPI_File_seek(mpi_removestarfile,mpi_removestarfile_ofst_total,MPI_SEEK_SET);
    if (THREEBODYBINARIES){
        MPI_File_seek(mpi_threebbfile,mpi_threebbfile_ofst_total,MPI_SEEK_SET);
        MPI_File_seek(mpi_threebbprobabilityfile,mpi_threebbprobabilityfile_ofst_total,MPI_SEEK_SET);
        MPI_File_seek(mpi_lightcollisionfile,mpi_lightcollisionfile_ofst_total,MPI_SEEK_SET);
        MPI_File_seek(mpi_threebbdebugfile,mpi_threebbdebugfile_ofst_total,MPI_SEEK_SET);
    }
    if (WRITE_BH_INFO) {
        MPI_File_seek(mpi_newbhfile,mpi_newbhfile_ofst_total,MPI_SEEK_SET);
        MPI_File_seek(mpi_bhmergerfile,mpi_bhmergerfile_ofst_total,MPI_SEEK_SET);
    }
	if(WRITE_PULSAR_INFO){
		 MPI_File_seek(mpi_pulsarfile,mpi_pulsarfile_ofst_total,MPI_SEEK_SET);
	}

	next_restart_t = CHECKPOINT_INTERVAL;
	NEXT_RESTART = RESTART_TCOUNT + 1;

	/*Bit of arrray-based housekeeping*/
	star_r[0] = ZERO; 
	star_r[clus.N_STAR + 1] = SF_INFINITY;
	cenma.m_new= star_m[0];
	star_m[0] = 0.0;

	rootprintf("******************************************************************************\n");
	rootprintf("Loading checkpoint %d at time %g in folder %s\n",RESTART_TCOUNT,TotalTime,restart_folder);
	rootprintf("******************************************************************************\n");

}
/*
void alloc_bin_buf()
{
#ifdef USE_MPI
	size_bin_buf = N_BIN_DIM / procs; //this is much more than what would be practically reqd. worst case scenario is N_BIN_DIM.

	//OPT: Modularize by declaring variables in main fn, and pass as argument to reqd functions.
	//the local binary array
	binary_buf = (binary_t *) calloc( size_bin_buf, sizeof(binary_t) );

	//array to store the indices of the local binary members in the original binary array
	num_bin_buf = (int *) calloc( size_bin_buf, sizeof(int) );
#endif
}
*/

/*
void mpi_merge_files()
{
#ifdef USE_MPI
	int i;
	char file_ext[64];
	MPI_Barrier(MPI_COMM_WORLD);

	// Total of 22 files. Variable no.of files for lagrad***.dat, core.dat
	// Using a simple brute force parallelization.
	if(0 % procs == myid)
		cat_and_rm_files("cmc.parsed");
	if(1 % procs == myid)
		cat_and_rm_files("lagrad.dat");
	if(2 % procs == myid)
		cat_and_rm_files("dyn.dat");
	if(3 % procs == myid)
		cat_and_rm_files("lagrad_10_info.dat");
	if(4 % procs == myid)
		cat_and_rm_files("avemass_lagrad.dat");
	if(5 % procs == myid)
		cat_and_rm_files("nostar_lagrad.dat");
	if(6 % procs == myid)
		cat_and_rm_files("rho_lagrad.dat");
	if(7 % procs == myid)
		cat_and_rm_files("ke_rad_lagrad.dat");
	if(8 % procs == myid)
		cat_and_rm_files("ke_tan_lagrad.dat");
	if(9 % procs == myid)
		cat_and_rm_files("v2_rad_lagrad.dat");
	if(10 % procs == myid)
		cat_and_rm_files("v2_tan_lagrad.dat");
	if(11 % procs == myid)
		cat_and_rm_files("centmass.dat");
	if(12 % procs == myid)
		cat_and_rm_files("log");
	if(13 % procs == myid)
		cat_and_rm_files("bin.dat");
	if(14 % procs == myid)
		cat_and_rm_files("binint.log");
	if(15 % procs == myid)
		cat_and_rm_files("esc.dat");
	if(16 % procs == myid)
		cat_and_rm_files("collision.log");
	if(17 % procs == myid)
		cat_and_rm_files("tidalcapture.log");
	if(18 % procs == myid)
		cat_and_rm_files("semergedisrupt.log");
	if(19 % procs == myid)
		cat_and_rm_files("removestar.log");
	if(20 % procs == myid)
		cat_and_rm_files("relaxation.dat");
	if(21 % procs == myid)
		for(i=0; i<NO_MASS_BINS-1; i++){
			sprintf(file_ext, "lagrad%d-%g-%g.dat", i, mass_bins[i], mass_bins[i+1]);
			cat_and_rm_files(file_ext);
		}
	if(22 % procs == myid)
		if (WRITE_EXTRA_CORE_INFO)
			cat_and_rm_files("core.dat");
	//MPI2: Is barrier needed? Probably not.
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}
*/
/*
void cat_and_rm_files(char* file_ext)
{
	char cmd[150];
	int i;

	for( i = 0; i < procs; ++i )
	{
		sprintf(cmd, "cat %s%d.%s >> %s.%s", outprefix_bak, i, file_ext, outprefix_bak, file_ext);
		dprintf("command is %s\n", cmd);
		system( cmd );
		sprintf(cmd, "rm %s%d.%s", outprefix_bak, i, file_ext);
		system( cmd );
	}
	dprintf("MPI Files merging: lag file = %s.%s\n", outprefix_bak, file_ext);
}

void cat_and_rm_files(char* file_ext)
{
   //char *filename_buf, temp[150];
	char* filename_buf = (char*)malloc(40 * procs * sizeof(char)); //40 is assumed to be the worst case length of a filename including outprefix.
	char* cmd = (char*)malloc(40 * procs * sizeof(char));
   int i;

   sprintf(filename_buf, "%s0.%s ", outprefix_bak, file_ext);
   for( i = 1; i < procs; ++i )
      sprintf( filename_buf, "%s%s%d.%s ", filename_buf, outprefix_bak, i, file_ext );

   sprintf( cmd, "cat %s> %s.%s", filename_buf, outprefix_bak, file_ext);
   system( cmd );
   dprintf("MPI Files merging: file = %s\n", filename_buf);

   for( i = 0; i < procs; ++i )
   {
      sprintf( cmd, "rm %s%d.%s", outprefix_bak, i, file_ext);
      system( cmd );
   }
	//free(filename_buf);
	//free(cmd);
}
*/
/*
void save_root_files()
{
#ifdef USE_MPI
	int i;
	char file_ext[64];
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid==0)
	{
		save_root_files_helper("cmc.parsed");
		save_root_files_helper("lagrad.dat");
		save_root_files_helper("dyn.dat");
		save_root_files_helper("lagrad_10_info.dat");
		save_root_files_helper("avemass_lagrad.dat");
		save_root_files_helper("nostar_lagrad.dat");
		save_root_files_helper("rho_lagrad.dat");
		save_root_files_helper("ke_rad_lagrad.dat");
		save_root_files_helper("ke_tan_lagrad.dat");
		save_root_files_helper("v2_rad_lagrad.dat");
		save_root_files_helper("v2_tan_lagrad.dat");
		save_root_files_helper("centmass.dat");
		save_root_files_helper("log");
		save_root_files_helper("bin.dat");
		save_root_files_helper("binint.log");
		save_root_files_helper("esc.dat");
		save_root_files_helper("collision.log");
		save_root_files_helper("tidalcapture.log		save_root_files_helper("semergedisrupt.log");
		save_root_files_helper("removestar.log");
		save_root_files_helper("relaxation.dat");
		for(i=0; i<NO_MASS_BINS-1; i++){
			sprintf(file_ext, "lagrad%d-%g-%g.dat", i, mass_bins[i], mass_bins[i+1]);
			save_root_files_helper(file_ext);
		}
		if (WRITE_EXTRA_CORE_INFO)
			save_root_files_helper("core.dat");
		//MPI2: Is barrier needed? Probably not.
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}
*/
/*
void save_root_files_helper(char* file_ext)
{
   char cmd[150];

   sprintf(cmd, "cp %s0.%s %s.%s", outprefix_bak, file_ext, outprefix_bak, file_ext);
	system( cmd );
*/
/*
   int i;
   for( i = 0; i < procs; ++i )
   {
      sprintf( cmd, "rm %s%d.%s &", outprefix_bak, i, file_ext);
      system( cmd );
   }
}
*/
/*
void rm_files()
{
#ifdef USE_MPI
   int i;
   char file_ext[64];

   rm_files_helper("cmc.parsed");
   rm_files_helper("lagrad.dat");
   rm_files_helper("dyn.dat");
   rm_files_helper("lagrad_10_info.dat");
   rm_files_helper("avemass_lagrad.dat");
   rm_files_helper("nostar_lagrad.dat");
   rm_files_helper("rho_lagrad.dat");
   rm_files_helper("ke_rad_lagrad.dat");
   rm_files_helper("ke_tan_lagrad.dat");
   rm_files_helper("v2_rad_lagrad.dat");
   rm_files_helper("v2_tan_lagrad.dat");
   rm_files_helper("centmass.dat");
   rm_files_helper("log");
   rm_files_helper("bin.dat");
   rm_files_helper("binint.log");
   rm_files_helper("esc.dat");
   rm_files_helper("collision.log");
   rm_files_helper("tidalcapture.log");
   rm_files_helper("semergedisrupt.log");
   rm_files_helper("removestar.log");
   rm_files_helper("relaxation.dat");
   for(i=0; i<NO_MASS_BINS-1; i++){
      sprintf(file_ext, "lagrad%d-%g-%g.dat", i, mass_bins[i], mass_bins[i+1]);
      rm_files_helper(file_ext);
   }
   if (WRITE_EXTRA_CORE_INFO)
      rm_files_helper("core.dat");

   //MPI2: Is barrier needed? Probably not.
   MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void rm_files_helper(char* file_ext)
{
   char cmd[150];

   sprintf(cmd, "rm %s%d.%s", outprefix_bak, myid, file_ext);
   system( cmd );
}
*/
/*
void print_small_output()
{
#ifdef USE_MPI
	int i;
	char filenm1[150], filenm2[150];

	if(myid==0)
	{
		if (tcount == 1)
		{
			sprintf(filenm1, "%s.%s", outprefix_bak, "small.lagrad.dat");
			sprintf(filenm2, "%s.%s", outprefix_bak, "small.log.dat");

			fp_lagrad = fopen(filenm1, "w");
			fp_log = fopen(filenm2, "w");
			fprintf(fp_lagrad, "# Lagrange radii [code units]\n");
			fprintf(fp_log, "# TotalTime\tDt\ttcount\tclus.N_MAX\tMtotal\tEtotal.tot\n");
		}

		fprintf(fp_lagrad, "%ld\t%.16e\t", tcount, TotalTime);
		for (i=0; i<MASS_PC_COUNT; i++)
			fprintf(fp_lagrad, "%e\t", mass_r[i]);
		fprintf(fp_lagrad, "\n");
		
		fprintf(fp_log, "%.8g\t%.8g\t%ld\t%ld\t%.8g\t%.8g\n", TotalTime, Dt, tcount, clus.N_MAX, Mtotal, Etotal.tot);
	}
#endif
}
*/
