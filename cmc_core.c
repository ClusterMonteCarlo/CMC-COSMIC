/* Various subroutines that caluculate quantities in the core including its 
 * extent according to various definitions */
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief ?
*
* @param k ?
* @param set ?
* @param len ?
*
* @return ?
*/
static inline int is_member(int k, int* set, int len) {
  int member;
  int i;

  /* could use bisection */
  member= 0;
  for (i=0; i<len; i++) {
    if (k==set[i]) {
      member=1;
      break;
    }
  }
  return(member);
}

/** 
 * @brief Calculates density estimators including only stars with certain types.
 * 
 *
 * @param n_points the number of points (stars) over which the local density 
 * is estimated (Casertano & Hut '85 suggest n_points=6)
 *
 * @param startypes Array containing the star types that should be considered.
 *
 * @param len The length of the startype array.
 * 
 * @return The array with the estimated local densities. The array is 
 * malloc'ed so you have to free it yourself.
 */
struct densities density_estimators(int n_points, int *startypes, int len) {
  long i, nave, ibuf;
  double m, m_tot, mrho, Vrj;
  struct densities rhoj;
  long jmin, jmax, j;

  /* calculate the total mass of all stars excluding remnants */
  m_tot= 0.;
  for (i=1; i< clus.N_MAX+1; i++) {
    if (is_member(star[i].se_k, startypes, len)) {
      m_tot+= star[i].m*madhoc;
    }
  }

  rhoj.idx= calloc(80, sizeof(long));

  /* .. and now the number of non-remnants within their half-mass radius */
  nave= 0; m= 0; i=1; ibuf=0;
  while (m< 0.5*m_tot) {
    if (is_member(star[i].se_k, startypes, len)) {
      m+= star[i].m*madhoc;
      rhoj.idx[nave]= i;
      nave++;
      ibuf++;
      if (ibuf==80) {
        rhoj.idx= (long *) realloc(rhoj.idx, sizeof(long)*(nave+80));
        ibuf= 0;
      }
    }
    i++;
  }


  rhoj.rho = (double *) calloc(nave, sizeof(double));
  rhoj.nave= nave;

  for (i=0; i<nave; i++) {
    jmin= MAX(i-n_points/2, 0);
    jmax= MIN(jmin + n_points, nave-1);
    mrho= 0.;
    /* this is equivalent to their J-1 factor for the case of equal masses,
       and seems like a good generalization for unequal masses */
    for (j=jmin+1; j<= jmax-1; j++) {
      mrho+= star[rhoj.idx[j]].m * madhoc;
    }
    Vrj = 4.0/3.0  * PI * (fb_cub(star[rhoj.idx[jmax]].r) - fb_cub(star[rhoj.idx[jmin]].r));

    rhoj.rho[i]= mrho/Vrj;
  }

  return(rhoj);
}

#ifdef USE_MPI
/**
 * @brief Calculates density estimators including only stars with certain types. (mpi version of density_estimators)
 *
 *
 * @param n_points the number of points (stars) over which the local density
 * is estimated (Casertano & Hut '85 suggest n_points=6)
 *
 * @param startypes Array containing the star types that should be considered.
 *
 * @param len The length of the startype array.
 *
 * @return The array with the estimated local densities. The array is
 * malloc'ed so you have to free it yourself.
 */
struct densities mpi_density_estimators(int n_points, int *startypes, int len) {
  long i, nave, ibuf;
  double m, m_tot, mrho, Vrj;
  struct densities rhoj;
  long jmin, jmax, j;
  long g_i;
  double m_cum;

  /* calculate the total mass of all stars excluding remnants */
  //MPI: First, find the local total mass.
  m_tot= 0.;
  for (i=1; i<= clus.N_MAX_NEW; i++) {
    if (is_member(star[i].se_k, startypes, len)) {
        g_i = get_global_idx(i);
        m_tot+= star_m[g_i]*madhoc;
    }
  }

  //MPI: find cumulative total mass in processors before this one. Also, find global total mass.
  m_cum=0.0;
  double tmpTimeStart = timeStartSimple();
  MPI_Exscan(&m_tot, &m_cum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  timeEndSimple(tmpTimeStart, &t_comm);


  long* rhoj_idx= calloc(80, sizeof(long));

  /* .. and now the number of non-remnants within their half-mass radius */
  //MPI: Fill up the local rho idx array
  m=m_cum;
  nave= 0; i=1; ibuf=0;
  while (m< 0.5*m_tot && i<=clus.N_MAX_NEW) {
    if (is_member(star[i].se_k, startypes, len)) {
        g_i = get_global_idx(i);
        m+= star_m[g_i]*madhoc;
        rhoj_idx[nave]= i;
        nave++;
        ibuf++;
        if (ibuf==80) {
            rhoj_idx= (long *) realloc(rhoj_idx, sizeof(long)*(nave+80));
            ibuf= 0;
        }
    }
    i++;
  }

  //MPI: Now, the idx array has to be corrected based on the last index of the previous processor.
  long idx_cum = 0;
  long idx_cum_s = 0;
  if(nave > 0)
      idx_cum_s = rhoj_idx[nave-1];

  tmpTimeStart = timeStartSimple();
  MPI_Exscan(&idx_cum_s, &idx_cum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  timeEndSimple(tmpTimeStart, &t_comm);

  for(i=0; i<nave; i++)
  {
      rhoj_idx[i] += idx_cum;
  }

  //MPI: The simplest way I could think of parallelizing this is by collecting the entire rho idx array. This is definitely inefficient and not scalable, but as of now we need a working solution. So, what follows is basically this - collecting chunks of the idx array from all the processors, putting them together and sending a copy to all nodes.
  long* nave_all_long = (long*) malloc(procs * sizeof(long));
  int* nave_all = (int*) malloc(procs * sizeof(long));
  int* nave_disp = (int*) calloc(procs, sizeof(int));

  //MPI: Gathering the size of the idx array from all processors .
  tmpTimeStart = timeStartSimple();
  MPI_Allgather(&nave, 1, MPI_LONG, nave_all_long, 1, MPI_LONG, MPI_COMM_WORLD);
  timeEndSimple(tmpTimeStart, &t_comm);

  //MPI: Calculating quantities required for the mpi call - displacement and count
  long total_nave = 0;
  int cum_disp=0;
  for(i=0; i<procs; i++)
  {
      nave_all[i] = nave_all_long[i];
      nave_disp[i] = cum_disp;
      cum_disp += nave_all[i];
  }
  total_nave = cum_disp;

  //MPI: Gathering all the idx chunks and sending the entire array to all processors.
  rhoj.idx = (long*) malloc(total_nave * sizeof(long));
  tmpTimeStart = timeStartSimple();
  MPI_Allgatherv(rhoj_idx, nave, MPI_LONG, rhoj.idx, nave_all, nave_disp, MPI_LONG, MPI_COMM_WORLD);
  timeEndSimple(tmpTimeStart, &t_comm);


  //MPI: Now we can go back and do the computations independently on each node.
  rhoj.rho = (double *) calloc(total_nave, sizeof(double));
  rhoj.nave= total_nave;

  for (i=0; i<rhoj.nave; i++) {
    jmin= MAX(i-n_points/2, 0);
    jmax= MIN(jmin + n_points, rhoj.nave-1);
    mrho= 0.;
    /* this is equivalent to their J-1 factor for the case of equal masses,
       and seems like a good generalization for unequal masses */
    for (j=jmin+1; j<= jmax-1; j++) {
      mrho+= star_m[rhoj.idx[j]] * madhoc;
    }
    Vrj = 4.0/3.0  * PI * (fb_cub(star_r[rhoj.idx[jmax]]) - fb_cub(star_r[rhoj.idx[jmin]]));

    rhoj.rho[i]= mrho/Vrj;
  }

  return(rhoj);
}

/**
* @brief calculate core quantities using density weighted averages (note that in Casertano & Hut (1985) only rho and rc are analyzed and tested) (mpi version of core_properties)
*
* @param rhoj ?
*
* @return ?
*/
struct core_t mpi_core_properties(struct densities rhoj) {
  double rhojsum;
  double rhoj2sum;
  struct core_t core;
  long idx, i;

  rhojsum = 0.0;
  rhoj2sum = 0.0;
  core.rho = 0.0;
  core.v_rms = 0.0;
  core.r = 0.0;
  core.m_ave = 0.0;
  for (i=0; i<rhoj.nave; i++) {
    idx= rhoj.idx[i];
    rhojsum += rhoj.rho[i];
    rhoj2sum += sqr(rhoj.rho[i]);
    core.rho += sqr(rhoj.rho[i]);
    core.r += rhoj.rho[i] * star_r[idx];
    core.m_ave += rhoj.rho[i] * star_m[idx] * madhoc;
    if( idx >= mpiBegin && idx <= mpiEnd )
        core.v_rms += rhoj.rho[i] * (sqr(star[get_local_idx(idx)].vr) + sqr(star[get_local_idx(idx)].vt));
  }

  //MPI: Reducing only vrms, the other values dont have to be reduced.
  double v_rms = core.v_rms;
  double tmpTimeStart = timeStartSimple();
  MPI_Reduce(&core.v_rms, &v_rms, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  timeEndSimple(tmpTimeStart, &t_comm);
  core.v_rms = v_rms;

  core.rho /= rhojsum;
  /* correction for inherent bias in estimator */
  core.rho *= 4.0/5.0;
  core.v_rms /= rhojsum;
  core.v_rms = sqrt(core.v_rms);
  core.r /= rhojsum;
  core.m_ave /= rhojsum;

  /* quantities derived from averages */
  core.n = core.rho / core.m_ave;
  core.r_spitzer = sqrt(3.0 * sqr(core.v_rms) / (4.0 * PI * core.rho));

  core.N = 4.0 / 3.0 * PI * cub(core.r) * core.n;
  /* core relaxation time, Spitzer (1987) eq. (2-62) */
  core.Trc = 0.065 * cub(core.v_rms) / (core.rho * core.m_ave);

  return(core);
}


#endif

/**
* @brief calculate core quantities using density weighted averages (note that in Casertano & Hut (1985) only rho and rc are analyzed and tested)
*
* @param rhoj ?
*
* @return ?
*/
struct core_t core_properties(struct densities rhoj) {
  double rhojsum;
  double rhoj2sum;
  struct core_t core;
  long idx, i;

  rhojsum = 0.0;
  rhoj2sum = 0.0;
  core.rho = 0.0;
  core.v_rms = 0.0;
  core.r = 0.0;
  core.m_ave = 0.0;
  for (i=0; i<rhoj.nave; i++) {
    idx= rhoj.idx[i];
    rhojsum += rhoj.rho[i];
    rhoj2sum += sqr(rhoj.rho[i]);
    core.rho += sqr(rhoj.rho[i]);
    core.v_rms += rhoj.rho[i] * (sqr(star[idx].vr) + sqr(star[idx].vt));
    core.r += rhoj.rho[i] * star[idx].r;
    core.m_ave += rhoj.rho[i] * star[idx].m * madhoc;
  }
  core.rho /= rhojsum;
  /* correction for inherent bias in estimator */
  core.rho *= 4.0/5.0;
  core.v_rms /= rhojsum;
  core.v_rms = sqrt(core.v_rms);
  core.r /= rhojsum;
  core.m_ave /= rhojsum;

  /* quantities derived from averages */
  core.n = core.rho / core.m_ave;
  core.r_spitzer = sqrt(3.0 * sqr(core.v_rms) / (4.0 * PI * core.rho));

  core.N = 4.0 / 3.0 * PI * cub(core.r) * core.n;
  /* core relaxation time, Spitzer (1987) eq. (2-62) */
  core.Trc = 0.065 * cub(core.v_rms) / (core.rho * core.m_ave);

  return(core);
}

/** 
 * @brief Appends the header information for a certain core to the core file.
 * 
 * No newline character will be added, so you have to do it yourself once 
 * all core headers are written.
 *
 * @param cfile       The core file.
 * @param tag         Shorthand for the type of core.
 * @param core_number The position of the core within the file (used to 
 *                    calculate the column numbers)
 */
void append_core_header(FILE *cfile, char *tag, int core_number) {
  int column;

  column= core_number*8+2;
  if (core_number==0) {
    rootfprintf(cfile, "# ");
    rootfprintf(cfile, "1:time");
  }
  rootfprintf(cfile, " %i:rho_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:v_rms_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:rc_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:r_spitzer_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:m_ave_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:n_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:N_%s", column, tag);
  column++;
  rootfprintf(cfile, " %i:Trc_%s ", column, tag);
}

/** 
 * @brief Writes out all core data to a file.
 *
 * No newline character will be added, so you have to do it yourself once
 * all core data is written.
 * 
 * @param cfile The core file (open file handle).
 * @param core  The struct that holds all core data.
 */

void write_core_data(FILE *cfile, struct core_t core) {
  rootfprintf(cfile, "%.12g %.12g %.12g %.12g %.12g %.12g %.12g %li %.12g ", TotalTime, core.rho, core.v_rms, core.r, 
      core.r_spitzer, core.m_ave, core.n, core.N, core.Trc);
};

/* Some convenience functions */
struct core_t no_remnants_core(int n_points) {
  struct densities rhoj;
  struct core_t core;
  int types[11]= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

#ifdef USE_MPI
  rhoj= mpi_density_estimators(n_points, types, 11);
  core= mpi_core_properties(rhoj);
#else
  rhoj= density_estimators(n_points, types, 11);
  core= core_properties(rhoj);
#endif
  free(rhoj.idx);
  free(rhoj.rho);
  return(core);
}

