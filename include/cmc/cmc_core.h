/* vi: set filetype=c.doxygen: */

/** 
 * @brief Struct holding the indices and estimates of the local density.
 */
struct densities {
/**
* @brief ?
*/
  long *idx;
/**
* @brief ?
*/
  double *rho;
/**
* @brief ?
*/
  long nave;
};

/** 
 * @brief Struct holding the properties of the cluster core.
 */
struct core_t {
/**
* @brief ?
*/
  double rho;
/**
* @brief ?
*/
  double v_rms;
/**
* @brief ?
*/
  double r;
/**
* @brief ?
*/
  double r_spitzer;
/**
* @brief ?
*/
  double m_ave;
/**
* @brief ?
*/
  double n;
/**
* @brief ?
*/
  long N;
/**
* @brief ?
*/
  double Trc;
};

void append_core_header(FILE *cfile, char *tag, int core_number);
void write_core_data(FILE *cfile, struct core_t);
struct core_t core_properties(struct densities rhoj);
struct densities density_estimators(int n_point, int *startypes, int len);
struct core_t no_remnants_core(int n_points);

