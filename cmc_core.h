
/** 
 * @brief Struct holding the indices and estimates of the local density.
 */
struct densities {
  long *idx;
  double *rho;
  long nave;
};

/** 
 * @brief Struct holding the properties of the cluster core.
 */
struct core_t {
  double rho;
  double v_rms;
  double r;
  double r_spitzer;
  double m_ave;
  double n;
  long N;
  double Trc;
};

void append_core_header(FILE *cfile, char *tag, int core_number);
void write_core_data(FILE *cfile, struct core_t);
struct core_t core_properties(struct densities rhoj);
struct densities density_estimators(int n_point, int *startypes, int len);
struct core_t no_remnants_core(int n_points);

