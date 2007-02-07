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


struct Search_Grid search_grid_initialize(double power_law_exponent);
  /* This is somewhat trivial at the moment but makes the usage of the 
   * Search_Grid structure clearer. */
  struct Search_Grid grid;
  int i;

  for (i=0; i<1000; i++) {
    grid.star_index[i]= 0;
  };
  grid.power_law_exponent= power_law_exponent;

  return(grid);
};

void search_grid_update(struct Search_Grid *grid) {
  int i, r_index;
  double coeff, r_max, r_min;
  
  r_max= star[clus.N_MAX+1].r;
  r_min= star[1].r;

  grid->interpol_coeff= (r_max-r_min)/pow(1000., grid->power_law_exponent);
  r_index=0;
  r_at_index= search_grid_get_r(grid, r_index);
  for (i=1; i<clus.N_MAX+1; i++) {
    if (star[i].r <= r_at_index) {
      grid->radius[r_index]= i;
    } else {
      do {
        r_index++;
        r_at_index= search_grid_get_r(grid, r_index);
        grid->radius[r_index]= i;
      } while (star[i].r> r_at_index);
    };
  };
};

double search_grid_get_r(struct Search_Grid *grid, int index) {
  if (index> 999) {
    printf("Index is too large! It must not be larger than 999, but got %i", index);
    exit(1);
  };
  return(star[1].r+ grid->interpol_coeff*pow(index+1., grid->power_law_exponent));
};

int search_grid_get_interval(struct Search_Grid *grid, double r) {
  double r_to_n, n;

  n= grid->power_law_exponent;
  r_to_n= (r-star[1].r)/grid->interpol_coeff;
  
  return(floor(pow(r_to_n, 1./n)));
};

