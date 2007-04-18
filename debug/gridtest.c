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


void main (void) {
  double r[100000000];
  long i;
  struct Search_Grid *grid;
  struct Interval search_interval;
  gsl_rng *rng;
  const gsl_rng_type *rng_type=gsl_rng_mt19937;

  rng= gsl_rng_alloc(rng_type);

  /* some initialization stuff needed to read the fits file */
  MEGA_YEAR= 1e-3;
  SOLAR_MASS_DYN= 1;
  GAMMA=0.01;
  N_STAR_DIM= 100000;
  N_STAR_DIM = (long) floor(1.5 * ((double) N_STAR_DIM));
  printf("#NSTAR_DIM is: %li\n", N_STAR_DIM);
  star = calloc(N_STAR_DIM, sizeof(star_t));
  clus.N_MAX= 100000;
  clus.N_STAR=100000;
  /* end of weirdness */
  read_fits_file_data("pltest");
  i=1;
  printf("#length and mass are %f %f\n", units.l, units.m);
  grid= search_grid_initialize(0.5, 0.6, 10, 0.97);
  search_grid_update(grid);
  printf("#grid power law is %f\n", grid->power_law_exponent);
/*  
  for (i=0; i<grid.length; i++) {
    r[i]= search_grid_get_r(&grid, i);
    printf("%li, %f, %li\n", i, r[i], grid.radius[i]);
  }; 

  printf("\n");
*/
/*  for (i=0; i< grid.length; i++) {
    double r;

    r= search_grid_get_r(&grid, i);
    printf("%li, %f, %li\n",i, r, i-search_grid_get_grid_index(&grid, r));
  };*/
/*  for (i=0; i<grid.length; i++) {
    printf("%li %li %li\n", i, grid.radius[i], grid.radius[i+1]-grid.radius[i]);
  };
 
  printf("\n");
*//*
  for (i=0; i<100000; i++) {
    printf("%li, %g\n", i, star[i].r);
  };
*//*
  for (i=0; i<clus.N_MAX/100+1; i++) {
    printf("%li %f\n", i, star[i*100].r);
  };
*/
/* for (i=0; i<1000000; i++) {
    long min, max, last_bin, ind2;
    double ind;
    //r[i]= gsl_rng_uniform(rng)*star[(long)(clus.N_STAR)].r+star[1].r;
    r[i]= i*0.01+star[1].r;
    search_interval= search_grid_get_interval(grid, r[i]);
    ind= search_grid_get_grid_indexf(grid, r[i]);
    ind2= search_grid_get_grid_index(grid, r[i]);
    last_bin= grid.radius[(grid.length)-1];
    min= search_interval.min;
    max= search_interval.max;
    if (star[min].r>r[i] || star[max].r<=r[i]) {
      printf("Harrgh! %f is not between %f and %f (%li, %li, %li, %f)\n", r[i], star[(int)min].r, \
          star[(int)max].r, min, max, search_grid_get_grid_index(grid, r[i]), ind);
      printf("radius[grid_index-1]= %li, radius[grid_index]= %li\n", grid.radius[ind2-1], grid.radius[ind2]);
      printf("radii: %f, %f\n", star[grid.radius[ind2-1]].r, star[grid.radius[ind2]].r);
      //printf("radius[index] and radius[index-1] %li, %li\n", grid.radius[(long)ind], grid.radius[(long)ind-1]);
    };
  }; */
  total_bisections= 0;
  search_interval.min=1;
  search_interval.max= clus.N_MAX+1;
  for (i=0; i<100000000; i++) {
    long min, max, last_bin, ind2;
    double ind;
    //r[i]= gsl_rng_uniform(rng)*star[(long)(clus.N_STAR)].r+star[1].r;
    r[i]= i*star[clus.N_MAX-1].r*0.95/1.e8+star[1].r;
    search_interval= search_grid_get_interval(grid, r[i]);
    min= search_interval.min;
    max= search_interval.max;
    ind2= FindZero_r(min, max, r[i]);
  };
    printf("Number of bisections is %li\n", total_bisections);
    search_grid_free(grid);
};


