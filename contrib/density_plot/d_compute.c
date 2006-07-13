#include <ruby.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

VALUE compute_class;

static VALUE c_do(VALUE self, VALUE j, VALUE p)
{
  // Allocate memory for pointers, ruby arrays
  double *rad;
  double *rho;
  double *rhoL;
  // Get internal variables and malloc
  int kount = NUM2INT(j);
  double rad_m = NUM2DBL(rb_iv_get(self, "@rRad_m"));
  double rratio = NUM2DBL(rb_iv_get(self, "@rRad_ratio"));
  double unit_l = NUM2DBL(rb_iv_get(self, "@unit_l"));
  rad = malloc(j*sizeof(double));
  rho = malloc(j*sizeof(double));
  rhoL = malloc(j*sizeof(double));
  double  zminus, z;
  double  sigma, sigmaL;
  // File for writing (temp.dat, to be renamed later)
  FILE* write;
  write = fopen("temp.dat", "w");
  int k;
  // Assign variables with loop
  VALUE local;
  local = rb_iv_get(self, "@rad");
  for(k=0;k<=kount;k++){
    rad[k] = NUM2DBL(rb_ary_entry(local,k));
     }
 local = rb_iv_get(self, "@rho");
  for(k=0;k<kount;k++){
    rho[k] = NUM2DBL(rb_ary_entry(local,k));
  }
  local = rb_iv_get(self, "@rhoL");
  for(k=0;k<kount;k++){
    rhoL[k] = NUM2DBL(rb_ary_entry(local, k));
  }
  int i;
  double rr;
  double reset = kount;
  // Perform integration
  int pkount = NUM2INT(p);
  for(i=0;i<=pkount;i++){
    kount = reset;
    sigma = 0.0;
    sigmaL = 0.0;
    rr = rad_m * pow(10,(((double) i)/(double) pkount)*log10(rratio));
    while (kount >= 0 && (rad[kount] >= rr)){
      z = sqrt(rad[kount]*rad[kount] - rr*rr);
      if (kount == 0 || (rad[kount-1] < rr)){
	    zminus = 0.0;
	    } else
	    {
	      zminus = sqrt(rad[kount-1]*rad[kount-1] - rr*rr);
	    }
    sigma += (rho[kount] * (z - zminus));
    sigmaL += (rhoL[kount] * (z - zminus));
    kount--;
    }
    sigma = sigma * 2.0;
    sigmaL = sigmaL * 2.0;
    double pi = 3.14159265;
    // Convert to flux?
    sigmaL = sigmaL / (4*pi);
    // Convert to parsecs^2 using unit_l
    sigmaL = sigmaL / (unit_l * unit_l);
    // How many parsecs are in an arcsecond^2?
    double sightarea = (pi*10./648000.)*(pi*10./648000.);
    // Convert to flux per arcsecond^2
    sigmaL = sigmaL * sightarea;
    // Find magnitude
    double sunmag = 4.83;
    double magnitude = -2.5*log10(sigmaL) + sunmag;
    if (sigma > 0.0){
    fprintf(write, "%24.16e %24.16e %24.16e\n", rr, sigma, magnitude);
    } else
      {
      }
  }
   return 0;
}
 
void Init_D_compute(){
  compute_class = rb_define_class("D_compute", rb_cObject);
  rb_define_method(compute_class, "do", c_do, 2);
  }
