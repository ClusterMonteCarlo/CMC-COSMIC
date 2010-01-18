#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "bse_wrap.h"
#include "popsynth.h"

void printoutputheader(void);
void printoutput(double time, binary_t *binarray);

int main(int argc, char *argv[]) {
  unsigned long int seed=235UL;
  binary_t *binarray;
  binary_t *lilbinarray;
  binary_t dummybinaryt;
  int i, j, idctr, k, active, tag, krecd;
  double X, norm, mmin, mmax, index, mass;
  gsl_rng *rng;
  const gsl_rng_type *rng_type=gsl_rng_mt19937;
  double z=0.03, *zpars, vs[3];
  double amin, amax, targettphysf, tphysf, tphysfbse, dtp, rlprimovera;
  int numprocs, myid;
  MPI_Status stat; 
  /* MPI derived datatype */
  int blockcounts[4];
  MPI_Datatype oldtypes[4];
  MPI_Aint offsets[4];
  MPI_Datatype mpibinaryt;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  /* set up MPI derived datatype */
  MPI_Address(&dummybinaryt.id1, &offsets[0]);
  MPI_Address(&dummybinaryt.a, &offsets[1]);
  MPI_Address(&dummybinaryt.bse_kw, &offsets[2]);
  MPI_Address(&dummybinaryt.bse_mass0, &offsets[3]);
  
  for (i=3; i>=0; i--) {
    offsets[i] -= offsets[0];
  }

  oldtypes[0] = MPI_LONG;
  blockcounts[0] = 2;
  oldtypes[1] = MPI_DOUBLE;
  blockcounts[1] = 2;
  oldtypes[2] = MPI_INT;
  blockcounts[2] = 2;
  oldtypes[3] = MPI_DOUBLE;
  blockcounts[3] = 32;

  MPI_Type_struct(4, blockcounts, offsets, oldtypes, &mpibinaryt);
  MPI_Type_commit(&mpibinaryt);
  
  lilbinarray = (binary_t *) malloc(BATCHSIZE * sizeof(binary_t));

  /* initialize stellar evolution */
  zpars = (double *) malloc(20 * sizeof(double));
  
  /* evolution parameters */
  bse_set_neta(0.5);
  bse_set_bwind(0.0);
  bse_set_hewind(1.0);
  bse_set_alpha1(0.5);
  bse_set_lambda(1.0);
  bse_set_ceflag(0);
  bse_set_tflag(1);
  bse_set_ifflag(0);
  bse_set_wdflag(1);
  bse_set_bhflag(0);
  bse_set_nsflag(1);
  bse_set_mxns(3.0);
  /* must fix this for MPI, since each process should have a different seed */
  bse_set_idum(29769);
  bse_set_pts1(0.05);
  bse_set_pts2(0.01);
  bse_set_pts3(0.02);
  bse_set_sigma(265.0);
  bse_set_beta(0.125);
  bse_set_xi(1.0);
  bse_set_acc2(1.5);
  bse_set_epsnov(0.001);
  bse_set_eddfac(2.0);
  bse_set_gamma(-1.0);
  
  bse_zcnsts(&z, zpars);
  bse_instar();

  if (myid == 0) {
    /* initialize GSL rng */
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(rng, seed);
    
    binarray = (binary_t *) malloc(NBIN * sizeof(binary_t));

    idctr = 0;
    for (i=0; i<NBIN; i++) {
      /* set arbitrary but sensible IDs */
      binarray[i].id1 = idctr++;
      binarray[i].id2 = idctr++;
      
      /* zero out some vars */
      binarray[i].bse_ospin[0] = 0.0;
      binarray[i].bse_ospin[1] = 0.0;
      binarray[i].bse_epoch[0] = 0.0;
      binarray[i].bse_epoch[1] = 0.0;
      binarray[i].bse_tphys = 0.0;
      
      /* set primary mass from power law */
      X = gsl_rng_uniform(rng);
      /* comparison with Tassos */
      /* mmin = 6.0; */
      /* mmax = 150.0; */
      /* comparison with Ivanova, et al. (2006) to form CVs */
      mmin = 0.5; /* mmin of secondary is 0.2 */
      mmax = 10.0;
      index = -2.35;
      norm = pow(mmax/mmin, index+1.0) - 1.0;
      mass = mmin*pow(norm*X+1, 1.0/(index+1.0));
      
      /* set secondary mass from uniform distribution in mass ratio */
      binarray[i].bse_mass0[0] = mass;
      binarray[i].bse_mass0[1] = gsl_rng_uniform(rng) * (mass - 0.08) + 0.08;
      
      binarray[i].bse_mass[0] = binarray[i].bse_mass0[0];
      binarray[i].bse_mass[1] = binarray[i].bse_mass0[1];
      
      binarray[i].bse_kw[0] = 1;
      binarray[i].bse_kw[1] = 1;
      
      /* evolve each star slightly to set radius */
      dtp = 0.0;
      binarray[i].bse_tphys = 0.0;
      tphysf = 1.0;
      bse_evolv1(&(binarray[i].bse_kw[0]), &(binarray[i].bse_mass0[0]), 
		 &(binarray[i].bse_mass[0]), &(binarray[i].bse_radius[0]), &(binarray[i].bse_lum[0]), 
		 &(binarray[i].bse_massc[0]), &(binarray[i].bse_radc[0]), 
		 &(binarray[i].bse_menv[0]), &(binarray[i].bse_renv[0]), 
		 &(binarray[i].bse_ospin[0]), &(binarray[i].bse_epoch[0]), 
		 &(binarray[i].bse_tms[0]), 
		 &binarray[i].bse_tphys, 
		 &tphysf, &dtp, &z, zpars, vs);
      
      dtp = 0.0;
      binarray[i].bse_tphys = 0.0;
      tphysf = 1.0;
      bse_evolv1(&(binarray[i].bse_kw[1]), &(binarray[i].bse_mass0[1]), 
		 &(binarray[i].bse_mass[1]), &(binarray[i].bse_radius[1]), &(binarray[i].bse_lum[1]), 
		 &(binarray[i].bse_massc[1]), &(binarray[i].bse_radc[1]), 
		 &(binarray[i].bse_menv[1]), &(binarray[i].bse_renv[1]), 
		 &(binarray[i].bse_ospin[1]), &(binarray[i].bse_epoch[1]), 
		 &(binarray[i].bse_tms[1]), 
		 &binarray[i].bse_tphys, 
		 &tphysf, &dtp, &z, zpars, vs);
      
      /* set semimajor axis uniform in log from primary radius at half Roche lobe to 10^5 R_sun, 
	 eccentricity thermal */
      rlprimovera = 0.46224 * pow(binarray[i].bse_mass[0]/(binarray[i].bse_mass[0]+binarray[i].bse_mass[1]), 1.0/3.0);
      /* comparison with Tassos */
      /* amin = 2.0 * binarray[i].bse_radius[0] / rlprimovera * RSUN / AU; */
      /* amax = 1.0e5 * RSUN / AU; */
      /* set semimajor axis uniform in log from 1d to 10^4 d orbital period */
      amin = pow((binarray[i].bse_mass[0]+binarray[i].bse_mass[1])*(1.0/365.25)*(1.0/365.25), 1.0/3.0);
      amax = pow((binarray[i].bse_mass[0]+binarray[i].bse_mass[1])*(1.0e4/365.25)*(1.0e4/365.25), 1.0/3.0);
      binarray[i].a = pow(10.0, gsl_rng_uniform(rng)*(log10(amax)-log10(amin))+log10(amin));
      binarray[i].bse_tb = 365.25 * sqrt(pow(binarray[i].a, 3.0)/(binarray[i].bse_mass0[0]+binarray[i].bse_mass0[1]));
      /* binarray[i].e = sqrt(gsl_rng_uniform(rng)); */
      binarray[i].e = 0.0;
      
      /* save initial state */
      binarray[i].m1init = binarray[i].bse_mass[0];
      binarray[i].m2init = binarray[i].bse_mass[1];
      binarray[i].ainit = binarray[i].a;
      binarray[i].einit = binarray[i].e;
    }
  }

  if (myid == 0) {
    printoutputheader();
  }
  
  /* distribute and collect data */
  tphysf = 5.0;
  /* targettphysf = 10.0; */
  targettphysf = 13.0e3;
  while (tphysf <= targettphysf) {
    if (myid == 0) {
      fprintf(stderr, "time=%g\n", tphysf);
      k = 0;
      active = 0;
      for (i=1; i<numprocs; i++) {
	MPI_Send(&k, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
	MPI_Send(&(binarray[k]), BATCHSIZE, mpibinaryt, i, 1, MPI_COMM_WORLD);
	/* fprintf(stderr, "master sent initial data: time=%g k=%d\n", tphysf, k); */
	k += BATCHSIZE;
	active++;
      }
      while (active) {
	MPI_Recv(&krecd, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
	MPI_Recv(&(binarray[krecd]), BATCHSIZE, mpibinaryt, stat.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
	if (k < NBIN) {
	  MPI_Send(&k, 1, MPI_INT, stat.MPI_SOURCE, 1, MPI_COMM_WORLD);
	  MPI_Send(&(binarray[k]), BATCHSIZE, mpibinaryt, stat.MPI_SOURCE, 1, MPI_COMM_WORLD);
	  /* fprintf(stderr, "master sent more data: time=%g k=%d\n", tphysf, k); */
	  k += BATCHSIZE;
	} else {
	  MPI_Send(&k, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
	  MPI_Send(&(binarray[0]), BATCHSIZE, mpibinaryt, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
	  active--;
	}
	printoutput(tphysf, &(binarray[krecd]));
      }
    } else {
      tag = 1;
      while (tag) {
	MPI_Recv(&krecd, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
	MPI_Recv(lilbinarray, BATCHSIZE, mpibinaryt, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
	tag = stat.MPI_TAG;

	/* do evolution if necessary */
	if (tag) {
	  /* fprintf(stderr, "myid=%d: Got some data, doing some work. time=%g\n", myid, tphysf); */
	  for (i=0; i<BATCHSIZE; i++) {
	    /* fprintf(stderr, "myid=%d, i=%d\n", myid, i); */
	    tphysfbse = tphysf;
	    dtp = 0.0;
	    bse_evolv2(&(lilbinarray[i].bse_kw[0]), &(lilbinarray[i].bse_mass0[0]), &(lilbinarray[i].bse_mass[0]), 
		       &(lilbinarray[i].bse_radius[0]), &(lilbinarray[i].bse_lum[0]), 
		       &(lilbinarray[i].bse_massc[0]), &(lilbinarray[i].bse_radc[0]), 
		       &(lilbinarray[i].bse_menv[0]), &(lilbinarray[i].bse_renv[0]), 
		       &(lilbinarray[i].bse_ospin[0]), &(lilbinarray[i].bse_epoch[0]), 
		       &(lilbinarray[i].bse_tms[0]), 
		       &lilbinarray[i].bse_tphys, &tphysfbse, &dtp, 
		       &z, zpars, 
		       &lilbinarray[i].bse_tb, &lilbinarray[i].e, vs);
	    
	    /* extract some binary info from BSE's bcm array */
	    j = 1;
	    while (bse_get_bcm(j, 1) >= 0.0) {
	      j++;
	    }
	    j--;
	    if (j >= 1) {
	      /* if (fabs((lilbinarray[i].bse_tphys - bse_get_bcm(j,1))/lilbinarray[i].bse_tphys) >= 1.0e-6) { */
/* 		fprintf(stderr, "lilbinarray[kb].bse_tphys=%g bcmtime=%g\n", lilbinarray[i].bse_tphys, bse_get_bcm(j,1)); */
/* 	      } */
	      lilbinarray[i].bse_bcm_dmdt[0] = bse_get_bcm(j, 14);
	      lilbinarray[i].bse_bcm_dmdt[1] = bse_get_bcm(j, 28);
	      lilbinarray[i].bse_bcm_radrol[0] = bse_get_bcm(j, 15);
	      lilbinarray[i].bse_bcm_radrol[1] = bse_get_bcm(j, 29);
	    } else {
	      fprintf(stderr, "Could not extract BSE bcm info!  Input dtp not exactly equal to tphysf-tphys?\n");
	      exit(1);
	    }
	  }
	  
	  /* send completed data back to master */
	  MPI_Send(&krecd, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	  MPI_Send(lilbinarray, BATCHSIZE, mpibinaryt, 0, 1, MPI_COMM_WORLD);
	} else {
	  /* fprintf(stderr, "myid=%d: Got quit signal, quitting. time=%g\n", myid, tphysf); */
	}
      }
    }
    tphysf += 5.0;
  }
  
  free(lilbinarray);
  
  if (myid == 0) {
    free(binarray);
    gsl_rng_free(rng);
  }

  MPI_Finalize();

  return(0);
}

void printoutputheader(void)
{
  fprintf(stdout, "#0:time #1:id1 #2:id2 #3:m1 #4:m2 #5:k1 #6:k2 #7:tb #8:e #9:rad1 #10:rad2 #11:lum1 #12:lum2 #13:massc1 #14:massc2 #15:radc1 #16:radc2 #17:menv1 #18:menv2 #19:renv1 #20:renv2 #21:ospin1 #22:ospin2 #23:tms1 #24:tms2 #25:dmdt1 #26:dmdt2 #27:radrol1 #28:radrol2 #29:m1init #30:m2init #31:ainit #32:einit\n");
}

void printoutput(double time, binary_t *binarray)
{
  double tb, dmdt0, dmdt1, radrol0, radrol1;
  int i, kw0, kw1;
  /* FILE *ofp; */
  /* char ofpname[1024]; */

/*   sprintf(ofpname, "popsynth.out.dat.gz"); */
/*   if ((ofp = (FILE *) gzopen(ofpname, "wb")) == NULL) { */
/*     fprintf(stderr, "cannot create output file %s\n", ofpname); */
/*     exit(1); */
/*   } */
  
  for (i=0; i<BATCHSIZE; i++) {
    tb = binarray[i].bse_tb;
    kw0 = binarray[i].bse_kw[0];
    kw1 = binarray[i].bse_kw[1];
    dmdt0 = binarray[i].bse_bcm_dmdt[0];
    dmdt1 = binarray[i].bse_bcm_dmdt[1];
    radrol0 = binarray[i].bse_bcm_radrol[0];
    radrol1 = binarray[i].bse_bcm_radrol[1];
    if ( tb > 0.0 && 
	 ((kw0 >= 10 && kw0 <= 12 && kw1 >= 10 && kw1 <= 12 && dmdt0 > 0.0 && radrol1 >= 0.99) || 
	  (kw1 >= 10 && kw1 <= 12 && kw0 >= 10 && kw0 <= 12 && dmdt1 > 0.0 && radrol0 >= 0.99)) ) {
      fprintf(stdout, "%g %ld %ld %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
	      time,
	      binarray[i].id1, binarray[i].id2,
	      binarray[i].bse_mass[0], binarray[i].bse_mass[1],
	      binarray[i].bse_kw[0], binarray[i].bse_kw[1],
	      binarray[i].bse_tb, binarray[i].e,
	      binarray[i].bse_radius[0], binarray[i].bse_radius[1],
	      binarray[i].bse_lum[0], binarray[i].bse_lum[1],
	      binarray[i].bse_massc[0], binarray[i].bse_massc[1],
	      binarray[i].bse_radc[0], binarray[i].bse_radc[1],
	      binarray[i].bse_menv[0], binarray[i].bse_menv[1],
	      binarray[i].bse_renv[0], binarray[i].bse_renv[1],
	      binarray[i].bse_ospin[0], binarray[i].bse_ospin[1],
	      binarray[i].bse_tms[0], binarray[i].bse_tms[1],
	      binarray[i].bse_bcm_dmdt[0], binarray[i].bse_bcm_dmdt[1],
	      binarray[i].bse_bcm_radrol[0], binarray[i].bse_bcm_radrol[1],
	      binarray[i].m1init, binarray[i].m2init, binarray[i].ainit, binarray[i].einit);
    }
  }
  
  /* fclose(ofp); */
}
