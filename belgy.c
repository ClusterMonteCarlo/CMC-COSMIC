#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cmc.h"
#include "cmc_vars.h"
#include "sinbin.h"

#ifdef SE

/* star types */
#define S_MASS_MS_STAR                    (0)
#define L_MASS_MS_STAR                    (1)
#define HERTZ_GAP_STAR                    (2)
#define FIRST_GIANT_BRANCH_STAR           (3)
#define CORE_He_BURNING_STAR              (4)
#define EARLY_ASYM_GIANT_BRANCH_STAR      (5)
#define THERM_PULS_ASYM_GIANT_BRANCH_STAR (6)
#define MS_NAKED_He_STAR                  (7)
#define HERTZ_GAP_NAKED_He_STAR           (8)
#define GIANT_BRANCH_NAKED_He_STAR        (9)
#define He_WD                             (10)
#define CO_WD                             (11)
#define ONe_WD                            (12)
#define NEUTRON_STAR                      (13)
#define BLACK_HOLE                        (14)
#define MASSLESS_REMNANT                  (15)
#define NAKED_CO_CORE                     (-1)

/* MEGA_YEAR is 10^6 years in dynamic units.
 * The time written in lagrange radii file, "Totaltime"
 * is in FP units, the exact conversion depends on,
 * among other things, total mass and size of the cluster
 * so it cannot be calculated automatically and is read
 * from the input file.
 */

/* the following is for a Plummer model trh(0)~20Myr */
/* #define MEGA_YEAR 0.004655 */
double dynamic2stellar_time(double time_dyn){
	return time_dyn/MEGA_YEAR;
}
double stellar2dynamic_time(double time_ste){
	return time_ste*MEGA_YEAR;
}

/* SOLAR_MASS_DYN is the mass of Sun in dynamic units. 
 * stellar mass is given in units of solar mass, 
 * dynamic mass is given in strange units, basically it looks like
 * sum m = N, ie. <m> = 1.0, 
 * this transformation loses information so the conversion
 * factor has to be determined at the beginning. it is read from 
 * the input file as SOLAR_MASS_DYN = 1.0/<m> ie. it is what value a 1 solar
 * mass star would assume. In late versions of IMGE this will go also into
 * the FITS file produced by massspec. 
 */

/* the following is a 0.5-100 Salpeter IMF, 65536 stars */
/* #define SOLAR_MASS_DYN 0.6142742970 */
double dynamic2stellar_mass(double mass_dyn){
	return mass_dyn/SOLAR_MASS_DYN;
}
double stellar2dynamic_mass(double mass_ste){
	return mass_ste*SOLAR_MASS_DYN;
}

void stellar_evolution_init(void){
	long int i, k;
	double frac;

	/* the following function calls need to be done for
	 * Chris' stellar evolution code */
	/* -- begin black_magic -- */
	M_hook=M_hookf();    
	M_HeF=M_HeFf();
	M_FGB=M_FGBf();
	coef_aa();
	coef_bb(); 
	/* --  end  black_magic -- */
	
	if((fp0=fopen("error.dat","w"))==NULL) {    /* errors and warnings */
		printf("error: can't open file fp0\n");
		exit_cleanly(-2);
	}

	for(k=1; k<=clus.N_MAX; k++){
		i = k;
		/* converting mass into solar mass unit */
		star[i].mass = dynamic2stellar_mass(star[i].m);
		/* setting the type */
		if(star[i].mass <= 0.7){
			star[i].k = S_MASS_MS_STAR;
		} else {
			star[i].k = L_MASS_MS_STAR;
		}
		/* setting the rest of the variables */
		star[i].mzams = star[i].m0 = star[i].mass;
		star[i].tbeg = star[i].tvir = 0.0;
		star[i].tend = star[i].tbeg + 1e-11;
		star[i].mc = star[i].mcHe = star[i].mcCO = 0.0;
		star[i].dt = star[i].mpre = star[i].tstart = frac = 0.0;
		star[i].kpre = star[i].k;
		star[i].flag = 0;
		star[i].lum = star[i].rad = 1.0;
	}
	
	/* evolving stars a little bit to set luminosity and radius */
	for(k=1; k<=clus.N_MAX; k++){
		i = k;
		frac = 0.0;
		singl(&(star[i].mzams), &(star[i].m0), &(star[i].mass), 
			&(star[i].k), &(star[i].tbeg), &(star[i].tvir),
			&(star[i].tend), &(star[i].lum), &(star[i].rad),
			&(star[i].mc), &(star[i].mcHe), &(star[i].mcCO),
			&(star[i].flag), &(star[i].dt), &(star[i].mpre),
			&(star[i].kpre), &(star[i].tstart), &frac);
	}
}

/* An accelaration scheme for SE:
 * the SE routine is called after dynamical evolution
 * tdesired is the time reached by dynamical evolution
 * star[i].tend is the time up to which star have already evolved
 *        .dt is the time scale over which SE stuff will change ~1%
 * therefore, if (tdesired-star[i].tend) < 0.5*star[i].dt, no
 * need to evolve, one can skip the star.
 */

void do_stellar_evolution(void){
	long int i, k;
	double frac, tdesired;
	
	tdesired = dynamic2stellar_time(TotalTime);
//	if(star[1].tend > tdesired){
//		printf(" PROBABLY NO  STELLAR EVOLUTION\n");
//	} else {
//		printf(" YES STELLAR EVOLUTION\n");
//	}
	for(k=1; k<=clus.N_MAX; k++){
		i = k;
		if ((tdesired-star[i].tend) < 0.5*star[i].dt) continue;
//		printf("Evolving star no:%8ld, mass = %e \r", 
//				i, star[i].mass);
		while (star[i].tend<tdesired) {
			star[i].tbeg = star[i].tend;
			star[i].tend = tdesired;
			frac = 0.0;
			singl(&(star[i].mzams), &(star[i].m0), &(star[i].mass), 
				&(star[i].k), &(star[i].tbeg), &(star[i].tvir),
				&(star[i].tend), &(star[i].lum), &(star[i].rad),
				&(star[i].mc), &(star[i].mcHe), &(star[i].mcCO),
				&(star[i].flag), &(star[i].dt), &(star[i].mpre),
				&(star[i].kpre), &(star[i].tstart), &frac);
		}
		star[i].m = stellar2dynamic_mass(star[i].mass);
	}
	printf("\n");
}

void write_stellar_data(void){
	long i, k;
	FILE *stel_file;
	char filename[1024];

	se_file_counter++;
	sprintf(filename, "%s_stellar_info.%05d", outprefix, se_file_counter);
	stel_file = fopen(filename, "w");
	if (stel_file==NULL){
		fprintf(stderr,
			"file cannot be opened to write stellar info\n");
		return;
	}
	fprintf(stel_file, "# time (Myr): %e\n", 
			dynamic2stellar_time(TotalTime));
	fprintf(stel_file, "# time (FP):  %e\n", TotalTime);
	fprintf(stel_file,
	       "#  id        mass        radius     luminosity  type\n");
	fprintf(stel_file,
	       "#======= ============ ============ ============ ====\n");
	for(k=1; k<=clus.N_MAX; k++){
		i = k;
		fprintf(stel_file, "%08ld ", k);
		fprintf(stel_file, "%e ", star[k].mass);
		fprintf(stel_file, "%e ", star[k].rad);
		fprintf(stel_file, "%e ", star[k].lum);
		fprintf(stel_file, "%2d ", star[k].k);
		fprintf(stel_file, "\n");
	}
	fclose(stel_file);
}

#endif
