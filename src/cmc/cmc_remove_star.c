/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cmc.h"
#include "cmc_vars.h"


/**
* @brief Sourav: checks whether a star should die or not before performing dynamics
*
* @param time ?
* @param k index of star
*
* @return ?
*/
int remove_old_star(double time, long k)
{
	long removed_id=0, knew, removed_bin_id, removed_component_id, left_component_id;	
	double mass=0.0, birth_time=0.0, life_time=0.0,bin_mass1, bin_mass2, 
		bin_mass_removed, bin_mass_left, removed_birth_time, removed_life_time, left_single_mass ;

	if (star[k].binind){
		if ((time-binary[star[k].binind].createtime_m1)>binary[star[k].binind].lifetime_m1){
			fprintf (removestarfile, "In the first loop:star[k].binind= %ld star[k].id= %ld\n", star[k].binind, star[k].id);
			/* create two stars for the binary components */
			knew = create_star(k, 0);
			/*knewp = create_star();*/
			
			/* and set the stars' properties */
			star_r[get_global_idx(knew)] = star_r[get_global_idx(k)];
			star[knew].vr = star[k].vr;
			/*star[knewp].vr = star[k].vr;*/
			
			star[knew].vt = star[k].vt;
			/*star[knewp].vt = star[k].vt;*/
			
			star_m[get_global_idx(knew)] = binary[star[k].binind].m2;
			star_phi[get_global_idx(knew)] = star_phi[get_global_idx(k)];

			set_star_EJ(knew);
			/*set_star_EJ(knewp);*/
			
			set_star_news(knew);
			/*set_star_news(knewp);*/
			
			set_star_olds(knew);
			/*set_star_olds(knewp);*/
			
			/* mark stars as interacted so they don't undergo E_CONS mode stuff */
			star[knew].interacted = 1;
			/*star[knewp].interacted = 1;*/
			
			star[knew].Eint = binary[star[k].binind].Eint2;
			/*star[knewp].Eint = binary[star[k].binind].Eint2;*/
			
			//debug stuff
			//star[knew].id = binary[star[k].binind].id2;
			//star[knew].id = star_get_id_new();
			star[knew].id = star_get_merger_id_new(binary[star[k].binind].id1, binary[star[k].binind].id2);
			star[knew].binind = 0;
	
			/*star[knewp].id = binary[star[k].binind].id2;*/
			
			star[knew].rad = binary[star[k].binind].rad2;
			/*star[knewp].rad = binary[star[k].binind].rad2;*/

			/*some variables we would like to log for bookkeeping */
			removed_id = star[k].id;	
			removed_bin_id = star[k].binind;
			removed_component_id = binary[star[k].binind].id1;
			left_component_id = binary[star[k].binind].id2;
			bin_mass1 = binary[star[k].binind].m1;
			bin_mass2 = binary[star[k].binind].m2;
			DMrejuv+= binary[star[k].binind].m1 * madhoc;  //Storing DMrejuv in Nbody units since comparing with mtotal later 
			bin_mass_removed = binary[star[k].binind].m1;
			bin_mass_left = binary[star[k].binind].m2;
			removed_birth_time = binary[star[k].binind].createtime_m1 ;
			removed_life_time = binary[star[k].binind].lifetime_m1;
			left_single_mass = star_m[get_global_idx(knew)];
			
			/* destroy this binary */
			destroy_obj(k);
			/*destroy_obj(knew);*/

			/*Now print all that stuff in the removed star log file*/
			fprintf (removestarfile, "binary destroyed: %g %ld %ld %ld %ld %g %g %g %g %g %g %g %g\n",
					time,
					removed_id,
					removed_bin_id,
					removed_component_id,
					left_component_id,
					bin_mass1*units.mstar/MSUN,
					bin_mass2*units.mstar/MSUN,
					bin_mass_removed*units.mstar/MSUN,
					bin_mass_left*units.mstar/MSUN,
					left_single_mass*units.mstar/MSUN,
					(time-removed_birth_time)*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					removed_birth_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					removed_life_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9);
			
		}
		else if ((time-binary[star[k].binind].createtime_m2)>binary[star[k].binind].lifetime_m2){
			fprintf (removestarfile, "In the second loop:star[k].binind= %ld star[k].id= %ld\n", star[k].binind, star[k].id);
			/* create two stars for the binary components */
			knew = create_star(k, 0);
			
			star_r[get_global_idx(knew)] = star_r[get_global_idx(k)];
			star[knew].vr = star[k].vr;
			/*star[knewp].vr = star[k].vr;*/
			
			star[knew].vt = star[k].vt;
			/*star[knewp].vt = star[k].vt;*/
			
			star_m[get_global_idx(knew)] = binary[star[k].binind].m1;
			star_phi[get_global_idx(knew)] = star_phi[get_global_idx(k)];
			
			set_star_EJ(knew);
			/*set_star_EJ(knewp);*/
			
			set_star_news(knew);
			/*set_star_news(knewp);*/
			
			set_star_olds(knew);
			/*set_star_olds(knewp);*/
			
			/* mark stars as interacted so they don't undergo E_CONS mode stuff */
			star[knew].interacted = 1;
			/*star[knewp].interacted = 1;*/
			
			star[knew].Eint = binary[star[k].binind].Eint1;
			/*star[knewp].Eint = binary[star[k].binind].Eint2;*/
			
			//debug Stuff
			//star[knew].id = binary[star[k].binind].id1;
			//star[knew].id = star_get_id_new();
			star[knew].id = star_get_merger_id_new(binary[star[k].binind].id1, binary[star[k].binind].id2);
			star[knew].binind = 0;

			/*star[knewp].id = binary[star[k].binind].id2;*/
			
			star[knew].rad = binary[star[k].binind].rad1;
			/*star[knewp].rad = binary[star[k].binind].rad2;*/
			
			/*some variables we would like to log for bookkeeping */
			removed_id = star[k].id;	
			removed_bin_id = star[k].binind;
			removed_component_id = binary[star[k].binind].id2;
			left_component_id = binary[star[k].binind].id1;
			bin_mass1 = binary[star[k].binind].m1;
			bin_mass2 = binary[star[k].binind].m2;
			DMrejuv+= binary[star[k].binind].m2 * madhoc; //Storing DMrejuv in Nbody units since comparing with mtotal later
			bin_mass_removed = binary[star[k].binind].m2;
			bin_mass_left = binary[star[k].binind].m1;
			removed_birth_time = binary[star[k].binind].createtime_m2 ;
			removed_life_time = binary[star[k].binind].lifetime_m2;
			left_single_mass = star_m[get_global_idx(knew)];
			
			/* destroy this binary */

			//debug stuff
			//destroy_obj(star[k].binind);
			destroy_obj(k);
			
			/*Now print all that stuff in the removed star log file*/
			fprintf (removestarfile, "binary destroyed: %g %ld %ld %ld %ld %g %g %g %g %g %g %g %g\n",
					time,
					removed_id,
					removed_bin_id,
					removed_component_id,
					left_component_id,
					bin_mass1*units.mstar/MSUN,
					bin_mass2*units.mstar/MSUN,
					bin_mass_removed*units.mstar/MSUN,
					bin_mass_left*units.mstar/MSUN,
					left_single_mass*units.mstar/MSUN,
					(time-removed_birth_time)*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					removed_birth_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					removed_life_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9);
			
		
		}
	
	}

	else { 
		if ((time-star[k].createtime)>star[k].lifetime){
			removed_id = star[k].id;
	 		mass = star_m[get_global_idx(k)];
			DMrejuv+= star_m[get_global_idx(k)] * madhoc; //Storing DMrejuv in Nbody units since comparing with mtotal later
			birth_time = star[k].createtime;
			life_time = star[k].lifetime;
			destroy_obj(k);
			fprintf (removestarfile, "single_destroyed: %g %ld %g %g %g %g\n",
					time,
					removed_id,
					mass*units.mstar/MSUN,
					(time-birth_time)*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					birth_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9,
					life_time*units.t*clus.N_STAR/log(GAMMA*clus.N_STAR)/YEAR/1.0e9);
		}
	}
	
	return(0);
}
