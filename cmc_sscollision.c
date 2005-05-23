/* -*- linux-c -*- */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

void sscollision_do(long k, long kp, double rcm, double vcm[4])
{
	long knew;

	/* create new star */
	knew = create_star();
	
	/* conserve mass and momentum in collision */
	star[knew].r = rcm;
	star[knew].vr = vcm[3];
	star[knew].vt = sqrt(sqr(vcm[1])+sqr(vcm[2]));
	star[knew].m = star[k].m + star[kp].m;
	star[knew].phi = potential(star[knew].r);
	set_star_EJ(knew);
	/* Calculate internal energy assuming potential doesn't change during collision.
	   Should there be a factor of 1/2 for the potential energy? */
	star[knew].Eint = star[k].Eint + star[kp].Eint 
		+ 0.5 * star[k].m * madhoc * (sqr(star[k].vr) + sqr(star[k].vt)) 
		+ 0.5 * star[kp].m * madhoc * (sqr(star[kp].vr) + sqr(star[kp].vt))
		- 0.5 * star[knew].m * madhoc * (sqr(star[knew].vr) + sqr(star[knew].vt))
		+ 0.5 * star[k].m * madhoc * star[k].phi
		+ 0.5 * star[kp].m * madhoc * star[kp].phi
		- 0.5 * star[knew].m * madhoc * star[knew].phi;
	set_star_news(knew);
	set_star_olds(knew);
	star[knew].interacted = 1;
	star[knew].id = star_get_id_new();
	star[knew].rad = r_of_m(star[knew].m);

	/* log collision */
	fprintf(collisionfile, "t=%g single-single idm=%ld(mm=%g) id1=%ld(m1=%g):id2=%ld(m2=%g) (r=%g)\n", 
		TotalTime, 
		star[knew].id, star[knew].m * units.mstar / FB_CONST_MSUN, 
		star[k].id, star[k].m * units.mstar / FB_CONST_MSUN, 
		star[kp].id, star[kp].m * units.mstar / FB_CONST_MSUN,
		star[knew].r);
	
	/* destroy two progenitors */
	destroy_obj(k);
	destroy_obj(kp);
}
