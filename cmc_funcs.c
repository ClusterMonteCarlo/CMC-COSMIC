/* -*- linux-c -*- */
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

/* binary - single interaction cross section */
double bin_single_sigma(double y)
{
	return(1.0 / (SQR(SQR(1.0 + y)) * sqrt(y)));
}

/* the mass--radius relationship, from Freitag, et al. (2004) */
/* M is expected to be input in the same units as star[].m */
/* radius is output in code (N-body) units */
double r_of_m(double M)
{
	if (M == 0.0) {
		/* special case of massless particles */
		return(0.0);
	} else if (M / clus.N_STAR * units.m / MSUN < 0.1) {
		return(0.1 * RSUN / units.l);
	} else if (M / clus.N_STAR * units.m / MSUN < 1.0) {
		return(RSUN / units.l * (M / clus.N_STAR * units.m / MSUN));
	} else if (M / clus.N_STAR * units.m / MSUN < 120.0) {
		return(RSUN / units.l * pow(M / clus.N_STAR * units.m / MSUN, 0.57));
	} else {
		return(1.6 * RSUN / units.l * pow(M / clus.N_STAR * units.m / MSUN, 0.47));
	}
}
