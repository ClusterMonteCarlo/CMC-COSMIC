/* -*- linux-c -*- */
/* cluster.h

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define FB_TIDALTOL 1.0e-3

/* number of stars in cluster */
#define FB_N 20

/* mass of each star in cluster */
#define FB_M (1.0 * FB_CONST_MSUN)

/* radius of each star in cluster */
#define FB_R (1.0 * FB_CONST_RSUN)

/* central velocity dispersion */
#define FB_SIGMA (10.0e5)

/* radius at which to truncate Plummer model */
#define FB_RMAX (10.0 * FB_CONST_PARSEC)

#define FB_DT 1.0 /* approximate output dt */
#define FB_TSTOP 1.0e7 /* in units of t_dyn */
#define FB_TPHYSSTOP (1.0e10 * FB_CONST_YR) /* in cgs */
#define FB_TCPUSTOP 28800.0 /* in seconds */

#define FB_ABSACC 1.0e-9 /* absolute accuracy of integrator */
#define FB_RELACC 1.0e-9 /* relative accuracy of integrator */
#define FB_NCOUNT 500 /* number of timesteps between calls to classify() */

#define FB_KS 0

#define FB_FEXP 3.0 /* expansion factor of merger product */

#define FB_SEED 0UL
#define FB_DEBUG 0

void print_usage(FILE *stream);
void calc_units(fb_hier_t hier, fb_units_t *units);
double fv(double v, void *params);
double vf(double f);
