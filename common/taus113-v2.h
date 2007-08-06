/* taus113-v2.h
 * Copyright (C) 2002 Atakan Gurkan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Declarations useful for using taus113-v2.c */

#ifndef _TAUS113_V2_H
#define _TAUS113_V2_H 1

struct rng_t113_state {
	unsigned long z1, z2, z3, z4;
};

void set_rng_t113(struct rng_t113_state st);
void get_rng_t113(struct rng_t113_state *st);
unsigned long rng_t113_int(void); 
double rng_t113_dbl(void);
void reset_rng_t113(unsigned long int s);

#endif /* _TAUS113_V2_H */
