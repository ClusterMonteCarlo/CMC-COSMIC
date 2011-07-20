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

/* New rng functions without private state variable */
unsigned long rng_t113_int_new(struct rng_t113_state *state); 
double rng_t113_dbl_new(struct rng_t113_state *state);
void reset_rng_t113_new(unsigned long int s, struct rng_t113_state *state);

/* Jump Polynomials Stuff */
struct state rng_t113_jump( struct rng_t113_state s );
struct rng_113_state rng_t113_next_state( struct rng_t113_state s );
//void rng_states_generate( struct state *h_states, unsigned int seed );

#endif /* _TAUS113_V2_H */
