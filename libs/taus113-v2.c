/* taus113-v2.c
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

/* This is a maximally equidistributed combined, collision free 
   Tausworthe generator, with a period ~2^{113}. The sequence is,

   x_n = (z1_n ^ z2_n ^ z3_n ^ z4_n)  

   b = (((z1_n <<  6) ^ z1_n) >> 13);
   z1_{n+1} = (((z1_n & 4294967294) << 18) ^ b);
   b = (((z2_n <<  2) ^ z2_n) >> 27);
   z2_{n+1} = (((z2_n & 4294967288) <<  2) ^ b);
   b = (((z3_n << 13) ^ z3_n) >> 21);
   z3_{n+1} = (((z3_n & 4294967280) <<  7) ^ b);
   b = (((z4_n <<  3)  ^ z4_n) >> 12);
   z4_{n+1} = (((z4_n & 4294967168) << 13) ^ b);

   computed modulo 2^32. In the formulas above '^' means exclusive-or 
   (C-notation), not exponentiation. 
   The algorithm is for 32-bit integers, hence a bitmask is used to clear 
   all but least significant 32 bits, after left shifts, to make the code 
   work on architectures where integers are 64-bit.

   The generator is initialized with 
   zi = (69069 * z{i+1}) MOD 2^32 where z0 is the seed provided
   During initialization a check is done to make sure that the initial seeds 
   have a required number of their most significant bits set.
   After this, the state is passed through the RNG 10 times to ensure the
   state satisfies a recurrence relation.

   References:
   P. L'Ecuyer, "Tables of Maximally-Equidistributed Combined LFSR Generators",
   Mathematics of Computation, 68, 225 (1999), 261--269.
     http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
   P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe Generators", 
   Mathematics of Computation, 65, 213 (1996), 203--213.
     http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
   the online version of the latter contains corrections to the print version.
*/
#include <stdio.h>

#define MASK 0xffffffffUL
#define LCG(n) ((69069UL * n) & 0xffffffffUL)

struct rng_t113_state {
	unsigned long z1, z2, z3, z4;
};

static struct rng_t113_state state;

void set_rng_t113(struct rng_t113_state st){
	state.z1 = st.z1;
	state.z2 = st.z2;
	state.z3 = st.z3;
	state.z4 = st.z4;
}

void get_rng_t113(struct rng_t113_state *st){
	st->z1 = state.z1;
	st->z2 = state.z2;
	st->z3 = state.z3;
	st->z4 = state.z4;
}

unsigned long rng_t113_int() {
	unsigned long b;

	b = ((((state.z1 <<  6) &MASK) ^ state.z1) >> 13);
	state.z1 = ((((state.z1 & 4294967294UL) << 18) &MASK) ^ b);
	b = ((((state.z2 <<  2) &MASK) ^ state.z2) >> 27);
	state.z2 = ((((state.z2 & 4294967288UL) <<  2) &MASK) ^ b);
	b = ((((state.z3 << 13) &MASK) ^ state.z3) >> 21);
	state.z3 = ((((state.z3 & 4294967280UL) <<  7) &MASK) ^ b);
	b = ((((state.z4 <<  3) &MASK) ^ state.z4) >> 12);
	state.z4 = ((((state.z4 & 4294967168UL) << 13) &MASK) ^ b);
  	return (state.z1 ^ state.z2 ^ state.z3 ^ state.z4);
}

double rng_t113_dbl() {
	return rng_t113_int() / 4294967296.0 ;
}

void reset_rng_t113(unsigned long int s) {

	if (s == 0) s = 1UL;	/* default seed is 1 */

	state.z1 = LCG (s);
	if (state.z1 < 2UL) state.z1 += 2UL;
	state.z2 = LCG (state.z1);
	if (state.z2 < 8UL) state.z2 += 8UL;
	state.z3 = LCG (state.z2);
	if (state.z3 < 16UL) state.z3 += 16UL;
	state.z4 = LCG (state.z3);
	if (state.z4 < 128UL) state.z4 += 128UL;

	/* Calling RNG ten times to satify recurrence condition */
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); 
	return;
}

/*************************************************************/
/*************************************************************/
/* New rng functions without implicit private state variable */
/*************************************************************/
/*************************************************************/
unsigned long rng_t113_int_new(struct rng_t113_state *state) {
	unsigned long b;

	b = ((((state->z1 <<  6) &MASK) ^ state->z1) >> 13);
	state->z1 = ((((state->z1 & 4294967294UL) << 18) &MASK) ^ b);
	b = ((((state->z2 <<  2) &MASK) ^ state->z2) >> 27);
	state->z2 = ((((state->z2 & 4294967288UL) <<  2) &MASK) ^ b);
	b = ((((state->z3 << 13) &MASK) ^ state->z3) >> 21);
	state->z3 = ((((state->z3 & 4294967280UL) <<  7) &MASK) ^ b);
	b = ((((state->z4 <<  3) &MASK) ^ state->z4) >> 12);
	state->z4 = ((((state->z4 & 4294967168UL) << 13) &MASK) ^ b);
  	return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4);
}

double rng_t113_dbl_new(struct rng_t113_state *state) {
	return rng_t113_int_new(state) / 4294967296.0 ;
}

void reset_rng_t113_new(unsigned long int s, struct rng_t113_state *state) {

	if (s == 0) s = 1UL;	/* default seed is 1 */

	state->z1 = LCG (s);
	if (state->z1 < 2UL) state->z1 += 2UL;
	state->z2 = LCG (state->z1);
	if (state->z2 < 8UL) state->z2 += 8UL;
	state->z3 = LCG (state->z2);
	if (state->z3 < 16UL) state->z3 += 16UL;
	state->z4 = LCG (state->z3);
	if (state->z4 < 128UL) state->z4 += 128UL;

	/* Calling RNG ten times to satify recurrence condition */
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); 
	return;
}
