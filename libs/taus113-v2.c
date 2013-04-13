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
#include <stdlib.h>
#include <string.h>
#include "../common/taus113-v2.h"

#define MASK 0xffffffffUL
#define LCG(n) ((69069UL * n) & 0xffffffffUL)

/**
* @brief Global rng state
*/
static struct rng_t113_state state_static;
/**
* @brief Jump displacement corresponding to 2^20
*/
unsigned int JPoly_2_20[4] = {0x0c382e31, 0x1b040425, 0x0b49a509, 0x0173f6b0};
/**
* @brief Jump displacement corresponding to 2^80
*/
unsigned int JPoly_2_80[4] = {0x487cf69c, 0x00be6310, 0x04bfe2bb, 0x000824f9};

/**
* @brief Sets the global rng state variables to the ones of the input state
*
* @param st Input state
*/
void set_rng_t113(struct rng_t113_state st){
	state_static.z[0] = st.z[0];
	state_static.z[1] = st.z[1];
	state_static.z[2] = st.z[2];
	state_static.z[3] = st.z[3];
}

/**
* @brief Sets the state variabls of the input state to the ones of the global rng
*
* @param st Input state
*/
void get_rng_t113(struct rng_t113_state *st){
	st->z[0] = state_static.z[0];
	st->z[1] = state_static.z[1];
	st->z[2] = state_static.z[2];
	st->z[3] = state_static.z[3];
}

/**
* @brief Returns a random unsigned int value using the global rng state
*
* @return unsigned int random value
*/
unsigned long rng_t113_int() {
	unsigned long b;

	b = ((((state_static.z[0] <<  6) &MASK) ^ state_static.z[0]) >> 13);
	state_static.z[0] = ((((state_static.z[0] & 4294967294UL) << 18) &MASK) ^ b);
	b = ((((state_static.z[1] <<  2) &MASK) ^ state_static.z[1]) >> 27);
	state_static.z[1] = ((((state_static.z[1] & 4294967288UL) <<  2) &MASK) ^ b);
	b = ((((state_static.z[2] << 13) &MASK) ^ state_static.z[2]) >> 21);
	state_static.z[2] = ((((state_static.z[2] & 4294967280UL) <<  7) &MASK) ^ b);
	b = ((((state_static.z[3] <<  3) &MASK) ^ state_static.z[3]) >> 12);
	state_static.z[3] = ((((state_static.z[3] & 4294967168UL) << 13) &MASK) ^ b);
  	return (state_static.z[0] ^ state_static.z[1] ^ state_static.z[2] ^ state_static.z[3]);
}

/**
* @brief Returns a random double value between 0 and 1 using the global rng state
*
* @return double precision random number between 0 and 1
*/
double rng_t113_dbl() {
	return rng_t113_int() / 4294967296.0 ;
}

/**
* @brief Resets global rng using the given seed.
*
* @param s seed
*/
void reset_rng_t113(unsigned long int s) {

	if (s == 0) s = 1UL;	/* default seed is 1 */

	state_static.z[0] = LCG (s);
	if (state_static.z[0] < 2UL) state_static.z[0] += 2UL;
	state_static.z[1] = LCG (state_static.z[0]);
	if (state_static.z[1] < 8UL) state_static.z[1] += 8UL;
	state_static.z[2] = LCG (state_static.z[1]);
	if (state_static.z[2] < 16UL) state_static.z[2] += 16UL;
	state_static.z[3] = LCG (state_static.z[2]);
	if (state_static.z[3] < 128UL) state_static.z[3] += 128UL;

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
/**
* @brief Returns a random unsigned int value
*
* @param state given rng state
*
* @return unsigned int random value
*/
unsigned long rng_t113_int_new(struct rng_t113_state *state) {
	unsigned long b;

	b = ((((state->z[0] <<  6) &MASK) ^ state->z[0]) >> 13);
	state->z[0] = ((((state->z[0] & 4294967294UL) << 18) &MASK) ^ b);
	b = ((((state->z[1] <<  2) &MASK) ^ state->z[1]) >> 27);
	state->z[1] = ((((state->z[1] & 4294967288UL) <<  2) &MASK) ^ b);
	b = ((((state->z[2] << 13) &MASK) ^ state->z[2]) >> 21);
	state->z[2] = ((((state->z[2] & 4294967280UL) <<  7) &MASK) ^ b);
	b = ((((state->z[3] <<  3) &MASK) ^ state->z[3]) >> 12);
	state->z[3] = ((((state->z[3] & 4294967168UL) << 13) &MASK) ^ b);
  	return (state->z[0] ^ state->z[1] ^ state->z[2] ^ state->z[3]);
}

/**
* @brief Returns a random double value between 0 and 1
*
* @param state given rng state
*
* @return double precision random number between 0 and 1
*/
double rng_t113_dbl_new(struct rng_t113_state *state) {
	return rng_t113_int_new(state) / 4294967296.0 ;
}

/**
* @brief Resets given rng using the given seed.
*
* @param s seed
* @param state given rng state
*/
void reset_rng_t113_new(unsigned long int s, struct rng_t113_state *state) {

	if (s == 0) s = 1UL;	/* default seed is 1 */

	state->z[0] = LCG (s);
	if (state->z[0] < 2UL) state->z[0] += 2UL;
	state->z[1] = LCG (state->z[0]);
	if (state->z[1] < 8UL) state->z[1] += 8UL;
	state->z[2] = LCG (state->z[1]);
	if (state->z[2] < 16UL) state->z[2] += 16UL;
	state->z[3] = LCG (state->z[2]);
	if (state->z[3] < 128UL) state->z[3] += 128UL;

	/* Calling RNG ten times to satify recurrence condition */
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); rng_t113_int_new(state); rng_t113_int_new(state); 
	rng_t113_int_new(state); 
	return;
}

//========================================================================
// Functions for jump polynomials
//========================================================================
/**
* @brief Returns the next state of the given rng state
*
* @param state given rng state
*
* @return the next state
*/
struct rng_t113_state
rng_t113_next_state(struct rng_t113_state state) {
  unsigned long b;

  b = ((((state.z[0] <<  6) &MASK) ^ state.z[0]) >> 13);
  state.z[0] = ((((state.z[0] & 4294967294UL) << 18) &MASK) ^ b);
  b = ((((state.z[1] <<  2) &MASK) ^ state.z[1]) >> 27);
  state.z[1] = ((((state.z[1] & 4294967288UL) <<  2) &MASK) ^ b);
  b = ((((state.z[2] << 13) &MASK) ^ state.z[2]) >> 21);
  state.z[2] = ((((state.z[2] & 4294967280UL) <<  7) &MASK) ^ b);
  b = ((((state.z[3] <<  3) &MASK) ^ state.z[3]) >> 12);
  state.z[3] = ((((state.z[3] & 4294967168UL) << 13) &MASK) ^ b);

  return (state);
}

/**
* @brief Implementation from "Testing, Selection, and Implementation of Random Number Generators" by Joseph C. Collins, Army Research Laboratory.
Given an state, jumps an amount of states corresponding to the jump displacement jpoly, and returns the new state.
*
* @param st given rng state
* @param jpoly jump displacement
*
* @return rng state after jumping
*/
struct rng_t113_state
rng_t113_jump(struct rng_t113_state st, unsigned int *jpoly) {
  unsigned int testbit;
  int i, j;
  struct rng_t113_state jumpstate= {{0, 0, 0, 0}};

  testbit=1;
  for (i=0; i< 32; i++) {
    for (j=0; j<4; j++) {
      if (jpoly[j]&testbit) {
        jumpstate.z[j]^=st.z[j];
      }
    }
    testbit= testbit<<1;
    st= rng_t113_next_state(st);
  }
  return (jumpstate);
}

void testjumping(void) {
  struct rng_t113_state s, jumps;
  unsigned int i;

  reset_rng_t113_new( 123456, &s );
  for (i=0; i<4; i++) {
    printf("z[%i]= %li\n", i,s.z[i]);
  }
  jumps= s;
  for (i=0; i<4; i++) {
    printf("jz[%i]= %li\n", i,jumps.z[i]);
  }
  for (i=0; i< 1048576; i++) {
    jumps= rng_t113_next_state(jumps);
  }
  printf("The random number after 2^20 iterations: %li\n", rng_t113_int_new(&jumps));

  jumps= rng_t113_jump(s, JPoly_2_20);
  printf("The random number after a jump: %li\n", rng_t113_int_new(&jumps));
  printf("Compare for yourself ... :)\n");
}
