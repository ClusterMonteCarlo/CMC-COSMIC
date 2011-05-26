#include "cmc_mpi.h"
#include "cmc.h"
#include "cmc_vars.h"

void mpiInitBcastGlobArrays()
{
	strcpy(funcName, __FUNCTION__);
	int i;
	star_r = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_m = (double *) malloc(N_STAR_DIM * sizeof(double));
	star_phi = (double *) malloc(N_STAR_DIM * sizeof(double));

	if(myid==0) {
		for(i=0; i<=N_STAR_DIM; i++) {
			star_r[i] = star[i].r;
			star_m[i] = star[i].m;
		}
	}
	MPI_Bcast(star_m, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(star_r, N_STAR_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

//Function to find start and end indices for each processor for a loop over N.
//There are problems if 2nd and 3rd arguments are long ints instead of ints. Be careful.
void mpiFindIndices( long N, int* mpiStart, int* mpiEnd )
{
	long chunkSize =  N / procs;
   if ( myid < N % procs )
   {
      *mpiStart = myid * chunkSize + myid + 1; //+1 since for loops go from 1 to N
		*mpiEnd = *mpiStart + chunkSize + 1 - 1;
      //*mpiLen = chunkSize + 1;
   } else {
      *mpiStart = myid * chunkSize + N % procs + 1; //+1 since for loops go from 1 to N
		*mpiEnd =  *mpiStart + chunkSize - 1;
      //*mpiLen = chunkSize;
   }
}

//Function to find start and end indices for each processor for a loop over N. The difference is that the number of elements each proc works on is an even number.
void mpiFindIndicesEven( long N, int* mpiStart, int* mpiEnd )
{
	N = N/2;
	int temp;
	long chunkSize = N / procs;
   if ( myid < N % procs )
   {
		temp = myid * chunkSize + myid;
		*mpiStart = temp * 2 + 1; //+1 since for loops go from 1 to N
		*mpiEnd = ( temp + chunkSize + 1 ) * 2 + 1 - 1;
   } else {
		temp = myid * chunkSize + N % procs;
		*mpiStart = temp * 2 + 1; //+1 since for loops go from 1 to N
		*mpiEnd =  ( temp + chunkSize ) * 2 + 1 - 1;
   }
}

//Function to find start index (displacement) and length for each processor for a loop over N.
void mpiFindDispAndLen( long N, int* mpiDisp, int* mpiLen )
{
	long chunkSize = N / procs;
	int i;
	for( i = 0; i < procs; i++ )
	{
		if ( i < N % procs )
		{
			*mpiDisp = i * chunkSize + i + 1; //+1 since for loops go from 1 to N
			*mpiLen = chunkSize + 1;
		} else {
			*mpiDisp = i * chunkSize + N % procs + 1; //+1 since for loops go from 1 to N
			*mpiLen =  chunkSize;
		}
		mpiDisp++; //Pointer increments
		mpiLen++;
	}
}

/*
Random number functions modified for parallel implementation. For original, see libs/taus113-v2.c
*/

/*
struct rng_t113_state {
	unsigned long z1, z2, z3, z4;
};

//static struct rng_t113_state state;


unsigned long rng_t113_int(struct rng_t113_state *st) {
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

double rng_t113_dbl(struct rng_t113_state *st) {
	return rng_t113_int(st) / 4294967296.0 ;
}

void reset_rng_t113(unsigned long int s, struct rng_t113_state *st) {

	if (s == 0) s = 1UL;	// default seed is 1 

	state.z1 = LCG (s);
	if (state.z1 < 2UL) state.z1 += 2UL;
	state.z2 = LCG (state.z1);
	if (state.z2 < 8UL) state.z2 += 8UL;
	state.z3 = LCG (state.z2);
	if (state.z3 < 16UL) state.z3 += 16UL;
	state.z4 = LCG (state.z3);
	if (state.z4 < 128UL) state.z4 += 128UL;

	// Calling RNG ten times to satify recurrence condition 
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); rng_t113_int(); rng_t113_int(); 
	rng_t113_int(); 
	return;
}
*/
