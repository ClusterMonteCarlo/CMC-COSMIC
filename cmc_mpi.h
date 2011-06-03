/* Contains all MPI related variables and functions */

/*
#ifndef USE_MPI
#define USE_MPI
#endif
*/
#include <mpi.h>
#include <stdio.h>

/* Global variables required for any code using MPI */
int myid, procs;

void mpiInitBcastGlobArrays();
void mpiFindIndices( long N, int* mpiOff, int* mpiLen );
void mpiFindIndicesSpecial( long N, int* mpiOff, int* mpiLen );
void mpiFindIndicesSpecial2( long N, int i, int* mpiStart, int* mpiEnd );
void mpiFindIndicesEven( long N, int* mpiStart, int* mpiEnd );
void mpiFindIndicesEven2( long N, int i, int* mpiStart, int* mpiEnd );
void mpiFindDispAndLen( long N, int* mpiDisp, int* mpiLen );
void mpiFindDispAndLenSpecial( long N, int* mpiDisp, int* mpiLen );
