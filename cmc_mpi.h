/* Contains all MPI related variables and functions */

/*
#ifndef USE_MPI
#define USE_MPI
#endif
*/
#include <mpi.h>

/* Global variables required for any code using MPI */
int myid, procs;

void mpiFindIndices( long N, int* mpiOff, int* mpiLen );
void mpiFindIndicesEven( long N, int* mpiStart, int* mpiEnd );
void mpiFindDispAndLen( long N, int* mpiDisp, int* mpiLen );
