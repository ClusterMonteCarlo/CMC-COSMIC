/* vi: set filetype=c.doxygen: */

/* Contains all MPI related variables and functions */
#include <mpi.h>
#include <stdio.h>

void mpiFindIndices( long N, int blkSize, int i, int* mpiStart, int* mpiEnd );
void mpiFindIndicesSpecial( long N, int* mpiStart, int* mpiEnd );
void mpiFindIndicesSpecialGen( long N, int i, int* mpiStart, int* mpiEnd );
void mpiFindIndicesEven( long N, int* mpiStart, int* mpiEnd );
void mpiFindIndicesEvenGen( long N, int i, int* mpiStart, int* mpiEnd );
void mpiFindDispAndLen( long N, int* mpiDisp, int* mpiLen );
void mpiFindDispAndLenCustom( long N, int blkSize, int* mpiDisp, int* mpiLen );
void mpiFindIndicesCustom( long N, int blkSize, int i, int* mpiStart, int* mpiEnd );
MPI_Comm inv_comm_create(int procs, MPI_Comm old_comm);
