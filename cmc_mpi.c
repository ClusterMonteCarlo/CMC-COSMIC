#include "cmc_mpi.h"

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

