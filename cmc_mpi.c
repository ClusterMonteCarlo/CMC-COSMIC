#include "cmc_mpi.h"
#include "cmc.h"
#include "cmc_vars.h"


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

//Same as mpiFindIndicesEven, but has input i instead of myid.
void mpiFindIndicesEven2( long N, int i, int* mpiStart, int* mpiEnd )
{
	N = N/2;
	int temp;
	long chunkSize = N / procs;
   if ( i < N % procs )
   {
		temp = i * chunkSize + i;
		*mpiStart = temp * 2 + 1; //+1 since for loops go from 1 to N
		*mpiEnd = ( temp + chunkSize + 1 ) * 2 + 1 - 1;
   } else {
		temp = i * chunkSize + N % procs;
		*mpiStart = temp * 2 + 1; //+1 since for loops go from 1 to N
		*mpiEnd =  ( temp + chunkSize ) * 2 + 1 - 1;
   }
}

//Function specially made for our code as dynamics_apply() takes 2 stars at a time. This function divides particles in pairs if the no.of stars is even, but if the no.of stars is odd, it gives one star to the last processor (since only this will be skipped by the serial code, and hence will help comparison) and divides the rest in pairs.
void mpiFindIndicesSpecial( long N, int* mpiStart, int* mpiEnd )
{
	if( N % 2 == 0 )
	{
		mpiFindIndicesEven2( N, myid, mpiStart, mpiEnd );
	} else {
		mpiFindIndicesEven2( N-1, myid, mpiStart, mpiEnd );
		if(myid == procs-1)
			*mpiEnd += 1;
	}
}

//Function specially made for our code so that the distribution of stars can be in pairs.
void mpiFindIndicesSpecial2( long N, int i, int* mpiStart, int* mpiEnd )
{
	if( N % 2 == 0 )
	{
		mpiFindIndicesEven2( N, i, mpiStart, mpiEnd );
	} else {
		mpiFindIndicesEven2( N-1, i, mpiStart, mpiEnd );
		if(i == procs-1)
			*mpiEnd += 1;
	}
}

//Function specially made for our code as dynamics_apply() takes 2 stars at a time.
void mpiFindDispAndLenSpecial( long N, int* mpiDisp, int* mpiLen )
{
	int i;
	for( i = 0; i < procs; i++ )
	{
		mpiFindIndicesSpecial2( N, i, &mpiDisp[i], &mpiLen[i] );
		mpiLen[i] = mpiLen[i] - mpiDisp[i] + 1; 
	}
}



//Function to find start index (displacement) and length for each processor for a loop over N. Not used anymore.
/*
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
*/

/* Future Work
void mpiReduceAndBcast( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
{
}
*/
