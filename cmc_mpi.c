#include "cmc_mpi.h"
#include "cmc.h"
#include "cmc_vars.h"

//Function to find start and end indices for each processor for a loop over N. Makes sure the number of elements each proc works on is an even number. Is a generic function that does not use myid, instead it should be passed through parameter i.
void mpiFindIndicesEvenGen( long N, int i, int* mpiStart, int* mpiEnd )
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
		mpiFindIndicesEvenGen( N, myid, mpiStart, mpiEnd );
	} else {
		mpiFindIndicesEvenGen( N-1, myid, mpiStart, mpiEnd );
		if(myid == procs-1)
			*mpiEnd += 1;
	}
}

//Same as mpiFindIndicesSpecial, but does not use myid, instead allows it to be passed through parameter i.
void mpiFindIndicesSpecialGen( long N, int i, int* mpiStart, int* mpiEnd )
{
	if( N % 2 == 0 )
	{
		mpiFindIndicesEvenGen( N, i, mpiStart, mpiEnd );
	} else {
		mpiFindIndicesEvenGen( N-1, i, mpiStart, mpiEnd );
		if(i == procs-1)
			*mpiEnd += 1;
	}
}

/**************************************************************************
Function to find start and end indices for each processor for a loop over N.
If blkSize is 1, N is divided equally amond procs. If blkSize > 0, N*blkSize
stars are divided among procs making sure each contains multiples of blkSize.
**************************************************************************/
void mpiFindIndices( long N, int blkSize, int i, int* mpiStart, int* mpiEnd )
{
	long chunkSize =  ( N / procs ) * blkSize;
   if ( i < N % procs )
   {
      *mpiStart = i * chunkSize + i * blkSize + 1; //+1 since for loops go from 1 to N
		*mpiEnd = *mpiStart + chunkSize + 1 * blkSize - 1;
   } else {
      *mpiStart = i * chunkSize + ( N % procs ) * blkSize + 1; //+1 since for loops go from 1 to N
		*mpiEnd =  *mpiStart + chunkSize - 1;
   }
}

void mpiFindIndicesCustom( long N, int blkSize, int i, int* mpiStart, int* mpiEnd )
{
	int blocks, remain;

	blocks = N / blkSize;
	remain = N % blkSize;

	mpiFindIndices( blocks, blkSize, i, mpiStart, mpiEnd );

	if(i == procs-1)
		*mpiEnd += remain;
}

void mpiFindDispAndLenSpecial( long N, int blkSize, int* mpiDisp, int* mpiLen )
{
	int i;
	for( i = 0; i < procs; i++ )
	{
		mpiFindIndicesCustom( N, blkSize, i, &mpiDisp[i], &mpiLen[i] );
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

//Function to find start and end indices for each processor for a loop over N. The difference is that the number of elements each proc works on is an even number.
/*
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
*/

/* Future Work
void mpiReduceAndBcast( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
{
}
*/
