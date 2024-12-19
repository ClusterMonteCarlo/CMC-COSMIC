/* vi: set filetype=c.doxygen: */
#include "cmc_mpi.h"
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief Function to find start and end indices for each processor for a loop over N. Makes sure the number of elements each proc works on is an even number. Is a generic function that does not use myid, instead it should be passed through parameter i.
*
* @param N data size for which data partitioning scheme needs to be found
* @param i parameters for processor i that is needed
* @param mpiStart start index of processor i in the global dataset
* @param mpiEnd end index of processor i in the global dataset
*/
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

/**
* @brief Function specially made for our code as dynamics_apply() takes 2 stars at a time. This function divides particles in pairs if the no.of stars is even, but if the no.of stars is odd, it gives one star to the last processor (since only this will be skipped by the serial code, and hence will help comparison) and divides the rest in pairs.
*
* @param N data size for which data partitioning scheme needs to be found
* @param mpiStart array for start indices
* @param mpiEnd array for end indices
*/
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

/**
* @brief Same as mpiFindIndicesSpecial, but does not use myid, instead allows it to be passed through parameter i.
*
* @param N data size for which data partitioning scheme needs to be found
* @param i parameters for processor i that is needed
* @param mpiStart start index of processor i in the global dataset
* @param mpiEnd end index of processor i in the global dataset
*/
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

/**
* @brief
Function to find start and end indices for the given processor i.
If blkSize is 1, N is divided equally among processors. If blkSize > 0, N*blkSize
stars are divided among procs in such a way that each contain multiples of blkSize.
*
* @param N the number of blocks of size blkSize to be divided among processors
* @param blkSize size of blocks
* @param i processor id
* @param mpiStart indices
* @param mpiEnd
*/
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

/**
* @brief Populates the mpiStart and mpiEnd values given the initial data size N, and the factor blkSize which data in each processor is required to be a multiple of (except the last processor). Followin is the data partitioning scheme - Each processor gets data that is a multiple of blkSize, and the remainder after division goes to the last processor. For more details, refer to: http://adsabs.harvard.edu/abs/2013ApJS..204...15P.
*
* @param N data size that needs to be partitioned
* @param blkSize the factor blkSize which data in each processor is required to be a multiple of (except the last one).
* @param mpiStart start index of the data elements in the global array belonging to processor i.
* @param mpiEnd array end index of the data elements in the global array belonging to processor i.
*/
void mpiFindIndicesCustom( long N, int blkSize, int i, int* mpiStart, int* mpiEnd )
{
	int blocks, remain;

	blocks = N / blkSize;
	remain = N % blkSize;

	mpiFindIndices( blocks, blkSize, i, mpiStart, mpiEnd );

	if(i == procs-1)
		*mpiEnd += remain;
}

/**
* @brief finds displacement and count for all processors
*
* @param N data size
* @param blkSize block size which each chunk has to be a multiple of
* @param mpiDisp displacement array
* @param mpiLen count array
*/
void mpiFindDispAndLenCustom( long N, int blkSize, int* mpiDisp, int* mpiLen )
{
	int i;
	for( i = 0; i < procs; i++ )
	{
		mpiFindIndicesCustom( N, blkSize, i, &mpiDisp[i], &mpiLen[i] );
		mpiLen[i] = mpiLen[i] - mpiDisp[i] + 1; 
	}
}


/**
* @brief Creates a new communicator with inverse order of processes.
*
* @param procs number of processors
* @param old_comm old communicator
*
* @return new communicator with inverse order of ranks
*/
MPI_Comm inv_comm_create(int procs, MPI_Comm old_comm)
{
	int i, *newranks;
	MPI_Group orig_group, new_group;
	MPI_Comm  new_comm;

	newranks = (int *) malloc (procs * sizeof(int));
	for (i=0; i<procs; i++)
		newranks[i] = procs - 1 - i;

	MPI_Comm_group(old_comm, &orig_group);
	MPI_Group_incl(orig_group, procs, newranks, &new_group);
	MPI_Comm_create(old_comm, new_group, &new_comm);

	free (newranks);
	return new_comm;
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

void mpiAllocFileBuffers(void){
    int error = 0;
    
    // Initialize all pointers to NULL
    mpi_logfile_buf = mpi_escfile_buf = mpi_binintfile_buf = mpi_collisionfile_buf = 
    mpi_pulsarfile_buf = mpi_morepulsarfile_buf = mpi_newnsfile_buf = mpi_morecollfile_buf = 
    mpi_triplefile_buf = mpi_tidalcapturefile_buf = mpi_tdefile_buf = mpi_semergedisruptfile_buf = 
    mpi_removestarfile_buf = mpi_relaxationfile_buf = mpi_bhsummaryfile_buf = 
    mpi_escbhsummaryfile_buf = mpi_newbhfile_buf = mpi_bhmergerfile_buf = mpi_threebbfile_buf = 
    mpi_threebbprobabilityfile_buf = mpi_lightcollisionfile_buf = mpi_threebbdebugfile_buf = NULL;

    mpi_logfile_wrbuf = mpi_escfile_wrbuf = mpi_binintfile_wrbuf = mpi_collisionfile_wrbuf = 
    mpi_pulsarfile_wrbuf = mpi_morepulsarfile_wrbuf = mpi_newnsfile_wrbuf = mpi_morecollfile_wrbuf = 
    mpi_triplefile_wrbuf = mpi_tidalcapturefile_wrbuf = mpi_tdefile_wrbuf = mpi_semergedisruptfile_wrbuf = 
    mpi_removestarfile_wrbuf = mpi_relaxationfile_wrbuf = mpi_bhsummaryfile_wrbuf = 
    mpi_escbhsummaryfile_wrbuf = mpi_newbhfile_wrbuf = mpi_bhmergerfile_wrbuf = mpi_threebbfile_wrbuf = 
    mpi_threebbprobabilityfile_wrbuf = mpi_lightcollisionfile_wrbuf = mpi_threebbdebugfile_wrbuf = NULL;

    // Regular buffers
    if ((mpi_logfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_escfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_binintfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_collisionfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_pulsarfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_morepulsarfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_newnsfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_morecollfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_triplefile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_tidalcapturefile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_tdefile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_semergedisruptfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_removestarfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_relaxationfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_bhsummaryfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_escbhsummaryfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_newbhfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_bhmergerfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbprobabilityfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_lightcollisionfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbdebugfile_buf = (char *)malloc(STR_BUF_LEN)) == NULL) error = 1;

    // Write buffers
    if ((mpi_logfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_escfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_binintfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_collisionfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_pulsarfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_morepulsarfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_newnsfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_morecollfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_triplefile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_tidalcapturefile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_tdefile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_semergedisruptfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_removestarfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_relaxationfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_bhsummaryfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_escbhsummaryfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_newbhfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_bhmergerfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbprobabilityfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_lightcollisionfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;
    if ((mpi_threebbdebugfile_wrbuf = (char *)malloc(STR_WRBUF_LEN)) == NULL) error = 1;

    if (error) {
        eprintf("Can't allocate buffer for mpi files\n");
        exit_cleanly(-1,__FUNCTION__);
    }
}
/**
 * @brief Frees all MPI buffer memory
 */
void mpiFreeFileBuffers(void) {
    // Regular buffers
    free(mpi_logfile_buf);
    free(mpi_escfile_buf);
    free(mpi_binintfile_buf);
    free(mpi_collisionfile_buf);
    free(mpi_pulsarfile_buf);
    free(mpi_morepulsarfile_buf);
    free(mpi_newnsfile_buf);
    free(mpi_morecollfile_buf);
    free(mpi_triplefile_buf);
    free(mpi_tidalcapturefile_buf);
    free(mpi_tdefile_buf);
    free(mpi_semergedisruptfile_buf);
    free(mpi_removestarfile_buf);
    free(mpi_relaxationfile_buf);
    free(mpi_bhsummaryfile_buf);
    free(mpi_escbhsummaryfile_buf);
    free(mpi_newbhfile_buf);
    free(mpi_bhmergerfile_buf);
    free(mpi_threebbfile_buf);
    free(mpi_threebbprobabilityfile_buf);
    free(mpi_lightcollisionfile_buf);
    free(mpi_threebbdebugfile_buf);

    // Write buffers
    free(mpi_logfile_wrbuf);
    free(mpi_escfile_wrbuf);
    free(mpi_binintfile_wrbuf);
    free(mpi_collisionfile_wrbuf);
    free(mpi_pulsarfile_wrbuf);
    free(mpi_morepulsarfile_wrbuf);
    free(mpi_newnsfile_wrbuf);
    free(mpi_morecollfile_wrbuf);
    free(mpi_triplefile_wrbuf);
    free(mpi_tidalcapturefile_wrbuf);
    free(mpi_tdefile_wrbuf);
    free(mpi_semergedisruptfile_wrbuf);
    free(mpi_removestarfile_wrbuf);
    free(mpi_relaxationfile_wrbuf);
    free(mpi_bhsummaryfile_wrbuf);
    free(mpi_escbhsummaryfile_wrbuf);
    free(mpi_newbhfile_wrbuf);
    free(mpi_bhmergerfile_wrbuf);
    free(mpi_threebbfile_wrbuf);
    free(mpi_threebbprobabilityfile_wrbuf);
    free(mpi_lightcollisionfile_wrbuf);
    free(mpi_threebbdebugfile_wrbuf);
}
