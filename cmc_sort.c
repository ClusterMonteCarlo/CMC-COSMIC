/* -*- linux-c -*- */

#include "cmc.h"
#include "cmc_vars.h"
#include <math.h>
#include <time.h>

#if (defined(USE_THREADS) || defined(USE_THREADS_SORT))

#define GENSORT_NAME                 qsortgs
#define GENSORT_TYPE                 star_t
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            a.r
#define GENSORT_COMPAREKEYS(k1,k2)   k1 < k2

#define GENSORT_USEPOINTERS

#include "common/gensort.h"


#include <math.h>
#include <pthread.h>

struct thread_data_sort {
      int		thread_id;
      star_t	*a;
      int		N;
};

#define NUM_THREADS_SORT	2

static struct thread_data_sort thread_data_sort_array[NUM_THREADS_SORT];

void *pqsorts(void *threadarg) {
	int taskid, N;
	star_t *a;
	struct thread_data_sort *my_data;

	my_data = (struct thread_data_sort *) threadarg;

	taskid = my_data->thread_id;
	N = my_data->N;
	a = my_data->a;
      
	qsortgs(a, N);

	pthread_exit(NULL);

	return(NULL);
}

void find_hoare(star_t *a, int N, int k){
	int i, j, l, r;
	double x;
	star_t w;

	l=0; r=N-1;
	while(l<r){
		x=a[k].r; i=l; j=r;
		while(j>=i){
			while(a[i].r<x) i++;
			while(x<a[j].r) j--;
			if(i<=j){
				w = a[i]; a[i++]=a[j]; a[j--]=w;
			}
		}
		if(j<k) l=i;
		if(k<i) r=j;
	}
}

void modifind_zabrodsky(star_t *a, int N, int k){
	int i, j, l, r;
	double x;
	star_t w;

	l=0; r=N-1;
	while(l<r){
		x=a[k].r; i=l; j=r;
		while (j>=k && k>=i){
			while(a[i].r<x) i++;
			while(x<a[j].r) j--;
			w = a[i]; a[i++]=a[j]; a[j--]=w;
		}
		if (j<k) l=i;
		if (k<i) r=j;
	}
}

#define SFRMAX(a,b) ((a)>(b)) ? (a) : (b)
#define SFRMIN(a,b) ((a)>(b)) ? (b) : (a)

void select_floyd_rivest(star_t *a, int l, int r, int k){
	int n, i, j, s, sd, ll, rr;
	double z, t;
	star_t w;
	
	while(l<r){
		if ( (r-l)>600 ){
			n = r-l+1;
			i = k-l+1;
			z = log(n);
			s = 0.5 * exp(2*z/3);
			sd = copysign(0.5 * sqrt(z*s*(n-s*1.0)/n), i-n/2.0);
			ll = SFRMAX(l, k-i*s*1.0/n+sd);
			rr = SFRMIN(r, k+(n-i)*s*1.0/n+sd);
			select_floyd_rivest(a, ll, rr, k);
		}
		t = a[k].r;
		i = l;
		j = r;
		w = a[l]; a[l]=a[k]; a[k]=w;
		if (a[r].r>t) {
			w = a[r]; a[r]=a[l]; a[l]=w;
		}
		while (i<j) {
			w = a[i]; a[i++]=a[j]; a[j--]=w;
			while(a[i].r<t) i++;
			while(t<a[j].r) j--;
		}
		if (a[l].r == t) {
			w = a[l]; a[l]=a[j]; a[j]=w;
		} else {
			w = a[++j]; a[j]=a[r]; a[r]=w;
		}
		if(j<=k) l = j+1;
		if(k<=j) r = j-1;
	}
}
#undef SFRMAX
#undef SFRMIN

void check_sort(star_t *s, long N){
	long i;
	for(i=1; i<N; i++){
		if(s[i].r>s[i+1].r){
			eprintf("sort not successful\n");
			eprintf("i=%6ld, N=%6ld\n", i, N);
			eprintf("s[i].r=%f, s[i+1].r=%f\n", s[i].r, s[i+1].r);
			exit(12);
		}
	}
}

/**
* @brief quicksort
*
* @param s pointer to array of stars
* @param N data size
*/
void qsorts(star_t *s, long N){
	strcpy(funcName, __FUNCTION__);
	pthread_t threads[NUM_THREADS_SORT];
	pthread_attr_t attr;
	int rc, t;

	pthread_attr_init(&attr);
#if defined(USE_FIND)
	find_hoare(s, N, N/2);
#elif defined(USE_MODIFIND)
	modifind_zabrodsky(s, N, N/2);
#else
	select_floyd_rivest(s, 0, N-1, N/2);
#endif
	/* launch the thread */
	for(t = 1; t<NUM_THREADS_SORT; t++){
		thread_data_sort_array[t].thread_id = t;
		thread_data_sort_array[t].a = s;
		thread_data_sort_array[t].N = N/2;
		rc = pthread_create(&threads[t], NULL, pqsorts, 
					(void *) &thread_data_sort_array[t]);
		if (rc) {
			eprintf("return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	/* do some sorting in main thread too */
	qsortgs(s+N/2+1, N/2);
	/* let the threads join */
	pthread_attr_destroy(&attr);
	for(t=1;t < NUM_THREADS_SORT;t++) {
		rc = pthread_join(threads[t], NULL);
		if (rc) {
			eprintf("return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}
	//check_sort(s, N);
}

#else

#define GENSORT_NAME                 qsorts
#define GENSORT_TYPE                 star_t
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            a.r
#define GENSORT_COMPAREKEYS(k1,k2)   k1 < k2

#define GENSORT_USEPOINTERS

#include "common/gensort.h"
#endif



#ifdef USE_MPI


// ************************************************************** //
// NEWER VERSION OF SAMPLE SORT WITH CLEANER LOAD BALANCING FUNCTION
// ************************************************************** //
//typedef star_t type;

/**
* @brief returns the key of a given data element
*
* @param buf input data element
*
* @return key of the data element
*/
keyType getKey(type* buf) { return (buf->r); }

/**
* @brief comparison function for sort.
*
* @param a first key
* @param b second key
*
* @return returns 1 if a > b, -1 if a < b, 0 if equal.
*/
int compare_keyType (const void * a, const void * b) 
{
	if( *((keyType*)a) < *((keyType*)b) ) return -1;
	if( *((keyType*)a) > *((keyType*)b) ) return 1;
	return 0;
}

/**
* @brief comparison function for sort
*
* @param a first data element
* @param b second data element
*
* @return returns 1 if a > b, -1 if a < b, 0 if equal.
*/
int compare_type (const void * a, const void * b)
{
	if( getKey((type*)a) < getKey((type*)b) ) return -1;
	if( getKey((type*)a) > getKey((type*)b) ) return 1;
	return 0;
}


/**
* @brief given the total number of data elements, number or processors, computes the expected count based on the data partitioning scheme
*
* @param expected_count array where the expected count will be store after computation
* @param N total number of data elements
* @param numprocs number of processors
*/
void find_expected_count( int* expected_count, int N, int numprocs )
{
	findLimits( N, MIN_CHUNK_SIZE );
	int i;
	for(i=0; i<procs; i++)
		expected_count[i] = End[i] - Start[i] + 1;
/*
	for(i=0; i<numprocs; i++)
	{
		expected_count[i] = N / numprocs;
		if (i < N % numprocs) expected_count[i]++;
	}
*/
}

/**
* @brief returns n_samples samples from the given N data elements in array buf
*
* @param buf array that needs to be sampled
* @param sample_array array in which samples will be stored
* @param N number of data elements in buf
* @param n_samples number of samples required
*
* @return 
*/
keyType* sample(type *buf, keyType *sample_array, int N, int n_samples)
{
	//srand ( time(NULL) );
	if(n_samples > N)
		eprintf("Oversampling occurred! id = %d num samples = %d local no. of stars = %d", myid, n_samples, N);

	//MPI: For now, we do regular sampling. In future we might want to explore random sampling.
	int i;
	for(i=1; i<=n_samples; i++)
		// RANDOM
		// if(i==1) srand ( time(NULL) );
		//buf[ rand() % N + 1 ].key;
		// UNIFORM
		sample_array[i-1] = getKey(&buf[ i * N/(n_samples+1) ]);

	return sample_array;
}

/**
* @brief binary search on buf for r
*
* @param buf input array on which search is done
* @param r target value
* @param kmin left start index
* @param kmax right start index
*
* @return index k such that buf[k] < r < buf[k+1]
*/
int binary_search( type* buf, keyType r, int kmin, int kmax )
{
   int ktry;

   if (getKey(&buf[kmin]) > r  || getKey(&buf[kmax]) < r)
      dprintf("r=%f is outside kmin=%f kmax=%f!!\n", r, getKey(&buf[kmin]), getKey(&buf[kmax]));
   if (getKey(&buf[kmin]) > r)
      return kmin;
   if (getKey(&buf[kmax]) < r)
      return kmax+1;

   do {
      ktry = (kmin+kmax+1)/2;
      if (getKey(&buf[ktry])<r)
      {
         kmin = ktry;
      } else {
         kmax = ktry-1;
      }
   } while (kmax!=kmin);

   return kmin+1;
}

/**
* @brief Stripped stars have their r values set to SF_INFINITY. The first step of the parallel sort is a local sort by each processor. During this, these stars are pushed to the end of the local array. This routine removes all of these stars starting from the end of the sorted array.
*
* @param buf array sorted by r
* @param local_N number of stars in the array
*/
void remove_stripped_stars(type* buf, int* local_N)
{
	// Fixing local_N to account for lost stars.
	int i = (*local_N)-1; //Starting from last star
	while(getKey(&buf[i]) == SF_INFINITY)
	{
		dprintf("stripped star found and removed.\n");
		//zero_star(i); //Check with Stefan if this is needed.
		(*local_N)--;
		i--;
	}
}

// NEWER VERSION OF SAMPLE SORT WITH CLEANER LOAD BALANCING FUNCTION
/**
* @brief Parallel sample sort. Following are the steps:
* 1. Local sort - each processor sorts the local data using sequential sort
* 2. Sampling - each processor samples s keys from the local dataset and sends them to the root node. The root node sorts these s*p keys, and picks p regularly sampled keys/points among these, called splitter points which are broadcasted to all processors. These points are supposed to divide the entire dataset into p buckets with approximately equaly number of points in them. This extent to which this division will be equal is a factor of the number of samples s contributed by each processor. If s=p, and regular sampling is used, the load balance will be ideal i.e. the maximum number of points ending up in a processor would be at most ~ 2*N/p. If s < p, it'll be worse, and is s > p, it'll be better.
* 3. Each processor places these splitter points into their sorted local array using binary search.
* 4. All processors have an all-to-all communication and exchange data
* 5. Each processor sorts the received chunks of data which completes the sort
* 6. Since the number of points ending up on each processor is non-deterministic, an optional phase is to exchange data between processors so that the number of data points on each processor is in accordance with out data partitioning scheme.
* @param buf the local data set (star) which is a part of the entire data set which is divided among many processors which is to be sorted in parallel
* @param local_N number of local data points
* @param dataType MPI datatype for the star data structure
* @param b_buf binary local data set
* @param b_dataType  MPI datatype for the binary data structure
* @param commgroup MPI communication group
* @param n_samples number of samples to be contributed by each processor
*
* @return total number of stars that were sorted
*/
int sample_sort( type	      *buf,
                 int	      *local_N,
                 MPI_Datatype dataType,
                 binary_t     *b_buf,
                 MPI_Datatype b_dataType,
                 MPI_Comm     commgroup,
                 int          n_samples )
{

	double tmpTimeStart = timeStartSimple();

	int    			i;
	int 				procs, myid;
	int  	  		 	dataSize;
	int				global_N;
	int  			  	*send_index;
	int  	  			*send_count, *recv_count;
	int				total_recv_count;
	int				max_alloc_outbuf_size;
	int 				*actual_count;
	int 				*expected_count;
	keyType			*sampleKeyArray_local;
	keyType			*sampleKeyArray_all;
	keyType			*splitterArray;
	type				*resultBuf;

	MPI_Comm_size(commgroup, &procs);
	MPI_Comm_rank(commgroup, &myid);
	MPI_Type_size(dataType, &dataSize);

	/* local in-place sort */
	double tmpTimeStart2 = timeStartSimple();
	qsort( buf, *local_N, sizeof(type), compare_type );
	timeEndSimple(tmpTimeStart2, &t_sort_lsort1);

	/* some stars are destroyed during the timestep, and their r values are set to infinity. Using this, we here remove them, and fix local_N to account for these lost stars. */
	tmpTimeStart2 = timeStartSimple();
	remove_stripped_stars(buf, local_N);

	/* find total number of elements to be sorted in parallel */
	double tmpTimeStart3 = timeStartSimple();
	MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* find the expected number of elements each processor should have as per the data partitioninig scheme */
	expected_count = (int*) malloc(procs * sizeof(int));
	find_expected_count( expected_count, global_N, procs );
	timeEndSimple(tmpTimeStart2, &t_sort_oth);

	/* Picking samples from the local data set */
	tmpTimeStart2 = timeStartSimple();
	sampleKeyArray_local = (keyType*) malloc(n_samples * sizeof(keyType));
	sample(buf, sampleKeyArray_local, *local_N, n_samples);

	/* root node gathers the samples from all nodes */
	sampleKeyArray_all = (keyType*) malloc(procs * n_samples * sizeof(keyType));
	tmpTimeStart3 = timeStartSimple();
	MPI_Gather(sampleKeyArray_local, n_samples * sizeof(keyType), MPI_BYTE, sampleKeyArray_all, n_samples * sizeof(keyType), MPI_BYTE, 0, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* procs-1 numbers are enough for p buckets (1 bucket/processor) */
	splitterArray = (keyType*) malloc((procs-1) * sizeof(keyType));

	/* Sorting the collected samples and determining splitters */
	if(myid==0)
	{
		qsort( sampleKeyArray_all, procs*n_samples, sizeof(keyType), compare_keyType );

		for(i=0; i<procs-1; i++)
			splitterArray[i] = sampleKeyArray_all[ (i+1) * n_samples - 1 ];
	}

	/* sending back splitters to all nodes */
	tmpTimeStart3 = timeStartSimple();
	MPI_Bcast(splitterArray, (procs-1) * sizeof(keyType), MPI_BYTE, 0, MPI_COMM_WORLD);
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* find the offset index for each send using binary search on splitter array */
	send_index = (int *) calloc(procs, sizeof(int));
	send_index[0] = 0;

	for (i=1; i<procs; i++)
		send_index[i] = binary_search(buf, splitterArray[i-1], 0, (*local_N)-1);

	free(sampleKeyArray_all);
	free(sampleKeyArray_local);
	free(splitterArray);
	timeEndSimple(tmpTimeStart2, &t_sort_splitters);

	tmpTimeStart2 = timeStartSimple();
	/* find send/recv count for each send/recv */
	send_count = (int*) malloc(procs * sizeof(int));
	recv_count = (int*) malloc(procs * sizeof(int));
	for (i=0; i<procs-1; i++)
		send_count[i] = send_index[i+1] - send_index[i];

	send_count[procs-1] = (*local_N) - send_index[procs-1]; // + 1;

	/* exchange the send counts to know how many stars are to be received */
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* calculate total number of stars to be received on this node */
	total_recv_count = 0;
	for (i=0; i<procs; i++) total_recv_count += recv_count[i];

	/* calculate the receive displacements for the mpi communication */
	int* recv_displ = (int*) malloc(procs * sizeof(int));
	recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		recv_displ[i] = recv_displ[i-1] + recv_count[i-1];

	/* find the actual counts on each processor - to be used later too */
	actual_count = (int*) malloc(procs * sizeof(int));
	tmpTimeStart3 = timeStartSimple();
	MPI_Allgather( &total_recv_count, 1, MPI_INT, actual_count, 1, MPI_INT, commgroup );
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* finding the maximum size of resultBuf to be allocated */
	max_alloc_outbuf_size = 0;
	for(i=0; i<procs; i++)
		max_alloc_outbuf_size = (actual_count[i] > max_alloc_outbuf_size) ? actual_count[i] : max_alloc_outbuf_size;

	/* allocate recv buffer */
	resultBuf = (type*) malloc(actual_count[myid] * sizeof(type)); 

	/* all to all communication */
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoallv(buf, send_count, send_index, dataType, resultBuf, recv_count, recv_displ, dataType, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	dprintf("new no.of stars in proc %d = %d\n", myid, total_recv_count);




	/***** Binary data *****/
	/* Now, we also have to send the binaries to the appropriate processor. For each bucket, we go over the single stars, and using their binind values, pack together all binaries that need to go to the appropriate processor. Once this is done for all buckets in each processor, we do an all to all communication and we're done. We also make sure the binind values are set appropriately. */
	int* b_send_index = (int*) malloc(procs * sizeof(int));
	int* b_send_count = (int*) malloc(procs * sizeof(int));
	int* b_recv_count = (int*) malloc(procs * sizeof(int));
/* buffer to pack binaries such that the binaries of each bucket are packet together and in sequence*/
	binary_t* b_tmp_buf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	int j, k=0;
	//Iterate over each bucket
	for(i=0; i<procs;i++)
	{
		b_send_index[i] = k;
		//Go over the single stars in the bucket
		for(j=send_index[i]; j<send_index[i]+send_count[i]; j++)
		{
			//If it is a binary, copy it to the buffer
			if(buf[j].binind > 0)
			{
				//MPI: binind-1 is used since binary+1 is passed as input to this sorting function.
				memcpy(&b_tmp_buf[k], &b_buf[buf[j].binind-1], sizeof(binary_t));
				k++;
			}
		}
		//Set the send count for this bucket
		b_send_count[i] = k - b_send_index[i];
	}

	//MPI: Set binary array to zeros, if not might cause problems when new stars are created in the next timestep. So it's best to wipe out the older data.
	memset (b_buf, 0, (N_BIN_DIM_OPT-1) * sizeof(binary_t));

	/* exchange send counts to know how many binaries to receive from each node */
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	/* total receive count for this node */
	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];

	/* receive displacement for all to all communication */
	int* b_recv_displ = (int*) malloc(procs * sizeof(int));
	b_recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		b_recv_displ[i] = b_recv_displ[i-1] + b_recv_count[i-1];

	binary_t* b_resultBuf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	/* all to all communication */
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoallv(b_tmp_buf, b_send_count, b_send_index, b_dataType, b_resultBuf+1, b_recv_count, b_recv_displ, b_dataType, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	//MPI: Before we do local sort, we need to fix the binary addressing i.e. the binind values.
	//MPI: k starts from 1 because binind has to be > 0 for binaries. and also the 0th element in the binary array is not to be used.
	k=1; 
	int kprev;
	for(i=0; i<procs; i++)
	{
		kprev=k;
		for(j=recv_displ[i]; j<recv_displ[i]+recv_count[i]; j++)
		{
			if(resultBuf[j].binind > 0)
			{
				resultBuf[j].binind = k;
				k++;
			}
		}
		//MPI: Checks to make sure the right number of binaries are recd.
		if(k-kprev!=b_recv_count[i]) eprintf("mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-kprev, b_recv_count[i]);
	}
	//MPI: Additional checks.
	if(k-1!=b_total_recv_count)
		eprintf("Binary numbers mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-1, b_total_recv_count);

	/***** End binary data *****/
	timeEndSimple(tmpTimeStart2, &t_sort_a2a);




	/* merge chunks recieved and local sort */
	tmpTimeStart2 = timeStartSimple();
	qsort(resultBuf, total_recv_count, sizeof(type), compare_type);
	timeEndSimple(tmpTimeStart2, &t_sort_lsort2);
	timeEndSimple(tmpTimeStart, &t_sort_only);


	/* exchange stars between processors to stay consistent with data partitioning scheme */
	tmpTimeStart = timeStartSimple();
	load_balance(resultBuf, buf, b_resultBuf, b_buf, expected_count, actual_count, myid, procs, dataType, b_dataType, commgroup);
	timeEndSimple(tmpTimeStart, &t_sort_lb);

	tmpTimeStart = timeStartSimple();
	tmpTimeStart2 = timeStartSimple();
	MPI_Barrier(commgroup);

	*local_N = actual_count[myid];

	free(expected_count);
	free(actual_count);
	free(send_index);
	free(send_count);
	free(recv_count);
	free(recv_displ);
	free(resultBuf);

	free(b_resultBuf);
	free(b_send_index);
	free(b_send_count);
	free(b_recv_count);
	free(b_recv_displ);
	free(b_tmp_buf);
	timeEndSimple(tmpTimeStart, &t_sort_oth);
	timeEndSimple(tmpTimeStart2, &t_sort_only);

	return global_N;
}

/**
* @brief after sample sort, exchanges data between processors to make sure the number of stars in each processor is in accordance with the data partitioning scheme
*/
void load_balance( 	type 				*inbuf,
							type 				*outbuf,
							binary_t			*b_inbuf,
							binary_t			*b_outbuf,
							int 				*expected_count, 
							int				*actual_count, 
							int 				myid, 
							int 				procs, 
							MPI_Datatype 	dataType, 	
							MPI_Datatype 	b_dataType, 	
							MPI_Comm			commgroup	)
{

	int 	i, local_count, total_recv_count;
	int 	*expected_cum_count, *actual_cum_count, *splitter;
	int  	*send_index;
	int  	*send_count, *recv_count;

	local_count = actual_count[myid];

	expected_cum_count = (int*) malloc(procs * sizeof(int));
	actual_cum_count = (int*) malloc(procs * sizeof(int));
	splitter = (int*) malloc(procs * sizeof(int));
	expected_cum_count[0] = 0;
	actual_cum_count[0] = 0;
	splitter[0] = 0;

	//Finding the expected and actual cumulative counts until the processors before this processor
	for(i=1; i<procs; i++)
	{
		expected_cum_count[i] = expected_count[i-1] + expected_cum_count[i-1];
		actual_cum_count[i] = actual_count[i-1] + actual_cum_count[i-1];
		//If you imagine an array putting together all the data, then this would be the splitter that would ideally divide this array into near-equal pieces (or as per the data partitioning scheme)
		splitter[i] = expected_cum_count[i];
	}

	send_index = (int *) calloc(procs, sizeof(int));
	send_count = (int*) calloc(procs, sizeof(int));
	recv_count = (int*) malloc(procs * sizeof(int));

	//check if splitter of 0 can be in any other proc that 0.
	for(i=1; i<procs; i++)
	{
		if(local_count == 0) break;

		//Based on the relative position between the each splitter and actual cumulative count for this processor, we decide which chunk of data goes to which processor. For a not so bad distribution of stars among processors, this will essentially be communication between neighbors, but in worse cases, this might need much more than that - an all to all communication, which is what we do here.
		if(splitter[i] >= actual_cum_count[myid] && splitter[i] < actual_cum_count[myid] + local_count)
		{
			send_index[i] = splitter[i] - actual_cum_count[myid];
			send_count[i-1] = send_index[i] - send_index[i-1];
		}

		if (splitter[i] >= actual_cum_count[myid] + local_count)
			break;
	}

	//figure out send and receive parameters
	send_count[i-1] = local_count - send_index[i-1];

	double tmpTimeStart3 = timeStartSimple();
	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	total_recv_count = 0;
	for (i=0; i<procs; i++) total_recv_count += recv_count[i];
	dprintf("after load-balancing stars in proc %d = %d expected_count = %d\n", myid, total_recv_count, expected_count[myid]);

	int* recv_displ = (int*) malloc(procs * sizeof(int));
	recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		recv_displ[i] = recv_displ[i-1] + recv_count[i-1];

	//all to all communication
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoallv(inbuf, send_count, send_index, dataType, outbuf, recv_count, recv_displ, dataType, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	actual_count[myid] = total_recv_count;




	/***** Binary data *****/
	//binary data is handled similar to th sorting routine
	int* b_send_index = (int*) malloc(procs * sizeof(int));
	int* b_send_count = (int*) malloc(procs * sizeof(int));
	int* b_recv_count = (int*) malloc(procs * sizeof(int));
	binary_t* b_tmp_buf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	int j, k=0;
	for(i=0; i<procs;i++)
	{
		b_send_index[i] = k;
		for(j=send_index[i]; j<send_index[i]+send_count[i]; j++)
			if(inbuf[j].binind > 0)
			{
				memcpy(&b_tmp_buf[k], &b_inbuf[inbuf[j].binind], sizeof(binary_t));
				k++;
			}
		b_send_count[i] = k - b_send_index[i];
	}
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];
	N_b_local = b_total_recv_count;

	int* b_recv_displ = (int*) malloc(procs * sizeof(int));
	b_recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		b_recv_displ[i] = b_recv_displ[i-1] + b_recv_count[i-1];

	//MPI: All to All Communication
	tmpTimeStart3 = timeStartSimple();
	MPI_Alltoallv(b_tmp_buf, b_send_count, b_send_index, b_dataType, b_outbuf, b_recv_count, b_recv_displ, b_dataType, commgroup);
	timeEndSimple(tmpTimeStart3, &t_comm);

	//MPI: k starts from 1 because binind has to be > 0 for binaries. and also the 0th element in the binary array is not to be used.
	k=1;
	int kprev;
	for(i=0; i<procs; i++)
	{
		kprev=k;
		for(j=recv_displ[i]; j<recv_displ[i]+recv_count[i]; j++)
		{
			if(outbuf[j].binind > 0)
			{
				outbuf[j].binind = k;
				k++;
			}
		}
		if(k-kprev!=b_recv_count[i]) eprintf("mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-kprev, b_recv_count[i]);
	}

	if(k-1!=b_total_recv_count)
		eprintf("Binary numbers mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-1, b_total_recv_count);

	free(b_send_index);
	free(b_send_count);
	free(b_recv_count);
	free(b_recv_displ);
	free(b_tmp_buf);

	/***** End binary data *****/




	free(splitter);
	free(send_index);
	free(send_count);
	free(recv_count);
	free(recv_displ);
	free(actual_cum_count);
	free(expected_cum_count);
}
#endif
