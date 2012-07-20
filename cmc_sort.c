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

keyType getKey(type* buf) { return (buf->r); }

int compare_keyType (const void * a, const void * b) 
{
	if( *((keyType*)a) < *((keyType*)b) ) return -1;
	if( *((keyType*)a) > *((keyType*)b) ) return 1;
	return 0;
}

int compare_type (const void * a, const void * b)
{
	if( getKey((type*)a) < getKey((type*)b) ) return -1;
	if( getKey((type*)a) > getKey((type*)b) ) return 1;
	return 0;
}


void find_expected_count( int* expected_count, int N, int numprocs )
{
	findLimits( N, 20 );
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

keyType* sample(type *buf, keyType *sample_array, int N, int n_samples)
{
	//srand ( time(NULL) );
	if(n_samples > N)
		eprintf("Oversampling occurred! id = %d num samples = %d local no. of stars = %d", myid, n_samples, N);

	int i;
	for(i=1; i<=n_samples; i++)
		// RANDOM
		// if(i==1) srand ( time(NULL) );
		//buf[ rand() % N + 1 ].key;
		// UNIFORM
		sample_array[i-1] = getKey(&buf[ i * N/(n_samples+1) ]);

	return sample_array;
}

//Re-write using compare fn.
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
int sample_sort( type			*buf,
                 int	         *local_N,
                 MPI_Datatype dataType,
					  binary_t		*b_buf,
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

	double tmpTimeStart2 = timeStartSimple();
	/* local in-place sort */
	qsort( buf, *local_N, sizeof(type), compare_type );

	// Fixing local_N to account for lost stars.
	remove_stripped_stars(buf, local_N);

	// Find global total number of elements to be sorted in parallel
	MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);

	expected_count = (int*) malloc(procs * sizeof(int));
	find_expected_count( expected_count, global_N, procs );
	timeEndSimple(tmpTimeStart2, &t_sort1);

	tmpTimeStart2 = timeStartSimple();
	/* Picking random/uniform samples and sending to root */
	sampleKeyArray_local = (keyType*) malloc(n_samples * sizeof(keyType));
	sample(buf, sampleKeyArray_local, *local_N, n_samples);

	/* Sample arrays to collect samples from each node and send to root */
	sampleKeyArray_all = (keyType*) malloc(procs * n_samples * sizeof(keyType));
	MPI_Gather(sampleKeyArray_local, n_samples * sizeof(keyType), MPI_BYTE, sampleKeyArray_all, n_samples * sizeof(keyType), MPI_BYTE, 0, commgroup);

	/* procs-1 numbers are enough to determine which elements belong to which bucket/proc. */
	splitterArray = (keyType*) malloc((procs-1) * sizeof(keyType));

	/* Sorting the collected samples, and sending back splitters */
	if(myid==0)
	{
		qsort( sampleKeyArray_all, procs*n_samples, sizeof(keyType), compare_keyType );

		for(i=0; i<procs-1; i++)
			splitterArray[i] = sampleKeyArray_all[ (i+1) * n_samples - 1 ];
	}

	MPI_Bcast(splitterArray, (procs-1) * sizeof(keyType), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* find the offset index for each send using binary search on splitter array */
	send_index = (int *) calloc(procs, sizeof(int));
	send_index[0] = 0;

	for (i=1; i<procs; i++)
		send_index[i] = binary_search(buf, splitterArray[i-1], 0, (*local_N)-1);

	free(sampleKeyArray_all);
	free(sampleKeyArray_local);
	free(splitterArray);

	/* find send/recv count for each send/recv */
	send_count = (int*) malloc(procs * sizeof(int));
	recv_count = (int*) malloc(procs * sizeof(int));
	for (i=0; i<procs-1; i++)
		send_count[i] = send_index[i+1] - send_index[i];

	send_count[procs-1] = (*local_N) - send_index[procs-1]; // + 1;

	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);

	total_recv_count = 0;
	for (i=0; i<procs; i++) total_recv_count += recv_count[i];

	int* recv_displ = (int*) malloc(procs * sizeof(int));
	recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		recv_displ[i] = recv_displ[i-1] + recv_count[i-1];

	actual_count = (int*) malloc(procs * sizeof(int));
	MPI_Allgather( &total_recv_count, 1, MPI_INT, actual_count, 1, MPI_INT, commgroup );
	timeEndSimple(tmpTimeStart2, &t_sort2);
	tmpTimeStart2 = timeStartSimple();

	//Finding the maximum size of resultBuf to be allocated.
	max_alloc_outbuf_size = 0;
	for(i=0; i<procs; i++)
		max_alloc_outbuf_size = (actual_count[i] > max_alloc_outbuf_size) ? actual_count[i] : max_alloc_outbuf_size;

	/* allocate recv buffer */
	//MPI3: Hopefully allocating actual_count[myid] elements would be sufficient.
	resultBuf = (type*) malloc(actual_count[myid] * sizeof(type)); 

	//MPI3: All to all communication
	MPI_Alltoallv(buf, send_count, send_index, dataType, resultBuf, recv_count, recv_displ, dataType, commgroup);

	dprintf("new no.of stars in proc %d = %d\n", myid, total_recv_count);




	/***** Binary data *****/
	int* b_send_index = (int*) malloc(procs * sizeof(int));
	int* b_send_count = (int*) malloc(procs * sizeof(int));
	int* b_recv_count = (int*) malloc(procs * sizeof(int));
	binary_t* b_tmp_buf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	int j, k=0;
	for(i=0; i<procs;i++)
	{
		b_send_index[i] = k;
		for(j=send_index[i]; j<send_index[i]+send_count[i]; j++)
		{
			if(buf[j].binind > 0)
			{
				//MPI3: binind-1 is used as binary+1 is passed as input to this sorting function.
				memcpy(&b_tmp_buf[k], &b_buf[buf[j].binind-1], sizeof(binary_t));
				k++;
			}
		}
		b_send_count[i] = k - b_send_index[i];
	}

	//MPI3: Set binary array to zeros, if not might cause problems when new stars are created in the next timestep. So it's best to wipe out the older data.
	memset (b_buf, 0, (N_BIN_DIM_OPT-1) * sizeof(binary_t));

	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);

	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];

	int* b_recv_displ = (int*) malloc(procs * sizeof(int));
	b_recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		b_recv_displ[i] = b_recv_displ[i-1] + b_recv_count[i-1];

	binary_t* b_resultBuf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	//MPI3: All to All Communication
	MPI_Alltoallv(b_tmp_buf, b_send_count, b_send_index, b_dataType, b_resultBuf+1, b_recv_count, b_recv_displ, b_dataType, commgroup);

	//MPI3: Before we do local sort, we need to fix the binary addressing.
	//Test thoroughly. Before alltoall use the id value of star elements to one of the ids of binaries, and check post sort - for testing.
	//MPI3: k starts from 1 because binind has to be > 0 for binaries. and also the 0th element in the binary array is not to be used.
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
		//MPI3: Checks to make sure the right number of binaries are recd.
		if(k-kprev!=b_recv_count[i]) eprintf("mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-kprev, b_recv_count[i]);
	}
	//MPI3: Additional checks.
	if(k-1!=b_total_recv_count)
		eprintf("Binary numbers mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-1, b_total_recv_count);

	/***** End binary data *****/




	timeEndSimple(tmpTimeStart2, &t_sort3);
	tmpTimeStart2 = timeStartSimple();
	qsort(resultBuf, total_recv_count, sizeof(type), compare_type);
	timeEndSimple(tmpTimeStart2, &t_sort4);
	timeEndSimple(tmpTimeStart, &t_sort_only);



	tmpTimeStart = timeStartSimple();
	load_balance(resultBuf, buf, b_resultBuf, b_buf, expected_count, actual_count, myid, procs, dataType, b_dataType, commgroup);
	timeEndSimple(tmpTimeStart, &t_load_bal);

	MPI_Barrier(commgroup);

	tmpTimeStart = timeStartSimple();
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
	timeEndSimple(tmpTimeStart, &t_sort_only);

	return global_N;
}


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

	for(i=1; i<procs; i++)
	{
		expected_cum_count[i] = expected_count[i-1] + expected_cum_count[i-1];
		actual_cum_count[i] = actual_count[i-1] + actual_cum_count[i-1];
		splitter[i] = expected_cum_count[i];
	}

	send_index = (int *) calloc(procs, sizeof(int));
	send_count = (int*) calloc(procs, sizeof(int));
	recv_count = (int*) malloc(procs * sizeof(int));

	//check if splitter of 0 can be in any other proc that 0.
	for(i=1; i<procs; i++)
	{
		if(local_count == 0) break;

		if(splitter[i] >= actual_cum_count[myid] && splitter[i] < actual_cum_count[myid] + local_count)
		{
			send_index[i] = splitter[i] - actual_cum_count[myid];
			send_count[i-1] = send_index[i] - send_index[i-1];
		}

		if (splitter[i] >= actual_cum_count[myid] + local_count)
			break;
	}

	send_count[i-1] = local_count - send_index[i-1];

	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);

	total_recv_count = 0;
	for (i=0; i<procs; i++) total_recv_count += recv_count[i];
	dprintf("after load-balancing stars in proc %d = %d expected_count = %d\n", myid, total_recv_count, expected_count[myid]);

	int* recv_displ = (int*) malloc(procs * sizeof(int));
	recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		recv_displ[i] = recv_displ[i-1] + recv_count[i-1];

	MPI_Alltoallv(inbuf, send_count, send_index, dataType, outbuf, recv_count, recv_displ, dataType, commgroup);

	actual_count[myid] = total_recv_count;




	/***** Binary data *****/
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
	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);

	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];
	N_b_local = b_total_recv_count;

	int* b_recv_displ = (int*) malloc(procs * sizeof(int));
	b_recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		b_recv_displ[i] = b_recv_displ[i-1] + b_recv_count[i-1];

	//MPI3: All to All Communication
	MPI_Alltoallv(b_tmp_buf, b_send_count, b_send_index, b_dataType, b_outbuf, b_recv_count, b_recv_displ, b_dataType, commgroup);

	//MPI3: k starts from 1 because binind has to be > 0 for binaries. and also the 0th element in the binary array is not to be used.
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
