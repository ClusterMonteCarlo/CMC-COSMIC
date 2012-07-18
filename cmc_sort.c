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
double getKey_old(star_t* starData) { return (starData->r); }

int compare (const void * a, const void * b) { return ( *(double*)a > *(double*)b ) ? 1 : -1; }

int binary_search_old( star_t* starData, double r, int kmin, int kmax )
{
	long ktry;

	if (starData[kmin].r > r  || starData[kmax].r < r)
		printf("proc %d: r=%g is outside kmin=%g kmax=%g!!\n", myid, r, starData[kmin].r, starData[kmax].r);
	if (starData[kmin].r > r)
		return kmin;
	if (starData[kmax].r < r)
		return kmax+1;

	do {      
		ktry = (kmin+kmax+1)/2;
		if (starData[ktry].r<r)
		{
			kmin = ktry;
		} else {
			kmax = ktry-1;
		}
	} while (kmax!=kmin);

	return kmin+1;
}

/* Tasks:
1. Implement uniform/random sampling option.
2. Post-sort Communication to match distribution requirement - done
3. Optimize for load-balancing - my method.
4. Use binning during sampling (Stefan's suggestion)
*/
int sample_sort_old( star_t        *starData,
                  long          *local_N, //Pass clus.N_MAX_NEW
                  MPI_Datatype  eType,
                  MPI_Comm      commgroup,
                  int           num_samples )
{

	int         i;
	int         eSize;
	int 			global_N;
	int         *send_index;
	int         *send_count, *recv_count;
	int			total_recv_count;
	int 			s_count_back, r_count_back, s_count_fwd, r_count_fwd;
	int 			*actual_count;
	int 			*actual_cum_count;
	int 			*expected_cum_count;
	double		*sampleKeyArray_local;
	double		*sampleKeyArray_all;
	double		*splitterArray;
	star_t      *resultBuf;
	star_t      *recvBuf;
	MPI_Status  *sendStatus, *recvStatus;
	MPI_Request *sendReq, *recvReq;
	MPI_Status	stat_back, stat_fwd;
#ifdef EXPERIMENTAL
	int         ideal_count;
#endif

	MPI_Type_size(eType, &eSize);

	/* local in-place sort */
	qsorts(starData+1, *local_N);

	dprintf("before stripping id = %d localN = %ld\n", myid, *local_N);
	// Fixing local_N to account for lost stars.
	i = *local_N;
	while(star[i].r == SF_INFINITY)
	{
		dmpiprintf("stripped star found and removed. i = %d id = %ld from proc. %d, local_N = %ld \n", i, star[i].id, myid, *local_N);
		zero_star(i); //Check with Stefan if this is needed.
		(*local_N)--;
		i--;
	}
	dprintf("after stripping id = %d localN = %ld\n", myid, *local_N);

	//find global total number of elements to be sorted in parallel
	MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);
	rootprintf("A total of %d stars to be sorted in parallel\n", global_N);

	//mpiFindIndicesCustom( global_N, 20, myid, &mpiBegin, &mpiEnd );
	findLimits( global_N, 20 );

	int *expected_count = (int*) malloc(procs * sizeof(int));
	for(i=0; i<procs; i++)
		expected_count[i] = End[i] - Start[i] + 1;

	//MPI3: If we use time as seed, reproducibility will become a problem. Using the timestep count as random seed instead to have both some randomness and reproducability..
	//srand ( time(NULL) );
	if(tcount==1)
		srand ( tcount );

	/* Sample arrays to collect samples from each node and send to root */
	sampleKeyArray_local = (double*) malloc(num_samples * sizeof(double));
	sampleKeyArray_all = (double*) malloc(procs * num_samples * sizeof(double));
	/* procs-1 numbers are enough to determine which elements belong ti which bucket/proc. */
	splitterArray = (double*) malloc((procs-1) * sizeof(double));

	/* Picking random samples and sending to root */
	for(i=0; i<num_samples; i++)
		sampleKeyArray_local[i] = starData[ rand() % *local_N + 1 ].r;
/*
	MPI_Datatype sampleType;
	MPI_Type_contiguous( sizeof(sampleKeyArray_local), MPI_BYTE, &sampleType );
	MPI_Type_commit( &sampleType );
*/
	MPI_Gather(sampleKeyArray_local, num_samples, MPI_DOUBLE, sampleKeyArray_all, num_samples, MPI_DOUBLE, 0, commgroup);
	//MPI_Gather(sampleKeyArray_local, num_samples, sampleType, sampleKeyArray_all, num_samples, sampleType, 0, commgroup);
	//MPI_Type_free(sampleType);

	/* Sorting the collected samples, and sending back splitters */
	if(myid==0)
	{
		qsort( sampleKeyArray_all, procs*num_samples, sizeof(double), compare );

		//dprintf("Splitter array is: ");
		for(i=0; i<procs-1; i++)
		{
			splitterArray[i] = sampleKeyArray_all[ (i+1) * num_samples - 1 ];
			//dprintf("%g ", splitterArray[i]);
		}
		//dprintf("\n\n");
	}

	MPI_Bcast(splitterArray, procs-1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	free(sampleKeyArray_all);
	free(sampleKeyArray_local);

	/* find the offset index for each send using binary search on splitter array */
	send_index = (int *) calloc(procs, sizeof(int));
	send_index[0] = 1; //0 is sentinel, so set to 1

	for (i=1; i<procs; i++)
		send_index[i] = binary_search_old(starData, splitterArray[i-1], 1, *local_N);

	/* find send/recv count for each send/recv */
	send_count = (int*) malloc(procs * sizeof(int));
	recv_count = (int*) malloc(procs * sizeof(int));
	for (i=0; i<procs-1; i++)
	{
		send_count[i] = send_index[i+1] - send_index[i];
		if(send_index[i] == -1)
		{
			send_index[i] = 0;
			send_count[i] = 0;
		}
	}
	send_count[i] = *local_N - send_index[i] + 1;
	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);

#ifdef EXPERIMENTAL
	/* Experimental: To improve load balancing */

	ideal_count = global_N / procs;

	//find total recv number
	total_recv_count = recv_count[0];
	for (i=1; i<procs; i++) total_recv_count += recv_count[i];

	/* End Experimental */
#endif

	/* allocate recv buffer */
	total_recv_count = recv_count[0];
	for (i=1; i<procs; i++) total_recv_count += recv_count[i];

	// Watch out! If the no.of stars after sort is very uneven, there is a possibility of it exceeding the allocated memore. In that case load_bal might have to be increased.
	double load_bal = 4;
	resultBuf = (star_t*) malloc((int)floor((double)expected_count[myid] * load_bal) * sizeof(star_t)); //The multiplying factor load_bal is for collecting data from neighbors fix load imbalance due to parallel sorting.
	//resultBuf = (star_t*) malloc(total_recv_count * sizeof(star_t));
	//dprintf("Allocating memory to receive %d stars on proc %d\n", (int)floor((double)expected_count[myid]*load_bal), myid);
	//dprintf("proc %d: local_N = %ld\n", myid, *local_N);

	dmpiprintf("Proc %d will receive: %d \t Ideal: %d\n", myid, send_count[myid]+total_recv_count, expected_count[myid]);
	MPI_Barrier(commgroup);

	//for(i=0; i<procs; i++)
		//dprintf("proc %d:\tsend_count[%d] = %d\tsend_index[%d] = %d\trecv_count[%d] = %d \n", myid, i, send_count[i], i, send_index[i], i, recv_count[i]);

	/* allocate asynchronous I/O request and wait status */
	sendStatus = (MPI_Status*)  malloc(2*procs* sizeof(MPI_Status));
	sendReq    = (MPI_Request*) malloc(2*procs* sizeof(MPI_Request));
	recvStatus = sendStatus + procs;
	recvReq    = sendReq    + procs;

	/* asynchronous sends */
	for (i=0; i<procs; i++) {
		sendReq[i] = MPI_REQUEST_NULL;
		if (myid != i && send_count[i] > 0)
			MPI_Isend(starData + send_index[i],
					send_count[i],
					eType,
					i,
					0,
					commgroup,
					&sendReq[i]);
	}

	/* asynchronous recv */
	recvBuf = resultBuf;
	for (i=0; i<procs; i++) {
		recvReq[i] = MPI_REQUEST_NULL;
		if (myid == i) {
			memcpy(recvBuf, starData+send_index[i], send_count[i]*eSize);
		}
		else if (recv_count[i] > 0) {
			MPI_Irecv(recvBuf,
					recv_count[i],
					eType,
					i,
					0,
					commgroup,
					&recvReq[i]);
		}
		recvBuf += recv_count[i];// * eSize;
	}

	/* wait till all asynchronous I/O done */
	MPI_Waitall(2*procs, sendReq, sendStatus);

	//dprintf("New no.of stars in proc %d = %d\n", myid, total_recv_count);

	/* merge sort might be faster */
	qsorts(resultBuf, total_recv_count);

	/* Post-sort communication for load balancing */
	// Variables with _back refer too backward communication, and fwd to forward communication.
	s_count_back = 0; r_count_back = 0; s_count_fwd = 0; r_count_fwd = 0;
	actual_count = (int*) malloc(procs * sizeof(int));
	actual_cum_count = (int*) malloc(procs * sizeof(int));
	expected_cum_count = (int*) malloc(procs * sizeof(int));

	MPI_Allgather( &total_recv_count, 1, MPI_INT, actual_count, 1, MPI_INT, commgroup );

	int tmp1 = 0, tmp2 = 0;
	for(i=0; i<procs; i++)
	{
		expected_cum_count[i] = tmp1 + expected_count[i];
		actual_cum_count[i] = tmp2 + actual_count[i];
		tmp1 = expected_cum_count[i];
		tmp2 = actual_cum_count[i];
	}

	dmpiprintf("IDCHECK %d\n", myid);
	if(myid != 0)
	{
		s_count_back = expected_cum_count[myid-1] - actual_cum_count[myid-1];
		r_count_fwd = -s_count_back;
	}
	s_count_back = (s_count_back > 0) ? s_count_back : 0;
	r_count_fwd = (r_count_fwd > 0) ? r_count_fwd : 0;

	if(myid != procs-1)
	{
		r_count_back = expected_cum_count[myid] - actual_cum_count[myid];
		s_count_fwd = -r_count_back;
	}
	r_count_back = (r_count_back > 0) ? r_count_back : 0;
	s_count_fwd = (s_count_fwd > 0) ? s_count_fwd : 0;

	// Backward Communication
	if( s_count_back > 0 )
		MPI_Send(resultBuf, s_count_back * sizeof(star_t), MPI_BYTE, myid - 1, 0, commgroup);

	if( r_count_back > 0 )
		MPI_Recv(resultBuf + total_recv_count, r_count_back * sizeof(star_t), MPI_BYTE, myid + 1, 0, commgroup, &stat_back);

	// Forward Communication
	if( s_count_fwd > 0 )
		MPI_Send(resultBuf + total_recv_count - s_count_fwd, s_count_fwd * sizeof(star_t), MPI_BYTE, myid + 1, 0, commgroup);

	if( r_count_fwd > 0 )
		MPI_Recv(star + 1, r_count_fwd * sizeof(star_t), MPI_BYTE, myid - 1, 0, commgroup, &stat_fwd);

	//int count_back, count_fwd;
	//MPI_Get_count(&stat_back, MPI_BYTE, &count_back);
	//MPI_Get_count(&stat_fwd, MPI_BYTE, &count_fwd);
	//dprintf("%d sent_fwd = %d sent_back = %d recd_fwd = %ld recd_back = %ld\n", myid, s_count_fwd, s_count_back, count_back/sizeof(star_t), count_fwd/sizeof(star_t));
	total_recv_count += r_count_fwd + r_count_back - s_count_fwd - s_count_back;

	//dprintf("%d sent_fwd = %d sent_back = %d recd_fwd = %d recd_back = %d\n", myid, s_count_fwd, s_count_back, r_count_fwd, r_count_back);
	//dprintf("%d total_recv_count = %d\n", myid, total_recv_count);

	//MPI3: Copy everything back to star array - make sure enough memory is allocated for star in each node. Retain sentinel - 0th star.
	for (i=0; i<total_recv_count; i++)
		starData[i+1+r_count_fwd] = resultBuf[i + s_count_back];

	*local_N = total_recv_count;

	MPI_Barrier(commgroup);
	free(send_index);
	free(send_count);
	free(recv_count);
	free(sendStatus);
	free(sendReq);
	free(splitterArray);

	free(resultBuf);
	free(actual_count);
	free(actual_cum_count);
	free(expected_count);
	free(expected_cum_count);

	return global_N;
}


// ************************************************************** //
// NEWER VERSION OF SAMPLE SORT WITH CLEANER LOAD BALANCING FUNCTION
// ************************************************************** //
//typedef star_t type;

keyType getKey(type* buf) { return (buf->r); }

int compare_keyType (const void * a, const void * b) 
{
	//return ( *((keyType*)a) - *((keyType*)b) );
	if( *((keyType*)a) < *((keyType*)b) ) return -1;
	if( *((keyType*)a) > *((keyType*)b) ) return 1;
	return 0;
}

int compare_type (const void * a, const void * b)
{
	//return ( getKey((type*)a) - getKey((type*)b) );
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
      //printf("proc %d: r=%g is outside kmin=%g kmax=%g!!\n", myid, r, getKey(&buf[kmin]), getKey(buf[kmax]));
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
	type  	    	*recvBuf;
	MPI_Status  	*sendStatus, *recvStatus;
	MPI_Request 	*sendReq, *recvReq;

	MPI_Comm_size(commgroup, &procs);
	MPI_Comm_rank(commgroup, &myid);
	MPI_Type_size(dataType, &dataSize);

	double tmpTimeStart2 = timeStartSimple();
	/* local in-place sort */
	qsort( buf, *local_N, sizeof(type), compare_type );

	//printf("before stripping id = %d localN = %d\n", myid, *local_N);
	// Fixing local_N to account for lost stars.
	remove_stripped_stars(buf, local_N);
	//printf("after stripping id = %d localN = %d\n", myid, *local_N);

	// Find global total number of elements to be sorted in parallel
	MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);
	//printf("A total of %d stars to be sorted in parallel\n", global_N);

	//mpiFindIndicesCustom( global_N, 20, myid, &mpiBegin, &mpiEnd );
	//findLimits( global_N, 20 );
	expected_count = (int*) malloc(procs * sizeof(int));
	find_expected_count( expected_count, global_N, procs );
/*
	if(myid==0)
		for(i=0; i<procs; i++)
			printf("expected_count[%d] = %d\n",i, expected_count[i]);
*/
	timeEndSimple(tmpTimeStart2, &t_sort1);
	tmpTimeStart2 = timeStartSimple();

	/* Picking random/uniform samples and sending to root */
	sampleKeyArray_local = (keyType*) malloc(n_samples * sizeof(keyType));
	sample(buf, sampleKeyArray_local, *local_N, n_samples);
/*
	if(myid ==3)
		for(i=0; i<n_samples; i++)
			printf("%f\n", sampleKeyArray_local[i]);	
*/
	/* Sample arrays to collect samples from each node and send to root */
	sampleKeyArray_all = (keyType*) malloc(procs * n_samples * sizeof(keyType));
	MPI_Gather(sampleKeyArray_local, n_samples * sizeof(keyType), MPI_BYTE, sampleKeyArray_all, n_samples * sizeof(keyType), MPI_BYTE, 0, commgroup);
/*
	if(myid ==0)
		for(i=0; i<procs*n_samples; i++)
			printf("%f\n", sampleKeyArray_all[i]);	
*/
	/* procs-1 numbers are enough to determine which elements belong to which bucket/proc. */
	splitterArray = (keyType*) malloc((procs-1) * sizeof(keyType));

	/* Sorting the collected samples, and sending back splitters */
	if(myid==0)
	{
		qsort( sampleKeyArray_all, procs*n_samples, sizeof(keyType), compare_keyType );
/*
		for(i=0; i<procs*n_samples; i++)
			printf("%f\n", sampleKeyArray_all[i]);	
*/
		//printf("Splitter array is: ");
		for(i=0; i<procs-1; i++)
		{
			splitterArray[i] = sampleKeyArray_all[ (i+1) * n_samples - 1 ];
         //printf("%f ", splitterArray[i]);
		}
      //printf("\n");

	}

	MPI_Bcast(splitterArray, (procs-1) * sizeof(keyType), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* find the offset index for each send using binary search on splitter array */
	send_index = (int *) calloc(procs, sizeof(int));
	send_index[0] = 0;

	for (i=1; i<procs; i++)
	{
		//printf("%d kmin = %d rkmin = %g kmax = %d rkmax = %g\n", myid, 0, getKey(&buf[0]), (*local_N)-1, getKey(&buf[(*local_N)-1]) );
		send_index[i] = binary_search(buf, splitterArray[i-1], 0, (*local_N)-1);
	}

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

	/***** Binary info *****/
	int* b_send_index = (int*) malloc(procs * sizeof(int));
	int* b_send_count = (int*) malloc(procs * sizeof(int));
	int* b_recv_count = (int*) malloc(procs * sizeof(int));
	binary_t* buf_bin = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));
	binary_t *b_recvBuf;
	MPI_Status *b_sendStatus, *b_recvStatus;
	MPI_Request	*b_sendReq, *b_recvReq;
	int b_DataSize;
	MPI_Type_size(b_dataType, &b_DataSize);

	int j, k=0;
	for(i=0; i<procs;i++)
	{
		b_send_index[i] = k;
		for(j=send_index[i]; j<send_index[i]+send_count[i]; j++)
		{
			if(buf[j].binind > 0/* && binary[buf[j].binind].inuse*/)
			{
				memcpy(&buf_bin[k], &binary[buf[j].binind], b_DataSize);
				k++;
			}
		}
		b_send_count[i] = k - b_send_index[i];
	}

	//Set binary array to zeros, if not might cause problems when new stars are created in the next timestep. So it's best to wipe out the older data.
	memset (binary, 0, N_BIN_DIM_OPT * sizeof(binary_t));

	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);

	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];
	/***** End binary info *****/


	actual_count = (int*) malloc(procs * sizeof(int));
	MPI_Allgather( &total_recv_count, 1, MPI_INT, actual_count, 1, MPI_INT, commgroup );
	timeEndSimple(tmpTimeStart2, &t_sort2);
	tmpTimeStart2 = timeStartSimple();

	//Finding the maximum size of resultBuf to be allocated.
	max_alloc_outbuf_size = 0;
	for(i=0; i<procs; i++)
		max_alloc_outbuf_size = (actual_count[i] > max_alloc_outbuf_size) ? actual_count[i] : max_alloc_outbuf_size;

	/* allocate recv buffer */
	//MPI3: May be one could just use actual_count[myid] instead of max_alloc...?
	resultBuf = (type*) malloc(max_alloc_outbuf_size * sizeof(type)); 

/*
	for(i=0; i<procs; i++)
		printf("proc %d:\tsend_count[%d] = %d\tsend_index[%d] = %d\trecv_count[%d] = %d \n", myid, i, b_send_count[i], i, b_send_index[i], i, b_recv_count[i]);
*/

	/* allocate asynchronous I/O request and wait status */
	sendStatus = (MPI_Status*)  malloc(2*procs* sizeof(MPI_Status));
	sendReq    = (MPI_Request*) malloc(2*procs* sizeof(MPI_Request));
	recvStatus = sendStatus + procs;
	recvReq    = sendReq    + procs;
	b_sendStatus = (MPI_Status*)  malloc(2*procs* sizeof(MPI_Status));
	b_sendReq    = (MPI_Request*) malloc(2*procs* sizeof(MPI_Request));
	b_recvStatus = b_sendStatus + procs;
	b_recvReq    = b_sendReq    + procs;
	binary_t* b_resultBuf = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));

	/* asynchronous sends */
	for (i=0; i<procs; i++) {
		sendReq[i] = MPI_REQUEST_NULL;
		if (myid != i && send_count[i] > 0)
			MPI_Isend(buf + send_index[i],// *dataSize,
					send_count[i],
					dataType,
					i,
					0,
					commgroup,
					&sendReq[i]);

		b_sendReq[i] = MPI_REQUEST_NULL;
		if (myid != i && b_send_count[i] > 0)
			MPI_Isend(buf_bin + b_send_index[i],// *dataSize,
					b_send_count[i],
					b_dataType,
					i,
					0,
					commgroup,
					&b_sendReq[i]);
	}

	/* asynchronous recv */
	recvBuf = resultBuf;
	b_recvBuf = b_resultBuf+1;
	for (i=0; i<procs; i++) {
		recvReq[i] = MPI_REQUEST_NULL;
		b_recvReq[i] = MPI_REQUEST_NULL;
		if (myid == i) {
			memcpy(recvBuf, buf+send_index[i], send_count[i]*dataSize);
			memcpy(b_recvBuf, buf_bin+b_send_index[i], b_send_count[i]*b_DataSize);
		} else {
			if (recv_count[i] > 0) {
				MPI_Irecv(recvBuf,
						recv_count[i],
						dataType,
						i,
						0,
						commgroup,
						&recvReq[i]);
			}
			if (b_recv_count[i] > 0) {
				MPI_Irecv(b_recvBuf,
						b_recv_count[i],
						b_dataType,
						i,
						0,
						commgroup,
						&b_recvReq[i]);
			}

		}
		recvBuf += recv_count[i];
		b_recvBuf += b_recv_count[i];
	}

	/* wait till all asynchronous i/o done */
	MPI_Waitall(2*procs, sendReq, sendStatus);
	MPI_Waitall(2*procs, b_sendReq, b_sendStatus);

	dprintf("new no.of stars in proc %d = %d\n", myid, total_recv_count);

	//MPI3: Before we do local sort, we need to fix the binary addressing.
	//Add additional checks!! Test thoroughly. Before alltoall use the id value of star elements to one of the ids of binaries, and check post sort - for testing.
	k=1; 
	int kprev;
//	for(i=0;i<total_recv_count;i++)
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
		if(k-kprev!=b_recv_count[i]) eprintf("mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-kprev, b_recv_count[i]);
	}

	if(k-1!=b_total_recv_count)
		eprintf("Binary numbers mismatch in proc %d j = %d recv_cnt = %d\n", myid, k-1, b_total_recv_count);

	timeEndSimple(tmpTimeStart2, &t_sort3);
	tmpTimeStart2 = timeStartSimple();
	qsort(resultBuf, total_recv_count, sizeof(type), compare_type);
	timeEndSimple(tmpTimeStart2, &t_sort4);
	timeEndSimple(tmpTimeStart, &t_sort_only);
/*
#ifdef USE_MPI
		strcpy(filename, "test_out_par");
		strcpy(tempstr, filename);
		sprintf(num, "%d", myid);
		strcat(tempstr, num);
		strcat(tempstr, ".dat");
		for( i = 0; i < procs; i++ )
		{
			if(myid == i)
			{
				//printf("id=%d\tbegin=%d\tend=\%d\n", myid, mpiBegin, mpiEnd);
				ftest = fopen( tempstr, "w" );
				for( j = 0; j <= total_recv_count; j++ )
				{
					//if(star[j].binind>0)
						//fprintf(ftest, "%ld\t%.18g\n", j, binary[star[j].binind].a);
					if(resultBuf[j].binind>0)
						fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\t%ld\n", get_global_idx(j+1), resultBuf[j].r, resultBuf[j].id, b_resultBuf[resultBuf[j].binind].id1, b_resultBuf[resultBuf[j].binind].id2, resultBuf[j].binind);
					else
						fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\n", get_global_idx(j+1), resultBuf[j].r, resultBuf[j].id, resultBuf[j].binind);
				}
				fclose(ftest);
			}
			MPI_Barrier(commgroup);
		}
		if(myid==0)
		{
			char process_str[30];
			sprintf(process_str, "./process.sh %d", procs);
			system(process_str);
		}
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
		{
			if(star[i].binind>0)
				fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\n", i, star[i].r, star[i].id, binary[star[i].binind].id1, binary[star[i].binind].id2);
			else
				fprintf(ftest, "%ld\t%.18g\t%ld\n", i, star[i].r, star[i].id);
		}
		fclose(ftest);
#endif
*/
	tmpTimeStart = timeStartSimple();
	load_balance(resultBuf, buf, b_resultBuf, binary, expected_count, actual_count, myid, procs, dataType, b_dataType, commgroup);
	timeEndSimple(tmpTimeStart, &t_load_bal);

	MPI_Barrier(commgroup);

	tmpTimeStart = timeStartSimple();
	*local_N = actual_count[myid];
	//N_b_local=b_total_recv_count;
/*
#ifdef USE_MPI
		strcpy(filename, "test_out_par");
		strcpy(tempstr, filename);
		sprintf(num, "%d", myid);
		strcat(tempstr, num);
		strcat(tempstr, ".dat");
		for( i = 0; i < procs; i++ )
		{
			if(myid == i)
			{
				//printf("id=%d\tbegin=%d\tend=\%d\n", myid, mpiBegin, mpiEnd);
				ftest = fopen( tempstr, "w" );
				for( j = 0; j < actual_count[i]; j++ )
				{
					//if(star[j].binind>0)
						//fprintf(ftest, "%ld\t%.18g\n", j, binary[star[j].binind].a);
					if(buf[j].binind>0)
						fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\t%ld\n", get_global_idx(j+1), buf[j].r, buf[j].id, binary[buf[j].binind].id1, binary[buf[j].binind].id2, buf[j].binind);
					else
						fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\n", get_global_idx(j+1), buf[j].r, buf[j].id, buf[j].binind);
				}
				fclose(ftest);
			}
			MPI_Barrier(commgroup);
		}
		if(myid==0)
		{
			char process_str[30];
			sprintf(process_str, "./process.sh %d", procs);
			system(process_str);
		}
#else
		strcpy(tempstr, "test_out_ser.dat");
		ftest = fopen( tempstr, "w" );
		for( i = 1; i <= clus.N_MAX; i++ )
		{
			if(star[i].binind>0)
				fprintf(ftest, "%ld\t%.18g\t%ld\t%ld\t%ld\n", i, star[i].r, star[i].id, binary[star[i].binind].id1, binary[star[i].binind].id2);
			else
				fprintf(ftest, "%ld\t%.18g\t%ld\n", i, star[i].r, star[i].id);
		}
		fclose(ftest);
#endif
*/

	free(expected_count);
	free(send_index);
	free(send_count);
	free(recv_count);
	free(b_send_index);
	free(b_send_count);
	free(b_recv_count);
	free(buf_bin);
	free(resultBuf);
	free(sendStatus);
	free(sendReq);
	free(b_resultBuf);
	free(b_sendStatus);
	free(b_sendReq);
	free(actual_count);
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

	int 	i, local_count, dataSize, total_recv_count;
	int 	*expected_cum_count, *actual_cum_count, *splitter;
	int  	*send_index;
	int  	*send_count, *recv_count;
	type 	*recvBuf;

	local_count = actual_count[myid];

	//TESTING
/*
	int global_N, sum=0;
	MPI_Allreduce(&local_count, &global_N, 1, MPI_INT, MPI_SUM, commgroup);

	for(i=0; i<procs-1; i++)
	{
		expected_count[i] = rand() % (global_N / procs) + 1;
		sum += expected_count[i];

	}
	expected_count[procs-1] = global_N - sum;
*/
	//END TESTING

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
		//printf("id = %d, expected_cum_count[%d] = %d\n", myid, i, expected_cum_count[i]);
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
			//send_count[i] = (send_index[i] + expected_count[i] > local_count) ? expected_count[i] : local_count;
		}

		if (splitter[i] >= actual_cum_count[myid] + local_count)
		{
			break;
			//send_index[i] = total_count - 1;
			//send_count[i-1] = local_count - send_index[i-1];
			//if(myid==0)
			//	printf("id = %di i = %d s_i = %d s_c = %d l_c = %d\n", myid, i, send_index[i-1], send_count[i-1], local_count);
			//if( (send_index[i-1] + send_count[i-1]) == local_count ) break;	
		}
	}

	send_count[i-1] = local_count - send_index[i-1];

	//for (i=0; i<procs-1; i++)
	//	send_count[i] = send_index[i+1] - send_index[i];

	//send_count[procs-1] = local_count - (send_index[procs-1] + send_count[proc-1]); // + 1;
	MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);

	total_recv_count = 0;
	for (i=0; i<procs; i++) total_recv_count += recv_count[i];
	dprintf("after load-balancing stars in proc %d = %d expected_count = %d\n", myid, total_recv_count, expected_count[myid]);

	int* recv_displ = (int*) malloc(procs * sizeof(int));
	recv_displ[0] = 0;
	for(i=1; i<procs; i++)
		recv_displ[i] = recv_displ[i-1] + recv_count[i-1];

	/***** Binary info *****/
	int* b_send_index = (int*) malloc(procs * sizeof(int));
	int* b_send_count = (int*) malloc(procs * sizeof(int));
	int* b_recv_count = (int*) malloc(procs * sizeof(int));
	binary_t* buf_bin = (binary_t*) malloc(N_BIN_DIM_OPT * sizeof(binary_t));
	binary_t *b_recvBuf;
	MPI_Status *b_sendStatus, *b_recvStatus;
	MPI_Request	*b_sendReq, *b_recvReq;

	int j, k=0;
	for(i=0; i<procs;i++)
	{
		b_send_index[i] = k;
		for(j=send_index[i]; j<send_index[i]+send_count[i]; j++)
			if(inbuf[j].binind > 0 && b_inbuf[inbuf[j].binind].inuse)
			{
				memcpy(&buf_bin[k], &b_inbuf[inbuf[j].binind], sizeof(binary_t));
				k++;
			}
		b_send_count[i] = k - b_send_index[i];
	}
	MPI_Alltoall(b_send_count, 1, MPI_INT, b_recv_count, 1, MPI_INT, commgroup);

	int b_total_recv_count = 0;
	for (i=0; i<procs; i++) b_total_recv_count += b_recv_count[i];
	N_b_local = b_total_recv_count;
	/***** End binary info *****/

	//for(i=0; i<procs; i++)
	//	printf("proc %d:\tsend_index[%d] = %d\tsend_count[%d] = %d\trecv_count[%d] = %d \n", myid, i, send_index[i], i, send_count[i], i, recv_count[i]);

	/* allocate asynchronous I/O request and wait status */
	MPI_Status 	*sendStatus = (MPI_Status*)  malloc(2*procs* sizeof(MPI_Status));
	MPI_Request *sendReq    = (MPI_Request*) malloc(2*procs* sizeof(MPI_Request));
	MPI_Status 	*recvStatus = sendStatus + procs;
	MPI_Request *recvReq    = sendReq    + procs;
	MPI_Type_size(dataType, &dataSize);
	b_sendStatus = (MPI_Status*)  malloc(2*procs* sizeof(MPI_Status));
	b_sendReq    = (MPI_Request*) malloc(2*procs* sizeof(MPI_Request));
	b_recvStatus = b_sendStatus + procs;
	b_recvReq    = b_sendReq    + procs;
	int b_DataSize;
	MPI_Type_size(b_dataType, &b_DataSize);

	/* asynchronous sends */
	for (i=0; i<procs; i++) {
		sendReq[i] = MPI_REQUEST_NULL;
		if (myid != i && send_count[i] > 0)
			MPI_Isend(inbuf + send_index[i],// *dataSize,
					send_count[i],
					dataType,
					i,
					0,
					commgroup,
					&sendReq[i]);

		b_sendReq[i] = MPI_REQUEST_NULL;
		if (myid != i && b_send_count[i] > 0)
			MPI_Isend(buf_bin + b_send_index[i],// *dataSize,
					b_send_count[i],
					b_dataType,
					i,
					0,
					commgroup,
					&b_sendReq[i]);
	}

	/* asynchronous recv */
	recvBuf = outbuf;
	b_recvBuf = b_outbuf+1;
	for (i=0; i<procs; i++) {
		recvReq[i] = MPI_REQUEST_NULL;
		b_recvReq[i] = MPI_REQUEST_NULL;
		if (myid == i) {
			memcpy(recvBuf, inbuf+send_index[i], send_count[i]*dataSize);
			memcpy(b_recvBuf, buf_bin+b_send_index[i], b_send_count[i]*b_DataSize);
		} else {
			if (recv_count[i] > 0) {
				MPI_Irecv(recvBuf,
						recv_count[i],
						dataType,
						i,
						0,
						commgroup,
						&recvReq[i]);
			}
			if (b_recv_count[i] > 0) {
				MPI_Irecv(b_recvBuf,
						b_recv_count[i],
						b_dataType,
						i,
						0,
						commgroup,
						&b_recvReq[i]);
			}

		}
		recvBuf += recv_count[i];
		b_recvBuf += b_recv_count[i];
	}
	/* wait till all asynchronous i/o done */
	MPI_Waitall(2*procs, sendReq, sendStatus);
	MPI_Waitall(2*procs, b_sendReq, b_sendStatus);


	k=1; 
	int kprev;
//	for(i=0;i<total_recv_count;i++)
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

	actual_count[myid] = total_recv_count;

	free(splitter);
	free(send_index);
	free(send_count);
	free(recv_count);
	free(sendStatus);
	free(sendReq);
	free(b_send_index);
	free(b_send_count);
	free(b_recv_count);
	free(b_sendStatus);
	free(b_sendReq);
	free(buf_bin);
	free(recv_displ);
	free(actual_cum_count);
	free(expected_cum_count);
}

#endif
/*
	int j;
	strcpy(filename, "test_out_par");
	strcpy(tempstr, filename);
	sprintf(num, "%d", myid);
	strcat(tempstr, num);
	strcat(tempstr, ".dat");
	ftest = fopen( tempstr, "w" );
	for( j = 1; j <= total_recv_count; j++ )
		fprintf(ftest, "%.18g\n", starData[j].r);
	fclose(ftest);
	if(myid==0)
		system("./process.sh");
*/
/*
      strcpy(filename, "test_keys");
      strcpy(tempstr, filename);
      sprintf(num, "%d", myid);
      strcat(tempstr, num);
      strcat(tempstr, ".dat");
      for( i = 0; i < procs; i++ )
      {
         if(myid == i)
         {
            ftest = fopen( tempstr, "w" );
				for( i = 0; i < num_samples; i++ )
					fprintf(ftest, "%d\t%.18g\t\n", i, sampleKeyArray_local[i]);
            fclose(ftest);
         }
      }

		 strcpy(tempstr, "test_presort.dat");
		 ftest = fopen( tempstr, "w" );
		 for( i = 0; i < procs*num_samples; i++ )
			 fprintf(ftest, "%d\t%.18g\t\n", i, sampleKeyArray_all[i]);
		 fclose(ftest);


		 strcpy(tempstr, "test_postsort.dat");
		 ftest = fopen( tempstr, "w" );
		 for( i = 0; i < procs*num_samples; i++ )
			 fprintf(ftest, "%d\t%.18g\t\n", i, sampleKeyArray_all[i]);
		 fclose(ftest);

	 // return resultBuf;

		int j;
		strcpy(filename, "test_out_par");
		strcpy(tempstr, filename);
		sprintf(num, "%d", myid);
		strcat(tempstr, num);
		strcat(tempstr, ".dat");
		for( i = 0; i < procs; i++ )
		{
			if(myid == i)
			{
				//printf("Start[i]=%d\tend=\%d\n", Start[i], End[i]);
				ftest = fopen( tempstr, "w" );
				for( j = 0; j < total_recv_count; j++ )
				//for( j = mpiBegin; j <= mpiEnd; j++ )
				//for( j = 1; j <= clus.N_MAX; j++ )
					fprintf(ftest, "%.18g\n", resultBuf[j].r);
					//fprintf(ftest, "%ld\t%ld\t\n", mpiBegin+j-1, star[j].id);
					//fprintf(ftest, "%ld\t%.18g\t\n", j, star_m[j]);
				fclose(ftest);
			}
		}
		if(myid==0)
			system("./process.sh");
MPI_Barrier(MPI_COMM_WORLD);

*/
