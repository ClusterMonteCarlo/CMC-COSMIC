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
double getKey(star_t* starData) { return (starData->r); }

int compare (const void * a, const void * b) { return ( *(double*)a > *(double*)b ) ? 1 : -1; }

int binary_search( star_t* starData, double r, int kmin, int kmax )
{
	long ktry;

	if (starData[kmin].r > r  || starData[kmax].r < r) {
		eprintf("r is outside kmin kmax!!\n");
		return -1;
	};      

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
int sample_sort( star_t        *starData,
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
	int			total_recv_num;
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

	dprintf("proc %d: local_N = %ld\n", myid, *local_N);

	MPI_Type_size(eType, &eSize);

	/* local in-place sort */
	qsorts(starData+1, *local_N);

	// Fixing local_N to account for lost stars.
	i = *local_N;
	while(star[i].r == SF_INFINITY)
	{
		dprintf("stripped star found and removed.\n");
		zero_star(i); //Check with Stefan if this is needed.
		(*local_N)--;
		i--;
	}

	//find global total number of elements to be sorted in parallel
	MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);

	//mpiFindIndicesCustom( global_N, 20, myid, &mpiBegin, &mpiEnd );
	findLimits( global_N, 20 );

	int *expected_count = (int*) malloc(procs * sizeof(int));
	for(i=0; i<procs; i++)
		expected_count[i] = End[i] - Start[i] + 1;

	//MPI3: If we use time as seed, reproducibility will become a problem. Using the timestep count as random seed instead to have both some randomness and reproducability..
	//srand ( time(NULL) );
	srand ( tcount );

	/* Sample arrays to collect samples from each node and send to root */
	sampleKeyArray_local = (double*) malloc(num_samples * sizeof(double));
	sampleKeyArray_all = (double*) malloc(procs * num_samples * sizeof(double));
	/* procs-1 numbers are enough to determine which elements belong ti which bucket/proc. */
	splitterArray = (double*) malloc((procs-1) * sizeof(double));

	/* Picking random samples and sending to root */
	for(i=0; i<num_samples; i++)
		sampleKeyArray_local[i] = starData[ rand() % *local_N + 1 ].r;

	MPI_Datatype sampleType;
	MPI_Type_contiguous( sizeof(sampleKeyArray_local), MPI_BYTE, &sampleType );
	MPI_Type_commit( &sampleType );

	MPI_Gather(sampleKeyArray_local, num_samples, sampleType, sampleKeyArray_all, num_samples, sampleType, 0, commgroup);

	/* Sorting the collected samples, and sending back splitters */
	if(myid==0)
	{
		qsort( sampleKeyArray_all, procs*num_samples, sizeof(double), compare );

		for(i=0; i<procs-1; i++)
			splitterArray[i] = sampleKeyArray_all[ (i+1) * num_samples - 1 ];
	}

	MPI_Bcast(splitterArray, procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/* find the offset index for each send using binary search on splitter array */
	send_index = (int *) calloc(procs, sizeof(int));
	send_index[0] = 1; //0 is sentinel, so set to 1

	for (i=1; i<procs; i++)
		send_index[i] = binary_search(starData, splitterArray[i-1], 1, *local_N);

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
	total_recv_num = recv_count[0];
	for (i=1; i<procs; i++) total_recv_num += recv_count[i];

	/* End Experimental */
#endif

	/* allocate recv buffer */
	total_recv_num = recv_count[0];
	for (i=1; i<procs; i++) total_recv_num += recv_count[i];

	// Watch out! If the no.of stars after sort is very uneven, there is a possibility of it exceeding the allocated memore. In that case load_bal might have to be increased.
	double load_bal = 2;
	resultBuf = (star_t*) malloc((int)floor((double)expected_count[myid] * load_bal) * sizeof(star_t)); //The multiplying factor load_bal is for collecting data from neighbors fix load imbalance due to parallel sorting.
	//resultBuf = (star_t*) malloc(total_recv_num * sizeof(star_t));
	//dprintf("Allocating memory to receive %d stars on proc %d\n", (int)floor((double)expected_count[myid]*load_bal), myid);
	dprintf("Allocating memory to receive %d stars on proc %d\n", total_recv_num, myid);

	for(i=0; i<procs; i++)
		dprintf("proc %d:\tsend_count[%d] = %d\tsend_index[%d] = %d\trecv_count[%d] = %d \n", myid, i, send_count[i], i, send_index[i], i, recv_count[i]);

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

	dprintf("New no.of stars in proc %d = %d\n", myid, total_recv_num);

	free(send_index);
	free(send_count);
	free(recv_count);
	free(sendStatus);
	free(sendReq);

	/* merge sort might be faster */
	qsorts(resultBuf, total_recv_num);

	/* Post-sort communication for load balancing */
	// Variables with _back refer too backward communication, and fwd to forward communication.
	s_count_back = 0; r_count_back = 0; s_count_fwd = 0; r_count_fwd = 0;
	actual_count = (int*) malloc(procs * sizeof(int));
	actual_cum_count = (int*) malloc(procs * sizeof(int));
	expected_cum_count = (int*) malloc(procs * sizeof(int));

	MPI_Allgather( &total_recv_num, 1, MPI_INT, actual_count, 1, MPI_INT, commgroup );

	int tmp1 = 0, tmp2 = 0;
	for(i=0; i<procs; i++)
	{
		expected_cum_count[i] = tmp1 + expected_count[i];
		actual_cum_count[i] = tmp2 + actual_count[i];
		tmp1 = expected_cum_count[i];
		tmp2 = actual_cum_count[i];
	}

	dprintf("actual\t%d\t%d\t%d\t%d\t%d\t%d\n", myid, total_recv_num, actual_count[0], actual_count[1], actual_count[2], actual_count[3] );
	if(myid==0)
		dprintf("\nideal\t%d\t%d\t%d\t%d\t%d\t%d\n", myid, global_N, expected_count[0], expected_count[1], expected_count[2], expected_count[3]);

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
		MPI_Recv(resultBuf + total_recv_num, r_count_back * sizeof(star_t), MPI_BYTE, myid + 1, 0, commgroup, &stat_back);

	// Forward Communication
	if( s_count_fwd > 0 )
		MPI_Send(resultBuf + total_recv_num - s_count_fwd, s_count_fwd * sizeof(star_t), MPI_BYTE, myid + 1, 0, commgroup);

	if( r_count_fwd > 0 )
		MPI_Recv(star + 1, r_count_fwd * sizeof(star_t), MPI_BYTE, myid - 1, 0, commgroup, &stat_fwd);

	//int count_back, count_fwd;
	//MPI_Get_count(&stat_back, MPI_BYTE, &count_back);
	//MPI_Get_count(&stat_fwd, MPI_BYTE, &count_fwd);
	//dprintf("%d sent_fwd = %d sent_back = %d recd_fwd = %ld recd_back = %ld\n", myid, s_count_fwd, s_count_back, count_back/sizeof(star_t), count_fwd/sizeof(star_t));
	total_recv_num += r_count_fwd + r_count_back - s_count_fwd - s_count_back;

	dprintf("%d sent_fwd = %d sent_back = %d recd_fwd = %d recd_back = %d\n", myid, s_count_fwd, s_count_back, r_count_fwd, r_count_back);
	dprintf("%d total_recv_num = %d\n", myid, total_recv_num);

	//MPI3: Copy everything back to star array - make sure enough memory is allocated for star in each node. Retain sentinel - 0th star.
	for (i=0; i<total_recv_num; i++)
		starData[i+1+r_count_fwd] = resultBuf[i + s_count_back];

	free(resultBuf);
	free(actual_count);
	free(actual_cum_count);
	free(expected_cum_count);

	*local_N = total_recv_num;

/*
	int j;
	strcpy(filename, "test_out_par");
	strcpy(tempstr, filename);
	sprintf(num, "%d", myid);
	strcat(tempstr, num);
	strcat(tempstr, ".dat");
	ftest = fopen( tempstr, "w" );
	for( j = 1; j <= total_recv_num; j++ )
		fprintf(ftest, "%.18g\n", starData[j].r);
	fclose(ftest);
	if(myid==0)
		system("./process.sh");
*/
	return global_N;
}

#endif
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
				for( j = 0; j < total_recv_num; j++ )
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
