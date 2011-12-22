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










#define _DEBUG 1

double getKey(star_t* starData) { return (starData->r); }

int compare (const void * a, const void * b) { return ( *(double*)a > *(double*)b ) ? 1 : -1; }

int binary_search( star_t* starData, double r, int kmin, int kmax )
{
   long ktry;
             
   if (starData[kmin].r > r  || starData[kmax].r < r) {
      printf("r is outside kmin kmax!!\n");
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
   
   return kmin;
}

#ifdef USE_MPI
void bucket_sort( star_t         *starData,
                  long          *local_N, //Pass clus.N_MAX_NEW
                  MPI_Datatype  eType,
                  MPI_Comm      commgroup,
                  int           num_buckets,
                  int           num_samples )
//                  double      (*convert2double)(const void*)
{

    int         i, j;
    int         eSize;
    double      interval;
    double      max, min;
    double      global_max, global_min;
    int         global_N;
//    int         myid, numprocs;
    int        *bucket_count, *global_bucket_count;
    int        *send_index;
    int         average;
    int         count;
    double      upper;
    int         *send_count;
    int         *recv_count;
    int          total_recv_num;
    star_t        *resultBuf;
    star_t        *recvBuf;
    MPI_Status  *sendStatus, *recvStatus;
    MPI_Request *sendReq,    *recvReq;

//    MPI_Comm_size(commgroup, &numprocs);
//    MPI_Comm_rank(commgroup, &myid);

//local_N = clus.N_MAX_NEW;

    if (_DEBUG) printf("local_N = %ld\n",  *local_N);
    /* check num_buckets, should not be smaller than numprocs */
    if (num_buckets < procs)
        num_buckets = procs;

    MPI_Type_size(eType, &eSize);

	 //MPI3: This needs to be fixed to work for our case
	 /* local in-place sort */
    //qsort(buf, *local_N, eSize, compare);

//MPI3: N_MAX_NEW is the current number of stars. So this should work.
qsorts(starData+1, *local_N);

    /* find local and global max, min */
    min = getKey(starData+1);
	//MPI3: No.of stars lost must also be included while passing local_N.
	if( myid == 1 ) (*local_N)--;
    max = getKey(starData+*local_N);

    if (_DEBUG) printf("myid = %d max = %f min = %f\n", myid, max, min);
    MPI_Allreduce(&min, &global_min, 1, MPI_DOUBLE, MPI_MIN, commgroup);
    MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, commgroup);
    if (_DEBUG) printf("myid = %d max = %f min = %f\n", myid, global_max, global_min);

    /* range of each bucket */
/*
    interval = (global_max - global_min) / num_buckets;
    if (_DEBUG) printf("interval = %f \n", interval);
*/

    /* find local bucket count */
/*
    bucket_count = (int *) calloc(num_buckets, sizeof(int));
    j = 0;
    upper = global_min + interval;
    for (i=1; i<*local_N; i++) {
        if (getKey(starData+i) <= upper)
            bucket_count[j]++;
        else {
            while (getKey(starData+i) > upper) {
                upper += interval;
                j++;
            }
            bucket_count[j]++;
        }
    }
*/

    /* find global accumulated buchet numbers */
/*
    global_bucket_count = (int *) calloc(num_buckets, sizeof(int));
    MPI_Allreduce(bucket_count, global_bucket_count, num_buckets,
                  MPI_INT, MPI_SUM, commgroup);
*/

    /* find global total number of elements to be sorted in parallel */
    MPI_Allreduce(local_N, &global_N, 1, MPI_INT, MPI_SUM, commgroup);
    if (_DEBUG) printf("global_N = %d \n", global_N);

    /* find approximate equal numbers of sorted sub-list for each proc */
/*
    average = global_N/procs;
    if (global_N % procs > 0)
        average++;
*/

	//MPI3: If we use time as seed, reproducability will become a problem. Using some random seed.
	 //srand ( time(NULL) );
	 srand ( 23143 );

	 double* sampleKeyArray_local = (double*) malloc(num_samples * sizeof(double));
	 double* sampleKeyArray_all = (double*) malloc(procs * num_samples * sizeof(double));
	//MPI3: procs-1 numbers are enough to determine which elements belong ti which bucket/proc.
	 double* splitterArray = (double*) malloc((procs-1) * sizeof(double));

	 for(i=0; i<num_samples; i++)
		 sampleKeyArray_local[i] = starData[ rand() % *local_N + 1 ].r;

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

	 MPI_Datatype sampleType;
	 MPI_Type_contiguous( sizeof(sampleKeyArray_local), MPI_BYTE, &sampleType );
	 MPI_Type_commit( &sampleType );

	 MPI_Gather(sampleKeyArray_local, num_samples, sampleType, sampleKeyArray_all, num_samples, sampleType, 0, commgroup);

	 if(myid==0)
	 {
		 strcpy(tempstr, "test_presort.dat");
		 ftest = fopen( tempstr, "w" );
		 for( i = 0; i < procs*num_samples; i++ )
			 fprintf(ftest, "%d\t%.18g\t\n", i, sampleKeyArray_all[i]);
		 fclose(ftest);

		 qsort( sampleKeyArray_all, procs*num_samples, sizeof(double), compare );


		 strcpy(tempstr, "test_postsort.dat");
		 ftest = fopen( tempstr, "w" );
		 for( i = 0; i < procs*num_samples; i++ )
			 fprintf(ftest, "%d\t%.18g\t\n", i, sampleKeyArray_all[i]);
		 fclose(ftest);

		 for(i=0; i<procs-1; i++)
			 splitterArray[i] = sampleKeyArray_all[ (i+1) * num_samples - 1 ];
	 }

	 MPI_Bcast(splitterArray, procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    /* find the offset index for each send */
    send_index = (int *) calloc(procs, sizeof(int));
    send_index[0] = 1; //0 is sentinel, so set to 1  /* local array index -> offset for each send */

    for (i=1; i<procs; i++) {
		send_index[i] = binary_search(starData, splitterArray[i-1], 1, *local_N);

/*
        if (average < count + global_bucket_count[i]) {
            send_index[j] += send_index[j-1] + bucket_count[i];
            count = 0;
            j++;
        }
        else {
            count += global_bucket_count[i];
            send_index[j] += bucket_count[i];
        }
*/
    }

/*
    if (_DEBUG) printf("average = %d \n", average);
    free(bucket_count);
    free(global_bucket_count);
*/

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
    send_count[i] = *local_N - send_index[i];
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, commgroup);


    /* allocate recv buffer */
    total_recv_num = recv_count[0];
    for (i=1; i<procs; i++) total_recv_num += recv_count[i];
    resultBuf = (star_t*) malloc(total_recv_num * sizeof(star_t));

	 if (_DEBUG) 
		 for(i=0; i<num_buckets; i++)
		 {
			 printf("%d send_count[%d] = %d\t", myid, i, send_count[i]);
			 printf("%d send_index[%d] = %d\t", myid, i, send_index[i]);
			 printf("%d recv_count[%d] = %d\t", myid, i, recv_count[i]);
		 }
	printf("\n");

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

printf("%d\n\n", total_recv_num);


    free(send_index);
    free(send_count);
    free(recv_count);
    free(sendStatus);
    free(sendReq);

    /* merge sort will be faster */
    //qsort(resultBuf, total_recv_num, eSize, compare);

//MPI3: N_MAX_NEW is the current number of stars. So this should work.
qsorts(resultBuf, total_recv_num);

//MPI3: Copy everything back to star array - make sure enough memory is allocated for star in each node. Retain sentinel - 0th star.
for (i=0; i<total_recv_num; i++)
	star[i+1] = resultBuf[i];

//MPI3: Free result buffer.
free(resultBuf);

    //free(buf);
    *local_N = total_recv_num;
//    return resultBuf;
}
#endif
