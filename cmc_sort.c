/* -*- linux-c -*- */

#include "cmc.h"

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
