
#ifndef _CUDA_MAIN_
#define _CUDA_EXTERN_ extern
#else
#define _CUDA_EXTERN_
#endif

#define GRID_DIM	16		// 16 is good
#define BLOCK_DIM 	128		//128 is good
#define THREADS		(GRID_DIM*BLOCK_DIM)
#define MAX_LOOPS	24		//based upon binary search

_CUDA_EXTERN_ long  *h_kmin;
_CUDA_EXTERN_ long  *h_kmax;
_CUDA_EXTERN_ long  *h_ktemp;

void cuInit();
void cuCopyDataToDevice();
void cuCalculateKs();
void cuCopyDataToHost();
void cuCleanUp();
