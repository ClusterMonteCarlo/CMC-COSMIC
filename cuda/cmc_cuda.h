#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda/cuda.h>

#ifndef _CUDA_MAIN_
#define _CUDA_EXTERN_ extern
#else
#define _CUDA_EXTERN_
#endif

/*=================================================================
  A brief discussion:
  The optimal sizing found by running various kernel sizes was 
  (28,256).  I bet this changes with N though, so soon there should
  be a way to find the optimal on the fly and store it for the duration
  
  All the _CUDA_EXTERN_ are host variables that need to be copied over
  to the __device__ variables that live in cmc_cuda.cu.  They are simply
  made to make the transfer easier.

  As for the __cplusplus, it wouldn't compile without that... not sure
  exactly why yet.

  
====================================================================*/
#define GRID_DIM	28	
#define BLOCK_DIM       256	
#define THREADS		(GRID_DIM*BLOCK_DIM)
#define MAX_LOOPS	30		//based upon the size, N

_CUDA_EXTERN_ double *m;
_CUDA_EXTERN_ double *r;
_CUDA_EXTERN_ double *phi;
_CUDA_EXTERN_ double *cE;
_CUDA_EXTERN_ double *cJ;
_CUDA_EXTERN_ double *cf;

_CUDA_EXTERN_ long  *h_kmin;
_CUDA_EXTERN_ long  *h_kmax;
_CUDA_EXTERN_ long  *h_ktemp;


#ifdef __cplusplus 
extern "C" void cuInitialize();
extern "C" void cuCopyDataToDevice();
extern "C" void cuCalculateKs();
extern "C" void cuCopyDataToHost();
extern "C" void cuCleanUp();
#else
void cuInitialize();
void cuCopyDataToDevice();
void cuCalculateKs();
void cuCopyDataToHost();
void cuCleanUp();
#endif


//====================================================================
// Unforunately cutil is written in cpp so I will duplicate them here
//====================================================================
#ifdef __cplusplus
extern "C" {
#endif
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                                 \
     cudaError err = call;                                                   \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0);

#  define CUDA_SAFE_CALL( call) do {                                         \
    CUDA_SAFE_CALL_NO_SYNC(call);                                            \
    cudaError err = cudaThreadSynchronize();                                 \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0);

#if __DEVICE_EMULATION__

#  define CUT_DEVICE_INIT()

#else
#  define CUT_DEVICE_INIT() do {                                             \
    int deviceCount;                                                         \
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&deviceCount));                \
    if (deviceCount == 0) {                                                  \
        fprintf(stderr, "There is no device.\n");                            \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    int dev;                                                                 \
    for (dev = 0; dev < deviceCount; ++dev) {                                \
        cudaDeviceProp deviceProp;                                           \
        CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&deviceProp, dev));   \
        if (deviceProp.major >= 1 && deviceProp.minor >= 3)                  \
            break;                                                           \
    }                                                                        \
    if (dev == deviceCount) {                                                \
        fprintf(stderr, "There is no device supporting CUDA.\n");            \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    else                                                                     \
        CUDA_SAFE_CALL(cudaSetDevice(dev));                                  \
} while (0);

#endif

#  define CUT_DEVICE_INIT_DRV(cuDevice) do {                                 \
    cuDevice = 0;                                                            \
    int deviceCount = 0;                                                     \
    CUresult err = cuInit(0);                                                \
    if (CUDA_SUCCESS == err)                                                 \
        CU_SAFE_CALL_NO_SYNC(cuDeviceGetCount(&deviceCount));                \
    if (deviceCount == 0) {                                                  \
        fprintf(stderr, "There is no device.\n");                            \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    int dev;                                                                 \
    for (dev = 0; dev < deviceCount; ++dev) {                                \
        int major, minor;                                                    \
        CU_SAFE_CALL_NO_SYNC(cuDeviceComputeCapability(&major, &minor, dev));\
        if (major >= 1)                                                      \
            break;                                                           \
    }                                                                        \
    if (dev == deviceCount) {                                                \
        fprintf(stderr, "There is no device supporting CUDA.\n");            \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    else                                                                     \
        CU_SAFE_CALL_NO_SYNC(cuDeviceGet(&cuDevice, dev));                   \
} while (0);

#ifdef __cplusplus
}
#endif
