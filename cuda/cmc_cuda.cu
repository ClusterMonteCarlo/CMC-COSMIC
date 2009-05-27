#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define _CUDA_MAIN_

#include "cmc_cuda.h" 
#include "../cmc.h"
#include "../cmc_vars.h"

__device__ double *cu_m;
__device__ double *cu_r;
__device__ double *cu_phi;
__device__ double *cu_E;
__device__ double *cu_J;

__device__ long  *cu_kmin;
__device__ long  *cu_kmax;
__device__ long  *cu_ktemp;

long totalStars;

__device__ int sign(double a) {
    return (a < 0 ? -1 : 1);
}


//=========================================================
// Some inline helper functions to keep it clean
//=========================================================
__device__ double cuSQR(double x){
    return x*x;
}

__device__ double cuPHI_S(double rad, long nstar, 
			  double cu_m, double cu_r) //8 flop
{	
    double a = cu_m/nstar;
    if (rad >= cu_r)
	return a/rad;
    return a/cu_r;
}

__device__ double cuFunction_Q(long j, long k, double E, double J, 
			       long nstar, double *cu_m, double *cu_r, double *cu_phi)
{
    double a = cu_r[k];
    return 2.0*(E-(cu_phi[k]+cuPHI_S(a, nstar, cu_m[j],cu_r[j])))-((J/a)*(J/a));
}

#ifndef EXPERIMENTAL
__device__ double cuFUNC(long j, long k, double E, double J, long nstar,  //15 flop
			 double *cu_m, double *cu_r, double *cu_phi)
{
    return 2.0 * cuSQR(cu_r[k]) 
	* (E - cu_phi[k] + cuPHI_S(cu_r[k], nstar, cu_m[j], cu_r[j])) - cuSQR(J);
}
#else
__device__ double cuFUNC(long j, long k, double E, double J, long nstar,  // 21 flop
			 double *cu_m, double *cu_r, double *cu_phi)
{
    return 2.0 * (E - cu_phi[k] + cuPHI_S(cu_r[k], nstar, cu_m[j], cu_r[j])) - cuSQR(J/cu_r[k]);
}
#endif


//========================================================================
// The equivalent of FindZero_Q
//========================================================================
__global__ void cuFindZero_Q(long start, long nstar, double *cu_m, double *cu_r, double *cu_phi, 
			     double *cu_E, double *cu_J, long *cu_kmin, long *cu_kmax, long *cu_ktemp)
{
    
    long n = nstar;
    long si = threadIdx.x + blockDim.x * blockIdx.x + start;
    long j = si;
    
    // should load these values to shared memory? check size of it
    // also, since we dont change E, J can they go into constant memory -> faster?
    // scratch that, no constant memory, maybe texture memory space
    double E, J, t;
    long ktry, kmin, kmax, sa, sb, k1, k2;
    E = cu_E[j];
    J = cu_J[j];
    
    //---------------------------------------------
    // search for kmin first 
    //---------------------------------------------
    kmin = 0;
    kmax = si;//cu_ktemp[si];	
    
    t = cuFUNC(j, kmin, E, J, n, cu_m, cu_r, cu_phi); 
    sa = sign(t);
    t = cuFUNC(j, kmax, E, J, n, cu_m, cu_r, cu_phi);
    sb = sign(t);
    
#pragma unroll 25
    //20 loops should be enough for at least
    for (long i=0; i<=MAX_LOOPS; i++)
	{	
	    ktry = (kmin+kmax+1)/2;
	    t=cuFUNC(j, ktry, E, J, n, cu_m, cu_r, cu_phi);
	    kmin = (1 + sa*sign(t))/2 * (ktry - kmin) + kmin;
	    kmax = (1 + sb*sign(t))/2 * (ktry - kmax - 1) + kmax;
	}
    k1 = kmax;
    
    // ---
    // try #2 on calculations
    // 3 (si) + 9 + 9 (t,sa,sb) + (3+9+13)*30 (loops)
    // 9 + 9 + (3+9+13)*30
    // = 1539 * 5e5 = 0.7695 GLOPs
    // 12.8 GFLOPS
    // ---

    //-------------------------------------------------------------------
    // 44 flop for initial t, sa, sb
    // MAX_LOOPS * (6 + 21 + 19 + 1 (i++)) = 1380 flop + 30flop
    // + 4 extra that were elsewhere
    // 
    // total flops for all of FindZero_Q = (44+1380+30+4)*2 = 2916
    // total flops performed for the search space = 2916 * N = 1.458GFLOP
    // total GFLOPS = 1.454 GFLOP / 0.060sec = 24.6GFLOPS
    //-------------------------------------------------------------------


    //---------------------------------------------
    // search for kmax next 
    //---------------------------------------------
    kmin = si;//cu_ktemp[si];
    kmax = nstar;
    
    t=cuFUNC(j, kmin, E, J, n, cu_m, cu_r, cu_phi);
    sa = sign(t);
    t=cuFUNC(j, kmax, E, J, n, cu_m, cu_r,  cu_phi)-(1e-8);
    sb = sign(t);
    
#pragma unroll 25
    //24 loops should be enough for at least 10^7 stars
    for (long i=0; i<=MAX_LOOPS; i++)
	{	
	    ktry = (kmin+kmax+1)/2;
	    t=cuFUNC(j, ktry, E, J, n, cu_m, cu_r, cu_phi);
	    kmin = (1 + sa*sign(t))/2 * (ktry - kmin) + kmin;
	    kmax = (1 + sb*sign(t))/2 * (ktry - kmax - 1) + kmax;
	}
    k2 = kmax;
    
    cu_kmin[j] = k1;
    cu_kmax[j] = k2;
}


//========================================================================
// prepares the ktemps for later
//========================================================================
__global__ void cuFindKtemps(long start, long nstar, double *cu_m, double *cu_r, 
			     double *cu_phi, double *cu_E, double *cu_J, long *cu_ktemp)
{
    long si = threadIdx.x + blockDim.x * blockIdx.x + start;
    
    //double E, J, Qtemp;
    long ktemp = si;
    //E = cu_E[si];
    //J = cu_J[si];
    
    /*    Qtemp = cuFunction_Q(si, ktemp, E, J, nstar, cu_m, cu_r, cu_phi);
    if (Qtemp < 0.0) 
	{
	    ktemp = -1;
	    do {
		ktemp++;
		Qtemp = cuFunction_Q(si, ktemp, E, J, nstar, cu_m, cu_r, cu_phi);
	    } while (Qtemp < 0.0 && ktemp <= nstar);		
	    if (ktemp >= nstar) 
		ktemp = si;
		}*/
    
    cu_ktemp[si] = ktemp;
}

//========================================================================
// Start up CUDA - needs to be called before anything else
//========================================================================
void cuInitialize()
{
    printf("------------------ USING CUDA -------------------\n");
    printf("\tInitializing devices...\n");
    
    CUT_DEVICE_INIT();
   
    printf("\tFreeing memory for host...\n");
    totalStars = clus.N_STAR*2; //there can be no more than twice the number of initial stars... i hope
    
    m = (double*)malloc(sizeof(double)*totalStars);
    r = (double*)malloc(sizeof(double)*totalStars);
    phi = (double*)malloc(sizeof(double)*totalStars);
    cE = (double*)malloc(sizeof(double)*totalStars);
    cJ = (double*)malloc(sizeof(double)*totalStars);
    
    h_kmin = (long*)malloc(sizeof(long)*totalStars);
    h_kmax = (long*)malloc(sizeof(long)*totalStars);
    h_ktemp = (long*)malloc(sizeof(long)*totalStars);
    
    printf("\tInitializing device memory for 2*%i objects (%ikB)...\n", totalStars/2, totalStars*15/8000);
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_m, sizeof(double)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_r, sizeof(double)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_phi, sizeof(double)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_E, sizeof(double)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_J, sizeof(double)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_kmin, sizeof(long)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_kmax, sizeof(long)*totalStars) );
    CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_ktemp, sizeof(long)*totalStars) );
    cudaError_t er = cudaGetLastError();
    printf("\tCUDA's comments regarding init: %s\n", cudaGetErrorString(er));
    printf("-------------------------------------------------\n");
}

//========================================================================
//copies all of the star info to the device at once
//hopefully this isn't a very big factor, who knows
//========================================================================
void cuCopyDataToDevice()
{
    for (long si = 1; si <= clus.N_MAX_NEW; si++)
	{
	    cE[si] = star[si].E + PHI_S(star[si].r, si);
	    cJ[si] = star[si].J;
	    m[si] = star[si].m;
	    r[si] = star[si].r;
	    phi[si] = star[si].phi;
	}
    
    CUDA_SAFE_CALL( cudaMemcpy( cu_m, m, sizeof(double)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_r, r, sizeof(double)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_phi, phi, sizeof(double)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_E, cE, sizeof(double)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_J, cJ, sizeof(double)*totalStars, cudaMemcpyHostToDevice) );
}

//========================================================================
//copies all of the final products back for processing
//========================================================================
void cuCopyDataToHost()
{
    CUDA_SAFE_CALL( cudaMemcpy( h_kmin, cu_kmin, sizeof(long)*totalStars, cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy( h_kmax, cu_kmax, sizeof(long)*totalStars, cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy( h_ktemp, cu_ktemp, sizeof(long)*totalStars, cudaMemcpyDeviceToHost) );
}


//========================================================================
// Host frontend for finding the k-values
//========================================================================
void cuCalculateKs()
{
    printf("------------------ USING CUDA -------------------\n");
    
    //clock_t start = clock();
    cuCopyDataToDevice();
    //clock_t end = clock();
    //printf("\tCUDA memcpy time: %f sec\n", double(end-start)/(CLOCKS_PER_SEC));
    
    for (long start = 1; start <= clus.N_MAX_NEW; start += THREADS)
    	cuFindKtemps<<< GRID_DIM, BLOCK_DIM >>>(start, clus.N_MAX, cu_m, cu_r, 
    						cu_phi, cu_E, cu_J, cu_ktemp);
    cudaThreadSynchronize();
    
    /*    int grids[6] = {8, 16, 20, 24, 28, 32};
    int blocks[7] = {32, 64, 128, 192, 256, 296, 384}; 
    for (int i=0; i<6; i++) {
	for (int j=0; j<7; j++){
	    clock_t begin = clock();
	    for (int x = 0; x<10; x++){
    */
    for (long start = 1; start <= clus.N_MAX_NEW; start += THREADS) //(grids[i]*blocks[j]))
	cuFindZero_Q<<< GRID_DIM, BLOCK_DIM >>>(start, clus.N_MAX+1, cu_m, cu_r, 
						cu_phi, cu_E, cu_J, cu_kmin, 
						cu_kmax, cu_ktemp);
    cudaThreadSynchronize();
		/*		
		
	    }
	    clock_t last = clock();
	    printf("\t%i x %i : %f sec", grids[i], blocks[j], double(last-begin)/(CLOCKS_PER_SEC*10));
	    cudaError er = cudaGetLastError();
	    printf("\t||\tCUDA says: %s\n", cudaGetErrorString(er));
	}
	}*/
    
    cuCopyDataToHost();
    cudaError er = cudaGetLastError();
    printf("\tCUDA says: %s\n", cudaGetErrorString(er));

    printf("-------------------------------------------------\n");
}

//========================================================================
// clean up all the CUDA crud
//========================================================================
void cuCleanUp()
{
    printf("------------------ USING CUDA -------------------\n");
    printf("\tFreeing device memory...\n");
    
    free(r);
    free(m);
    free(phi);
    free(cE);
    free(cJ);
    free(h_kmin);
    free(h_kmax);
    free(h_ktemp);
    
    CUDA_SAFE_CALL( cudaFree(cu_m) );
    CUDA_SAFE_CALL( cudaFree(cu_r) );
    CUDA_SAFE_CALL( cudaFree(cu_phi) );
    CUDA_SAFE_CALL( cudaFree(cu_J) );
    CUDA_SAFE_CALL( cudaFree(cu_E) );
    CUDA_SAFE_CALL( cudaFree(cu_kmin) );
    CUDA_SAFE_CALL( cudaFree(cu_kmax) );
    CUDA_SAFE_CALL( cudaFree(cu_ktemp) );
    
    printf("\tReleasing devices...\n");
    printf("-------------------------------------------------\n");
}

