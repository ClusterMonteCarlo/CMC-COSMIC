#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cutil.h>

#define _CUDA_MAIN_

#include "cuda.h"
#include "../cmc.h"
#include "../cmc_vars.h"

__device__ float *cu_m;
__device__ float *cu_m_;

__device__ float *cu_r;
__device__ float *cu_r_;

__device__ float *cu_phi;
__device__ float *cu_phi_;

__device__ float *cu_E;
__device__ float *cu_E_;

__device__ float *cu_J;
__device__ float *cu_J_;

__device__ long  *cu_kmin;
__device__ long  *cu_kmax;
__device__ long  *cu_ktemp;

long totalStars;

__device__ int sign(float a)
{
	return (int)(fabsf(a)/a);
}

__device__ int sign(float a, float aa)
{
	if (fabsf(a) > 0.0000001)
		return sign(a);
	return sign(aa);
}


//-------------------------------------------------------
// Found in CUDA examples / yanked from dsfun90 directly
//-------------------------------------------------------
__device__ inline void D2F_subtract(float &p0, float &p1, float a0, 
	float a1, float b0, float b1)
{
	float t1 = a0 - b0;
	float e = t1 - a0;
	float t2 = ((-b0 - e) + (a0 - (t1 - e))) + a1 - b1;

	p0 = e = t1 + t2;
	p1 = t2 - (e - t1);
}

__device__ inline void D2F_multiply(float &p0, float &p1,
	float a0, float a1, float b0, float b1)
{
	float cona = a0 * 8193.0f;
	float conb = b0 * 8193.0f;
	float sa1 = cona - (cona - a0);
	float sb1 = conb - (conb - b0);
	float sa2 = a0 - sa1;
	float sb2 = b0 - sb1;

	float c11 = a0 * b0;
	float c21 = (((sa1 * sb1 - c11) + sa1 * sb2) + sa2 * sb1) + sa2 * sb2;

	float c2 = a0 * b1 + a1 * b0;
	
	float t1 = c11 + c2;
	float e = t1 - c11;
	float t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + a1 * b1;

	p0 = e = t1 + t2;
	p1 = t2 - (e - t1);
}

__device__ inline void D2F_divide(float &p, float &pp, 
	float a, float aa, float b, float bb)
{
	float CONST = 4097.0f; //could be one of (2049, 4097, 8193)
	float s1 = a / b;
	float cona = CONST * s1;
	float conb = CONST * b;
	float a1 = cona - (cona - s1);
	float b1 = conb - (conb - b);
	float a2 = s1 - a1;
	float b2 = b - b1;

	float c11 = s1 * b;
	float c21 = (((a1*b1 - c11) + a1*b2) + a2*b1) + a2*b2;
	float c2 = s1 * bb;

	float t1 = c11 + c2;
	float e = t1 - c11;
	float t2 = ((c2-e)+ (c11-(t1-e))) + c21;

	float t12 = t1 + t2;
	float t22 = t2 - (t12 - t1);

	float t11 = a - t12;
	e = t11 - a;
	float t21 = ((-t12-e)+(a-(t11-e))) + aa - t22;
	float s2 = (t11 + t21) / b;
	p = s1 + s2;
	pp = s2 - (p - s1);
}

__device__ void cuPHI_S(float &a, float &aa, float rad, float rad_,
	long nstar, float cu_m, float cu_m_, float cu_r, float cu_r_)
{	
	float n = (float)nstar/2;
	float nn = 0.0f;//(float)(nstar - n);

	float a0, a1, d0, d1;
	D2F_multiply(d0, d1, cu_r, cu_r_, n, nn);
	D2F_divide(a0, a1, cu_m, cu_m_, d0, d1);

	float b0, b1;
	D2F_multiply(d0, d1, rad, rad_, n, nn);
	D2F_divide(b0, b1, cu_m, cu_m_, d0, d1);

	if (rad>=cu_r){	a = b0; aa = b1; }
	else { a = a0; aa = a1; }
}

__device__ float cuPHI_S(float rad, float rad_,
	long nstar, float cu_m, float cu_m_, float cu_r, float cu_r_)
{	
	float n = (float)nstar;
	float nn = (float)(nstar - n);

	float a0, a1, d0, d1;
	D2F_multiply(d0, d1, cu_r, cu_r_, n, nn);
	D2F_divide(a0, a1, cu_m, cu_m_, d0, d1);

	float b0, b1;
	D2F_multiply(d0, d1, rad, rad_, n, nn);
	D2F_divide(b0, b1, cu_m, cu_m_, d0, d1);

	if (rad>=cu_r)
		return b0;
	return a0;
}

__device__ float cuFunction_Q(long j, long k, float E, float E_, float J, float J_,
	long nstar, float *cu_m, float *cu_m_, float *cu_r, float *cu_r_, float *cu_phi, float *cu_phi_)
{
	float a = cu_r_[k];
	float b = cu_r[k];
	float e = a/b;
	float d = (1 - e + e*e - e*e*e)/b;
	return 2.0*(E-(cu_phi[k]+cuPHI_S(cu_r[k], cu_r_[k], nstar, 
			cu_m[j], cu_m_[j], cu_r[j], cu_r_[j])))-((J*d)*(J*d));
}

__device__ void cuFUNC(float &a0, float &a1, long j, long k, float E, float E_, float J, float J_,
	long nstar, float *cu_m, float *cu_m_, float *cu_r, float *cu_r_, float *cu_phi, float *cu_phi_)
{
	float a, aa, b, bb, l, ll, m, mm, f, ff, r, rr, p, pp, r2, rr2;

	D2F_subtract(a, aa, E, E_, cu_phi[k], cu_phi_[k]);
	cuPHI_S(p, pp, cu_r[k], cu_r_[k], nstar, cu_m[j], cu_m_[j], cu_r[j], cu_r_[j]);

	D2F_subtract(b, bb, a, aa, p, pp);
	D2F_multiply(l, ll, J, J_, J, J_);
	//D2F_multiply(r, rr, cu_r[k], cu_r_[k], cu_r[k], cu_r_[k]);
	r = cu_r[k]*cu_r[k]; rr =0.0f;
	D2F_multiply(r2, rr2, r, rr, 2.0, 0.0);
	D2F_multiply(m, mm, b, bb, r2, rr2);
	D2F_subtract(f, ff, m, mm, l, ll);

	a0 = f;
	a1 = ff;
}

__device__ float cuFUNC(long j, long k, float E, float E_, float J, float J_,
	long nstar, float *cu_m, float *cu_m_, float *cu_r, float *cu_r_, float *cu_phi, float *cu_phi_)
{
	float a, aa, b, bb, l, ll, m, mm, f, ff, r, rr, p, pp, r2, rr2;

	D2F_subtract(a, aa, E, E_, cu_phi[k], cu_phi_[k]);
	cuPHI_S(p, pp, cu_r[k], cu_r_[k], nstar, cu_m[j], cu_m_[j], cu_r[j], cu_r_[j]);

	D2F_subtract(b, bb, a, aa, p, pp);
	D2F_multiply(l, ll, J, J_, J, J_);
	//D2F_multiply(r, rr, cu_r[k], cu_r_[k], cu_r[k], cu_r_[k]);
	r = cu_r[k]*cu_r[k]; rr =0.0f;
	D2F_multiply(r2, rr2, r, rr, 2.0, 0.0);
	D2F_multiply(m, mm, b, bb, r2, rr2);
	D2F_subtract(f, ff, m, mm, l, ll);

	return f;
}

//========================================================================
//start up CUDA - needs to be called before anything else
//========================================================================
void cuInit()
{
	printf("------------------ USING CUDA -------------------\n");
	printf("\tInitializing devices...\n");
	CUT_DEVICE_INIT();

	dfile  = fopen("Qdiff.csv", "w+");
	dfile2 = fopen("kdiff.csv", "w+");
	fprintf(dfile, "Star Index, Device kmin, Host kmin, Device kmax, Host kmax\n");
	fprintf(dfile2, "Star Index, Host kmin, Device kmin, Host kmax, Device kmax\n");

	printf("\tCreating timer...\n");
	cutCreateTimer(&timer);
	totalStars = clus.N_STAR*2; //there can be no more than twice the number of initial stars... i hope

	m = (float*)malloc(sizeof(float)*totalStars);
	m_ = (float*)malloc(sizeof(float)*totalStars);
	r = (float*)malloc(sizeof(float)*totalStars);
	r_ = (float*)malloc(sizeof(float)*totalStars);
    phi = (float*)malloc(sizeof(float)*totalStars);
    phi_ = (float*)malloc(sizeof(float)*totalStars);
	cE = (float*)malloc(sizeof(float)*totalStars);
	cE_ = (float*)malloc(sizeof(float)*totalStars);
	cJ = (float*)malloc(sizeof(float)*totalStars);
	cJ_ = (float*)malloc(sizeof(float)*totalStars);

	h_kmin = (long*)malloc(sizeof(long)*totalStars);
	h_kmax = (long*)malloc(sizeof(long)*totalStars);
	h_ktemp = (long*)malloc(sizeof(long)*totalStars);

	printf("\tInitializing device memory for 2*%i objects (%ikB)...\n", totalStars/2, totalStars*15/8000);
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_m, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_m_, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_r, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_r_, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_phi, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_phi_, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_E, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_E_, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_J, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_J_, sizeof(float)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_kmin, sizeof(long)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_kmax, sizeof(long)*totalStars) );
	CUDA_SAFE_CALL( cudaMalloc( (void**) &cu_ktemp, sizeof(long)*totalStars) );
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
		float E = star[si].E + PHI_S(star[si].r, si);
		cE[si] = (float)(E);
		cE_[si] = (float)(E - cE[si]);

		cJ[si] = (float)star[si].J;
		cJ_[si] = (float)(star[si].J - cJ[si]);

		m[si] = (float)star[si].m;
		m_[si] = (float)(star[si].m - m[si]);

		r[si] = (float)star[si].r;
		r_[si] = (float)(star[si].r - r[si]);

		phi[si] = (float)star[si].phi;
		phi_[si] = (float)(star[si].phi - phi[si]);
	}

	CUDA_SAFE_CALL( cudaMemcpy( cu_m, m, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( cu_m_, m_, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_r, r, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_r_, r, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_phi, phi, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy( cu_phi_, phi_, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( cu_E, cE, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( cu_E_, cE_, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( cu_J, cJ, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( cu_J_, cJ_, sizeof(float)*totalStars, cudaMemcpyHostToDevice) );
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


//TODO - rewrite linear search aspect of calc_orbit_rs and find all ktemps before hand
//if Qtemp is still negative after linear search??? then use negative ktemps and test for that
//later
//========================================================================
// The equivalent of FindZero_Q
//========================================================================
__global__ void cuFindZero_Q(long start, long nstar, float *cu_m, float *cu_m_, 
		float *cu_r, float *cu_r_, float *cu_phi, float *cu_phi_, float *cu_E, 
		float *cu_E_, float *cu_J, float *cu_J_, long *cu_kmin, long *cu_kmax, long *cu_ktemp)
{

	long n = nstar;
	long si = threadIdx.x + blockDim.x * blockIdx.x + start;
	long j = si;

	float E, J, t, tt, E_, J_; 
	long ktry, kmin, kmax, sa, sb, k1, k2;
	E = cu_E[j];  	E_ = cu_E_[j];
	J = cu_J[j];	J_ = cu_J_[j];

/*
//=============================================================
	if(cuFUNC(j, kmin, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)<
		cuFUNC(j, kmax, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)){
		do {
			ktry = (kmin+kmax+1)/2;
			if (cuFUNC(j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)<0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	} else {
		do {
			ktry = (kmin+kmax+1)/2;
			if (cuFUNC(j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)>0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	}

	k1 = kmin;

	kmin = cu_ktemp[si];
	kmax = n;

	if(cuFUNC(j, kmin, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)<
		cuFUNC(j, kmax, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)){
		do {
			ktry = (kmin+kmax+1)/2;
			if (cuFUNC(j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)<0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	} else {
		do {
			ktry = (kmin+kmax+1)/2;
			if (cuFUNC(j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_)>0){
				kmin = ktry;
			} else {
				kmax = ktry-1;
			}
		} while (kmax!=kmin);
	}

	k2 = kmin;
//============================================================
*/
	
	//---------------------------------------------
	// search for kmin first 
	//---------------------------------------------
	kmin = 0;
	kmax = cu_ktemp[si];	

	cuFUNC(t, tt, j, kmin, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
	sa = sign(t, tt);
	cuFUNC(t, tt, j, kmax, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
	sb = sign(t, tt);

	//20 loops should be enough for at least 
	for (int i=0; i<=MAX_LOOPS; i++)
	{	
		ktry = (kmin+kmax+1)/2;
		cuFUNC(t, tt, j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
		kmin = (1 + sa*sign(t, tt))/2 * (ktry - kmin) + kmin;
		kmax = (1 + sb*sign(t, tt))/2 * (ktry - kmax - 1) + kmax;
	}
	k1 = kmax;

	//---------------------------------------------
	// search for kmax next 
	//---------------------------------------------
	kmin = cu_ktemp[si];
	kmax = n;
	
	cuFUNC(t, tt, j, kmin, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
	sa = sign(t, tt);
	cuFUNC(t, tt, j, kmax, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
	sb = sign(t, tt);

	//24 loops should be enough for at least 10^7 stars
	for (int i=0; i<=MAX_LOOPS; i++)
	{	
		ktry = (kmin+kmax+1)/2;
		cuFUNC(t, tt, j, ktry, E, E_, J, J_, n, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
		kmin = (1 + sa*sign(t, tt))/2 * (ktry - kmin) + kmin;
		kmax = (1 + sb*sign(t, tt))/2 * (ktry - kmax - 1) + kmax;
	}
	k2 = kmin;

	cu_kmin[j] = k1;
	cu_kmax[j] = k2;
}


//========================================================================
// prepares the ktemps for later
//========================================================================
__global__ void cuFindKtemps(long start, long nstar, float *cu_m, float *cu_m_, 
		float *cu_r, float *cu_r_, float *cu_phi, float *cu_phi_, float *cu_E, 
		float *cu_E_, float *cu_J, float *cu_J_, long *cu_ktemp)
{
	long si = threadIdx.x + blockDim.x * blockIdx.x + start;

//	float E, J, Qtemp, E_, J_;
	long ktemp = si;
//	E = cu_E[si];  	E_ = cu_E_[si];
//	J = cu_J[si];	J_ = cu_J_[si];
   	
/*	Qtemp = cuFunction_Q(si, ktemp, E, E_, J, J_, nstar, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
	if (Qtemp < 0.0) 
	{
		ktemp = -1;
		do {
			ktemp++;
			Qtemp = cuFunction_Q(si, ktemp, E, E_, J, J_, nstar, cu_m, cu_m_, cu_r, cu_r_, cu_phi, cu_phi_);
		} while (Qtemp < 0.0 && ktemp <= nstar);		
		if (ktemp >= nstar) 
			ktemp = si;
	}*/

	cu_ktemp[si] = ktemp;
}

//========================================================================
// Host frontend for finding the k-values
//========================================================================
void cuCalculateKs()
{
	printf("------------------ USING CUDA -------------------\n");

	cutResetTimer(timer);
	cutStartTimer(timer);
	cuCopyDataToDevice();

	for (long start = 1; start <= clus.N_MAX_NEW; start += THREADS)
		cuFindKtemps<<< GRID_DIM, BLOCK_DIM >>>(start, clus.N_MAX, cu_m, cu_m_, cu_r, 
				cu_r_, cu_phi, cu_phi_, cu_E, cu_E_, cu_J, cu_J_, cu_ktemp);

	cudaThreadSynchronize();

	for (long start = 1; start <= clus.N_MAX_NEW; start += THREADS)
		cuFindZero_Q<<< GRID_DIM, BLOCK_DIM >>>(start, clus.N_MAX, cu_m, cu_m_, cu_r, 
				cu_r_, cu_phi, cu_phi_, cu_E, cu_E_, cu_J, cu_J_, cu_kmin, cu_kmax, cu_ktemp);

	cudaThreadSynchronize();
	cuCopyDataToHost();
	cutStopTimer(timer);
	float time = cutGetTimerValue(timer);
	printf("\tCUDA device time: %f\n", time);
	printf("-------------------------------------------------\n");
}

//========================================================================
// clean up all the CUDA crud
//========================================================================
void cuCleanUp(int argc, char **argv)
{
	printf("------------------ USING CUDA -------------------\n");
	printf("\tFreeing device memory...\n");

    free(r);
	free(r_);
    free(m);
    free(phi);
	free(cE);
	free(cJ);
	free(h_kmin);
	free(h_kmax);
	free(h_ktemp);

    CUDA_SAFE_CALL( cudaFree(cu_m) );
    CUDA_SAFE_CALL( cudaFree(cu_r) );
	CUDA_SAFE_CALL( cudaFree(cu_r_) );
    CUDA_SAFE_CALL( cudaFree(cu_phi) );
    CUDA_SAFE_CALL( cudaFree(cu_J) );
    CUDA_SAFE_CALL( cudaFree(cu_E) );
	CUDA_SAFE_CALL( cudaFree(cu_kmin) );
	CUDA_SAFE_CALL( cudaFree(cu_kmax) );
	CUDA_SAFE_CALL( cudaFree(cu_ktemp) );

	printf("\tReleasing devices...\n");
    CUT_EXIT(0, 0);
	printf("-------------------------------------------------\n");
}

