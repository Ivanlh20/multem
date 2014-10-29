#ifndef hgeneralGPU_H
#define hgeneralGPU_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

__device__ unsigned int stRetCount = 0;
__device__ double stTotSumd = 0;
__device__ double2 stTotSumc;

/*------------------------------------------------------------*/
/* PURPOSE: Atomic addition double values					  */
/*------------------------------------------------------------*/
__device__ inline void atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
__device__ inline double bessi0GPU(double x)
{
	 double ax,ans;
	 double y;

	 if ((ax=fabs(x)) < 3.75){
		 y=x/3.75,y=y*y;
		 ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	 }else{
		 y=3.75/ax;
		 ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
		 +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
		 +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
		 +y*0.392377e-2))))))));
	 }
	 return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1. */
/*------------------------------------------------------------*/
__device__ inline double bessi1GPU(double x)
{
	 double ax,ans;
	 double y;

	 if ((ax=fabs(x)) < 3.75){
		 y=x/3.75,y=y*y;
		 ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
		 +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	 }else{
		 y=3.75/ax;
		 ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
		 -y*0.420059e-2));
		 ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
		 +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		 ans *= (exp(ax)/sqrt(ax));
	 }
	 return x < 0.0 ? -ans : ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.  */
/*------------------------------------------------------------*/
__device__ inline double bessk0GPU(double x)
{
	 double y,ans;

	 if (x <= 2.0){
		 y=x*x/4.0;
		 ans=(-log(x/2.0)*bessi0GPU(x))+(-0.57721566+y*(0.42278420
		 +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
		 +y*(0.10750e-3+y*0.74e-5))))));
	 }else{
		 y=2.0/x;
		 ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
		 +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
		 +y*(-0.251540e-2+y*0.53208e-3))))));
	 }
	 return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
__device__ inline double bessk1GPU(double x)
{
	 double y,ans;

	 if (x <= 2.0){
		 y=x*x/4.0;
		 ans=(log(x/2.0)*bessi1GPU(x))+(1.0/x)*(1.0+y*(0.15443144
		 +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
		 +y*(-0.110404e-2+y*(-0.4686e-4)))))));
	 }else{
		 y=2.0/x;
		 ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
		 +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
		 +y*(0.325614e-2+y*(-0.68245e-3)))))));
	 }
	 return ans;
}

/************************************************************************************/
__device__ inline int unrolledBinarySearch128(double &x0, volatile double *x){
	int i0 = 0, ie = stnR-1;
	int im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //64
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //32
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //16
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //8
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //4
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //2
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //1
	
	return i0;
}

/************************************************************************************/
__device__ inline void reduceBlockDouble128(volatile double *M1, volatile double *M2, volatile double *M3, volatile double *M4, int tid)
{ 
	if (tid < 64){
		M1[tid] += M1[tid + 64];
		M2[tid] += M2[tid + 64];
		M3[tid] += M3[tid + 64];
		M4[tid] += M4[tid + 64];
	}
	__syncthreads();

	if (tid < 32){
		M1[tid] += M1[tid + 32];
		M2[tid] += M2[tid + 32];
		M3[tid] += M3[tid + 32];
		M4[tid] += M4[tid + 32];

		M1[tid] += M1[tid + 16];
		M2[tid] += M2[tid + 16];
		M3[tid] += M3[tid + 16];
		M4[tid] += M4[tid + 16];

		M1[tid] += M1[tid + 8];
		M2[tid] += M2[tid + 8];
		M3[tid] += M3[tid + 8];
		M4[tid] += M4[tid + 8];

		M1[tid] += M1[tid + 4];
		M2[tid] += M2[tid + 4];
		M3[tid] += M3[tid + 4];
		M4[tid] += M4[tid + 4];

		M1[tid] += M1[tid + 2];
		M2[tid] += M2[tid + 2];
		M3[tid] += M3[tid + 2];
		M4[tid] += M4[tid + 2];

		M1[tid] += M1[tid + 1];
		M2[tid] += M2[tid + 1];
		M3[tid] += M3[tid + 1];
		M4[tid] += M4[tid + 1];
	}

	__syncthreads();
}

/************************************************************************************/
// reduction double
__device__ inline void k_reduceBlockDouble256(volatile double *M, int tid)
{
	if (tid < 128)
		M[tid] += M[tid + 128];
	__syncthreads();

	if (tid < 64)
		M[tid] += M[tid + 64];
	__syncthreads();

	if (tid < 32){
		M[tid] += M[tid + 32];
		M[tid] += M[tid + 16];
		M[tid] += M[tid + 8];
		M[tid] += M[tid + 4];
		M[tid] += M[tid + 2];
		M[tid] += M[tid + 1];
	}
}

// reduction double
__device__ inline void k_reduceBlockDouble256(volatile double *M1, volatile double *M2, int tid)
{
	if (tid < 128){
		M1[tid] += M1[tid + 128];
		M2[tid] += M2[tid + 128];
	}
	__syncthreads();

	if (tid < 64){
		M1[tid] += M1[tid + 64];
		M2[tid] += M2[tid + 64];
	}
	__syncthreads();

	if (tid < 32){
		M1[tid] += M1[tid + 32];
		M2[tid] += M2[tid + 32];
		M1[tid] += M1[tid + 16];
		M2[tid] += M2[tid + 16];
		M1[tid] += M1[tid + 8];
		M2[tid] += M2[tid + 8];
		M1[tid] += M1[tid + 4];
		M2[tid] += M2[tid + 4];
		M1[tid] += M1[tid + 2];
		M2[tid] += M2[tid + 2];
		M1[tid] += M1[tid + 1];
		M2[tid] += M2[tid + 1];
	}
}

/************************************************************************************/
// reduction complex
__device__ inline void k_reduceBlockComplex256(volatile double2 *M, int tid)
{
	if (tid < 128){
		M[tid].x += M[tid + 128].x;
		M[tid].y += M[tid + 128].y;
	}
	__syncthreads();

	if (tid < 64){
		M[tid].x += M[tid + 64].x;
		M[tid].y += M[tid + 64].y;
	}
	__syncthreads();

	if (tid < 32){
		M[tid].x += M[tid + 32].x;
		M[tid].y += M[tid + 32].y;
		M[tid].x += M[tid + 16].x;
		M[tid].y += M[tid + 16].y;
		M[tid].x += M[tid + 8].x;
		M[tid].y += M[tid + 8].y;
		M[tid].x += M[tid + 4].x;
		M[tid].y += M[tid + 4].y;
		M[tid].x += M[tid + 2].x;
		M[tid].y += M[tid + 2].y;
		M[tid].x += M[tid + 1].x;
		M[tid].y += M[tid + 1].y;
	}
}

// reduction complex
__device__ inline void k_reduceBlockComplex256(volatile double2 *M1, volatile double2 *M2, int tid)
{
	if (tid < 128){
		M1[tid].x += M1[tid + 128].x;
		M1[tid].y += M1[tid + 128].y;
		M2[tid].x += M2[tid + 128].x;
		M2[tid].y += M2[tid + 128].y;
	}
	__syncthreads();

	if (tid < 64){
		M1[tid].x += M1[tid + 64].x;
		M1[tid].y += M1[tid + 64].y;
		M2[tid].x += M2[tid + 64].x;
		M2[tid].y += M2[tid + 64].y;
	}
	__syncthreads();

	if (tid < 32){
		M1[tid].x += M1[tid + 32].x;
		M1[tid].y += M1[tid + 32].y;
		M2[tid].x += M2[tid + 32].x;
		M2[tid].y += M2[tid + 32].y;
		M1[tid].x += M1[tid + 16].x;
		M1[tid].y += M1[tid + 16].y;
		M2[tid].x += M2[tid + 16].x;
		M2[tid].y += M2[tid + 16].y;
		M1[tid].x += M1[tid + 8].x;
		M1[tid].y += M1[tid + 8].y;
		M2[tid].x += M2[tid + 8].x;
		M2[tid].y += M2[tid + 8].y;
		M1[tid].x += M1[tid + 4].x;
		M1[tid].y += M1[tid + 4].y;
		M2[tid].x += M2[tid + 4].x;
		M2[tid].y += M2[tid + 4].y;
		M1[tid].x += M1[tid + 2].x;
		M1[tid].y += M1[tid + 2].y;
		M2[tid].x += M2[tid + 2].x;
		M2[tid].y += M2[tid + 2].y;
		M1[tid].x += M1[tid + 1].x;
		M1[tid].y += M1[tid + 1].y;
		M2[tid].x += M2[tid + 1].x;
		M2[tid].y += M2[tid + 1].y;
	}
}

/************************************************************************************/

void f_ScaleVectorD(dim3 grid, dim3 threads, int n, double f, double *&Vg);

void f_ScaleVectorC(dim3 grid, dim3 threads, int n, double f, double2 *&Vg);

void f_SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *&Vg1, double *&Vg2);

void f_SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *&Vg);

void f_SetValueVectorC(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *&Vg);

void f_SetValueVectorCD(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *&Vcg, double v, double *&Vdg);

void f_GetRealImagVectorC(dim3 grid, dim3 threads, int n, double2 *&Vg, sComplex &Vgri);

void f_SetRealImagVectorC(dim3 grid, dim3 threads, int n, sComplex &Vgri, double2 *&Vg);

void f_AddwDtoD(dim3 grid, dim3 threads, bool add, int n, double w, double *&Vig, double *&Vog);

void f_AddwM2CtoD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double *&Vg2);

void f_AddwCtoC(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vig, double2 *&Vog);

void f_AddwCwM2CtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double2 *&aVg, double *&aVg2);

void f_AddwCwDtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double *&Vg2, double2 *&aVg, double *&aVg2);

void f_fft2ShiftD(dim3 grid, dim3 threads, int nxh, int nyh, double *&Mg);

void f_fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, double2 *&Mg);

void f_fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, sComplex &Mg);

void f_setRetCount(unsigned int ih);

double f_getTotSumd();

// Sum over all the elements
double f_SumMd(dim3 grid, dim3 threads, int nxy, double *&M, double *&Mb);

// Sum module square over all the elements
double f_SumMC2(dim3 grid, dim3 threads, int nxy, double2 *&M, double *&Mb);

// Sum module square over a circular detector
void f_SumMDet(dim3 grid, dim3 threads, sGP &GP, double *&M, double gmin2, double gmax2, double *&Mp, int iSum, double *&SumM);

// Sum module square over a circular detector
void f_SumM1M2Det(dim3 grid, dim3 threads, sGP &GP, double *&M1, double *&M2, double gmin2, double gmax2, double *&M1p, double *&M2p, int iSum, double *&SumM1, double *&SumM2);

/****************************************************************************/
/****************************************************************************/
void f_ReadQuadratureGPU(int typ, int nQGPU, sQ1 &QGPU);

void f_ReadQuadratureCPUGPU(int typ, int nQ, sQ1 &QCPU, sQ1 &QGPU);

/****************************************************************************/
/****************************************************************************/
void f_sCoefPar_cudaFree(sCoefPar &CoefPar);

void f_sCoefPar_cudaInit(sCoefPar &CoefPar);

void f_sCoefPar_cudaMalloc(int nCoefPar, sCoefPar &CoefPar);

/****************************************************************************/
/****************************************************************************/

void f_sciVn_cudaFree(sciVn &ciVn);

void f_sciVn_cudaInit(sciVn &ciVn);

void f_sciVn_cudaMalloc(int nciVn, sciVn &ciVn);

/****************************************************************************/
/****************************************************************************/

void f_sDetCir_cudaFree(sDetCir &DetCir);

void f_sDetCir_cudaInit(sDetCir &DetCir);

void f_sDetCir_cudaMalloc(int nDetCir, sDetCir &DetCir);

/****************************************************************************/
/****************************************************************************/

void Cd_2_Ch(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, sComplex &Psih);

void CdDd_2_ChDh(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih);

void CdDdDd_2_ChDhDh(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h);

/****************************************************************************/
/****************************************************************************/

// Phase multiplication
void f_PhaseMul(dim3 grid, dim3 threads, sGP &GP, double gxu, double gyu, sERg &ERg, double2 *&Psi);

// AntiAliasing
void f_BandwidthLimited2D(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double2 *Psi);

// Transmission function WPO
void f_TransmissionWPO(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double2 *&Trans);

// Transmission function
void f_Transmission1(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double2 *&Trans);

// Transmission function
void f_Transmission2(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double *&V1, double *&V2, eSlicePos SlicePos, double2 *&Trans);

// Trasmit
void f_Transmit(dim3 grid, dim3 threads, sGP &GP, double2 *&Trans, double2 *&Psi);

// Propagate
void f_Propagate(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sProp &Prop, double2 *&Psi);

#endif