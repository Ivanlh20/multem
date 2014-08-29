#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hgeneralGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>

/*******************************General kernels********************************/
// Scale Double vector
__global__ void ScaleVectorD(int n, double f, double *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i] *= f;
		i += ii;
	}
}

// Scale Complex vector
__global__ void ScaleVectorC(int n, double f, double2 *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i].x *= f;
		Vg[i].y *= f;
		i += ii;
	}
}

// Set value to 2 Double vector:
__global__ void SetValueVectorD(int n, double v, double *Vg1, double *Vg2)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg1[i] = v;
		Vg2[i] = v;
		i += ii;
	}
}

// Set value to Double vector:
__global__ void SetValueVectorD(int n, double v, double *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i] = v;
		i += ii;
	}
}

// Set value to Complex vector
__global__ void SetValueVectorC(int n, double vr, double vi, double2 *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){	
		Vg[i].x = vr;
		Vg[i].y = vi;
		i += ii;
	}
}

// Set value to Complex vector and Double vector
__global__ void SetValueVectorCD(int n, double vr, double vi, double2 *Vcg, double v, double *Vdg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){	
		Vcg[i].x = vr;
		Vcg[i].y = vi;
		Vdg[i] = v;
		i += ii;
	}
}

// Get Real and Imaginary part of a Complex vector
__global__ void GetRealImagVectorC(int n, double2 *Vg, double *Vrg, double *Vig)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vrg[i] = Vg[i].x;
		Vig[i] = Vg[i].y;
		i += ii;
	}
}

// Set Real and Imaginary part of a Complex vector
__global__ void SetRealImagVectorC(int n, double *Vrg, double *Vig, double2 *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i].x = Vrg[i];
		Vg[i].y = Vig[i];
		i += ii;
	}
}

// Add Mi to Mo and multiply by a weighted factor
template <bool add>
__global__ void AddwDtoD(int n, double w, double *Vig, double *Vog)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;

	while (i < n){
		if (add)
			Vog[i] += w*Vig[i];
		else
			Vog[i] = w*Vig[i];
		i += ii;
	}
}

// Add M*M' to M2 and multiply by a weighted factor
template <bool add>
__global__ void AddwM2CtoD(int n, double w, double2 *Vg, double *Vg2)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double x, y;

	while (i < n){
		x = Vg[i].x; y = Vg[i].y;
		if (add)
			Vg2[i] += w*(x*x + y*y);
		else
			Vg2[i] = w*(x*x + y*y);
		i += ii;
	}
}

// Add Ci to Co and multiply by a weighted factor
template <bool add>
__global__ void AddwCtoC(int n, double w, double2 *Vig, double2 *Vog)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;

	while (i < n){
		if (add){
			Vog[i].x += w*Vig[i].x;
			Vog[i].y += w*Vig[i].y;
		}else{
			Vog[i].x = w*Vig[i].x;
			Vog[i].y = w*Vig[i].y;
		}
		i += ii;
	}
}

// Add Mc to aMc and aMc2
template <bool add>
__global__ void AddwCwM2CtoCD(int n, double w, double2 *Vg, double2 *Vga, double *Vg2a)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double x, y;

	while (i < n){
		x = Vg[i].x; y = Vg[i].y;
		if (add){
			Vga[i].x += w*x; Vga[i].y += w*y;
			Vg2a[i] += w*(x*x + y*y);
		}else{
			Vga[i].x = w*x; Vga[i].y = w*y;
			Vg2a[i] = w*(x*x + y*y);
		}
		i += ii;
	}
}

// Add Mc to aMc and MMca2
template <bool add>
__global__ void AddwCwDtoCD(int n, double w, double2 *Vg, double *Vg2, double2 *Vga, double *Vg2a)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double x, y;

	while (i < n){
		x = Vg[i].x; y = Vg[i].y;
		if (add){
			Vga[i].x += w*x; Vga[i].y += w*y;
			Vg2a[i] += w*Vg2[i];
		}else{
			Vga[i].x = w*x; Vga[i].y += w*y;
			Vg2a[i] = w*Vg2[i];
		}
		i += ii;
	}
}

// Shift Double matrix respect to (nxh, nyh)
__global__ void fft2ShiftD(int nxh, int nyh, double *Mg)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < nxh)&&(j < nyh)){
		int ny = 2*nyh, k1, k2;
		double z;

		k1 = j + i*ny; k2 = (nyh+j)+(nxh+i)*ny;
		z = Mg[k1];
		Mg[k1] = Mg[k2];
		Mg[k2] = z;

		k1 = (nyh+j) + i*ny; k2 = j+(nxh+i)*ny;
		z = Mg[k1];
		Mg[k1] = Mg[k2];
		Mg[k2] = z;
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void fft2ShiftC(int nxh, int nyh, double2 *Mg)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < nxh)&&(j < nyh)){
		int ny = 2*nyh, k1, k2;
		double zr, zi;

		k1 = j + i*ny; k2 = (nyh+j)+(nxh+i)*ny;
		zr = Mg[k1].x; zi = Mg[k1].y;
		Mg[k1].x = Mg[k2].x; Mg[k1].y = Mg[k2].y;
		Mg[k2].x = zr; Mg[k2].y = zi;

		k1 = (nyh+j) + i*ny; k2 = j+(nxh+i)*ny;
		zr = Mg[k1].x; zi = Mg[k1].y;
		Mg[k1].x = Mg[k2].x; Mg[k1].y = Mg[k2].y;
		Mg[k2].x = zr; Mg[k2].y = zi;
	}
}

/****************************************************************************/
/****************************************************************************/

// reduction double
__device__ void reduceBlockDoubleMatrix(volatile double *M, double SumT, const unsigned int tid){

    M[tid] = SumT;
    __syncthreads();

    if (thrnxy >= 512){
        if (tid < 256)
            M[tid] = SumT = SumT + M[tid + 256];
        __syncthreads();
    }

    if (thrnxy >= 256){
        if (tid < 128)
            M[tid] = SumT = SumT + M[tid + 128];
        __syncthreads();
    }

    if (thrnxy >= 128){
        if (tid <  64)
            M[tid] = SumT = SumT + M[tid +  64];
        __syncthreads();
    }

    if (tid < 32){
        if (thrnxy >=  64)
            M[tid] = SumT = SumT + M[tid + 32];

        if (thrnxy >=  32)
            M[tid] = SumT = SumT + M[tid + 16];

        if (thrnxy >=  16)
            M[tid] = SumT = SumT + M[tid +  8];

        if (thrnxy >=   8)
            M[tid] = SumT = SumT + M[tid +  4];

        if (thrnxy >=   4)
            M[tid] = SumT = SumT + M[tid +  2];

        if (thrnxy >=   2)
            M[tid] = SumT = SumT + M[tid +  1];
    }
}

// reduction complex
__device__ void reduceBlockComplexMatrix(volatile double2 *M, double2 SumT, const unsigned int tid){

    M[tid].x = SumT.x; 
    M[tid].y = SumT.y;
    __syncthreads();

    if (thrnxy >= 512){
        if (tid < 256){
            M[tid].x = SumT.x = SumT.x + M[tid + 256].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 256].y;
		}
        __syncthreads();
    }

    if (thrnxy >= 256){
        if (tid < 128){
            M[tid].x = SumT.x = SumT.x + M[tid + 128].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 128].y;
		}

        __syncthreads();
    }

    if (thrnxy >= 128){
        if (tid <  64){
            M[tid].x = SumT.x = SumT.x + M[tid + 64].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 64].y;
		}
        __syncthreads();
    }

    if (tid < 32){
        if (thrnxy >=  64){
            M[tid].x = SumT.x = SumT.x + M[tid + 32].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 32].y;
		}

        if (thrnxy >=  32){
            M[tid].x = SumT.x = SumT.x + M[tid + 16].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 16].y;
		}

        if (thrnxy >=  16){
            M[tid].x = SumT.x = SumT.x + M[tid + 8].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 8].y;
		}

        if (thrnxy >=   8){
            M[tid].x = SumT.x = SumT.x + M[tid + 4].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 4].y;
		}

        if (thrnxy >=   4){
            M[tid].x = SumT.x = SumT.x + M[tid + 2].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 2].y;
		}

        if (thrnxy >=   2)
			{
            M[tid].x = SumT.x = SumT.x + M[tid + 1].x;
            M[tid].y = SumT.y = SumT.y + M[tid + 1].y;
		}
    }
}

// Sum over all the elements
__global__ void SumAd(unsigned int nxy, bool nIsPow2, const double *A, double *TB){    
    __shared__ double M[thrnxy];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(thrnxy*2) + threadIdx.x;
    unsigned int gridSize = thrnxy*2*gridDim.x;
    double SumT = 0;

    while (i < nxy){
        SumT += A[i];
        if (nIsPow2 || i+thrnxy < nxy)
            SumT += A[i+thrnxy];
        i += gridSize;
    }
	
    reduceBlockDoubleMatrix(M, SumT, tid);

    if (tid == 0) TB[blockIdx.x] = M[0];
	
    if (gridDim.x > 1) {
        __shared__ bool amLast;

        __threadfence();

        if (tid==0){
            unsigned int ticket = atomicInc(&stRetCount, gridDim.x);
            amLast = (ticket == gridDim.x-1);
        }
        __syncthreads();

        if (amLast){
            i = tid;
            SumT = 0;

            while (i < gridDim.x){
                SumT += TB[i];
                i += thrnxy;
            }

            reduceBlockDoubleMatrix(M, SumT, tid);

            if (tid==0){
                stTotSumd = M[0];
                stRetCount = 0;
            }
        }
    }
}

// Sum module square over all the elements
__global__ void SumAc2(unsigned int nxy, bool nIsPow2, const double2 *A, double *TB){   
    __shared__ double M[thrnxy];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(thrnxy*2) + threadIdx.x;
    unsigned int gridSize = thrnxy*2*gridDim.x;
    double SumT = 0, x, y;

    while (i < nxy){
		x = A[i].x; y = A[i].y;
        SumT += x*x + y*y;
        if (nIsPow2 || i+thrnxy < nxy){
			x = A[i+thrnxy].x; y = A[i+thrnxy].y;
            SumT += x*x + y*y;
		}
        i += gridSize;
    }
	
    reduceBlockDoubleMatrix(M, SumT, tid);

    if (tid == 0) TB[blockIdx.x] = M[0];
	
    if (gridDim.x > 1) {
        __shared__ bool amLast;

        __threadfence();

        if (tid==0){
            unsigned int ticket = atomicInc(&stRetCount, gridDim.x);
            amLast = (ticket == gridDim.x-1);
        }
        __syncthreads();

        if (amLast){
            i = tid;
            SumT = 0;

            while (i < gridDim.x){
                SumT += TB[i];
                i += thrnxy;
            }

            reduceBlockDoubleMatrix(M, SumT, tid);

            if (tid==0){
                stTotSumd = M[0];
                stRetCount = 0;
            }
        }
    }
}

// Sum module square over all the elements with index
__global__ void SumAc2Ind(int nind, bool nIsPow2, int *ind, const double2 *A, double *TB){
    __shared__ double M[thrnxy];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(thrnxy*2) + threadIdx.x, j;
    unsigned int gridSize = thrnxy*2*gridDim.x;
    double SumT = 0, x, y;

    while (i < nind){
		j = ind[i];
		x = A[j].x; y = A[j].y;
		SumT += x*x + y*y; 
        if (nIsPow2 || i+thrnxy < nind){
			j = ind[i+thrnxy];
			x = A[j].x; y = A[j].y;
			SumT += x*x+y*y; 
		}
        i += gridSize;
    }
	
    reduceBlockDoubleMatrix(M, SumT, tid);

    if (tid == 0) TB[blockIdx.x] = M[0];
	
    if (gridDim.x > 1) {
        __shared__ bool amLast;

        __threadfence();

        if (tid==0){
            unsigned int ticket = atomicInc(&stRetCount, gridDim.x);
            amLast = (ticket == gridDim.x-1);
        }
        __syncthreads();

        if (amLast){
            i = tid;            
            SumT = 0;

            while (i < gridDim.x){
                SumT += TB[i];
                i += thrnxy;
            }

            reduceBlockDoubleMatrix(M, SumT, tid);

            if (tid==0){
                stTotSumd = M[0];
                stRetCount = 0;
            }
        }
    }
}

/**************************************************************************/

void ScaleVectorD(dim3 grid, dim3 threads, int n, double f, double *Vg){
    ScaleVectorD<<<grid, threads>>>(n, f, Vg);
}

void ScaleVectorC(dim3 grid, dim3 threads, int n, double f, double2 *Vg){
    ScaleVectorC<<<grid, threads>>>(n, f, Vg);
}

void SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *Vg1, double *Vg2){
    SetValueVectorD<<<grid, threads>>>(n, v, Vg1, Vg2);
}

void SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *Vg){
    SetValueVectorD<<<grid, threads>>>(n, v, Vg);
}

void SetValueVectorC(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *Vg){
    SetValueVectorC<<<grid, threads>>>(n, vr, vi, Vg);
}

void SetValueVectorCD(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *Vcg, double v, double *Vdg){
    SetValueVectorCD<<<grid, threads>>>(n, vr, vi, Vcg, v, Vdg);
}

void GetRealImagVectorC(dim3 grid, dim3 threads, int n, double2 *Vg, double *Vrg, double *Vig){
    GetRealImagVectorC<<<grid, threads>>>(n, Vg, Vrg, Vig);
}

void SetRealImagVectorC(dim3 grid, dim3 threads, int n, double *Vrg, double *Vig, double2 *Vg){
    SetRealImagVectorC<<<grid, threads>>>(n, Vrg, Vig, Vg);
}

void AddwDtoD(dim3 grid, dim3 threads, bool add, int n, double w, double *Vig, double *Vog){
    if(add)
		AddwDtoD<true><<<grid, threads>>>(n, w, Vig, Vog);
	else
		AddwDtoD<false><<<grid, threads>>>(n, w, Vig, Vog);
}

void AddwM2CtoD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double *Vg2){
    if(add)
		AddwM2CtoD<true><<<grid, threads>>>(n, w, Vg, Vg2);
	else
		AddwM2CtoD<false><<<grid, threads>>>(n, w, Vg, Vg2);
}

void AddwCtoC(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vig, double2 *Vog){
    if(add)
		AddwCtoC<true><<<grid, threads>>>(n, w, Vig, Vog);
	else
		AddwCtoC<false><<<grid, threads>>>(n, w, Vig, Vog);
}

void AddwCwM2CtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double2 *Vga, double *Vg2a){
    if(add)
		AddwCwM2CtoCD<true><<<grid, threads>>>(n, w, Vg, Vga, Vg2a);
	else
		AddwCwM2CtoCD<false><<<grid, threads>>>(n, w, Vg, Vga, Vg2a);
}

void AddwCwDtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double *Vg2, double2 *Vga, double *Vg2a){
    if(add)
		AddwCwDtoCD<true><<<grid, threads>>>(n, w, Vg, Vg2, Vga, Vg2a);
	else
		AddwCwDtoCD<false><<<grid, threads>>>(n, w, Vg, Vg2, Vga, Vg2a);
}

void fft2ShiftD(dim3 grid, dim3 threads, int nxh, int nyh, double *Mg){
    fft2ShiftD<<<grid, threads>>>(nxh, nyh, Mg);
}

void fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, double2 *Mg){
    fft2ShiftC<<<grid, threads>>>(nxh, nyh, Mg);
}

/******************************************************************************/
void SetRetCount(unsigned int ih) {
	cudaMemcpyToSymbol(stRetCount, &ih, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

double getTotSumd(){
	double r;
	cudaMemcpyFromSymbol(&r, stTotSumd, cSizeofRD, 0, cudaMemcpyDeviceToHost);
	return r;
}

// Sum over all the elements
double SumAd(dim3 grid, dim3 threads, unsigned int nxy, const double *A, double *TB){
	bool nIsPow2 = (nxy&(nxy-1))==0;
	SetRetCount();
	SumAd<<<grid, threads>>>(nxy, nIsPow2, A, TB);
	return getTotSumd();
}

// Sum module square over all the elements
double SumAc2(dim3 grid, dim3 threads, unsigned int nxy, const double2 *A, double *TB){
	bool nIsPow2 = (nxy&(nxy-1))==0;
	SetRetCount();
	SumAc2<<<grid, threads>>>(nxy, nIsPow2, A, TB);
	return getTotSumd();
}

// Sum module square over all the elements with index
double SumAc2Ind(dim3 grid, dim3 threads, int nind, int *ind, const double2 *A, double *TB){
	bool nIsPow2 = (nind&(nind-1))==0;
	SetRetCount();
	SumAc2Ind<<<grid, threads>>>(nind, nIsPow2, ind, A, TB);
	return getTotSumd();
}

/******************************************************************************/
void ReadQuadratureGPU(int typ, int nQGPU, sQ1 &QGPU){
	sQ1 QCPU;
	QCPU.x = new double[nQGPU]; 
	QCPU.w = new double[nQGPU];
	cQuadrature Quad;
	Quad.ReadQuadrature(typ, nQGPU, QCPU);
	cudaMalloc((void**)&(QGPU.x), nQGPU*cSizeofRD);
	cudaMemcpy(QGPU.x, QCPU.x, nQGPU*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(QGPU.w), nQGPU*cSizeofRD);
	cudaMemcpy(QGPU.w, QCPU.w, nQGPU*cSizeofRD, cudaMemcpyHostToDevice);
	delete [] QCPU.x;
	delete [] QCPU.w;
}

void ReadQuadratureCPUGPU(int typ, int nQ, sQ1 &QCPU, sQ1 &QGPU){
	QCPU.x = new double[nQ]; 
	QCPU.w = new double[nQ];
	cQuadrature Quad;
	Quad.ReadQuadrature(typ, nQ, QCPU);
	cudaMalloc((void**)&(QGPU.x), nQ*cSizeofRD);
	cudaMemcpy(QGPU.x, QCPU.x, nQ*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(QGPU.w), nQ*cSizeofRD);
	cudaMemcpy(QGPU.w, QCPU.w, nQ*cSizeofRD, cudaMemcpyHostToDevice);
}