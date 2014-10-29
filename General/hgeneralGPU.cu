#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hgeneralGPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

/****************************************************************************/
__global__ void k_ScaleVectorD(int n, double f, double *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i] *= f;
		i += ii;
	}
}

// Scale Complex vector
__global__ void k_ScaleVectorC(int n, double f, double2 *Vg)
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
__global__ void k_SetValueVectorD(int n, double v, double *Vg1, double *Vg2)
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
__global__ void k_SetValueVectorD(int n, double v, double *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i] = v;
		i += ii;
	}
}

// Set value to Complex vector
__global__ void k_SetValueVectorC(int n, double vr, double vi, double2 *Vg)
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
__global__ void k_SetValueVectorCD(int n, double vr, double vi, double2 *Vcg, double v, double *Vdg)
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
__global__ void k_GetRealImagVectorC(int n, double2 *Vg, sComplex Vgri)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vgri.real[i] = Vg[i].x;
		Vgri.imag[i] = Vg[i].y;
		i += ii;
	}
}

// Set Real and Imaginary part of a Complex vector
__global__ void k_SetRealImagVectorC(int n, sComplex Vgri, double2 *Vg)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	while (i < n){
		Vg[i].x = Vgri.real[i];
		Vg[i].y = Vgri.imag[i];
		i += ii;
	}
}

// Add Mi to Mo and multiply by a weighted factor
template <bool add>
__global__ void k_AddwDtoD(int n, double w, double *Vig, double *Vog)
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
__global__ void k_AddwM2CtoD(int n, double w, double2 *Vg, double *Vg2)
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
__global__ void k_AddwCtoC(int n, double w, double2 *Vig, double2 *Vog)
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
__global__ void k_AddwCwM2CtoCD(int n, double w, double2 *Vg, double2 *aVg, double *aVg2)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double x, y;

	while (i < n){
		x = Vg[i].x; y = Vg[i].y;
		if (add){
			aVg[i].x += w*x; aVg[i].y += w*y;
			aVg2[i] += w*(x*x + y*y);
		}else{
			aVg[i].x = w*x; aVg[i].y = w*y;
			aVg2[i] = w*(x*x + y*y);
		}
		i += ii;
	}
}

// Add Mc to aMc and MMca2
template <bool add>
__global__ void k_AddwCwDtoCD(int n, double w, double2 *Vg, double *Vg2, double2 *aVg, double *aVg2)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double x, y;

	while (i < n){
		x = Vg[i].x; y = Vg[i].y;
		if (add){
			aVg[i].x += w*x; aVg[i].y += w*y;
			aVg2[i] += w*Vg2[i];
		}else{
			aVg[i].x = w*x; aVg[i].y = w*y;
			aVg2[i] = w*Vg2[i];
		}
		i += ii;
	}
}

// Shift Double matrix respect to (nxh, nyh)
__global__ void k_fft2ShiftD(int nxh, int nyh, double *Mg)
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
__global__ void k_fft2ShiftC(int nxh, int nyh, double2 *Mg)
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

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2ShiftC(int nxh, int nyh, sComplex Mg)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < nxh)&&(j < nyh)){
		int ny = 2*nyh, k1, k2;
		double zr, zi;

		k1 = j + i*ny; k2 = (nyh+j)+(nxh+i)*ny;
		zr = Mg.real[k1]; zi = Mg.imag[k1];
		Mg.real[k1] = Mg.real[k2]; Mg.imag[k1] = Mg.imag[k2];
		Mg.real[k2] = zr; Mg.imag[k2] = zi;

		k1 = (nyh+j) + i*ny; k2 = j+(nxh+i)*ny;
		zr = Mg.real[k1]; zi = Mg.imag[k1];
		Mg.real[k1] = Mg.real[k2]; Mg.imag[k1] = Mg.imag[k2];
		Mg.real[k2] = zr; Mg.imag[k2] = zi;
	}
}

// Sum over all the elements
__global__ void k_SumMd(int nxy, bool nIsPow2, double *M, double *Mp){ 
	__shared__ double Ms[thrnxy];
	unsigned int tid = threadIdx.x, gridSize = 2*thrnxy*gridDim.x;
	unsigned int i =  threadIdx.x + 2*thrnxy*blockIdx.x;
	double SumT = 0;

	while (i < nxy){
		SumT += M[i];
		if (nIsPow2 || i+thrnxy < nxy)
			SumT += M[i+thrnxy];
		i += gridSize;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid == 0) Mp[blockIdx.x] = Ms[0];
	
	if (gridDim.x > 1){
		__shared__ bool amLast;

		__threadfence();

		if (tid == 0){
			unsigned int ticket = atomicInc(&stRetCount, gridDim.x);
			amLast = (ticket == gridDim.x-1);
			}
			__syncthreads();

			if (amLast){
				i = tid;
				SumT = 0;
				while (i < gridDim.x){
					SumT += Mp[i];
					i += thrnxy;
				}
				Ms[tid] = SumT;
				__syncthreads();

				k_reduceBlockDouble256(Ms, tid);

				if (tid == 0){
					stTotSumd = Ms[0];
					stRetCount = 0;
				}
			}
	}
}

// Sum module square over all the elements
__global__ void k_SumMC2(int nxy, bool nIsPow2, double2 *M, double *Mp){ 
	__shared__ double Ms[thrnxy];
	unsigned int tid = threadIdx.x, gridSize = 2*thrnxy*gridDim.x;
	unsigned int i =  threadIdx.x + 2*thrnxy*blockIdx.x;
	double x, y, SumT = 0;

	while (i < nxy){
		x = M[i].x; y = M[i].y;
		SumT += x*x + y*y;
		if (nIsPow2 || i+thrnxy < nxy){
			x = M[i+thrnxy].x; y = M[i+thrnxy].y;
			SumT += x*x + y*y;
		}
		i += gridSize;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid == 0) Mp[blockIdx.x] = Ms[0];
	
	if (gridDim.x > 1){
		__shared__ bool amLast;

		__threadfence();

		if (tid == 0){
			unsigned int ticket = atomicInc(&stRetCount, gridDim.x);
			amLast = (ticket == gridDim.x-1);
		}
		__syncthreads();

		if (amLast){
			i = tid;
			SumT = 0;
			while (i < gridDim.x){
				SumT += Mp[i];
				i += thrnxy;
			}
			Ms[tid] = SumT;
			__syncthreads();

			k_reduceBlockDouble256(Ms, tid);

			if (tid == 0){
				stTotSumd = Ms[0];
				stRetCount = 0;
			}
		}
	}
}

/************************************************************************************/
/************************************************************************************/
// Sum module square over all the elements
__global__ void k_SumMDet(sGP GP, double *M, double gmin2, double gmax2, double *Mp){ 
	__shared__ double Ms[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double gx, gy, g2, b;
	double SumT=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			b = ((gmin2<=g2)&&(g2<=gmax2))?1.0:0.0;
			SumT += b*M[ixy];
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		Mp[bid] = Ms[0];
}

// Sum module square over all the elements
__global__ void k_SumM(int n, double *Mp, int iSum, double *SumM){ 
	__shared__ double Ms[thrnxy];

	int tid = threadIdx.x;
	int ix =  threadIdx.x + blockIdx.x*blockDim.x;
	int gridSize = blockDim.x*gridDim.x;
	double SumT=0;

	while (ix < n){
		SumT += Mp[ix];
		ix += gridSize;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		SumM[iSum] = Ms[0];
}

/************************************************************************************/
/************************************************************************************/
// Sum module square over all the elements
__global__ void k_SumM1M2Det(sGP GP, double *M1, double *M2, double gmin2, double gmax2, double *M1p, double *M2p){ 
	__shared__ double Ms1[thrnxny*thrnxny];
	__shared__ double Ms2[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double gx, gy, g2, t;
	double SumT1=0, SumT2=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			t = ((gmin2<=g2)&&(g2<=gmax2))?1.0:0.0;
			SumT1 += t*M1[ixy];
			SumT2 += t*M2[ixy];	
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms1[tid] = SumT1;
	Ms2[tid] = SumT2;
	__syncthreads();

	k_reduceBlockDouble256(Ms1, Ms2, tid);

	if (tid==0){
		M1p[bid] = Ms1[0];
		M2p[bid] = Ms2[0];
	}
}

// Sum module square over all the elements
__global__ void k_SumM1M2(int n, double *M1p, double *M2p, int iSum, double *SumM1, double *SumM2){ 
	__shared__ double Ms1[thrnxy];
	__shared__ double Ms2[thrnxy];

	int tid = threadIdx.x;
	int ix =  threadIdx.x + blockIdx.x*blockDim.x;
	int gridSize = blockDim.x*gridDim.x;
	double SumT1=0, SumT2=0;

	while (ix < n){
		SumT1 += M1p[ix];
		SumT2 += M2p[ix];
		ix += gridSize;
	}
	Ms1[tid] = SumT1;
	Ms2[tid] = SumT2;
	__syncthreads();

	k_reduceBlockDouble256(Ms1, Ms2, tid);

	if (tid==0){
		SumM1[iSum] = Ms1[0];
		SumM2[iSum] = Ms2[0];
	}
}

/************************************************************************************/
/************************************************************************************/
// Phase
__global__ void k_Phase(sGP GP, double gxu, double gyu, sERg ERgg){
	int ixy = threadIdx.x + blockIdx.x*blockDim.x;
	double theta, R;

	if (ixy < GP.nx){
		R = IsRS(ixy,GP.nxh)*GP.dRx;
		theta = c2Pi*R*gxu;
		sincos(theta, &(ERgg.x[ixy].y), &(ERgg.x[ixy].x));
	}

	if (ixy < GP.ny){
		R = IsRS(ixy,GP.nyh)*GP.dRy;
		theta = c2Pi*R*gyu;
		sincos(theta, &(ERgg.y[ixy].y), &(ERgg.y[ixy].x));
	}
}

// Phase multiplication
__global__ void k_PhaseMul(sGP GP, sERg ERgg, double2 *Psig){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int ixy = i*GP.ny+j;
		double z1r = ERgg.x[i].x, z1i = ERgg.x[i].y;
		double z2r = ERgg.y[j].x, z2i = ERgg.y[j].y;
		double z3r = z1r*z2r-z1i*z2i, z3i = z1i*z2r+z1r*z2i;
		z2r = Psig[ixy].x;		z2i = Psig[ixy].y; 
		Psig[ixy].x = z3r*z2r-z3i*z2i; 
		Psig[ixy].y = z3i*z2r+z3r*z2i; 
	}
}

// Calculated transmission function: Exp(i*sigma*Vpg)
__global__ void k_TransmissionWPO(sGP GP, double f, double *Vpg, double2 *Transg){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		Transg[ixy].x = 1.0;
		Transg[ixy].y = f*Vpg[ixy];
	}
}

// Calculated transmission function: Exp(i*sigma*Vpg)
__global__ void k_Transmission1(sGP GP, double f, double *Vpg, double2 *Transg){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double theta = f*Vpg[ixy];
		sincos(theta, &(Transg[ixy].y) , &(Transg[ixy].x));
	}
}

// Calculated transmission function
__global__ void k_Transmission2(sGP GP, eSlicePos SlicePo, double f, double *Vpg, double *zVpg, double *zVpog, double2 *Transg){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double VR = Vpg[ixy], zVp = zVpg[ixy];
		double theta;
		switch (SlicePo){
			case eSPFirst: // initial slice
				theta = f*(VR-zVp);
				zVpog[ixy] = zVp;
				break;
			case eSPMedium: // intermediate slice
				theta = f*(VR-zVp+zVpog[ixy]);
				zVpog[ixy] = zVp;
				break;
			case eSPLast: // last slice
				theta = f*zVpog[ixy];
				break;
		}	
		sincos(theta, &(Transg[ixy].y) , &(Transg[ixy].x));
	}
}

// Anti-Aliasing, scale with cut-off (2/3)gmax
__global__ void k_BandwidthLimited2D(sGP GP, double2 *Transg){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if (g2 < GP.gmaxl2){
			Transg[ixy].x *= GP.inxy;
			Transg[ixy].y *= GP.inxy;
		}else{
			Transg[ixy].x = 0.0; 
			Transg[ixy].y = 0.0; 
		}
	}
}

// Element by element multiplication
__global__ void k_Transmit(sGP GP, double2 *Transg, double2 *Psig){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double z1r = Transg[ixy].x, z1i = Transg[ixy].y;
		double z2r = Psig[ixy].x, z2i = Psig[ixy].y;
		Psig[ixy].x = z1r*z2r-z1i*z2i;
		Psig[ixy].y = z1i*z2r+z1r*z2i;
	}
}

// Build propagator function
__global__ void k_Propagator(sGP GP, double gxu, double gyu, double scale, sProp Propg){
	int ixy = threadIdx.x + blockIdx.x*blockDim.x;
	double theta, g;

	if (ixy < GP.nx){
		g = IsFS(ixy,GP.nxh)*GP.dgx;
		theta = scale*(g+gxu)*(g+gxu);
		sincos(theta, &(Propg.x[ixy].y), &(Propg.x[ixy].x));
	}
	if (ixy < GP.ny){
		g = IsFS(ixy,GP.nyh)*GP.dgy;
		theta = scale*(g+gyu)*(g+gyu);
		sincos(theta, &(Propg.y[ixy].y), &(Propg.y[ixy].x));
	}
}

// Propagate, scale with cut-off (2/3)gmax
__global__ void k_Propagate(sGP GP, sProp Propg, double2 *Psig){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if (g2 < GP.gmaxl2){
			double z1r = Propg.x[ix].x, z1i = Propg.x[ix].y;
			double z2r = Propg.y[iy].x, z2i = Propg.y[iy].y;
			double z3r = z1r*z2r-z1i*z2i, z3i = z1i*z2r+z1r*z2i;
			z2r = Psig[ixy].x;		z2i = Psig[ixy].y; 
			Psig[ixy].x = (z3r*z2r-z3i*z2i)*GP.inxy;
			Psig[ixy].y = (z3i*z2r+z3r*z2i)*GP.inxy; 
		}else{
 			Psig[ixy].x = 0.0; 
			Psig[ixy].y = 0.0; 
		}
	}
}

/************************************************************************************/
/************************************************************************************/
void f_ScaleVectorD(dim3 grid, dim3 threads, int n, double f, double *&Vg){
 k_ScaleVectorD<<<grid, threads>>>(n, f, Vg);
}

void f_ScaleVectorC(dim3 grid, dim3 threads, int n, double f, double2 *&Vg){
 k_ScaleVectorC<<<grid, threads>>>(n, f, Vg);
}

void f_SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *&Vg1, double *&Vg2){
 k_SetValueVectorD<<<grid, threads>>>(n, v, Vg1, Vg2);
}

void f_SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *&Vg){
 k_SetValueVectorD<<<grid, threads>>>(n, v, Vg);
}

void f_SetValueVectorC(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *&Vg){
 k_SetValueVectorC<<<grid, threads>>>(n, vr, vi, Vg);
}

void f_SetValueVectorCD(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *&Vcg, double v, double *&Vdg){
 k_SetValueVectorCD<<<grid, threads>>>(n, vr, vi, Vcg, v, Vdg);
}

void f_GetRealImagVectorC(dim3 grid, dim3 threads, int n, double2 *&Vg, sComplex &Vgri){
 k_GetRealImagVectorC<<<grid, threads>>>(n, Vg, Vgri);
}

void f_SetRealImagVectorC(dim3 grid, dim3 threads, int n, sComplex &Vgri, double2 *&Vg){
 k_SetRealImagVectorC<<<grid, threads>>>(n, Vgri, Vg);
}

void f_AddwDtoD(dim3 grid, dim3 threads, bool add, int n, double w, double *&Vig, double *&Vog){
 if(add)
		k_AddwDtoD<true><<<grid, threads>>>(n, w, Vig, Vog);
	else
		k_AddwDtoD<false><<<grid, threads>>>(n, w, Vig, Vog);
}

void f_AddwM2CtoD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double *&Vg2){
 if(add)
		k_AddwM2CtoD<true><<<grid, threads>>>(n, w, Vg, Vg2);
	else
		k_AddwM2CtoD<false><<<grid, threads>>>(n, w, Vg, Vg2);
}

void f_AddwCtoC(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vig, double2 *&Vog){
 if(add)
		k_AddwCtoC<true><<<grid, threads>>>(n, w, Vig, Vog);
	else
		k_AddwCtoC<false><<<grid, threads>>>(n, w, Vig, Vog);
}

void f_AddwCwM2CtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double2 *&aVg, double *&aVg2){
 if(add)
		k_AddwCwM2CtoCD<true><<<grid, threads>>>(n, w, Vg, aVg, aVg2);
	else
		k_AddwCwM2CtoCD<false><<<grid, threads>>>(n, w, Vg, aVg, aVg2);
}

void f_AddwCwDtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *&Vg, double *&Vg2, double2 *&aVg, double *&aVg2){
 if(add)
		k_AddwCwDtoCD<true><<<grid, threads>>>(n, w, Vg, Vg2, aVg, aVg2);
	else
		k_AddwCwDtoCD<false><<<grid, threads>>>(n, w, Vg, Vg2, aVg, aVg2);
}

void f_fft2ShiftD(dim3 grid, dim3 threads, int nxh, int nyh, double *&Mg){
 k_fft2ShiftD<<<grid, threads>>>(nxh, nyh, Mg);
}

void f_fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, double2 *&Mg){
 k_fft2ShiftC<<<grid, threads>>>(nxh, nyh, Mg);
}

void f_fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, sComplex &Mg){
 k_fft2ShiftC<<<grid, threads>>>(nxh, nyh, Mg);
}

void f_setRetCount(int ih){
	cudaMemcpyToSymbol(stRetCount, &ih, sizeof(int), 0, cudaMemcpyHostToDevice);
}

double f_getTotSumd(){
	double r;
	cudaMemcpyFromSymbol(&r, stTotSumd, cSizeofRD, 0, cudaMemcpyDeviceToHost);
	return r;
}

// Sum over all the elements
double f_SumMd(dim3 grid, dim3 threads, int nxy, double *&M, double *&Mb){
	bool nIsPow2 = (nxy&(nxy-1)) == 0;
	f_setRetCount(0);
	k_SumMd<<<grid, threads>>>(nxy, nIsPow2, M, Mb);
	return f_getTotSumd();
}

// Sum module square over all the elements
double f_SumMC2(dim3 grid, dim3 threads, int nxy, double2 *&M, double *&Mb){
	bool nIsPow2 = (nxy&(nxy-1)) == 0;
	f_setRetCount(0);
	k_SumMC2<<<grid, threads>>>(nxy, nIsPow2, M, Mb);
	return f_getTotSumd();
}

// Sum all the elements
void f_SumMDet(dim3 grid, dim3 threads, sGP &GP, double *&M, double gmin2, double gmax2, double *&Mp, int iSum, double *&SumM){
	k_SumMDet<<<grid, threads>>>(GP, M, gmin2, gmax2, Mp);
	k_SumM<<<1, thrnxy>>>(grid.x*grid.y, Mp, iSum, SumM);
}

// Sum all the elements
void f_SumM1M2Det(dim3 grid, dim3 threads, sGP &GP, double *&M1, double *&M2, double gmin2, double gmax2, double *&M1p, double *&M2p, int iSum, double *&SumM1, double *&SumM2){
	k_SumM1M2Det<<<grid, threads>>>(GP, M1, M2, gmin2, gmax2, M1p, M2p);
	k_SumM1M2<<<1, thrnxy>>>(grid.x*grid.y, M1p, M2p, iSum, SumM1, SumM2);
}

/****************************************************************************/
/****************************************************************************/

void f_ReadQuadratureGPU(int typ, int nQGPU, sQ1 &QGPU){
	if(nQGPU<=0) return;

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

void f_ReadQuadratureCPUGPU(int typ, int nQ, sQ1 &QCPU, sQ1 &QGPU){
	if(nQ<=0) return;

	QCPU.x = new double[nQ]; 
	QCPU.w = new double[nQ];
	cQuadrature Quad;
	Quad.ReadQuadrature(typ, nQ, QCPU);
	cudaMalloc((void**)&(QGPU.x), nQ*cSizeofRD);
	cudaMemcpy(QGPU.x, QCPU.x, nQ*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(QGPU.w), nQ*cSizeofRD);
	cudaMemcpy(QGPU.w, QCPU.w, nQ*cSizeofRD, cudaMemcpyHostToDevice);
}

/****************************************************************************/
/****************************************************************************/

void f_sCoefPar_cudaFree(sCoefPar &CoefPar){
	cudaFreen(CoefPar.cl);
	cudaFreen(CoefPar.cnl);
}

void f_sCoefPar_cudaInit(sCoefPar &CoefPar){
	CoefPar.cl = 0;
	CoefPar.cnl = 0;
}

void f_sCoefPar_cudaMalloc(int nCoefPar, sCoefPar &CoefPar){
	if(nCoefPar<=0) return;

	cudaMalloc((void**)&(CoefPar.cl), nCoefPar*cSizeofRD);
	cudaMalloc((void**)&(CoefPar.cnl), nCoefPar*cSizeofRD);
}

/****************************************************************************/
/****************************************************************************/

void f_sciVn_cudaFree(sciVn &ciVn){
	cudaFreen(ciVn.c0);
	cudaFreen(ciVn.c1);
	cudaFreen(ciVn.c2);
	cudaFreen(ciVn.c3);
}

void f_sciVn_cudaInit(sciVn &ciVn){
	ciVn.c0 = 0;
	ciVn.c1 = 0;
	ciVn.c2 = 0;
	ciVn.c3 = 0;
}

void f_sciVn_cudaMalloc(int nciVn, sciVn &ciVn){
	if(nciVn<=0) return;

	cudaMalloc((void**)&(ciVn.c0), nciVn*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c1), nciVn*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c2), nciVn*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c3), nciVn*cSizeofRD);
}

/****************************************************************************/
/****************************************************************************/

void f_sDetCir_cudaFree(sDetCir &DetCir){
	cudaFreen(DetCir.g2min);
	cudaFreen(DetCir.g2max);
}

void f_sDetCir_cudaInit(sDetCir &DetCir){
	DetCir.g2min = 0;
	DetCir.g2max = 0;
}

void f_sDetCir_cudaMalloc(int nDetCir, sDetCir &DetCir){
	if(nDetCir<=0) return;

	cudaMalloc((void**)&(DetCir.g2min), nDetCir*cSizeofRD);
	cudaMalloc((void**)&(DetCir.g2max), nDetCir*cSizeofRD);
}

/****************************************************************************/
/****************************************************************************/

// From Device To Host
void Cd_2_Ch(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, sComplex &Psih){	
	dim3 grid(MIN(64, (nxy+thrnxy-1)/thrnxy), 1, 1);
	dim3 threads(thrnxy, 1, 1);
	sComplex Psi;

	Psi.real = Psi_real; Psi.imag = Psi_Imag;
	// Get Real and Imaginary part to complex matrix
	f_GetRealImagVectorC(grid, threads, nxy, Psid, Psi);
	// Copy real part of the wave function to the host
	cudaMemcpy(Psih.real, Psi.real, nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	// Copy imaginary part of the wave function to the host
	cudaMemcpy(Psih.imag, Psi.imag, nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void CdDd_2_ChDh(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih){
	Cd_2_Ch(nxy, Psi_real, Psi_Imag, Psid, Psih);
	// Copy wave function squared to the host
	cudaMemcpy(M2Psih, M2Psid, nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void CdDdDd_2_ChDhDh(int nxy, double *&Psi_real, double *&Psi_Imag, double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h){	
	Cd_2_Ch(nxy, Psi_real, Psi_Imag, Psid, Psih);
	// Copy wave function squared to the host
	cudaMemcpy(M2Psi1h, M2Psi1d, nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(M2Psi2h, M2Psi2d, nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

/****************************************************************************/
/****************************************************************************/

// Phase multiplication
void f_PhaseMul(dim3 grid, dim3 threads, sGP &GP, double gxu, double gyu, sERg &ERg, double2 *&Psi){
	// exp(2*1i*pi(rx*gxu+ry*gyu)
	k_Phase<<<grid, threads>>>(GP, gxu, gyu, ERg);
	// Phase multiplication
	k_PhaseMul<<<grid, threads>>>(GP, ERg, Psi);
}

// AntiAliasing
void f_BandwidthLimited2D(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double2 *Psi){
	// Forward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
	// AntiAliasing, scale and bandlimited
	k_BandwidthLimited2D<<<grid, threads>>>(GP, Psi);
	// Backward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
}

// Transmission function WPO
void f_TransmissionWPO(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double2 *&Trans){
	// Transmission
	k_TransmissionWPO<<<grid, threads>>>(GP, fPot, V0, Trans);
	// AntiAliasing
	f_BandwidthLimited2D(grid, threads, PlanPsi, GP, Trans);
}

// Transmission function
void f_Transmission1(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double2 *&Trans){
	// Transmission
	k_Transmission1<<<grid, threads>>>(GP, fPot, V0, Trans);
	// AntiAliasing
	f_BandwidthLimited2D(grid, threads, PlanPsi, GP, Trans);
}

// Transmission function
void f_Transmission2(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double fPot, double *&V0, double *&V1, double *&V2, eSlicePos SlicePos, double2 *&Trans){
	// Transmission
	k_Transmission2<<<grid, threads>>>(GP, SlicePos, fPot, V0, V1, V2, Trans);
	// AntiAliasing
	f_BandwidthLimited2D(grid, threads, PlanPsi, GP, Trans);
}

// Trasmit
void f_Transmit(dim3 grid, dim3 threads, sGP &GP, double2 *&Trans, double2 *&Psi){
	k_Transmit<<<grid, threads>>>(GP, Trans, Psi);
}

// Propagate
void f_Propagate(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sProp &Prop, double2 *&Psi){
	// Forward propagator
	k_Propagator<<<grid, threads>>>(GP, gxu, gyu, -cPi*lambda*z, Prop);
	// Forward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
	// Propagate, scale and bandlimited the wave function
	k_Propagate<<<grid, threads>>>(GP, Prop, Psi);
	// Backward fft2
	if(Space == eSReal)
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
}