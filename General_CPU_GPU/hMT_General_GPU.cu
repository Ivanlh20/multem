/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

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

/***************************************************************************/
/***************************************************************************/

__global__ void k_Scale_MD(sGP GP, double w, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MD_io[ixy] *= w;
	}
}

__global__ void k_Scale_MD(sGP GP, double w, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MD1_io[ixy] *= w;
		MD2_io[ixy] *= w;
	}
}

__global__ void k_Scale_MC(sGP GP, double w, double2 * __restrict MC_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC_io[ixy].x *= w;
		MC_io[ixy].y *= w;
	}
}

__global__ void k_Scale_MC(sGP GP, double w, double2 * __restrict MC1_io, double2 * __restrict MC2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC1_io[ixy].x *= w;
		MC1_io[ixy].y *= w;
		MC2_io[ixy].x *= w;
		MC2_io[ixy].y *= w;
	}
}

void f_Scale_MD(sGP &GP, double w, double *&MD_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Scale_MD<<<Bnxny, Tnxny>>>(GP, w, MD_io);
}

void f_Scale_MD(sGP &GP, double w, double *&MD1_io, double *&MD2_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Scale_MD<<<Bnxny, Tnxny>>>(GP, w, MD1_io, MD2_io);
}

void f_Scale_MC(sGP &GP, double w, double2 *&MC_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Scale_MC<<<Bnxny, Tnxny>>>(GP, w, MC_io);
}

void f_Scale_MC(sGP &GP, double w, double2 *&MC1_io, double2 *&MC2_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Scale_MC<<<Bnxny, Tnxny>>>(GP, w, MC1_io, MC2_io);
}

/***************************************************************************/
/***************************************************************************/

// Set value to Double vector:
__global__ void k_Set_MD(sGP GP, double M, double * __restrict MD_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MD_o[ixy] = M;
	}
}

// Set value to 2 Double vector:
__global__ void k_Set_MD(sGP GP, double M, double * __restrict MD1_o, double * __restrict MD2_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MD1_o[ixy] = M;
		MD2_o[ixy] = M;
	}
}

// Set value to Complex vector
__global__ void k_Set_MC(sGP GP, double Mr, double Mi, double2 * __restrict MC_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC_o[ixy].x = Mr;
		MC_o[ixy].y = Mi;
	}
}

// Set value to 2 Double vector:
__global__ void k_Set_MC(sGP GP, double Mr, double Mi, double2 * __restrict MC1_o, double2 * __restrict MC2_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC1_o[ixy].x = Mr;
		MC1_o[ixy].y = Mi;
		MC2_o[ixy].x = Mr;
		MC2_o[ixy].y = Mi;
	}
}

// Set value to Complex vector and Double vector
__global__ void k_Set_MC_MD(sGP GP, double Mr, double Mi, double2 * __restrict MC_o, double M, double * __restrict MD_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC_o[ixy].x = Mr;
		MC_o[ixy].y = Mi;
		MD_o[ixy] = M;
	}
}

// Set Real and Imaginary part of a Complex vector
__global__ void k_Set_MC(sGP GP, const sComplex MC_i, double2 * __restrict MC_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC_o[ixy].x = MC_i.real[ixy];
		MC_o[ixy].y = MC_i.imag[ixy];
	}
}

// Get Real and Imaginary part of a Complex vector
__global__ void k_Get_MC(sGP GP, const double2 * __restrict MC_i, sComplex MC_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		MC_o.real[ixy] = MC_i[ixy].x;
		MC_o.imag[ixy] = MC_i[ixy].y;
	}
}

void f_Set_MD(sGP &GP, double M, double *&MD_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MD<<<Bnxny, Tnxny>>>(GP, M, MD_o);
}

void f_Set_MD(sGP &GP, double M, double *&MD1_o, double *&MD2_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MD<<<Bnxny, Tnxny>>>(GP, M, MD1_o, MD2_o);
}

void f_Set_MC(sGP &GP, double Mr, double Mi, double2 *&MC_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MC<<<Bnxny, Tnxny>>>(GP, Mr, Mi, MC_o);
}

void f_Set_MC(sGP &GP, double Mr, double Mi, double2 *&MC1_o, double2 *&MC2_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MC<<<Bnxny, Tnxny>>>(GP, Mr, Mi, MC1_o, MC2_o);
}

void f_Set_MC_MD(sGP &GP, double Mr, double Mi, double2 *&MC_o, double M, double *&MD_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MC_MD<<<Bnxny, Tnxny>>>(GP, Mr, Mi, MC_o, M, MD_o);
}

void f_Set_MC(sGP &GP, sComplex &MC_i, double2 *&MC_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Set_MC<<<Bnxny, Tnxny>>>(GP, MC_i, MC_o);
}

void f_Get_MC(sGP &GP, double2 *&MC_i, sComplex &MC_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Get_MC<<<Bnxny, Tnxny>>>(GP, MC_i, MC_o);
}

/***************************************************************************/
/***************************************************************************/

template <bool add>
__global__ void k_Add_wMD(sGP GP, double w, const double * __restrict MD_i, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x = MD_i[ixy];
		if (add)
			MD_io[ixy] += x = w*x;
		else
			MD_io[ixy] = x = w*x;
	}
}

template <bool add>
__global__ void k_Add_wMD(sGP GP, double w, const double * __restrict MD1_i, const double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x1 = MD1_i[ixy], x2 = MD2_i[ixy];
		if (add){
			MD1_io[ixy] += x1 = w*x1;
			MD2_io[ixy] += x2 = w*x2;
		}else{
			MD1_io[ixy] = x1 = w*x1;
			MD2_io[ixy] = x2 = w*x2;
		}
	}
}

template <bool add>
__global__ void k_Add_wMC(sGP GP, double w, const double2 * __restrict MC_i, double2 * __restrict MC_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x = MC_i[ixy].x, y = MC_i[ixy].y;
		if (add){
			MC_io[ixy].x += x = w*x;
			MC_io[ixy].y += y = w*y;
		}else{
			MC_io[ixy].x = x = w*x;
			MC_io[ixy].y = y = w*y;
		}
	}
}

template <bool add>
__global__ void k_Add_wMC(sGP GP, double w, const double2 * __restrict MC1_i, const double2 * __restrict MC2_i, double2 * __restrict MC1_io, double2 * __restrict MC2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x1 = MC1_i[ixy].x, y1 = MC1_i[ixy].y;
		double x2 = MC2_i[ixy].x, y2 = MC2_i[ixy].y;
		if (add){
			MC1_io[ixy].x += x1 = w*x1;
			MC1_io[ixy].y += y1 = w*y1;
			MC2_io[ixy].x += x2 = w*x2;
			MC2_io[ixy].y += y2 = w*y2;
		}else{
			MC1_io[ixy].x = x1 = w*x1;
			MC1_io[ixy].y = y1 = w*y1;
			MC2_io[ixy].x = x2 = w*x2;
			MC2_io[ixy].y = y2 = w*y2;
		}
	}
}

template <bool add>
__global__ void k_Add_wMC2(sGP GP, double w, const double2 * __restrict MC_i, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x = MC_i[ixy].x, y = MC_i[ixy].y, z = x*x+y*y;
		if (add)
			MD_io[ixy] += z = w*z;
		else
			MD_io[ixy] = z = w*z;
	}
}

template <bool add>
__global__ void k_Add_wMC_wMD(sGP GP, double w, const double2 * __restrict MC_i, double2 * __restrict MC_io, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x = MC_i[ixy].x, y = MC_i[ixy].y, z = x*x+y*y;
		if (add){
			MC_io[ixy].x += x = w*x; 
			MC_io[ixy].y += y = w*y;
			MD_io[ixy] += z = w*z;
		}else{
			MC_io[ixy].x = x = w*x; 
			MC_io[ixy].y = y = w*y;
			MD_io[ixy] = z = w*z;
		}
	}
}

template <bool add>
__global__ void k_Add_wMC_wMD(sGP GP, double w, const double2 * __restrict MC_i, const double * __restrict MD_i, double2 * __restrict MC_io, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double x = MC_i[ixy].x, y = MC_i[ixy].y, z = MD_i[ixy];
		if (add){
			MC_io[ixy].x += x = w*x; 
			MC_io[ixy].y += y = w*y;
			MD_io[ixy] += z = w*z;
		}else{
			MC_io[ixy].x = x = w*x; 
			MC_io[ixy].y = y = w*y;
			MD_io[ixy] = z = w*z;
		}
	}
}

void f_Add_wMD(bool add, sGP &GP, double w, double *&MD_i, double *&MD_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMD<true><<<Bnxny, Tnxny>>>(GP, w, MD_i, MD_io);
	else
		k_Add_wMD<false><<<Bnxny, Tnxny>>>(GP, w, MD_i, MD_io);
}

void f_Add_wMD(bool add, sGP &GP, double w, double *&MD1_i, double *&MD2_i, double *&MD1_io, double *&MD2_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMD<true><<<Bnxny, Tnxny>>>(GP, w, MD1_i, MD2_i, MD1_io, MD2_io);
	else
		k_Add_wMD<false><<<Bnxny, Tnxny>>>(GP, w, MD1_i, MD2_i, MD1_io, MD2_io);
}

void f_Add_wMC(bool add, sGP &GP, double w, double2 *&MC_i, double2 *&MC_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMC<true><<<Bnxny, Tnxny>>>(GP, w, MC_i, MC_io);
	else
		k_Add_wMC<false><<<Bnxny, Tnxny>>>(GP, w, MC_i, MC_io);
}

void f_Add_wMC(bool add, sGP &GP, double w, double2 *&MC1_i, double2 *&MC2_i, double2 *&MC1_io, double2 *&MC2_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMC<true><<<Bnxny, Tnxny>>>(GP, w, MC1_i, MC2_i, MC1_io, MC2_io);
	else
		k_Add_wMC<false><<<Bnxny, Tnxny>>>(GP, w, MC1_i, MC2_i, MC1_io, MC2_io);
}

void f_Add_wMC2(bool add, sGP &GP, double w, double2 *&MC_i, double *&MD_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMC2<true><<<Bnxny, Tnxny>>>(GP, w, MC_i, MD_io);
	else
		k_Add_wMC2<false><<<Bnxny, Tnxny>>>(GP, w, MC_i, MD_io);
}

void f_Add_wMC_wMD(bool add, sGP &GP, double w, double2 *&MC_i, double2 *&MC_io, double *&MD_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMC_wMD<true><<<Bnxny, Tnxny>>>(GP, w, MC_i, MC_io, MD_io);
	else
		k_Add_wMC_wMD<false><<<Bnxny, Tnxny>>>(GP, w, MC_i, MC_io, MD_io);
}

void f_Add_wMC_wMD(bool add, sGP &GP, double w, double2 *&MC_i, double *&MD_i, double2 *&MC_io, double *&MD_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	if(add)
		k_Add_wMC_wMD<true><<<Bnxny, Tnxny>>>(GP, w, MC_i, MD_i, MC_io, MD_io);
	else
		k_Add_wMC_wMD<false><<<Bnxny, Tnxny>>>(GP, w, MC_i, MD_i, MC_io, MD_io);
}

/***************************************************************************/
/***************************************************************************/
__device__ inline void ExcVal(double &v1, double &v2)
{
	 double v = v1;
	 v1 = v2;
	 v2 = v; 
}

__device__ inline void ExcVal(double2 &v1, double2 &v2)
{
	 double vx = v1.x, vy = v1.y;
	 v1.x = v2.x; v1.y = v2.y;
	 v2.x = vx; v2.y = vy;
}

// Shift Double matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MD(sGP GP, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MD_io[ixy], MD_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MD_io[ixy], MD_io[ixys]);
	}
}

// Shift Double matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MD(sGP GP, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MD1_io[ixy], MD1_io[ixys]);
		ExcVal(MD2_io[ixy], MD2_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MD1_io[ixy], MD1_io[ixys]);
		ExcVal(MD2_io[ixy], MD2_io[ixys]);
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MC(sGP GP, double2 * __restrict MC_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MC_io[ixy], MC_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MC_io[ixy], MC_io[ixys]);
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MC(sGP GP, double2 * __restrict MC1_io, double2 * __restrict MC2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MC1_io[ixy], MC1_io[ixys]);
		ExcVal(MC2_io[ixy], MC2_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MC1_io[ixy], MC1_io[ixys]);
		ExcVal(MC2_io[ixy], MC2_io[ixys]);
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MC_MD(sGP GP, double2 * __restrict MC_io, double * __restrict MD_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MC_io[ixy], MC_io[ixys]);
		ExcVal(MD_io[ixy], MD_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MC_io[ixy], MC_io[ixys]);
		ExcVal(MD_io[ixy], MD_io[ixys]);
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MC_MD(sGP GP, double2 * __restrict MC_io, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		ExcVal(MC_io[ixy], MC_io[ixys]);
		ExcVal(MD1_io[ixy], MD1_io[ixys]);
		ExcVal(MD2_io[ixy], MD2_io[ixys]);

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		ExcVal(MC_io[ixy], MC_io[ixys]);
		ExcVal(MD1_io[ixy], MD1_io[ixys]);
		ExcVal(MD2_io[ixy], MD2_io[ixys]);
	}
}

// Shift Complex matrix respect to (nxh, nyh)
__global__ void k_fft2Shift_MC(sGP GP, sComplex MC_io)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nxh)&&(iy < GP.nyh)){
		int ixy, ixys;
		double zr, zi;

		ixy = ix*GP.ny + iy; ixys = (GP.nxh+ix)*GP.ny+(GP.nyh+iy);
		zr = MC_io.real[ixy]; zi = MC_io.imag[ixy];
		MC_io.real[ixy] = MC_io.real[ixys]; MC_io.imag[ixy] = MC_io.imag[ixys];
		MC_io.real[ixys] = zr; MC_io.imag[ixys] = zi;

		ixy = (GP.nyh+iy) + ix*GP.ny; ixys = (GP.nxh+ix)*GP.ny+iy;
		zr = MC_io.real[ixy]; zi = MC_io.imag[ixy];
		MC_io.real[ixy] = MC_io.real[ixys]; MC_io.imag[ixy] = MC_io.imag[ixys];
		MC_io.real[ixys] = zr; MC_io.imag[ixys] = zi;
	}
}

void f_fft2Shift_MD(sGP &GP, double *&MD_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MD<<<Bnxhnyh, Tnxhnyh>>>(GP, MD_io);
}

void f_fft2Shift_MD(sGP &GP, double *&MD1_io, double *&MD2_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MD<<<Bnxhnyh, Tnxhnyh>>>(GP, MD1_io, MD2_io);
}

void f_fft2Shift_MC(sGP &GP, double2 *&MC_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MC<<<Bnxhnyh, Tnxhnyh>>>(GP, MC_io);
}

void f_fft2Shift_MC(sGP &GP, double2 *&MC1_io, double2 *&MC2_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MC<<<Bnxhnyh, Tnxhnyh>>>(GP, MC1_io, MC2_io);
}

void f_fft2Shift_MC_MD(sGP &GP, double2 *&MC_io, double *&MD_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MC_MD<<<Bnxhnyh, Tnxhnyh>>>(GP, MC_io, MD_io);
}

void f_fft2Shift_MC_MD(sGP &GP, double2 *&MC_io, double *&MD1_io, double *&MD2_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MC_MD<<<Bnxhnyh, Tnxhnyh>>>(GP, MC_io, MD1_io, MD2_io);
}

void f_fft2Shift_MC(sGP &GP, sComplex &MC_io){
	dim3 Bnxhnyh, Tnxhnyh;
	f_get_BTnxhnyh(GP, Bnxhnyh, Tnxhnyh);
	k_fft2Shift_MC<<<Bnxhnyh, Tnxhnyh>>>(GP, MC_io);
}

/***************************************************************************/
/***************************************************************************/

// Sum module square over all the elements
__global__ void k_Sum_MD(sGP GP, const double * __restrict MD_i, double * __restrict MDp_o){ 
	__shared__ double Ms[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double SumT=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			SumT += MD_i[ixy];
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		MDp_o[bid] = Ms[0];
}

// Sum module square over all the elements
__global__ void k_Sum_MC2(sGP GP, const double2 * __restrict MC_i, double * __restrict MDp_o){ 
	__shared__ double Ms[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double x, y;
	double SumT=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			x = MC_i[ixy].x; 
			y = MC_i[ixy].y;
			SumT += x*x+y*y;
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		MDp_o[bid] = Ms[0];
}

// Sum over the detector
__global__ void k_Sum_MD_Det(sGP GP, const double * __restrict MD_i, double gmin2, double gmax2, double * __restrict MDp_o){ 
	__shared__ double Ms[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double gx, gy, g2;
	double SumT=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			if((gmin2<=g2)&&(g2<=gmax2))
				SumT += MD_i[ixy];
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		MDp_o[bid] = Ms[0];
}

// Sum over the detector
__global__ void k_Sum_MD_Det(sGP GP, const double * __restrict MD1_i, const double * __restrict MD2_i, double gmin2, double gmax2, double * __restrict MD1p_o, double * __restrict MD2p_o){ 
	__shared__ double Ms1[thrnxny*thrnxny];
	__shared__ double Ms2[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double gx, gy, g2;
	double SumT1=0, SumT2=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			if((gmin2<=g2)&&(g2<=gmax2)){
				SumT1 += MD1_i[ixy];
				SumT2 += MD2_i[ixy];
			}
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms1[tid] = SumT1;
	Ms2[tid] = SumT2;
	__syncthreads();

	k_reduceBlockDouble256(Ms1, Ms2, tid);

	if (tid==0){
		MD1p_o[bid] = Ms1[0];
		MD2p_o[bid] = Ms2[0];
	}
}

// Sum over the detector
__global__ void k_Sum_MC_Det(sGP GP, const double2 * __restrict MC_i, double gmin2, double gmax2, double * __restrict MDp_o){ 
	__shared__ double Ms[thrnxny*thrnxny];

	int tid = threadIdx.x + threadIdx.y*blockDim.x;
	int bid = blockIdx.x + blockIdx.y*gridDim.x;
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int gridSizex = blockDim.x*gridDim.x;
	int gridSizey = blockDim.y*gridDim.y;
	int ix, iy, ixy;
	double gx, gy, g2;
	double x, y;
	double SumT=0;

	ix = ix0;
	while (ix < GP.nx){
		iy = iy0;
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		while (iy < GP.ny){
			ixy = ix*GP.ny+iy;	
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			if((gmin2<=g2)&&(g2<=gmax2)){
				x = MC_i[ixy].x; 
				y = MC_i[ixy].y;
				SumT += x*x+y*y;
			}
			iy += gridSizey;
		}
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		MDp_o[bid] = Ms[0];
}

// Sum module square over all the elements
__global__ void k_Sum_MDp(int n, double * MDp_i, int iSum, double * Sum_MD_o){ 
	__shared__ double Ms[thrnxy];

	int tid = threadIdx.x;
	int ix = threadIdx.x + blockIdx.x*blockDim.x;
	int gridSizex = blockDim.x*gridDim.x;
	double SumT=0;

	while (ix < n){
		SumT += MDp_i[ix];
		ix += gridSizex;
	}
	Ms[tid] = SumT;
	__syncthreads();

	k_reduceBlockDouble256(Ms, tid);

	if (tid==0)
		Sum_MD_o[iSum] = Ms[0];
}

// Sum module square over all the elements
__global__ void k_Sum_MDp(int n, double * MD1p_i, const double * MD2p_i, int iSum, double * Sum_MD1_o, double * Sum_MD2_o){ 
	__shared__ double Ms1[thrnxy];
	__shared__ double Ms2[thrnxy];

	int tid = threadIdx.x;
	int ix = threadIdx.x + blockIdx.x*blockDim.x;
	int gridSizex = blockDim.x*gridDim.x;
	double SumT1=0, SumT2=0;

	while (ix < n){
		SumT1 += MD1p_i[ix];
		SumT2 += MD2p_i[ix];
		ix += gridSizex;
	}
	Ms1[tid] = SumT1;
	Ms2[tid] = SumT2;
	__syncthreads();

	k_reduceBlockDouble256(Ms1, Ms2, tid);

	if (tid==0){
		Sum_MD1_o[iSum] = Ms1[0];
		Sum_MD2_o[iSum] = Ms2[0];
	}
}

double f_Sum_MD(sGP &GP, double *&MD_i, double *&MDp_i){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	Bnxny.x = MIN(Bnxny.x, 32); Bnxny.y = MIN(Bnxny.y, 32);
	k_Sum_MD<<<Bnxny, Tnxny>>>(GP, MD_i, MDp_i);
	k_Sum_MDp<<<1, thrnxy>>>(Bnxny.x*Bnxny.y, MDp_i, 0, MDp_i);
	double sum;
	cudaMemcpy(&sum, MDp_i, cSizeofRD, cudaMemcpyDeviceToHost);
	return sum;
}

double f_Sum_MC2(sGP &GP, double2 *&MC_i, double *&MDp_i){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	Bnxny.x = MIN(Bnxny.x, 32); Bnxny.y = MIN(Bnxny.y, 32);
	k_Sum_MC2<<<Bnxny, Tnxny>>>(GP, MC_i, MDp_i);
	k_Sum_MDp<<<1, thrnxy>>>(Bnxny.x*Bnxny.y, MDp_i, 0, MDp_i);
	double sum;
	cudaMemcpy(&sum, MDp_i, cSizeofRD, cudaMemcpyDeviceToHost);
	return sum;
}

void f_Sum_MD_Det(sGP &GP, double *&MD_i, double gmin2, double gmax2, double *&MDp_i, int iSum, double *&Sum_MD_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	Bnxny.x = MIN(Bnxny.x, 32); Bnxny.y = MIN(Bnxny.y, 32);
	k_Sum_MD_Det<<<Bnxny, Tnxny>>>(GP, MD_i, gmin2, gmax2, MDp_i);
	k_Sum_MDp<<<1, thrnxy>>>(Bnxny.x*Bnxny.y, MDp_i, iSum, Sum_MD_o);
}

void f_Sum_MD_Det(sGP &GP, double *&MD1_i, double *&MD2_i, double gmin2, double gmax2, double *&MD1p_i, double *&MD2p_i, int iSum, double *&Sum_MD1_o, double *&Sum_MD2_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	Bnxny.x = MIN(Bnxny.x, 32); Bnxny.y = MIN(Bnxny.y, 32);	
	k_Sum_MD_Det<<<Bnxny, Tnxny>>>(GP, MD1_i, MD2_i, gmin2, gmax2, MD1p_i, MD2p_i);
	k_Sum_MDp<<<1, thrnxy>>>(Bnxny.x*Bnxny.y, MD1p_i, MD2p_i, iSum, Sum_MD1_o, Sum_MD2_o);
}

void f_Sum_MC_Det(sGP &GP, double2 *&MC_i, double gmin2, double gmax2, double *&MDp_i, int iSum, double *&Sum_MD_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	Bnxny.x = MIN(Bnxny.x, 32); Bnxny.y = MIN(Bnxny.y, 32);
	k_Sum_MC_Det<<<Bnxny, Tnxny>>>(GP, MC_i, gmin2, gmax2, MDp_i);
	k_Sum_MDp<<<1, thrnxy>>>(Bnxny.x*Bnxny.y, MDp_i, iSum, Sum_MD_o);
}

/***************************************************************************/
/***************************************************************************/

// Phase
__global__ void k_Phase(sGP GP, double gxu, double gyu, sACD ExpR_x_o, sACD ExpR_y_o){
	int ix = threadIdx.x + blockIdx.x*blockDim.x, iy = ix;
	double theta, R, x, y;

	if (ix < GP.nx){
		R = IsRS(ix,GP.nxh)*GP.dRx;
		theta = c2Pi*R*gxu;
		sincos(theta, &y, &x);
		ExpR_x_o.x[ix] = x;
		ExpR_x_o.y[ix] = y;
	}

	if (iy < GP.ny){
		R = IsRS(iy,GP.nyh)*GP.dRy;
		theta = c2Pi*R*gxu;
		sincos(theta, &y, &x);
		ExpR_y_o.x[iy] = x;
		ExpR_y_o.y[iy] = y;
	}
}

// Phase multiplication
__global__ void k_PhaseMul(sGP GP, const sACD ExpR_x_i, const sACD ExpR_y_i, double2 * __restrict Psi_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double z1r = ExpR_x_i.x[ix], z1i = ExpR_x_i.y[ix];
		double z2r = ExpR_y_i.x[iy], z2i = ExpR_y_i.y[iy];
		double z3r = z1r*z2r-z1i*z2i, z3i = z1i*z2r+z1r*z2i;
		z2r = Psi_io[ixy].x;		z2i = Psi_io[ixy].y; 
		Psi_io[ixy].x = z1r = z3r*z2r-z3i*z2i; 
		Psi_io[ixy].y = z1i = z3i*z2r+z3r*z2i; 
	}
}

// Anti-Aliasing, scale with cut-off (2/3)gmax
__global__ void k_BandwidthLimit2D(sGP GP, double2 * __restrict M_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if (g2 < GP.gmaxl2){
			M_io[ixy].x *= GP.inxy; 
			M_io[ixy].y *= GP.inxy; 
		}else{
 			M_io[ixy].x = 0.0; 
			M_io[ixy].y = 0.0; 
		}
	}
}

// Element by element multiplication
__global__ void k_Transmit(sGP GP, const double2 * __restrict Trans_i, double2 * __restrict Psi_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double z1r = Trans_i[ixy].x, z1i = Trans_i[ixy].y;
		double z2r = Psi_io[ixy].x, z2i = Psi_io[ixy].y;
		double z3r = z1r*z2r-z1i*z2i, z3i = z1i*z2r+z1r*z2i;
		Psi_io[ixy].x = z3r;
		Psi_io[ixy].y = z3i;
	}
}

// Build propagator function
__global__ void k_Propagator(sGP GP, double gxu, double gyu, double scale, sACD Prop_x_o, sACD Prop_y_o){
	int ix = threadIdx.x + blockIdx.x*blockDim.x, iy = ix;
	double theta, gx, gy, x, y;

	if (ix < GP.nx){
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		theta = scale*(gx+gxu)*(gx+gxu);
		sincos(theta, &y, &x);
		Prop_x_o.x[ix] = x;
		Prop_x_o.y[ix] = y;
	}

	if (iy < GP.ny){
		gy = IsFS(iy,GP.nyh)*GP.dgy;
		theta = scale*(gy+gyu)*(gy+gyu);
		sincos(theta, &y, &x);
		Prop_y_o.x[iy] = x;
		Prop_y_o.y[iy] = y;
	}
}

// Propagate, scale with cut-off (2/3)gmax
__global__ void k_Propagate(sGP GP, const sACD Prop_x_i, const sACD Prop_y_i, double2 * __restrict Psi_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if (!GP.BWL||(g2 < GP.gmaxl2)){
			double z1r = Prop_x_i.x[ix], z1i = Prop_x_i.y[ix];
			double z2r = Prop_y_i.x[iy], z2i = Prop_y_i.y[iy];
			double z3r = z1r*z2r-z1i*z2i, z3i = z1i*z2r+z1r*z2i;
			z2r = Psi_io[ixy].x;		z2i = Psi_io[ixy].y; 
			Psi_io[ixy].x = z1r = (z3r*z2r-z3i*z2i)*GP.inxy;
			Psi_io[ixy].y = z1i = (z3i*z2r+z3r*z2i)*GP.inxy; 
		}else{
 			Psi_io[ixy].x = 0.0; 
			Psi_io[ixy].y = 0.0; 
		}
	}
}

// Probe in Fourier space
__global__ void k_Probe_FS(sGP GP, sLens Lens, double x, double y, double2 * __restrict fPsi_o){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){	
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;
		if ((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2)){
			double chi = x*gx + y*gy + g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
			if ((Lens.m!=0)||(Lens.cmfa2!=0)||(Lens.cmfa3!=0)){
				double g = sqrt(g2);
				double phi = (g==0)?0:acos(gy/g);
				phi = (gy<0)?cPi+phi:phi;
				chi += Lens.m*phi + Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
			}
			double chix, chiy;
			sincos(chi, &chiy, &chix);		
			fPsi_o[ixy].x = chix; 
			fPsi_o[ixy].y = chiy;	
		}else{
 			fPsi_o[ixy].x = 0.0; 
			fPsi_o[ixy].y = 0.0; 
		}
	}
}

// Apply Coherent transfer function
__global__ void k_Apply_CTF(sGP GP, sLens Lens, double gxu, double gyu, double2 * __restrict fPsi_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;
		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){
			g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
			double chi = g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
			if ((Lens.cmfa2!=0)||(Lens.cmfa3!=0)){
				double g = sqrt(g2);
				double phi = (g==0)?0:acos(gy/g);
				phi = (gy<0)?cPi+phi:phi;
				chi += Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
			}
			double Psix = fPsi_io[ixy].x, Psiy = fPsi_io[ixy].y; 
			double chix, chiy;
			sincos(chi, &chiy, &chix);			
			fPsi_io[ixy].x = gx = chix*Psix-chiy*Psiy; 
			fPsi_io[ixy].y = gy = chix*Psiy+chiy*Psix;
		}
		else{
 			fPsi_io[ixy].x = 0.0; 
			fPsi_io[ixy].y = 0.0; 
		}
	}
}

// Apply Coherent transfer function
__global__ void k_Apply_CTF(sGP GP, sLens Lens, double gxu, double gyu, const double2 * __restrict fPsi_i, double2 * __restrict fPsi_o){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = ((ix<GP.nxh)?ix:(ix-GP.nx))*GP.dgx;
		double gy = ((iy<GP.nyh)?iy:(iy-GP.ny))*GP.dgy;
		double g2 = gx*gx + gy*gy;
		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){
			g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
			double chi = g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
			if ((Lens.cmfa2!=0)||(Lens.cmfa3!=0)){
				double g = sqrt(g2);
				double phi = (g==0)?0:acos(gy/g);
				phi = (gy<0)?cPi+phi:phi;
				chi += Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
			}
			double Psix = fPsi_i[ixy].x, Psiy = fPsi_i[ixy].y; 
			double chix, chiy;
			sincos(chi, &chiy, &chix);			
			fPsi_o[ixy].x = gx = chix*Psix-chiy*Psiy; 
			fPsi_o[ixy].y = gy = chix*Psiy+chiy*Psix;
		}
		else{
 			fPsi_o[ixy].x = 0.0; 
			fPsi_o[ixy].y = 0.0; 
		}
	}
}

// Partially coherent transfer function, linear image model and weak phase object
__global__ void k_Apply_PCTF(sGP GP, sLens Lens, double2 * __restrict fPsi_io){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){			
			double chi = g2*(Lens.cCs3*g2+Lens.cf);
			double c = cPi*Lens.beta*Lens.sf;
			double u = 1.0 + c*c*g2;

			c = cPi*Lens.sf*Lens.lambda*g2;
			double sie = 0.25*c*c;
			c = cPi*Lens.beta*(Lens.Cs3*Lens.lambda2*g2-Lens.f);
			double tie = c*c*g2;
			double sti = exp(-(sie+tie)/u);

			double Psix = fPsi_io[ixy].x, Psiy = fPsi_io[ixy].y; 
			double chix, chiy;
			sincos(chi, &chiy, &chix);		 	

			fPsi_io[ixy].x = gx = sti*(chix*Psix-chiy*Psiy); 
			fPsi_io[ixy].y = gy = sti*(chix*Psiy+chiy*Psix);
		}else{
 			fPsi_io[ixy].x = 0.0; 
			fPsi_io[ixy].y = 0.0; 
		}
	}
}

// Partially coherent transfer function, linear image model and weak phase object
__global__ void k_Apply_PCTF(sGP GP, sLens Lens, double2 * __restrict fPsi_i, double2 * __restrict fPsi_o){
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < GP.nx)&&(iy < GP.ny)){
		int ixy = ix*GP.ny+iy;
		double gx = IsFS(ix,GP.nxh)*GP.dgx;
		double gy = IsFS(iy,GP.nyh)*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){			
			double chi = g2*(Lens.cCs3*g2+Lens.cf);
			double c = cPi*Lens.beta*Lens.sf;
			double u = 1.0 + c*c*g2;

			c = cPi*Lens.sf*Lens.lambda*g2;
			double sie = 0.25*c*c;
			c = cPi*Lens.beta*(Lens.Cs3*Lens.lambda2*g2-Lens.f);
			double tie = c*c*g2;
			double sti = exp(-(sie+tie)/u);

			double Psix = fPsi_i[ixy].x, Psiy = fPsi_i[ixy].y; 
			double chix, chiy;
			sincos(chi, &chiy, &chix);		 	

			fPsi_o[ixy].x = gx = sti*(chix*Psix-chiy*Psiy); 
			fPsi_o[ixy].y = gy = sti*(chix*Psiy+chiy*Psix);
		}else{
 			fPsi_o[ixy].x = 0.0; 
			fPsi_o[ixy].y = 0.0; 
		}
	}
}

void f_PhaseMul(sGP &GP, double gxu, double gyu, sACD &ExpRg_x_o, sACD &ExpRg_y_o, double2 *&Psi_io){
	dim3 Bmnxny, Tmnxny;
	f_get_BTmnxny(GP, Bmnxny, Tmnxny);	
	// exp(2*1i*pi(rx*gxu+ry*gyu)
	k_Phase<<<Bmnxny, Tmnxny>>>(GP, gxu, gyu, ExpRg_x_o, ExpRg_y_o);
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	// Phase multiplication
	k_PhaseMul<<<Bnxny, Tnxny>>>(GP, ExpRg_x_o, ExpRg_y_o, Psi_io);
}

void f_BandwidthLimit2D(cufftHandle &PlanPsi, sGP &GP, double2 *&MC_io){
	if(!GP.BWL) return;

	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);

	// Forward fft2
	cufftExecZ2Z(PlanPsi, MC_io, MC_io, CUFFT_FORWARD);
	// AntiAliasing, scale and bandlimited
	k_BandwidthLimit2D<<<Bnxny, Tnxny>>>(GP, MC_io);		
	// Backward fft2
	cufftExecZ2Z(PlanPsi, MC_io, MC_io, CUFFT_INVERSE);
}

void f_Transmit(sGP &GP, double2 *&Trans_i, double2 *&Psi_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);	
	k_Transmit<<<Bnxny, Tnxny>>>(GP, Trans_i, Psi_io);
}

void f_Propagate(cufftHandle &PlanPsi, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sACD &Prop_x_o, sACD &Prop_y_o, double2 *&Psi_io){
	// Forward fft2
	cufftExecZ2Z(PlanPsi, Psi_io, Psi_io, CUFFT_FORWARD);

	dim3 Bmnxny, Tmnxny;
	f_get_BTmnxny(GP, Bmnxny, Tmnxny);	
	// Forward propagator
	k_Propagator<<<Bmnxny, Tmnxny>>>(GP, gxu, gyu, -cPi*lambda*z, Prop_x_o, Prop_y_o);
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	// Propagate, scale and bandlimited the wave function
	k_Propagate<<<Bnxny, Tnxny>>>(GP, Prop_x_o, Prop_y_o, Psi_io);

	// Backward fft2
	if(Space==eSReal)
		cufftExecZ2Z(PlanPsi, Psi_io, Psi_io, CUFFT_INVERSE);
}

void f_Probe_FS(sGP &GP, sLens Lens, double x, double y, double2 *&fPsi_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Probe_FS<<<Bnxny, Tnxny>>>(GP, Lens, x, y, fPsi_o);
}

void f_Apply_CTF(sGP &GP, sLens &Lens, double gxu, double gyu, double2 *&fPsi_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Apply_CTF<<<Bnxny, Tnxny>>>(GP, Lens, gxu, gyu, fPsi_io);
}

void f_Apply_CTF(sGP &GP, sLens &Lens, double gxu, double gyu, double2 *&fPsi_i, double2 *&fPsi_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Apply_CTF<<<Bnxny, Tnxny>>>(GP, Lens, gxu, gyu, fPsi_i, fPsi_o);
}

void f_Apply_PCTF(sGP &GP, sLens &Lens, double2 *&fPsi_io){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Apply_PCTF<<<Bnxny, Tnxny>>>(GP, Lens, fPsi_io);
}

void f_Apply_PCTF(sGP &GP, sLens &Lens, double2 *&fPsi_i, double2 *&fPsi_o){
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);
	k_Apply_PCTF<<<Bnxny, Tnxny>>>(GP, Lens, fPsi_i, fPsi_o);
}

/***************************************************************************/
/***************************************************************************/

// From Device To Host
void f_Copy_MCd(sGP &GP, double2 *&MC_d_i, double *&MCr_d_i, double *&MCi_d_i, sComplex &MC_h_o){	
	sComplex MC_d;
	MC_d.real = MCr_d_i; MC_d.imag = MCi_d_i;
	// Get Real and Imaginary part to complex matrix
	f_Get_MC(GP, MC_d_i, MC_d);
	// Copy real part of the wave function to the host
	cudaMemcpy(MC_h_o.real, MC_d.real, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	// Copy imaginary part of the wave function to the host
	cudaMemcpy(MC_h_o.imag, MC_d.imag, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void f_Copy_MCd_MDd(sGP &GP, double2 *&MC_d_i, double *&MD_d_i, double *&MCr_d_i, double *&MCi_d_i, sComplex &MC_h_o, double *&MD_h_o){
	f_Copy_MCd(GP, MC_d_i, MCr_d_i, MCi_d_i, MC_h_o);
	// Copy wave function squared to the host
	cudaMemcpy(MD_h_o, MD_d_i, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void f_Copy_MCd_MDd(sGP &GP, double2 *&MC_d_i, double *&MD1_d_i, double *&MD2_d_i, double *&MCr_d_i, double *&MCi_d_i, sComplex &MC_h_o, double *&MD1_h_o, double *&MD2_h_o){	
	f_Copy_MCd(GP, MC_d_i, MCr_d_i, MCi_d_i, MC_h_o);
	// Copy wave function squared to the host
	cudaMemcpy(MD1_h_o, MD1_d_i, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(MD2_h_o, MD2_d_i, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

/***************************************************************************/
/***************************************************************************/

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

/***************************************************************************/
/***************************************************************************/

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

/***************************************************************************/
/***************************************************************************/

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

/***************************************************************************/
/***************************************************************************/

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

/***************************************************************************/
/***************************************************************************/

void f_sACD_cudaFree(sACD &ACD){
	cudaFreen(ACD.x);
	cudaFreen(ACD.y);
}

void f_sACD_cudaInit(sACD &ACD){
	ACD.x = 0;
	ACD.y = 0;
}

void f_sACD_cudaMalloc(int nACD, sACD &ACD){
	if(nACD<=0) return;

	cudaMalloc((void**)&(ACD.x), nACD*cSizeofRD);
	cudaMalloc((void**)&(ACD.y), nACD*cSizeofRD);
}

/***************************************************************************/
/***************************************************************************/

void f_GPU_Sync_CPU(int &iSynCPU, int &cSynCPU){
	if(cSynCPU<0) cSynCPU = ccSynCPU;

	iSynCPU++;
	if(iSynCPU%cSynCPU==0)
		cudaDeviceSynchronize();
}