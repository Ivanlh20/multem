/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef DEVICE_FUNCTIONS_H
#define DEVICE_FUNCTIONS_H

#include <type_traits>
#include <algorithm>

#include "math.cuh"
#include "types.hpp"

#include "host_device_functions.cuh"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_atomic_functions.h>
#include <cufft.h>

#include <thrust/swap.h>
#include <thrust/complex.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/for_each.h>
#include <thrust/fill.h>

namespace multem
{
	namespace device_detail
	{
		template<class TGrid>
		GridBT get_grid_nxny(const TGrid &grid, const dim3 Blk_max = dim3(0, 0, 0))
		{
			GridBT gridBT;
			gridBT.Blk = dim3((grid.ny+c_thrnxny-1)/c_thrnxny, (grid.nx+c_thrnxny-1)/c_thrnxny);
			gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

			if(Blk_max.x != 0)
			{
				gridBT.Blk.x = min(Blk_max.x, gridBT.Blk.x);
			}
			if(Blk_max.y != 0)
			{
				gridBT.Blk.y = min(Blk_max.y, gridBT.Blk.y);
			}
			if(Blk_max.z != 0)
			{
				gridBT.Blk.z = min(Blk_max.z, gridBT.Blk.z);
			}

			return gridBT;
		}

		template<class T>
		__device__ __forceinline__ 
		void atomicAdd(T *address, T val)
		{
		}

		template<>
		__device__ __forceinline__ 
		void atomicAdd<float>(float *address, float val)
		{
			::atomicAdd(address, val);
		}

		template<>
		__device__ __forceinline__ 
		void atomicAdd<double>(double *address, double val)
		{
			unsigned long long int* address_as_ull = (unsigned long long int*)address;
			unsigned long long int old = *address_as_ull, assumed;

			do {
				assumed = old;
				old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
			} while (assumed != old);
		}

		__device__ __forceinline__ 
			double __shfl_down(double var, unsigned int srcLane, int width=32)
		{
			int2 a = *reinterpret_cast<int2*>(&var);
			a.x = __shfl_down(a.x, srcLane, width);
			a.y = __shfl_down(a.y, srcLane, width);
			return *reinterpret_cast<double*>(&a);
		}

		// Linear projected potential: V and zV
		template<ePotential_Type potential_type, class T> 
		__global__ void linear_Vz(rQ1<T> qz, Atom_Vp<T> atom_Vp)
		{	
			__shared__ T V0s[c_nqz];
			__shared__ T dV0s[c_nqz];

			int ix = threadIdx.x, iR = blockIdx.x;

			T x = qz.x[ix];
			T w = qz.w[ix];
			T R2 = atom_Vp.R2[iR];
			T V, dVir;

			T a = (atom_Vp.split)?(-atom_Vp.z0h):(atom_Vp.zeh-atom_Vp.z0h);
			T b = (atom_Vp.split)?(atom_Vp.z0h):(atom_Vp.zeh+atom_Vp.z0h);
			T z = a*x + b;
			T r = sqrt(z*z + R2);
			Vr_dVrir<potential_type, T>(r, atom_Vp.cl, atom_Vp.cnl, a*w, V, dVir);

			V0s[ix] = V; 
			dV0s[ix] = dVir;

			if (atom_Vp.split)
			{
				a = atom_Vp.zeh;
				b = atom_Vp.zeh;
				z = a*x + b; 
				r = sqrt(z*z + R2);
				Vr_dVrir<potential_type, T>(r, atom_Vp.cl, atom_Vp.cnl, a*w, V, dVir);

				V = V0s[ix] += V; 
				dVir = dV0s[ix] += dVir;
			}

			__syncthreads();

			if(ix < 64)
			{
				V0s[ix] = V = V + V0s[ix + 64];
				dV0s[ix] = dVir = dVir + dV0s[ix + 64];
			}
			__syncthreads();

			if(ix < 32)
			{
				V0s[ix] = V = V + V0s[ix + 32];
				dV0s[ix] = dVir = dVir + dV0s[ix + 32];
			}
			__syncthreads();

			if(ix < 16)
			{
				V0s[ix] = V = V + V0s[ix + 16];
				dV0s[ix] = dVir = dVir + dV0s[ix + 16];
			}
			__syncthreads();

			if(ix < 8)
			{
				V0s[ix] = V = V + V0s[ix + 8];
				dV0s[ix] = dVir = dVir + dV0s[ix + 8];
			}
			__syncthreads();

			if(ix < 4)
			{
				V0s[ix] = V = V + V0s[ix + 4];
				dV0s[ix] = dVir = dVir + dV0s[ix + 4];
			}
			__syncthreads();

			if(ix < 2)
			{
				V0s[ix] = V = V + V0s[ix + 2];
				dV0s[ix] = dVir = dVir + dV0s[ix + 2];
			}
			__syncthreads();

			if(ix < 1)
			{
				V0s[ix] = V = V + V0s[ix + 1];
				dV0s[ix] = dVir = dVir + dV0s[ix + 1];
			}
			__syncthreads();

			if(ix == 0 )
			{
				atom_Vp.c0[iR] = V; 		// V0
				atom_Vp.c1[iR] = 0.5*dVir; 	// dR2V0
			}
		}

		// Get Local interpolation coefficients
		template<class T> 
		__global__ void cubic_poly_coef(Atom_Vp<T> atom_Vp)
		{
			int iR = threadIdx.x;

			if (iR < c_nR-1)
			{
				T dR2 = 1.0/(atom_Vp.R2[iR+1]-atom_Vp.R2[iR]);
				T V = atom_Vp.c0[iR]; 
				T Vn = atom_Vp.c0[iR+1];
				T dV = atom_Vp.c1[iR]; 
				T dVn = atom_Vp.c1[iR+1];
				T m = (Vn-V)*dR2; 
				T n = dV+dVn;
				atom_Vp.c0[iR] = V-atom_Vp.c0[c_nR-1];
				atom_Vp.c2[iR] = (3.0*m-n-dV)*dR2;
				atom_Vp.c3[iR] = (n-2.0*m)*dR2*dR2;
			}
		}

		// Cubic polynomial evaluation
		template<class T> 
		__global__ void eval_cubic_poly(Grid<T> grid, Atom_Vp<T> atom_Vp, rVector<T> V0g)
		{
			int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;

			if ((ix0 < atom_Vp.ixn)&&(iy0 < atom_Vp.iyn))
			{
				int ix = ix0 + atom_Vp.ix0;
				int iy = iy0 + atom_Vp.iy0;

				T Rx = grid.Rx(ix) - atom_Vp.x;
				T Ry = grid.Ry(iy) - atom_Vp.y;
				T R2 = Rx*Rx + Ry*Ry;
				if (R2 < atom_Vp.R_max2)
				{
					R2 = max(R2, atom_Vp.R_min2);

					ix -= static_cast<int>(floor(grid.Rx(ix)/grid.lx))*grid.nx;
					iy -= static_cast<int>(floor(grid.Ry(iy)/grid.ly))*grid.ny;

					ix = grid.iRx_shift(ix);
					iy = grid.iRy_shift(iy);
					int ixy = grid.ind_col(ix, iy);

					ix = host_device_detail::unrolledBinarySearch_c_nR<T>(R2, atom_Vp.R2);

					T dx = R2 - atom_Vp.R2[ix]; 
					T dx2 = dx*dx;
					T V = atom_Vp.occ*(atom_Vp.c0[ix] + atom_Vp.c1[ix]*dx + atom_Vp.c2[ix]*dx2 + atom_Vp.c3[ix]*dx2*dx);
						
					atomicAdd<T>(&(V0g.V[ixy]), V);
				}
			}
		}
		
		// Shift matrix respect to (nxh, nyh)
		template<class TGrid, class T>
		__global__ void fft2_shift(TGrid grid, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nxh)&&(iy < grid.nyh))
			{
				int ixy = grid.ind_col(ix, iy); 
				int ixy_shift = grid.ind_col(grid.nxh+ix, grid.nyh+iy);
				thrust::swap(M_io.V[ixy], M_io.V[ixy_shift]);

				ixy = grid.ind_col(ix, grid.nyh+iy); 
				ixy_shift = grid.ind_col(grid.nxh+ix, iy);
				thrust::swap(M_io.V[ixy], M_io.V[ixy_shift]);
			}
		}

		// sum over the detector
		template<class TGrid, class T>
		__global__ void sum_over_Det(TGrid grid, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, rVector<T> M_i, rVector<T> Mp_o)
		{ 
			using value_type_r = Value_type<TGrid>;

			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;

			T sum = 0;

			int ix = ix0;
			while (ix < grid.nx)
			{
				int iy = iy0;
				while (iy < grid.ny)
				{
					int ixy = grid.ind_col(ix, iy); 	
					value_type_r g2 = grid.g2_shift(ix, iy);
					if((g2_min <= g2)&&(g2 <= g2_max))
					{
						sum += M_i.V[ixy];
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}
			Mshare[tid] = sum;

			__syncthreads();

			if(tid < 128)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 128];
			}
			__syncthreads();

			if(tid < 64)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 64];
			}
			__syncthreads();

			if(tid < 32)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 32];
			}
			__syncthreads();

			if(tid < 16)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 16];
			}
			__syncthreads();

			if(tid < 8)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 8];
			}
			__syncthreads();

			if(tid < 4)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 4];
			}
			__syncthreads();

			if(tid < 2)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 2];
			}
			__syncthreads();

			if(tid < 1)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 1];
			}
			__syncthreads();

			if(tid == 0 )
			{
				Mp_o.V[bid] = sum;
			}
		}

		// sum over the detector
		template<class TGrid, class T>
		__global__ void sum_square_over_Det(TGrid grid, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, rVector<T> M_i, rVector<Value_type<TGrid>> Mp_o)
		{ 
			using value_type_r = Value_type<TGrid>;

			__shared__ value_type_r Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;

			value_type_r sum = 0;

			int ix = ix0;
			while (ix < grid.nx)
			{
				int iy = iy0;
				while (iy < grid.ny)
				{
					int ixy = grid.ind_col(ix, iy); 	
					value_type_r g2 = grid.g2_shift(ix, iy);
					if((g2_min <= g2)&&(g2 <= g2_max))
					{
						sum += thrust::norm(M_i.V[ixy]);
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}
			Mshare[tid] = sum;

			__syncthreads();

			if(tid < 128)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 128];
			}
			__syncthreads();

			if(tid < 64)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 64];
			}
			__syncthreads();

			if(tid < 32)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 32];
			}
			__syncthreads();

			if(tid < 16)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 16];
			}
			__syncthreads();

			if(tid < 8)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 8];
			}
			__syncthreads();

			if(tid < 4)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 4];
			}
			__syncthreads();

			if(tid < 2)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 2];
			}
			__syncthreads();

			if(tid < 1)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 1];
			}
			__syncthreads();

			if(tid == 0 )
			{
				Mp_o.V[bid] = sum;
			}
		}

		// Anti-Aliasing, scale with cut-off (2/3)g_max
		template<class TGrid, class T>
		__global__ void bandwidth_limit(TGrid grid, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, T w, rVector<T> M_io)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy); 
				value_type_r g2 = grid.g2_shift(ix, iy);

				if((g2_min <= g2)&&(g2 <= g2_max))
				{
					M_io.V[ixy] *= w;
				}
				else
				{
 					M_io.V[ixy] = 0;
				}
			}
		}

		// Phase multiplication
		template<class TGrid, class T>
		__global__ void phase_mul(TGrid grid, rVector<T> V_x_i, 
		rVector<T> V_y_i, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				psi_o.V[ixy] = psi_i.V[ixy]*V_x_i.V[ix]*V_y_i.V[iy]; 
			}
		}

		// Propagate, scale with cut-off (2/3)g_max
		template<class TGrid, class T>
		__global__ void propagator_mul(TGrid grid, rVector<T> V_x_i, 
		rVector<T> V_y_i, rVector<T> psi_i, rVector<T> psi_o)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r g2 = grid.g2_shift(ix, iy);

				if((!grid.bwl)||(g2 < grid.gl2_max))
				{
					psi_o.V[ixy] = static_cast<T>(grid.inxy)*psi_i.V[ixy]*V_x_i.V[ix]*V_y_i.V[iy];
				}
				else
				{
 					psi_o.V[ixy] = static_cast<T>(0);
				}
			}
		}

		// Probe in Fourier space
		template<class TGrid, class T>
		__global__ void probe(TGrid grid, Lens<Value_type<TGrid>> lens, 
		Value_type<TGrid> x, Value_type<TGrid> y, rVector<T> fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{	
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;

				if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
				{
					value_type_r chi = x*gx + y*gy + g2*(lens.cCs5*g2*g2+lens.cCs3*g2+lens.cf);
					if(nonZero(lens.m)||nonZero(lens.cmfa2)||nonZero(lens.cmfa3))
					{
						value_type_r g = sqrt(g2);
						value_type_r phi = atan2(gy, gx);
						chi += lens.m*phi + lens.cmfa2*g2*sin(2*(phi-lens.afa2)) + lens.cmfa3*g*g2*sin(3*(phi-lens.afa3)); 			
					}	
					fPsi_o.V[ixy] = thrust::euler(chi); 
				}
				else
				{
 					fPsi_o.V[ixy] = static_cast<T>(0);
				}
			}
		}

		// Apply Coherent transfer function
		template<class TGrid, class T>
		__global__ void apply_CTF(TGrid grid, Lens<Value_type<TGrid>> lens, 
		Value_type<TGrid> gxu, Value_type<TGrid> gyu, rVector<T> fPsi_i, rVector<T> fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{	
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;

				if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
				{
					g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
					value_type_r chi = g2*(lens.cCs5*g2*g2+lens.cCs3*g2+lens.cf);
					if(nonZero(lens.cmfa2)||nonZero(lens.cmfa3))
					{
						value_type_r g = sqrt(g2);
						value_type_r phi = atan2(gy, gx);
						chi += lens.cmfa2*g2*sin(2*(phi-lens.afa2)) + lens.cmfa3*g*g2*sin(3*(phi-lens.afa3)); 			
					}
					fPsi_o.V[ixy] = fPsi_i.V[ixy]*thrust::euler(chi);
				}
				else
				{
 					fPsi_o.V[ixy] = static_cast<T>(0);
				}
			}
		}

		// Partially coherent transfer function, linear image model and weak phase_components object
		template<class TGrid, class T>
		__global__ void apply_PCTF(TGrid grid, Lens<Value_type<TGrid>> lens, rVector<T> fPsi_i, rVector<T> fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;
			const value_type_r c_Pi = 3.141592653589793238463;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r g2 = grid.g2_shift(ix, iy);

				if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
				{			
					value_type_r chi = g2*(lens.cCs3*g2+lens.cf);
					value_type_r c = c_Pi*lens.beta*lens.sf;
					value_type_r u = 1.0 + c*c*g2;

					c = c_Pi*lens.sf*lens.lambda*g2;
					value_type_r sie = 0.25*c*c;
					c = c_Pi*lens.beta*(lens.Cs3*lens.lambda2*g2-lens.f);
					value_type_r tie = c*c*g2;
					value_type_r sti = exp(-(sie+tie)/u);

					fPsi_o.V[ixy] = fPsi_i.V[ixy]*thrust::polar(sti, chi);
				}
				else
				{
					fPsi_o.V[ixy] = static_cast<T>(0);
				}
			}
		}

		// sum over the detector
		template<class T>
		__global__ void norm_factor_lorentz(Grid<T> grid, T gc2, T ge2, rVector<T> Mp_o)
		{ 
			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;

			T sum = 0;

			int ix = ix0;
			while (ix < grid.nx)
			{
				int iy = iy0;
				while (iy < grid.ny)
				{
					T g2 = grid.g2_shift(ix, iy);
					if(g2 < gc2)
					{
						sum += 1.0/(g2 + ge2);
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}
			Mshare[tid] = sum;

			__syncthreads();

			if(tid < 128)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 128];
			}
			__syncthreads();

			if(tid < 64)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 64];
			}
			__syncthreads();

			if(tid < 32)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 32];
			}
			__syncthreads();

			if(tid < 16)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 16];
			}
			__syncthreads();

			if(tid < 8)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 8];
			}
			__syncthreads();

			if(tid < 4)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 4];
			}
			__syncthreads();

			if(tid < 2)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 2];
			}
			__syncthreads();

			if(tid < 1)
			{
				Mshare[tid] = sum = sum + Mshare[tid + 1];
			}
			__syncthreads();

			if(tid == 0 )
			{
				Mp_o.V[bid] = sum;
			}
		}

		template<class TGrid>
		Value_type<TGrid> Lorentz_factor(TGrid &grid, EELS<Value_type<TGrid>> &eels)
		{
			device_vector<Value_type<TGrid>> sum_v(c_thrnxny*c_thrnxny);

			auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

			norm_factor_lorentz<TGrid::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels.gc2, eels.ge2, sum_v);

			return sqrt(eels.occ)/multem::sum(grid, sum_v);
		}

		template<class TGrid, class T>
		__global__ void kernel_xyz(TGrid grid, EELS<Value_type<TGrid>> eels,
		rVector<T> k_x, rVector<T> k_y, rVector<T> k_z)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = eels.factor/(g2 + eels.ge2);
					k_x.V[ixy] = static_cast<T>(gx*lorentz)*pos;
					k_y.V[ixy] = static_cast<T>(gy*lorentz)*pos;
					k_z.V[ixy] = static_cast<T>(eels.ge*lorentz)*pos;
				}
				else
				{
					k_x.V[ixy] = static_cast<T>(0);
					k_y.V[ixy] = static_cast<T>(0);
					k_z.V[ixy] = static_cast<T>(0);
				}
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_x(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_x)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = eels.factor/(g2 + eels.ge2);
					k_x.V[ixy] = static_cast<T>(gx*lorentz)*pos;
				}
				else
				{
					k_x.V[ixy] = static_cast<T>(0);
				}
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_y(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_y)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = eels.factor/(g2 + eels.ge2);
					k_y.V[ixy] = static_cast<T>(gy*lorentz)*pos;
				}
				else
				{
					k_y.V[ixy] = static_cast<T>(0);
				}
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_z(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_z)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = eels.factor/(g2 + eels.ge2);
					k_z.V[ixy] = static_cast<T>(eels.ge*lorentz)*pos;
				}
				else
				{
					k_z.V[ixy] = static_cast<T>(0);
				}
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_mn1(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_mn1)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					const value_type_r c_i2i2 = 0.70710678118654746; 

					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
					T k_x = T(gx*lorentz, 0);
					T k_y = T(0, gy*lorentz);
					k_mn1.V[ixy] = (k_x - k_y)*pos;
				}
				else
				{
					k_mn1.V[ixy] = static_cast<T>(0);
				}
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_mp1(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_mp1)
		{
			using value_type_r = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				int ixy = grid.ind_col(ix, iy);
				value_type_r gx = grid.gx_shift(ix);
				value_type_r gy = grid.gy_shift(iy);
				value_type_r g2 = gx*gx + gy*gy;
				
				if(g2 < eels.gc2)
				{
					const value_type_r c_i2i2 = 0.70710678118654746; 

					T pos = thrust::euler(eels.x*gx + eels.y*gy);
					value_type_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
					T k_x = T(gx*lorentz, 0);
					T k_y = T(0, gy*lorentz);
					k_mp1.V[ixy] = (k_x + k_y)*pos;
				}
				else
				{
					k_mp1.V[ixy] = static_cast<T>(0);
				}
			}
		}

	} // namespace device_detail

	template<class TQ1>
	enable_if_Device<TQ1, void>
	get_cubic_poly_coef_Vz(ePotential_Type potential_type, TQ1 &qz,
	Stream<Value_type<TQ1>, e_Device> &stream, Vector<Atom_Vp<Value_type<TQ1>>, e_Host> &atom_Vp)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					device_detail::linear_Vz<ePT_Doyle_0_4, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
				case ePT_Peng_0_4:
					device_detail::linear_Vz<ePT_Peng_0_4, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
				case ePT_Peng_0_12:
					device_detail::linear_Vz<ePT_Peng_0_12, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
				case ePT_Kirkland_0_12:
					device_detail::linear_Vz<ePT_Kirkland_0_12, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
				case ePT_Weickenmeier_0_12:
					device_detail::linear_Vz<ePT_Weickenmeier_0_12, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
				case ePT_Lobato_0_12:
					device_detail::linear_Vz<ePT_Lobato_0_12, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
					break;
			}
			device_detail::cubic_poly_coef<typename TQ1::value_type><<<dim3(1), dim3(c_nR), 0, stream[istream]>>>(atom_Vp[istream]);
		}
	}

	template<class TGrid, class TVector_r>
	enable_if_Device<TVector_r, void>
	eval_cubic_poly(TGrid &grid, Stream<Value_type<TGrid>, e_Device> &stream, 
	Vector<Atom_Vp<Value_type<TGrid>>, e_Host> &atom_Vp, TVector_r &V0)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			GridBT gridBT = atom_Vp[istream].get_eval_cubic_poly_gridBT();

			device_detail::eval_cubic_poly<typename TVector_r::value_type><<<gridBT.Blk, gridBT.Thr, 0, stream[istream]>>>(grid, atom_Vp[istream], V0);
		}
	}

	template<class TGrid, class TVector>
	enable_if_Device<TVector, void>
	fft2_shift(const TGrid &grid, TVector &M_io)
	{
		GridBT gridBT;
		gridBT.Blk = dim3((grid.nyh+c_thrnxny-1)/c_thrnxny, (grid.nxh+c_thrnxny-1)/c_thrnxny);
		gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

		device_detail::fft2_shift<TGrid, TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, M_io); 	
	}

	template<class TGrid, class TVector>
	enable_if_Device<TVector, Value_type<TVector>>
	sum_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector &M_i)
	{
		TVector sum_v(c_thrnxny*c_thrnxny);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_over_Det<TGrid, TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_v);

		return multem::sum(grid, sum_v);
	}

	template<class TGrid, class TVector_r>
	enable_if_Device<TVector_r, Value_type<TGrid>>
	sum_square_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector_r &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, TVector_r::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_t);
		return sum(grid, sum_t);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	bandwidth_limit(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, Value_type<TVector_c> w, TVector_c &M_io)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::bandwidth_limit<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), w, M_io); 
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	phase_mul(const TGrid &grid, TVector_c &exp_x_i, TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::phase_mul<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, exp_x_i, exp_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	propagator_mul(const TGrid &grid, TVector_c &prop_x_i, TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::propagator_mul<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, prop_x_i, prop_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	probe(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_o)
	{
		using value_type_r = Value_type<TGrid>;

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::probe<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, x, y, fPsi_o);

		value_type_r total = sum_square(grid, fPsi_o);
		scale(fPsi_o, sqrt(value_type_r(grid.nxy())/total));
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	apply_CTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_CTF<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, gxu, gyu, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	apply_PCTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_PCTF<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_xyz(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_xyz<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x, k_y, k_z);

		fft2.inverse(k_x);
		fft2.inverse(k_y);
		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_x(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_x)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_x<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x);

		fft2.inverse(k_x);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_y(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_y)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_y<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_y);

		fft2.inverse(k_y);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_z(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_z<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_z);

		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_mn1(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_mn1)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mn1<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mn1);

		fft2.inverse(k_mn1);
	}

	template<class TGrid, class TVector_c>
	enable_if_Device<TVector_c, void>
	kernel_mp1(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_mp1)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mp1<TGrid, TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mp1);

		fft2.inverse(k_mp1);
	}

	template<class TVector_c>
	enable_if_Device<TVector_c, void>
	scomplex_to_complex(m_matrix_c &V_i, TVector_c &V_o)
	{
		using value_type_c = Value_type<TVector_c>;
		
		Vector<value_type_c, e_Host> V_h(V_i.size);

		for(auto ixy = 0; ixy < V_i.size; ixy++)
		{
			V_h[ixy] = value_type_c(V_i.real[ixy], V_i.imag[ixy]);
		}

		multem::assign(V_h, V_o);
	}

	template<class TGrid, class TVector_i, class TVector_o>
	typename std::enable_if<is_Device<TVector_i>::value && !is_m_matrix_rc<TVector_o>::value, void>::type
	to_host(TGrid &grid, TVector_i &M_i, TVector_o &M_o)
	{
		multem::device_synchronize<multem::e_Device>();

		thrust::copy(M_i.begin(), M_i.end(), M_o.begin());
	}

	template<class TGrid, class TVector, class Tm_matrix_r>
	typename std::enable_if<is_Device<TVector>::value && is_m_matrix_r<Tm_matrix_r>::value, void>::type
	to_host(TGrid &grid, TVector &M_i, Tm_matrix_r &M_o)
	{
		multem::device_synchronize<multem::e_Device>();

		thrust::copy(M_i.begin(), M_i.end(), M_o.real);
	}

	template<class TGrid, class TVector, class Tm_matrix_c>
	typename std::enable_if<is_Device<TVector>::value && is_m_matrix_c<Tm_matrix_c>::value, void>::type
	to_host(TGrid &grid, TVector &M_i, Tm_matrix_c &M_o)
	{
		multem::device_synchronize<multem::e_Device>();

		Vector<Value_type<TVector>, e_Host> M_t(M_i.begin(), M_i.end());
		for(auto i = 0; i<M_t.size(); i++)
		{
			M_o.real[i] = M_t[i].real();
			M_o.imag[i] = M_t[i].imag();
		}
	}

	inline
	bool is_gpu_available()
	{
		int device_count = 0;
		cudaError_t error_id = cudaGetDeviceCount(&device_count);

		if ((error_id != cudaSuccess)||(device_count == 0))
		{
			return false;
		}
		return true;
	}
} // namespace multem

#endif