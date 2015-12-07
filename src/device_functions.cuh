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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef DEVICE_FUNCTIONS_H
#define DEVICE_FUNCTIONS_H

#include <thread>
#include <type_traits>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft2.cuh"

#include "host_device_functions.cuh"
#include "host_functions.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_atomic_functions.h>
#include <cufft.h>

#define reduce_array_256(tid, sum, Mshare)						\
	{															\
		Mshare[tid] = sum;										\
		__syncthreads();										\
		if(tid < 128)											\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 128];		\
		}														\
		__syncthreads();										\
		if(tid < 64)											\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 64];			\
		}														\
		__syncthreads();										\
		if(tid < 32)											\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 32];			\
		}														\
		__syncthreads();										\
		if(tid < 16)											\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 16];			\
		}														\
		__syncthreads();										\
		if(tid < 8)												\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 8];			\
		}														\
		__syncthreads();										\
		if(tid < 4)												\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 4];			\
		}														\
		__syncthreads();										\
		if(tid < 2)												\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 2];			\
		}														\
		__syncthreads();										\
		if(tid < 1)												\
		{														\
			Mshare[tid] = sum = sum + Mshare[tid + 1];			\
		}														\
		__syncthreads();										\
	}

namespace multem
{
	namespace device_detail
	{

		template<class TGrid>
		GridBT get_grid_nxny(TGrid &grid, const dim3 Blk_max = dim3(0, 0, 0))
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
			double __shfl_down(double var, unsigned int srcLane, int width =32)
		{
			int2 a = *reinterpret_cast<int2*>(&var);
			a.x = __shfl_down(a.x, srcLane, width);
			a.y = __shfl_down(a.y, srcLane, width);
			return *reinterpret_cast<double*>(&a);
		}

		// modify
		template<class T>
		__global__ void atom_cost_function(Grid<T> grid, Atom_Ip<T> atom_Ip, rVector<T> M_i, rVector<T> Mp_o)
		{ 
			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;

			T sum = 0;

			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
			while (ix0 < atom_Ip.ixn)
			{
				int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
				while (iy0 < atom_Ip.iyn)
				{
					int ix = ix0 + atom_Ip.ix0;
					int iy = iy0 + atom_Ip.iy0;

					T R2 = grid.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid.ind_col(ix, iy);
						T M = M_i[ixy];
						T V = host_device_detail::eval_cubic_poly<T>(R2, atom_Ip);
						sum += (V-2*M)*V;
					}
					iy0 += blockDim.y*gridDim.y;
				}
				ix0 += blockDim.x*gridDim.x;
			}
			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[tid] = sum;
			}
		}

		template<class T>
		__global__ void subtract_atom(Grid<T> grid, Atom_Ip<T> atom_Ip, rVector<T> M_i)
		{ 
			int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
			while (ix0 < atom_Ip.ixn)
			{
				int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
				while (iy0 < atom_Ip.iyn)
				{
					int ix = ix0 + atom_Ip.ix0;
					int iy = iy0 + atom_Ip.iy0;

					T R2 = grid.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid.ind_col(ix, iy);
						T V = host_device_detail::eval_cubic_poly<T>(R2, atom_Ip);

						atomicAdd<T>(&(M_i.V[ixy]), V);
					}
					iy0 += blockDim.y*gridDim.y;
				}
				ix0 += blockDim.x*gridDim.x;
			}
		}

		// Linear projected potential: V and zV
		template<ePotential_Type potential_type, int charge, class T> 
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
			Vr_dVrir<potential_type, charge, T>(r, atom_Vp.cl, atom_Vp.cnl, a*w, V, dVir);

			V0s[ix] = V; 
			dV0s[ix] = dVir;

			if (atom_Vp.split)
			{
				a = atom_Vp.zeh;
				b = atom_Vp.zeh;
				z = a*x + b; 
				r = sqrt(z*z + R2);
				Vr_dVrir<potential_type, charge, T>(r, atom_Vp.cl, atom_Vp.cnl, a*w, V, dVir);

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
				atom_Vp.c0[iR] = V; 		// V_0
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
				host_device_detail::cubic_poly_coef(iR, atom_Vp);
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

				T R2 = grid.R2(ix, iy, atom_Vp.x, atom_Vp.y);
				if (R2 < atom_Vp.R2_max)
				{
					int ixy;
					T V = host_device_detail::eval_cubic_poly(ix, iy, grid, R2, atom_Vp, ixy);

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
				host_device_detail::fft2_shift(ix, iy, grid, M_io);
			}
		}

		// sum over the detector
		template<class TGrid, class T>
		__global__ void sum_over_Det(TGrid grid, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, rVector<T> M_i, rVector<T> Mp_o)
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
					host_device_detail::sum_over_Det(ix, iy, grid, g2_min, g2_max, M_i, sum);
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}
			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
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
					host_device_detail::sum_square_over_Det(ix, iy, grid, g2_min, g2_max, M_i, sum);
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		// sum over the detector
		template<class TGrid, class T>
		__global__ void sum_square_over_Det(TGrid grid, rVector<Value_type<TGrid>> S_i, rVector<T> M_i, rVector<Value_type<TGrid>> Mp_o)
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
					host_device_detail::sum_square_over_Det(ix, iy, grid, S_i, M_i, sum);
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		// Anti-Aliasing, scale with cut-off (2/3)g_max
		template<class TGrid, class T>
		__global__ void bandwidth_limit(TGrid grid, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, T w, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::bandwidth_limit(ix, iy, grid, g2_min, g2_max, w, M_io);
			}
		}

		// Phase multiplication
		template<class TGrid, class T>
		__global__ void phase_multiplication(TGrid grid, rVector<T> exp_x_i, 
		rVector<T> exp_y_i, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::phase_multiplication(ix, iy, grid, exp_x_i, exp_y_i, psi_i, psi_o);
			}
		}

		// Propagate, scale with cut-off (2/3)g_max
		template<class TGrid, class T>
		__global__ void propagator_multiplication(TGrid grid, rVector<T> prop_x_i, 
		rVector<T> prop_y_i, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::propagator_multiplication(ix, iy, grid, prop_x_i, prop_y_i, psi_i, psi_o);
			}
		}

		// shift function in Fourier space
		template<class TGrid, class T>
		__global__ void phase_factor_2D(TGrid grid, Value_type<TGrid> x, Value_type<TGrid> y, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::phase_factor_2D(ix, iy, grid, x, y, psi_i, psi_o);
			}
		}

		// Convergent incident wave in Fourier space
		template<class TGrid, class T>
		__global__ void probe(TGrid grid, Lens<Value_type<TGrid>> lens, 
		Value_type<TGrid> x, Value_type<TGrid> y, rVector<T> fPsi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{	
				host_device_detail::probe(ix, iy, grid, lens, x, y, fPsi_o);
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
				host_device_detail::apply_CTF(ix, iy, grid, lens, gxu, gyu, fPsi_i, fPsi_o);
			}
		}

		// Partially coherent transfer function, linear image model and weak phase_components object
		template<class TGrid, class T>
		__global__ void apply_PCTF(TGrid grid, Lens<Value_type<TGrid>> lens, rVector<T> fPsi_i, rVector<T> fPsi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::apply_PCTF(ix, iy, grid, lens, fPsi_i, fPsi_o);
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
					host_device_detail::Lorentz_factor(ix, iy, grid, gc2, ge2, sum);
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		template<class TGrid>
		Value_type<TGrid> Lorentz_factor(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels)
		{
			device_vector<Value_type<TGrid>> sum_v(c_thrnxny*c_thrnxny);

			auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

			norm_factor_lorentz<Value_type<TGrid>><<<gridBT.Blk, gridBT.Thr>>>(grid, eels.gc2, eels.ge2, sum_v);

			return sqrt(eels.occ)/thrust::reduce(sum_v.begin(), sum_v.end());
		}

		template<class TGrid, class T>
		__global__ void kernel_xyz(TGrid grid, EELS<Value_type<TGrid>> eels, 
		rVector<T> k_x, rVector<T> k_y, rVector<T> k_z)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_xyz(ix, iy, grid, eels, k_x, k_y, k_z);
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_x(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_x)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_x(ix, iy, grid, eels, k_x);
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_y(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_y)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_y(ix, iy, grid, eels, k_y);
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_z(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_z)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_z(ix, iy, grid, eels, k_z);
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_mn1(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_mn1)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_mn1(ix, iy, grid, eels, k_mn1);
			}
		}

		template<class TGrid, class T>
		__global__ void kernel_mp1(TGrid grid, EELS<Value_type<TGrid>> eels, rVector<T> k_mp1)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid.nx)&&(iy < grid.ny))
			{
				host_device_detail::kernel_mp1(ix, iy, grid, eels, k_mp1);
			}
		}

	} // namespace device_detail

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	assign(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h =nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template<class TVector>
	enable_if_device_vector<TVector, void>
	fill(Stream<e_device> &stream, TVector &M_io, Value_type<TVector> value_i)
	{
		thrust::fill(M_io.begin(), M_io.end(), value_i);
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::scale<value_type>(w_i));
	}

	template<class TVector>
	enable_if_device_vector<TVector, void>
	scale(Stream<e_device> &stream, Value_type<TVector> w_i, TVector &M_io)
	{
		scale(stream, w_i, M_io, M_io);
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	square(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square<value_type>());
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	square_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square_scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add<value_type>());
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		add(stream, M_i, M_io, M_io);
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_device> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, 
	Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_scale_i<value_type>(w1_i, w2_i));
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_i<value_type>());
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_square<value_type>());
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_device> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_scale_i<value_type>(w1_i, w2_i));
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_square_scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::multiply<value_type>());
	}

	template<class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		multiply(stream, M_i, M_io, M_io);
	}

	template<class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	sum(Stream<e_device> &stream, TVector &M_i)
	{
		return thrust::reduce(M_i.begin(), M_i.end());
	}

	template<class TVector>
	enable_if_device_vector<TVector, Value_type_r<TVector>>
	sum_square(Stream<e_device> &stream, TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;
		return thrust::transform_reduce(M_i.begin(), M_i.end(), 
		functor::square<value_type_r>(), value_type_r(0), functor::add<value_type_r>());
	}

	template<class TVector>
	enable_if_device_vector<TVector, Value_type_r<TVector>>
	mean(Stream<e_device> &stream, TVector &M_i)
	{
		return sum(stream, M_i)/M_i.size();
	}

	template<class TVector>
	enable_if_device_vector<TVector, void>
	mean_std(Stream<e_device> &stream, TVector &M_i, Value_type_r<TVector> &x_mean, Value_type_r<TVector> &x_std)
	{
		using value_type_r = Value_type_r<TVector>;

		x_mean = mean(stream, M_i);
		x_std = thrust::transform_reduce(M_i.begin(), M_i.begin(), functor::square_dif<value_type_r>(x_mean), value_type_r(0), functor::add<value_type_r>());

		x_std = sqrt(x_std/M_i.size());
	}

	template<class TVector>
	enable_if_device_vector<TVector, Value_type_r<TVector>>
	std(Stream<e_device> &stream, TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;

		value_type_r x_mean , x_std;
		mean_std(stream, M_i, x_mean, x_std);
		return x_std;
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	phase_factor_1D(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> x, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		thrust::counting_iterator<int> first(0);
		thrust::counting_iterator<int> last = first + grid.nx;

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, fPsi_i.begin(), fPsi_o.begin())), 
						 thrust::make_zip_iterator(thrust::make_tuple(last, fPsi_i.end(), fPsi_o.end())), 
						 functor::phase_factor<TGrid>(grid, x));
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	phase_factor_2D(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		using value_type_r = Value_type<TGrid>;

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::phase_factor_2D<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, x, y, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, Value_type<TVector_r>>
	atom_cost_function(TGrid &grid, const Atom_Ip<Value_type<TGrid>> &atom_Ip, TVector_r &M_i)
	{
		TVector_r sum_v(1);

		GridBT gridBT;
		gridBT.Blk = dim3();
		gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

		device_detail::atom_cost_function<typename TVector_r::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, atom_Ip, M_i, sum_v);
		return sum_v[0];
	}

	template<class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, void>
	subtract_atom(Stream<e_device> &stream, TGrid &grid, Vector<Atom_Ip<Value_type<TGrid>>, e_host> &atom_Ip, TVector_r &M_i)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			GridBT gridBT;
			gridBT.Blk = dim3();
			gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

			device_detail::subtract_atom<typename TVector_r::value_type><<<gridBT.Blk, gridBT.Thr, 0, stream[istream]>>>(grid, atom_Ip[istream], M_i);
		}
	}

	template<class TQ1>
	enable_if_device<TQ1, void>
	get_cubic_poly_coef_Vz(Stream<e_device> &stream, ePotential_Type potential_type, 
	TQ1 &qz, Vector<Atom_Vp<Value_type<TQ1>>, e_host> &atom_Vp)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			if(atom_Vp[istream].charge==0)
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						device_detail::linear_Vz<ePT_Doyle_0_4, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Peng_0_4:
						device_detail::linear_Vz<ePT_Peng_0_4, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Peng_0_12:
						device_detail::linear_Vz<ePT_Peng_0_12, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Kirkland_0_12:
						device_detail::linear_Vz<ePT_Kirkland_0_12, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Weickenmeier_0_12:
						device_detail::linear_Vz<ePT_Weickenmeier_0_12, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Lobato_0_12:
						device_detail::linear_Vz<ePT_Lobato_0_12, 0, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
				}
			}
			else
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						device_detail::linear_Vz<ePT_Doyle_0_4, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Peng_0_4:
						device_detail::linear_Vz<ePT_Peng_0_4, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Peng_0_12:
						device_detail::linear_Vz<ePT_Peng_0_12, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Kirkland_0_12:
						device_detail::linear_Vz<ePT_Kirkland_0_12, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Weickenmeier_0_12:
						device_detail::linear_Vz<ePT_Weickenmeier_0_12, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
					case ePT_Lobato_0_12:
						device_detail::linear_Vz<ePT_Lobato_0_12, 1, typename TQ1::value_type><<<dim3(c_nR), dim3(c_nqz), 0, stream[istream]>>>(qz, atom_Vp[istream]);
						break;
				}		
			}
			device_detail::cubic_poly_coef<typename TQ1::value_type><<<dim3(1), dim3(c_nR), 0, stream[istream]>>>(atom_Vp[istream]);
		}
	}

	template<class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, void>
	eval_cubic_poly(Stream<e_device> &stream, TGrid &grid, 
	Vector<Atom_Vp<Value_type<TGrid>>, e_host> &atom_Vp, TVector_r &V_0)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			GridBT gridBT = atom_Vp[istream].get_eval_cubic_poly_gridBT();

			device_detail::eval_cubic_poly<typename TVector_r::value_type><<<gridBT.Blk, gridBT.Thr, 0, stream[istream]>>>(grid, atom_Vp[istream], V_0);
		}
	}

	template<class TGrid, class TVector>
	enable_if_device_vector<TVector, void>
	fft2_shift(Stream<e_device> &stream, TGrid &grid, TVector &M_io)
	{
		GridBT gridBT;
		gridBT.Blk = dim3((grid.nyh+c_thrnxny-1)/c_thrnxny, (grid.nxh+c_thrnxny-1)/c_thrnxny);
		gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

		device_detail::fft2_shift<TGrid, typename TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, M_io); 	
	}

	template<class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	sum_over_Det(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		TVector sum_v(c_thrnxny*c_thrnxny);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_over_Det<TGrid, typename TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_v);
		return multem::sum(stream, sum_v);
	}

	template<class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny, 0);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, typename TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_t);
		return sum(stream, sum_t);
	}

	template<class TGrid, class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_device> &stream, TGrid &grid, TVector_1 &S_i, TVector_2 &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny, 0);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, typename TVector_2::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, S_i, M_i, sum_t);
		return sum(stream, sum_t);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	bandwidth_limit(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, Value_type<TVector_c> w, TVector_c &M_io)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::bandwidth_limit<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), w, M_io); 
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	phase_components(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &V_x_o, TVector_c &V_y_o)
	{
		thrust::counting_iterator<int> first(0);
		thrust::counting_iterator<int> last = first + grid.nx_ny_max();

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x_o.begin(), V_y_o.begin())), 
						 thrust::make_zip_iterator(thrust::make_tuple(last, V_x_o.end(), V_y_o.end())), 
						 functor::phase_components<TGrid>(grid, c_2Pi, gxu, gyu));
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	phase_multiplication(Stream<e_device> &stream, TGrid &grid, TVector_c &exp_x_i, TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::phase_multiplication<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, exp_x_i, exp_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	propagator_components(Stream<e_device> &stream, TGrid &grid, Value_type<TGrid> gxu, Value_type<TGrid> gyu, Value_type<TGrid> w, TVector_c &V_x_o, TVector_c &V_y_o)
	{
		thrust::counting_iterator<int> first(0);
		thrust::counting_iterator<int> last = first + grid.nx_ny_max();

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x_o.begin(), V_y_o.begin())), 
						 thrust::make_zip_iterator(thrust::make_tuple(last, V_x_o.end(), V_y_o.end())), 
						 functor::propagator_components<TGrid>(grid, w, gxu, gyu));
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	propagator_multiplication(Stream<e_device> &stream, TGrid &grid, TVector_c &prop_x_i, TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::propagator_multiplication<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, prop_x_i, prop_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	transmission_function(Stream<e_device> &stream, TGrid &grid, eElec_Spec_Int_Model elec_spec_int_model, Value_type<TGrid> w, TVector_1 &V0_i, TVector_2 &Trans_o)
	{	
		using value_type_r = Value_type<TGrid>;

		thrust::transform(V0_i.begin(), V0_i.end(), Trans_o.begin(), 
			functor::transmission_function<value_type_r>(w, elec_spec_int_model));
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	probe(Stream<e_device> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_o)
	{
		using value_type_r = Value_type<TGrid>;

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::probe<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, x, y, fPsi_o);

		auto total = sum_square(stream, fPsi_o);
		scale(stream, sqrt(1.0/total), fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_CTF(Stream<e_device> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_CTF<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, gxu, gyu, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_PCTF(Stream<e_device> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_PCTF<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_xyz(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_xyz<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x, k_y, k_z);

		fft2.inverse(k_x);
		fft2.inverse(k_y);
		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_x(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_x)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_x<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x);

	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_y(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_y)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_y<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_y);

		fft2.inverse(k_y);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_z(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_z<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_z);

		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mn1(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_mn1)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mn1<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mn1);

		fft2.inverse(k_mn1);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mp1(Stream<e_device> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_device> &fft2, TVector_c &k_mp1)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mp1<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mp1);

		fft2.inverse(k_mp1);
	}

	/***************************************************************************/
	/***************************************************************************/
	template<class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector_and_host_vector<TVector_i, TVector_o>::value 
	&& is_complex<Value_type<TVector_o>>::value && !std::is_same<Value_type<TVector_i>, Value_type<TVector_o>>::value, void>::type
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o, 
	Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		// copy data from host to host
		copy_to_host(stream, *M_i_h, M_o);
	}

	template<class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector_and_host_vector<TVector_i, TVector_o>::value 
	&& (!is_complex<Value_type<TVector_o>>::value || std::is_same<Value_type<TVector_i>, Value_type<TVector_o>>::value), void>::type
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o, 
	Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template<class TVector_i, class TVector_o>
	enable_if_device_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		// add and scale
		multem::add_scale_to_host(stream, w_i, *M_i_h, M_o);
	}

	template<class TVector_i, class TVector_o>
	enable_if_device_vector_and_host_vector<TVector_i, TVector_o, void>
	add_square_scale_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		multem::add_square_scale_to_host(stream, w_i, *M_i_h, M_o);
	}

	template<class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_device_vector_and_host_vector<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i, 
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_host> *psi_i_h =nullptr)
	{
		Vector<Value_type<TVector_c_i>, e_host> M_h;
		psi_i_h = (psi_i_h == nullptr)?&M_h:psi_i_h;

		// data transfer from GPU to CPU
		psi_i_h->assign(psi_i.begin(), psi_i.end());

		multem::add_scale_m2psi_psi_to_host(stream, w_i, *psi_i_h, m2psi_o, psi_o);
	}

	/***************************************************************************/
	/***************************************************************************/

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

	inline
	int number_of_gpu_available()
	{
		int device_count = 0;
		cudaError_t error_id = cudaGetDeviceCount(&device_count);

		if (error_id != cudaSuccess)
		{
			device_count = 0;
		}
		return (device_count>0)?device_count:0;
	}
} // namespace multem

#endif