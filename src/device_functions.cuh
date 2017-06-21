/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "fft.cuh"
#include "blas.cuh"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

#include "host_device_functions.cuh"
#include "host_functions.hpp"

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

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
	__device__ __forceinline__
	double atomicAdd(double *address, double val)
	{
		unsigned long long int* address_as_ull = (unsigned long long int*)address;
		unsigned long long int old = *address_as_ull, assumed;

		do
		{
			assumed = old;
			old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
		} while (assumed != old);

		return __longlong_as_double(old);
	}
#endif

namespace mt
{
	namespace device_detail
	{
		// modify
		template <class T>
		__global__ void atom_cost_function(Grid_2d<T> grid_2d, Atom_Sa<T> atom_Ip, rVector<T> M_i, rVector<T> Mp_o)
		{ 
			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;

			T sum = 0;

			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;
			while (ix_0 < atom_Ip.ixn)
			{
				int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
				while (iy_0 < atom_Ip.iyn)
				{
					int ix = ix_0 + atom_Ip.ix_0;
					int iy = iy_0 + atom_Ip.iy_0;

					T R2 = grid_2d.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid_2d.ind_col(ix, iy);
						T M = M_i[ixy];
						const T V = host_device_detail::eval_cubic_poly(R2, atom_Ip);
						sum += (V-2*M)*V;
					}
					iy_0 += blockDim.x*gridDim.x;
				}
				ix_0 += blockDim.y*gridDim.y;
			}
			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[tid] = sum;
			}
		}

		template <class T>
		__global__ void subtract_atom(Grid_2d<T> grid_2d, Atom_Sa<T> atom_Ip, rVector<T> M_i)
		{ 
			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;
			while (ix_0 < atom_Ip.ixn)
			{
				int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
				while (iy_0 < atom_Ip.iyn)
				{
					int ix = ix_0 + atom_Ip.ix_0;
					int iy = iy_0 + atom_Ip.iy_0;

					T R2 = grid_2d.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid_2d.ind_col(ix, iy);
						const T V = host_device_detail::eval_cubic_poly(R2, atom_Ip);

						atomicAdd(&(M_i.V[ixy]), V);
						// atomicAdd<T>(&(M_i.V[ixy]), V);
					}
					iy_0 += blockDim.x*gridDim.x;
				}
				ix_0 += blockDim.y*gridDim.y;
			}
		}

		// Linear projected potential: V and zV
		template <ePotential_Type potential_type, int charge, class TAtom> 
		__global__ void linear_Vz(rQ1<Value_type<TAtom>> qz, TAtom atom)
		{	
			using T = Value_type<TAtom>;

			__shared__ T V0s[c_nqz];
			__shared__ T dV0s[c_nqz];

			int ix = threadIdx.x, iR = blockIdx.x;

			T x = qz.x[ix];
			T w = qz.w[ix];
			T R2 = atom.R2[iR];
			T V, dVir;

			T a = (atom.split)?(-atom.z0h):(atom.zeh-atom.z0h);
			T b = (atom.split)?(atom.z0h):(atom.zeh+atom.z0h);
			T z = a*x + b;
			T r = sqrt(z*z + R2);
			Vr_dVrir<potential_type, charge, T>(r, atom.cl, atom.cnl, a*w, V, dVir);

			V0s[ix] = V; 
			dV0s[ix] = dVir;

			if (atom.split)
			{
				a = atom.zeh;
				b = atom.zeh;
				z = a*x + b; 
				r = sqrt(z*z + R2);
				Vr_dVrir<potential_type, charge, T>(r, atom.cl, atom.cnl, a*w, V, dVir);

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
				dVir = 0.5*dVir;

				auto R2_tap = atom.R2_tap;
				auto tap_cf = atom.tap_cf;
				host_device_detail::apply_tapering(R2_tap, tap_cf, R2, V, dVir);
				atom.c0[iR] = V;			// V_0
				atom.c1[iR] = dVir; 		// dR2V0
			}
		}

		// Get Local interpolation coefficients
		template <class TAtom> 
		__global__ void cubic_poly_coef(TAtom atom)
		{
			int iR = threadIdx.x;

			if (iR < c_nR-1)
			{
				host_device_detail::cubic_poly_coef(iR, atom);
			}
		}

		// Cubic polynomial evaluation
		template <class T> 
		__global__ void eval_cubic_poly(Grid_2d<T> grid_2d, rVector<Atom_Vp<T>> atoms, rVector<T> M_o)
		{
			auto atom = atoms.V[blockIdx.x];

			int ix_t = threadIdx.y;
			while (ix_t < atom.nx)
			{
				int ix = atom.ix_0+ix_t;

				int iy_t = threadIdx.x;
				while (iy_t < atom.ny)
				{
					int iy = atom.iy_0+iy_t;
					const auto R2 = grid_2d.R2(ix, iy, atom.x, atom.y);

					if (R2 < atom.R2_max)
					{
						const T V = atom.occ*host_device_detail::eval_cubic_poly(R2, atom);
						const int ixy = grid_2d.ind_col_pbc_shift(ix, iy);

						atomicAdd(&(M_o.V[ixy]), V);
					}

					iy_t += blockDim.x;
				}
				ix_t += blockDim.y;
			}
		}

		// Gaussian evaluation
		template <class T> 
		__global__ void gauss_eval(Grid_2d<T> grid_2d, rVector<Gauss_Sp<T>> atoms, rVector<T> M_o)
		{
			const auto gauss = atoms.V[blockIdx.x];

			int ix_t = threadIdx.y;
			while (ix_t < gauss.nx)
			{
				int ix = gauss.ix_0+ix_t;

				int iy_t = threadIdx.x;
				while (iy_t < gauss.ny)
				{
					int iy = gauss.iy_0+iy_t;
					auto R2 = grid_2d.R2(ix, iy, gauss.x, gauss.y);

					if (R2 < gauss.R2_max)
					{
						int ixy = grid_2d.ind_col_pbc(ix, iy);
						atomicAdd(&(M_o.V[ixy]), gauss(R2));
					}

					iy_t += blockDim.x;
				}
				ix_t += blockDim.y;
			}
		}

		// Shift matrix respect to nxh
		template <class TGrid, class T>
		__global__ void fft1_shift(TGrid grid_1d, rVector<T> M_io)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nxh)
			{
				host_device_detail::fft1_shift(ix, grid_1d, M_io);
			}
		}

		// Shift matrix respect to nyh
		template <class TGrid, class T>
		__global__ void fft2_sft_bc(TGrid grid_2d, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.nyh))
			{
				host_device_detail::fft2_sft_bc(ix, iy, grid_2d, M_io);
			}
		}

		// Shift matrix respect to (nxh, nyh)
		template <class TGrid, class T>
		__global__ void fft2_shift(TGrid grid_2d, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nxh)&&(iy < grid_2d.nyh))
			{
				host_device_detail::fft2_shift(ix, iy, grid_2d, M_io);
			}
		}

		// sum over the detector
		template <class TGrid, class T>
		__global__ void sum_over_Det(TGrid grid_2d, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, rVector<T> M_i, rVector<T> Mp_o)
		{ 
			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;

			T sum = 0;

			int ix = ix_0;
			while (ix < grid_2d.nx)
			{
				int iy = iy_0;
				while (iy < grid_2d.ny)
				{
					host_device_detail::sum_over_Det(ix, iy, grid_2d, g2_min, g2_max, M_i, sum);
					iy += blockDim.x*gridDim.x;
				}
				ix += blockDim.y*gridDim.y;
			}
			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		// sum over the detector
		template <class TGrid, class T>
		__global__ void sum_square_over_Det(TGrid grid_2d, Value_type<TGrid> g2_min, Value_type<TGrid> g2_max, rVector<T> M_i, rVector<Value_type<TGrid>> Mp_o)
		{ 
			using T_r = Value_type<TGrid>;

			__shared__ T_r Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;

			T_r sum = 0;

			int ix = ix_0;
			while (ix < grid_2d.nx)
			{
				int iy = iy_0;
				while (iy < grid_2d.ny)
				{
					host_device_detail::sum_square_over_Det(ix, iy, grid_2d, g2_min, g2_max, M_i, sum);
					iy += blockDim.x*gridDim.x;
				}
				ix += blockDim.y*gridDim.y;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		// sum over the detector
		template <class TGrid, class T>
		__global__ void sum_square_over_Det(TGrid grid_2d, rVector<Value_type<TGrid>> S_i, rVector<T> M_i, rVector<Value_type<TGrid>> Mp_o)
		{ 
			using T_r = Value_type<TGrid>;

			__shared__ T_r Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;

			T_r sum = 0;

			int ix = ix_0;
			while (ix < grid_2d.nx)
			{
				int iy = iy_0;
				while (iy < grid_2d.ny)
				{
					host_device_detail::sum_square_over_Det(ix, iy, grid_2d, S_i, M_i, sum);
					iy += blockDim.x*gridDim.x;
				}
				ix += blockDim.y*gridDim.y;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		// Anti-Aliasing, scale with cut-off (2/3)g_max
		template <class TGrid, class T>
		__global__ void bandwidth_limit(TGrid grid_2d, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::bandwidth_limit(ix, iy, grid_2d, M_io);
			}
		}

		template <class TGrid, class T>
		__global__ void hard_aperture(TGrid grid_2d, Value_type<TGrid> g2_max, Value_type<TGrid> w, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::hard_aperture(ix, iy, grid_2d, g2_max, w, M_io);
			}
		}

		// Propagate
		template <class TGrid, class T>
		__global__ void propagate(TGrid grid_2d, Value_type<TGrid> w, 
		Value_type<TGrid> gxu, Value_type<TGrid> gyu, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::propagate(ix, iy, grid_2d, w, gxu, gyu, psi_i, psi_o);
			}
		}

		/***********************************************************************/
		// phase factor 1d
		template <class TGrid, class T>
		__global__ void exp_r_factor_1d(TGrid grid_1d, Value_type<TGrid> gx, rVector<T> psi_i, rVector<T> psi_o)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::exp_r_factor_1d(ix, grid_1d, gx, psi_i, psi_o);
			}
		}

		// phase factor 2d by col
		template <class TGrid, class T>
		__global__ void exp_r_factor_2d_bc(TGrid grid_2d, Value_type<TGrid> alpha, 
		rVector<Value_type<TGrid>> gy, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::exp_r_factor_2d_bc(ix, iy, grid_2d, alpha, gy, psi_i, psi_o);
			}
		}

		// phase factor 2d
		template <class TGrid, class T>
		__global__ void exp_r_factor_2d(TGrid grid_2d, Value_type<TGrid> gx, 
		Value_type<TGrid> gy, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::exp_r_factor_2d(ix, iy, grid_2d, gx, gy, psi_i, psi_o);
			}
		}

		// phase factor 2d
		template <class TGrid, class T>
		__global__ void mul_exp_r_factor_2d(TGrid grid_2d, rVector<Value_type<TGrid>> gx, 
		rVector<Value_type<TGrid>> gy, rVector<T> psi_i, rVector<T> psi_o)
		{
			using Tr = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::mul_exp_r_factor_2d(ix, iy, grid_2d, gx, gy, psi_i, psi_o);
			}
		}

		/***********************************************************************/
		// phase factor 1d
		template <class TGrid, class T>
		__global__ void exp_g_factor_1d(TGrid grid_1d, Value_type<TGrid> Rx, rVector<T> psi_i, rVector<T> psi_o)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::exp_g_factor_1d(ix, grid_1d, Rx, psi_i, psi_o);
			}
		}

		// phase factor 2d by col
		template <class TGrid, class T>
		__global__ void exp_g_factor_2d_bc(TGrid grid_2d, Value_type<TGrid> alpha, 
		rVector<Value_type<TGrid>> Ry, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::exp_g_factor_2d_bc(ix, iy, grid_2d, alpha, Ry, psi_i, psi_o);
			}
		}

		// phase factor 2d
		template <class TGrid, class T>
		__global__ void exp_g_factor_2d(TGrid grid_2d, Value_type<TGrid> Rx, 
		Value_type<TGrid> Ry, rVector<T> psi_i, rVector<T> psi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::exp_g_factor_2d(ix, iy, grid_2d, Rx, Ry, psi_i, psi_o);
			}
		}

		// phase factor 2d
		template <class TGrid, class T>
		__global__ void mul_exp_g_factor_2d(TGrid grid_2d, rVector<Value_type<TGrid>> Rx, 
		rVector<Value_type<TGrid>> Ry, rVector<T> psi_i, rVector<T> psi_o)
		{
			using Tr = Value_type<TGrid>;

			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::mul_exp_g_factor_2d(ix, iy, grid_2d, Rx, Ry, psi_i, psi_o);
			}
		}

		/***********************************************************************/
		// Convergent incident wave in Fourier space
		template <class TGrid, class T>
		__global__ void probe(TGrid grid_2d, Lens<Value_type<TGrid>> lens, 
		Value_type<TGrid> x, Value_type<TGrid> y, Value_type<TGrid> gxu, Value_type<TGrid> gyu, rVector<T> fPsi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{	
				host_device_detail::probe(ix, iy, grid_2d, lens, x, y, gxu, gyu, fPsi_o);
			}
		}

		// Apply Coherent transfer function
		template <class TGrid, class T>
		__global__ void apply_CTF(TGrid grid_2d, Lens<Value_type<TGrid>> lens, 
		Value_type<TGrid> gxu, Value_type<TGrid> gyu, rVector<T> fPsi_i, rVector<T> fPsi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{	
				host_device_detail::apply_CTF(ix, iy, grid_2d, lens, gxu, gyu, fPsi_i, fPsi_o);
			}
		}

		// Partially coherent transfer function, linear image model and weak phase_components object
		template <class TGrid, class T>
		__global__ void apply_PCTF(TGrid grid_2d, Lens<Value_type<TGrid>> lens, rVector<T> fPsi_i, rVector<T> fPsi_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::apply_PCTF(ix, iy, grid_2d, lens, fPsi_i, fPsi_o);
			}
		}

		// sum over the detector
		template <class T>
		__global__ void norm_factor_lorentz(Grid_2d<T> grid_2d, T gc2, T ge2, rVector<T> Mp_o)
		{ 
			__shared__ T Mshare[c_thrnxny*c_thrnxny];

			int tid = threadIdx.x + threadIdx.y*blockDim.x;
			int bid = blockIdx.x + blockIdx.y*gridDim.x;

			int iy_0 = threadIdx.x + blockIdx.x*blockDim.x;
			int ix_0 = threadIdx.y + blockIdx.y*blockDim.y;

			T sum = 0;

			int ix = ix_0;
			while (ix < grid_2d.nx)
			{
				int iy = iy_0;
				while (iy < grid_2d.ny)
				{
					host_device_detail::Lorentz_factor(ix, iy, grid_2d, gc2, ge2, sum);
					iy += blockDim.x*gridDim.x;
				}
				ix += blockDim.y*gridDim.y;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0)
			{
				Mp_o[bid] = sum;
			}
		}

		template <class TGrid>
		Value_type<TGrid> Lorentz_factor(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels)
		{
			device_vector<Value_type<TGrid>> sum_v(c_thrnxny*c_thrnxny);

			auto grid_bt = grid_2d.cuda_grid(dim3(c_thrnxny, c_thrnxny));

			norm_factor_lorentz<Value_type<TGrid>><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels.gc2, eels.ge2, sum_v);

			return sqrt(eels.occ)/thrust::reduce(sum_v.begin(), sum_v.end());
		}

		template <class TGrid, class T>
		__global__ void kernel_xyz(TGrid grid_2d, EELS<Value_type<TGrid>> eels, 
		rVector<T> k_x, rVector<T> k_y, rVector<T> k_z)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_xyz(ix, iy, grid_2d, eels, k_x, k_y, k_z);
			}
		}

		template <class TGrid, class T>
		__global__ void kernel_x(TGrid grid_2d, EELS<Value_type<TGrid>> eels, rVector<T> k_x)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_x(ix, iy, grid_2d, eels, k_x);
			}
		}

		template <class TGrid, class T>
		__global__ void kernel_y(TGrid grid_2d, EELS<Value_type<TGrid>> eels, rVector<T> k_y)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_y(ix, iy, grid_2d, eels, k_y);
			}
		}

		template <class TGrid, class T>
		__global__ void kernel_z(TGrid grid_2d, EELS<Value_type<TGrid>> eels, rVector<T> k_z)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_z(ix, iy, grid_2d, eels, k_z);
			}
		}

		template <class TGrid, class T>
		__global__ void kernel_mn1(TGrid grid_2d, EELS<Value_type<TGrid>> eels, rVector<T> k_mn1)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_mn1(ix, iy, grid_2d, eels, k_mn1);
			}
		}

		template <class TGrid, class T>
		__global__ void kernel_mp1(TGrid grid_2d, EELS<Value_type<TGrid>> eels, rVector<T> k_mp1)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::kernel_mp1(ix, iy, grid_2d, eels, k_mp1);
			}
		}

		// trs
		template <class T>
		__global__ void trs(int ncols, int nrows, rVector<T> M_i, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < ncols)&&(iy < nrows))
			{
				host_device_detail::trs(ix, iy, ncols, nrows, M_i, M_o);
			}
		}

		// vector column x matrix
		template <class TGrid, class T>
		__global__ void vector_col_x_matrix(TGrid grid_2d, rVector<Value_type_r<T>> fg, 
		rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::vector_col_x_matrix(ix, iy, grid_2d, fg, M_io);
			}
		}

		// Gaussian convolution
		template <class TGrid, class T>
		__global__ void gauss_cv_1d(TGrid grid_1d, Value_type<TGrid> alpha, 
		rVector<T> M_io)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::gauss_cv_1d(ix, grid_1d, alpha, M_io);
			}
		}

		// Gaussian convolution
		template <class TGrid, class T>
		__global__ void gauss_cv_2d(TGrid grid_2d, 
		Value_type<TGrid> alpha, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::gauss_cv_2d(ix, iy, grid_2d, alpha, M_io);
			}
		}

		// Gaussian convolution
		template <class TGrid, class T>
		__global__ void gauss_dcv_1d(TGrid grid_1d, Value_type<TGrid> alpha, 
		Value_type<TGrid> PSNR, rVector<T> M_io)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::gauss_dcv_1d(ix, grid_1d, alpha, PSNR, M_io);
			}
		}

		// Gaussian deconvolution
		template <class TGrid, class T>
		__global__ void gauss_dcv_2d(TGrid grid_2d, 
		Value_type<TGrid> alpha, Value_type<TGrid> PSNR, rVector<T> M_io)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::gauss_dcv_2d(ix, iy, grid_2d, alpha, PSNR, M_io);
			}
		}

		// pcf_1d preprocessing
		template <class TGrid, class T>
		__global__ void pcf_1d_pp(int ix_s, TGrid grid_1d, Butterworth_1d<Value_type<TGrid>> bw_1d, 
		rVector<Value_type_r<T>> M_i, rVector<T> M_o)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::pcf_1d_pp(ix, ix_s, grid_1d, bw_1d, M_i, M_o);
			}
		}

		// gaussian 1d envelope
		template <class TGrid, class T>
		__global__ void pcf_1d_gaussian(TGrid grid_1d, Gauss_1d<Value_type<TGrid>> gs_1d, 
		rVector<T> M_1, rVector<T> M_2, rVector<T> pcf)
		{
			int ix = threadIdx.x + blockIdx.x*blockDim.x;

			if(ix < grid_1d.nx)
			{
				host_device_detail::pcf_1d_gaussian(ix, grid_1d, gs_1d, M_1, M_2, pcf);
			}
		}

		// pcf_2d by col preprocessing
		template <class TGrid, class T>
		__global__ void pcf_2d_bc_pp(int iy_s, TGrid grid_2d_i, TGrid grid_2d_o, 
		rVector<Value_type_r<T>> M_i, rVector<Value_type_r<T>> fh, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_i.nx)&&(iy < grid_2d_i.ny))
			{
				host_device_detail::pcf_2d_bc_pp(ix, iy, iy_s, grid_2d_i, grid_2d_o, M_i, fh, M_o);
			}
		}

		// gaussian 2d by col envelope
		template <class TGrid, class T>
		__global__ void pcf_2d_bc_gaussian(TGrid grid_2d, rVector<T> M_1, rVector<T> M_2, 
		rVector<Value_type_r<T>> fg, rVector<T> pcf)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::pcf_2d_bc_gaussian(ix, iy, grid_2d, M_1, M_2, fg, pcf);
			}
		}

		// gaussian 2d by col envelope
		template <class TGrid, class T>
		__global__ void pcf_2d_bc_assign_real(int iy_s, TGrid grid_2d_i, TGrid grid_2d_o, 
		rVector<T> M_i, rVector<Value_type_r<T>> M_o, bool b_pos)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_i.nx)&&(iy < grid_2d_i.ny))
			{
				int ixy_i = grid_2d_o.ind_col(ix, iy+iy_s);
				int ixy_o = grid_2d_i.ind_col(ix, iy);
				auto v = M_i[ixy_i].real();
				M_o[ixy_o] = (!b_pos || (v>0))?v:0;
			}
		}

		// pcf_2d preprocessing
		template <class TGrid, class T>
		__global__ void pcf_2d_pp(TGrid grid_2d, Butterworth_2d<Value_type<TGrid>> bw_2d, 
		rVector<Value_type_r<T>> M_i, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::pcf_2d_pp(ix, iy, grid_2d, bw_2d, M_i, M_o);
			}
		}

		// gaussian 2d envelope
		template <class TGrid, class T>
		__global__ void pcf_2d_gaussian(TGrid grid_2d, Gauss_2d<Value_type<TGrid>> gs_2d, 
		rVector<T> M_1, rVector<T> M_2, rVector<T> pcf)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d.nx)&&(iy < grid_2d.ny))
			{
				host_device_detail::pcf_2d_gaussian(ix, iy, grid_2d, gs_2d, M_1, M_2, pcf);
			}
		}

		// scaling xy
		template <class TGrid, class T>
		__global__ void sc_2d(TGrid grid_2d_i, rVector<T> M_i, T fxy, TGrid grid_2d_o, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_o.nx)&&(iy < grid_2d_o.ny))
			{
				host_device_detail::sc_2d(ix, iy, grid_2d_i, M_i, fxy, grid_2d_o, M_o);
			}
		}

		// rotate xy
		template <class TGrid, class T>
		__global__ void rot_2d(TGrid grid_2d_i, rVector<T> M_i, T theta, r2d<T> p0, 
		T bg, TGrid grid_2d_o, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_o.nx)&&(iy < grid_2d_o.ny))
			{
				host_device_detail::rot_2d(ix, iy, grid_2d_i, M_i, theta, p0, bg, grid_2d_o, M_o);
			}
		}

		// rotate, scale and shift xy
		template <class TGrid, class T>
		__global__ void rot_sca_sft_2d(TGrid grid_2d_i, rVector<T> M_i, T theta, r2d<T> p0, 
		T fxy, r2d<T> ps, T bg, TGrid grid_2d_o, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_o.nx)&&(iy < grid_2d_o.ny))
			{
				host_device_detail::rot_sca_sft_2d(ix, iy, grid_2d_i, M_i, theta, p0, fxy, ps, bg, grid_2d_o, M_o);
			}
		}

		// x-shear and y-scaling
		template <class TGrid, class T>
		__global__ void shx_scy(TGrid grid_2d_i, rVector<T> M_i, T fx, T fy, T bg, TGrid grid_2d_o, rVector<T> M_o)
		{
			int iy = threadIdx.x + blockIdx.x*blockDim.x;
			int ix = threadIdx.y + blockIdx.y*blockDim.y;

			if((ix < grid_2d_o.nx)&&(iy < grid_2d_o.ny))
			{
				host_device_detail::shx_scy(ix, iy, grid_2d_i, M_i, fx, fy, bg, grid_2d_o, M_o);
			}
		}

	} // namespace device_detail

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template <class TVector>
	enable_if_device_vector<TVector, void>
	fill(Stream<e_device> &stream, TVector &M_io, Value_type<TVector> value_i)
	{
		thrust::fill(thrust::device, M_io.begin(), M_io.end(), value_i);
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::scale<value_type>(w_i));
	}

	template <class TVector>
	enable_if_device_vector<TVector, void>
	scale(Stream<e_device> &stream, Value_type<TVector> w_i, TVector &M_io)
	{
		scale(stream, w_i, M_io, M_io);
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	square(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square<value_type>());
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	square_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square_scale<value_type>(w_i));
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add<value_type>());
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		add(stream, M_i, M_io, M_io);
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_device> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, 
	Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_scale_i<value_type>(w1_i, w2_i));
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_scale<value_type>(w_i));
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_i<value_type>());
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_square<value_type>());
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_device> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_scale_i<value_type>(w1_i, w2_i));
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_device> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), functor::add_square_scale<value_type>(w_i));
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_device> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::multiply<value_type>());
	}

	template <class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_device> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		multiply(stream, M_i, M_io, M_io);
	}

	template <class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	sum(Stream<e_device> &stream, TVector &M_i)
	{
		return thrust::reduce(M_i.begin(), M_i.end());
	}

	template <class TVector>
	enable_if_device_vector<TVector, Value_type_r<TVector>>
	sum_square(Stream<e_device> &stream, TVector &M_i)
	{
		using T_r = Value_type_r<TVector>;
		return thrust::transform_reduce(M_i.begin(), M_i.end(), 
		functor::square<T_r>(), T_r(0), functor::add<T_r>());
	}

	template <class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	mean(Stream<e_device> &stream, TVector &M_i)
	{
		return sum(stream, M_i)/M_i.size();
	}

	template <class TVector>
	enable_if_device_vector<TVector, void>
	mean_var(Stream<e_device> &stream, TVector &M_i, Value_type<TVector> &x_mean, Value_type_r<TVector> &x_var)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		x_mean = mean(stream, M_i);
		x_var = thrust::transform_reduce(M_i.begin(), M_i.begin(), functor::square_dif<T, T_r>(x_mean), T_r(0), functor::add<T_r>());

		x_var = x_var/M_i.size();
	}

	template <class TVector>
	enable_if_device_vector<TVector, Value_type_r<TVector>>
	variance(Stream<e_device> &stream, TVector &M_i)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		T x_mean;
		T_r x_var;
		mean_var(stream, M_i, x_mean, x_var);
		return x_var;
	}

	/***********************************************************************/
	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_r_factor_1d(TGrid &grid_1d, Value_type<TGrid> gx, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_1d.cuda_grid();

		device_detail::exp_r_factor_1d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d, gx, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_r, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_r_factor_2d_bc(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> alpha, 
	TVector_r &gy, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::exp_r_factor_2d_bc<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, alpha, gy, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_r_factor_2d(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> gx, Value_type<TGrid> gy, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::exp_r_factor_2d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, gx, gy, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	mul_exp_r_factor_2d(Stream<e_device> &stream, TGrid &grid_2d, Vector<Value_type<TGrid>, e_host> &gxh, 
	Vector<Value_type<TGrid>, e_host> &gyh, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		Vector<Value_type<TGrid>, e_device> gx = gxh;
		Vector<Value_type<TGrid>, e_device> gy = gyh;
		device_detail::mul_exp_r_factor_2d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, gx, gy, fPsi_i, fPsi_o);
	}

	/***********************************************************************/
	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_g_factor_1d(TGrid &grid_1d, Value_type<TGrid> x, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_1d.cuda_grid();

		device_detail::exp_g_factor_1d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d, x, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_r, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_g_factor_2d_bc(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> alpha, 
	TVector_r &y, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::exp_g_factor_2d_bc<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, alpha, y, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	exp_g_factor_2d(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::exp_g_factor_2d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, x, y, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	mul_exp_g_factor_2d(Stream<e_device> &stream, TGrid &grid_2d, Vector<Value_type<TGrid>, e_host> &xh, 
	Vector<Value_type<TGrid>, e_host> &yh, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		Vector<Value_type<TGrid>, e_device> x = xh;
		Vector<Value_type<TGrid>, e_device> y = yh;
		device_detail::mul_exp_g_factor_2d<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, x, y, fPsi_i, fPsi_o);
	}

	/***********************************************************************/
	template <class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, Value_type<TVector_r>>
	atom_cost_function(TGrid &grid_2d, const Atom_Sa<Value_type<TGrid>> &atom_Ip, TVector_r &M_i)
	{
		TVector_r sum_v(1);

		Grid_BT grid_bt;
		grid_bt.Blk = dim3();
		grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);

		device_detail::atom_cost_function<typename TVector_r::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, atom_Ip, M_i, sum_v);
		return sum_v[0];
	}

	template <class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, void>
	subtract_atom(Stream<e_device> &stream, TGrid &grid_2d, Vector<Atom_Sa<Value_type<TGrid>>, e_host> &atom_Ip, TVector_r &M_i)
	{
		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			Grid_BT grid_bt;
			grid_bt.Blk = dim3();
			grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);

			device_detail::subtract_atom<typename TVector_r::value_type><<<grid_bt.Blk, grid_bt.Thr, 0, stream[istream]>>>(grid_2d, atom_Ip[istream], M_i);
		}
	}

	// Linear projected potential: V and zV
	template <class TQ1, class TVAtom>
	enable_if_device<TQ1, void>
	linear_Vz(Stream<e_device> &stream, ePotential_Type potential_type, TQ1 &qz, TVAtom &vatom)
	{	
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		auto str_linear_Vz = [](cudaStream_t &stream, const ePotential_Type &potential_type, TQ1 &qz, TAtom &atom)
		{
			if(atom.charge == 0)
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						device_detail::linear_Vz<ePT_Doyle_0_4, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Peng_0_4:
						device_detail::linear_Vz<ePT_Peng_0_4, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Peng_0_12:
						device_detail::linear_Vz<ePT_Peng_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Kirkland_0_12:
						device_detail::linear_Vz<ePT_Kirkland_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Weickenmeier_0_12:
						device_detail::linear_Vz<ePT_Weickenmeier_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Lobato_0_12:
						device_detail::linear_Vz<ePT_Lobato_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
				}
			}
			else
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						device_detail::linear_Vz<ePT_Doyle_0_4, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Peng_0_4:
						device_detail::linear_Vz<ePT_Peng_0_4, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Peng_0_12:
						device_detail::linear_Vz<ePT_Peng_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Kirkland_0_12:
						device_detail::linear_Vz<ePT_Kirkland_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Weickenmeier_0_12:
						device_detail::linear_Vz<ePT_Weickenmeier_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
					case ePT_Lobato_0_12:
						device_detail::linear_Vz<ePT_Lobato_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
						break;
				}		
			}
		};

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			str_linear_Vz(stream[istream], potential_type, qz, vatom[istream]);
		}
	}

	// Get Local interpolation coefficients
	template <class TVAtom> 
	enable_if_device<typename TVAtom::value_type, void>
	cubic_poly_coef(Stream<e_device> &stream, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			device_detail::cubic_poly_coef<TAtom><<<dim3(1), dim3(c_nR), 0, stream[istream]>>>(vatom[istream]);
		}
	}

	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, void>
	fft1_shift(TGrid &grid_1d, TVector &M_io)
	{
		auto grid_bt = grid_1d.cuda_grid_h();

		device_detail::fft1_shift<TGrid, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d, M_io); 	
	}

	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, void>
	fft2_sft_bc(Stream<e_device> &stream, TGrid &grid_2d, TVector &M_io)
	{
		auto grid_bt = grid_2d.cuda_grid(dim3((grid_2d.nyh+c_thrnxny-1)/c_thrnxny, (grid_2d.nx+c_thrnxny-1)/c_thrnxny));

		device_detail::fft2_sft_bc<TGrid, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_io); 	
	}

	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, void>
	fft2_shift(Stream<e_device> &stream, TGrid &grid_2d, TVector &M_io)
	{
		auto grid_bt = grid_2d.cuda_grid_h();

		device_detail::fft2_shift<TGrid, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_io); 	
	}

	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	sum_over_Det(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		TVector sum_v(c_thrnxny*c_thrnxny);

		auto grid_bt = grid_2d.cuda_grid(dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_over_Det<TGrid, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, pow(g_min, 2), pow(g_max, 2), M_i, sum_v);
		return mt::sum(stream, sum_v);
	}

	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny, 0);

		auto grid_bt = grid_2d.cuda_grid(dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, pow(g_min, 2), pow(g_max, 2), M_i, sum_t);
		return sum(stream, sum_t);
	}

	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_device> &stream, TGrid &grid_2d, TVector_1 &S_i, TVector_2 &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny, 0);

		auto grid_bt = grid_2d.cuda_grid(dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, typename TVector_2::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, S_i, M_i, sum_t);
		return sum(stream, sum_t);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	bandwidth_limit(Stream<e_device> &stream, TGrid &grid_2d, TVector_c &M_io)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::bandwidth_limit<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_io); 
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	hard_aperture(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> g_max, Value_type<TGrid> w, TVector_c &M_io)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::hard_aperture<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, pow(g_max, 2), w, M_io); 
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	propagate(Stream<e_device> &stream, TGrid &grid_2d, Value_type<TGrid> w, 
	Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();

		device_detail::propagate<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, w, gxu, gyu, psi_i, psi_o);
	}

	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_device_vector_and_device_vector<TVector_1, TVector_2, void>
	transmission_function(Stream<e_device> &stream, TGrid &grid_2d, eElec_Spec_Int_Model elec_spec_int_model, 
	Value_type<TGrid> w, TVector_1 &V0_i, TVector_2 &Trans_o)
	{	
		using T_r = Value_type<TGrid>;

		thrust::transform(V0_i.begin(), V0_i.end(), Trans_o.begin(), 
		functor::transmission_function<T_r>(w, elec_spec_int_model));
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	probe(Stream<e_device> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, 
	Value_type<TGrid> y, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_o)
	{
		using T_r = Value_type<TGrid>;

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::probe<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, lens, x, y, gxu, gyu, fPsi_o);

		auto total = sum_square(stream, fPsi_o);
		scale(stream, sqrt(1.0/total), fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_CTF(Stream<e_device> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();
		device_detail::apply_CTF<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, lens, gxu, gyu, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_PCTF(Stream<e_device> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto grid_bt = grid_2d.cuda_grid();
		device_detail::apply_PCTF<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, lens, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_xyz(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_xyz<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_x, k_y, k_z);

		fft_2d.inverse(k_x);
		fft_2d.inverse(k_y);
		fft_2d.inverse(k_z);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_x(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_x)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_x<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_x);

	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_y(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_y)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_y<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_y);

		fft_2d.inverse(k_y);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_z(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_z<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_z);

		fft_2d.inverse(k_z);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mn1(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_mn1)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_mn1<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_mn1);

		fft_2d.inverse(k_mn1);
	}

	template <class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mp1(Stream<e_device> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_device> &fft_2d, TVector_c &k_mp1)
	{
		eels.factor = device_detail::Lorentz_factor(stream, grid_2d, eels);

		auto grid_bt = grid_2d.cuda_grid();

		device_detail::kernel_mp1<TGrid, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, eels, k_mp1);

		fft_2d.inverse(k_mp1);
	}

	/***************************************************************************/
	/***************************************************************************/
	template <class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector_and_host_vector<TVector_i, TVector_o>::value 
	&& is_complex<Value_type<TVector_o>>::value && !std::is_same<Value_type<TVector_i>, Value_type<TVector_o>>::value, void>::type
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o, 
	Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		// copy data from host to host
		copy_to_host(stream, *M_i_h, M_o);
	}

	template <class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector_and_host_vector<TVector_i, TVector_o>::value 
	&& (!is_complex<Value_type<TVector_o>>::value || std::is_same<Value_type<TVector_i>, Value_type<TVector_o>>::value), void>::type
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o, 
	Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template <class TVector_i, class TVector_o>
	enable_if_device_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		// add and scale
		mt::add_scale_to_host(stream, w_i, *M_i_h, M_o);
	}

	template <class TVector_i, class TVector_o>
	enable_if_device_vector_and_host_vector<TVector_i, TVector_o, void>
	add_square_scale_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		Vector<Value_type<TVector_i>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// data transfer from GPU to CPU
		M_i_h->assign(M_i.begin(), M_i.end());

		mt::add_square_scale_to_host(stream, w_i, *M_i_h, M_o);
	}

	template <class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_device_vector_and_host_vector<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i, 
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_host> *psi_i_h = nullptr)
	{
		Vector<Value_type<TVector_c_i>, e_host> M_h;
		psi_i_h = (psi_i_h == nullptr)?&M_h:psi_i_h;

		// data transfer from GPU to CPU
		psi_i_h->assign(psi_i.begin(), psi_i.end());

		mt::add_scale_m2psi_psi_to_host(stream, w_i, *psi_i_h, m2psi_o, psi_o);
	}

	/***************************************************************************/
	/***************************************************************************/

	template <class TVector>
	enable_if_device_vector<TVector, void>
	trs(Stream<e_device> &stream, const int &nrows, const int &ncols, TVector &M)
	{
		Grid_BT grid_bt;
		grid_bt.Blk = dim3((nrows+c_thrnxny-1)/c_thrnxny, (ncols+c_thrnxny-1)/c_thrnxny);
		grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);

		TVector M_t(nrows*ncols);
		device_detail::trs<typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(ncols, nrows, M, M_t);
		M = M_t;
	}


	/*******************************************************************/
	// find peak position 1d
	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	fit_max_pos_1d(TGrid &grid_1d, TVector &Im, Value_type<TVector> p_i, 
	Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		Vector<T, e_host> Im_h = Im;
		return fit_max_pos_1d(grid_1d, Im_h, p_i, sigma_i, radius_i);
	}

	// find peak position 2d
	template <class TGrid, class TVector>
	enable_if_device_vector<TVector, r2d<Value_type<TVector>>>
	fit_max_pos_2d(TGrid &grid_2d, TVector &Im, r2d<Value_type<TVector>> p_i, 
	Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		Vector<T, e_host> Im_h = Im;
		return fit_max_pos_2d(grid_2d, Im_h, p_i, sigma_i, radius_i);
	}

	inline
	bool is_gpu_available()
	{
		bool is_available = false;
		try
		{
			int device_count = 0;
			cudaError_t error_id = cudaGetDeviceCount(&device_count);

			is_available = !((error_id != cudaSuccess)||(device_count == 0));
		}
		catch(...)
		{
			is_available = false;
		}

		return is_available;
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

} // namespace mt

#endif
