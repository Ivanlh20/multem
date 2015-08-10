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
#include "traits.cuh"
#include "fft2.cuh"
#include "stream.cuh"

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
				if (R2 < atom_Vp.R_max2)
				{
					int ixy;
					T V = host_device_detail::eval_cubic_poly(ix, iy, R2, grid, atom_Vp, ixy);

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
					value_type_r g2 = grid.g2_shift(ix, iy);
					if((g2_min <= g2)&&(g2 <= g2_max))
					{
						int ixy = grid.ind_col(ix, iy); 
						sum += M_i[ixy];
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}
			reduce_array_256(tid, sum, Mshare);

			if(tid == 0 )
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
					value_type_r g2 = grid.g2_shift(ix, iy);
					if((g2_min <= g2)&&(g2 <= g2_max))
					{
						int ixy = grid.ind_col(ix, iy);
						sum += thrust::norm(M_i[ixy]);
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0 )
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

		// Probe in Fourier space
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
					T g2 = grid.g2_shift(ix, iy);
					if(g2 < gc2)
					{
						sum += 1.0/(g2 + ge2);
					}
					iy += blockDim.y*gridDim.y;
				}
				ix += blockDim.x*gridDim.x;
			}

			reduce_array_256(tid, sum, Mshare);

			if(tid == 0 )
			{
				Mp_o[bid] = sum;
			}
		}

		template<class TGrid>
		Value_type<TGrid> Lorentz_factor(TGrid &grid, EELS<Value_type<TGrid>> &eels)
		{
			device_vector<Value_type<TGrid>> sum_v(c_thrnxny*c_thrnxny);

			auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

			norm_factor_lorentz<Value_type<TGrid>><<<gridBT.Blk, gridBT.Thr>>>(grid, eels.gc2, eels.ge2, sum_v);

			return sqrt(eels.occ)/multem::sum(grid, sum_v);
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
	enable_if_device_vector<TVector_r, void>
	eval_cubic_poly(TGrid &grid, Stream<Value_type<TGrid>, e_Device> &stream, 
	Vector<Atom_Vp<Value_type<TGrid>>, e_Host> &atom_Vp, TVector_r &V_0)
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
	fft2_shift(const TGrid &grid, TVector &M_io)
	{
		GridBT gridBT;
		gridBT.Blk = dim3((grid.nyh+c_thrnxny-1)/c_thrnxny, (grid.nxh+c_thrnxny-1)/c_thrnxny);
		gridBT.Thr = dim3(c_thrnxny, c_thrnxny);

		device_detail::fft2_shift<TGrid, typename TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, M_io); 	
	}

	template<class TGrid, class TVector>
	enable_if_device_vector<TVector, Value_type<TVector>>
	sum_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector &M_i)
	{
		TVector sum_v(c_thrnxny*c_thrnxny);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_over_Det<TGrid, typename TVector::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_v);

		return multem::sum(grid, sum_v);
	}

	template<class TGrid, class TVector_r>
	enable_if_device_vector<TVector_r, Value_type<TGrid>>
	sum_square_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector_r &M_i)
	{
		device_vector<Value_type<TGrid>> sum_t(c_thrnxny*c_thrnxny);

		auto gridBT = device_detail::get_grid_nxny(grid, dim3(c_thrnxny, c_thrnxny));

		device_detail::sum_square_over_Det<TGrid, typename TVector_r::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), M_i, sum_t);
		return sum(grid, sum_t);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	bandwidth_limit(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, Value_type<TVector_c> w, TVector_c &M_io)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::bandwidth_limit<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, pow(g_min, 2), pow(g_max, 2), w, M_io); 
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	phase_multiplication(const TGrid &grid, TVector_c &exp_x_i, TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::phase_multiplication<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, exp_x_i, exp_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	propagator_multiplication(const TGrid &grid, TVector_c &prop_x_i, TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::propagator_multiplication<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, prop_x_i, prop_y_i, psi_i, psi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	probe(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_o)
	{
		using value_type_r = Value_type<TGrid>;

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::probe<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, x, y, fPsi_o);

		value_type_r total = sum_square(grid, fPsi_o);
		scale(fPsi_o, sqrt(value_type_r(grid.nxy())/total));
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_CTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_CTF<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, gxu, gyu, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	apply_PCTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto gridBT = device_detail::get_grid_nxny(grid);
		device_detail::apply_PCTF<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, lens, fPsi_i, fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_xyz(TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_xyz<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x, k_y, k_z);

		fft2.inverse(k_x);
		fft2.inverse(k_y);
		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_x(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_x)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_x<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_x);

		fft2.inverse(k_x);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_y(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_y)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_y<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_y);

		fft2.inverse(k_y);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_z(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_z)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_z<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_z);

		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mn1(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_mn1)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mn1<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mn1);

		fft2.inverse(k_mn1);
	}

	template<class TGrid, class TVector_c>
	enable_if_device_vector<TVector_c, void>
	kernel_mp1(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Device> &fft2, TVector_c &k_mp1)
	{
		eels.factor = device_detail::Lorentz_factor(grid, eels);

		auto gridBT = device_detail::get_grid_nxny(grid);

		device_detail::kernel_mp1<TGrid, typename TVector_c::value_type><<<gridBT.Blk, gridBT.Thr>>>(grid, eels, k_mp1);

		fft2.inverse(k_mp1);
	}

	template<class TVector_o>
	enable_if_device_vector<TVector_o, void>
	assign(rmatrix_c &M_i, TVector_o &M_o, Vector<Value_type<TVector_o>, e_Host> *M_i_h=nullptr)
	{
		Vector<Value_type<TVector_o>, e_Host> M_h;
		if(M_i_h==nullptr)
		{
			M_i_h = &M_h;
		}

		M_i_h->resize(M_i.size);

		for(auto ixy = 0; ixy < M_i_h->size(); ixy++)
		{
			(*M_i_h)[ixy] = M_i[ixy];
		}
		multem::assign(*M_i_h, M_o);
	}

	/***************************************************************************/
	/***************************************************************************/

	template<class TGrid, class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector<TVector_i>::value && !is_rmatrix_c<TVector_o>::value, void>::type
	copy_to_host(const TGrid &grid, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		thrust::copy(M_i.begin(), M_i.end(), M_o.begin());
	}

	template<class TGrid, class TVector_i, class TVector_o>
	typename std::enable_if<is_device_vector<TVector_i>::value && is_rmatrix_c<TVector_o>::value, void>::type
	copy_to_host(const TGrid &grid, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		Vector<Value_type<TVector_i>, e_Host> M_h;
		if(M_i_h==nullptr)
		{
			M_i_h = &M_h;
		}

		M_i_h->assign(M_i.begin(), M_i.end());
		for(auto ixy = 0; ixy < M_i_h->size(); ixy++)
		{
			M_o[ixy] = (*M_i_h)[ixy];
		}
	}

	template<class TGrid, class TVector_i, class TVector_o>
	enable_if_device_vector<TVector_i, void>
	add_scale_to_host(TGrid &grid, Value_type<TVector_i> w_i, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		Vector<Value_type<TVector_i>, e_Host> M_h;
		if(M_i_h==nullptr)
		{
			M_i_h = &M_h;
		}

		M_i_h->assign(M_i.begin(), M_i.end());

		for(auto ixy = 0; ixy < M_i_h->size(); ixy++)
		{
			M_o[ixy] += w_i*((*M_i_h)[ixy]);
		}
	}

	template<class TGrid, class TVector_i, class TVector_o>
	enable_if_device_vector<TVector_i, void>
	add_square_scale_to_host(TGrid &grid, Value_type<TGrid> w_i, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		Vector<Value_type<TVector_i>, e_Host> M_h;
		if(M_i_h==nullptr)
		{
			M_i_h = &M_h;
		}

		M_i_h->assign(M_i.begin(), M_i.end());

		multem::add_square_scale(w_i, *M_i_h, M_o);
	}

	template<class TGrid, class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_device_vector<TVector_c_i, void>
	add_scale_m2psi_psi(TGrid &grid, Value_type<TGrid> w_i, TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_Host> *psi_i_h=nullptr)
	{
		Vector<Value_type<TVector_c_i>, e_Host> M_h;
		if(psi_i_h==nullptr)
		{
			psi_i_h = &M_h;
		}

		psi_i_h->assign(psi_i.begin(), psi_i.end());

		Value_type<TVector_c_i> w_i_c = w_i;
		for(auto i = 0; i < psi_i_h->size(); i++)
		{
			m2psi_o[i] += w_i*thrust::norm((*psi_i_h)[i]);
			psi_o[i] += w_i_c*((*psi_i_h)[i]);
		}
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
} // namespace multem

#endif