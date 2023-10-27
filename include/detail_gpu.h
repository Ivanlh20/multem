/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version of the License, or
* (at your option) any later version.
*
* Multem is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

// #if defined __CUDACC__ && !defined GPU_DETAIL_H

#pragma once

#include "const_enum.h"
#include "math_mt.h"
#include "fcn_cos_tap.h"
#include "fcn_gauss.h"
#include "fcn_exp.h"
#include "fcn_fermi.h"
#include "fcn_butwth.h"
#include "fcn_hann.h"
#include "type_traits_gen.h"
#include "stream_gpu.h"
#include "detail_cgpu.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cooperative_groups.h>

namespace gpu_cg = cooperative_groups;

namespace mt
{
	/* shared memory arrays */
	namespace detail_gpu
	{
		// Utility class used to avoid linker errors with extern
		// unsized shared memory arrays with templated type
		template <class T>
		struct SharedMemory
		{
			__device__ inline operator T *()
			{
				extern __shared__ int __smem[];
				return (T *)__smem;
			}

			__device__ inline operator const T *() const
			{
				extern __shared__ int __smem[];
				return (T *)__smem;
			}
		};

		// specialize for double to avoid unaligned memory
		// access compile errors
		template<>
		struct SharedMemory<double>
		{
			__device__ inline operator double *()
			{
				extern __shared__ double __smem_d[];
				return (double *)__smem_d;
			}

			__device__ inline operator const double *() const
			{
				extern __shared__ double __smem_d[];
				return (double *)__smem_d;
			}
		};
	}

	/* macros reduce */
	namespace detail_gpu
	{
		#define GPU_BLK_REDUCE(sum, tid, sdata, bid, gdata_o)								\
		{																					\
			/* each thread puts its local sum into shared memory */							\
			sdata[tid] = sum;																\
			cta.sync();																		\
																							\
			/* do reduction in shared mem */												\
			if ((blocksize >= 512) && (tid < 256))											\
			{																				\
				sdata[tid] = sum = sum + sdata[tid + 256];									\
			}																				\
			cta.sync();																		\
																							\
			if ((blocksize >= 256) && (tid < 128))											\
			{																				\
				sdata[tid] = sum = sum + sdata[tid + 128];									\
			}																				\
			cta.sync();																		\
																							\
			if ((blocksize >= 128) && (tid < 64))											\
			{																				\
				sdata[tid] = sum = sum + sdata[tid + 64];									\
			}																				\
			cta.sync();																		\
																							\
			gpu_cg::thread_block_tile<32> tile32 = gpu_cg::tiled_partition<32>(cta);		\
																							\
			if (cta.thread_rank() < 32)														\
			{																				\
				/* Fetch final intermediate sum from 2nd warp */							\
				if (blocksize >= 64)														\
				{																			\
					sum += sdata[tid + 32];													\
				}																			\
				/* Reduce final warp using shuffle */										\
				for (int offset = tile32.size()/2; offset > 0; offset /= 2) 				\
				{																			\
						sum += tile32.shfl_down(sum, offset);									\
				}																			\
			}																				\
																							\
			if (cta.thread_rank() == 0)														\
			{																				\
				gdata_o[bid] = sum;															\
			}																				\
		}
	}

	/* atomicAdd - double */
	namespace detail_gpu
	{
	#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
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
	}

	/* remainder summation */
	namespace detail_gpu
	{
		template <class T>
		__global__ void fcn_sum_rem(pVctr_gpu_32<T> mx_i)
		{ 
			gpu_cg::thread_block cta = gpu_cg::this_thread_block();
			T *sdata = SharedMemory<T>();

			dt_uint32 tid = threadIdx.x;
			dt_uint32 blocksize = blockDim.x;

			KS<T> sum = mx_i[tid];

			GPU_BLK_REDUCE(sum, tid, sdata, 0, mx_i);
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_sum_rem_gpu(TVctr& mx_i)
		{
			using T = Value_type<TVctr>;
			int threads = mx_i.size();
			int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
			fcn_sum_rem<T><<<dim3(1), dim3(threads), smemSize>>>(mx_i.ptr_32());
			T sum = 0;
			cudaMemcpy(&sum, mx_i.data(), sizeof(T), cudaMemcpyDeviceToHost);

			return sum;
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_sum_rem_cpu(TVctr& mx_i)
		{
			using T = Value_type<TVctr>;

			Vctr_cpu<T> mx(mx_i);
			KS<T> sum(0);

			for (dt_int32 ik=0; ik<mx.m_size; ik++)
			{
				sum += mx[ik];
			}

			return sum;
		}
	}

	/* vector functions */
	namespace detail_gpu
	{
		// These were implemented using the thrust library: functors
	}

	/* add - assign - crop - norm_2 - fftsft */
	namespace detail_gpu
	{
		/* shift matrix respect to nx_h */
		template <class T>
		__global__ void fcn_fftsft_1d(iGrid_1d igrid, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_1DC(igrid.nx_h, detail_cgpu::fcn_fftsft_1d(ix, igrid, mx_io.m_data));
		}

		/* shift matrix respect to ny_h */
		template <class T>
		__global__ void fcn_fftsft_bc_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny_h, detail_cgpu::fcn_fftsft_bc_2d(ix, iy, igrid, mx_io.m_data));
		}

		/* shift matrix respect to (nx_h, ny_h) */
		template <class T>
		__global__ void fcn_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_fftsft_2d(ix, iy, igrid, mx_io.m_data));
		}

 		/* shift matrix respect to (nx_h, ny_h) */
		template <class T>
		__global__ void fcn_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_fftsft_2d(ix, iy, igrid, mx_i.m_data, mx_o.m_data));
		}

		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class T>
		__global__ void fcn_fftsft_3d(iGrid_3d igrid, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_3DC(igrid.nx_h, igrid.ny_h, igrid.nz_h, detail_cgpu::fcn_fftsft_3d(ix, iy, iz, igrid, mx_io.m_data));
		}

 		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class T>
		__global__ void fcn_fftsft_3d(iGrid_3d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_3DC(igrid.nx_h, igrid.ny_h, igrid.nz_h, detail_cgpu::fcn_fftsft_3d(ix, iy, iz, igrid, mx_i.m_data, mx_o.m_data));
		}

		/***************************************************************************************/
 		/* add, scale and shift */
		template <class T>
		__global__ void fcn_add_sc_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, T w, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_add_sc_fftsft_2d(ix, iy, igrid, mx_i.m_data, w, mx_o.m_data));
		}

 		/* add, scale, square and shift */
		template <class T, class U>
		__global__ void fcn_add_sc_norm_2_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, U w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_add_sc_norm_2_fftsft_2d(ix, iy, igrid, mx_i.m_data, w, mx_o.m_data));
		}

		/***************************************************************************************/
 		/* assign and crop */
		template <class T>
		__global__ void fcn_assign_crop_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, iRegion_Rect_2d& iregion, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_assign_crop_2d(ix, iy, igrid, mx_i.m_data, iregion, mx_o.m_data));
		}

 		/* assign, crop and shift */
		template <class T>
		__global__ void fcn_assign_crop_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, iRegion_Rect_2d& iregion, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_assign_crop_fftsft_2d(ix, iy, igrid, mx_i.m_data, iregion, mx_o.m_data));
		}
 
 		/* add, scale, crop and shift */
		template <class T>
		__global__ void fcn_add_sc_crop_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, iRegion_Rect_2d& iregion, T w, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_add_sc_crop_fftsft_2d(ix, iy, igrid, mx_i.m_data, iregion, w, mx_o.m_data));
		}

 		/* add, scale, square, crop and shift */
		template <class T, class U>
		__global__ void fcn_add_sc_norm_2_crop_fftsft_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, iRegion_Rect_2d& iregion, U w, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx_h, igrid.ny_h, detail_cgpu::fcn_add_sc_norm_2_crop_fftsft_2d(ix, iy, igrid, mx_i.m_data, iregion, w, mx_o.m_data));
		}
	}

	/* transpose - element wise matrix op vector */
	namespace detail_gpu
	{
		/* transpose */
		template <class T>
		__global__ void fcn_trs_2d(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_trs_2d(ix, iy, igrid, mx_i.m_data, mx_o.m_data));
		}

		/***************************************************************************************/
		/* element wise addition: matrix + vector row */
		template <class T>
		__global__ void fcn_ew_add_mx_vctr_row(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_add_mx_vctr_row(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}

		/* element wise addition: matrix + vector col */
		template <class T>
		__global__ void fcn_ew_add_mx_vctr_col(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_add_mx_vctr_col(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}

		/***************************************************************************************/
		/* element wise subtraction: matrix - vector row */
		template <class T>
		__global__ void fcn_ew_sub_mx_vctr_row(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_sub_mx_vctr_row(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}

		/* element wise subtraction: matrix - vector col */
		template <class T>
		__global__ void fcn_ew_sub_mx_vctr_col(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_sub_mx_vctr_col(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}

		/***************************************************************************************/
		/* element wise multiplication matrix X vector row */
		template <class T>
		__global__ void fcn_ew_mult_mx_vctr_row(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_mult_mx_vctr_row(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}

		/* element wise multiplication matrix X vector col */
		template <class T>
		__global__ void fcn_ew_mult_mx_vctr_col(iGrid_2d igrid, pVctr_gpu_32<T> vctr, pVctr_gpu_32<T> mx_io)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_ew_mult_mx_vctr_col(ix, iy, igrid, vctr.m_data, mx_io.m_data));
		}
	}

	/* aperture functions */
	namespace detail_gpu
	{
		template <class T, class U>
		__global__ void fcn_fermi_aperture(Grid_2d<T> grid, T g2_cut, T alpha, T w, pVctr_gpu_32<U> mx_io)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_fermi_aperture(ix, iy, grid, g2_cut, alpha, w, mx_io.m_data));
		}

		template <class T, class U>
		__global__ void fcn_hard_aperture(Grid_2d<T> grid, T g2_cut, T w, pVctr_gpu_32<U> mx_io)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_hard_aperture(ix, iy, grid, g2_cut, w, mx_io.m_data));
		}
	}

	/* phase shifts real space*/
	namespace detail_gpu
	{
		// phase factor 1d
		template <class T, class U>
		__global__ void fcn_rs_exp_factor_1d(Grid_1d<T> grid, pVctr_gpu_32<U> psi_i, T gx, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_1DC(grid.nx, detail_cgpu::fcn_rs_exp_factor_1d(ix, grid, psi_i.m_data, gx, w, psi_o.m_data));
		}

		// phase factor 2d by col
		template <class T, class U>
		__global__ void fcn_rs_exp_factor_2d_bc(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, T alpha, pVctr_gpu_32<T> gy, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_rs_exp_factor_2d_bc(ix, iy, grid, psi_i.m_data, alpha, gy, w, psi_o.m_data));
		}

		// phase factor 2d
		template <class T, class U>
		__global__ void fcn_rs_exp_factor_2d(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, R_2d<T> g, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_rs_exp_factor_2d(ix, iy, grid, psi_i.m_data, g, w, psi_o.m_data));
		}

		// phase factor 2d multipositions
		template <class T, class U>
		__global__ void fcn_rs_mul_exp_factor_2d(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, pVctr_gpu_32<R_2d<T>> g, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_rs_mul_exp_factor_2d(ix, iy, grid, psi_i.m_data, g.m_data, g.m_size, w, psi_o.m_data));
		}
	}

	/* phase shifts fourier space*/
	namespace detail_gpu
	{
		// phase factor 1d
		template <class T, class U>
		__global__ void fcn_fs_exp_factor_1d(Grid_1d<T> grid, pVctr_gpu_32<U> psi_i, T rx, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_1DC(grid.nx, detail_cgpu::fcn_fs_exp_factor_1d(ix, grid, psi_i.m_data, rx, w, psi_o.m_data));
		}

		// phase factor 2d by col
		template <class T, class U>
		__global__ void fcn_fs_exp_factor_2d_bc(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, T alpha, pVctr_gpu_32<T> ry, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_fs_exp_factor_2d_bc(ix, iy, grid, psi_i.m_data, alpha, ry, w, psi_o.m_data));
		}

		// phase factor 2d
		template <class T, class U>
		__global__ void fcn_fs_exp_factor_2d(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, R_2d<T> r, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_fs_exp_factor_2d(ix, iy, grid, psi_i.m_data, r, w, psi_o.m_data));
		}

		// phase factor 2d multi-positions
		template <class T, class U>
		__global__ void fcn_fs_mul_exp_factor_2d(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, pVctr_gpu_32<R_2d<T>> r, T w, pVctr_gpu_32<U> psi_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_fs_mul_exp_factor_2d(ix, iy, grid, psi_i.m_data, r.m_data, r.m_size, w, psi_o.m_data));
		}
	}

	/* gradient */
	namespace detail_gpu
	{
		template <class T>
		__global__ void fcn_grad_x(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> dmx_x)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_grad_x(ix, iy, igrid, mx_i.m_data, dmx_x.m_data));
		}

		template <class T>
		__global__ void fcn_grad_y(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> dmx_y)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_grad_y(ix, iy, igrid, mx_i.m_data, dmx_y.m_data));
		}

		template <class T>
		__global__ void fcn_grad(iGrid_2d igrid, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> dmx_x, pVctr_gpu_32<T> dmx_y)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_grad(ix, iy, igrid, mx_i.m_data, dmx_x.m_data, dmx_y.m_data));
		}
	}

	/* function multiplication fourier space */
	namespace detail_gpu
	{
		#define FCN_MULT_FS_FCN_GPU_DET(FN, FCN, POW, DIM)																						\
		template <class T, class U>																												\
		__global__ void fcn_mult_fs_##FN##_##DIM##d(Grid_##DIM##d<T> grid, FCN<T> fcn, T w, pVctr_gpu_32<U> mx_io)								\
		{																																		\
			FOR_LOOP_NDC(DIM, F_SHAPE_ND(grid, DIM), detail_cgpu::fcn_mult_fs_fcn_g##POW##_##DIM##d(IDX_ND(DIM), grid, fcn, w, mx_io.m_data))	\
		}

		/* gaussian */
		FCN_MULT_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 1);		// fcn_mult_fs_gauss_1d
		FCN_MULT_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 2);		// fcn_mult_fs_gauss_2d
		FCN_MULT_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 3);		// fcn_mult_fs_gauss_3d

		/* exponential */
		FCN_MULT_FS_FCN_GPU_DET(exp, Fcn_Exp, , 1);				// fcn_mult_fs_exp_1d
		FCN_MULT_FS_FCN_GPU_DET(exp, Fcn_Exp, , 2);				// fcn_mult_fs_exp_2d
		FCN_MULT_FS_FCN_GPU_DET(exp, Fcn_Exp, , 3);				// fcn_mult_fs_exp_3d

		/* fermi */
		FCN_MULT_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 1);		// fcn_mult_fs_fermi_1d
		FCN_MULT_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 2);		// fcn_mult_fs_fermi_2d
		FCN_MULT_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 3);		// fcn_mult_fs_fermi_3d

		/* butterworth */
		FCN_MULT_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 1);		// fcn_mult_fs_butwth_1d
		FCN_MULT_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 2);		// fcn_mult_fs_butwth_2d
		FCN_MULT_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 3);		// fcn_mult_fs_butwth_3d

		/* hann */
		FCN_MULT_FS_FCN_GPU_DET(hann, Fcn_Hann, , 1);			// fcn_mult_fs_hann_1d
		FCN_MULT_FS_FCN_GPU_DET(hann, Fcn_Hann, , 2);			// fcn_mult_fs_hann_2d
		FCN_MULT_FS_FCN_GPU_DET(hann, Fcn_Hann, , 3);			// fcn_mult_fs_hann_3d
	}

	/* deconvolution */
	namespace detail_gpu
	{
		#define FCN_DCV_FS_FCN_GPU_DET(FN, FCN, POW, DIM)																								\
		template <class T, class U>																														\
		__global__ void fcn_dcv_fs_##FN##_##DIM##d(Grid_##DIM##d<T> grid, FCN<T> fcn, T psnr, T w, pVctr_gpu_32<U> mx_io)								\
		{																																				\
			FOR_LOOP_NDC(DIM, F_SHAPE_ND(grid, DIM), detail_cgpu::fcn_dcv_fs_fcn_g##POW##_##DIM##d(IDX_ND(DIM), grid, fcn, psnr, w, mx_io.m_data))		\
		}

		/* gaussian */
		FCN_DCV_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 1);			// fcn_dcv_fs_gauss_1d
		FCN_DCV_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 2);			// fcn_dcv_fs_gauss_2d
		FCN_DCV_FS_FCN_GPU_DET(gauss, Fcn_Gauss, 2, 3);			// fcn_dcv_fs_gauss_3d

		/* exponential */
		FCN_DCV_FS_FCN_GPU_DET(exp, Fcn_Exp, , 1);				// fcn_dcv_fs_exp_1d
		FCN_DCV_FS_FCN_GPU_DET(exp, Fcn_Exp, , 2);				// fcn_dcv_fs_exp_2d
		FCN_DCV_FS_FCN_GPU_DET(exp, Fcn_Exp, , 3);				// fcn_dcv_fs_exp_3d

		/* fermi */
		FCN_DCV_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 1);			// fcn_dcv_fs_fermi_1d
		FCN_DCV_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 2);			// fcn_dcv_fs_fermi_2d
		FCN_DCV_FS_FCN_GPU_DET(fermi, Fcn_Fermi, 2, 3);			// fcn_dcv_fs_fermi_3d

		/* butterworth */
		FCN_DCV_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 1);		// fcn_dcv_fs_butwth_1d
		FCN_DCV_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 2);		// fcn_dcv_fs_butwth_2d
		FCN_DCV_FS_FCN_GPU_DET(butwth, Fcn_Butwth, 2, 3);		// fcn_dcv_fs_butwth_3d

		/* hann */
		FCN_DCV_FS_FCN_GPU_DET(hann, Fcn_Hann, , 1);			// fcn_dcv_fs_hann_1d
		FCN_DCV_FS_FCN_GPU_DET(hann, Fcn_Hann, , 2);			// fcn_dcv_fs_hann_2d
		FCN_DCV_FS_FCN_GPU_DET(hann, Fcn_Hann, , 3);			// fcn_dcv_fs_hann_3d
	}

	/* window functions */
	namespace detail_gpu
	{
		#define FCN_WD_FCN_GPU_DET(FN, FCN, DIM)																									\
		template <dt_bool sft, class T, class U>																									\
		__global__ void fcn_wd_##FN##_##DIM##d(Grid_##DIM##d<T> grid, FCN<T> fcn, pVctr_gpu_32<U> mx_o)												\
		{																																			\
			if (sft)																																\
			{																																		\
				FOR_LOOP_NDC(DIM, F_SHAPE_ND(grid, DIM), detail_cgpu::fcn_wd_fcn_r2_sft_##DIM##d(IDX_ND(DIM), grid, fcn, mx_o.m_data))				\
			}																																		\
			else																																	\
			{																																		\
				FOR_LOOP_NDC(DIM, F_SHAPE_ND(grid, DIM), detail_cgpu::fcn_wd_fcn_r2_##DIM##d(IDX_ND(DIM), grid, fcn, mx_o.m_data))					\
			}																																		\
		}

			FCN_WD_FCN_GPU_DET(gauss, Wd_Gauss_1d, 1);		// fcn_wd_gauss_1d
			FCN_WD_FCN_GPU_DET(gauss, Wd_Gauss_2d, 2);		// fcn_wd_gauss_2d
			FCN_WD_FCN_GPU_DET(gauss, Wd_Gauss_3d, 3);		// fcn_wd_gauss_3d

			FCN_WD_FCN_GPU_DET(exp, Wd_Exp_1d, 1);			// fcn_wd_exp_1d
			FCN_WD_FCN_GPU_DET(exp, Wd_Exp_2d, 2);			// fcn_wd_exp_2d
			FCN_WD_FCN_GPU_DET(exp, Wd_Exp_3d, 3);			// fcn_wd_exp_3d

			FCN_WD_FCN_GPU_DET(fermi, Wd_Fermi_1d, 1);		// fcn_wd_fermi_1d
			FCN_WD_FCN_GPU_DET(fermi, Wd_Fermi_2d, 2);		// fcn_wd_fermi_2d
			FCN_WD_FCN_GPU_DET(fermi, Wd_Fermi_3d, 3);		// fcn_wd_fermi_3d

			FCN_WD_FCN_GPU_DET(butwth, Wd_Butwth_1d, 1);	// fcn_wd_butwth_1d
			FCN_WD_FCN_GPU_DET(butwth, Wd_Butwth_2d, 2);	// fcn_wd_butwth_2d
			FCN_WD_FCN_GPU_DET(butwth, Wd_Butwth_3d, 3);	// fcn_wd_butwth_3d

			FCN_WD_FCN_GPU_DET(hann, Wd_Hann_1d, 1);		// fcn_wd_hann_1d
			FCN_WD_FCN_GPU_DET(hann, Wd_Hann_2d, 2);		// fcn_wd_hann_2d
			FCN_WD_FCN_GPU_DET(hann, Wd_Hann_3d, 3);		// fcn_wd_hann_3d
	}

	/* phase correlation */
	namespace detail_gpu
	{
		/****************** pcf data processing real space *******************/
		template <class T, class U>
		__global__ void fcn_rs_pcf_1d_dp(Grid_1d<T> grid, pVctr_gpu_32<T> mx_i, 
			Wd_Butwth_1d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_1DC(grid.nx, detail_cgpu::fcn_rs_pcf_1d_dp(ix, grid, mx_i.m_data, wd, w, mx_o.m_data));
		}

		template <class T, class U>
		__global__ void fcn_rs_pcf_2d_dp(Grid_2d<T> grid, pVctr_gpu_32<T> mx_i, 
			Wd_Butwth_2d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_rs_pcf_2d_dp(ix, iy, grid, mx_i.m_data, wd, w, mx_o.m_data));
		}

		template <class T, class U>
		__global__ void fcn_rs_pcf_3d_dp(Grid_3d<T> grid, pVctr_gpu_32<T> mx_i, 
			Wd_Butwth_3d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_3DC(grid.nx, grid.ny, grid.nz, detail_cgpu::fcn_rs_pcf_2d_dp(ix, iy, iz, grid, mx_i.m_data, wd, w, mx_o.m_data));
		}

		/***************** pcf data processing fourier space *****************/
		template <class T, class U>
		__global__ void fcn_fs_pcf_1d_dp(Grid_1d<T> grid, pVctr_gpu_32<T> mx_1i, pVctr_gpu_32<T> mx_2i, 
			Wd_Gauss_1d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_1DC(grid.nx, detail_cgpu::fcn_fs_pcf_1d_dp(ix, grid, mx_1i.m_data, mx_2i.m_data, wd, w, mx_o.m_data));
		}

		template <class T, class U>
		__global__ void fcn_fs_pcf_2d_dp(Grid_2d<T> grid, pVctr_gpu_32<T> mx_1i, pVctr_gpu_32<T> mx_2i, 
			Wd_Gauss_2d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_2DC(grid.nx, grid.ny, detail_cgpu::fcn_fs_pcf_2d_dp(ix, iy, grid, mx_1i.m_data, mx_2i.m_data, wd, w, mx_o.m_data));
		}

		template <class T, class U>
		__global__ void fcn_fs_pcf_3d_dp(Grid_3d<T> grid, pVctr_gpu_32<T> mx_1i, pVctr_gpu_32<T> mx_2i, 
			Wd_Gauss_3d<T> wd, T w, pVctr_gpu_32<U> mx_o)
		{
			FOR_LOOP_3DC(grid.nx, grid.ny, grid.nz, detail_cgpu::fcn_fs_pcf_2d_dp(ix, iy, iz, grid, mx_1i.m_data, mx_2i.m_data, wd, w, mx_o.m_data));
		}
	}

	/* optical flow */
	namespace detail_gpu
	{
		template <class T>
		__global__ void fcn_opt_flow(iGrid_2d igrid, pVctr_gpu_32<T> mx_s, pVctr_gpu_32<T> mx_m, 
			T alpha, pVctr_gpu_32<T> dphi_x, pVctr_gpu_32<T> dphi_y)
		{
			FOR_LOOP_2DC(igrid.nx, igrid.ny, detail_cgpu::fcn_opt_flow(ix, iy, igrid, mx_s.m_data, mx_m.m_data, alpha, dphi_x.m_data, dphi_y.m_data));
		}
	}

	/* bilinear interpolation */
	namespace detail_gpu
	{
		template <class T>
		__global__ void fcn_intrpl_bl_rg_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, 
			pVctr_gpu_32<T> vrx, pVctr_gpu_32<T> vry, T bg, pVctr_gpu_32<T> mx_o)
		{
			FOR_LOOP_1DC(vrx.m_size, detail_cgpu::fcn_intrpl_bl_rg_2d(ix, grid_i, mx_i.m_data, vrx.m_data, vry.m_data, bg, mx_o.m_data));
		}
	}

	// namespace detail_gpu
	// {
// 		// scan noise xy
	// 	template <class T>
	// 	__global__ void sd_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> dx_i, pVctr_gpu_32<T> dy_i, T bg, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		 (dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_i.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_i.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::sd_2d(ix, iy, grid_i, mx_i, dx_i, dy_i, bg, mx_o);
	// 			}
	// 		}
	// 	}

	// 	template <class T>
	// 	__global__ void sd_nr_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, pVctr_gpu_32<T> dx_i, pVctr_gpu_32<T> dy_i, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_i.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_i.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::sd_nr_2d(ix, iy, grid_i, mx_i, dx_i, dy_i, mx_o);
	// 			}
	// 		}
	// 	}

	// 	// scaling xy
	// 	template <class T>
	// 	__global__ void sc_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, T fxy, Grid_2d<T> grid_o, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_o.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_o.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::sc_2d(ix, iy, grid_i, mx_i, fxy, grid_o, mx_o);
	// 			}
	// 		}
	// 	}

	// 	// rotate, scale and shift xy
	// 	template <class T>
	// 	__global__ void rot_sca_sft_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, T theta, R_2d<T> p0, 
	// 	T fx, T fy, R_2d<T> ps, T bg, Grid_2d<T> grid_o, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_o.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_o.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::rot_sca_sft_2d(ix, iy, grid_i, mx_i, theta, p0, fx, fy, ps, bg, grid_o, mx_o);
	// 			}
	// 		}
	// 	}

	// 	// 2d general affine transformation
	// 	template <class T>
	// 	__global__ void at_2d(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, Mx_2x2<T> A, R_2d<T> txy, T bg, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_i.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_i.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::at_2d(ix, iy, grid_i, mx_i, A, txy, bg, mx_o);
	// 			}
	// 		}
	// 	}

// 		// x-shear and y-scaling
	// 	template <class T>
	// 	__global__ void shx_scy(Grid_2d<T> grid_i, pVctr_gpu_32<T> mx_i, T fx, T fy, T bg, Grid_2d<T> grid_o, pVctr_gpu_32<T> mx_o)
	// 	{
	// 		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < grid_o.nx; ix += blockDim.y*gridDim.y)
	// 		{
	// 			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < grid_o.ny; iy += blockDim.x*gridDim.x)
	// 			{
	// 				detail_cgpu::shx_scy(ix, iy, grid_i, mx_i, fx, fy, bg, grid_o, mx_o);
	// 			}
	// 		}
	// 	}

	// }
}
