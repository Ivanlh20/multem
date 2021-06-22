/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

// #if defined __CUDACC__ && !defined GPU_DETAIL_MT_H
#ifndef GPU_DETAIL_MT_H
	#define GPU_DETAIL_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "type_traits_mt.cuh"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	#include "quad_data.cuh"
	#include "cgpu_detail_mt.cuh"
	#include "gpu_detail.cuh"

	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cufft.h>
	#include <cooperative_groups.h>

	namespace gpu_cg = cooperative_groups;

	namespace mt
	{
		/* atomic functions */
		namespace gpu_detail_mt
		{
			/************************** pFcn_clnl_3_x<T> ****************************/
			// evaluate function
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, pLNL_Coef_gpu<T> coef, pVctr_gpu_32<T> fx, pFcn_clnl_3_1<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fx[ix] = fcn(x[ix], coef.cl, coef.cnl);
				}
			}

			// evaluate two functions
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, pLNL_Coef_gpu<T> coef, pVctr_gpu_32<T> fx1, pVctr_gpu_32<T> fx2, pFcn_clnl_3_2<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fcn(x[ix], coef.cl, coef.cnl, fx1[ix], fx2[ix]);
				}
			}

			/************************** pFcn_clnl_4_x<T> ****************************/
			// evaluate function
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, T c, pLNL_Coef_gpu<T> coef, pVctr_gpu_32<T> fx, pFcn_clnl_4_1<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fx[ix] = fcn(x[ix], c, coef.cl, coef.cnl);
				}
			}

			// evaluate two functions
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, T c, pLNL_Coef_gpu<T> coef, pVctr_gpu_32<T> fx1, pVctr_gpu_32<T> fx2, pFcn_clnl_4_2<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fcn(x[ix], c, coef.cl, coef.cnl, fx1[ix], fx2[ix]);
				}
			}

			/************************** pFcn_clnl_6_x<T> ****************************/
			// evaluate function
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, pLNL_Coef_gpu<T> coef, pQuad_Coef_1d_gpu<T> quad, pVctr_gpu_32<T> fx, pFcn_clnl_6_1<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fx[ix] = fcn(x[ix], coef.cl, coef.cnl, quad.m_size, quad.x, quad.w);
				}
			}

			// evaluate two functions
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, T z_0, T z_e, pLNL_Coef_gpu<T> coef, pQuad_Coef_1d_gpu<T> quad, pVctr_gpu_32<T> fx1, pVctr_gpu_32<T> fx2, pFcn_clnl_6_2<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fcn(x[ix], coef.cl, coef.cnl, quad.m_size, quad.x, quad.w, fx1[ix], fx2[ix]);
				}
			}

			/************************** pFcn_clnl_8_x<T> ****************************/
			// evaluate function
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, T z_0, T z_e, pLNL_Coef_gpu<T> coef, pQuad_Coef_1d_gpu<T> quad, pVctr_gpu_32<T> fx, pFcn_clnl_8_1<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fx[ix] = fcn(x[ix], z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w);
				}
			}

			// evaluate two functions
			template <class T>
			__global__ void fcn_eval_fcn_coef_lnl(pVctr_gpu_32<T> x, T z_0, T z_e, pLNL_Coef_gpu<T> coef, pQuad_Coef_1d_gpu<T> quad, pVctr_gpu_32<T> fx1, pVctr_gpu_32<T> fx2, pFcn_clnl_8_2<T> fcn)
			{
				FOR_IX_1DC(x.m_size)
				{
					fcn(x[ix], z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w, fx1[ix], fx2[ix]);
				}
			}
		} // gpu_detail

		/* detector integration */
		namespace gpu_detail_mt
		{
			/* integration over a detector ring */
			template <class T, class U>
			__global__ void fcn_int_det_ring(Grid_2d<T> grid, T g2_min, T g2_max, pVctr_gpu_32<U> mx_i, pVctr_gpu_32<U> mx_o)
			{ 
				gpu_cg::thread_block cta = gpu_cg::this_thread_block();
				U *sdata = gpu_detail::SharedMemory<U>();

				dt_uint32 tid = threadIdx.x + threadIdx.y*blockDim.x;
				dt_uint32 bid = blockIdx.x + blockIdx.y*gridDim.x;
				dt_uint32 blocksize = blockDim.x*blockDim.y;

				KS<U> sum = U(0);

				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_int_det_ring(ix, iy, grid, g2_min, g2_max, mx_i.m_data, sum));

				GPU_BLK_REDUCE(sum, tid, sdata, bid, mx_o);
			}

			/* norm_2 integration over a detector ring */
			template <class T, class U>
			__global__ void fcn_int_det_ring_norm_2(Grid_2d<T> grid, T g2_min, T g2_max, pVctr_gpu_32<U> mx_i, pVctr_gpu_32<Value_type_r<U>> mx_o)
			{ 
				using U_r = Value_type_r<U>;
				gpu_cg::thread_block cta = gpu_cg::this_thread_block();
				U_r *sdata = gpu_detail::SharedMemory<U_r>();

				dt_uint32 tid = threadIdx.x + threadIdx.y*blockDim.x;
				dt_uint32 bid = blockIdx.x + blockIdx.y*gridDim.x;
				dt_uint32 blocksize = blockDim.x*blockDim.y;

				KS<U_r> sum = U_r(0);

				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_int_det_ring_norm_2(ix, iy, grid, g2_min, g2_max, mx_i.m_data, sum));

				GPU_BLK_REDUCE(sum, tid, sdata, bid, mx_o);
			}

			/* integration over a detector with sensitivity */
			template <class T, class U>
			__global__ void fcn_int_det_sen(Grid_2d<T> grid, pVctr_gpu_32<U> sen_i, pVctr_gpu_32<U> mx_i, pVctr_gpu_32<U> mx_o)
			{ 
				gpu_cg::thread_block cta = gpu_cg::this_thread_block();
				U *sdata = gpu_detail::SharedMemory<U>();

				dt_uint32 tid = threadIdx.x + threadIdx.y*blockDim.x;
				dt_uint32 bid = blockIdx.x + blockIdx.y*gridDim.x;
				dt_uint32 blocksize = blockDim.x*blockDim.y;

				KS<U> sum = U(0);

				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_int_det_sen(ix, iy, grid, sen_i.m_data, mx_i.m_data, sum));

				GPU_BLK_REDUCE(sum, tid, sdata, bid, mx_o);
			}

			/* norm_2 integration over a detector with sensitivity */
			template <class T, class U>
			__global__ void fcn_int_det_sen_norm_2(Grid_2d<T> grid, pVctr_gpu_32<Value_type_r<U>> sen_i, pVctr_gpu_32<U> mx_i, pVctr_gpu_32<Value_type_r<U>> mx_o)
			{ 
				using U_r = Value_type_r<U>;
				gpu_cg::thread_block cta = gpu_cg::this_thread_block();
				U_r *sdata = gpu_detail::SharedMemory<U_r>();

				dt_uint32 tid = threadIdx.x + threadIdx.y*blockDim.x;
				dt_uint32 bid = blockIdx.x + blockIdx.y*gridDim.x;
				dt_uint32 blocksize = blockDim.x*blockDim.y;

				KS<U_r> sum = U_r(0);

				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_int_det_sen_norm_2(ix, iy, grid, sen_i.m_data, mx_i.m_data, sum));

				GPU_BLK_REDUCE(sum, tid, sdata, bid, mx_o);
			}
		}

		/* wave propagation */
		namespace gpu_detail_mt
		{
			/* propagate */
			template <class T, class U>
			__global__ void fcn_fs_propagate(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, R_2d<T> g_0, T w_g2, T w, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate(ix, iy, grid, psi_i.m_data, g_0, w_g2, w, psi_o.m_data));
			}

			/* propagate and bandwith limit using a fermi aperture */
			template <class T, class U>
			__global__ void fcn_fs_propagate_bw_f(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T alpha, T w, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate_bw_f(ix, iy, grid, psi_i.m_data, g_0, w_g2, g2_cut, alpha, w, psi_o.m_data));
			}

			/* propagate and bandwith limit using a hard aperture */
			template <class T, class U>
			__global__ void fcn_fs_propagate_bw_h(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T w, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate_bw_h(ix, iy, grid, psi_i.m_data, g_0, w_g2, g2_cut, w, psi_o.m_data));
			}
		}

		/* probe - ctf - pctf */
		namespace gpu_detail_mt
		{
			/* Convergent incident wave in Fourier space */
			template <dt_bool bb_phi, class T, class U>
			__global__ void fcn_fs_probe(Grid_2d<T> grid, Lens<T> lens, R_2d<T> r, R_2d<T> gu, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_probe<bb_phi>(ix, iy, grid, lens, r, gu, psi_o.m_data));
			}

			/* apply coherent transfer function */
			template <dt_bool bb_phi, class T, class U>
			__global__ void fcn_fs_apply_ctf(Grid_2d<T> grid, pVctr_gpu_32<U> psi_i, Lens<T> lens, R_2d<T> gu, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_apply_ctf<bb_phi>(ix, iy, grid, psi_i.m_data, lens, gu, psi_o.m_data));
			}

			/* apply partial coherent transfer function */
			// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 611
			// tp_inc_iehwgd: 08, spt_inc_theta_c: 611
			template <dt_bool bb_phi, class T, class U>
			__global__ void fcn_fs_apply_pctf(iGrid_2d grid, pVctr_gpu_32<U> psi_i, Lens<T> lens, pVctr_gpu_32<U> psi_o)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_apply_pctf<bb_phi>(ix, iy, grid, psi_i.m_data, lens, psi_o.m_data));
			}
		}

		/* transmission function */
		namespace gpu_detail_mt
		{
			template <eElec_Spec_Int_Mod esim, class T, class U>
			__global__ void fcn_trans_fcn(iGrid_1d igrid, pVctr_gpu_32<T> vzp_i, const T& w, pVctr_gpu_32<U> tfcn_o)
			{
				FOR_LOOP_1DC(igrid.nx, cgpu_detail_mt::fcn_trans_fcn<esim>(ix, igrid, vzp_i.m_data, w, tfcn_o.m_data));
			}
		}

		/* eels */
		namespace gpu_detail_mt
		{
			template <class T>
			__global__ void fcn_eels_lorentz_norm_factor(Grid_2d<T> grid, T gc2, T ge2, pVctr_gpu_32<T> mx_o)
			{ 
				gpu_cg::thread_block cta = gpu_cg::this_thread_block();
				T *sdata = gpu_detail::SharedMemory<T>();

				dt_uint32 tid = threadIdx.x + threadIdx.y*blockDim.x;
				dt_uint32 bid = blockIdx.x + blockIdx.y*gridDim.x;
				dt_uint32 blocksize = blockDim.x*blockDim.y;

				KS<T> sum = T(0);

				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_lorentz_norm_factor(ix, iy, grid, gc2, ge2, sum));

				GPU_BLK_REDUCE(sum, tid, sdata, bid, mx_o);
			}

			template <class T, class U>
			__global__ void fcn_eels_w_xyz(Grid_2d<T> grid, EELS<T> eels, 
				pVctr_gpu_32<U> w_x, pVctr_gpu_32<U> w_y, pVctr_gpu_32<U> w_z)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_xyz(ix, iy, grid, eels, w_x.m_data, w_y.m_data, w_z.m_data));
			}

			template <class T, class U>
			__global__ void fcn_eels_w_x(Grid_2d<T> grid, EELS<T> eels, pVctr_gpu_32<U> w_x)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_x(ix, iy, grid, eels, w_x.m_data));
			}

			template <class T, class U>
			__global__ void fcn_eels_w_y(Grid_2d<T> grid, EELS<T> eels, pVctr_gpu_32<U> w_y)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_y(ix, iy, grid, eels, w_y.m_data));
			}

			template <class T, class U>
			__global__ void fcn_eels_w_z(Grid_2d<T> grid, EELS<T> eels, pVctr_gpu_32<U> w_z)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_z(ix, iy, grid, eels, w_z.m_data));
			}

			template <class T, class U>
			__global__ void fcn_eels_w_mn1(Grid_2d<T> grid, EELS<T> eels, pVctr_gpu_32<U> w_mn1)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_mn1(ix, iy, grid, eels, w_mn1.m_data));
			}

			template <class T, class U>
			__global__ void fcn_eels_w_mp1(Grid_2d<T> grid, EELS<T> eels, pVctr_gpu_32<U> w_mp1)
			{
				FOR_LOOP_2DC(grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_mp1(ix, iy, grid, eels, w_mp1.m_data));
			}
		}
	}

#endif
