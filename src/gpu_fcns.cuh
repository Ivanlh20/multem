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

// #if defined __CUDACC__ && !defined GPU_FCNS_H
#ifndef GPU_FCNS_H
	#define GPU_FCNS_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "types.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	#include "cgpu_fcns.cuh"
	#include "grid.cuh"
	#include "gpu_detail.cuh"

	/* vector functions */
	namespace mt
	{
		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_assign_real(TVctr_c& mx_i, TVctr_r& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::assign_real<T>());
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_assign_max_real(TVctr_c& mx_i, Value_type<TVctr_r> M_v, TVctr_r& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::assign_max_real<T>(M_v, M_v));
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_and_rvctr<TVctr_c, TVctr_r, void>
		fcn_assign_real(TVctr_c& mx_i, Value_type<TVctr_r> M_v, TVctr_r& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			if (M_v>0)
			{
				fcn_assign_max_real(mx_i, T(0), mx_o, pstream);
			}
			else
			{
				fcn_assign_real(mx_i, mx_o, pstream);
			}
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_assign_abs_real(TVctr_c& mx_i, TVctr_r& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::assign_abs_real<T>());
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_c, TVctr_r, void>
		fcn_assign_abs(TVctr_c& mx_i, TVctr_r& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::assign_abs<T>());
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fill(TVctr& mx_io, Value_type<TVctr> val, Stream_gpu* pstream = nullptr)
		{
			thrust::fill(mx_io.begin(), mx_io.end(), val);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_scale(const Value_type<TVctr_2>& w_i, TVctr_1& mx_i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::scale<T>(w_i));
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_scale(const Value_type<TVctr>& w_i, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			thrust::transform(mx_io.begin(), mx_io.end(), mx_io.begin(), cgpu_fctr::scale<T>(w_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_norm_2(TVctr_1& mx_i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::norm_2<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_scale_norm_2(Value_type<TVctr_2> w_i, TVctr_1& mx_i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_o.begin(), cgpu_fctr::scale_norm_2<T>(w_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_sub(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::sub<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::add<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add(TVctr_1& mx_i, TVctr_2& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_io.begin(), mx_io.begin(), cgpu_fctr::add<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_scale(Value_type<TVctr_1> w1_i, TVctr_1& mx_1i, 
		Value_type<TVctr_1> w2_i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::add_scale_i<T>(w1_i, w2_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_scale(Value_type<TVctr_1> w_i, TVctr_1& mx_i, TVctr_2& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_io.begin(), mx_io.begin(), cgpu_fctr::add_scale<T>(w_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_norm_2(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::add_norm_2_i<T>());

		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_norm_2(TVctr_1& mx_i, TVctr_2& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_io.begin(), mx_io.begin(), cgpu_fctr::add_norm_2<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_scale_norm_2(Value_type<TVctr_2> w1_i, TVctr_1& mx_1i, 
		Value_type<TVctr_2> w2_i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::add_scale_norm_2_i<T>(w1_i, w2_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_scale_norm_2(Value_type<TVctr_2> w_i, TVctr_1& mx_i, TVctr_2& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_io.begin(), mx_io.begin(), cgpu_fctr::add_scale_norm_2<T>(w_i));
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_mult(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_1i.begin(), mx_1i.end(), mx_2i.begin(), mx_o.begin(), cgpu_fctr::mult<T>());
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_mult(TVctr_1& mx_i, TVctr_2& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			thrust::transform(mx_i.begin(), mx_i.end(), mx_io.begin(), mx_io.begin(), cgpu_fctr::mult<T>());
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_div_sft(Value_type<TVctr> mx_sft, Value_type<TVctr> mx_div, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			thrust::transform(mx_io.begin(), mx_io.end(), mx_io.begin(), cgpu_fctr::div_sft<T>(mx_sft, mx_div));
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_sum(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			return thrust::reduce(mx_i.begin(), mx_i.end());
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_norm_2(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::norm_2<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_norm_2_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::norm_2_sft<T, T_r>(x_sft), T_r(0), cgpu_fctr::add<T_r>());
		}

 		template <class TVctr>
		enable_if_cvctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_max_real(TVctr& mx_i, Value_type_r<TVctr> v_min, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::assign_max_real<T_r>(v_min, T_r(0)), T_r(0), cgpu_fctr::add<T_r>());
		}

 		template <class TVctr>
		enable_if_cvctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs_real(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::assign_abs_real<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
		}

 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::assign_abs<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
		}

 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;
			return thrust::transform_reduce(mx_i.begin(), mx_i.end(), 
			cgpu_fctr::abs_sft<T, T_r>(), T_r(0), cgpu_fctr::add<T_r>());
		}

		/***************************************************************************************/
		/* calculate mean */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_mean(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum(mx_i, pstream)/T(mx_i.size());
		}

		/* calculate mean of the norm square */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_norm_2(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_norm_2(mx_i, pstream)/T_r(mx_i.size());
		}

		/* calculate mean of the shifted norm square */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_norm_2_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_norm_2_sft(mx_i, x_sft, pstream)/T_r(mx_i.size());
		}

 		template <class TVctr>
		enable_if_cvctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_max_real(TVctr& mx_i, Value_type_r<TVctr> v_min, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum_max_real(mx_i, v_min, pstream)/T(mx_i.size());
		}

 		template <class TVctr>
		enable_if_cvctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs_real(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum_abs_real(mx_i, pstream)/T(mx_i.size());
		}

 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_abs(mx_i, pstream)/T_r(mx_i.size());
		}

 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_gpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_abs_sft(mx_i, x_sft, pstream)/T_r(mx_i.size());
		}

		/***************************************************************************************/
		/* calculate mean and variance */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_mean_var(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_var, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			x_mean = fcn_mean(mx_i, pstream);
			x_var = fcn_mean_norm_2_sft(mx_i, x_mean, pstream);
		}

		/* calculate mean and standard deviation */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_mean_std(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_std, Stream_gpu* pstream = nullptr)
		{
			fcn_mean_var(mx_i, x_mean, x_std, pstream);
			x_std = sqrt(x_std);
		}

		// mean and first absolute moment
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_mean_moment_1a(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_ma1, Stream_gpu* pstream = nullptr)
		{
			x_mean = fcn_mean(mx_i, pstream);
			x_ma1 = fcn_mean_abs_sft(mx_i, x_mean, pstream);
		}

		/* variance */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_variance(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			T x_mean;
			T_r x_var;
			fcn_mean_var(mx_i, x_mean, x_var, pstream);

			return x_var;
		}

		/* standard deviation */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_std(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			auto x_var = fcn_variance(mx_i, pstream);

			return sqrt(x_var);
		}

		// first absolute moment
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_moment_1a(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			T x_mean;
			T_r x_ma1;
			fcn_mean_moment_1a(mx_i, x_mean, x_ma1, pstream);

			return x_ma1;
		}

		/***************************************************************************************/
		/* minimum element of a vector */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_min_element(const TVctr& x, Stream_gpu* pstream = nullptr)
		{
			return *thrust::min_element(x.begin(), x.end());
		}

		/* maximum element of a vector */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_max_element(const TVctr& x, Stream_gpu* pstream = nullptr)
		{
			return *thrust::max_element(x.begin(), x.end());
		}

		/* minimum and maximum element of a vector */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_minmax_element(const TVctr& x, Value_type<TVctr>& x_min, Value_type<TVctr>& x_max, Stream_gpu* pstream = nullptr)
		{
			auto x_min_max = thrust::minmax_element(x.begin(), x.end());
			x_min = *(x_min_max.first);
			x_max = *(x_min_max.second);
		}
	}

	/* add - assign - crop - norm_2 - fftsft */
	namespace mt
	{
		/* shift matrix respect to nx_h */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_1d(TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_1d();

			gpu_detail::fcn_fftsft_1d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_io.ptr_32());
		}

		/* shift matrix respect to ny_h */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_bc_2d(TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			gpu_detail::fcn_fftsft_bc_2d<T><<<igrid.d_grid_d0_h(), igrid.d_blk()>>>(igrid, mx_io.ptr_32());
		}

		/* shift matrix respect to (nx_h, ny_h) */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_2d(TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			gpu_detail::fcn_fftsft_2d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_io.ptr_32());
		}

		/* shift matrix respect to (nx_h, ny_h) */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_2d(TVctr& mx_i, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_fftsft_2d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), mx_o.ptr_32());
		}

		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_3d( TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_3d();

			gpu_detail::fcn_fftsft_3d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_io.ptr_32());
		}

		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fftsft_3d(TVctr& mx_i, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_3d();

			gpu_detail::fcn_fftsft_3d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), mx_o.ptr_32());
		}

		/***************************************************************************************/
		/* add, scale and shift */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_add_sc_fftsft_2d(TVctr& mx_i, Value_type<TVctr> w, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_add_sc_fftsft_2d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), w, mx_o.ptr_32());
		} 

		/* add, scale, square and shift */
 		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_sc_norm_2_fftsft_2d(TVctr_1& mx_i, Value_type<TVctr_2> w, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_1>;
			using U = Value_type<TVctr_2>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_add_sc_norm_2_fftsft_2d<T, U><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), w, mx_o.ptr_32());
		}

		/***************************************************************************************/
		/* assign and crop */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_assign_crop_2d(TVctr& mx_i, TVctr& mx_o, iRegion_Rect_2d& iregion, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_assign_crop_2d<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), iregion, mx_o.ptr_32());
		}

		/* assign, crop and shift */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_assign_crop_fftsft_2d(TVctr& mx_i, iRegion_Rect_2d& iregion, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_assign_crop_fftsft_2d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), iregion, mx_o.ptr_32());
		}

 		/* add, scale, crop and shift */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_add_sc_crop_fftsft_2d(TVctr& mx_i, iRegion_Rect_2d& iregion, Value_type<TVctr> w, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_add_sc_crop_fftsft_2d<T><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), iregion, w, mx_o.ptr_32());
		}

		/* add, scale, square, crop and shift */
 		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		fcn_add_sc_norm_2_crop_fftsft_2d(TVctr_1& mx_i, iRegion_Rect_2d& iregion, Value_type<TVctr_2> w, TVctr_2& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_1>;
			using U = Value_type<TVctr_2>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_add_sc_norm_2_crop_fftsft_2d<T, U><<<igrid.d_grid_h(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), iregion, w, mx_o.ptr_32());
		}
	}

	/* transpose - element wise matrix op vector */
	namespace mt
	{
		/* transpose */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, Vctr_gpu<Value_type<TVctr>>>
		fcn_trs_2d(TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			Vctr_gpu<T> mx_t(mx_i.shape_2d_trs());
			gpu_detail::fcn_trs_2d<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), mx_t.ptr_32());

			return mx_t;
		}

		/* transpose */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_trs_2d(TVctr& mx_i, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			if (&mx_i != &mx_o)
			{
				gpu_detail::fcn_trs_2d<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), mx_o.ptr_32());
			}
			else
			{
				Vctr_gpu<T> mx_t(mx_i.shape_2d_trs());
				gpu_detail::fcn_trs_2d<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), mx_t.ptr_32());
				mx_t.cpy_to_gpu_ptr(mx_o.m_data, mx_o.size());
			}
		}

		/* element wise addition: matrix + vector*/
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_ew_add_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				gpu_detail::fcn_ew_add_mx_vctr_col<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				gpu_detail::fcn_ew_add_mx_vctr_row<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
		}

		/* element wise subtraction: matrix - vector */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_ew_sub_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				gpu_detail::fcn_ew_sub_mx_vctr_col<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				gpu_detail::fcn_ew_sub_mx_vctr_row<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
		}

		/* element wise multiplication matrix X vector */
 		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_ew_mult_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				gpu_detail::fcn_ew_mult_mx_vctr_col<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				gpu_detail::fcn_ew_mult_mx_vctr_row<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vctr, mx_io.ptr_32());
			}
		}
	}

	/* aperture functions */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fermi_aperture(Grid_2d<T>& grid, T g2_cut, T alpha, T w, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_fermi_aperture<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, g2_cut, alpha, w, mx_io.ptr_32());
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_hard_aperture(Grid_2d<T>& grid, T g2_cut, T w, TVctr& mx_io, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_hard_aperture<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, g2_cut, w, mx_io.ptr_32());
		}
	}

	/* phase shifts real space*/
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_rs_exp_factor_1d(Grid_1d<T>& grid, TVctr& psi_i, T gx, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_rs_exp_factor_1d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), gx, w, psi_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_rs_exp_factor_2d_bc(Grid_2d<T>& grid, TVctr_c& psi_i, T alpha, TVctr_r &gy, T w, TVctr_c& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_rs_exp_factor_2d_bc<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, alpha, psi_i.ptr_32(), gy.ptr_32(), w, psi_o.ptr_32());
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_rs_exp_factor_2d(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_rs_exp_factor_2d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), g, w, psi_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_rs_mul_exp_factor_2d(Grid_2d<T>& grid, TVctr_c& psi_i, TVctr_r& g, T w, TVctr_c& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_rs_mul_exp_factor_2d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), g.ptr_32(), w, psi_o.ptr_32());
		}
	}

	/* phase shifts fourier space */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_exp_factor_1d(Grid_1d<T>& grid, TVctr& psi_i, T rx, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_fs_exp_factor_1d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), rx, w, psi_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_fs_exp_factor_2d_bc(Grid_2d<T>& grid, TVctr_c& psi_i, T alpha, TVctr_r &ry, T w, TVctr_c& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_fs_exp_factor_2d_bc<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, alpha, psi_i.ptr_32(), ry.ptr_32(), w, psi_o.ptr_32());
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_exp_factor_2d(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> r, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_fs_exp_factor_2d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), r, w, psi_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_fs_mul_exp_factor_2d(Grid_2d<T>& grid, TVctr_c& psi_i, TVctr_r& r, T w, TVctr_c& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_fs_mul_exp_factor_2d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), r.ptr_32(), w, psi_o.ptr_32());
		}
	}

	/* gradient */
	namespace mt
	{
		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_grad_x(TVctr& mx_i, TVctr& dmx_x, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_grad_x<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), dmx_x.ptr_32());
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_grad_y(TVctr& mx_i, TVctr& dmx_y, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_grad_y<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), dmx_y.ptr_32());
		}

		template <class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_grad(TVctr& mx_i, TVctr& dmx_x, TVctr& dmx_y, Stream_gpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			gpu_detail::fcn_grad<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_i.ptr_32(), dmx_x.ptr_32(), dmx_y.ptr_32());
		}
	}

	/* function multiplication fourier space */
	namespace mt
	{
		#define FCN_MULT_FS_FCN_GPU(FN, FCN, DIM)																							\
		template <class T, class TVctr_c>																									\
		enable_if_vctr_gpu<TVctr_c, void>																									\
		fcn_mult_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FCN<T>& fcn, T w, TVctr_c& mx_io, Stream_gpu* pstream = nullptr)				\
		{																																	\
			using U = Value_type<TVctr_c>;																									\
																																			\
			gpu_detail::fcn_mult_fs_##FN##_##DIM##d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, fcn, w, mx_io.ptr_32());					\
		}

		/* gaussian */
		FCN_MULT_FS_FCN_GPU(gauss, Fcn_Gauss, 1);		// fcn_mult_fs_gauss_1d
		FCN_MULT_FS_FCN_GPU(gauss, Fcn_Gauss, 2);		// fcn_mult_fs_gauss_2d
		FCN_MULT_FS_FCN_GPU(gauss, Fcn_Gauss, 3);		// fcn_mult_fs_gauss_3d

		/* exponential */
		FCN_MULT_FS_FCN_GPU(exp, Fcn_Exp, 1);			// fcn_mult_fs_exp_1d
		FCN_MULT_FS_FCN_GPU(exp, Fcn_Exp, 2);			// fcn_mult_fs_exp_2d
		FCN_MULT_FS_FCN_GPU(exp, Fcn_Exp, 3);			// fcn_mult_fs_exp_3d

		/* fermi */
		FCN_MULT_FS_FCN_GPU(fermi, Fcn_Fermi, 1);		// fcn_mult_fs_fermi_1d
		FCN_MULT_FS_FCN_GPU(fermi, Fcn_Fermi, 2);		// fcn_mult_fs_fermi_2d
		FCN_MULT_FS_FCN_GPU(fermi, Fcn_Fermi, 3);		// fcn_mult_fs_fermi_3d

		/* butterworth */
		FCN_MULT_FS_FCN_GPU(butwth, Fcn_Butwth, 1);		// fcn_mult_fs_butwth_1d
		FCN_MULT_FS_FCN_GPU(butwth, Fcn_Butwth, 2);		// fcn_mult_fs_butwth_2d
		FCN_MULT_FS_FCN_GPU(butwth, Fcn_Butwth, 3);		// fcn_mult_fs_butwth_3d

		/* hann */
		FCN_MULT_FS_FCN_GPU(hann, Fcn_Hann, 1);			// fcn_mult_fs_hann_1d
		FCN_MULT_FS_FCN_GPU(hann, Fcn_Hann, 2);			// fcn_mult_fs_hann_2d
		FCN_MULT_FS_FCN_GPU(hann, Fcn_Hann, 3);			// fcn_mult_fs_hann_3d
	}
	
	/* convolution */
	namespace mt
	{
		#define FCN_CV_FS_FCN_GPU(FN, FCN, DIM)																									\
		template <class T, class TVctr_c>																										\
		enable_if_cvctr_gpu<TVctr_c, void>																										\
		fcn_cv_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FFT_gpu<T>& fft, FCN<T>& fcn, TVctr_c& mx_io, Stream_gpu* pstream = nullptr)			\
		{																																		\
			using U = Value_type<TVctr_c>;																										\
																																				\
			fcn_fftsft_##DIM##d(mx_io, pstream);																								\
			fft.forward(mx_io);																													\
																																				\
			const T w = grid.isize_r();																											\
																																				\
			fcn_mult_fs_##FN##_##DIM##d(grid, fcn, w, mx_io, pstream);																			\
																																				\
			fft.inverse(mx_io);																													\
			fcn_fftsft_##DIM##d(mx_io, pstream);																								\
		}

		/* gaussian convolution */
		FCN_CV_FS_FCN_GPU(gauss, Fcn_Gauss, 1);			// fcn_cv_fs_gauss_1d
		FCN_CV_FS_FCN_GPU(gauss, Fcn_Gauss, 2);			// fcn_cv_fs_gauss_2d
		FCN_CV_FS_FCN_GPU(gauss, Fcn_Gauss, 3);			// fcn_cv_fs_gauss_3d

		/* exponential convolution */
		FCN_CV_FS_FCN_GPU(exp, Fcn_Exp, 1);				// fcn_cv_fs_exp_1d
		FCN_CV_FS_FCN_GPU(exp, Fcn_Exp, 2);				// fcn_cv_fs_exp_2d
		FCN_CV_FS_FCN_GPU(exp, Fcn_Exp, 3);				// fcn_cv_fs_exp_3d

		/* fermi convolution */
		FCN_CV_FS_FCN_GPU(fermi, Fcn_Fermi, 1);			// fcn_cv_fs_fermi_1d
		FCN_CV_FS_FCN_GPU(fermi, Fcn_Fermi, 2);			// fcn_cv_fs_fermi_2d
		FCN_CV_FS_FCN_GPU(fermi, Fcn_Fermi, 3);			// fcn_cv_fs_fermi_3d

		/* butterworth convolution */
		FCN_CV_FS_FCN_GPU(butwth, Fcn_Butwth, 1);		// fcn_cv_fs_butwth_1d
		FCN_CV_FS_FCN_GPU(butwth, Fcn_Butwth, 2);		// fcn_cv_fs_butwth_2d
		FCN_CV_FS_FCN_GPU(butwth, Fcn_Butwth, 3);		// fcn_cv_fs_butwth_3d

		/* hann convolution */
		FCN_CV_FS_FCN_GPU(hann, Fcn_Hann, 1);			// fcn_cv_fs_hann_1d
		FCN_CV_FS_FCN_GPU(hann, Fcn_Hann, 2);			// fcn_cv_fs_hann_2d
		FCN_CV_FS_FCN_GPU(hann, Fcn_Hann, 3);			// fcn_cv_fs_hann_3d
	}

	/* deconvolution */
	namespace mt
	{
		#define FCN_DCV_FS_FCN_GPU(FN, FCN, DIM)																									\
		template <class T, class TVctr_c>																											\
		enable_if_cvctr_gpu<TVctr_c, void>																											\
		fcn_dcv_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FFT_gpu<T>& fft, FCN<T>& fcn, T psnr, TVctr_c& mx_io, Stream_gpu* pstream = nullptr)		\
		{																																			\
			using U = Value_type<TVctr_c>;																											\
																																					\
			fcn_fftsft_##DIM##d(mx_io, pstream);																									\
			fft.forward(mx_io);																														\
																																					\
			const T w = grid.isize_r();																												\
																																					\
			gpu_detail::fcn_dcv_fs_##FN##_##DIM##d<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, fcn, psnr, w, mx_io.ptr_32());						\
																																					\
			fft.inverse(mx_io);																														\
			fcn_fftsft_##DIM##d(mx_io, pstream);																									\
		}

		/* gaussian convolution */
		FCN_DCV_FS_FCN_GPU(gauss, Fcn_Gauss, 1);			// fcn_dcv_fs_gauss_1d
		FCN_DCV_FS_FCN_GPU(gauss, Fcn_Gauss, 2);			// fcn_dcv_fs_gauss_2d
		FCN_DCV_FS_FCN_GPU(gauss, Fcn_Gauss, 3);			// fcn_dcv_fs_gauss_3d

		/* exponential convolution */
		FCN_DCV_FS_FCN_GPU(exp, Fcn_Exp, 1);				// fcn_dcv_fs_exp_1d
		FCN_DCV_FS_FCN_GPU(exp, Fcn_Exp, 2);				// fcn_dcv_fs_exp_2d
		FCN_DCV_FS_FCN_GPU(exp, Fcn_Exp, 3);				// fcn_dcv_fs_exp_3d

		/* fermi convolution */
		FCN_DCV_FS_FCN_GPU(fermi, Fcn_Fermi, 1);			// fcn_dcv_fs_fermi_1d
		FCN_DCV_FS_FCN_GPU(fermi, Fcn_Fermi, 2);			// fcn_dcv_fs_fermi_2d
		FCN_DCV_FS_FCN_GPU(fermi, Fcn_Fermi, 3);			// fcn_dcv_fs_fermi_3d

		/* butterworth convolution */
		FCN_DCV_FS_FCN_GPU(butwth, Fcn_Butwth, 1);			// fcn_dcv_fs_butwth_1d
		FCN_DCV_FS_FCN_GPU(butwth, Fcn_Butwth, 2);			// fcn_dcv_fs_butwth_2d
		FCN_DCV_FS_FCN_GPU(butwth, Fcn_Butwth, 3);			// fcn_dcv_fs_butwth_3d

		/* hann convolution */
		FCN_DCV_FS_FCN_GPU(hann, Fcn_Hann, 1);				// fcn_dcv_fs_hann_1d
		FCN_DCV_FS_FCN_GPU(hann, Fcn_Hann, 2);				// fcn_dcv_fs_hann_2d
		FCN_DCV_FS_FCN_GPU(hann, Fcn_Hann, 3);				// fcn_dcv_fs_hann_3d
	}
	
	/* window functions */
	namespace mt
	{
		#define FCN_WD_FCN_GPU(FN, FCN, DIM)																						\
		template <class T, class TVctr>																								\
		enable_if_vctr_gpu<TVctr, void>																								\
		fcn_wd_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FCN<T>& fcn, dt_bool sft, TVctr& mx_o, Stream_gpu* pstream = nullptr)		\
		{																															\
			using U = Value_type<TVctr>;																							\
																																	\
			if (sft)																												\
			{																														\
				gpu_detail::fcn_wd_##FN##_##DIM##d<true, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, fcn, mx_o.ptr_32());			\
			}																														\
			else																													\
			{																														\
				gpu_detail::fcn_wd_##FN##_##DIM##d<false, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, fcn, mx_o.ptr_32());			\
			}																														\
		}

		 FCN_WD_FCN_GPU(gauss, Wd_Gauss_1d, 1);		// fcn_wd_gauss_1d
		 FCN_WD_FCN_GPU(gauss, Wd_Gauss_2d, 2);		// fcn_wd_gauss_2d
		 FCN_WD_FCN_GPU(gauss, Wd_Gauss_3d, 3);		// fcn_wd_gauss_3d

		 FCN_WD_FCN_GPU(exp, Wd_Exp_1d, 1);			// fcn_wd_exp_1d
		 FCN_WD_FCN_GPU(exp, Wd_Exp_2d, 2);			// fcn_wd_exp_2d
		 FCN_WD_FCN_GPU(exp, Wd_Exp_3d, 3);			// fcn_wd_exp_3d

		 FCN_WD_FCN_GPU(fermi, Wd_Fermi_1d, 1);		// fcn_wd_fermi_1d
		 FCN_WD_FCN_GPU(fermi, Wd_Fermi_2d, 2);		// fcn_wd_fermi_2d
		 FCN_WD_FCN_GPU(fermi, Wd_Fermi_3d, 3);		// fcn_wd_fermi_3d

		 FCN_WD_FCN_GPU(butwth, Wd_Butwth_1d, 1);	// fcn_wd_butwth_1d
		 FCN_WD_FCN_GPU(butwth, Wd_Butwth_2d, 2);	// fcn_wd_butwth_2d
		 FCN_WD_FCN_GPU(butwth, Wd_Butwth_3d, 3);	// fcn_wd_butwth_3d

		 FCN_WD_FCN_GPU(hann, Wd_Hann_1d, 1);		// fcn_wd_hann_1d
		 FCN_WD_FCN_GPU(hann, Wd_Hann_2d, 2);		// fcn_wd_hann_2d
		 FCN_WD_FCN_GPU(hann, Wd_Hann_3d, 3);		// fcn_wd_hann_3d
	}

	/* phase correlation */
	namespace mt
	{
		/****************** pcf data processing real space *******************/
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_1d_dp(Grid_1d<T>& grid, TVctr_r& mx_i, Wd_Butwth_1d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_rs_pcf_1d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_2d_dp(Grid_2d<T>& grid, TVctr_r& mx_i, Wd_Butwth_2d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_rs_pcf_2d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}
		
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_3d_dp(Grid_3d<T>& grid, TVctr_r& mx_i, Wd_Butwth_3d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_rs_pcf_3d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}

		/***************** pcf data processing fourier space *****************/
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_1d_dp(Grid_1d<T>& grid, TVctr_r& mx_i, Wd_Gauss_1d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_fs_pcf_1d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_2d_dp(Grid_2d<T>& grid, TVctr_r& mx_i, Wd_Gauss_2d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_fs_pcf_2d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}
		
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_3d_dp(Grid_3d<T>& grid, TVctr_r& mx_i, Wd_Gauss_3d<T>& fcn, T w, TVctr_c& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			gpu_detail::fcn_fs_pcf_3d_dp<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, mx_i.ptr_32(), fcn, w, mx_o.ptr_32());
		}
	}

	/* optical flow */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_opt_flow(TVctr& mx_s, TVctr& mx_m, T alpha, TVctr& dphi_x, TVctr& dphi_y, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto igrid = mx_s.igrid_2d();

			gpu_detail::fcn_opt_flow<T><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, mx_s.ptr_32(), mx_m.ptr_32(), alpha, dphi_x.ptr_32(), dphi_y.ptr_32());
		}
	}

	/* bilinear interpolation */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_intrpl_bl_rg_2d(Grid_2d<T>& grid, TVctr& mx_i, TVctr& vx, TVctr& vy, T bg, TVctr& mx_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail::fcn_intrpl_bl_rg_2d<T><<<grid.d_grid_size(), grid.d_blk_size()>>>(grid, mx_i.ptr_32(), vx.ptr_32(), vy.ptr_32(), bg, mx_o.ptr_32());
		}
	}

	namespace mt
	{
	// 	/***************************************************************************************/
	// 	template <class TGrid, class TVctr_r>
	// 	enable_if_vctr_gpu<TVctr_r, Value_type<TVctr_r>>
	// 	atom_cost_function(Grid_2d<T>& grid, const Atom_Sa<Value_type<TGrid>>& atom_Ip, TVctr_r& mx_i)
	// 	{
	// 		TVctr_r sum_v(1);

	// 		D_Grid_Blk d_grid_blk;
	// 		d_grid_blk.grid = dim3();
	// 		d_grid_blk.blk = dim3(c_thr_2d_x, c_thr_2d_y);

	// 		gpu_detail::atom_cost_function<typename TVctr_r::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, atom_Ip, mx_i.ptr_32(), sum_v);
	// 		return sum_v[0];
	// 	}

	// 	template <class TGrid, class TVctr_r>
	// 	enable_if_vctr_gpu<TVctr_r, void>
	// 	subtract_atom(Grid_2d<T>& grid, Vctr<Atom_Sa<Value_type<TGrid>>, edev_cpu>& atom_Ip, TVctr_r& mx_i)
	// 	{
	// 		if (stream.n_stream_act<= 0)
	// 		{
	// 			return;
	// 		}

	// 		for(auto istm = 0; istm < stream.n_stream_act; istm++)
	// 		{
	// 			D_Grid_Blk d_grid_blk;
	// 			d_grid_blk.grid = dim3();
	// 			d_grid_blk.blk = dim3(c_thr_2d_x, c_thr_2d_y);

	// 			gpu_detail::subtract_atom<typename TVctr_r::value_type><<<d_grid_blk.grid, d_grid_blk.blk, 0, stream[istm]>>>(grid_2d, atom_Ip[istm], mx_i);
	// 		}
	// 	}

	// 	// Linear projected potential: V and zV
	// 	template <class TQ1, class TVAtom>
	// 	enable_if_device<TQ1, void>
	// 	linear_Vz(ePot_Parm_Typ pot_parm_typ, TQ1 &qz, TVAtom &vatom)
	// 	{	
	// 		using TAtom = Value_type<TVAtom>;

	// 		if (stream.n_stream_act<= 0)
	// 		{
	// 			return;
	// 		}

	// 		auto str_linear_Vz = [](cudaStream_t &stream, const ePot_Parm_Typ& pot_parm_typ, TQ1 &qz, TAtom &atom)
	// 		{
	// 			if (atom.charge == 0)
	// 			{
	// 				switch(pot_parm_typ)
	// 				{
	// 					case ePPT_doyle_0_4:
	// 						gpu_detail::linear_Vz<ePPT_doyle_0_4, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_peng_0_4:
	// 						gpu_detail::linear_Vz<ePPT_peng_0_4, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_peng_0_12:
	// 						gpu_detail::linear_Vz<ePPT_peng_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_kirkland_0_12:
	// 						gpu_detail::linear_Vz<ePPT_kirkland_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_weickenmeier_0_12:
	// 						gpu_detail::linear_Vz<ePPT_weickenmeier_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_lobato_0_12:
	// 						gpu_detail::linear_Vz<ePPT_lobato_0_12, 0, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 				}
	// 			}
	// 			else
	// 			{
	// 				switch(pot_parm_typ)
	// 				{
	// 					case ePPT_doyle_0_4:
	// 						gpu_detail::linear_Vz<ePPT_doyle_0_4, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_peng_0_4:
	// 						gpu_detail::linear_Vz<ePPT_peng_0_4, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_peng_0_12:
	// 						gpu_detail::linear_Vz<ePPT_peng_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_kirkland_0_12:
	// 						gpu_detail::linear_Vz<ePPT_kirkland_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_weickenmeier_0_12:
	// 						gpu_detail::linear_Vz<ePPT_weickenmeier_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 					case ePPT_lobato_0_12:
	// 						gpu_detail::linear_Vz<ePPT_lobato_0_12, 1, TAtom><<<dim3(c_nR), dim3(c_nqz), 0, stream>>>(qz, atom);
	// 						break;
	// 				}		
	// 			}
	// 		};

	// 		for(auto istm = 0; istm < stream.n_stream_act; istm++)
	// 		{
	// 			str_linear_Vz(stream[istm], pot_parm_typ, qz, vatom[istm]);
	// 		}
	// 	}

	// 	// Get local interpolation coefficients
	// 	template <class TVAtom> 
	// 	enable_if_device<typename TVAtom::value_type, void>
	// 	fcn_vd_2_coef_poly3(TVAtom &vatom)
	// 	{
	// 		using TAtom = Value_type<TVAtom>;

	// 		if (stream.n_stream_act<= 0)
	// 		{
	// 			return;
	// 		}

	// 		for(auto istm = 0; istm < stream.n_stream_act; istm++)
	// 		{
	// 			gpu_detail::fcn_vd_2_coef_poly3<TAtom><<<dim3(1), dim3(c_nR), 0, stream[istm]>>>(vatom[istm]);
	// 		}
	// 	}
	// 
	// 	template <class TVctr_i, class TVctr_o>
	// 	typename std::enable_if<is_vctr_gpu_and_vctr_cpu<TVctr_i, TVctr_o>::value 
	// 	&& is_cfloat<Value_type<TVctr_o>>::value && !std::is_same<Value_type<TVctr_i>, Value_type<TVctr_o>>::value, void>::type
	// 	cpy_to_host(Stream<edev_cpu>& stream, TVctr_i &mx_i, TVctr_o &mx_o, 
	// 	Vctr<Value_type<TVctr_i>, edev_cpu> *M_i_h = nullptr)
	// 	{
	// 		Vctr<Value_type<TVctr_i>, edev_cpu> M_h;
	// 		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

	// 		// data transfer from GPU to CPU
	// 		M_i_h->assign(mx_i.begin(), mx_i.end());

	// 		// copy data from CPU to CPU
	// 		cpy_to_host(stream, *M_i_h, mx_o.ptr_32());
	// 	}

	// 	template <class TVctr_i, class TVctr_o>
	// 	typename std::enable_if<is_vctr_gpu_and_vctr_cpu<TVctr_i, TVctr_o>::value 
	// 	&& (!is_cfloat<Value_type<TVctr_o>>::value || std::is_same<Value_type<TVctr_i>, Value_type<TVctr_o>>::value), void>::type
	// 	cpy_to_host(Stream<edev_cpu>& stream, TVctr_i &mx_i, TVctr_o &mx_o, 
	// 	Vctr<Value_type<TVctr_i>, edev_cpu> *M_i_h = nullptr)
	// 	{
	// 		mx_o.assign(mx_i.begin(), mx_i.end());
	// 	}

	// 	template <class TVctr_i, class TVctr_o>
	// 	enable_if_vctr_gpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	// 	add_sc_to_host(Stream<edev_cpu>& stream, Value_type<TVctr_i> w_i, 
	// 	TVctr_i &mx_i, TVctr_o &mx_o, Vctr<Value_type<TVctr_i>, edev_cpu> *M_i_h = nullptr)
	// 	{
	// 		Vctr<Value_type<TVctr_i>, edev_cpu> M_h;
	// 		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

	// 		// data transfer from GPU to CPU
	// 		M_i_h->assign(mx_i.begin(), mx_i.end());

	// 		// add and scale
	// 		mt::add_sc_to_host(stream, w_i, *M_i_h, mx_o.ptr_32());
	// 	}

	// 	template <class TVctr_i, class TVctr_o>
	// 	enable_if_vctr_gpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	// 	add_sc_square_to_host(Stream<edev_cpu>& stream, Value_type<TVctr_o> w_i, 
	// 	TVctr_i &mx_i, TVctr_o &mx_o, Vctr<Value_type<TVctr_i>, edev_cpu> *M_i_h = nullptr)
	// 	{
	// 		Vctr<Value_type<TVctr_i>, edev_cpu> M_h;
	// 		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

	// 		// data transfer from GPU to CPU
	// 		M_i_h->assign(mx_i.begin(), mx_i.end());

	// 		mt::add_sc_square_to_host(stream, w_i, *M_i_h, mx_o.ptr_32());
	// 	}

	// 	template <class TVctr_c_i, class TVctr_r_o, class TVctr_c_o>
	// 	enable_if_vctr_gpu_and_vctr_cpu<TVctr_c_i, TVctr_c_o, void>
	// 	add_sc_m2psi_psi_to_host(Stream<edev_cpu>& stream, Value_type<TVctr_r_o> w_i, 
	// 	TVctr_c_i &psi_i, TVctr_r_o &m2psi_o, TVctr_c_o &psi_o, Vctr<Value_type<TVctr_c_i>, edev_cpu> *psi_i_h = nullptr)
	// 	{
	// 		Vctr<Value_type<TVctr_c_i>, edev_cpu> M_h;
	// 		psi_i_h = (psi_i_h == nullptr)?&M_h:psi_i_h;

	// 		// data transfer from GPU to CPU
	// 		psi_i_h->assign(psi_i.begin(), psi_i.end());

	// 		mt::add_sc_m2psi_psi_to_host(stream, w_i, *psi_i_h, m2psi_o, psi_o);
	// 	}
	// 	/***************************************************************************************/
	// 	// find peak position 1d
	// 	template <class TGrid, class TVctr>
	// 	enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
	// 	fit_max_pos_1d(TGrid& grid_1d, TVctr& Im, Value_type<TVctr> p_i, 
	// 	Value_type<TVctr> sigma_i, Value_type<TVctr> radius_i)
	// 	{
	// 		using T = Value_type<TVctr>;

	// 		Vctr<T, edev_cpu> Im_h = Im;
	// 		return fit_max_pos_1d(grid_1d, Im_h, p_i, sigma_i, radius_i);
	// 	}

	// 	inline
	// 	D_Grid_Blk d_grid_blk(dt_int32 nx, dt_int32 ny, dt_int32 gpu_device, dt_int32 n_thread)
	// 	{
 // 			dt_int32 n_SMs;
	// 		cudaDeviceGetAttribute(&n_SMs, cudaDevAttrMultiProcessorCount, gpu_device);

 // 			dt_int32 n_thr_SMs;
	// 		cudaDeviceGetAttribute(&n_thr_SMs, cudaDevAttrMaxThreadsPerMultiProcessor, gpu_device);

 // 			dt_int32 n_thr_blk;
	// 		cudaDeviceGetAttribute(&n_thr_blk, cudaDevAttrMaxThreadsPerBlock, gpu_device);

	// 		dt_int32 tnr_nxy = c_thr_2d_x*c_thr_2d_y;
	// 		dt_int32 blk_max = 8*max(1, n_thr_SMs/(tnr_nxy*n_thread))*n_SMs;
	// 		dt_int32 blk_x = (ny+c_thr_2d_x-1)/c_thr_2d_x;
	// 		dt_int32 blk_y = (nx+c_thr_2d_y-1)/c_thr_2d_y;
	// 		dt_float64 f = sqrt(dt_float64(blk_max)/dt_float64(blk_x*blk_y));
	// 		blk_x = max(c_blk_x, static_cast<dt_int32>(f*blk_x));
	// 		blk_y = max(c_blk_y, static_cast<dt_int32>(f*blk_y));

	// 		D_Grid_Blk d_grid_blk;
	// 		d_grid_blk.grid = dim3(blk_x, blk_y);
	// 		d_grid_blk.blk = dim3(c_thr_2d_x, c_thr_2d_y);

	// 		return d_grid_blk;
	// 	}
	}

#endif