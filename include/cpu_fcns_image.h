/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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
#pragma once

#include <thread>
#include <algorithm>
#include <deque>

#include "math_mt.h"
#include "type_traits_gen.h"
#include "vctr_cpu.h"
#include "stream_cpu.h"
#include "r_2d.h"
#include "cgpu_functors.h"

namespace mt
{
	/* threshold */
	namespace cpu_fcn_image
	{
		/* otsu threshold*/
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_otsu_thr(const TVctr& vctr, const dt_uint32& n_bins);

		/* matrix binarization */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_binarize_mx(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr);

		/* maximum threshold */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_threshold_max(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr);

		/* maximum threshold */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_threshold_min(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr);
	}

	/* morphological operations */
	namespace cpu_fcn_image
	{
		/* gray value data dilation */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_dilate(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr);
		
		/* gray value data erosion */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_erode(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr);

		/* gray value data opening */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_open(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr);

		/* gray value data closing */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_close(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr);

		/* gray value datatophat */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_tophat(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr);
	}

	/* filter wiener */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);
		
		/* 2d by col */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_2d_bc(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);
	}

	/* filter mean */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_mean_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_mean_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);
	}

	/* filter median */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		/* 2d by col */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_2d_bc(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr);
	}

	/* filter median by position */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr, class TVctr_idx>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr, TVctr_idx, void>
		fltr_median_1d_pos(TVctr& mx_i, dt_int32 n_kr, TVctr_idx&& vctr_idx, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		// 2d for specific points
		template <class TVctr, class TVctr_idx>
		enable_if_vctr_cpu_r_2d_and_vctr_cpu<TVctr_idx, TVctr, void>
		fltr_median_2d_pos(TVctr& mx_i, dt_int32 n_kr, TVctr_idx&& vctr_idx, TVctr& mx_o, Stream_cpu* pstream = nullptr);
	}

	/* filter poisson noise */
	namespace cpu_fcn_image
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_1d(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_2d_bc(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr);

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_2d(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr);
	}


	/* peak signal to noise ratio */
	namespace cpu_fcn_image
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_get_psnr(TVctr& mx_i, TVctr& mx_r, Stream_cpu* pstream = nullptr);
	}

	/* others */
	namespace cpu_fcn_image
	{
		//// scale_image
		//template <class TVctr>
		//TVctr scale_image_mean(Stream<edev_cpu>& stream, dt_int32 ny_i, dt_int32 nx_i, TVctr& Im_i, Value_type<TVctr> shrink_factor, dt_int32& ny_o, dt_int32& nx_o);
	}
}

#include "../src/cpu_fcns_image.inl"