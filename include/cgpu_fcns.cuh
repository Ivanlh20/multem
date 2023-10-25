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

#ifndef CGPU_FCNS_H
	#define CGPU_FCNS_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include <type_traits>
	#include <algorithm>

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "type_traits_gen.h"
	#include "cgpu_vctr.cuh"

	#include <thrust/complex.h>
	#include <thrust/swap.h>
	#include <thrust/extrema.h>
	#include <thrust/transform.h>
	#include <thrust/reduce.h>
	#include <thrust/transform_reduce.h>
	#include <thrust/functional.h>
	#include <thrust/sort.h>
	#include <thrust/iterator/zip_iterator.h>
	#include <thrust/tuple.h>
	#include <thrust/execution_policy.h>

	namespace mt
	{
		// snr using squared difference
 		template <class TVctr>
		Value_type_r<TVctr> fcn_snr_ma2(TVctr& mx_i, TVctr& M_r)
		{
			using T = Value_type<TVctr>;

			auto ma1_r = fcn_variance(M_r);

			TVctr M_n(mx_i.size());
			for(auto ik = 0; ik<M_n.size(); ik++)
			{
				M_n[ik] = mx_i[ik] - M_r[ik];
			}

			auto ma1_n = fcn_variance(M_n);

			return sqrt(ma1_r/ma1_n);

		}

		// snr using absolute difference
 		template <class TVctr>
		Value_type_r<TVctr> fcn_snr_ma1(TVctr& mx_i, TVctr& M_r)
		{
			using T = Value_type<TVctr>;

			auto ma1_r = fcn_moment_1a(M_r);

			TVctr M_n(mx_i.size());
			for(auto ik = 0; ik<M_n.size(); ik++)
			{
				M_n[ik] = mx_i[ik] - M_r[ik];
			}

			auto ma1_n = fcn_moment_1a(M_n);

			return ma1_r/ma1_n;

		}
	}

#endif
