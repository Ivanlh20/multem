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

#include <type_traits>
#include <complex>

#include "const_enum.h"
#include "math_mt.h"
#include "type_traits_gen.h"
#include "fcns_cgpu_gen.h"
#include "r_2d.h"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <math.h>
	#include <thrust/complex.h>
#endif

/* R_3d */
namespace mt
{
	template <class T>
	using R_3d = R_xtd<T, edim_3>;

 	template <class T>
	class R_xtd<T, edim_3>
	{
	public:
		using value_type = T;

		T x;
		T y;
		T z;

		/************************************* constructors ************************************/
		/* Default constructor */
		CGPU_EXEC
		R_xtd();

		/* converting constructor by initializer list */
		template <class U>
		CPU_EXEC
		R_xtd(const dt_init_list<U>& list);

		/* constructor by pointer */
		template <class U, class = enable_if_real<U>>
		CPU_EXEC
		R_xtd(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0);

		/* constructor by pointer */
		template <class U, class = enable_if_real<U>>
		CPU_EXEC
		R_xtd(U* v);

		/* constructor by pointer */
		template <class U, class = enable_if_real<U>>
		CPU_EXEC
		R_xtd(U* v, const dt_int32& n_v);

		/* constructor a R_3d with an x part */
		CGPU_EXEC
		R_xtd(const T& x);

		/* converting constructor: R_3d with an x part */
		template <class U>
		CGPU_EXEC
		R_xtd(const U& x);

		/* constructor a R_3d from its x and y parts */
		CGPU_EXEC
		R_xtd(const T& x, const T& y);

		/* converting constructor: R_3d from its x and y parts */
		template <class U>
		CGPU_EXEC
		R_xtd(const U& x, const U& y);

		/* constructor a R_3d from its x, y and z parts */
		CGPU_EXEC
		R_xtd(const T& x, const T& y, const T& z);

		/* converting constructor: R_3d from its x, y and z parts */
		template <class U>
		CGPU_EXEC
		R_xtd(const U& x, const U& y, const U& z);

		/* copy constructor */
		CGPU_EXEC
		R_xtd(const R_xtd<T, edim_3>& r);

		/* Move constructor */
		CGPU_EXEC
		R_xtd(R_xtd<T, edim_3>&& r);

		/* converting copy constructor */
		template <class U>
		CGPU_EXEC
		R_xtd(const R_xtd<U, edim_3>& r);

		/******************************** assignment operators *********************************/
		/* Assign x and set y to 0 */
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const T& x);

		/* Convert and assign x and set y to 0 */
		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const U& x);

		/* copy assignment operator */
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const R_xtd<T, edim_3>& r);

		/* Move assignment operator */
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(R_xtd<T, edim_3>&& r);

		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const R_xtd<U, edim_3>& r);

		/* Assignment operator from R_2d */
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const R_2d<T>& r);

		/* converting assignment operator from R_2d */
		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const R_2d<U>& r);

#ifdef __CUDACC__
		/* Assignment operator from thrust::complex */
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const thrust::complex<T>& r);

		/* converting assignment operator from thrust::complex */
		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator=(const thrust::complex<U>& r);
#endif

		/* Assignment operator from std::complex */
		CPU_EXEC
		R_xtd<T, edim_3>& operator=(const std::complex<T>& r);

		/* converting assignment operator from std::complex */
		template <class U> 
		CPU_EXEC
		R_xtd<T, edim_3>& operator=(const std::complex<U>& r);

		/******************* Compound assignment operators *******************/

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator+=(const R_xtd<U, edim_3>& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator+=(const U& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator-=(const R_xtd<U, edim_3>& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator-=(const U& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator*=(const R_xtd<U, edim_3>& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator*=(const U& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator/=(const R_xtd<U, edim_3>& r);

		template <class U>
		CGPU_EXEC
		R_xtd<T, edim_3>& operator/=(const U& r);

		/************************** Other operators **************************/
		CGPU_EXEC
		dt_bool fcn_is_zero() const;

		CGPU_EXEC
		dt_bool is_nzero() const;

		CGPU_EXEC
		dt_bool is_one() const;

		CGPU_EXEC
		R_2d<T> floor() const;

		CGPU_EXEC
		R_2d<T> ceil() const;

		CGPU_EXEC
		R_2d<T> round() const;

		CGPU_EXEC
		T min() const;

		CGPU_EXEC
		T max() const;

		CGPU_EXEC
		T norm_2() const;

		CGPU_EXEC
		T norm() const;

		CGPU_EXEC
		void normalize();

		CGPU_EXEC
		void fill(const T& v);
	};
}

#include "../src/r_3d.inl"