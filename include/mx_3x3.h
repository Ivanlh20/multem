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

#include "const_enum.h"
#include "math_mt.h"
#include "type_traits_gen.h"
#include "r_3d.h"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <math.h>
	#include <thrust/complex.h>
#endif

namespace mt
{
	/* 3x3 matrix */
	template <class T>
	class Mx_3x3
	{
	public:
		using value_type = T;

		/************************************* constructors ************************************/
		/* Default constructor */
		CGPU_EXEC
		Mx_3x3();

		/* constructor */
		CGPU_EXEC
		Mx_3x3(const T& c_11, const T& c_21, const T& c_31, const T& c_12, 
		const T& c_22, const T& c_32, const T& c_13, const T& c_23, const T& c_33);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Mx_3x3(const U& c_11, const U& c_21, const U& c_31, const U& c_12, 
		const U& c_22, const U& c_32, const U& c_13, const U& c_23, const U& c_33);


		/* Copy constructor */
		CGPU_EXEC
		Mx_3x3(const Mx_3x3<T>& mx);

		/* Move constructor */
		template <class U>
		CGPU_EXEC
		Mx_3x3(Mx_3x3<U>&& mx);

		/* converting copy constructor */
		template <class U>
		CGPU_EXEC
		Mx_3x3(const Mx_3x3<U>& mx);

		template <class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_3x3(U* v);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Mx_3x3<T>& operator=(const Mx_3x3<T>& mx);

		/* Move assignment operator */
		CGPU_EXEC
		Mx_3x3<T>& operator=(Mx_3x3<T>&& mx);

		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator=(const Mx_3x3<U>& mx);

		/******************* Compound assignment operators *******************/
		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator+=(const Mx_3x3<U>& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator+=(const U& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator-=(const Mx_3x3<U>& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator-=(const U& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator*=(const Mx_3x3<U>& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator*=(const U& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator/=(const Mx_3x3<U>& mx);

		template <class U>
		CGPU_EXEC
		Mx_3x3<T>& operator/=(const U& mx);

		/************************** Other operators **************************/
		dt_int32 size() const;

		CGPU_EXEC_INL
		T& operator[](const dt_int32 idx);

		CGPU_EXEC_INL
		const T& operator[](const dt_int32 idx) const;

		CGPU_EXEC_INL
		T& operator()(const dt_int32 ir, const dt_int32 ic);

		CGPU_EXEC_INL
		const T& operator()(const dt_int32 ir, const dt_int32 ic) const;

		CGPU_EXEC
		dt_bool fcn_is_zero() const;

		CGPU_EXEC
		dt_bool is_nzero() const;

		CGPU_EXEC
		dt_bool is_id_mx() const;

		CGPU_EXEC
		T min() const;

		CGPU_EXEC
		T max() const;

		CGPU_EXEC
		T det() const;

		CGPU_EXEC
		Mx_3x3<T> inv() const;

	private:
		T m_data[9];
	};
}

#include "../src/mx_3x3.inl"