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

#ifndef MX_2x2_H
	#define MX_2x2_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <math.h>
		#include <thrust/complex.h>
	#endif

	#include "r_2d.cuh"

	namespace mt
	{
		/***************************************************************************************/
		/******************************* forward declarations **********************************/
		/***************************************************************************************/
		template <class T>
		class Mx_2x2;

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_zero(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_nzero(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		dt_bool is_id_mx(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		T fmin(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		T fmax(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		T det(const Mx_2x2<T>& mx);

		template <class T>
		CGPU_EXEC
		Mx_2x2<T> inv(const Mx_2x2<T>& mx);

		/***************************************************************************************/
		/******************************* 2x2 matrix *****************************/
		/***************************************************************************************/
	
		template <class T>
		class Mx_2x2
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			/* Default constructor */
			CGPU_EXEC
			Mx_2x2(): m_data{0, 0, 0, 0} {}

			// ! constructor
			CGPU_EXEC
			Mx_2x2(const T& c_11, const T& c_21, const T& c_12, const T& c_22) 
			{
				m_data[0] = c_11;
				m_data[1] = c_21;
				m_data[2] = c_12;
				m_data[3] = c_22;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Mx_2x2(const U& c_11, const U& c_21, const U& c_12, const U& c_22) 
			{
				m_data[0] = T(c_11);
				m_data[1] = T(c_21);
				m_data[2] = T(c_12);
				m_data[3] = T(c_22);
			}

			/* Copy constructor */
			CGPU_EXEC
			Mx_2x2(const Mx_2x2<T>& mx)
			{
				*this = mx;
			}

			/* Move constructor */
			template <class U>
			CGPU_EXEC
			Mx_2x2(Mx_2x2<U>&& mx)
			{
				*this = std::move(mx);
			}

			/* converting copy constructor */
			template <class U>
			CGPU_EXEC
			Mx_2x2(const Mx_2x2<U>& mx)
			{
				*this = mx;
			}

			template <class U, class = enable_if_real<U>>
			CGPU_EXEC
			Mx_2x2(U* v) 
			{
				m_data[0] = T(v[0]);
				m_data[1] = T(v[1]);
				m_data[2] = T(v[2]);
				m_data[3] = T(v[3]);
			}

			/******************************** assignment operators *********************************/

			/* copy assignment operator */
			CGPU_EXEC
			Mx_2x2<T>& operator=(const Mx_2x2<T>& mx)
			{
				if (this != &mx)
				{
					m_data[0] = mx.m_data[0];
					m_data[1] = mx.m_data[1];
					m_data[2] = mx.m_data[2];
					m_data[3] = mx.m_data[3];
				}

				return *this;
			}

			/* Move assignment operator */
			CGPU_EXEC
			Mx_2x2<T>& operator=(Mx_2x2<T>&& mx)
			{
				if (this != &mx)
				{
					m_data[0] = mx.m_data[0];
					m_data[1] = mx.m_data[1];
					m_data[2] = mx.m_data[2];
					m_data[3] = mx.m_data[3];
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator=(const Mx_2x2<U>& mx)
			{
				m_data[0] = T(mx.m_data[0]);
				m_data[1] = T(mx.m_data[1]);
				m_data[2] = T(mx.m_data[2]);
				m_data[3] = T(mx.m_data[3]);

				return *this;
			}

			/******************* Compound assignment operators *******************/
			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator+=(const Mx_2x2<U>& mx)
			{
				*this = *this + mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator+=(const U& mx)
			{
				*this = *this + mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator-=(const Mx_2x2<U>& mx)
			{
				*this = *this - mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator-=(const U& mx)
			{
				*this = *this - mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator*=(const Mx_2x2<U>& mx)
			{
				*this = *this*mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator*=(const U& mx)
			{
				*this = *this*mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator/=(const Mx_2x2<U>& mx)
			{
				*this = *this/mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_2x2<T>& operator/=(const U& mx)
			{
				*this = *this/mx;

				return *this;
			}

			/************************** Other operators **************************/
			dt_int32 size() const
			{
				return 4;
			}

			CGPU_EXEC_INL
			T& operator[](const dt_int32 idx) { return m_data[idx]; }

			CGPU_EXEC_INL
			const T& operator[](const dt_int32 idx) const { return m_data[idx]; }

			CGPU_EXEC_INL
			T& operator()(const dt_int32 ir, const dt_int32 ic) { return m_data[ir - 1 + 2*(ic - 1)]; }

			CGPU_EXEC_INL
			const T& operator()(const dt_int32 ir, const dt_int32 ic) const { return m_data[ir - 1 + 2*(ic - 1)]; }


			CGPU_EXEC
			dt_bool fcn_is_zero() const
			{
				return mt::fcn_is_zero(*this);
			}

			CGPU_EXEC
			dt_bool is_nzero() const
			{
				return mt::fcn_is_nzero(*this);
			}

			CGPU_EXEC
			dt_bool is_id_mx() const
			{
				return mt::is_id_mx(*this);
			}

			CGPU_EXEC
			T min() const
			{
				return mt::fmin(*this);
			}

			CGPU_EXEC
			T max() const
			{
				return mt::fmax(*this);
			}

			CGPU_EXEC
			T det() const
			{
				return mt::det(*this);
			}

			CGPU_EXEC
			Mx_2x2<T> inv() const
			{
				return mt::inv(*this);
			}

		private:
			T m_data[4];
		};

		/******************** Unitary arithmetic operators ********************/
		template <class T>
		CGPU_EXEC
		Mx_2x2<T> operator-(const Mx_2x2<T>& mx)
		{
			return {-mx[0], -mx[1], -mx[2], -mx[3]};
		}

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_zero(const Mx_2x2<T>& mx)
		{
			return fcn_is_zero(mx[0], mx[1], mx[2], mx[3]);
		}

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_nzero(const Mx_2x2<T>& mx)
		{
			return fcn_is_nzero(mx[0], mx[1], mx[2], mx[3]);
		}

		template <class T>
		CGPU_EXEC
		dt_bool is_id_mx(const Mx_2x2<T>& mx)
		{
			return fcn_is_one(mx[0], mx[3]) && fcn_is_zero(mx[1], mx[2]);
		}

		template <class T>
		CGPU_EXEC
		T fmin(const Mx_2x2<T>& mx)
		{
			return ::fmin(mx[0], ::fmin(mx[1], ::fmin(mx[2], mx[3])));
		}

		template <class T>
		CGPU_EXEC
		T fmax(const Mx_2x2<T>& mx)
		{
			return ::fmax(mx[0], ::fmax(mx[1], ::fmax(mx[2], mx[3])));
		}

		template <class T>
		CGPU_EXEC
		T det(const Mx_2x2<T>& mx)
		{
			return mx(1, 1)*mx(2, 2) - mx(1, 2)*mx(2, 1);
		}

		template <class T>
		CGPU_EXEC
		Mx_2x2<T> inv(const Mx_2x2<T>& mx)
		{
			Mx_2x2<T> mx_inv{mx(2, 2), -mx(2, 1), -mx(1, 2), mx(1, 1)};

			return mx_inv/det(mx);
		}

		/************************ arithmetic operators ***********************/
		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> operator+(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_2x2<T> operator+(const Mx_2x2<T>& lhs, const U& rhs)
		{
			return {lhs[0] + rhs, lhs[1] + rhs, lhs[2] + rhs, lhs[3] + rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_2x2<U> operator+(const T& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs + rhs[0], lhs + rhs[1], lhs + rhs[2], lhs + rhs[3]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> operator-(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_2x2<T> operator-(const Mx_2x2<T>& lhs, const U& rhs)
		{
			return {lhs[0] - rhs, lhs[1] - rhs, lhs[2] - rhs, lhs[3] -rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_2x2<U> operator-(const T& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs - rhs[0], lhs - rhs[1], lhs - rhs[2], lhs - rhs[3]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> operator*(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			Mx_2x2<X> mx;

			mx[0] = lhs[0]*rhs[0] + lhs[2]*rhs[1];
			mx[1] = lhs[1]*rhs[0] + lhs[3]*rhs[1];
		
			mx[2] = lhs[0]*rhs[2] + lhs[2]*rhs[3];
			mx[3] = lhs[1]*rhs[2] + lhs[3]*rhs[3];

			return mx;
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_2x2<T> operator*(const Mx_2x2<T>& lhs, const U& rhs)
		{
			return {lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_2x2<U> operator*(const T& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs*rhs[0], lhs*rhs[1], lhs*rhs[2], lhs*rhs[3]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		R_2d<U> operator*(const Mx_2x2<T>& lhs, const R_2d<U>& rhs)
		{
			return {lhs[0]*rhs.x + lhs[2]*rhs.y, lhs[1]*rhs.x + lhs[3]*rhs.y};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		R_2d<T> operator*(const R_2d<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return {lhs.x*rhs[0] + lhs.y*rhs[1], lhs.x*rhs[2] + lhs.y*rhs[3]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> operator/(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return lhs*inv(rhs);
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_2x2<T> operator/(const Mx_2x2<T>& lhs, const U& rhs)
		{
			return {fcn_div(lhs[0], rhs), fcn_div(lhs[1], rhs), fcn_div(lhs[2], rhs), fcn_div(lhs[3], rhs)};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_2x2<U> operator/(const T& lhs, const Mx_2x2<U>& rhs)
		{
			return {fcn_div(lhs, rhs[0]), fcn_div(lhs, rhs[1]), fcn_div(lhs, rhs[2]), fcn_div(lhs, rhs[3])};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> fmin(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return {::fmin(lhs[0], rhs[0]), ::fmin(lhs[1], rhs[1]), ::fmin(lhs[2], rhs[2]), ::fmin(lhs[3], rhs[3])};
		}
	
		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_2x2<X> fmax(const Mx_2x2<T>& lhs, const Mx_2x2<U>& rhs)
		{
			return {::fmax(lhs[0], rhs[0]), ::fmax(lhs[1], rhs[1]), ::fmax(lhs[2], rhs[2]), ::fmax(lhs[3], rhs[3])};
		}

		/***************************************************************************************/
		/******************************** initializer list *************************************/
		/***************************************************************************************/
		using dt_init_list_mx_2x2_f64 = std::initializer_list<Mx_2x2<dt_float64>>;

		/************************** other operators **************************/
		template <class T>
		CGPU_EXEC 
		Mx_2x2<T> fcn_rot_mx_2d(const T& theta)
		{
			T cos_theta, sin_theta;
			sincos(theta, &sin_theta, &cos_theta);

			return {cos_theta, sin_theta, -sin_theta, cos_theta};
		}

	}
#endif