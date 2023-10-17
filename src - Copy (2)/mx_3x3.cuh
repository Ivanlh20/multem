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

#ifndef MX_3x3_H
	#define MX_3x3_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "r_3d.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <math.h>
		#include <thrust/complex.h>
	#endif

	namespace mt
	{
		 /***************************************************************************************/
		/******************************* forward declarations **********************************/
		 /***************************************************************************************/
		template <class T>
		class Mx_3x3;

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_zero(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_nzero(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		dt_bool is_id_mx(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		T fmin(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		T fmax(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		T det(const Mx_3x3<T>& mx);

		template <class T>
		CGPU_EXEC
		Mx_3x3<T> inv(const Mx_3x3<T>& mx);

		 /***************************************************************************************/
		/******************************* 3x3 matrix *****************************/
		 /***************************************************************************************/

		template <class T>
		class Mx_3x3
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			/* Default constructor */
			CGPU_EXEC
			Mx_3x3(): m_data{T(), T(), T(), T(), T(), T(), T(), T(), T()} {}

			// ! constructor
			CGPU_EXEC
			Mx_3x3(const T& c_11, const T& c_21, const T& c_31, const T& c_12, 
			const T& c_22, const T& c_32, const T& c_13, const T& c_23, const T& c_33) 
			{
				m_data[0] = c_11;
				m_data[1] = c_21;
				m_data[2] = c_31;
				m_data[3] = c_12;
				m_data[4] = c_22;
				m_data[5] = c_32;
				m_data[6] = c_13;
				m_data[7] = c_23;
				m_data[8] = c_33;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Mx_3x3(const U& c_11, const U& c_21, const U& c_31, const U& c_12, 
			const U& c_22, const U& c_32, const U& c_13, const U& c_23, const U& c_33) 
			{
				m_data[0] = T(c_11);
				m_data[1] = T(c_21);
				m_data[2] = T(c_31);
				m_data[3] = T(c_12);
				m_data[4] = T(c_22);
				m_data[5] = T(c_32);
				m_data[6] = T(c_13);
				m_data[7] = T(c_23);
				m_data[8] = T(c_33);
			}

			///* converting constructor */
			//template <class U>
			//CGPU_EXEC
			//Mx_3x3(dt_init_list<U>& list) 
			//{
			//	m_data[0] = T(ptr[0]);
			//	m_data[1] = T(ptr[1]);
			//	m_data[2] = T(ptr[2]);
			//	m_data[3] = T(ptr[3]);
			//	m_data[4] = T(ptr[4]);
			//	m_data[5] = T(ptr[5]);
			//	m_data[6] = T(ptr[6]);
			//	m_data[7] = T(ptr[7]);
			//	m_data[8] = T(ptr[8]);
			//}

			/* Copy constructor */
			CGPU_EXEC
			Mx_3x3(const Mx_3x3<T>& mx)
			{
				*this = mx;
			}

			/* Move constructor */
			template <class U>
			CGPU_EXEC
			Mx_3x3(Mx_3x3<U>&& mx)
			{
				*this = std::move(mx);
			}

			/* converting copy constructor */
			template <class U>
			CGPU_EXEC
			Mx_3x3(const Mx_3x3<U>& mx)
			{
				*this = mx;
			}

			template <class U, class = enable_if_real<U>>
			CGPU_EXEC
			Mx_3x3(U* v) 
			{
				m_data[0] = T(v[0]);
				m_data[1] = T(v[1]);
				m_data[2] = T(v[2]);
				m_data[3] = T(v[3]);
				m_data[4] = T(v[4]);
				m_data[5] = T(v[5]);
				m_data[6] = T(v[6]);
				m_data[7] = T(v[7]);
				m_data[8] = T(v[8]);
			}

			/******************************** assignment operators *********************************/

			/* copy assignment operator */
			CGPU_EXEC
			Mx_3x3<T>& operator=(const Mx_3x3<T>& mx)
			{
				if (this != &mx)
				{
					m_data[0] = mx.m_data[0];
					m_data[1] = mx.m_data[1];
					m_data[2] = mx.m_data[2];
					m_data[3] = mx.m_data[3];
					m_data[4] = mx.m_data[4];
					m_data[5] = mx.m_data[5];
					m_data[6] = mx.m_data[6];
					m_data[7] = mx.m_data[7];
					m_data[8] = mx.m_data[8];
				}

				return *this;
			}

			/* Move assignment operator */
			CGPU_EXEC
			Mx_3x3<T>& operator=(Mx_3x3<T>&& mx)
			{
				if (this != &mx)
				{
					m_data[0] = mx.m_data[0];
					m_data[1] = mx.m_data[1];
					m_data[2] = mx.m_data[2];
					m_data[3] = mx.m_data[3];
					m_data[4] = mx.m_data[4];
					m_data[5] = mx.m_data[5];
					m_data[6] = mx.m_data[6];
					m_data[7] = mx.m_data[7];
					m_data[8] = mx.m_data[8];
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator=(const Mx_3x3<U>& mx)
			{
				m_data[0] = T(mx.m_data[0]);
				m_data[1] = T(mx.m_data[1]);
				m_data[2] = T(mx.m_data[2]);
				m_data[3] = T(mx.m_data[3]);
				m_data[4] = T(mx.m_data[4]);
				m_data[5] = T(mx.m_data[5]);
				m_data[6] = T(mx.m_data[6]);
				m_data[7] = T(mx.m_data[7]);
				m_data[8] = T(mx.m_data[8]);

				return *this;
			}

			/******************* Compound assignment operators *******************/
			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator+=(const Mx_3x3<U>& mx)
			{
				*this = *this + mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator+=(const U& mx)
			{
				*this = *this + mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator-=(const Mx_3x3<U>& mx)
			{
				*this = *this - mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator-=(const U& mx)
			{
				*this = *this - mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator*=(const Mx_3x3<U>& mx)
			{
				*this = *this*mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator*=(const U& mx)
			{
				*this = *this*mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator/=(const Mx_3x3<U>& mx)
			{
				*this = *this/mx;

				return *this;
			}

			template <class U>
			CGPU_EXEC
			Mx_3x3<T>& operator/=(const U& mx)
			{
				*this = *this/mx;

				return *this;
			}

			/************************** Other operators **************************/
			dt_int32 size() const
			{
				return 9;
			}

			CGPU_EXEC_INL
			T& operator[](const dt_int32 idx) { return m_data[idx]; }

			CGPU_EXEC_INL
			const T& operator[](const dt_int32 idx) const { return m_data[idx]; }

			CGPU_EXEC_INL
			T& operator()(const dt_int32 ir, const dt_int32 ic) { return m_data[ir - 1 + 3*(ic - 1)]; }

			CGPU_EXEC_INL
			const T& operator()(const dt_int32 ir, const dt_int32 ic) const { return m_data[ir - 1 + 3*(ic - 1)]; }


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
			Mx_3x3<T> inv() const
			{
				return mt::inv(*this);
			}

		private:
			T m_data[9];
		};

		/******************** Unitary arithmetic operators ********************/
		template <class T>
		CGPU_EXEC
		Mx_3x3<T> operator-(const Mx_3x3<T>& mx)
		{
			return {-mx[0], -mx[1], -mx[2], -mx[3], -mx[4], -mx[5], -mx[6], -mx[7], -mx[8]};
		}

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_zero(const Mx_3x3<T>& mx)
		{
			return fcn_is_zero(mx[0], mx[1], mx[2], mx[3], mx[4], mx[5], mx[6], mx[7], mx[8]);
		}

		template <class T>
		CGPU_EXEC
		dt_bool fcn_is_nzero(const Mx_3x3<T>& mx)
		{
			return fcn_is_nzero(mx[0], mx[1], mx[2], mx[3], mx[4], mx[5], mx[6], mx[7], mx[8]);
		}

		template <class T>
		CGPU_EXEC
		dt_bool is_id_mx(const Mx_3x3<T>& mx)
		{
			return fcn_is_one(mx[0], mx[4], mx[8]) && fcn_is_zero(mx[1], mx[2], mx[3], mx[5], mx[6], mx[7]);
		}

		template <class T>
		CGPU_EXEC
		T fmin(const Mx_3x3<T>& mx)
		{
			return ::fmin(mx[0], ::fmin(mx[1], ::fmin(mx[2], ::fmin(mx[3], ::fmin(mx[4], ::fmin(mx[5], ::fmin(mx[6], ::fmin(mx[7], mx[8]))))))));
		}

		template <class T>
		CGPU_EXEC
		T fmax(const Mx_3x3<T>& mx)
		{
			return ::fmax(mx[0], ::fmax(mx[1], ::fmax(mx[2], ::fmax(mx[3], ::fmax(mx[4], ::fmax(mx[5], ::fmax(mx[6], ::fmax(mx[7], mx[8]))))))));
		}

		template <class T>
		CGPU_EXEC
		T det(const Mx_3x3<T>& mx)
		{
			T val = mx(1, 1)*mx(2, 2)*mx(3, 3) + mx(1, 2)*mx(2, 3)*mx(3, 1) + mx(1, 3)*mx(2, 1)*mx(3, 2);
			val -= mx(1, 3)*mx(2, 2)*mx(3, 1) + mx(1, 2)*mx(2, 1)*mx(3, 3) + mx(1, 1)*mx(2, 3)*mx(3, 2);

			return val;
		}

		template <class T>
		CGPU_EXEC
		Mx_3x3<T> inv(const Mx_3x3<T>& mx)
		{
			Mx_3x3<T> mx_inv;

			mx_inv(1, 1) = mx(2, 2)*mx(3, 3) - mx(2, 3)*mx(3, 2);
			mx_inv(2, 1) = mx(2, 3)*mx(3, 1) - mx(2, 1)*mx(3, 3);
			mx_inv(3, 1) = mx(2, 1)*mx(3, 2) - mx(2, 2)*mx(3, 1);

			mx_inv(1, 2) = mx(1, 3)*mx(3, 2) - mx(1, 2)*mx(3, 3);
			mx_inv(2, 2) = mx(1, 1)*mx(3, 3) - mx(1, 3)*mx(3, 1);
			mx_inv(3, 2) = mx(1, 2)*mx(3, 1) - mx(1, 1)*mx(3, 2);

			mx_inv(1, 3) = mx(1, 2)*mx(2, 3) - mx(1, 3)*mx(2, 2);
			mx_inv(3, 3) = mx(1, 1)*mx(2, 2) - mx(1, 2)*mx(2, 1);
			mx_inv(2, 3) = mx(1, 3)*mx(2, 1) - mx(1, 1)*mx(2, 3);

			return mx_inv/det(mx);
		}

		/************************ arithmetic operators ***********************/
		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> operator+(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3], 
				lhs[4] + rhs[4], lhs[5] + rhs[5], lhs[6] + rhs[6], lhs[7] + rhs[7], lhs[8] + rhs[8]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_3x3<T> operator+(const Mx_3x3<T>& lhs, const U& rhs)
		{
			return {lhs[0] + rhs, lhs[1] + rhs, lhs[2] + rhs, lhs[3] + rhs, 
				lhs[4] + rhs, lhs[5] + rhs, lhs[6] + rhs, lhs[7] + rhs, lhs[8] + rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_3x3<U> operator+(const T& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs + rhs[0], lhs + rhs[1], lhs + rhs[2], lhs + rhs[3], 
				lhs + rhs[4], lhs + rhs[5], lhs + rhs[6], lhs + rhs[7], lhs + rhs[8]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> operator-(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3], 
				lhs[4] - rhs[4], lhs[5] - rhs[5], lhs[6] - rhs[6], lhs[7] - rhs[7], lhs[8] - rhs[8]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_3x3<T> operator-(const Mx_3x3<T>& lhs, const U& rhs)
		{
			return {lhs[0] - rhs, lhs[1] - rhs, lhs[2] - rhs, lhs[3] - rhs, 
				lhs[4] - rhs, lhs[5] - rhs, lhs[6] - rhs, lhs[7] - rhs, lhs[8] - rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_3x3<U> operator-(const T& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs - rhs[0], lhs - rhs[1], lhs - rhs[2], lhs - rhs[3], 
				lhs - rhs[4], lhs - rhs[5], lhs - rhs[6], lhs - rhs[7], lhs - rhs[8]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> operator*(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			Mx_3x3<X> mx;

			mx[0] = lhs[0]*rhs[0] + lhs[3]*rhs[1] + lhs[6]*rhs[2];
			mx[1] = lhs[1]*rhs[0] + lhs[4]*rhs[1] + lhs[7]*rhs[2];
			mx[2] = lhs[2]*rhs[0] + lhs[5]*rhs[1] + lhs[8]*rhs[2];

			mx[3] = lhs[0]*rhs[3] + lhs[3]*rhs[4] + lhs[6]*rhs[5];
			mx[4] = lhs[1]*rhs[3] + lhs[4]*rhs[4] + lhs[7]*rhs[5];
			mx[5] = lhs[2]*rhs[3] + lhs[5]*rhs[4] + lhs[8]*rhs[5];

			mx[6] = lhs[0]*rhs[6] + lhs[3]*rhs[7] + lhs[6]*rhs[8];
			mx[7] = lhs[1]*rhs[6] + lhs[4]*rhs[7] + lhs[7]*rhs[8];
			mx[8] = lhs[2]*rhs[6] + lhs[5]*rhs[7] + lhs[8]*rhs[8];

			return mx;
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_3x3<T> operator*(const Mx_3x3<T>& lhs, const U& rhs)
		{
			return {lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs, 
				lhs[4]*rhs, lhs[5]*rhs, lhs[6]*rhs, lhs[7]*rhs, lhs[8]*rhs};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_3x3<U> operator*(const T& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs*rhs[0], lhs*rhs[1], lhs*rhs[2], lhs*rhs[3], 
				lhs*rhs[4], lhs*rhs[5], lhs*rhs[6], lhs*rhs[7], lhs*rhs[8]};
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		R_3d<U> operator*(const Mx_3x3<T>& lhs, const R_3d<U>& rhs)
		{
			return {lhs[0]*rhs.x + lhs[3]*rhs.y + lhs[6]*rhs.z, 
				lhs[1]*rhs.x + lhs[4]*rhs.y + lhs[7]*rhs.z, 
				lhs[2]*rhs.x + lhs[5]*rhs.y + lhs[8]*rhs.z};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		R_3d<T> operator*(const R_3d<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return {lhs.x*rhs[0] + lhs.y*rhs[1] + lhs.z*rhs[2], 
				lhs.x*rhs[3] + lhs.y*rhs[4] + lhs.z*rhs[5], 
				lhs.x*rhs[6] + lhs.y*rhs[7] + lhs.z*rhs[8]};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> operator/(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return lhs*inv(rhs);
		}

		template <class T, class U, class = enable_if_real<U>>
		CGPU_EXEC
		Mx_3x3<T> operator/(const Mx_3x3<T>& lhs, const U& rhs)
		{
			return {fcn_div(lhs[0], rhs), fcn_div(lhs[1], rhs), fcn_div(lhs[2], rhs), fcn_div(lhs[3], rhs), 
				fcn_div(lhs[4], rhs), fcn_div(lhs[5], rhs), fcn_div(lhs[6], rhs), fcn_div(lhs[7], rhs), fcn_div(lhs[8], rhs)};
		}

		template <class T, class U, class = enable_if_real<T>>
		CGPU_EXEC
		Mx_3x3<U> operator/(const T& lhs, const Mx_3x3<U>& rhs)
		{
			return {fcn_div(lhs, rhs[0]), fcn_div(lhs, rhs[1]), fcn_div(lhs, rhs[2]), fcn_div(lhs, rhs[3]), 
				fcn_div(lhs, rhs[4]), fcn_div(lhs, rhs[5]), fcn_div(lhs, rhs[6]), fcn_div(lhs, rhs[7]), fcn_div(lhs, rhs[8])};
		}

		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> fmin(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return {::fmin(lhs[0], rhs[0]), ::fmin(lhs[1], rhs[1]), ::fmin(lhs[2], rhs[2]), ::fmin(lhs[3], rhs[3]), 
				::fmin(lhs[4], rhs[4]), ::fmin(lhs[5], rhs[5]), ::fmin(lhs[6], rhs[6]), ::fmin(lhs[7], rhs[7]), ::fmin(lhs[8], rhs[8])};
		}
	
		template <class T, class U, class X = sel_lg_type<T, U>>
		CGPU_EXEC
		Mx_3x3<X> fmax(const Mx_3x3<T>& lhs, const Mx_3x3<U>& rhs)
		{
			return {::fmax(lhs[0], rhs[0]), ::fmax(lhs[1], rhs[1]), ::fmax(lhs[2], rhs[2]), ::fmax(lhs[3], rhs[3]), 
				::fmax(lhs[4], rhs[4]), ::fmax(lhs[5], rhs[5]), ::fmax(lhs[6], rhs[6]), ::fmax(lhs[7], rhs[7]), ::fmax(lhs[8], rhs[8])};
		}


		/***************************************************************************************/
		/******************************** initializer list *************************************/
		/***************************************************************************************/
		using dt_init_list_mx_3x3_f64 = std::initializer_list<Mx_3x3<dt_float64>>;

		/************************** other operators **************************/
		template <class T>
		CGPU_EXEC 
		Mx_3x3<T> fcn_rot_mx_3d(const T& theta, const R_3d<T>& u0)
		{
			T alpha = T(1)-cos(theta);
			alpha = fcn_is_zero(alpha)?0:alpha;
			T beta = sin(theta);
			beta = fcn_is_zero(beta)?0:beta;

			return {T(1) + alpha*(u0.x*u0.x-T(1)), u0.y*u0.x*alpha + u0.z*beta, u0.z*u0.x*alpha - u0.y*beta, 
			u0.x*u0.y*alpha - u0.z*beta, T(1) + alpha*(u0.y*u0.y-1), u0.z*u0.y*alpha + u0.x*beta, 
			u0.x*u0.z*alpha + u0.y*beta, u0.y*u0.z*alpha - u0.x*beta, T(1) + alpha*(u0.z*u0.z-1)};
		}
	}
#endif