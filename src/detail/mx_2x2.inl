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

#include "mx_2x2.h"

/* forward declarations */
namespace mt
{
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
}

/* 2x2 matrix */
namespace mt
{
	/************************************* constructors ************************************/
	/* Default constructor */
	template <class T>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(): m_data{0, 0, 0, 0} {}

	/* constructor */
	template <class T>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(const T& c_11, const T& c_21, const T& c_12, const T& c_22) 
	{
		m_data[0] = c_11;
		m_data[1] = c_21;
		m_data[2] = c_12;
		m_data[3] = c_22;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(const U& c_11, const U& c_21, const U& c_12, const U& c_22) 
	{
		m_data[0] = T(c_11);
		m_data[1] = T(c_21);
		m_data[2] = T(c_12);
		m_data[3] = T(c_22);
	}

	/* Copy constructor */
	template <class T>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(const Mx_2x2<T>& mx)
	{
		*this = mx;
	}

	/* Move constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(Mx_2x2<U>&& mx)
	{
		*this = std::move(mx);
	}

	/* converting copy constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(const Mx_2x2<U>& mx)
	{
		*this = mx;
	}

	template <class T>
	template <class U, class>
	CGPU_EXEC
	Mx_2x2<T>::Mx_2x2(U* v) 
	{
		m_data[0] = T(v[0]);
		m_data[1] = T(v[1]);
		m_data[2] = T(v[2]);
		m_data[3] = T(v[3]);
	}

	/******************************** assignment operators *********************************/

	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator=(const Mx_2x2<T>& mx)
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
	template <class T>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator=(Mx_2x2<T>&& mx)
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
	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator=(const Mx_2x2<U>& mx)
	{
		m_data[0] = T(mx.m_data[0]);
		m_data[1] = T(mx.m_data[1]);
		m_data[2] = T(mx.m_data[2]);
		m_data[3] = T(mx.m_data[3]);

		return *this;
	}

	/******************* Compound assignment operators *******************/
	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator+=(const Mx_2x2<U>& mx)
	{
		*this = *this + mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator+=(const U& mx)
	{
		*this = *this + mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator-=(const Mx_2x2<U>& mx)
	{
		*this = *this - mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator-=(const U& mx)
	{
		*this = *this - mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator*=(const Mx_2x2<U>& mx)
	{
		*this = *this*mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator*=(const U& mx)
	{
		*this = *this*mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator/=(const Mx_2x2<U>& mx)
	{
		*this = *this/mx;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	Mx_2x2<T>& Mx_2x2<T>::operator/=(const U& mx)
	{
		*this = *this/mx;

		return *this;
	}

	/************************** Other operators **************************/
	template <class T>
	dt_int32 Mx_2x2<T>::size() const
	{
		return 4;
	}

	template <class T>
	CGPU_EXEC_INL
	T& Mx_2x2<T>::operator[](const dt_int32 idx) { return m_data[idx]; }

	template <class T>
	CGPU_EXEC_INL
	const T& Mx_2x2<T>::operator[](const dt_int32 idx) const { return m_data[idx]; }

	template <class T>
	CGPU_EXEC_INL
	T& Mx_2x2<T>::operator()(const dt_int32 ir, const dt_int32 ic) { return m_data[ir - 1 + 2*(ic - 1)]; }

	template <class T>
	CGPU_EXEC_INL
	const T& Mx_2x2<T>::operator()(const dt_int32 ir, const dt_int32 ic) const { return m_data[ir - 1 + 2*(ic - 1)]; }

	template <class T>
	CGPU_EXEC
	dt_bool Mx_2x2<T>::fcn_is_zero() const
	{
		return mt::fcn_is_zero(*this);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Mx_2x2<T>::is_nzero() const
	{
		return mt::fcn_is_nzero(*this);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Mx_2x2<T>::is_id_mx() const
	{
		return mt::is_id_mx(*this);
	}

	template <class T>
	CGPU_EXEC
	T Mx_2x2<T>::min() const
	{
		return mt::fmin(*this);
	}

	template <class T>
	CGPU_EXEC
	T Mx_2x2<T>::max() const
	{
		return mt::fmax(*this);
	}

	template <class T>
	CGPU_EXEC
	T Mx_2x2<T>::det() const
	{
		return mt::det(*this);
	}

	template <class T>
	CGPU_EXEC
	Mx_2x2<T> Mx_2x2<T>::inv() const
	{
		return mt::inv(*this);
	}
}

/* unitary operators */
namespace mt
{
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
}

/* binary operators */
namespace mt
{
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
}

/* traits */
namespace mt
{
	/***************************************************************************************/
	/******************************** initializer list *************************************/
	/***************************************************************************************/
	using dt_init_list_mx_2x2_f64 = std::initializer_list<Mx_2x2<dt_float64>>;
}

/* other operators */
namespace mt
{
	template <class T>
	CGPU_EXEC 
	Mx_2x2<T> fcn_rot_mx_2d(const T& theta)
	{
		T cos_theta, sin_theta;
		sincos(theta, &sin_theta, &cos_theta);

		return {cos_theta, sin_theta, -sin_theta, cos_theta};
	}
}