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

#include "kahan_sum.h"

/* Kahan summation */
namespace mt
{
	/************************************* constructors ************************************/
	/* Default constructor */
	template <class T>
	CGPU_EXEC
	KS<T>::KS(): m_sum(0), m_error(0) {}

	/* Copy constructor */
	template <class T>
	CGPU_EXEC
	KS<T>::KS(const T& sum, T error): m_sum(sum), m_error(error) {}

	/* Copy constructor */
	template <class T>
	CGPU_EXEC
	KS<T>::KS(const KS<T>& ks)
	{
		*this = ks;
	}

	/* Move constructor */
	template <class T>
	CGPU_EXEC
	KS<T>::KS(KS<T>&& ks)
	{
		*this = std::move(ks);
	}

	/* converting copy constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	KS<T>::KS(const KS<U>& ks)
	{
		*this = ks;
	}

	/******************************** assignment operators *********************************/
	/* constructor from r */
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator=(const T& r)
	{
		m_sum = r;
		m_error = 0;

		return *this;
	}

	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator=(const KS<T>& ks)
	{
		m_sum = ks.m_sum;
		m_error = ks.m_error;

		return *this;
	}

	/* Move assignment operator */
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator=(KS<T>&& ks)
	{
		m_sum = std::move(ks.m_sum);
		m_error = std::move(ks.m_error);

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	KS<T>& KS<T>::operator=(const KS<U>& ks)
	{
		m_sum = T(ks.m_sum);
		m_error = T(ks.m_error);

		return *this;
	}

	// user define conversion
	template <class T>
	CGPU_EXEC
	T KS<T>::operator()()
	{
		return m_sum;
	}
			

	// user define conversion
	template <class T>
	CGPU_EXEC
	KS<T>::operator T()
	{
		return m_sum;
	}

	/******************* Compound assignment operators *******************/
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator+=(T r)
	{
		this->add(r);

		return *this;
	}
		
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator+=(const KS<T>& r)
	{
		this->add(r.m_sum);

		return *this;
	}
		
	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator-=(T r)
	{
		this->add(-r);

		return *this;
	}	

	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator-=(const KS<T>& r)
	{
		this->add(-r.m_sum);

		return *this;
	}

	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator/=(T r)
	{
		*this = *this/r;

		return *this;
	}

	template <class T>
	CGPU_EXEC
	KS<T>& KS<T>::operator/=(const KS<T>& r)
	{
		*this = *this/r;

		return *this;
	}

	template <class T>
	CGPU_EXEC
	void KS<T>::add(T r)
	{
		r = r - m_error;
		T t = m_sum + r;			// m_sum is big, r small, so low-order digits of r are lost.
		m_error = (t-m_sum)-r;		// (t - m_sum) cancels the high-order part of r;subtracting r recovers negative (low part of r)
		m_sum = t;					// Algebraically, m_error should always be zero. Beware overly-aggressive optimizing compilers!

		// Next time around, the lost low part will be added to y in a fresh attempt
	}
}

/* unitary operators */
namespace mt
{
	template <class T>
	CGPU_EXEC
	KS<T> operator-(const KS<T>& r)
	{
		return {-r.m_sum, r.m_error};
	}
}

/* binary operators */
namespace mt
{
	template <class T>
	CGPU_EXEC
	KS<T> operator+(const KS<T>& lhs, const KS<T>& rhs)
	{
		KS<T> ks(lhs);
		ks += rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	KS<T> operator+(const KS<T>& lhs, const T& rhs)
	{
		KS<T> ks(lhs);
		ks += rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	KS<T> operator+(const T& lhs, const KS<T>& rhs)
	{
		KS<T> ks(lhs);
		ks += rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	KS<T> operator-(const KS<T>& lhs, const KS<T>& rhs)
	{
		KS<T> ks(lhs);
		ks -= rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	KS<T> operator-(const KS<T>& lhs, const T& rhs)
	{
		KS<T> ks(lhs);
		ks -= rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	KS<T> operator-(const T& lhs, const KS<T>& rhs)
	{
		KS<T> ks(lhs);
		ks -= rhs;
		return ks;
	}

	template <class T>
	CGPU_EXEC
	T operator/(const KS<T>& lhs, const KS<T>& rhs)
	{
		return lhs.m_sum/rhs.m_sum;
	}

	template <class T>
	CGPU_EXEC
	T operator/(const KS<T>& lhs, const T& rhs)
	{
		return lhs.m_sum/rhs;
	}

	template <class T>
	CGPU_EXEC
	T operator/(const T& lhs, const KS<T>& rhs)
	{
		return lhs/rhs.m_sum;
	}
}