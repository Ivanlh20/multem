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

#ifndef KAHAN_SUM_H
	#define KAHAN_SUM_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "macros.cuh"

	namespace mt
	{
		/***************************************************************************************/
		/**************************** Kahan summation algorithm ********************************/
		/* https:// en.wikipedia.org/wiki/Kahan_summation_algorithmnnn */
		/***************************************************************************************/
		template <class T>
		class KS
		{
		public:
			using value_type = T;
		
			T m_sum;
			T m_error;

			/************************************* constructors ************************************/
			/* Default constructor */
			CGPU_EXEC
			KS(): m_sum(0), m_error(0) {}

			/* Copy constructor */
			CGPU_EXEC
			KS(const T& sum, T error=0): m_sum(sum), m_error(error) {}

			/* Copy constructor */
			CGPU_EXEC
			KS(const KS<T>& ks)
			{
				*this = ks;
			}

			/* Move constructor */
			CGPU_EXEC
			KS(KS<T>&& ks)
			{
				*this = std::move(ks);
			}

			/* converting copy constructor */
			template <class U>
			CGPU_EXEC
			KS(const KS<U>& ks)
			{
				*this = ks;
			}

			/******************************** assignment operators *********************************/
			// ! constructor from r
			CGPU_EXEC
			KS<T>& operator=(const T& r)
			{
				m_sum = r;
				m_error = 0;

				return *this;
			}

			/* copy assignment operator */
			CGPU_EXEC
			KS<T>& operator=(const KS<T>& ks)
			{
				m_sum = ks.m_sum;
				m_error = ks.m_error;

				return *this;
			}

			/* Move assignment operator */
			CGPU_EXEC
			KS<T>& operator=(KS<T>&& ks)
			{
				m_sum = std::move(ks.m_sum);
				m_error = std::move(ks.m_error);

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			KS<T>& operator=(const KS<U>& ks)
			{
				m_sum = T(ks.m_sum);
				m_error = T(ks.m_error);

				return *this;
			}

			// user define conversion
			CGPU_EXEC
			T operator()()
			{
				return m_sum;
			}
			

			// user define conversion
			CGPU_EXEC
			operator T()
			{
				return m_sum;
			}

			/******************* Compound assignment operators *******************/
			CGPU_EXEC
			KS<T>& operator+=(T r)
			{
				this->add(r);

				return *this;
			}
		
			CGPU_EXEC
			KS<T>& operator+=(const KS<T>& r)
			{
				this->add(r.m_sum);

				return *this;
			}
		
			CGPU_EXEC
			KS<T>& operator-=(T r)
			{
				this->add(-r);

				return *this;
			}	

			CGPU_EXEC
			KS<T>& operator-=(const KS<T>& r)
			{
				this->add(-r.m_sum);

				return *this;
			}

			CGPU_EXEC
			KS<T>& operator/=(T r)
			{
				*this = *this/r;

				return *this;
			}

			CGPU_EXEC
			KS<T>& operator/=(const KS<T>& r)
			{
				*this = *this/r;

				return *this;
			}

			CGPU_EXEC
			void add(T r)
			{
				r = r - m_error;
				T t = m_sum + r;			// m_sum is big, r small, so low-order digits of r are lost.
				m_error = (t-m_sum)-r;		// (t - m_sum) cancels the high-order part of r;subtracting r recovers negative (low part of r)
				m_sum = t;					// Algebraically, m_error should always be zero. Beware overly-aggressive optimizing compilers!

				// Next time around, the lost low part will be added to y in a fresh attempt
			}
		};

		/******************** Unitary arithmetic operators ********************/
		template <class T>
		CGPU_EXEC
		KS<T> operator-(const KS<T>& r)
		{
			return {-r.m_sum, r.m_error};
		}

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
#endif