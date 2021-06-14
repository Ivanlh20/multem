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

#ifndef SHAPE_H
	#define SHAPE_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <stdint.h>
	#include <algorithm>

	#include "macros.cuh"

	template <class T>
	class dt_shape_st
	{
	public:
		/************************************* constructors ************************************/
		dt_shape_st(): dim_max(4) 
		{
			for (int32_t ik=0; ik<dim_max; ik++)
			{
				m_data[ik] = T(1);
			}
		}

		template <class U>
		dt_shape_st(const U& s0): dt_shape_st()
		{
			m_data[0] = T(s0);
		}

		template <class U, class V>
		dt_shape_st(const U& s0, const V& s1): dt_shape_st()
		{
			m_data[0] = T(s0);
			m_data[1] = T(s1);
		}

		template <class U, class V, class X>
		dt_shape_st(const U& s0, const V& s1, const X& s2): dt_shape_st()
		{
			m_data[0] = T(s0);
			m_data[1] = T(s1);
			m_data[2] = T(s2);
		}

		template <class U, class V, class X, class Y>
		dt_shape_st(const U& s0, const V& s1, const X& s2, const Y& s3): dt_shape_st()
		{
			m_data[0] = T(s0);
			m_data[1] = T(s1);
			m_data[2] = T(s2);
			m_data[3] = T(s3);
		}

		/* copy constructor */
		dt_shape_st(const dt_shape_st<T>& shape): dt_shape_st()
		{
			*this = shape;
		}

		/* converting constructor */
		template <class U>
		dt_shape_st(const dt_shape_st<U>& shape): dt_shape_st()
		{
			*this = shape;
		}

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		dt_shape_st<T>& operator=(const dt_shape_st<T>& shape)
		{
			for (int32_t ik=0; ik<dim_max; ik++)
			{
				m_data[ik] = shape.m_data[ik];
			}

			return *this;
		}

		/* converting assignment operator */
		template <class U>
		dt_shape_st<T>& operator=(const dt_shape_st<U>& shape)
		{
			for (int32_t ik=0; ik<dim_max; ik++)
			{
				m_data[ik] = T(shape.m_data[ik]);
			}

			return *this;
		}

		T& operator[](const int32_t& iy)
		{ 
			return m_data[iy];
		}

		const T& operator[](const int32_t& iy) const 
		{ 
			return m_data[iy];
		}

		T* data()
		{
			return m_data;
		}

		const T* data() const
		{
			return m_data;
		}

		T size()
		{
			return T(dim_max);
		}
		
		T shape_size()
		{
			T ss = 1;
			for (int32_t ik=0; ik<dim_max; ik++) 
			{
				ss = ss*std::max(T(1), m_data[ik]);
			}

			return ss;
		}	

		void swap(const int32_t& ind_0, const int32_t& ind_e)
		{
			std::swap(m_data[ind_0], m_data[ind_e]);
		}

		T dim() const
		{
			int32_t ik_0 = dim_max - 1;
			int32_t ic = 0;
			for (int32_t ik=ik_0; ik>1; ik--) 
			{
				if (m_data[ik]<2)
				{
					ic++;
				}
				else
					break;
			}

			return T(dim_max-ic);
		}

		T m_data[4];

	private:
		const int32_t dim_max;
	};

	using dt_shape = dt_shape_st<uint32_t>;

	using dt_shape_64 = dt_shape_st<uint64_t>;
#endif