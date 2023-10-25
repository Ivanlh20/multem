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

#include <stdint.h>
#include <algorithm>

#include "macros.h"

/* template definition */
template <class T> class dt_shape_st;


/* derived class */
using dt_shape = dt_shape_st<uint32_t>;

using dt_shape_64 = dt_shape_st<uint64_t>;


template <class T>
class dt_shape_st
{
public:
	/************************************* constructors ************************************/
	dt_shape_st();

	template <class U>
	dt_shape_st(const U& s0);

	template <class U, class V>
	dt_shape_st(const U& s0, const V& s1);

	template <class U, class V, class X>
	dt_shape_st(const U& s0, const V& s1, const X& s2);

	template <class U, class V, class X, class Y>
	dt_shape_st(const U& s0, const V& s1, const X& s2, const Y& s3);

	/* copy constructor */
	dt_shape_st(const dt_shape_st<T>& shape);

	/* converting constructor */
	template <class U>
	dt_shape_st(const dt_shape_st<U>& shape);

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	dt_shape_st<T>& operator=(const dt_shape_st<T>& shape);

	/* converting assignment operator */
	template <class U>
	dt_shape_st<T>& operator=(const dt_shape_st<U>& shape);

	T& operator[](const int32_t& iy);

	const T& operator[](const int32_t& iy) const;

	T* data();

	const T* data() const;

	T size();
		
	T shape_size();	

	void swap(const int32_t& ind_0, const int32_t& ind_e);

	T dim() const;

	T m_data[4];

private:
	const int32_t dim_max;
};

#include "../src/shape_t.inl"