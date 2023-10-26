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
#include "vctr_cpu.h"
#include "vctr_gpu.h"

/* template definition */
namespace mt
{
#ifndef QUAD_COEF_2D
	#define QUAD_COEF_2D
	template <class T, eDev Dev> class Quad_Coef_2d;
	
	template <class T, eDev Dev> class pQuad_Coef_2d;
#endif
}

/* class two dimensional quadrature coefficients */
namespace mt
{
	template <class T>
	using Quad_Coef_2d_cpu = Quad_Coef_2d<T, edev_cpu>;

	template <class T>
	using Quad_Coef_2d_gpu = Quad_Coef_2d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class Quad_Coef_2d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> x;
		mutable Vctr<T, Dev> y;
		mutable Vctr<T, Dev> w;

		/************************************* constructors ************************************/
		Quad_Coef_2d() = default;

		Quad_Coef_2d(const Vctr_cpu<T>& x, const Vctr_cpu<T>& y, const Vctr_cpu<T>& w);

		Quad_Coef_2d(const dt_init_list_f64& x, const dt_init_list_f64& y, const dt_init_list_f64& w);

		Quad_Coef_2d(const size_type& new_size);

		Quad_Coef_2d(const size_type& new_size, const T& value);

		/* copy constructor */
		Quad_Coef_2d(const Quad_Coef_2d<T, Dev>& coef_quad_2d);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		Quad_Coef_2d(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Quad_Coef_2d<T, Dev>& operator=(const Quad_Coef_2d<T, Dev>& coef_quad_2d);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		Quad_Coef_2d<T, Dev>& operator=(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d);

		template <class U, eDev Dev_u> 
		void assign(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d);

		/**************** user define conversion operators *******************/
		pQuad_Coef_2d<T, Dev> ptr() const;

		/* user define conversion for pointer Vctr */
		operator pQuad_Coef_2d<T, Dev>() const;

		void fill(const T& val_x, const T& val_y, const T& val_w);

		size_type size() const;

		void clear();

		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};
}

/* pointer class two dimensional quadrature coefficients */
namespace mt
{
	template <class T>
	using pQuad_Coef_2d_cpu = pQuad_Coef_2d<T, edev_cpu>;

	template <class T>
	using pQuad_Coef_2d_gpu = pQuad_Coef_2d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class pQuad_Coef_2d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ x;
		T* __restrict__ y;
		T* __restrict__ w;
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pQuad_Coef_2d();

		CGPU_EXEC
		pQuad_Coef_2d(T* x, T* y, T* w, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pQuad_Coef_2d(const pQuad_Coef_2d<T, Dev>& pcoef_quad_2d);

		/* constructor from Quad_Coef_2d */
		explicit pQuad_Coef_2d(const Quad_Coef_2d<T, Dev>& coef_quad_2d);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pQuad_Coef_2d<T, Dev>& operator=(const pQuad_Coef_2d<T, Dev>& pcoef_quad_2d);

		/* Assignment operator: Quad_Coef_2d -> pQuad_Coef_2d */
		CPU_EXEC
		pQuad_Coef_2d<T, Dev>& operator=(Quad_Coef_2d<T, Dev>& coef_quad_2d);

		size_type size() const;
	};
}

#include "../src/quad_coef_2d.inl"