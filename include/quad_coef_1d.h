/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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
#ifndef QUAD_COEF_1D
	#define QUAD_COEF_1D
	template <class T, eDev Dev> class Quad_Coef_1d;

	template <class T, eDev Dev> class pQuad_Coef_1d;
#endif
}

/* class one dimensional quadrature coefficients */
namespace mt
{
	template <class T>
	using Quad_Coef_1d_cpu = Quad_Coef_1d<T, edev_cpu>;

	template <class T>
	using Quad_Coef_1d_gpu = Quad_Coef_1d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class Quad_Coef_1d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> x;
		mutable Vctr<T, Dev> w;

		/************************************* constructors ************************************/
		Quad_Coef_1d() = default;

		Quad_Coef_1d(const Vctr_cpu<T>& x, const Vctr_cpu<T>& w);

		Quad_Coef_1d(const dt_init_list_f64& x, const dt_init_list_f64& w);

		Quad_Coef_1d(const size_type& new_size);

		Quad_Coef_1d(const size_type& new_size, const T& value);

		/* copy constructor */
		Quad_Coef_1d(const Quad_Coef_1d<T, Dev>& coef_quad_1d);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		Quad_Coef_1d(const Quad_Coef_1d<U, Dev_u>& coef_quad_1d);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Quad_Coef_1d<T, Dev>& operator=(const Quad_Coef_1d<T, Dev>& coef_quad_1d);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		Quad_Coef_1d<T, Dev>& operator=(const Quad_Coef_1d<U, Dev_u>& coef_quad_1d);

		template <class U, eDev Dev_u> 
		void assign(const Quad_Coef_1d<U, Dev_u>& coef_quad_1d);

		/**************** user define conversion operators *******************/
		pQuad_Coef_1d<T, Dev> ptr() const;

		/* user define conversion for pointer Vctr */
		operator pQuad_Coef_1d<T, Dev>() const;

		void fill(const T& val_x, const T& val_w);

		size_type size() const;

		void clear();

		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};
}

/* pointer class one dimensional quadrature coefficients */
namespace mt
{
	template <class T>
	using pQuad_Coef_1d_cpu = pQuad_Coef_1d<T, edev_cpu>;

	template <class T>
	using pQuad_Coef_1d_gpu = pQuad_Coef_1d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class pQuad_Coef_1d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ x;
		T* __restrict__ w;
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pQuad_Coef_1d();

		CGPU_EXEC
		pQuad_Coef_1d(T* x, T* w, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pQuad_Coef_1d(const pQuad_Coef_1d<T, Dev>& pcoef_quad_1d);

		/* constructor from Quad_Coef_1d */
		explicit pQuad_Coef_1d(const Quad_Coef_1d<T, Dev>& coef_quad_1d);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pQuad_Coef_1d<T, Dev>& operator=(const pQuad_Coef_1d<T, Dev>& pcoef_quad_1d);

		/* Assignment operator: Quad_Coef_1d -> pQuad_Coef_1d */
		CPU_EXEC
		pQuad_Coef_1d<T, Dev>& operator=(const Quad_Coef_1d<T, Dev>& coef_quad_1d);

		size_type size() const;
	};
}

#include "../src/quad_coef_1d.inl"