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
#include "vctr_cpu.h"
#include "vctr_gpu.h"

/* template definition */
namespace mt
{
#ifndef POLY_COEF_2D_DEC
	#define POLY_COEF_2D_DEC
	template <class ST, eDev Dev> class Poly_Coef_2d;

	template <class ST, eDev Dev> class pPoly_Coef_2d;
#endif
}

// class quadratic polynomial coefficients
namespace mt
{
	template <class T>
	using Poly_Coef_2d_cpu = Poly_Coef_2d<T, edev_cpu>;

	template <class T>
	using Poly_Coef_2d_gpu = Poly_Coef_2d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class Poly_Coef_2d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> c0;
		mutable Vctr<T, Dev> c1;
		mutable Vctr<T, Dev> c2;

		/************************************* constructors ************************************/
		Poly_Coef_2d() {}

		Poly_Coef_2d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2);

		Poly_Coef_2d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2);

		Poly_Coef_2d(const size_type& new_size);

		Poly_Coef_2d(const size_type& new_size, const T& value);

		/* copy constructor */
		Poly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		Poly_Coef_2d(const Poly_Coef_2d<U, Dev_u>& coef_poly2);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Poly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<T, Dev>& coef_poly2);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		Poly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<U, Dev_u>& coef_poly2);

		template <class U, eDev Dev_u> 
		void assign(const Poly_Coef_2d<U, Dev_u>& coef_poly2);

		/**************** user define conversion operators *******************/
		pPoly_Coef_2d<T, Dev> ptr() const;

		/* user define conversion for pointer Vctr */
		operator pPoly_Coef_2d<T, Dev>() const;

		void fill(const T& val_c0, const T& val_c1, const T& val_c2);

		size_type size() const;

		void clear();

		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};
}


// pointer class quadratic polynomial coefficients
namespace mt
{
	template <class T>
	using pPoly_Coef_2d_cpu = pPoly_Coef_2d<T, edev_cpu>;

	template <class T>
	using pPoly_Coef_2d_gpu = pPoly_Coef_2d<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class pPoly_Coef_2d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ c0;
		T* __restrict__ c1;
		T* __restrict__ c2;
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pPoly_Coef_2d();

		CGPU_EXEC
		pPoly_Coef_2d(T* c0, T* c1, T* c2, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pPoly_Coef_2d(const pPoly_Coef_2d<T, Dev>& pcoef_poly2);

		/* constructor from Poly_Coef_2d */
		explicit pPoly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pPoly_Coef_2d<T, Dev>& operator=(const pPoly_Coef_2d<T, Dev>& pcoef_poly2);

		/* Assignment operator: Poly_Coef_2d -> pPoly_Coef_2d */
		CPU_EXEC
		pPoly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<T, Dev>& coef_poly2);

		size_type size() const;
	};
}

#include "../src/poly_coef_2d.inl"