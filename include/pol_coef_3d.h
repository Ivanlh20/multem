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

#include <type_traits>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
#include <utility>
#include <functional> 

#include "const_enum.h"
#include "type_traits_gen.h"
#include "math_mt.h"
#include "r_2d.h"
#include "r_3d.h"
#include "fcns_cgpu_gen.h"
#include "vctr_cpu.h"	
#include "vctr_gpu.h"

/***************************************************************************************/
/************************* linear and non-linear coefficients **************************/
/***************************************************************************************/
/* vector forward declaration */
namespace mt
{
#ifndef LNL_COEF_DEC
	#define LNL_COEF_DEC
	template <class ST, eDev Dev> class pLNL_Coef;

	template <class ST, eDev Dev> class LNL_Coef;
#endif
}

/* derived class */
namespace mt
{
	template <class T>
	using pLNL_Coef_cpu = pLNL_Coef<T, edev_cpu>;

	template <class T>
	using pLNL_Coef_gpu = pLNL_Coef<T, edev_gpu>;
}

namespace mt
{
	template <class T, eDev Dev>
	class pLNL_Coef
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ cl;		// linear coefficient
		T* __restrict__ cnl;	// non-linear coefficient
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pLNL_Coef();

		CGPU_EXEC
		pLNL_Coef(T* cl, T* cnl, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pLNL_Coef(const pLNL_Coef<T, Dev>& pcoef_lnl);

		/* constructor from LNL_Coef */
		explicit pLNL_Coef(const LNL_Coef<T, Dev>& coef_lnl);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pLNL_Coef<T, Dev>& operator=(const pLNL_Coef<T, Dev>& pcoef_lnl);

		/* Assignment operator: LNL_Coef -> pLNL_Coef */
		CPU_EXEC
		pLNL_Coef<T, Dev>& operator=(const LNL_Coef<T, Dev>& coef_lnl);

		size_type size() const;
	};
}

namespace mt
{
	namespace mt
	{
		template <class T>
		using LNL_Coef_cpu = LNL_Coef<T, edev_cpu>;

		template <class T>
		using LNL_Coef_gpu = LNL_Coef<T, edev_gpu>;
	}
}

namespace mt
{
	template <class T, eDev Dev>
	class LNL_Coef
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> cl;
		mutable Vctr<T, Dev> cnl;

		/************************************* constructors ************************************/
		LNL_Coef() {}

		LNL_Coef(const Vctr_cpu<T>& cl, const Vctr_cpu<T>& cnl);

		LNL_Coef(const dt_init_list_f64& cl, const dt_init_list_f64& cnl);

		LNL_Coef(const size_type& new_size);

		LNL_Coef(const size_type& new_size, const T& value);

		/* copy constructor */
		LNL_Coef(const LNL_Coef<T, Dev>& coef_lnl);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		LNL_Coef(const LNL_Coef<U, Dev_u>& coef_lnl);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		LNL_Coef<T, Dev>& operator=(const LNL_Coef<T, Dev>& coef_lnl);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		LNL_Coef<T, Dev>& operator=(const LNL_Coef<U, Dev_u>& coef_lnl);

		template <class U, eDev Dev_u> 
		void assign(const LNL_Coef<U, Dev_u>& coef_lnl);

		/**************** user define conversion operators *******************/
		pLNL_Coef<T, Dev> ptr() const;

		/* user define conversion for pointer */
		operator pLNL_Coef<T, Dev>() const;

		void fill(const T& val_l, const T& val_nl);

		size_type size() const;

		void clear();
		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};
}

/***************************************************************************************/
/*************************** linear polynomial coefficients ****************************/
/***************************************************************************************/
/* vector forward declaration */
namespace mt
{
#ifndef LNL_POLY_COEF_DEC
	#define LNL_POLY_COEF_DEC
	template <class ST, eDev Dev> class pPoly_Coef_1d;

	template <class ST, eDev Dev> class Poly_Coef_1d;
#endif
}

/* derived class */
namespace mt
{
	namespace mt
	{
		template <class T>
		using pPoly_Coef_1d_cpu = pPoly_Coef_1d<T, edev_cpu>;

		template <class T>
		using pPoly_Coef_1d_gpu = pPoly_Coef_1d<T, edev_gpu>;
	}
}

namespace mt
{	
	template <class T, eDev Dev> class Poly_Coef_1d;

	template <class T, eDev Dev>
	class pPoly_Coef_1d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ c0;
		T* __restrict__ c1;
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pPoly_Coef_1d();

		CGPU_EXEC
		pPoly_Coef_1d(T* c0, T* c1, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pPoly_Coef_1d(const pPoly_Coef_1d<T, Dev>& pcoef_poly1);

		/* constructor from Poly_Coef_1d */
		explicit pPoly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pPoly_Coef_1d<T, Dev>& operator=(const pPoly_Coef_1d<T, Dev>& pcoef_poly1);

		/* Assignment operator: Poly_Coef_1d -> pPoly_Coef_1d */
		CPU_EXEC
		pPoly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<T, Dev>& coef_poly1);

		size_type size() const;
	};

	/* derived class */
	namespace mt
	{
		namespace mt
		{
			template <class T>
			using Poly_Coef_1d_cpu = Poly_Coef_1d<T, edev_cpu>;

			template <class T>
			using Poly_Coef_1d_gpu = Poly_Coef_1d<T, edev_gpu>;
		}
	}

	template <class T, eDev Dev>
	class Poly_Coef_1d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> c0;
		mutable Vctr<T, Dev> c1;

		/************************************* constructors ************************************/
		Poly_Coef_1d() {}

		Poly_Coef_1d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1);

		Poly_Coef_1d(const dt_init_list_f64& c0, const dt_init_list_f64& c1);

		Poly_Coef_1d(const size_type& new_size);

		Poly_Coef_1d(const size_type& new_size, const T& value);

		/* copy constructor */
		Poly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		Poly_Coef_1d(const Poly_Coef_1d<U, Dev_u>& coef_poly1);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Poly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<T, Dev>& coef_poly1);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		Poly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<U, Dev_u>& coef_poly1);

		template <class U, eDev Dev_u> 
		void assign(const Poly_Coef_1d<U, Dev_u>& coef_poly1);

		/**************** user define conversion operators *******************/
		pPoly_Coef_1d<T, Dev> ptr() const;

		/* user define conversion for pointer Vctr */
		operator pPoly_Coef_1d<T, Dev>() const;

		void fill(const T& val_c0, const T& val_c1);

		size_type size() const;

		void clear();

		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};
}

/***************************************************************************************/
/************************* quadratic polynomial coefficients ***************************/
/***************************************************************************************/
namespace mt
{
	template <class T, eDev Dev> class Poly_Coef_2d;

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
		pPoly_Coef_2d(): c0(nullptr), c1(nullptr), c2(nullptr), m_size(0) {}

		CGPU_EXEC
		pPoly_Coef_2d(T* c0, T* c1, T* c2, const size_type& size): c0(c0), c1(c1), c2(c2), m_size(size) {}

		/* copy constructor */
		CGPU_EXEC
		pPoly_Coef_2d(const pPoly_Coef_2d<T, Dev>& pcoef_poly2)
		{
			*this = pcoef_poly2;
		}

		/* constructor from Poly_Coef_2d */
		explicit pPoly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2)
		{
			*this = coef_poly2;
		}

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pPoly_Coef_2d<T, Dev>& operator=(const pPoly_Coef_2d<T, Dev>& pcoef_poly2)
		{
			if (this != &pcoef_poly2)
			{
				c0 = pcoef_poly2.c0;
				c1 = pcoef_poly2.c1;
				c2 = pcoef_poly2.c2;
				m_size = pcoef_poly2.m_size;
			}

			return *this;
		}

		/* Assignment operator: Poly_Coef_2d -> pPoly_Coef_2d */
		CPU_EXEC
		pPoly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<T, Dev>& coef_poly2)
		{
			c0 = coef_poly2.c0.data();
			c1 = coef_poly2.c1.data();
			c2 = coef_poly2.c2.data();
			m_size = size_type(coef_poly2.size());

			return *this;
		}

		size_type size() const
		{
			return m_size;
		}
	};

	template <class T>
	using pPoly_Coef_2d_cpu = pPoly_Coef_2d<T, edev_cpu>;

	template <class T>
	using pPoly_Coef_2d_gpu = pPoly_Coef_2d<T, edev_gpu>;

	/***************************************************************************************/
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

	template <class T>
	using Poly_Coef_2d_cpu = Poly_Coef_2d<T, edev_cpu>;

	template <class T>
	using Poly_Coef_2d_gpu = Poly_Coef_2d<T, edev_gpu>;
}

/***************************************************************************************/
/************************** cubic interpolation coefficients ***************************/
/***************************************************************************************/
namespace mt
{
	template <class T, eDev Dev> class Poly_Coef_3d;

	template <class T, eDev Dev>
	class pPoly_Coef_3d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		T* __restrict__ c0;
		T* __restrict__ c1;
		T* __restrict__ c2;
		T* __restrict__ c3;
		size_type m_size;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pPoly_Coef_3d();

		CGPU_EXEC
		pPoly_Coef_3d(T* c0, T* c1, T* c2, T* c3, const size_type& size);

		/* copy constructor */
		CGPU_EXEC
		pPoly_Coef_3d(const pPoly_Coef_3d<T, Dev>& pcoef_poly3);

		/* constructor from Poly_Coef_3d */
		explicit pPoly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pPoly_Coef_3d<T, Dev>& operator=(const pPoly_Coef_3d<T, Dev>& pcoef_poly3);

		/* Assignment operator: Poly_Coef_3d -> pPoly_Coef_3d */
		CPU_EXEC
		pPoly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<T, Dev>& coef_poly3);

		size_type size() const;
	};

	template <class T>
	using pPoly_Coef_3d_cpu = pPoly_Coef_3d<T, edev_cpu>;

	template <class T>
	using pPoly_Coef_3d_gpu = pPoly_Coef_3d<T, edev_gpu>;

	/***************************************************************************************/
	template <class T, eDev Dev>
	class Poly_Coef_3d
	{
	public:
		using value_type = T;
		using size_type = dt_int32;
		static const eDev device = Dev;

		mutable Vctr<T, Dev> c0;
		mutable Vctr<T, Dev> c1;
		mutable Vctr<T, Dev> c2;
		mutable Vctr<T, Dev> c3;

		/************************************* constructors ************************************/
		Poly_Coef_3d() {}

		Poly_Coef_3d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2, const Vctr_cpu<T>& c3);

		Poly_Coef_3d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2, const dt_init_list_f64& c3);

		Poly_Coef_3d(const size_type& new_size);

		Poly_Coef_3d(const size_type& new_size, const T& value);

		/* copy constructor */
		Poly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3);

		/* converting constructor */
		template <class U, eDev Dev_u> 
		Poly_Coef_3d(const Poly_Coef_3d<U, Dev_u>& coef_poly3);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Poly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<T, Dev>& coef_poly3);

		/* converting assignment operator */
		template <class U, eDev Dev_u> 
		Poly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<U, Dev_u>& coef_poly3);

		template <class U, eDev Dev_u> 
		void assign(const Poly_Coef_3d<U, Dev_u>& coef_poly3);

		/**************** user define conversion operators *******************/
		pPoly_Coef_3d<T, Dev> ptr() const;

		/* user define conversion for pointer Vctr */
		operator pPoly_Coef_3d<T, Dev>() const;

		void fill(const T& val_c0, const T& val_c1, const T& val_c2, const T& val_c3);

		size_type size() const;

		void clear();

		void reserve(const size_type& new_size);

		void resize(const size_type& new_size);

		void resize(const size_type& new_size, const T& value);

		void shrink_to_fit();
	};

	template <class T>
	using Poly_Coef_3d_cpu = Poly_Coef_3d<T, edev_cpu>;

	template <class T>
	using Poly_Coef_3d_gpu = Poly_Coef_3d<T, edev_gpu>;
}