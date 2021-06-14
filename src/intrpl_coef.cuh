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

#ifndef INTRPL_COEF_H
	#define INTRPL_COEF_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include <type_traits>
	#include <vector>
	#include <algorithm>
	#include <thread>
	#include <mutex>
	#include <utility>
	#include <functional> 

	#include "const_enum.cuh"
	#include "type_traits_gen.cuh"
	#include "math.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "cgpu_vctr.cuh"	
	#include "grid.cuh"	

	/***************************************************************************************/
	/************************* linear and non-linear coefficients **************************/
	/***************************************************************************************/
	namespace mt
	{		
		template <class T, eDev Dev> class LNL_Coef;

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
			pLNL_Coef(): cl(nullptr), cnl(nullptr), m_size(0) {}

			CGPU_EXEC
			pLNL_Coef(T* cl, T* cnl, const size_type& size): cl(cl), cnl(cnl), m_size(size) {}

			/* copy constructor */
			CGPU_EXEC
			pLNL_Coef(const pLNL_Coef<T, Dev>& pcoef_lnl)
			{
				*this = pcoef_lnl;
			}

			// ! constructor from LNL_Coef
			explicit pLNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
			{
				*this = coef_lnl;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			pLNL_Coef<T, Dev>& operator=(const pLNL_Coef<T, Dev>& pcoef_lnl)
			{
				if (this != &pcoef_lnl)
				{
					cl = pcoef_lnl.cl;
					cnl = pcoef_lnl.cnl;
					m_size = pcoef_lnl.m_size;
				}

				return *this;
			}

			// ! Assignment operator: LNL_Coef -> pLNL_Coef
			CPU_EXEC
			pLNL_Coef<T, Dev>& operator=(const LNL_Coef<T, Dev>& coef_lnl)
			{
				cl = coef_lnl.cl.data();
				cnl = coef_lnl.cnl.data();
				m_size = size_type(coef_lnl.size());

				return *this;
			}

			size_type size() const
			{
				return m_size;
			}
		};

		template <class T>
		using pLNL_Coef_cpu = pLNL_Coef<T, edev_cpu>;

		template <class T>
		using pLNL_Coef_gpu = pLNL_Coef<T, edev_gpu>;

		/***************************************************************************************/
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

			LNL_Coef(const Vctr_cpu<T>& cl, const Vctr_cpu<T>& cnl): 
				cl(cl), cnl(cnl){}

			LNL_Coef(const dt_init_list_f64& cl, const dt_init_list_f64& cnl): 
				cl(cl), cnl(cnl) {}

			LNL_Coef(const size_type& new_size)
			{
				resize(new_size);
			}

			LNL_Coef(const size_type& new_size, const T& value)
			{
				resize(new_size, value);
			}

			/* copy constructor */
			LNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
			{
				*this = coef_lnl;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			LNL_Coef(const LNL_Coef<U, Dev_u>& coef_lnl)
			{
				*this = coef_lnl;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			LNL_Coef<T, Dev>& operator=(const LNL_Coef<T, Dev>& coef_lnl)
			{
				this->assign(coef_lnl);
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			LNL_Coef<T, Dev>& operator=(const LNL_Coef<U, Dev_u>& coef_lnl)
			{
				this->assign(coef_lnl);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const LNL_Coef<U, Dev_u>& coef_lnl)
			{ 
				cl.assign(coef_lnl.cl);
				cnl.assign(coef_lnl.cnl);
			}

			/**************** user define conversion operators *******************/
			pLNL_Coef<T, Dev> ptr() const
			{
				return pLNL_Coef<T, Dev>(*this);
			}

			// ! user define conversion for pointer
			operator pLNL_Coef<T, Dev>() const
			{
				return pLNL_Coef<T, Dev>(*this);
			}

			void fill(const T& val_l, const T& val_nl)
			{
				cl.fill(val_l);
				cnl.fill(val_nl);
			}

			size_type size() const
			{
				return size_type(cl.size());
			}

			void clear()
			{
				cl.clear();
				cnl.clear();
			}

			void reserve(const size_type& new_size)
			{
				cl.reserve(new_size);
				cnl.reserve(new_size);
			}

			void resize(const size_type& new_size)
			{
				cl.resize(new_size);
				cnl.resize(new_size);
			}

			void resize(const size_type& new_size, const T& value)
			{
				cl.resize(new_size, value);
				cnl.resize(new_size, value);
			}

			void shrink_to_fit()
			{
				cl.shrink_to_fit();
				cnl.shrink_to_fit();
			}
		};

		template <class T>
		using LNL_Coef_cpu = LNL_Coef<T, edev_cpu>;

		template <class T>
		using LNL_Coef_gpu = LNL_Coef<T, edev_gpu>;
	}

	/***************************************************************************************/
	/*************************** linear polynomial coefficients ****************************/
	/***************************************************************************************/
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
			pPoly_Coef_1d(): c0(nullptr), c1(nullptr), m_size(0) {}

			CGPU_EXEC
			pPoly_Coef_1d(T* c0, T* c1, const size_type& size): c0(c0), c1(c1), m_size(size) {}

			/* copy constructor */
			CGPU_EXEC
			pPoly_Coef_1d(const pPoly_Coef_1d<T, Dev>& pcoef_poly1)
			{
				*this = pcoef_poly1;
			}

			// ! constructor from Poly_Coef_1d
			explicit pPoly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1)
			{
				*this = coef_poly1;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			pPoly_Coef_1d<T, Dev>& operator=(const pPoly_Coef_1d<T, Dev>& pcoef_poly1)
			{
				if (this != &pcoef_poly1)
				{
					c0 = pcoef_poly1.c0;
					c1 = pcoef_poly1.c1;
					m_size = pcoef_poly1.m_size;
				}

				return *this;
			}

			// ! Assignment operator: Poly_Coef_1d -> pPoly_Coef_1d
			CPU_EXEC
			pPoly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<T, Dev>& coef_poly1)
			{
				c0 = coef_poly1.c0.data();
				c1 = coef_poly1.c1.data();
				m_size = size_type(coef_poly1.size());

				return *this;
			}

			size_type size() const
			{
				return m_size;
			}
		};

		template <class T>
		using pPoly_Coef_1d_cpu = pPoly_Coef_1d<T, edev_cpu>;

		template <class T>
		using pPoly_Coef_1d_gpu = pPoly_Coef_1d<T, edev_gpu>;

		/***************************************************************************************/
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

			Poly_Coef_1d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1): 
				c0(c0), c1(c1){}

			Poly_Coef_1d(const dt_init_list_f64& c0, const dt_init_list_f64& c1): 
				c0(c0), c1(c1) {}

			Poly_Coef_1d(const size_type& new_size)
			{
				resize(new_size);
			}

			Poly_Coef_1d(const size_type& new_size, const T& value)
			{
				resize(new_size, value);
			}

			/* copy constructor */
			Poly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1)
			{
				*this = coef_poly1;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			Poly_Coef_1d(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
			{
				*this = coef_poly1;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Poly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<T, Dev>& coef_poly1)
			{
				this->assign(coef_poly1);
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			Poly_Coef_1d<T, Dev>& operator=(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
			{
				this->assign(coef_poly1);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
			{ 
				c0.assign(coef_poly1.c0);
				c1.assign(coef_poly1.c1);
			}

			/**************** user define conversion operators *******************/
			pPoly_Coef_1d<T, Dev> ptr() const
			{
				return pPoly_Coef_1d<T, Dev>(*this);
			}

			// ! user define conversion for pointer Vctr
			operator pPoly_Coef_1d<T, Dev>() const
			{
				return pPoly_Coef_1d<T, Dev>(*this);
			}

			void fill(const T& val_c0, const T& val_c1)
			{
				c0.fill(val_c0);
				c1.fill(val_c1);
			}

			size_type size() const
			{
				return c0.size();
			}

			void clear()
			{
				c0.clear();
				c1.clear();
			}

			void reserve(const size_type& new_size)
			{
				c0.reserve(new_size);
				c1.reserve(new_size);
			}

			void resize(const size_type& new_size)
			{
				c0.resize(new_size);
				c1.resize(new_size);
			}

			void resize(const size_type& new_size, const T& value)
			{
				c0.resize(new_size, value);
				c1.resize(new_size, value);
			}

			void shrink_to_fit()
			{
				c0.shrink_to_fit();
				c1.shrink_to_fit();
			}
		};

		template <class T>
		using Poly_Coef_1d_cpu = Poly_Coef_1d<T, edev_cpu>;

		template <class T>
		using Poly_Coef_1d_gpu = Poly_Coef_1d<T, edev_gpu>;
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

			// ! constructor from Poly_Coef_2d
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

			// ! Assignment operator: Poly_Coef_2d -> pPoly_Coef_2d
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

			Poly_Coef_2d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2): 
				c0(c0), c1(c1), c2(c2){}

			Poly_Coef_2d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2): 
				c0(c0), c1(c1), c2(c2) {}

			Poly_Coef_2d(const size_type& new_size)
			{
				resize(new_size);
			}

			Poly_Coef_2d(const size_type& new_size, const T& value)
			{
				resize(new_size, value);
			}

			/* copy constructor */
			Poly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2)
			{
				*this = coef_poly2;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			Poly_Coef_2d(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
			{
				*this = coef_poly2;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Poly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<T, Dev>& coef_poly2)
			{
				this->assign(coef_poly2);
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			Poly_Coef_2d<T, Dev>& operator=(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
			{
				this->assign(coef_poly2);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
			{ 
				c0.assign(coef_poly2.c0);
				c1.assign(coef_poly2.c1);
				c2.assign(coef_poly2.c2);
			}

			/**************** user define conversion operators *******************/
			pPoly_Coef_2d<T, Dev> ptr() const
			{
				return pPoly_Coef_2d<T, Dev>(*this);
			}

			// ! user define conversion for pointer Vctr
			operator pPoly_Coef_2d<T, Dev>() const
			{
				return pPoly_Coef_2d<T, Dev>(*this);
			}

			void fill(const T& val_c0, const T& val_c1, const T& val_c2)
			{
				c0.fill(val_c0);
				c1.fill(val_c1);
				c2.fill(val_c2);
			}

			size_type size() const
			{
				return c0.size();
			}

			void clear()
			{
				c0.clear();
				c1.clear();
				c2.clear();
			}

			void reserve(const size_type& new_size)
			{
				c0.reserve(new_size);
				c1.reserve(new_size);
				c2.reserve(new_size);
			}

			void resize(const size_type& new_size)
			{
				c0.resize(new_size);
				c1.resize(new_size);
				c2.resize(new_size);
			}

			void resize(const size_type& new_size, const T& value)
			{
				c0.resize(new_size, value);
				c1.resize(new_size, value);
				c2.resize(new_size, value);
			}

			void shrink_to_fit()
			{
				c0.shrink_to_fit();
				c1.shrink_to_fit();
				c2.shrink_to_fit();
			}
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
			pPoly_Coef_3d(): c0(nullptr), c1(nullptr), c2(nullptr), c3(nullptr), m_size(0) {}

			CGPU_EXEC
			pPoly_Coef_3d(T* c0, T* c1, T* c2, T* c3, const size_type& size): c0(c0), c1(c1), c2(c2), c3(c3), m_size(size) {}

			/* copy constructor */
			CGPU_EXEC
			pPoly_Coef_3d(const pPoly_Coef_3d<T, Dev>& pcoef_poly3)
			{
				*this = pcoef_poly3;
			}

			// ! constructor from Poly_Coef_3d
			explicit pPoly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3)
			{
				*this = coef_poly3;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			pPoly_Coef_3d<T, Dev>& operator=(const pPoly_Coef_3d<T, Dev>& pcoef_poly3)
			{
				if (this != &pcoef_poly3)
				{
					c0 = pcoef_poly3.c0;
					c1 = pcoef_poly3.c1;
					c2 = pcoef_poly3.c2;
					c3 = pcoef_poly3.c3;
					m_size = pcoef_poly3.m_size;
				}

				return *this;
			}

			// ! Assignment operator: Poly_Coef_3d -> pPoly_Coef_3d
			CPU_EXEC
			pPoly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<T, Dev>& coef_poly3)
			{
				c0 = coef_poly3.c0.data();
				c1 = coef_poly3.c1.data();
				c2 = coef_poly3.c2.data();
				c3 = coef_poly3.c3.data();
				m_size = size_type(coef_poly3.size());

				return *this;
			}

			size_type size() const
			{
				return m_size;
			}
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

			Poly_Coef_3d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2, const Vctr_cpu<T>& c3): 
				c0(c0), c1(c1), c2(c2), c3(c3){}

			Poly_Coef_3d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2, const dt_init_list_f64& c3): 
				c0(c0), c1(c1), c2(c2), c3(c3) {}

			Poly_Coef_3d(const size_type& new_size)
			{
				resize(new_size);
			}

			Poly_Coef_3d(const size_type& new_size, const T& value)
			{
				resize(new_size, value);
			}

			/* copy constructor */
			Poly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3)
			{
				*this = coef_poly3;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			Poly_Coef_3d(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
			{
				*this = coef_poly3;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Poly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<T, Dev>& coef_poly3)
			{
				this->assign(coef_poly3);
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			Poly_Coef_3d<T, Dev>& operator=(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
			{
				this->assign(coef_poly3);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
			{ 
				c0.assign(coef_poly3.c0);
				c1.assign(coef_poly3.c1);
				c2.assign(coef_poly3.c2);
				c3.assign(coef_poly3.c3);
			}

			/**************** user define conversion operators *******************/
			pPoly_Coef_3d<T, Dev> ptr() const
			{
				return pPoly_Coef_3d<T, Dev>(*this);
			}

			// ! user define conversion for pointer Vctr
			operator pPoly_Coef_3d<T, Dev>() const
			{
				return pPoly_Coef_3d<T, Dev>(*this);
			}

			void fill(const T& val_c0, const T& val_c1, const T& val_c2, const T& val_c3)
			{
				c0.fill(val_c0);
				c1.fill(val_c1);
				c2.fill(val_c2);
				c3.fill(val_c3);
			}

			size_type size() const
			{
				return c0.size();
			}

			void clear()
			{
				c0.clear();
				c1.clear();
				c2.clear();
				c3.clear();
			}

			void reserve(const size_type& new_size)
			{
				c0.reserve(new_size);
				c1.reserve(new_size);
				c2.reserve(new_size);
				c3.reserve(new_size);
			}

			void resize(const size_type& new_size)
			{
				c0.resize(new_size);
				c1.resize(new_size);
				c2.resize(new_size);
				c3.resize(new_size);
			}

			void resize(const size_type& new_size, const T& value)
			{
				c0.resize(new_size, value);
				c1.resize(new_size, value);
				c2.resize(new_size, value);
				c3.resize(new_size, value);
			}

			void shrink_to_fit()
			{
				c0.shrink_to_fit();
				c1.shrink_to_fit();
				c2.shrink_to_fit();
				c3.shrink_to_fit();
			}
		};

		template <class T>
		using Poly_Coef_3d_cpu = Poly_Coef_3d<T, edev_cpu>;

		template <class T>
		using Poly_Coef_3d_gpu = Poly_Coef_3d<T, edev_gpu>;
	}
#endif 