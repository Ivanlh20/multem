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

/* template definition */
namespace mt
{
#ifndef LNL_COEF_DEC
	#define LNL_COEF_DEC
	template <class ST, eDev Dev> class LNL_Coef;
	
	template <class ST, eDev Dev> class pLNL_Coef;
#endif
}

// class linear and non-linear coefficients
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

// pointer class linear and non-linear coefficients
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

#include "../src/lnl_coef.inl"