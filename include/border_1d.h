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

#include "math_mt.h"
#include "const_enum.h"
#include "fcns_cgpu_gen.h"
#include "vctr_cpu.h"

/* template definition */
namespace mt
{
	template <class T, eDim Dim> class Border_Rect_xd;

	template <eDim Dim>
	using iBorder_Rect_xd = Border_Rect_xd<dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class T>
	using Border_Rect_1d = Border_Rect_xd<T, edim_1>;

	using iBorder_Rect_1d = Border_Rect_xd<dt_int32, edim_1>;

	using iBorder_Rect_1d_64 = Border_Rect_xd<dt_int64, edim_1>;
}	

/***************************************************************************************/
/******************************** rectangular border ***********************************/
/***************************************************************************************/
/* template specialization 1d */
namespace mt
{
	template <class T>
	class Border_Rect_xd<T, edim_1>
	{
	public:
		using value_type = T;
		using size_type = dt_int32;

		T bx_0;		// initial x position
		T bx_e;		// final x position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Border_Rect_xd();

		Border_Rect_xd(const T& bx_0, const T& bx_e);

		template <class U>
		Border_Rect_xd(const dt_init_list<U>& list);

		template <class U>
		Border_Rect_xd(const Vctr_cpu<U>& vctr);

		/* copy constructor */
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<T, edim_1>& border);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<U, edim_1>& border);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Border_Rect_xd<T, edim_1>& operator=(const Border_Rect_xd<T, edim_1>& border);		
			
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd<T, edim_1>& operator=(const Border_Rect_xd<U, edim_1>& border);

		template <class U> 
		CGPU_EXEC
		void assign(const Border_Rect_xd<U, edim_1>& border);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& bx_0, const U& bx_e);

		template <class U> 
		void set_in_data(const dt_init_list<U>& list);

		template <class U>
		void set_in_data(const Vctr_cpu<U>& vctr);

		CGPU_EXEC
		void clear();
			
		CGPU_EXEC
		T bx_sum() const;

		CGPU_EXEC
		T bx_min() const;

		CGPU_EXEC
		T bx_max() const;

		CGPU_EXEC
		T bs_x(const T& bs_x_i) const;

		CGPU_EXEC
		T bs(const T& bs_i) const;

		CGPU_EXEC
		T bs_min(const T& bs) const;

		CGPU_EXEC
		T bs_max(const T& bs) const;

		CGPU_EXEC
		T bs_x_h(const T& bs_x) const;

		CGPU_EXEC
		T bs_h(const T& bs) const;

		CGPU_EXEC
		T rx_c(const T& bs_x) const;

		CGPU_EXEC
		T r_c(const T& bs) const;

		CGPU_EXEC
		T radius_x(const T& bs_x) const;

		CGPU_EXEC
		T radius(const T& bs) const;

		CGPU_EXEC
		T radius_x_p(const T& bs_x, const T& p) const;

		CGPU_EXEC
		T radius_p(const T& bs, const T& p) const;

		void set_by_sft(const T& bs, const T& dr);

		void sft_bdr(const T& bs, const T& dr);
	};
}

#include "../src/border_1d.inl"