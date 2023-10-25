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

#include "r_2d.h"
#include "border_1d.h"

/* derived class */
namespace mt
{
	template <class T>
	using Border_Rect_2d = Border_Rect_xd<T, edim_2>;

	using iBorder_Rect_2d = Border_Rect_xd<dt_int32, edim_2>;

	using iBorder_Rect_2d_64 = Border_Rect_xd<dt_int64, edim_2>;
}

/* template specialization 2d */
namespace mt
{
	template <class T>
	class Border_Rect_xd<T, edim_2>: public Border_Rect_xd<T, edim_1>
	{		
	public:			
		using value_type = T;
		using size_type = dt_int32;

		T by_0;		// initial y position
		T by_e;		// final y position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Border_Rect_xd();

		Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e);

		template <class U> 
		Border_Rect_xd(const dt_init_list<U>& list);

		template <class U> 
		Border_Rect_xd(const Vctr_cpu<U>& vctr);

		/* copy constructor */
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<T, edim_2>& border);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<U, edim_2>& border);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Border_Rect_xd<T, edim_2>& operator=(const Border_Rect_xd<T, edim_2>& border);		
			
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd<T, edim_2>& operator=(const Border_Rect_xd<U, edim_2>& border);

		template <class U> 
		CGPU_EXEC
		void assign(const Border_Rect_xd<U, edim_2>& border);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e);

		template <class U> 
		void set_in_data(const dt_init_list<U>& list);

		template <class U>
		void set_in_data(const Vctr_cpu<U>& vctr);

		CGPU_EXEC
		void clear();
	
		CGPU_EXEC
		T by_sum() const;

		CGPU_EXEC
		T by_min() const;

		CGPU_EXEC
		T by_max() const;

		CGPU_EXEC
		T bs_y(const T& bs_y_i) const;

		CGPU_EXEC
		R_2d<T> bs(const R_2d<T>& bs_i) const;

		CGPU_EXEC
		T bs_min(const R_2d<T>& bs) const;

		CGPU_EXEC
		T bs_max(const R_2d<T>& bs) const;

		CGPU_EXEC
		T bs_y_h(const T& bs_y) const;

		CGPU_EXEC
		R_2d<T> bs_h(const R_2d<T>& bs) const;

		CGPU_EXEC
		T ry_c(const T& bs_y) const;

		CGPU_EXEC
		R_2d<T> r_c(const R_2d<T>& bs) const;

		CGPU_EXEC
		T radius_y(const T& bs_y) const;

		CGPU_EXEC
		R_2d<T> radius(const R_2d<T>& bs) const;

		CGPU_EXEC
		T radius_y_p(const T& bs_y, const T& p) const;

		CGPU_EXEC
		R_2d<T> radius_p(const R_2d<T>& bs, const T& p) const;

		void set_by_sft(const R_2d<T>& bs, const R_2d<T>& dr);

		void sft_bdr(const R_2d<T>& bs, const R_2d<T>& dr);
	};
}

#include "../src/border_2d.inl"