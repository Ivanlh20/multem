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

#include "r_3d.h"
#include "border_2d.h"

/* derived class */
namespace mt
{
	template <class T>
	using Border_Rect_3d = Border_Rect_xd<T, edim_3>;

	using iBorder_Rect_3d = Border_Rect_xd<dt_int32, edim_3>;

	using iBorder_Rect_3d_64 = Border_Rect_xd<dt_int64, edim_3>;
}

/* template specialization 3d */
namespace mt
{
	template <class T>
	class Border_Rect_xd<T, edim_3>: public Border_Rect_xd<T, edim_2>
	{		
	public:			
		using value_type = T;
		using size_type = dt_int32;

		T bz_0;		// initial z position
		T bz_e;		// final z position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Border_Rect_xd();

		Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e, const T& bz_0, const T& bz_e);

		template <class U> 
		Border_Rect_xd(const dt_init_list<U>& list);

		template <class U> 
		Border_Rect_xd(const Vctr_cpu<U>& vctr);

		/* copy constructor */
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<T, edim_3>& border);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd(const Border_Rect_xd<U, edim_3>& border);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Border_Rect_xd<T, edim_3>& operator=(const Border_Rect_xd<T, edim_3>& border);			
			
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Border_Rect_xd<T, edim_3>& operator=(const Border_Rect_xd<U, edim_3>& border);

		template <class U> 
		CGPU_EXEC
		void assign(const Border_Rect_xd<U, edim_3>& border);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e, const U& bz_0, const U& bz_e);

		template <class U> 
		void set_in_data(const dt_init_list<U>& list);

		template <class U>
		void set_in_data(const Vctr_cpu<U>& vctr);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		T bz_sum() const;

		CGPU_EXEC
		T bz_min() const;

		CGPU_EXEC
		T bz_max() const;

		CGPU_EXEC
		T bs_z(const T& bs_z_i) const;

		CGPU_EXEC
		R_3d<T> bs(const R_3d<T>& bs_i) const;

		CGPU_EXEC
		T bs_min(const R_3d<T>& bs) const;

		CGPU_EXEC
		T bs_max(const R_3d<T>& bs) const;

		CGPU_EXEC
		T bs_z_h(const T& bs_z) const;

		CGPU_EXEC
		R_3d<T> bs_h(const R_3d<T>& bs) const;

		CGPU_EXEC
		T rz_c(const T& bs_z) const;

		CGPU_EXEC
		R_3d<T> r_c(const R_3d<T>& bs) const;

		CGPU_EXEC
		T radius_z(const T& bs_z) const;

		CGPU_EXEC
		R_3d<T> radius(const R_3d<T>& bs) const;

		CGPU_EXEC
		T radius_z_p(const T& bs_z, const T& p) const;

		CGPU_EXEC
		R_3d<T> radius_p(const R_3d<T>& bs, const T& p) const;

		void set_by_sft(const R_3d<T>& bs, const R_3d<T>& dr);

		void sft_bdr(const R_3d<T>& bs, const R_3d<T>& dr);
	};
}

#include "../src/border_3d.inl"