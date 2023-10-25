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
#include "region_1d.h"

/* derived class */
namespace mt
{
	template <class T>
	using Region_Rect_2d = Region_Rect_xd<T, edim_2>;

	using iRegion_Rect_2d = Region_Rect_xd<dt_int32, edim_2>;

	using iRegion_Rect_2d_64 = Region_Rect_xd<dt_int64, edim_2>;
}

/* template specialization 2d */
namespace mt
{
	template <class T>
	class Region_Rect_xd<T, edim_2>: public Region_Rect_xd<T, edim_1>
	{		
	public:			
		using value_type = T;
		using size_type = dt_int32;

		T ry_0;		// initial y position
		T ry_e;		// final y position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Region_Rect_xd();

		Region_Rect_xd(const T& rx_0, const T& rx_e, const T& ry_0, const T& ry_e);

		template <class U>
		Region_Rect_xd(U *data, const size_type& n_data);

		/* copy constructor */
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<T, edim_2>& region);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<U, edim_2>& region);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Region_Rect_xd<T, edim_2>& operator=(const Region_Rect_xd<T, edim_2>& region);
		
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd<T, edim_2>& operator=(const Region_Rect_xd<U, edim_2>& region);

		template <class U> 
		CGPU_EXEC
		void assign(const Region_Rect_xd<U, edim_2>& region);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& rx_0, const U& rx_e, const U& ry_0, const U& ry_e);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		T bs_y() const;

		CGPU_EXEC
		R_2d<T> bs() const;

		CGPU_EXEC
		T bs_y_h() const;

		CGPU_EXEC
		R_2d<T> bs_h() const ;

		CGPU_EXEC
		T ry_c() const;
		
		CGPU_EXEC
		R_2d<T> r_c() const;

		CGPU_EXEC
		T radius_y() const;

		CGPU_EXEC
		R_2d<T> radius() const;

		CGPU_EXEC
		T radius_y_p(const T& p) const;

		CGPU_EXEC
		R_2d<T> radius_p(const T& p) const;

		void sft_region(const R_2d<T>& bs, const R_2d<T>& dr);

		/***************************************************************************************/
		CGPU_EXEC
		dt_bool chk_bound_y(const T& y) const;

		CGPU_EXEC
		dt_bool chk_bound(const R_2d<T>& r) const;

		CGPU_EXEC
		dt_bool chk_bound(const T& rx, const T& ry) const;

		CGPU_EXEC
		dt_bool chk_bound_y_eps(const T& y) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const R_2d<T>& r) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const T& rx, const T& ry) const;

		template <class U=T, class = enable_if_int<U>>
		CGPU_EXEC
		T sub_2_ind(const T& ix, const T& iy) const;
	};
}

#include "../src/region_2d.inl"