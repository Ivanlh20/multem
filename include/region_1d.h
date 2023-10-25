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

#include "math_mt.h"
#include "const_enum.h"
#include "type_traits_gen.h"
#include "fcns_cgpu_gen.h"
#include "vctr_cpu.h"

/* template definition */
namespace mt
{
	template <class T, eDim Dim> class Region_Rect_xd;

	template <eDim Dim>
	using iRegion_Rect_xd = Region_Rect_xd<dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class T>
	using Region_Rect_1d = Region_Rect_xd<T, edim_1>;

	using iRegion_Rect_1d = Region_Rect_xd<dt_int32, edim_1>;

	using iRegion_Rect_1d_64 = Region_Rect_xd<dt_int64, edim_1>;
}	

/* template specialization 1d */
namespace mt
{
	template <class T>
	class Region_Rect_xd<T, edim_1>
	{
	public:
		using value_type = T;
		using size_type = dt_int32;

		T rx_0;		// initial x position
		T rx_e;		// final x position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Region_Rect_xd();

		Region_Rect_xd(const T& rx_0, const T& rx_e);

		template <class U>
		Region_Rect_xd(U *data, const size_type& n_data);

		/* copy constructor */
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<T, edim_1>& region);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<U, edim_1>& region);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Region_Rect_xd<T, edim_1>& operator=(const Region_Rect_xd<T, edim_1>& region);		
		
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd<T, edim_1>& operator=(const Region_Rect_xd<U, edim_1>& region);

		template <class U> 
		CGPU_EXEC
		void assign(const Region_Rect_xd<U, edim_1>& region);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& rx_0, const U& rx_e);

		CGPU_EXEC
		void clear();
		
		CGPU_EXEC
		T bs_x() const;

		CGPU_EXEC
		T bs() const;

		CGPU_EXEC
		T bs_x_h() const;

		CGPU_EXEC
		T bs_h() const ;

		CGPU_EXEC
		T rx_c() const;
		
		CGPU_EXEC
		T r_c() const;

		CGPU_EXEC
		T radius_x() const;

		CGPU_EXEC
		T radius() const;

		CGPU_EXEC
		T radius_x_p(const T& p) const;

		CGPU_EXEC
		T radius_p(const T& p) const;

		void sft_region(const T& bs, const T& dr);

		/***************************************************************************************/
		CGPU_EXEC
		dt_bool chk_bound_x(const T& x) const;

		CGPU_EXEC
		dt_bool chk_bound(const T& r) const;

		CGPU_EXEC
		dt_bool chk_bound_x_eps(const T& x) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const T& r) const;

		template <class U=T, class = enable_if_int<U>>
		CGPU_EXEC
		T sub_2_ind(const T& ix) const;
	};
}

#include "../src/region_1d.inl"