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

#include "region_2d.h"

/* template specialization 2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_2>::Region_Rect_xd(): Region_Rect_xd<T, edim_1>(), ry_0(0), ry_e(0) {}

	template <class T>
	Region_Rect_xd<T, edim_2>::Region_Rect_xd(const T& rx_0, const T& rx_e, const T& ry_0, const T& ry_e)
	{
		set_in_data(rx_0, rx_e, ry_0, ry_e);
	}

	template <class T>
	template <class U>
	Region_Rect_xd<T, edim_2>::Region_Rect_xd(U *data, const size_type& n_data)
	{
		Region_Rect_xd<T, edim_1>(data, n_data);

		ry_0 = (n_data>2)?data[2]:T(0);
		ry_e = (n_data>3)?data[3]:T(0);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_2>::Region_Rect_xd(const Region_Rect_xd<T, edim_2>& region)
	{
		*this = region;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_2>::Region_Rect_xd(const Region_Rect_xd<U, edim_2>& region)
	{
		*this = region;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_2>& Region_Rect_xd<T, edim_2>::operator=(const Region_Rect_xd<T, edim_2>& region)
	{
		if (this != &region)
		{
			Region_Rect_xd<T, edim_1>::operator=(region);

			ry_0 = region.ry_0;
			ry_e = region.ry_e;
		}

		return *this;
	}			
	
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_2>& Region_Rect_xd<T, edim_2>::operator=(const Region_Rect_xd<U, edim_2>& region)
	{
		if (this != &region)
		{
			Region_Rect_xd<T, edim_1>::operator=(region);

			ry_0 = T(region.ry_0);
			ry_e = T(region.ry_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Region_Rect_xd<T, edim_2>::assign(const Region_Rect_xd<U, edim_2>& region)
	{
		*this = region;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Region_Rect_xd<T, edim_2>::set_in_data(const U& rx_0, const U& rx_e, const U& ry_0, const U& ry_e)
	{
		Region_Rect_xd<T, edim_1>::set_in_data(rx_0, rx_e);

		this->ry_0 = ry_0;
		this->ry_e = ry_e;
	}

	template <class T>
	CGPU_EXEC
	void Region_Rect_xd<T, edim_2>::clear()
	{
		Region_Rect_xd<T, edim_1>::clear();

		ry_0 = T(0);
		ry_e = T(0);
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::bs_y() const 
	{ 
		return fcn_max(ry_e-ry_0, T(0));
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Region_Rect_xd<T, edim_2>::bs() const 
	{ 
		return {this->bs_x(), bs_y()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::bs_y_h() const 
	{ 
		return bs_y()/T(2);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Region_Rect_xd<T, edim_2>::bs_h() const 
	{ 
		return {this->bs_x_h(), bs_y_h()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::ry_c() const 
	{ 
		return (ry_e+ry_0)/T(2);
	}
	
	template <class T>
	CGPU_EXEC
	R_2d<T> Region_Rect_xd<T, edim_2>::r_c() const 
	{ 
		return {this->rx_c(), ry_c()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::radius_y() const
	{
		return bs_y_h();
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Region_Rect_xd<T, edim_2>::radius() const
	{
		return {this->radius_x(), radius_y()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::radius_y_p(const T& p) const
	{
		return (T(1)-p)*radius_y();
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Region_Rect_xd<T, edim_2>::radius_p(const T& p) const
	{
		return {this->radius_x_p(p), radius_y_p(p)};
	}

	template <class T>
	void Region_Rect_xd<T, edim_2>::sft_region(const R_2d<T>& bs, const R_2d<T>& dr)
	{
		Region_Rect_xd<T, edim_1>::sft_region(bs.x, dr.x);

		ry_0 = fcn_set_bound(ry_0 + dr.y, T(0), bs.y);
		ry_e = fcn_set_bound(ry_e + dr.y, T(0), bs.y);
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound_y(const T& y) const
	{
		return fcn_chk_bound(y, ry_0, ry_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound(const R_2d<T>& r) const
	{
		return this->chk_bound_x(r.x) && chk_bound_y(r.y);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound(const T& rx, const T& ry)const
	{
		return this->chk_bound_x(rx) && chk_bound_y(ry);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound_y_eps(const T& y) const 
	{ 
		return fcn_chk_bound_eps(y, ry_0, ry_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound_eps(const R_2d<T>& r) const 
	{ 
		return this->chk_bound_x_eps(r.x) && chk_bound_y_eps(r.y);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_2>::chk_bound_eps(const T& rx, const T& ry)const
	{
		return this->chk_bound_x_eps(rx) && chk_bound_y_eps(ry);
	}

	template <class T>
	template <class U, class>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_2>::sub_2_ind(const T& ix, const T& iy) const 
	{ 
		return (iy-ry_0) + (ix-this->rx_0)*bs_y();
	}
}