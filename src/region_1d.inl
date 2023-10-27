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

#include "region_1d.h"

/* template specialization 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_1>::Region_Rect_xd(): rx_0(0), rx_e(0) {}

	template <class T>
	Region_Rect_xd<T, edim_1>::Region_Rect_xd(const T& rx_0, const T& rx_e)
	{
		set_in_data(rx_0, rx_e);
	}

	template <class T>
	template <class U>
	Region_Rect_xd<T, edim_1>::Region_Rect_xd(U *data, const size_type& n_data)
	{
		rx_0 = (n_data>0)?data[0]:T(0);
		rx_e = (n_data>1)?data[1]:T(0);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_1>::Region_Rect_xd(const Region_Rect_xd<T, edim_1>& region)
	{
		*this = region;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_1>::Region_Rect_xd(const Region_Rect_xd<U, edim_1>& region)
	{
		*this = region;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_1>& Region_Rect_xd<T, edim_1>::operator=(const Region_Rect_xd<T, edim_1>& region)
	{
		if (this != &region)
		{
			rx_0 = region.rx_0;
			rx_e = region.rx_e;
		}

		return *this;
	}			
	
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_1>& Region_Rect_xd<T, edim_1>::operator=(const Region_Rect_xd<U, edim_1>& region)
	{
		if (this != &region)
		{
			rx_0 = T(region.rx_0);
			rx_e = T(region.rx_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Region_Rect_xd<T, edim_1>::assign(const Region_Rect_xd<U, edim_1>& region)
	{
		*this = region;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Region_Rect_xd<T, edim_1>::set_in_data(const U& rx_0, const U& rx_e)
	{
		this->rx_0 = rx_0;
		this->rx_e = rx_e;
	}

	template <class T>
	CGPU_EXEC
	void Region_Rect_xd<T, edim_1>::clear()
	{
		rx_0 = T(0);
		rx_e = T(0);
	}
	
	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::bs_x() const 
	{ 
		return fcn_max(rx_e-rx_0, T(0));
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::bs() const 
	{ 
		return bs_x();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::bs_x_h() const 
	{ 
		return bs_x()/T(2);
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::bs_h() const 
	{ 
		return bs_x_h();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::rx_c() const 
	{ 
		return (rx_e+rx_0)/T(2);
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::r_c() const 
	{ 
		return rx_c();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::radius_x() const
	{
		return bs_x_h();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::radius() const
	{
		return radius_x();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::radius_x_p(const T& p) const
	{
		return (T(1)-p)*radius_x();
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::radius_p(const T& p) const
	{
		return radius_x_p(p);
	}

	template <class T>
	void Region_Rect_xd<T, edim_1>::sft_region(const T& bs, const T& dr)
	{
		rx_0 = fcn_set_bound(rx_0 + dr, T(0), bs);
		rx_e = fcn_set_bound(rx_e + dr, T(0), bs);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_1>::chk_bound_x(const T& x) const
	{
		return fcn_chk_bound(x, rx_0, rx_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_1>::chk_bound(const T& r) const
	{
		return chk_bound_x(r);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_1>::chk_bound_x_eps(const T& x) const 
	{ 
		return fcn_chk_bound_eps(x, rx_0, rx_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_1>::chk_bound_eps(const T& r) const 
	{ 
		return chk_bound_x_eps(r);
	}

	template <class T>
	template <class U, class>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_1>::sub_2_ind(const T& ix) const 
	{ 
		return ix-rx_0;
	}
}