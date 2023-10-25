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

#include "border_2d.h"

/* template specialization 2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(): Border_Rect_xd<T, edim_1>(), by_0(0), by_e(0) {}

	template <class T>
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e)
	{
		set_in_data(bx_0, bx_e, by_0, by_e);
	}

	template <class T>
	template <class U> 
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(const dt_init_list<U>& list)
	{
		set_in_data(list);
	}

	template <class T>
	template <class U> 
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(const Vctr_cpu<U>& vctr)
	{
		set_in_data(vctr);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(const Border_Rect_xd<T, edim_2>& border)
	{
		*this = border;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_2>::Border_Rect_xd(const Border_Rect_xd<U, edim_2>& border)
	{
		*this = border;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_2>& Border_Rect_xd<T, edim_2>::operator=(const Border_Rect_xd<T, edim_2>& border)
	{
		if (this != &border)
		{
			Border_Rect_xd<T, edim_1>::operator=(border);

			by_0 = border.by_0;
			by_e = border.by_e;
		}

		return *this;
	}			
			
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_2>& Border_Rect_xd<T, edim_2>::operator=(const Border_Rect_xd<U, edim_2>& border)
	{
		if (this != &border)
		{
			Border_Rect_xd<T, edim_1>::operator=(border);

			by_0 = T(border.by_0);
			by_e = T(border.by_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Border_Rect_xd<T, edim_2>::assign(const Border_Rect_xd<U, edim_2>& border)
	{
		*this = border;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_2>::set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e)
	{
		Border_Rect_xd<T, edim_1>::set_in_data(bx_0, bx_e);

		this->by_0 = by_0;
		this->by_e = by_e;
	}

	template <class T>
	template <class U> 
	void Border_Rect_xd<T, edim_2>::set_in_data(const dt_init_list<U>& list)
	{
		auto ptr = list.begin();

		Border_Rect_xd<T, edim_1>::set_in_data({ptr[0], ptr[1]});

		by_0 = T(ptr[2]);
		by_e = T(ptr[3]);
	}

	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_2>::set_in_data(const Vctr_cpu<U>& vctr)
	{
		const auto shape = vctr.shape();

		if (shape[0]==1)
		{
			if (shape[1]<2)
				set_in_data({vctr[0], vctr[0], vctr[0], vctr[0]});
			else if (shape[1]<3)
				set_in_data({vctr[0], vctr[0], vctr[1], vctr[1]});
			else if (shape[1]<4)
				set_in_data({vctr[0], vctr[1], vctr[2], T(0)});			
			else
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3]});	
		}
		else if (shape[0]==2)
		{
			if (shape[1]<2)
				set_in_data({vctr[0], vctr[1], T(0), T(0)});
			else
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3]});
		}
	}

	template <class T>
	CGPU_EXEC
	void Border_Rect_xd<T, edim_2>::clear()
	{
		Border_Rect_xd<T, edim_1>::clear();

		by_0 = T(0);
		by_e = T(0);
	}
	
	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::by_sum() const
	{
		return by_0 + by_e;
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::by_min() const
	{
		return fcn_min(by_0, by_e);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::by_max() const
	{
		return fcn_max(by_0, by_e);
	}
	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::bs_y(const T& bs_y_i) const 
	{ 
		return fcn_max(bs_y_i - by_sum(), T(0));
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Border_Rect_xd<T, edim_2>::bs(const R_2d<T>& bs_i) const 
	{ 
		return {this->bs_x(bs_i.x), bs_y(bs_i.y)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::bs_min(const R_2d<T>& bs) const 
	{ 
		return fmin(bs(bs));
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::bs_max(const R_2d<T>& bs) const 
	{ 
		return fmax(bs(bs));
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::bs_y_h(const T& bs_y) const 
	{ 
		return this->bs_y(bs_y)/T(2);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Border_Rect_xd<T, edim_2>::bs_h(const R_2d<T>& bs) const 
	{ 
		return {this->bs_x_h(bs.x), bs_y_h(bs.y)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::ry_c(const T& bs_y) const 
	{ 
		return by_0 + bs_y_h(bs_y);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Border_Rect_xd<T, edim_2>::r_c(const R_2d<T>& bs) const 
	{ 
		return {this->rx_c(bs.x), ry_c(bs.y)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::radius_y(const T& bs_y) const
	{
		return bs_y_h(bs_y);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Border_Rect_xd<T, edim_2>::radius(const R_2d<T>& bs) const
	{
		return {this->radius_x(bs.x), radius_y(bs.y)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_2>::radius_y_p(const T& bs_y, const T& p) const
	{
		return (T(1)-p)*radius_y(bs_y);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Border_Rect_xd<T, edim_2>::radius_p(const R_2d<T>& bs, const T& p) const
	{
		return {this->radius_x_p(bs.x, p), radius_y_p(bs.y, p)};
	}

	template <class T>
	void Border_Rect_xd<T, edim_2>::set_by_sft(const R_2d<T>& bs, const R_2d<T>& dr)
	{
		Border_Rect_xd<T, edim_1>::set_by_sft(bs.x, dr.x);

		by_0 = fcn_set_bound(dr.y, T(0), bs.y);
		by_e = fcn_set_bound(bs.y - dr.y, T(0), bs.y);
	}

	template <class T>
	void Border_Rect_xd<T, edim_2>::sft_bdr(const R_2d<T>& bs, const R_2d<T>& dr)
	{
		Border_Rect_xd<T, edim_1>::sft_bdr(bs.x, dr.x);

		by_0 = fcn_set_bound(by_0 + dr.y, T(0), bs.y);
		by_e = fcn_set_bound(by_e - dr.y, T(0), bs.y);
	}
}