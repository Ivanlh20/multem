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

#include "border_1d.h"

/* template specialization 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(): bx_0(0), bx_e(0) {}

	template <class T>
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(const T& bx_0, const T& bx_e)
	{
		set_in_data(bx_0, bx_e);
	}

	template <class T>
	template <class U>
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(const dt_init_list<U>& list)
	{
		set_in_data(list);
	}

	template <class T>
	template <class U>
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(const Vctr_cpu<U>& vctr)
	{
		set_in_data(vctr);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(const Border_Rect_xd<T, edim_1>& border)
	{
		*this = border;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_1>::Border_Rect_xd(const Border_Rect_xd<U, edim_1>& border)
	{
		*this = border;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_1>& Border_Rect_xd<T, edim_1>::operator=(const Border_Rect_xd<T, edim_1>& border)
	{
		if (this != &border)
		{
			bx_0 = border.bx_0;
			bx_e = border.bx_e;
		}

		return *this;
	}			
			
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_1>& Border_Rect_xd<T, edim_1>::operator=(const Border_Rect_xd<U, edim_1>& border)
	{
		if (this != &border)
		{
			bx_0 = T(border.bx_0);
			bx_e = T(border.bx_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Border_Rect_xd<T, edim_1>::assign(const Border_Rect_xd<U, edim_1>& border)
	{
		*this = border;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_1>::set_in_data(const U& bx_0, const U& bx_e)
	{
		this->bx_0 = bx_0;
		this->bx_e = bx_e;
	}

	template <class T>
	template <class U> 
	void Border_Rect_xd<T, edim_1>::set_in_data(const dt_init_list<U>& list)
	{
		auto ptr = list.begin();

		bx_0 = T(ptr[0]);
		bx_e = T(ptr[1]);
	}

	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_1>::set_in_data(const Vctr_cpu<U>& vctr)
	{
		if (vctr.size_32()==1)
		{
			set_in_data({vctr[0], vctr[0]});
		}
		else
		{
			set_in_data({vctr[0], vctr[1]});
		}
	}

	template <class T>
	CGPU_EXEC
	void Border_Rect_xd<T, edim_1>::clear()
	{
		bx_0 = T(0);
		bx_e = T(0);
	}
			
	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bx_sum() const
	{
		return bx_0 + bx_e;
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bx_min() const
	{
		return fcn_min(bx_0, bx_e);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bx_max() const
	{
		return fcn_max(bx_0, bx_e);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs_x(const T& bs_x_i) const 
	{ 
		return fcn_max(bs_x_i - bx_sum(), T(0));
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs(const T& bs_i) const 
	{ 
		return bs_x(bs_i);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs_min(const T& bs) const 
	{ 
		return bs(bs);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs_max(const T& bs) const 
	{ 
		return bs(bs);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs_x_h(const T& bs_x) const 
	{ 
		return this->bs_x(bs_x)/T(2);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::bs_h(const T& bs) const 
	{ 
		return bs_x_h(bs);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::rx_c(const T& bs_x) const 
	{ 
		return bx_0 + bs_x_h(bs_x);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::r_c(const T& bs) const 
	{ 
		return rx_c(bs);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::radius_x(const T& bs_x) const
	{
		return bs_x_h(bs_x);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::radius(const T& bs) const
	{
		return radius_x(bs);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::radius_x_p(const T& bs_x, const T& p) const
	{
		return (T(1)-p)*radius_x(bs_x);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_1>::radius_p(const T& bs, const T& p) const
	{
		return radius_x_p(bs, p);
	}

	template <class T>
	void Border_Rect_xd<T, edim_1>::set_by_sft(const T& bs, const T& dr)
	{
		bx_0 = fcn_set_bound(dr, T(0), bs);
		bx_e = fcn_set_bound(bs - dr, T(0), bs);
	}

	template <class T>
	void Border_Rect_xd<T, edim_1>::sft_bdr(const T& bs, const T& dr)
	{
		bx_0 = fcn_set_bound(bx_0 + dr, T(0), bs);
		bx_e = fcn_set_bound(bx_e - dr, T(0), bs);
	}
}
