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

#include "border_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(): Border_Rect_xd<T, edim_2>(), bz_0(0), bz_e(0) {}

	template <class T>
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e, const T& bz_0, const T& bz_e)
	{
		set_in_data(bx_0, bx_e, by_0, by_e, bz_0, bz_e);
	}

	template <class T>
	template <class U> 
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(const dt_init_list<U>& list)
	{
		set_in_data(list);
	}

	template <class T>
	template <class U> 
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(const Vctr_cpu<U>& vctr)
	{
		set_in_data(vctr);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(const Border_Rect_xd<T, edim_3>& border)
	{
		*this = border;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_3>::Border_Rect_xd(const Border_Rect_xd<U, edim_3>& border)
	{
		*this = border;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Border_Rect_xd<T, edim_3>& Border_Rect_xd<T, edim_3>::operator=(const Border_Rect_xd<T, edim_3>& border)
	{
		if (this != &border)
		{
			Border_Rect_xd<T, edim_2>::operator=(border);

			bz_0 = border.bz_0;
			bz_e = border.bz_e;
		}

		return *this;
	}			
			
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Border_Rect_xd<T, edim_3>& Border_Rect_xd<T, edim_3>::operator=(const Border_Rect_xd<U, edim_3>& border)
	{
		if (this != &border)
		{
			Border_Rect_xd<T, edim_2>::operator=(border);

			bz_0 = T(border.bz_0);
			bz_e = T(border.bz_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Border_Rect_xd<T, edim_3>::assign(const Border_Rect_xd<U, edim_3>& border)
	{
		*this = border;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_3>::set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e, const U& bz_0, const U& bz_e)
	{
		Border_Rect_xd<T, edim_2>::set_in_data(bx_0, bx_e, by_0, by_e);

		this->bz_0 = bz_0;
		this->bz_e = bz_e;
	}

	template <class T>
	template <class U> 
	void Border_Rect_xd<T, edim_3>::set_in_data(const dt_init_list<U>& list)
	{
		auto ptr = list.begin();

		Border_Rect_xd<T, edim_2>::set_in_data({ptr[0], ptr[1], ptr[2], ptr[3]});

		bz_0 = T(ptr[4]);
		bz_e = T(ptr[5]);
	}

	template <class T>
	template <class U>
	void Border_Rect_xd<T, edim_3>::set_in_data(const Vctr_cpu<U>& vctr)
	{
		const auto shape = vctr.shape();

		if (shape[0]==1)
		{
			if (shape[1]<2)
				set_in_data({vctr[0], vctr[0], vctr[0], vctr[0], T(0), T(0)});
			else if (shape[1]<3)
				set_in_data({vctr[0], vctr[0], vctr[1], vctr[1], T(0), T(0)});
			else if (shape[1]<4)
				set_in_data({vctr[0], vctr[0], vctr[1], vctr[1], vctr[2], vctr[2]});
			else if (shape[1]<5)
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], T(0), T(0)});		
			else if (shape[1]<6)
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], T(0)});			
			else
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], vctr[5]});
		}
		else if (shape[0]==2)
		{
			if (shape[1]<2)
				set_in_data({vctr[0], vctr[1], T(0), T(0), T(0), T(0)});
			else if (shape[1]<3)
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], T(0), T(0)});
			else
				set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], vctr[5]});
		}
	}

	template <class T>
	CGPU_EXEC
	void Border_Rect_xd<T, edim_3>::clear()
	{
		Border_Rect_xd<T, edim_2>::clear();

		bz_0 = T(0);
		bz_e = T(0);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bz_sum() const
	{
		return bz_0 + bz_e;
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bz_min() const
	{
		return fcn_min(bz_0, bz_e);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bz_max() const
	{
		return fcn_max(bz_0, bz_e);
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bs_z(const T& bs_z_i) const 
	{ 
		return fcn_max(bs_z_i - bz_sum(), T(0));
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Border_Rect_xd<T, edim_3>::bs(const R_3d<T>& bs_i) const 
	{ 
		return {this->bs_x(bs_i.x), this->bs_y(bs_i.y), bs_z(bs_i.z)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bs_min(const R_3d<T>& bs) const 
	{ 
		return fmin(bs(bs));
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bs_max(const R_3d<T>& bs) const 
	{ 
		return fmax(bs(bs));
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::bs_z_h(const T& bs_z) const 
	{ 
		return this->bs_z(bs_z)/T(2);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Border_Rect_xd<T, edim_3>::bs_h(const R_3d<T>& bs) const 
	{ 
		return {this->bs_x_h(bs.x), this->bs_y_h(bs.y), bs_z_h(bs.z)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::rz_c(const T& bs_z) const 
	{ 
		return bz_0 + bs_z_h(bs_z);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Border_Rect_xd<T, edim_3>::r_c(const R_3d<T>& bs) const 
	{ 
		return {this->rx_c(bs.x), this->ry_c(bs.y), rz_c(bs.z)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::radius_z(const T& bs_z) const
	{
		return bs_z_h(bs_z);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Border_Rect_xd<T, edim_3>::radius(const R_3d<T>& bs) const
	{
		return {this->radius_x(bs.x), this->radius_y(bs.y), radius_z(bs.z)};
	}

	template <class T>
	CGPU_EXEC
	T Border_Rect_xd<T, edim_3>::radius_z_p(const T& bs_z, const T& p) const
	{
		return (T(1)-p)*radius_z(bs_z);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Border_Rect_xd<T, edim_3>::radius_p(const R_3d<T>& bs, const T& p) const
	{
		return {this->radius_x_p(bs.x, p), this->radius_y_p(bs.y, p), radius_z_p(bs.z, p)};
	}

	template <class T>
	void Border_Rect_xd<T, edim_3>::set_by_sft(const R_3d<T>& bs, const R_3d<T>& dr)
	{
		Border_Rect_xd<T, edim_2>::set_bz_sft({bs.x, bs.y}, {dr.x, dr.y});

		bz_0 = fcn_set_bound(dr.z, T(0), bs.z);
		bz_e = fcn_set_bound(bs.z - dr.z, T(0), bs.z);
	}

	template <class T>
	void Border_Rect_xd<T, edim_3>::sft_bdr(const R_3d<T>& bs, const R_3d<T>& dr)
	{
		Border_Rect_xd<T, edim_2>::sft_bdr({bs.x, bs.y}, {dr.x, dr.y});

		bz_0 = fcn_set_bound(bz_0 + dr.z, T(0), bs.z);
		bz_e = fcn_set_bound(bz_e - dr.z, T(0), bs.z);
	}
}