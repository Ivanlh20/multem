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

#include "ithread_rect_2d.h"

/* template specialization 2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(): iThread_Rect_sxd<ST, edim_1>(), iy_0(0), iy_e(0) {}

	template <class ST>
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const dt_init_list_int32& list)
	{	
		auto ptr = list.begin();

		this->ix_0 = ptr[0];
		this->ix_e = ptr[1];
		this->iy_0 = ptr[2];
		this->iy_e = ptr[3];

		this->istm = 0;

		this->ind_0 = ST(0);
		this->ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e)
		: iThread_Rect_sxd<ST, edim_1>(ix_0, ix_e), iy_0(iy_0), iy_e(iy_e)
	{
		this->ind_0 = ST(0);
		this->ind_e = size();
	}
			
	template <class ST>
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const ST& nx, const ST& ny)
		: iThread_Rect_sxd<ST, edim_1>(nx), iy_0(0), iy_e(ny)
	{
		this->ind_0 = ST(0);
		this->ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& ind_0, const ST& ind_e)
		: iThread_Rect_sxd<ST, edim_1>(ix_0, ix_e, ind_0, ind_e), iy_0(iy_0), iy_e(iy_e) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_2>& ithread)
	{
		*this = ithread;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_2>::iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_2>& ithread)
	{
		*this = ithread;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_2>& iThread_Rect_sxd<ST, edim_2>::operator=(const iThread_Rect_sxd<ST, edim_2>& ithread)
	{
		if (this != &ithread)
		{
			iThread_Rect_sxd<ST, edim_1>::operator=(ithread);

			iy_0 = ithread.iy_0;
			iy_e = ithread.iy_e;
		}

		return *this;
	}
			
	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_2>& iThread_Rect_sxd<ST, edim_2>::operator=(const iThread_Rect_sxd<SU, edim_2>& ithread)
	{
		if (this != &ithread)
		{
			iThread_Rect_sxd<ST, edim_1>::operator=(ithread);

			iy_0 = ST(ithread.iy_0);
			iy_e = ST(ithread.iy_e);
		}

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_2>::assign(const iThread_Rect_sxd<SU, edim_2>& ithread)
	{
		*this = ithread;
	}

	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_2>::clear()
	{
		iThread_Rect_sxd<ST, edim_1>::clear();

		iy_0 = ST(0);
		iy_e = ST(0);
	}


	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_2>::ny() const 
	{ 
		return iy_e - iy_0;
	}

	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_2>::size() const 
	{ 
		return iThread_Rect_sxd<ST, edim_1>::size()*ny();
	}

	template <class ST>
	CGPU_EXEC
 	dt_bool iThread_Rect_sxd<ST, edim_2>::chk_bound_iy(const ST& iy) const
	{
		return fcn_chk_bound(iy, iy_0, iy_e);
	}

	template <class ST>
	CGPU_EXEC
	dt_bool iThread_Rect_sxd<ST, edim_2>::chk_bound(const ST& ix, const ST& iy) const
	{
		return this->chk_bound_ix(ix) && chk_bound_iy(iy);
	}

	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_2>::apply_bound_iy(const ST& iy_0, const ST& iy_e)
	{
		this->iy_0 = fcn_set_bound(this->iy_0, iy_0, iy_e);
		this->iy_e = fcn_set_bound(this->iy_e, iy_0, iy_e);
	}

	template <class ST>
 	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_2>::sub_2_ind(const ST& ix, const ST& iy) const 
	{ 
		return (iy-iy_0) + (ix-this->ix_0)*ny();
	}

	template <class ST>
	void iThread_Rect_sxd<ST, edim_2>::set_ind_asc_order()
	{ 
		iThread_Rect_sxd<ST, edim_1>::set_ind_asc_order();

		if (iy_0>iy_e)
		{
			std::swap(iy_0, iy_e);
		}
	}
}