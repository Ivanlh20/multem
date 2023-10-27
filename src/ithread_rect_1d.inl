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

#include "ithread_rect_1d.h"

/* template specialization 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(): ix_0(0), ix_e(0), ind_0(0), ind_e(0), istm(0) {}

	template <class ST>
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const dt_init_list_int32& list): istm(0)
	{	
		auto ptr = list.begin();

		ix_0 = ptr[0];
		ix_e = ptr[1];

		ind_0 = ST(0);
		ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e): ix_0(ix_0), ix_e(ix_e), istm(0)
	{
		ind_0 = ST(0);
		ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const ST& nx): ix_0(0), ix_e(nx), istm(0)
	{
		ind_0 = ST(0);
		ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& ind_0, const ST& ind_e)
		: ix_0(ix_0), ix_e(ix_e), ind_0(ind_0), ind_e(ind_e), istm(0) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_1>& ithread)
	{
		*this = ithread;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_1>::iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_1>& ithread)
	{
		*this = ithread;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_1>& iThread_Rect_sxd<ST, edim_1>::operator=(const iThread_Rect_sxd<ST, edim_1>& ithread)
	{
		if (this != &ithread)
		{
			ix_0 = ithread.ix_0;
			ix_e = ithread.ix_e;

			ind_0 = ithread.ind_0;
			ind_e = ithread.ind_e;

			istm = ithread.istm;
		}

		return *this;
	}			
			
	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_1>& iThread_Rect_sxd<ST, edim_1>::operator=(const iThread_Rect_sxd<SU, edim_1>& ithread)
	{
		if (this != &ithread)
		{
			ix_0 = ST(ithread.ix_0);
			ix_e = ST(ithread.ix_e);

			ind_0 = ST(ithread.ind_0);
			ind_e = ST(ithread.ind_e);

			istm = ithread.istm;
		}

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_1>::assign(const iThread_Rect_sxd<SU, edim_1>& ithread)
	{
		*this = ithread;
	}

	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_1>::clear()
	{
		ix_0 = ST(0);
		ix_e = ST(0);

		ind_0 = ST(0);
		ind_e = ST(0);

		istm = 0;
	}

	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_1>::nx() const 
	{ 
		return ix_e - ix_0;
	}

	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_1>::size() const 
	{ 
		return nx();
	}

	template <class ST>
	CGPU_EXEC
 	dt_bool iThread_Rect_sxd<ST, edim_1>::chk_bound_ix(const ST& ix) const
	{
		return fcn_chk_bound(ix, ix_0, ix_e);
	}

	template <class ST>
	CGPU_EXEC
	dt_bool iThread_Rect_sxd<ST, edim_1>::chk_bound(const ST& ix) const
	{
		return chk_bound_ix(ix);
	}

	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_1>::apply_bound_ix(const ST& ix_0, const ST& ix_e)
	{
		this->ix_0 = fcn_set_bound(this->ix_0, ix_0, ix_e);
		this->ix_e = fcn_set_bound(this->ix_e, ix_0, ix_e);
	}

 	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_1>::sub_2_ind(const ST& ix) const 
	{ 
		return ix-ix_0;
	}

	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_1>::set_ind(const ST& ind_0, const ST& ind_e)
	{ 
		this->ind_0 = ind_0;
		this->ind_e = ind_e;
	}

	template <class ST>
	void iThread_Rect_sxd<ST, edim_1>::set_ind_asc_order()
	{ 
		if (ix_0>ix_e)
		{
			std::swap(ix_0, ix_e);
		}
	}
}