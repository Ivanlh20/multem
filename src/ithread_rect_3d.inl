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

#include "ithread_rect_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(): iThread_Rect_sxd<ST, edim_2>(), iz_0(0), iz_e(0) {}

	template <class ST>
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const dt_init_list_int32& list)
	{	
		auto ptr = list.begin();

		this->ix_0 = ptr[0];
		this->ix_e = ptr[1];
		this->iy_0 = ptr[2];
		this->iy_e = ptr[3];
		this->iz_0 = ptr[4];
		this->iz_e = ptr[5];

		this->istm = 0;

		this->ind_0 = ST(0);
		this->ind_e = size();
	}
	
	template <class ST>
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& iz_0, const ST& iz_e)
		: iThread_Rect_sxd<ST, edim_2>(ix_0, ix_e, iy_0, iy_e), iz_0(iz_0), iz_e(iz_e)
	{
		this->ind_0 = ST(0);
		this->ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const ST& nx, const ST& ny, const ST& nz)
		: iThread_Rect_sxd<ST, edim_2>(nx, ny), iz_0(0), iz_e(nz)
	{
		this->ind_0 = ST(0);
		this->ind_e = size();
	}

	template <class ST>
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& iz_0, const ST& iz_e, const ST& ind_0, const ST& ind_e)
		: iThread_Rect_sxd<ST, edim_2>(ix_0, ix_e, iy_0, iy_e, ind_0, ind_e), iz_0(iz_0), iz_e(iz_e) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_3>& ithread)
	{
		*this = ithread;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_3>::iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_3>& ithread)
	{
		*this = ithread;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_3>& iThread_Rect_sxd<ST, edim_3>::operator=(const iThread_Rect_sxd<ST, edim_3>& ithread)
	{
		if (this != &ithread)
		{
			iThread_Rect_sxd<ST, edim_2>::operator=(ithread);

			iz_0 = ithread.iz_0;
			iz_e = ithread.iz_e;
		}

		return *this;
	}
			
	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iThread_Rect_sxd<ST, edim_3>& iThread_Rect_sxd<ST, edim_3>::operator=(const iThread_Rect_sxd<SU, edim_3>& ithread)
	{
		if (this != &ithread)
		{
			iThread_Rect_sxd<ST, edim_2>::operator=(ithread);

			iz_0 = ST(ithread.iz_0);
			iz_e = ST(ithread.iz_e);
		}

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_3>::assign(const iThread_Rect_sxd<SU, edim_3>& ithread)
	{
		*this = ithread;
	}

	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_3>::clear()
	{
		iThread_Rect_sxd<ST, edim_2>::clear();

		iz_0 = ST(0);
		iz_e = ST(0);
	}

	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_3>::nz() const 
	{ 
		return iz_e - iz_0;
	}

	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_3>::size() const 
	{ 
		return iThread_Rect_sxd<ST, edim_2>::size()*nz();
	}

	template <class ST>
	CGPU_EXEC
 	dt_bool iThread_Rect_sxd<ST, edim_3>::chk_bound_iz(const ST& iz) const
	{
		return fcn_chk_bound(iz, iz_0, iz_e);
	}

	template <class ST>
	CGPU_EXEC
	dt_bool iThread_Rect_sxd<ST, edim_3>::chk_bound(const ST& ix, const ST& iy, const ST& iz) const
	{
		return this->chk_bound_ix(ix) && this->chk_bound_iy(iy) && chk_bound_iz(iz);
	}

	template <class ST>
	CGPU_EXEC
	void iThread_Rect_sxd<ST, edim_3>::apply_bound_iz(const ST& iz_0, const ST& iz_e)
	{
		this->iz_0 = fcn_set_bound(this->iz_0, iz_0, iz_e);
		this->iz_e = fcn_set_bound(this->iz_e, iz_0, iz_e);
	}

 	template <class ST>
	CGPU_EXEC
	ST iThread_Rect_sxd<ST, edim_3>::sub_2_ind(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return (iy-this->iy_0) + (ix-this->ix_0)*this->ny_s() + (iz-iz_0)*this->ny()*this->nx();
	}

	template <class ST>
	void iThread_Rect_sxd<ST, edim_3>::set_ind_asc_order()
	{ 
		iThread_Rect_sxd<ST, edim_2>::set_ind_asc_order();

		if (iz_0>iz_e)
		{
			std::swap(iz_0, iz_e);
		}
	}
}