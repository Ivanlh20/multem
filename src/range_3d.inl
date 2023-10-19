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

#include "range_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_3>::Range_xd(): Range_xd<ST, edim_2>(), iz_0(0), nz(0) {}

	template <class ST>
	Range_xd<ST, edim_3>::Range_xd(const ST& ix_0, const ST& nx, const ST& iy_0, const ST& ny, const ST& iz_0, const ST& nz)
		: Range_xd<ST, edim_2>(ix_0, nx, iy_0, ny), iz_0(iz_0), nz(nz) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_3>::Range_xd(const Range_xd<ST, edim_3>& range)
	{
		*this = range;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	Range_xd<ST, edim_3>::Range_xd(const Range_xd<SU, edim_3>& range)
	{
		*this = range;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_3>&Range_xd<ST, edim_3>:: operator=(const Range_xd<ST, edim_3>& range)
	{
		if (this != &range)
		{
			Range_xd<ST, edim_2>::operator=(range);

			iz_0 = range.iz_0;
			nz = range.nz;
		}

		return *this;
	}
			
	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	Range_xd<ST, edim_3>& Range_xd<ST, edim_3>::operator=(const Range_xd<SU, edim_3>& range)
	{
		if (this != &range)
		{
			Range_xd<ST, edim_2>::operator=(range);

			iz_0 = ST(range.iz_0);
			nz = ST(range.nz);
		}

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void Range_xd<ST, edim_3>::assign(const Range_xd<SU, edim_3>& range)
	{
		*this = range;
	}

	/***************************************************************************************/
	template <class ST>
	template <class SU> 
	void Range_xd<ST, edim_3>::set_in_data(const SU& ix_0, const SU& nx, const SU& iy_0, const SU& ny, const SU& iz_0, const SU& nz)
	{
		Range_xd<ST, edim_2>::set_in_data(ix_0, nx, iy_0, ny);

		this->iz_0 = ST(iz_0);
		this->nz = ST(nz);
	}

	template <class ST>
	template <class T> 
	void Range_xd<ST, edim_3>::set_in_data(const R_3d<T>& r, const T& r_max, const Grid_3d<T>& grid)
	{
		grid.ix_0_ix_n(r.x, r_max, this->ix_0, this->nx);
		grid.iy_0_iy_n(r.y, r_max, this->iy_0, this->ny);
		grid.iz_0_iz_n(r.z, r_max, iz_0, nz);
	}

	template <class ST>
	CGPU_EXEC
	void Range_xd<ST, edim_3>::clear()
	{
		Range_xd<ST, edim_2>::clear();

		iz_0 = ST(0);
		nz = ST(0);
	}
}