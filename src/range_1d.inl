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

#include "range_1d.h"

/* template specialization 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_1>::Range_xd(): ix_0(0), nx(0) {}

	template <class ST>
	Range_xd<ST, edim_1>::Range_xd(const ST& ix_0, const ST& nx) : ix_0(ix_0), nx(nx) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_1>::Range_xd(const Range_xd<ST, edim_1>& range)
	{
		*this = range;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	Range_xd<ST, edim_1>::Range_xd(const Range_xd<SU, edim_1>& range)
	{
		*this = range;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	Range_xd<ST, edim_1>& Range_xd<ST, edim_1>::operator=(const Range_xd<ST, edim_1>& range)
	{
		if (this != &range)
		{
			ix_0 = range.ix_0;
			nx = range.nx;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	Range_xd<ST, edim_1>& Range_xd<ST, edim_1>::operator=(const Range_xd<SU, edim_1>& range)
	{
		if (this != &range)
		{
			ix_0 = ST(range.ix_0);
			nx = ST(range.nx);
		}

		return *this;
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	void Range_xd<ST, edim_1>::assign(const Range_xd<SU, edim_1>& range)
	{
		*this = range;
	}

	/***************************************************************************************/
	template <class ST>
	template <class SU>
	void Range_xd<ST, edim_1>::set_in_data(const SU& ix_0, const SU& nx)
	{
		this->ix_0 = ST(ix_0);
		this->nx = ST(nx);
	}

	template <class ST>
	template <class T>
	void Range_xd<ST, edim_1>::set_in_data(const T& r, const T& r_max, const Grid_1d<T>& grid)
	{
		grid.ix_0_ix_n(r, r_max, ix_0, nx);
	}

	template <class ST>
	CGPU_EXEC
	void Range_xd<ST, edim_1>::clear()
	{
		ix_0 = ST(0);
		nx = ST(0);
	}
}