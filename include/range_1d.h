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

#include "const_enum.h"

#include "grid_1d.h"

/* template definition */
namespace mt
{
	template <class ST, eDim Dim> class Range_xd;

	template <eDim Dim>
	using Range = Range_xd<dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class ST>
	using Range_1d_st = Range_xd<ST, edim_1>;

	using Range_1d = Range_xd<dt_int32, edim_1>;

	using Range_1d_64 = Range_xd<dt_int64, edim_1>;
}

/* template specialization 1d */
namespace mt
{
	template <class ST>
	class Range_xd<ST, edim_1>
	{
	public:
		using size_type = ST;

		ST ix_0;		// initial x index
		ST nx;			// x points

		/************************************* constructors ************************************/
		CGPU_EXEC
		Range_xd();

		Range_xd(const ST& ix_0, const ST& nx);

		/* copy constructor */
		CGPU_EXEC
		Range_xd(const Range_xd<ST, edim_1>& range);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		Range_xd(const Range_xd<SU, edim_1>& range);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Range_xd<ST, edim_1>& operator=(const Range_xd<ST, edim_1>& range);

		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		Range_xd<ST, edim_1>& operator=(const Range_xd<SU, edim_1>& range);

		template <class SU>
		CGPU_EXEC
		void assign(const Range_xd<SU, edim_1>& range);

		/***************************************************************************************/
		template <class SU>
		void set_in_data(const SU& ix_0, const SU& nx);

		template <class T>
		void set_in_data(const T& r, const T& r_max, const Grid_1d<T>& grid);

		CGPU_EXEC
		void clear();
	};
}

#include "../src/range_1d.inl"