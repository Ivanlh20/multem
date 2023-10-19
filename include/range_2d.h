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

#include "r_2d.h"
#include "grid_2d.h"

#include "range_1d.h"

/* derived class */
namespace mt
{
	template <class ST>
	using Range_2d_st = Range_xd<ST, edim_2>;

	using Range_2d = Range_xd<dt_int32, edim_2>;

	using Range_2d_64 = Range_xd<dt_int64, edim_2>;
}

/* template specialization 2d */
namespace mt
{
	template <class ST>
	class Range_xd<ST, edim_2>: public Range_xd<ST, edim_1>
	{
	public:
		using size_type = ST;

		ST iy_0;	// initial y index
		ST ny;		// y points

		/************************************* constructors ************************************/
		CGPU_EXEC
		Range_xd();

		Range_xd(const ST& ix_0, const ST& nx, const ST& iy_0, const ST& ny);

		/* copy constructor */
		CGPU_EXEC
		Range_xd(const Range_xd<ST, edim_2>& range);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		Range_xd(const Range_xd<SU, edim_2>& range);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Range_xd<ST, edim_2>& operator=(const Range_xd<ST, edim_2>& range);
			
		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		Range_xd<ST, edim_2>& operator=(const Range_xd<SU, edim_2>& range);

		template <class SU> 
		CGPU_EXEC
		void assign(const Range_xd<SU, edim_2>& range);

		/***************************************************************************************/
		template <class SU> 
		void set_in_data(const SU& ix_0, const SU& nx, const SU& iy_0, const SU& ny);

		template <class T> 
		void set_in_data(const R_2d<T>& r, const T& r_max, const Grid_2d<T>& grid);

		CGPU_EXEC
		void clear();
	};
}

#include "../src/range_2d.inl"