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

#include "r_3d.h"
#include "grid_3d.h"

#include "range_2d.h"

/* derived class */
namespace mt
{
	template <class ST>
	using Range_3d_st = Range_xd<ST, edim_3>;

	using Range_3d = Range_xd<dt_int32, edim_3>;

	using Range_3d_64 = Range_xd<dt_int64, edim_3>;	
}

/* template specialization 3d */
namespace mt
{
	template <class ST>
	class Range_xd<ST, edim_3>: public Range_xd<ST, edim_2>
	{
	public:
		using size_type = ST;

		ST iz_0;	// initial z index
		ST nz;		// z points

		/************************************* constructors ************************************/
		CGPU_EXEC
		Range_xd();

		Range_xd(const ST& ix_0, const ST& nx, const ST& iy_0, const ST& ny, const ST& iz_0, const ST& nz);

		/* copy constructor */
		CGPU_EXEC
		Range_xd(const Range_xd<ST, edim_3>& range);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		Range_xd(const Range_xd<SU, edim_3>& range);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Range_xd<ST, edim_3>& operator=(const Range_xd<ST, edim_3>& range);
			
		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		Range_xd<ST, edim_3>& operator=(const Range_xd<SU, edim_3>& range);

		template <class SU> 
		CGPU_EXEC
		void assign(const Range_xd<SU, edim_3>& range);

		/***************************************************************************************/
		template <class SU> 
		void set_in_data(const SU& ix_0, const SU& nx, const SU& iy_0, const SU& ny, const SU& iz_0, const SU& nz);

		template <class T> 
		void set_in_data(const R_3d<T>& r, const T& r_max, const Grid_3d<T>& grid);

		CGPU_EXEC
		void clear();
	};	
}

#include "../src/range_3d.inl"