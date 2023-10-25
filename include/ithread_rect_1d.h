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

/* forward declaration */
namespace mt
{
	template <class ST, eDim Dim> class iThread_Rect_sxd;

	template <eDim Dim>
	using iThread_Rect_xd = iThread_Rect_sxd<dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class ST>
	using iThread_Rect_1d_st = iThread_Rect_sxd<ST, edim_1>;

	using iThread_Rect_1d = iThread_Rect_sxd<dt_int32, edim_1>;

	using iThread_Rect_1d_64 = iThread_Rect_sxd<dt_int64, edim_1>;
}

/* template specialization 1d */
namespace mt
{
	template <class ST>
	class iThread_Rect_sxd<ST, edim_1>
	{
	public:
		using size_type = ST;

		ST ind_0;			// initial linear index
		ST ind_e;			// final linear index

		ST ix_0;			// x initial index
		ST ix_e;			// x final index

		dt_int32 istm;		// stream number

		/************************************* constructors ************************************/
		iThread_Rect_sxd();

		iThread_Rect_sxd(const dt_init_list_int32& list);

		iThread_Rect_sxd(const ST& ix_0, const ST& ix_e);

		iThread_Rect_sxd(const ST& nx);

		iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& ind_0, const ST& ind_e);

		/* copy constructor */
		CGPU_EXEC
		iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_1>& ithread);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_1>& ithread);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		iThread_Rect_sxd<ST, edim_1>& operator=(const iThread_Rect_sxd<ST, edim_1>& ithread);			
			
		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		iThread_Rect_sxd<ST, edim_1>& operator=(const iThread_Rect_sxd<SU, edim_1>& ithread);

		template <class SU> 
		CGPU_EXEC
		void assign(const iThread_Rect_sxd<SU, edim_1>& ithread);

		/***************************************************************************************/
		CGPU_EXEC
		void clear();

		CGPU_EXEC
		ST nx() const;

		CGPU_EXEC
		ST size() const;

		CGPU_EXEC
 		dt_bool chk_bound_ix(const ST& ix) const;

		CGPU_EXEC
		dt_bool chk_bound(const ST& ix) const;

		CGPU_EXEC
		void apply_bound_ix(const ST& ix_0, const ST& ix_e);

 		CGPU_EXEC
		ST sub_2_ind(const ST& ix) const;

		CGPU_EXEC
		void set_ind(const ST& ind_0, const ST& ind_e);

		void set_ind_asc_order();
	};
}

#include "../src/ithread_rect_1d.inl"