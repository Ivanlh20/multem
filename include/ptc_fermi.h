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

#include "macros.h"
#include "math_mt.h"
#include "grid_1d.h"
#include "grid_2d.h"
#include "grid_3d.h"
#include "range_1d.h"
#include "range_2d.h"
#include "range_3d.h"

#include "fcn_fermi.h"
#include "fcn_cos_tap.h"

/* template definition */
namespace mt
{
#ifndef PTC_ELEM_DEC
	#define PTC_ELEM_DEC
	template <class T, eDim Dim, eFcn_typ Fcn_typ> class Ptc_s_fcn_xd;
#endif
}

/* derived class */
namespace mt
{
	template <class T, eDim Dim>
	using Ptc_s_Fermi_xd = Ptc_s_fcn_xd<T, Dim, efcn_fermi>;

	template <class T>
	using Ptc_s_Fermi_1d = Ptc_s_fcn_xd<T, edim_1, efcn_fermi>;

	template <class T>
	using Ptc_s_Fermi_2d = Ptc_s_fcn_xd<T, edim_2, efcn_fermi>;

	template <class T>
	using Ptc_s_Fermi_3d = Ptc_s_fcn_xd<T, edim_3, efcn_fermi>;
}

namespace mt
{
	template <class T, eDim Dim>
	class Ptc_s_fcn_xd<T, Dim, efcn_fermi>: public Fcn_Elem<T, efcn_fermi>, public Fcn_Cos_Tap<T>, public Range<Dim>
	{
	public:
		using value_type = T;

		R_xd<T, Dim> r;
		T r_max_2;

		/************************************* constructors ************************************/
		CGPU_EXEC
		Ptc_s_fcn_xd();

		Ptc_s_fcn_xd(const R_xd<T, Dim>& r, const T& a, const T& alpha, const T& radius, const T& r_tap, const T& r_max);

		Ptc_s_fcn_xd(const R_xd<T, Dim>& r, const T& a, const T& alpha, const T& radius, const T& r_tap, const T& r_max, const Grid_xd<T, Dim>& grid);

		/* copy constructor */
		CGPU_EXEC
		Ptc_s_fcn_xd(const Ptc_s_fcn_xd<T, Dim, efcn_fermi>& ptc);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC 
		Ptc_s_fcn_xd<T, Dim, efcn_fermi>& operator=(const Ptc_s_fcn_xd<T, Dim, efcn_fermi>& ptc);

		CGPU_EXEC
		void assign(const Ptc_s_fcn_xd<T, Dim, efcn_fermi>& ptc);

		/***************************************************************************************/
		void set_in_data(const R_xd<T, Dim>& r, const T& a, const T& alpha, const T& radius, const T& r_tap, const T& r_max);

		void set_in_data(const R_xd<T, Dim>& r, const T& a, const T& alpha, const T& radius, const T& r_tap, const T& r_max, const Grid_xd<T, Dim>& grid);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		T eval_r2(const T& r2) const;
	};
}


#include "../src/ptc_fermi.inl"