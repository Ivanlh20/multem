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
* GNU General Public License for more details.Gauss_wd_
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

#pragma once

#include "const_enum.h"
#include "r_2d.h"
#include "r_3d.h"

/* template definition */
namespace mt
{
#ifndef FCN_ELEM_DEC
	#define FCN_ELEM_DEC
	template <class T, eFcn_typ Fcn_typ> class Fcn_Elem;
#endif
}

/* template definition */
namespace mt
{
	template <class T, eDim Dim, eFcn_typ Fcn_typ> class Wdb_fcn_xd;
}	

/* base class fcn window */
namespace mt
{
	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	class Wdb_fcn_xd: public Fcn_Elem<T, Fcn_typ>
	{
	public:
		using value_type = T;

		R_xd<T, Dim> r;
		T r_wd_2;
		T r_max_2;
		T sft;
		T sc;

		/************************************* constructors ************************************/
		CGPU_EXEC
		Wdb_fcn_xd();

		/* copy constructor */
		CGPU_EXEC
		Wdb_fcn_xd(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC 
		Wdb_fcn_xd<T, Dim, Fcn_typ>& operator=(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd);

		CGPU_EXEC
		void assign(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd);

		/***************************************************************************************/
		CGPU_EXEC
		void clear();

		CGPU_EXEC
		T eval_r2(const T& r2) const;

		CGPU_EXEC
		T eval_r2(const R_2d<T>& r2) const;

		CGPU_EXEC
		T eval_r2(const R_3d<T>& r2) const;

		void set_in_data(const R_xd<T, Dim>& r, const T& r_wd, const T& r_max);
	};
}

#include "../src/wd_fcn.inl"