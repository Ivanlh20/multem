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

#include "wd_fcn.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	Wdb_fcn_xd<T, Dim, Fcn_typ>::Wdb_fcn_xd(): Fcn_Elem<T, Fcn_typ>(), r(), r_wd_2(0), r_max_2(0), sft(0), sc(1) {}

	/* copy constructor */
	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	Wdb_fcn_xd<T, Dim, Fcn_typ>::Wdb_fcn_xd(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
	{
		*this = wd;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC 
	Wdb_fcn_xd<T, Dim, Fcn_typ>& Wdb_fcn_xd<T, Dim, Fcn_typ>::operator=(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
	{
		if (this != &wd)
		{
			Fcn_Elem<T, Fcn_typ>::operator=(wd);

			r = wd.r;
			r_wd_2 = wd.r_wd_2;
			r_max_2 = wd.r_max_2;
			sft = wd.sft;
			sc = wd.sc;
		}

		return *this;
	}

	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	void Wdb_fcn_xd<T, Dim, Fcn_typ>::assign(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
	{
		*this = wd;
	}

	/***************************************************************************************/
	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	void Wdb_fcn_xd<T, Dim, Fcn_typ>::clear()
	{
		Fcn_Elem<T, Fcn_typ>::clear();

		r = T(0);
		r_wd_2 = T(0);
		r_max_2 = T(0);
		sft = T(0);
		sc = T(1);
	}

	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	T Wdb_fcn_xd<T, Dim, Fcn_typ>::eval_r2(const T& r2) const
	{ 
		return (r2<r_wd_2)?T(1):(r2<r_max_2)?::fmax(T(0), (Fcn_Elem<T, Fcn_typ>::eval_r2(r2) - sft)/sc):T(0);
	}

	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	T Wdb_fcn_xd<T, Dim, Fcn_typ>::eval_r2(const R_2d<T>& r2) const
	{ 
		return eval_r2(r2.x)*eval_r2(r2.y);
	}

	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	CGPU_EXEC
	T Wdb_fcn_xd<T, Dim, Fcn_typ>::eval_r2(const R_3d<T>& r2) const
	{ 
		return eval_r2(r2.x)*eval_r2(r2.y)*eval_r2(r2.z);
	}

	template <class T, eDim Dim, eFcn_typ Fcn_typ>
	void Wdb_fcn_xd<T, Dim, Fcn_typ>::set_in_data(const R_xd<T, Dim>& r, const T& r_wd, const T& r_max)
	{
		this->r = r;
		r_wd_2 = ::square(r_wd);
		r_max_2 = ::square(r_max);
		sft = Fcn_Elem<T, Fcn_typ>::eval_r2(r_max_2);
		sc = Fcn_Elem<T, Fcn_typ>::eval_r2(r_wd_2) - sft;
	}
}