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
* GNU General Public License for more details.Gauss_wd_
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

#include "wd_gauss.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDim Dim>
	CGPU_EXEC
	Wd_fcn_xd<T, Dim, efcn_gauss>::Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_gauss>() {}

	template <class T, eDim Dim>
	Wd_fcn_xd<T, Dim, efcn_gauss>::Wd_fcn_xd(const R_xd<T, Dim>& r, const T& sigma, const T& r_wd, const T& r_max) 
	{
		set_in_data(r, sigma, r_wd, r_max);
	}

	/* copy constructor */
	template <class T, eDim Dim>
	CGPU_EXEC
	Wd_fcn_xd<T, Dim, efcn_gauss>::Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_gauss>& wd): Wdb_fcn_xd<T, Dim, efcn_gauss>(wd){}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDim Dim>
	CGPU_EXEC 
	Wd_fcn_xd<T, Dim, efcn_gauss>& Wd_fcn_xd<T, Dim, efcn_gauss>::operator=(const Wdb_fcn_xd<T, Dim, efcn_gauss>& wd)
	{
		Wdb_fcn_xd<T, Dim, efcn_gauss>::operator=(wd);

		return *this;
	}

	/***************************************************************************************/
	template <class T, eDim Dim>
	void Wd_fcn_xd<T, Dim, efcn_gauss>::set_in_data(const R_xd<T, Dim>& r, const T& sigma, const T& r_wd, const T& r_max)
	{
		Fcn_Elem<T, efcn_gauss>::set_in_data(T(1), sigma);
		Wdb_fcn_xd<T, Dim, efcn_gauss>::set_in_data(r, r_wd, r_max);
	}
}