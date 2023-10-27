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

#include "ptc_exp.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDim Dim>
	CGPU_EXEC
	Ptc_s_fcn_xd<T, Dim, efcn_exp>::Ptc_s_fcn_xd(): Fcn_Elem<T, efcn_exp>(), Fcn_Cos_Tap<T>(), Range<Dim>(), r(), r_max_2(0) {}

	template <class T, eDim Dim>
	Ptc_s_fcn_xd<T, Dim, efcn_exp>::Ptc_s_fcn_xd(const R_xd<T, Dim>& r, const T& a, const T& beta, const T& r_tap, const T& r_max) 
	{
		set_in_data(r, a, beta, r_tap, r_max);
	}

	template <class T, eDim Dim>
	Ptc_s_fcn_xd<T, Dim, efcn_exp>::Ptc_s_fcn_xd(const R_xd<T, Dim>& r, const T& a, const T& beta, const T& r_tap, const T& r_max, const Grid_xd<T, Dim>& grid)
	{
		set_in_data(r, a, beta, r_tap, r_max, grid);
	}

	/* copy constructor */
	template <class T, eDim Dim>
	CGPU_EXEC
	Ptc_s_fcn_xd<T, Dim, efcn_exp>::Ptc_s_fcn_xd(const Ptc_s_fcn_xd<T, Dim, efcn_exp>& ptc)
	{
		*this = ptc;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDim Dim>
	CGPU_EXEC 
	Ptc_s_fcn_xd<T, Dim, efcn_exp>& Ptc_s_fcn_xd<T, Dim, efcn_exp>::operator=(const Ptc_s_fcn_xd<T, Dim, efcn_exp>& ptc)
	{
		if (this != &ptc)
		{
			Fcn_Elem<T, efcn_exp>::operator=(ptc);
			Fcn_Cos_Tap<T>::operator=(ptc);
			Range<Dim>::operator=(ptc);

			r = ptc.r;
			r_max_2 = ptc.r_max_2;
		}

		return *this;
	}

	template <class T, eDim Dim>
	CGPU_EXEC
	void Ptc_s_fcn_xd<T, Dim, efcn_exp>::assign(const Ptc_s_fcn_xd<T, Dim, efcn_exp>& ptc)
	{
		*this = ptc;
	}

	/***************************************************************************************/
	template <class T, eDim Dim>
	void Ptc_s_fcn_xd<T, Dim, efcn_exp>::set_in_data(const R_xd<T, Dim>& r, const T& a, const T& beta, const T& r_tap, const T& r_max)
	{
		Fcn_Elem<T, efcn_exp>::set_in_data(a, beta);
		Fcn_Cos_Tap<T>::set_in_data(::square(r_tap), ::square(r_max));
		this->r = r;
		r_max_2 = ::square(r_max);
	}

	template <class T, eDim Dim>
	void Ptc_s_fcn_xd<T, Dim, efcn_exp>::set_in_data(const R_xd<T, Dim>& r, const T& a, const T& beta, const T& r_tap, const T& r_max, const Grid_xd<T, Dim>& grid)
	{
		set_in_data(r, a, beta, r_tap, r_max);

		Range<Dim>::set_in_data(r, r_max, grid);
	}

	template <class T, eDim Dim>
	CGPU_EXEC
	void Ptc_s_fcn_xd<T, Dim, efcn_exp>::clear()
	{
		Fcn_Elem<T, efcn_exp>::clear();
		Fcn_Cos_Tap<T>::clear();
		Range<Dim>::clear();

		r = T(0);
		r_max_2 = T(0);
	}

	template <class T, eDim Dim>
	CGPU_EXEC
	T Ptc_s_fcn_xd<T, Dim, efcn_exp>::eval_r2(const T& r2) const
	{ 
		return Fcn_Elem<T, efcn_exp>::eval_r2(r2)*Fcn_Cos_Tap<T>::eval_r(r2);
	}
}