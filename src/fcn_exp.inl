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

#include "fcn_exp.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_exp>::Fcn_Elem(): a_1(0), b_1(0) {}

	template <class T>
	Fcn_Elem<T, efcn_exp>::Fcn_Elem(const T& a, const T& beta)
	{
		set_in_data(a, beta);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_exp>::Fcn_Elem(const Fcn_Elem<T, efcn_exp>& parm)
	{
		*this = parm;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_exp>& Fcn_Elem<T, efcn_exp>::operator=(const Fcn_Elem<T, efcn_exp>& parm)
	{
		if (this != &parm)
		{
			a_1 = parm.a_1;
			b_1 = parm.b_1;
		}

		return *this;
	}

	/***************************************************************************************/
	template <class T>
	void Fcn_Elem<T, efcn_exp>::set_in_data(const T& a, const T& beta)
	{
		a_1 = a;
		b_1 = T(1)/(beta);
	}

	template <class T>
	CGPU_EXEC
	void Fcn_Elem<T, efcn_exp>::clear()
	{
		a_1 = T(0);
		b_1 = T(0);
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::eval_r(const T& r, const T& beta) const
	{ 
		return a_1*exp(-r/beta);
	}

	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::eval_r(const T& r) const
	{ 
		return a_1*exp(-b_1*r);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Fcn_Elem<T, efcn_exp>::eval_r(const R_2d<T>& r) const
	{ 
		return {eval_r(r.x), eval_r(r.y)};
	}

	template <class T>
	CGPU_EXEC
	R_3d<T>Fcn_Elem<T, efcn_exp>:: eval_r(const R_3d<T>& r) const
	{ 
		return {eval_r(r.x), eval_r(r.y), eval_r(r.z)};
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::eval_r2(const T& r2, const T& beta) const
	{ 
		return eval_r(::sqrt(r2), beta);
	}

	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::eval_r2(const T& r2) const
	{ 
		return eval_r(::sqrt(r2));
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Fcn_Elem<T, efcn_exp>::eval_r2(const R_2d<T>& r2) const
	{ 
		return eval_r(::sqrt(r2));
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Fcn_Elem<T, efcn_exp>::eval_r2(const R_3d<T>& r2) const
	{ 
		return eval_r(::sqrt(r2));
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::operator()(const T& r, const T& beta) const
	{ 
		return eval_r(r, beta);
	}

	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_exp>::operator()(const T& r) const
	{ 
		return eval_r(r);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Fcn_Elem<T, efcn_exp>::operator()(const R_2d<T>& r) const
	{ 
		return eval_r(r);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Fcn_Elem<T, efcn_exp>::operator()(const R_3d<T>& r) const
	{ 
		return eval_r(r);
	}
}
