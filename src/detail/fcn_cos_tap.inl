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

#include "fcn_cos_tap.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_cos_tap>::Fcn_Elem(): r_tap(0), r_max(0), coef_tap(0){}

		template <class T>
	Fcn_Elem<T, efcn_cos_tap>::Fcn_Elem(const T& r_tap, const T& r_max)
	{
		set_in_data(r_tap, r_max);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_cos_tap>::Fcn_Elem(const Fcn_Elem<T, efcn_cos_tap>& parm)
	{
		*this = parm;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Fcn_Elem<T, efcn_cos_tap>::Fcn_Elem(const Fcn_Elem<U, efcn_cos_tap>& parm)
	{
		*this = parm;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Fcn_Elem<T, efcn_cos_tap>& Fcn_Elem<T, efcn_cos_tap>::operator=(const Fcn_Elem<T, efcn_cos_tap>& parm)
	{
		if (this != &parm)
		{
			r_tap = parm.r_tap;
			r_max = parm.r_max;
			coef_tap = parm.coef_tap;
		}

		return *this;
	}
			
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Fcn_Elem<T, efcn_cos_tap>& Fcn_Elem<T, efcn_cos_tap>::operator=(const Fcn_Elem<U, efcn_cos_tap>& parm)
	{
		r_tap = T(parm.r_tap);
		r_max = T(parm.r_max);
		coef_tap = T(parm.coef_tap);

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	void Fcn_Elem<T, efcn_cos_tap>::assign(const Fcn_Elem<U, efcn_cos_tap>& parm)
	{
		*this = parm;
	}

	/***************************************************************************************/
	template <class T>
	void Fcn_Elem<T, efcn_cos_tap>::set_in_data(const T& r_tap, const T& r_max)
	{
		this->r_tap = r_tap;
		this->r_max = r_max;

		coef_tap = fcn_coef_tap(r_tap, r_max);
	}

	template <class T>
	CGPU_EXEC
	void Fcn_Elem<T, efcn_cos_tap>::clear()
	{
		r_tap = T(0);
		r_max = T(0);
		coef_tap = T(0);
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_cos_tap>::eval_r(const T& r, const T& r_tap, const T& r_max) const
	{ 
		const T coef_tap = fcn_coef_tap(r_tap, r_max);
		return fcn_cos_tap(r_tap, coef_tap, r);
	}

	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_cos_tap>::eval_r(const T& r) const
	{ 
		return fcn_cos_tap(r_tap, coef_tap, r);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Fcn_Elem<T, efcn_cos_tap>::eval_r(const R_2d<T>& r) const
	{ 
		return {eval_r(r.x), eval_r(r.y)};
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Fcn_Elem<T, efcn_cos_tap>::eval_r(const R_3d<T>& r) const
	{ 
		return {eval_r(r.x), eval_r(r.y), eval_r(r.z)};
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_cos_tap>::operator()(const T& r, const T& r_tap, const T& r_max) const
	{ 
		return eval_r(r, r_tap, r_max);
	}

	template <class T>
	CGPU_EXEC
	T Fcn_Elem<T, efcn_cos_tap>::operator()(const T& r) const
	{ 
		return eval_r(r);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> Fcn_Elem<T, efcn_cos_tap>::operator()(const R_2d<T>& r) const
	{ 
		return eval_r(r);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Fcn_Elem<T, efcn_cos_tap>::operator()(const R_3d<T>& r) const
	{ 
		return eval_r(r);
	}
}