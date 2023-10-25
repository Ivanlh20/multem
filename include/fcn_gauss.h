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
#include "math_mt.h"
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

/* derived class */
namespace mt
{
	template <class T>
	using Fcn_Gauss = Fcn_Elem<T, efcn_gauss>;
}

namespace mt
{
	template <class T>
	class Fcn_Elem<T, efcn_gauss>
	{
	public:
		using value_type = T;

		T a_1;
		T b_1;

		/************************************* constructors ************************************/
		CGPU_EXEC
		Fcn_Elem();

		Fcn_Elem(const T& sigma);

		Fcn_Elem(const T& a, const T& sigma);

		/* copy constructor */
		CGPU_EXEC
		Fcn_Elem(const Fcn_Elem<T, efcn_gauss>& parm);

		/* copy assignment operator */
		CGPU_EXEC
		Fcn_Elem<T, efcn_gauss>& operator=(const Fcn_Elem<T, efcn_gauss>& parm);

		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Fcn_Elem<T, efcn_gauss>& operator=(const Fcn_Elem<U, efcn_gauss>& parm);

		template <class U>
		CGPU_EXEC
		void assign(const Fcn_Elem<U, efcn_gauss>& parm);

		/***************************************************************************************/
		void set_in_data(const T& a, const T& sigma);

		CGPU_EXEC
		void clear();

		/***************************************************************************************/
		CGPU_EXEC
		T eval_r2(const T& r2, const T& sigma) const;

		CGPU_EXEC
		T eval_r2(const T& r2) const;

		CGPU_EXEC
		R_2d<T> eval_r2(const R_2d<T>& r2) const;

		CGPU_EXEC
		R_3d<T> eval_r2(const R_3d<T>& r2) const;

		/***************************************************************************************/
		CGPU_EXEC
		T eval_r(const T& r, const T& sigma) const;

		CGPU_EXEC
		T eval_r(const T& r) const;

		CGPU_EXEC
		R_2d<T> eval_r(const R_2d<T>& r) const;

		CGPU_EXEC
		R_3d<T> eval_r(const R_3d<T>& r) const;

		/***************************************************************************************/
		CGPU_EXEC
		T operator()(const T& r2, const T& sigma) const;

		CGPU_EXEC
		T operator()(const T& r2) const;

		CGPU_EXEC
		R_2d<T> operator()(const R_2d<T>& r2) const;

		CGPU_EXEC
		R_3d<T> operator()(const R_3d<T>& r2) const;
	};
}

#include "../src/fcn_gauss.inl"