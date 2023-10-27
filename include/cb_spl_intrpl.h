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

#pragma once

#include <vector>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "grid_2d.h"
#include "vctr_cpu.h"
#include "poly_coef_3d.h"

namespace mt
{
	// Cubic spline interpolation
	template <class T>
	class Cb_spl_Intrpl
	{
	public:
		Cb_spl_Intrpl();

		Cb_spl_Intrpl(const pVctr_cpu_64<T>& x, const pVctr_cpu_64<T>& y);

		void set_points(const pVctr_cpu_64<T>& x, const pVctr_cpu_64<T>& y);

		T operator() (T x) const;

		void operator() (const pVctr_cpu_64<T>& x, pVctr_cpu_64<T>& y);

		void get_coef_poly3(Poly_Coef_3d<T, edev_cpu>& ci_coef);

		template <class U, class TVctr> 
		void eval_radial_function(const Grid_2d<U>& grid_2d, U x, U y, TVctr& V);

		T x_min() const;

		T x_max() const;

	private:
		// interpolation parameters
		// f(x) = c3*(x-x_i)^3 + c2*(x-x_i)^2 + c1*(x-x_i) + c0
		Vctr_cpu<T> m_x;
		Vctr_cpu<T> m_c0;
		Vctr_cpu<T> m_c1;
		Vctr_cpu<T> m_c2;
		Vctr_cpu<T> m_c3;
	};
}

#include "../src/cb_spl_intrpl.inl"