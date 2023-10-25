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

#ifndef CB_SPL_INTRPL_H
	#define CB_SPL_INTRPL_H

	#ifdef _MSC_VER
	#pragma once
	#endif

	#include <vector>
	#include <cstdlib>
	#include <cassert>
	#include <algorithm>

	#include "type_traits_gen.h"
	#include "intrpl_coef.cuh"

	namespace mt
	{
		// Cubic spline interpolation
		template <class T>
		class Cb_spl_Intrpl
		{
		public:
			Cb_spl_Intrpl(){}

			Cb_spl_Intrpl(const pVctr_cpu_64<T>& x, const pVctr_cpu_64<T>& y)
			{
				set_points(x, y);
			}

			void set_points(const pVctr_cpu_64<T>& x, const pVctr_cpu_64<T>& y)
			{
				dt_int32 n = static_cast<dt_int32>(x.size());

				m_x.assign(x.begin(), x.end());
				m_c0.assign(y.begin(), y.end());
				m_c1.resize(n);
				m_c2.resize(n);
				m_c3.resize(n);

				Vctr<T, edev_cpu> m_h(n);
				Vctr<T, edev_cpu> cl(n);
				Vctr<T, edev_cpu> m_mu(n);
				Vctr<T, edev_cpu> m_z(n);

				n--;

				for(auto ik = 0; ik < n; ik++)
				{
					m_h[ik] = m_x[ik+1]-m_x[ik];
				}

				cl[0] = 1;
				m_mu[0] = 0;
				m_z[0] = 0;

				for(auto ik = 1; ik < n; ik++)
				{
					cl[ik] = T(2)*(m_x[ik+1]-m_x[ik-1])-m_h[ik-1]*m_mu[ik-1];
					m_mu[ik] = m_h[ik]/cl[ik];
					const auto alpha = T(3)*(m_c0[ik+1]-m_c0[ik])/m_h[ik]-T(3)*(m_c0[ik]-m_c0[ik-1])/m_h[ik-1];
					m_z[ik] = (alpha-m_h[ik-1]*m_z[ik-1])/cl[ik];
				}

				cl[n] = 1;
				m_z[n] = 0;
				m_c2[n] = 0;

				for(auto ik =n-1; ik >= 0; ik--)
				{
					m_c2[ik] = m_z[ik] - m_mu[ik]*m_c2[ik+1];
					m_c1[ik] = (m_c0[ik+1]-m_c0[ik])/m_h[ik]-m_h[ik]*(m_c2[ik+1]+2*m_c2[ik])/3;
					m_c3[ik] = (m_c2[ik+1]-m_c2[ik])/(T(3)*m_h[ik]);
				}
			}

			T operator() (T x) const
			{
				// find the closest point m_x[idx] < m_x, idx = 0 even if m_x<m_x[0]
				auto it = std::lower_bound(m_x.begin(), m_x.end(), x);
				const dt_int32 idx = std::max(dt_int32(it-m_x.begin())-1, 0);

				const auto h = x-m_x[idx];
				return ((m_c3[idx]*h + m_c2[idx])*h + m_c1[idx])*h + m_c0[idx];
			}

			void operator() (const pVctr_cpu_64<T>& x, pVctr_cpu_64<T>& y)
			{	
				for(auto ik = 0; ik<x.size(); ik++)
				{
					y[ik] = this->operator()(x[ik]);
				}
			}

			void get_coef_poly3(Poly_Coef_3d<T, edev_cpu>& ci_coef)
			{
				ci_coef.c0.assign(m_c0.begin(), m_c0.end());
				ci_coef.c1.assign(m_c1.begin(), m_c1.end());
				ci_coef.c2.assign(m_c2.begin(), m_c2.end());
				ci_coef.c3.assign(m_c3.begin(), m_c3.end());
			}

			template <class TGrid, class TVctr> 
			void eval_radial_function(const TGrid& grid_2d, Value_type<TGrid> x, Value_type<TGrid> y, TVctr& V)
			{	
				using X = Value_type<TGrid>;

				X r_max = x_max();
				X fr_max = (*this)(r_max);
				for(auto ix = 0; ix < grid_2d.nx; ix++)
				{
					for(auto iy = 0; iy < grid_2d.ny; iy++)
					{
						X r = grid_2d.r2(ix, iy, x, y);
						dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
						V[ixy] = (r < r_max)?(*this)(r):fr_max;
					}
				}
			}

			T x_min() const
			{
				return m_x.front();
			}

			T x_max() const
			{
				return m_x.back();
			}

		private:
			// interpolation parameters
			// f(x) = c3*(x-x_i)^3 + c2*(x-x_i)^2 + c1*(x-x_i) + c0
			Vctr<T, edev_cpu> m_x;
			Vctr<T, edev_cpu> m_c0;
			Vctr<T, edev_cpu> m_c1;
			Vctr<T, edev_cpu> m_c2;
			Vctr<T, edev_cpu> m_c3;
		};
	}

#endif
