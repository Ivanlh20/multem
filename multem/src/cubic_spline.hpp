/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef SPLINE_H
#define SPLINE_H

#ifdef _MSC_VER
#pragma once
#endif// _MSC_VER

#include <vector>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "traits.cuh"

namespace mt
{
	// Cubic_Spline interpolation
	template <class T>
	class Cubic_Spline
	{
	private:
		// interpolation parameters
		// f(x) = c3*(x-x_i)^3 + c2*(x-x_i)^2 + c1*(x-x_i) + c0
		Vector<T, e_host> m_x, m_c3, m_c2, m_c1, m_c0;
	public:
		template <class TVector>
		void set_points(const TVector &x, const TVector &y)
		{
			assert(x.size() == y.size());
			assert(x.size()>2);
			int	n =static_cast<int>(x.size());

			m_x.assign(x.begin(), x.end());
			m_c0.assign(y.begin(), y.end());
			m_c1.resize(n);
			m_c2.resize(n);
			m_c3.resize(n);

			Vector<T, e_host> m_h(n);
			Vector<T, e_host> m_l(n);
			Vector<T, e_host> m_mu(n);
			Vector<T, e_host> m_z(n);

			n--;

			for(auto i = 0; i < n; i++)
			{
				m_h[i] = m_x[i+1]-m_x[i];
			}

			m_l[0] = 1;
			m_mu[0] = 0;
			m_z[0] = 0;

			for(auto i = 1; i < n; i++)
			{
				m_l[i] = 2*(m_x[i+1]-m_x[i-1])-m_h[i-1]*m_mu[i-1];
				m_mu[i] = m_h[i]/m_l[i];
				auto alpha = 3*(m_c0[i+1]-m_c0[i])/m_h[i]-3*(m_c0[i]-m_c0[i-1])/m_h[i-1];
				m_z[i] = (alpha-m_h[i-1]*m_z[i-1])/m_l[i];
			}

			m_l[n] = 1;
			m_z[n] = 0;
			m_c2[n] = 0;

			for(auto i =n-1; i >= 0; i--)
			{
				m_c2[i] = m_z[i] - m_mu[i]*m_c2[i+1];
				m_c1[i] = (m_c0[i+1]-m_c0[i])/m_h[i]-m_h[i]*(m_c2[i+1]+2*m_c2[i])/3;
				m_c3[i] = (m_c2[i+1]-m_c2[i])/(3*m_h[i]);
			}
		}

		T operator() (T x) const
		{
			size_t n =m_x.size();
			// find the closest point m_x[idx] < x, idx = 0 even if x<m_x[0]
			// std::vector<T>::const_iterator it;
			auto it =std::lower_bound(m_x.begin(), m_x.end(), x);
			int idx =std::max(int(it-m_x.begin())-1, 0);

			T h =x-m_x[idx];
			return ((m_c3[idx]*h + m_c2[idx])*h + m_c1[idx])*h + m_c0[idx];
		}

		template <class TCI_Coef>
		void get_coeff(TCI_Coef &ci_coef)
		{
			ci_coef.c0.assign(m_c0.begin(), m_c0.end());
			ci_coef.c1.assign(m_c1.begin(), m_c1.end());
			ci_coef.c2.assign(m_c2.begin(), m_c2.end());
			ci_coef.c3.assign(m_c3.begin(), m_c3.end());
		}

		template <class TVector_1, class TVector_2> 
		void eval_function(TVector_1 &x, TVector_2 &y)
		{	
			for(auto i = 0; i<x.size(); i++)
			{
				y[i] = (*this)(x[i]);
			}
		}

		template <class TGrid, class TVector> 
		void eval_radial_function(const TGrid &grid_2d, Value_type<TGrid> x, Value_type<TGrid> y, TVector &V)
		{	
			using X = Value_type<TGrid>;

			X r_max = x_max();
			X fr_max = (*this)(r_max);
			for(auto ix = 0; ix < grid_2d.nx; ix++)
			{
				for(auto iy = 0; iy < grid_2d.ny; iy++)
				{
					X r = grid_2d.R2(ix, iy, x, y);
					int ixy = grid_2d.ind_col(ix, iy);
					V[ixy] = (r < r_max)?(*this)(r):fr_max;
				}
			}
		}

		double x_min() const
		{
			return m_x.front();
		}

		double x_max() const
		{
			return m_x.back();
		}
	};

} // namespace mt

#endif
