/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef SCAN_PAT_H
	#define SCAN_PAT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "r_2d.cuh"
	#include "cgpu_vctr.cuh"

	namespace mt
	{
		/************************* scanning pattern **************************/
		template <class T>
		struct Scan_Pat
		{
			public:
				using value_type = T;
				using size_type = dt_int32;

				eScan_Pat_Typ typ;		// 1: Line, 2: Area, 3: User_Define
				dt_bool pbc;			// periodic boundary conditions
				dt_bool spxs;			// square pixel size
				R_2d<dt_int32> nsp;		// number of sampling points
				R_2d<T> r_0;			// initial scanning position
				R_2d<T> r_e;			// final scanning position
				Vctr_r_2d_cpu<T> r;		// user define positions

				Scan_Pat(): typ(espt_line), pbc(false), spxs(true), nsp(0, 0), 
				r_0(0, 0), r_e(0, 0), dr(0, 0){};

				template <class U> 
				void assign(Scan_Pat<U>& scan)
				{
					if (this !=  &scan)
					{
						typ = scan.typ;
						pbc = scan.pbc;
						spxs = scan.spxs;
						nsp = scan.nsp;
						r_0 = scan.r_0;
						r_e = scan.r_e;

						dr = scan.dr;
						r = scan.r;
					}
				}

				void set_in_data(const eScan_Pat_Typ& typ, const dt_bool& pbc, const dt_bool& spxs, 
					const R_2d<dt_int32>& nsp, const R_2d<T>& r_0, const R_2d<T>& r_e, const Vctr_r_2d_cpu<T>& r)
				{
					this->typ = typ;
					this->pbc = pbc;
					this->spxs = spxs;
					this->nsp = nsp;
					this->r_0 = r_0;
					this->r_e = r_e;
					this->r = r;

					set_dep_var();
				}

				void set_dep_var()
				{
					if ((nsp.x <= 0) || (nsp.y <= 0))
					{
						clear();
						return;
					}

					if (is_scan_pat_line())
					{
						set_grid_scan_pat_line();
					}
					else if (is_scan_pat_area())
					{
						set_grid_scan_pat_area();
					}
					else if (is_scan_pat_user_def())
					{
						nsp.x = r.size();
						nsp.y = nsp.x;
					}
				}

				template <class U> 
				Scan_Pat<T>& operator=(Scan_Pat<U>& scan)
				{
					assign(scan);
					return *this;
				}

				size_type size() const
				{
					return nsp.x*nsp.y;
				}

				void clear()
				{
					typ = espt_line;
					pbc = false;
					spxs = true;
					nsp.x = 0;
					nsp.y = 0;
					r_0 = 0;
					r_e = 0;
					r.clear_shrink_to_fit();
					dr = 0;
				}

				R_2d<T> operator()(const dt_int32& ind) const 
				{
					if (is_scan_pat_user_def())
					{
						return r_ud[ind];
					}
					else
					{
						const dt_int32 ix = ind/nsp.y;
						const dt_int32 iy = ind - ix*nsp.y;

						return (r_0 + R_2d<T>(ix, iy)*dr);
					}
				}

				Vctr_cpu<T> rx_vctr() const
				{
					Vctr_cpu<T> x;
					x.reserve(nsp.x);
					if (is_scan_pat_user_def())
					{
						for(auto ix = 0; ix<nsp.x; ix++)
						{
							x.push_back(r[ix].x);
						}
					}
					else
					{
						for(auto ix = 0; ix<nsp.x; ix++)
						{
							x.push_back(r_0.x + T(ix)*dr.x);
						}
					}
					return x;
				}

				Vctr_cpu<T> ry_vctr() const
				{
					Vctr_cpu<T> y;
					y.reserve(nsp.y);
					if (is_scan_pat_user_def())
					{
						for(auto iy = 0; iy<nsp.y; iy++)
						{
							y.push_back(r[iy].y);
						}
					}
					else
					{
						for(auto iy = 0; iy<nsp.y; iy++)
						{
							y.push_back(r_0.y + T(iy)*dr.y);
						}
					}
					return y;
				}

				dt_bool is_scan_pat_line() const
				{
					return mt::is_scan_pat_line(typ);
				}

				dt_bool is_scan_pat_area() const
				{
					return mt::is_scan_pat_area(typ);
				}

 				dt_bool is_scan_pat_user_def() const
				{
					return mt::is_scan_pat_user_def(typ);
				}
			private:
				R_2d<T> dr;	// pixel size

				void set_grid_scan_pat_line()
				{
					const auto r_u = r_e-r_0;
					dr = r_u/((pbc)?nsp:(nsp-1));
				}

				void set_grid_scan_pat_area()
				{
					const auto r_u = r_e-r_0;
					dr = r_u/((pbc)?nsp:(nsp-1));

					if (spxs)
					{
						if (fabs(r_u.x)>fabs(r_u.y))
						{
							dr.y = std::copysign(dr.x, r_u.y);
							nsp.y = fcn_cfloor<dt_int32>(r_u.y/dr.y+Epsilon<T>::rel+0.5);
							nsp.y += (pbc)?0:1;
						}
						else
						{
							dr.x = std::copysign(dr.y, r_u.x);
							nsp.x = fcn_cfloor<dt_int32>(r_u.x/dr.x+Epsilon<T>::rel+0.5);
							nsp.x += (pbc)?0:1;
						}
					}
				}

		};
	}

#endif