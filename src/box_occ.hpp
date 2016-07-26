/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOX_OCC_H
#define BOX_OCC_H

#include <vector>
#include <deque>

#include "types.cuh"
#include "math.cuh"
#include "lin_alg_def.cuh"
#include "atom_data.hpp"

namespace mt
{
	class Box_Occ
	{
		public:
			Box_Occ(): r_min(0), r2_min(0), a_min(0), lx(0), ly(0), lz(0), nx(0), ny(0), nz(0), nxy(0){};

			void init()
			{
				for(auto i = 0; i<occ.size(); i++)
				{
					if(occ[i]>-2)
					{
						occ[i] = -1;
					}
				}
			}

			void set_input_data(double r_min_i, double lx_i, double ly_i, double lz_i)
			{	
				r_min = r_min_i;
				r2_min = pow(r_min, 2);
				a_min = r_min/mt::c_3i2;
				lx = lx_i;
				ly = ly_i;
				lz = lz_i;
				nx = static_cast<int>(ceil(lx/a_min));
				ny = static_cast<int>(ceil(ly/a_min));
				nz = static_cast<int>(ceil(lz/a_min));
				nxy = nx*ny;
				occ.resize(nxy*nz, -1);
			}

			inline
			int64_t xyz_2_ind(const int &ix, const int &iy, const int &iz)
			{
				return int64_t(iz)*nxy+ int64_t(iy)*nx + int64_t(ix);
			}

			template<class T>
			inline void set_occ(const r3d<T> &r, const int &iatom)
			{
				int ix = static_cast<int>(floor(r.x/a_min));
				int iy = static_cast<int>(floor(r.y/a_min));
				int iz = static_cast<int>(floor(r.z/a_min));
				auto ixyz = xyz_2_ind(ix, iy, iz);
				occ[ixyz] = iatom;
			}

			template<class T>
			inline bool check_r_min(const Atom_Data_Sa<T> &atoms, const r3d<T> &r_i)
			{
				int ix_i = static_cast<int>(floor(r_i.x/a_min));
				int iy_i = static_cast<int>(floor(r_i.y/a_min));
				int iz_i = static_cast<int>(floor(r_i.z/a_min));
				auto ixyz = xyz_2_ind(ix_i, iy_i, iz_i);

				if (occ[ixyz]>-1)
				{
					return false;
				}

				auto get_limits = [](const int &k, const int &nk, int &k_0, int &k_e)->void
				{
					k_0 = (k-1<0)?0:(k-1);
					k_e = (k+2>nk)?nk:(k+2);
				};

				int ix_0, ix_e;
				get_limits(ix_i, nx, ix_0, ix_e);
				int iy_0, iy_e;
				get_limits(iy_i, ny, iy_0, iy_e);
				int iz_0, iz_e;
				get_limits(iz_i, nz, iz_0, iz_e);

				for(auto iz = iz_0; iz<iz_e; iz++)
				{
					for(auto iy = iy_0; iy<iy_e; iy++)
					{
						for(auto ix = ix_0; ix<ix_e; ix++)
						{
							auto iatom = occ[xyz_2_ind(ix, iy, iz)];
							if((iatom>-1)&&(norm(atoms.r_n[iatom]-r_i)<r2_min))
							{
								return false;
							}
						}
					}
				}

				return true;
			}

		private:
			double r_min;
			double r2_min;
			double a_min;
			double lx;
			double ly;
			double lz;
			int nx;
			int ny;
			int nz;
			int64_t nxy;

			Vector<int, e_host> occ;
	};

	template<class T>
	class Neigh_2d
	{
		public:
			using value_type = T;
			using size_type = std::size_t;

			Neigh_2d(){};

			template<class TVector>
			Neigh_2d(Stream<e_host> &stream, TVector &x_i, TVector &y_i, Value_type<TVector> r_neigh_i)
			{
				operator()(stream, x_i, y_i, r_neigh_i);
			}

			size_type size() const {return neigh.size();};

			vector<int>& operator[](const int i){ return neigh[i]; }

			const vector<int>& operator[](const int i) const { return neigh[i]; }

			template<class TVector>
			void operator()(Stream<e_host> &stream, TVector &x_i, TVector &y_i, Value_type<TVector> r_neigh_i)
			{
				const int n_neigh = 16;
				r_neigh = r_neigh_i;

				T lx;
				x = shift_vector(x_i, lx);
				int nx = static_cast<int>(ceil(lx/r_neigh));

				T ly;
				y = shift_vector(y_i, ly);
				int ny = static_cast<int>(ceil(ly/r_neigh));

				// reserve memory
				vector<vector<int>> neigh_grid(nx*ny);
				for(auto i=0; i<neigh_grid.size(); i++)
				{
					neigh_grid[i].reserve(n_neigh);
				}

				// set grid neighbors
				for(auto i=0; i<x.size(); i++)
				{
					int ix = static_cast<int>(floor(x[i]/r_neigh));
					int iy = static_cast<int>(floor(y[i]/r_neigh));
					neigh_grid[ix*ny+iy].push_back(i);
				};

				auto get_neighbors_sort = [&](T x_i, T y_i, T radius)->vector<int> 
				{
					auto vd_neigh = reserve_vector<vector<int>>(n_neigh);
					auto r2d_neigh = reserve_vector<vector<T>>(n_neigh);

					T radius2 = radius*radius;

					int ix_i = static_cast<int>(floor(x_i/r_neigh));
					int iy_i = static_cast<int>(floor(y_i/r_neigh));
					int ix_0 = max(0, ix_i-1);
					int ix_e = min(nx, ix_i+2);
					int iy_0 = max(0, iy_i-1);
					int iy_e = min(ny, iy_i+2);

					for(auto ix=ix_0; ix<ix_e; ix++)
					{
						for(auto iy=iy_0; iy<iy_e; iy++)
						{
							vector<int> &v = neigh_grid[ix*ny+iy];
							for(auto iv=0; iv<v.size(); iv++)
							{
								auto idx = v[iv];
								auto rx = x[idx]-x_i;
								auto ry = y[idx]-y_i;
								T r2d = rx*rx+ry*ry; 
								if(r2d<radius2)
								{
									r2d_neigh.push_back(r2d);
									vd_neigh.push_back(idx);
								}
							}
						}
					}

					if(vd_neigh.size()>1)
					{
						auto first = thrust::make_zip_iterator(thrust::make_tuple(r2d_neigh.begin(), vd_neigh.begin()));
						auto last = thrust::make_zip_iterator(thrust::make_tuple(r2d_neigh.end(), vd_neigh.end()));
						thrust::sort(first, last);
					}

					vd_neigh.shrink_to_fit();
					return vd_neigh;
				};

				// get neighbors
				neigh.resize(x.size());

				auto thr_neighbors_sort = [&](const Range &range)
				{
					for(auto i=range.ixy_0; i<range.ixy_e; i++)
					{
						neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh);
					}
				};

				stream.set_n_act_stream(x.size());
				stream.set_grid(x.size(), 1);
				stream.exec(thr_neighbors_sort);

				// get average minimum radius
				T r_min = 0;
				int cr_min = 0;
				for(auto i=0; i<x.size(); i++)
				{
					if(neigh[i].size()>1)
					{
						int idx = neigh[i][1];
						auto rx = x[idx]-x[i];
						auto ry = y[idx]-y[i];
						r_min += sqrt(rx*rx+ry*ry);
						cr_min++;
					}
				}
				r_min /= cr_min;

				T r_neigh_n = ::fmax(2*r_min, 1.25*r_neigh);
				// correct outsiders
				for(auto i=0; i<x.size(); i++)
				{
					if(neigh[i].size()==1)
					{
						neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
					}
				}

				////check minimun distance
				//auto thr_neighbors_sort_n = [&](const Range &range)
				//{
				//	for(auto i=range.ixy_0; i<range.ixy_e; i++)
				//	{
				//		T r_neigh_n = 1.5*d_min(i);
				//		if(r_neigh_n>d_max(i))
				//		{
				//			neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
				//		}
				//	}
				//};

				//stream.set_n_act_stream(x.size());
				//stream.set_grid(x.size(), 1);
				//stream.exec(thr_neighbors_sort_n);
			}

			template<class TGrid, class TVector>
			void delete_points(TVector &x_i, TVector &y_i, TGrid &grid_i)
			{
				const int n_neigh = 16;
				T r_neigh = 1.6*grid_i.dR_min();

				T lx;
				x = shift_vector(x_i, lx);
				int nx = static_cast<int>(ceil(lx/r_neigh));

				T ly;
				y = shift_vector(y_i, ly);
				int ny = static_cast<int>(ceil(ly/r_neigh));

				// reserve memory
				vector<vector<int>> neigh_grid(nx*ny);
				for(auto i=0; i<neigh_grid.size(); i++)
				{
					neigh_grid[i].reserve(n_neigh);
				}

				// set grid neighbors
				for(auto i=0; i<x.size(); i++)
				{
					int ix = static_cast<int>(floor(x[i]/r_neigh));
					int iy = static_cast<int>(floor(y[i]/r_neigh));
					neigh_grid[ix*ny+iy].push_back(i);
				};

				auto get_neighbors = [&](T x_i, T y_i, T radius)->vector<int> 
				{
					auto vd_neigh = reserve_vector<vector<int>>(n_neigh);

					T radius2 = radius*radius;

					int ix_i = static_cast<int>(floor(x_i/r_neigh));
					int iy_i = static_cast<int>(floor(y_i/r_neigh));
					int ix_0 = max(0, ix_i-1);
					int ix_e = min(nx, ix_i+2);
					int iy_0 = max(0, iy_i-1);
					int iy_e = min(ny, iy_i+2);

					for(auto ix=ix_0; ix<ix_e; ix++)
					{
						for(auto iy=iy_0; iy<iy_e; iy++)
						{
							vector<int> &v = neigh_grid[ix*ny+iy];
							for(auto iv=0; iv<v.size(); iv++)
							{
								auto idx = v[iv];
								auto rx = x[idx]-x_i;
								auto ry = y[idx]-y_i;
								T r2d = rx*rx+ry*ry; 
								if(r2d<radius2)
								{
									vd_neigh.push_back(idx);
								}
							}
						}
					}
					return vd_neigh;
				};

				// get neighbors
				vector<bool> bb(x.size(), true);
				auto x_o = reserve_vector<TVector>(x.size());
				auto y_o = reserve_vector<TVector>(y.size());

				for(auto i=0; i<x.size(); i++)
				{
					if(bb[i])
					{
						T rx = x_i[i];
						T ry = y_i[i];
						auto neigh = get_neighbors(x[i], y[i], r_neigh);
						if(neigh.size()>1)
						{
							rx = 0;
							ry = 0;
							for(auto j=0; j<neigh.size(); j++)
							{
								int idx = neigh[j];
								rx += x_i[idx];
								ry += y_i[idx];
								if(j>0)
								{
									bb[idx] = false;
								}
							}
							rx /= neigh.size();
							ry /= neigh.size();
						}
						x_o.push_back(rx);
						y_o.push_back(ry);
					}
				}
				x_i = x_o;
				y_i = y_o;
			}

			int size(const int &idx)
			{
				if(idx>=neigh.size())
				{
					return 1;
				}

				return neigh[idx].size();
			}

			T d_min(const int &idx)
			{
				T d = r_neigh;
				if(idx>=neigh.size())
				{
					return d;
				}

				if(neigh[idx].size()>1)
				{
					int i = neigh[idx][1];
					auto rx = x[i]-x[idx];
					auto ry = y[i]-y[idx];
					d = sqrt(rx*rx+ry*ry);
				}
				return d;
			}

			T d_max(const int &idx)
			{
				T d = r_neigh;
				if(idx>=neigh.size())
				{
					return d;
				}

				if(neigh[idx].size()>1)
				{
					int i = neigh[idx].back();
					auto rx = x[i]-x[idx];
					auto ry = y[i]-y[idx];
					d = sqrt(rx*rx+ry*ry);
				}
				return d;
			}

			template<class TGrid, class TVector>
			T radius_min(const int &idx, TGrid &grid, TVector &Im)
			{
				T r_min = 2.0*grid.dR_min();
				if(idx>=neigh.size())
				{
					return r_min;
				}

				auto interp2 = [](const r2d<T> &p, TGrid &grid, TVector &Im)->T
				{
					auto ix = grid.lb_index_x(p.x);
					auto iy = grid.lb_index_y(p.y);

					T f11 = Im[grid.ind_col(ix, iy)];
					T f12 = Im[grid.ind_col(ix, iy+1)];
					T f21 = Im[grid.ind_col(ix+1, iy)];
					T f22 = Im[grid.ind_col(ix+1, iy+1)];

					T x1 = grid.Rx(ix);
					T x2 = grid.Rx(ix+1);
					T y1 = grid.Ry(iy);
					T y2 = grid.Ry(iy+1);
					
					T dx1 = p.x-x1;
					T dx2 = x2-p.x;
					T dy1 = p.y-y1;
					T dy2 = y2-p.y;

					T f = (dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1))/((x2-x1)*(y2-y1));
					return f;

				};

				if(neigh[idx].size()>1)
				{
					int i = neigh[idx][1];
					r2d<T> p1(x[idx], y[idx]);
					r2d<T> p2(x[i], y[i]);
					r2d<T> p12 = r2d<T>(x[i], y[i])-p1;
					auto u = normalized(p12);

					T r = p12.module();
					T dr = 0.75*grid.dR_min();
					int nr = static_cast<int>(ceil(r/dr+0.5));
					dr = r/nr;

					T f_min = Im[grid.ixy(p1.x, p1.y)];
					for(auto ir=1; ir<nr; ir++)
					{
						T d = ir*dr;
						auto p = p1 + d*u;
						T f = interp2(p, grid, Im);
						if(f_min>f)
						{
							f_min = f;
							r_min = d;
						}
					}
					r_min = ::fmax(r_min, 0.5*r);
				}
				return r_min;
			}

			template<class TGrid, class TVector>
			TVector radius_min(Stream<e_host> &stream, TGrid &grid, TVector &Im)
			{
				TVector radius(x.size());

				auto thr_radius_min = [&](const Range &range)
				{
					for(auto idx=range.ixy_0; idx<range.ixy_e; idx++)
					{
						radius[idx] = radius_min(idx, grid, Im);
					}
				};
				stream.set_n_act_stream(x.size());
				stream.set_grid(x.size(), 1);
				stream.exec(thr_radius_min);

				return radius;
			}

			template<class TVector>
			TVector select(const int &idx, TVector &x, Value_type<TVector> x_sf=0, Value_type<TVector> x_sc=1)
			{
				TVector v;
				if(idx>=neigh.size())
				{
					return v;
				}

				vector<int> &neigh_i = neigh[idx];
				v.reserve(neigh_i.size());

				for(const auto &in:neigh_i)
				{
					v.push_back((x[in]-x_sf)/x_sc);
				}
				return v;
			}

			template<class TVector>
			Atom_Data_Sp<Value_type<TVector>> select(const int &idx, TVector &a, TVector &sigma, TVector &x, TVector &y)
			{
				Atom_Data_Sp<Value_type<TVector>> atom;
				if(idx>=neigh.size())
				{
					return atom;
				}

				vector<int> &neigh_i = neigh[idx];
				atom.reserve(neigh_i.size());

				for(const auto &in:neigh_i)
				{
					atom.push_back(x[in], y[in], a[in], sigma[in]);
				}
				return atom;
			}

			template<class TVector>
			T min_val(const int &idx, TVector &x)
			{
				T x_min = 0;
				if(idx>=neigh.size())
				{
					return x_min;
				}

				vector<int> &neigh_i = neigh[idx];

				x_min = neigh_i[0];
				for(const auto &in:neigh_i)
				{
					x_min = ::fmin(x_min, x[in]);
				}

				return x_min;
			}

			template<class TVector>
			T mean_val(const int &idx, TVector &x)
			{
				T x_mean = 0;
				if(idx>=neigh.size())
				{
					return x_mean;
				}

				vector<int> &neigh_i = neigh[idx];

				for(const auto &in:neigh_i)
				{
					x_mean += x[in];
				}
				x_mean /= neigh_i.size();

				return x_mean;
			}
		private:
			vector<T> x;
			vector<T> y;
			vector<vector<int>> neigh;
			T r_neigh;

			template<class TVector>
			vector<T> shift_vector(TVector &x_i, T &l_m)
			{
				T d_e = 1e-03;
				auto x_minmax = std::minmax_element(x_i.begin(), x_i.end()); 
				T x_min = *(x_minmax.first)-d_e;
				l_m = (*(x_minmax.second)+d_e-x_min);

				auto x = reserve_vector<vector<T>>(x_i.size());
				for(auto i=0; i<x_i.size(); i++)
				{
					x.push_back(x_i[i]-x_min);
				}
				return x;
			}
	};

} // namespace mt
#endif