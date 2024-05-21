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

#ifndef BOX_OCC_H
#define BOX_OCC_H

#include <vector>
#include <deque>

#include "types.cuh"
#include "math.cuh"
#include "lin_alg_def.cuh"
#include "atomic_data_mt.hpp"
#include "stream.cuh"

namespace mt
{
	template <class T>
	class Box_Occ
	{
		public:
			Box_Occ(): d_min(0), d2_min(0), a_min(0), l_x(0), l_y(0), l_z(0), nx(0), ny(0), nz(0), nxy(0){};

			void init()
			{
				thrust::fill(occ.begin(), occ.end(), T(-1));
			}

			void set_input_data(T d_min_i, T lx_i, T ly_i, T lz_i)
			{	
				d_min = d_min_i;
				d2_min = pow(d_min, 2);
				a_min = d_min/mt::c_3i2;
				l_x = lx_i;
				l_y = ly_i;
				l_z = lz_i;
				nx = static_cast<int64_t>(ceil(l_x/a_min));
				ny = static_cast<int64_t>(ceil(l_y/a_min));
				nz = static_cast<int64_t>(ceil(l_z/a_min));
				nxy = nx*ny;
				occ.clear();
				occ.resize(nxy*nz, -1);
			}

			std::size_t xyz_2_ind(const int &ix, const int &iy, const int &iz) const
			{
				return (int64_t(iz)*nxy+ int64_t(iy)*nx + int64_t(ix));
			}

			int get_occ(const int &ix, const int &iy, const int &iz) const
			{
				return occ[xyz_2_ind(ix, iy, iz)];
			}

			void set_occ(const int &ix, const int &iy, const int &iz, const int &val)
			{
				occ[xyz_2_ind(ix, iy, iz)] = val;
			}

			void set_occ(const r3d<T> &r, const int &iatoms)
			{
				const int ix = static_cast<int>(floor(r.x/a_min));
				const int iy = static_cast<int>(floor(r.y/a_min));
				const int iz = static_cast<int>(floor(r.z/a_min));
				set_occ(ix, iy, iz, iatoms);
			}

			std::vector<int> range_loop(const int &k, const int &n_k)
			{
				if(k==0)
				{
					return std::vector<int>{k, k+1, k+2};
				}
				else if(k==1)
				{
					return std::vector<int>{k-1, k, k+1, k+2};
				}
				else if(k==n_k-1)
				{
					return std::vector<int>{k-2, k-1, k};
				}
				else if(k==n_k-2)
				{
					return std::vector<int>{k-2, k-1, k, k+1};
				}
				else
				{
					return std::vector<int>{k-2, k-1, k, k+1, k+2};
				}
			}

			// std::vector<int> range_loop_pbc(const int &k, const int &n_k) {
			// 	return {
			// 		(k - 2 + n_k) % n_k,
			// 		(k - 1 + n_k) % n_k,
			// 		k,
			// 		(k + 1) % n_k,
			// 		(k + 2) % n_k
			// 	};
			// }

			std::vector<int> range_loop_pbc(const int &k, const int &n_k) {
				return {
					(k - 1 + n_k) % n_k,
					k,
					(k + 1) % n_k,
				};
			}

			// std::vector<int> range_loop_pbc(const int &k, const int &n_k)
			// {
			// 	if(k==0)
			// 	{
			// 		return std::vector<int>{n_k-2, n_k-1, k, k+1, k+2};
			// 	}
			// 	else if(k==1)
			// 	{
			// 		return std::vector<int>{n_k-1, k-1, k, k+1, k+2};
			// 	}
			// 	else if(k==n_k-1)
			// 	{
			// 		return std::vector<int>{k-2, k-1, k, 0, 1};
			// 	}
			// 	else if(k==n_k-2)
			// 	{
			// 		return std::vector<int>{k-2, k-1, k, k+1, 0};
			// 	}
			// 	else
			// 	{
			// 		return std::vector<int>{k-2, k-1, k, k+1, k+2};
			// 	}
			// }	 

			template <class TAtom>
			bool check_r_min(TAtom &atoms, const r3d<T> &r, T &d2)
			{
				const int ix_i = static_cast<int>(floor(r.x/a_min));
				const int iy_i = static_cast<int>(floor(r.y/a_min));
				const int iz_i = static_cast<int>(floor(r.z/a_min));

				if (get_occ(ix_i, iy_i, iz_i)>-1)
				{
					return false;
				}

				auto ax = range_loop_pbc(ix_i, nx);
				auto ay = range_loop_pbc(iy_i, ny);
				auto az = range_loop(iz_i, nz);

				for(auto iz: az)
				{
					for(auto iy: ay)
					{
						for(auto ix: ax)
						{
							auto iatoms = get_occ(ix, iy, iz);
							if(iatoms>-1)
							{
								d2 = atoms.norm_pbc_xy(iatoms, r);
								if(d2<d2_min)
								{
									return false;
								}
							}
						}
					}
				}

				return true;
			}

			T d_min;
			T d2_min;
			T a_min;
			T l_x;
			T l_y;
			T l_z;
			int64_t nx;
			int64_t ny;
			int64_t nz;
			int64_t nxy;
			Vector<int, e_host> occ;
	};

	template <class T>
	class Neigh_2d
	{
		public:
			using value_type = T;
			using size_type = std::size_t;

			Neigh_2d(){};

			template <class TVector>
			Neigh_2d(Stream<e_host> &stream, TVector &x_i, TVector &y_i, Value_type<TVector> r_neigh_i)
			{
				operator()(stream, x_i, y_i, r_neigh_i);
			}

			size_type size() const {return neigh.size();};

			vector<int>& operator[](const int i){ return neigh[i]; }

			const vector<int>& operator[](const int i) const { return neigh[i]; }

			template <class TVector>
			void operator()(Stream<e_host> &stream, TVector &x_i, TVector &y_i, Value_type<TVector> r_neigh_i)
			{
				const int n_neigh = 16;
				r_neigh = r_neigh_i;

				T lx;
				x = sft_vector(x_i, lx);
				int nx = static_cast<int>(ceil(lx/r_neigh));

				T ly;
				y = sft_vector(y_i, ly);
				int ny = static_cast<int>(ceil(ly/r_neigh));

				// reserve memory
				vector<vector<int>> neigh_grid(nx*ny);
				for(auto i=0; i<neigh_grid.size(); i++)
				{
					neigh_grid[i].reserve(n_neigh);
				}

				// set grid_2d neighbors
				for(auto i=0; i<x.size(); i++)
				{
					int ix = static_cast<int>(floor(x[i]/r_neigh));
					int iy = static_cast<int>(floor(y[i]/r_neigh));
					neigh_grid[ix*ny+iy].push_back(i);
				};

				auto get_neighbors_sort = [&](T x_i, T y_i, T radius)->vector<int> 
				{
					vector<int> vd_neigh;
					vd_neigh.reserve(n_neigh);

					vector<T> r2d_neigh;
					r2d_neigh.reserve(n_neigh);

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

				auto thr_neighbors_sort = [&](const Range_2d &range)
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
					if(neigh[i].size() == 1)
					{
						neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
					}
				}

				// // check minimun distance
				// auto thr_neighbors_sort_n = [&](const Range_2d &range)
				// {
				// 	for(auto i=range.ixy_0; i<range.ixy_e; i++)
				// 	{
				// 		T r_neigh_n = 1.5*d_min(i);
				// 		if(r_neigh_n>d_max(i))
				// 		{
				// 			neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
				// 		}
				// 	}
				// };

				// stream.set_n_act_stream(x.size());
				// stream.set_grid(x.size(), 1);
				// stream.exec(thr_neighbors_sort_n);
			}

			template <class TGrid, class TVector>
			void delete_points(TVector &x_i, TVector &y_i, TGrid &grid_i)
			{
				const int n_neigh = 16;
				T r_neigh = 1.6*grid_i.dR_min();

				T lx;
				x = sft_vector(x_i, lx);
				int nx = static_cast<int>(ceil(lx/r_neigh));

				T ly;
				y = sft_vector(y_i, ly);
				int ny = static_cast<int>(ceil(ly/r_neigh));

				// reserve memory
				vector<vector<int>> neigh_grid(nx*ny);
				for(auto i=0; i<neigh_grid.size(); i++)
				{
					neigh_grid[i].reserve(n_neigh);
				}

				// set grid_2d neighbors
				for(auto i=0; i<x.size(); i++)
				{
					int ix = static_cast<int>(floor(x[i]/r_neigh));
					int iy = static_cast<int>(floor(y[i]/r_neigh));
					neigh_grid[ix*ny+iy].push_back(i);
				};

				auto get_neighbors = [&](T x_i, T y_i, T radius)->vector<int> 
				{
					vector<int> vd_neigh;
					vd_neigh.reserve(n_neigh);

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

				TVector x_o;
				x_o.reserve(x.size());

				TVector y_o;
				y_o.reserve(y.size());

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

			template <class TGrid, class TVector>
			T radius_min(const int &idx, TGrid &grid_2d, TVector &Im)
			{
				T r_min = 2.0*grid_2d.dR_min();
				if(idx>=neigh.size())
				{
					return r_min;
				}

				auto interp2 = [](const r2d<T> &p, TGrid &grid_2d, TVector &Im)->T
				{
					auto ix = grid_2d.lb_index_x(p.x);
					auto iy = grid_2d.lb_index_y(p.y);

					T f11 = Im[grid_2d.ind_col(ix, iy)];
					T f12 = Im[grid_2d.ind_col(ix, iy+1)];
					T f21 = Im[grid_2d.ind_col(ix+1, iy)];
					T f22 = Im[grid_2d.ind_col(ix+1, iy+1)];

					T x1 = grid_2d.Rx(ix);
					T x2 = grid_2d.Rx(ix+1);
					T y1 = grid_2d.Ry(iy);
					T y2 = grid_2d.Ry(iy+1);
					
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
					T dr = 0.75*grid_2d.dR_min();
					int nr = static_cast<int>(ceil(r/dr+0.5));
					dr = r/nr;

					T f_min = Im[grid_2d.ixy(p1.x, p1.y)];
					for(auto ir=1; ir<nr; ir++)
					{
						T d = ir*dr;
						auto p = p1 + d*u;
						T f = interp2(p, grid_2d, Im);
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

			template <class TGrid, class TVector>
			TVector radius_min(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im)
			{
				TVector radius(x.size());

				auto thr_radius_min = [&](const Range_2d &range)
				{
					for(auto idx=range.ixy_0; idx<range.ixy_e; idx++)
					{
						radius[idx] = radius_min(idx, grid_2d, Im);
					}
				};
				stream.set_n_act_stream(x.size());
				stream.set_grid(x.size(), 1);
				stream.exec(thr_radius_min);

				return radius;
			}

			template <class TVector>
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

			template <class TVector>
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

			template <class TVector>
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

			template <class TVector>
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

			template <class TVector>
			vector<T> sft_vector(TVector &x_i, T &l_m)
			{
				T d_e = 1e-03;
				auto x_minmax = std::minmax_element(x_i.begin(), x_i.end()); 
				T x_min = *(x_minmax.first)-d_e;
				l_m = (*(x_minmax.second)+d_e-x_min);

				vector<T> x;
				x.reserve(x_i.size());

				for(auto i=0; i<x_i.size(); i++)
				{
					x.push_back(x_i[i]-x_min);
				}
				return x;
			}
	};

} // namespace mt
#endif