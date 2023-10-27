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

#ifndef BOX_OCC_H
	#define BOX_OCC_H

	#include "fcns_cgpu_gen.h"
	#include "math_mt.h"
	#include "r_3d.h"
	#include "types.cuh"
	#include "particles.cuh"
	#include "cgpu_stream.cuh"

	namespace mt
	{
		template <class T>
		class Box_Occ_3d
		{
			public:
				Box_Occ_3d(): d2_min(), a_min(), r_0(), nx(), ny(), nz(){};

				Box_Occ_3d(T d_min, R_3d<T> bs, R_3d<T> r_0 = R_3d<T>())
				{
					set_in_data(d_min, bs, r_0);
				}

				void set_in_data(T d_min, R_3d<T> bs, R_3d<T> r_0 = R_3d<T>())
				{	
					d2_min = ::square(d_min);
					a_min = d_min/mt::c_3i2<T>;
					this->r_0 = r_0;

					nx = fcn_cceil<dt_int32>(bs.x/a_min);
					ny = fcn_cceil<dt_int32>(bs.y/a_min);
					nz = fcn_cceil<dt_int32>(bs.z/a_min);

					occ.clear();
					occ.shrink_to_fit();
					occ.resize({ny, nx, nz}, -1);
				}

				void set(R_3d<T> r, const dt_int32& val)
				{
					r = (r-r_0)/a_min;

					const auto ix = fcn_cfloor<dt_int32>(r.x);
					const auto iy = fcn_cfloor<dt_int32>(r.y);
					const auto iz = fcn_cfloor<dt_int32>(r.z);

					occ(iy, ix, iz) = val;
				}

				template <class TPtc_3d>
				dt_bool check_r_min(const TPtc_3d& atoms, R_3d<T> r)
				{
					const auto r_idx = (r-r_0)/a_min;

					const auto ix_m = fcn_cfloor<dt_int32>(r_idx.x);
					const auto iy_m = fcn_cfloor<dt_int32>(r_idx.y);
					const auto iz_m = fcn_cfloor<dt_int32>(r_idx.z);

					const auto bb_bound = fcn_chk_bound(ix_m, 0, nx) && fcn_chk_bound(iy_m, 0, ny) && fcn_chk_bound(iz_m, 0, nz);

					if (bb_bound && (occ(iy_m, ix_m, iz_m) > -1))
					{
						return false;
					}

					const auto ax = range_loop_pbc(ix_m, nx);
					const auto ay = range_loop_pbc(iy_m, ny);
					const auto az = range_loop(iz_m, nz);

					for(auto iz: az)
					{
						for(auto ix: ax)
						{
							for(auto iy: ay)
							{
								auto iatom = occ(iy, ix, iz);
								if (iatom>-1)
								{
									auto d2 = atoms.norm_2_pbc_xy(iatom, r);
									if (d2<d2_min)
									{
										return false;
									}
								}
							}
						}
					}

					return true;
				}

				R_3d<T> r_0;

				Vctr_cpu<dt_int32> occ;
			private:
				T a_min;
				T d2_min;
				dt_int32 nx;
				dt_int32 ny;
				dt_int32 nz;

				Vctr_std<dt_int32> range_loop(const dt_int32& k, const dt_int32& n_k)
				{
					if (k==0)
					{
						return {k, k+1, k+2};
					}
					else if (k==1)
					{
						return {k-1, k, k+1, k+2};
					}
					else if (k==n_k-1)
					{
						return {k-2, k-1, k};
					}
					else if (k==n_k-2)
					{
						return {k-2, k-1, k, k+1};
					}
					else
					{
						return {k-2, k-1, k, k+1, k+2};
					}
				}

				Vctr_std<dt_int32> range_loop_pbc(const dt_int32& k, const dt_int32& n_k)
				{
					if (k==0)
					{
						return {n_k-2, n_k-1, k, k+1, k+2};
					}
					else if (k==1)
					{
						return {n_k-1, k-1, k, k+1, k+2};
					}
					else if (k==n_k-1)
					{
						return {k-2, k-1, k, 0, 1};
					}
					else if (k==n_k-2)
					{
						return {k-2, k-1, k, k+1, 0};
					}
					else
					{
						return {k-2, k-1, k, k+1, k+2};
					}
				}	 
		};

		template <class T>
		class Neigh_2d
		{
			public:
				using value_type = T;
				using size_type = dt_uint64;

				Neigh_2d() {};

				template <class TVctr>
				Neigh_2d(Stream<edev_cpu>& stream, TVctr& x_i, TVctr& y_i, Value_type<TVctr> r_neigh_i)
				{
					operator()(stream, x_i, y_i, r_neigh_i);
				}

				size_type size() const {return neigh.size(); };

				Vctr_std<dt_int32>& operator[](const dt_int32 i){ return neigh[i]; }

				const Vctr_std<dt_int32>& operator[](const dt_int32 i) const { return neigh[i]; }

				template <class TVctr>
				void operator()(Stream<edev_cpu>& stream, TVctr& x_i, TVctr& y_i, Value_type<TVctr> r_neigh_i)
				{
					const dt_int32 n_neigh = 16;
					r_neigh = r_neigh_i;

					T bs_x;
					x = sft_vector(x_i, bs_x);
					dt_int32 nx = static_cast<dt_int32>(::ceil(bs_x/r_neigh));

					T bs_y;
					y = sft_vector(y_i, bs_y);
					dt_int32 ny = static_cast<dt_int32>(::ceil(bs_y/r_neigh));

					// reserve memory
					Vctr_std<Vctr_std<dt_int32>> neigh_grid(nx*ny);
					for(auto i=0; i<neigh_grid.size(); i++)
					{
						neigh_grid[i].reserve(n_neigh);
					}

					// set grid_2d neighbors
					for(auto i=0; i<x.size(); i++)
					{
						dt_int32 ix = static_cast<dt_int32>(::floor(x[i]/r_neigh));
						dt_int32 iy = static_cast<dt_int32>(::floor(y[i]/r_neigh));
						neigh_grid[ix*ny+iy].push_back(i);
					};

					auto get_neighbors_sort = [&](T x_i, T y_i, T radius)->Vctr_std<dt_int32> 
					{
						Vctr_std<dt_int32> vd_neigh;
						vd_neigh.reserve(n_neigh);

						Vctr_std<T> r2d_neigh;
						r2d_neigh.reserve(n_neigh);

						T radius2 = radius*radius;

						dt_int32 ix_i = static_cast<dt_int32>(::floor(x_i/r_neigh));
						dt_int32 iy_i = static_cast<dt_int32>(::floor(y_i/r_neigh));
						dt_int32 ix_0 = max(0, ix_i-1);
						dt_int32 ix_e = min(nx, ix_i+2);
						dt_int32 iy_0 = max(0, iy_i-1);
						dt_int32 iy_e = min(ny, iy_i+2);

						for(auto ix=ix_0; ix<ix_e; ix++)
						{
							for(auto iy=iy_0; iy<iy_e; iy++)
							{
								Vctr_std<dt_int32>& v = neigh_grid[ix*ny+iy];
								for(auto iv=0; iv<v.size(); iv++)
								{
									auto idx = v[iv];
									auto rx = x[idx]-x_i;
									auto ry = y[idx]-y_i;
									T R_2d = rx*rx+ry*ry;
									if (R_2d<radius2)
									{
										r2d_neigh.push_back(R_2d);
										vd_neigh.push_back(idx);
									}
								}
							}
						}

						if (vd_neigh.size()>1)
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

					auto thr_neighbors_sort = [&](const iThread_Rect_2d& range)
					{
						for(auto i=range.ind_0; i<range.ind_e; i++)
						{
							neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh);
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					// check out
					//stream.exec(thr_neighbors_sort);

					// get average minimum radius
					T r_min = 0;
					dt_int32 cr_min = 0;
					for(auto i=0; i<x.size(); i++)
					{
						if (neigh[i].size()>1)
						{
							dt_int32 idx = neigh[i][1];
							auto rx = x[idx]-x[i];
							auto ry = y[idx]-y[i];
							r_min += sqrt(rx*rx+ry*ry);
							cr_min++;
						}
					}
					r_min /= cr_min;

					T r_neigh_n = ::fmax(2*r_min, 1.25*r_neigh);
					// correct outliers
					for(auto i=0; i<x.size(); i++)
					{
						if (neigh[i].size() == 1)
						{
							neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
						}
					}

					// // check minimum distance
					// auto thr_neighbors_sort_n = [&](const iThread_Rect_2d& range)
					// {
					// 	for(auto i=range.ind_0; i<range.ind_e; i++)
					// 	{
					// 		T r_neigh_n = 1.5*d_min(i);
					// 		if (r_neigh_n>d_max(i))
					// 		{
					// 			neigh[i] = get_neighbors_sort(x[i], y[i], r_neigh_n);
					// 		}
					// 	}
					// };

					// stream.set_n_stream_act(x.size());
					// stream.set_grid(x.size(), 1);
					// stream.exec(thr_neighbors_sort_n);
				}

				template <class TGrid, class TVctr>
				void delete_points(TVctr& x_i, TVctr& y_i, TGrid& grid_i)
				{
					const dt_int32 n_neigh = 16;
					T r_neigh = 1.6*grid_i.dR_min();

					T bs_x;
					x = sft_vector(x_i, bs_x);
					dt_int32 nx = static_cast<dt_int32>(::ceil(bs_x/r_neigh));

					T bs_y;
					y = sft_vector(y_i, bs_y);
					dt_int32 ny = static_cast<dt_int32>(::ceil(bs_y/r_neigh));

					// reserve memory
					Vctr_std<Vctr_std<dt_int32>> neigh_grid(nx*ny);
					for(auto i=0; i<neigh_grid.size(); i++)
					{
						neigh_grid[i].reserve(n_neigh);
					}

					// set grid_2d neighbors
					for(auto i=0; i<x.size(); i++)
					{
						dt_int32 ix = static_cast<dt_int32>(::floor(x[i]/r_neigh));
						dt_int32 iy = static_cast<dt_int32>(::floor(y[i]/r_neigh));
						neigh_grid[ix*ny+iy].push_back(i);
					};

					auto get_neighbors = [&](T x_i, T y_i, T radius)->Vctr_std<dt_int32> 
					{
						Vctr_std<dt_int32> vd_neigh;
						vd_neigh.reserve(n_neigh);

						T radius2 = radius*radius;

						dt_int32 ix_i = static_cast<dt_int32>(::floor(x_i/r_neigh));
						dt_int32 iy_i = static_cast<dt_int32>(::floor(y_i/r_neigh));
						dt_int32 ix_0 = max(0, ix_i-1);
						dt_int32 ix_e = min(nx, ix_i+2);
						dt_int32 iy_0 = max(0, iy_i-1);
						dt_int32 iy_e = min(ny, iy_i+2);

						for(auto ix=ix_0; ix<ix_e; ix++)
						{
							for(auto iy=iy_0; iy<iy_e; iy++)
							{
								Vctr_std<dt_int32>& v = neigh_grid[ix*ny+iy];
								for(auto iv=0; iv<v.size(); iv++)
								{
									auto idx = v[iv];
									auto rx = x[idx]-x_i;
									auto ry = y[idx]-y_i;
									T R_2d = rx*rx+ry*ry;
									if (R_2d<radius2)
									{
										vd_neigh.push_back(idx);
									}
								}
							}
						}
						return vd_neigh;
					};

					// get neighbors
					Vctr_std<dt_bool> bb(x.size(), true);

					TVctr x_o;
					x_o.reserve(x.size());

					TVctr y_o;
					y_o.reserve(y.size());

					for(auto i=0; i<x.size(); i++)
					{
						if (bb[i])
						{
							T rx = x_i[i];
							T ry = y_i[i];
							auto neigh = get_neighbors(x[i], y[i], r_neigh);
							if (neigh.size()>1)
							{
								rx = 0;
								ry = 0;
								for(auto j=0; j<neigh.size(); j++)
								{
									dt_int32 idx = neigh[j];
									rx += x_i[idx];
									ry += y_i[idx];
									if (j>0)
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

				dt_int32 size(const dt_int32& idx)
				{
					if (idx>=neigh.size())
					{
						return 1;
					}

					return neigh[idx].size();
				}

				T d_min(const dt_int32& idx)
				{
					T d = r_neigh;
					if (idx>=neigh.size())
					{
						return d;
					}

					if (neigh[idx].size()>1)
					{
						dt_int32 i = neigh[idx][1];
						auto rx = x[i]-x[idx];
						auto ry = y[i]-y[idx];
						d = sqrt(rx*rx+ry*ry);
					}
					return d;
				}

				T d_max(const dt_int32& idx)
				{
					T d = r_neigh;
					if (idx>=neigh.size())
					{
						return d;
					}

					if (neigh[idx].size()>1)
					{
						dt_int32 i = neigh[idx].back();
						auto rx = x[i]-x[idx];
						auto ry = y[i]-y[idx];
						d = sqrt(rx*rx+ry*ry);
					}
					return d;
				}

				template <class TGrid, class TVctr>
				T radius_min(const dt_int32& idx, TGrid& grid_2d, TVctr& Im)
				{
					T r_min = 2.0*grid_2d.dr_min();
					if (idx>=neigh.size())
					{
						return r_min;
					}

					auto interp2 = [](const R_2d<T>& p, TGrid& grid_2d, TVctr& Im)->T
					{
						auto ix = grid_2d.rx_2_irx_bfds(p.x);
						auto iy = grid_2d.ry_2_iry_bfds(p.y);

						T f11 = Im[grid_2d.sub_2_ind(ix, iy)];
						T f12 = Im[grid_2d.sub_2_ind(ix, iy+1)];
						T f21 = Im[grid_2d.sub_2_ind(ix+1, iy)];
						T f22 = Im[grid_2d.sub_2_ind(ix+1, iy+1)];

						T x1 = grid_2d.rx(ix);
						T x2 = grid_2d.rx(ix+1);
						T y1 = grid_2d.ry(iy);
						T y2 = grid_2d.ry(iy+1);
					
						T dx1 = p.x-x1;
						T dx2 = x2-p.x;
						T dy1 = p.y-y1;
						T dy2 = y2-p.y;

						T f = (dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1))/((x2-x1)*(y2-y1));
						return f;

					};

					if (neigh[idx].size()>1)
					{
						dt_int32 i = neigh[idx][1];
						R_2d<T> p1(x[idx], y[idx]);
						R_2d<T> p2(x[i], y[i]);
						R_2d<T> p12 = R_2d<T>(x[i], y[i])-p1;
						auto u = normalize(p12);

						T r = p12.norm();
						T dr = 0.75*grid_2d.dr_min();
						dt_int32 nr = static_cast<dt_int32>(::ceil(r/dr+0.5));
						dr = r/nr;

						T f_min = Im[grid_2d.rv_2_ir_bfds(p1.x, p1.y)];
						for(auto ir=1; ir<nr; ir++)
						{
							T d = ir*dr;
							auto p = p1 + d*u;
							T f = interp2(p, grid_2d, Im);
							if (f_min>f)
							{
								f_min = f;
								r_min = d;
							}
						}
						r_min = ::fmax(r_min, 0.5*r);
					}
					return r_min;
				}

				template <class TGrid, class TVctr>
				TVctr radius_min(Stream<edev_cpu>& stream, TGrid& grid_2d, TVctr& Im)
				{
					TVctr radius(x.size());

					auto thr_radius_min = [&](const iThread_Rect_2d& range)
					{
						for(auto idx=range.ind_0; idx<range.ind_e; idx++)
						{
							radius[idx] = radius_min(idx, grid_2d, Im);
						}
					};
					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					// check out
					//stream.exec(thr_radius_min);

					return radius;
				}

				template <class TVctr>
				TVctr select(const dt_int32& idx, TVctr& x, Value_type<TVctr> x_sf=0, Value_type<TVctr> x_sc=1)
				{
					TVctr v;
					if (idx>=neigh.size())
					{
						return v;
					}

					Vctr_std<dt_int32>& neigh_i = neigh[idx];
					v.reserve(neigh_i.size());

					for(const auto &in:neigh_i)
					{
						v.push_back((x[in]-x_sf)/x_sc);
					}
					return v;
				}

				//template <class TVctr>
				//Ptc_gauss_2d<Value_type<TVctr>> select(const dt_int32& idx, TVctr& a, TVctr& sigma, TVctr& x, TVctr& y)
				//{
				//	Ptc_gauss_2d<Value_type<TVctr>> atom;
				//	if (idx>=neigh.size())
				//	{
				//		return atom;
				//	}

				//	Vctr_std<dt_int32>& neigh_i = neigh[idx];
				//	atom.reserve(neigh_i.size());

				//	for(const auto &in:neigh_i)
				//	{
				//		atom.push_back(x[in], y[in], a[in], sigma[in]);
				//	}
				//	return atom;
				//}

				template <class TVctr>
				T min_val(const dt_int32& idx, TVctr& x)
				{
					T x_min = 0;
					if (idx>=neigh.size())
					{
						return x_min;
					}

					Vctr_std<dt_int32>& neigh_i = neigh[idx];

					x_min = neigh_i[0];
					for(const auto &in:neigh_i)
					{
						x_min = ::fmin(x_min, x[in]);
					}

					return x_min;
				}

				template <class TVctr>
				T mean_val(const dt_int32& idx, TVctr& x)
				{
					T x_mean = 0;
					if (idx>=neigh.size())
					{
						return x_mean;
					}

					Vctr_std<dt_int32>& neigh_i = neigh[idx];

					for(const auto &in:neigh_i)
					{
						x_mean += x[in];
					}
					x_mean /= neigh_i.size();

					return x_mean;
				}
			private:
				Vctr_std<T> x;
				Vctr_std<T> y;
				Vctr_std<Vctr_std<dt_int32>> neigh;
				T r_neigh;

				template <class TVctr>
				Vctr_std<T> sft_vector(TVctr& x_i, T& l_m)
				{
					T d_e = 1e-03;
					auto x_minmax = fcn_minmax_element(x_i.begin(), x_i.end());
					T x_min = *(x_minmax.first)-d_e;
					l_m = (*(x_minmax.second)+d_e-x_min);

					Vctr_std<T> x;
					x.reserve(x_i.size());

					for(auto i=0; i<x_i.size(); i++)
					{
						x.push_back(x_i[i]-x_min);
					}
					return x;
				}
		};

	}
#endif