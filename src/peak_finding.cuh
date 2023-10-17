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

#ifndef PEAK_FINDING_H
	#define PEAK_FINDING_H

	#include <vector>
	#include <deque>

	#include <fftw3.h>
	#include "math_mt.h"
	#include "types.cuh"
	#include "type_traits_gen.h"
	#include "cgpu_fft.cuh"
	#include "cgpu_stream.cuh"
	#include "lapack.hpp"
	#include "box_occ.hpp"
	#include "cgpu_classes.cuh"
	#include "cpu_fcns_image.hpp"

	namespace mt
	{
		template <class T, eDev Dev>
		class Peak_Finding
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;
				using size_type = dt_uint64;

				Peak_Finding(): sigma(1.0), thres(0.5), niter(10), d_error(1e-3), 
				image_mean(0), image_std(1), radius_n(6*sigma), ref_mean(500), ref_std(80) {};

				Peak_Finding(Grid_2d<T>& grid_i, TVctr_r &image_i, T sigma_i, T thres_i, dt_int32 niter_i)
				{
					operator()(grid_i, image_i, sigma_i, thres_i, niter_i);
				}

				void cleanup()
				{
					fft_2d.cleanup();
				}

				void operator()(Grid_2d<T>& grid_i, TVctr_r &image_i, T sigma_i, T thres_i, dt_int32 niter_i)
				{
					grid_2d = grid_i;
					sigma = sigma_i;
					thres = thres_i;
					niter = niter_i;
					d_error = 1e-4;

					ref_mean = 500;
					ref_std = 80;
					radius_n = 6*sigma;

					stream.resize(4);
					fft_2d.create_plan_2d(grid_2d.ny, grid_2d.nx, stream.size());

					// image_mean and image_std are defined here
					image = scale_in_image(stream, image_i);
				}

				void find()
				{
					// denoise im
					dt_int32 nkr_w = max(1, static_cast<dt_int32>(::floor(sigma/2+0.5)));
					dt_int32 nkr_m = max(1, static_cast<dt_int32>(::floor(sigma/3+0.5)));
					image_den = fltr_poiss_nois_2d(stream, grid_2d, image, nkr_w, nkr_m);
		
					// get peak signal to noise ratio
					auto PSNR = fcn_get_psnr(stream, image, image_den);

					// deconvolute in image
					auto image_dcon = image;
					Gauss_Dcv_2d<T, Dev> gauss_dcv_2d(&stream, &fft_2d, grid_2d);
					gauss_dcv_2d(sigma, PSNR, image_dcon);

					thrust::for_each(image_dcon.begin(), image_dcon.end(), [&](T& v){ v = (v<0)?0:v; });

					// get background
					dt_int32 nkr = static_cast<dt_int32>(::round(5.0*sigma/grid_2d.dR_min()));
					auto image_thr = fcn_morp_op_open(stream, grid_2d.ny, grid_2d.nx, image_dcon, nkr);

					Gauss_Cv_2d<T, Dev> gauss_cv_2d(&stream, &fft_2d, grid_2d);
					gauss_cv_2d(sigma, image_thr);

					// background substraction
					thrust::transform(image_dcon.begin(), image_dcon.end(), image_thr.begin(), 
					image_thr.begin(), [](const T& a, const T&b){ return ::fmax(a-b, T(0)); });

					// fast peak finding
					fast_peak_finding(stream, image_thr, x, y);

					// deleting close neighbors
					neigh.delete_points(x, y, grid_2d);

					// set output data
					A.resize(x.size());
					std::fill(A.begin(), A.end(), T(0));

					A_min.resize(x.size());
					std::fill(A_min.begin(), A_min.end(), T(0));

					A_max.resize(x.size());
					std::fill(A_max.begin(), A_max.end(), T(0));

					S.resize(x.size());
					std::fill(S.begin(), S.end(), T(0));

					S_min.resize(x.size());
					std::fill(S_min.begin(), S_min.end(), T(0));

					S_max.resize(x.size());
					std::fill(S_max.begin(), S_max.end(), T(0));

					bg = 0;

					dt_int32 ibor = max(10, grid_2d.r_2_ir_cds_dr_min(3*sigma));
					range_ct.ix_0 = ibor;
					range_ct.ix_e = grid_2d.nx-ibor;
					range_ct.iy_0 = ibor;
					range_ct.iy_e = grid_2d.ny-ibor;
				}

				void fit()
				{
					// get neighbors
					neigh(stream, x, y, radius_n);

					// get image mask
					mask = get_mask(stream, neigh, x, y);

					// set bg values
					T bg_min, bg_max;
					set_init_bg_values(stream, mask, image_den, x, y, bg, bg_min, bg_max);

					// set A values
					set_init_A_values(stream, image_den, neigh, x, y, A, A_min, A_max, bg);

					// get first approximations
					S = fit_b_0(stream, mask, image, x, y, A, sigma, bg, 15, d_error);
					A = fit_a_0(stream, mask, image, x, y, A, A_min, A_max, S, bg);
				
					for(auto ipk = 0; ipk < x.size(); ipk++)
					{
						S_min[ipk] = ::fmax(0.5*grid_2d.dR_min(), 0.25*S[ipk]);
						S_max[ipk] = ::fmin(0.9*neigh.d_min(ipk), 1.5*S[ipk]);
					}
				
					region.set(stream, grid_2d, image, neigh, x, y, 1, 1);
					for(auto it = 0; it < 1; it++)
					{
	 				fit_a_b_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 20, d_error);
					}

					// region.set(stream, grid_2d, image, neigh, x, y, 1, 1.0);
					// for(auto it = 0; it < 3; it++)
					// {
					// 	fit_a_b_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 20, d_error);
					// }
					// bg = fit_c_0(stream, image, x, y, A, S, bg, bg_min, bg_max);
			
					for(auto it = 0; it < niter; it++)
					{
	 					region.set(stream, grid_2d, image, neigh, x, y, 1, 0.4);
	 					fit_x_y_0(stream, region, x, y, A, S, bg, 10, d_error);

						region.set(stream, grid_2d, image, neigh, x, y, 1, 1);
						fit_a_b_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 20, d_error);
					}

					// region.set(stream, grid_2d, image, neigh, x, y, 1, 1.0);
					// fit_a_b_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 10, d_error);
					// fit_a_b_x_y_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 10, d_error);
					// for(auto it = 0; it < 3; it++)
					// {
					// 	region.set(stream, grid_2d, image, neigh, x, y, 1, 1.0);
					// 	fit_a_b_x_y_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 10, d_error);
					// }
					// neigh(stream, x, y, radius_n);
					// region.set(stream, grid_2d, image, neigh, x, y, 1, 1.0);
					// fit_a_b_x_y(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 10, d_error);
					// fit_a_b_x_y_0(stream, region, x, y, A, A_min, A_max, S, S_min, S_max, bg, 20, d_error);
				}

				void get_data(TVctr_r &x_o, TVctr_r &y_o, TVctr_r &A_o, TVctr_r &S_o, T& bg_o)
				{
					dt_int32 npeaks = x.size();

					x_o.resize(npeaks);
					y_o.resize(npeaks);
					A_o.resize(npeaks);
					S_o.resize(npeaks);

					T sf = ref_std/image_std;
					for(auto idx=0; idx<x.size(); idx++)
					{
						x_o[idx] = x[idx];
						y_o[idx] = y[idx];
						A_o[idx] = A[idx]*sf;
						S_o[idx] = S[idx];
					}

					bg_o = (bg-image_mean)*sf + ref_mean;
				}

			private:
				struct Region_Rad_2d
				{
					public:
						using size_type = dt_uint64;

						TVctr_r Rx;
						TVctr_r Ry;
						TVctr_r R2;
						TVctr_r Ixy;

						T R_max;
						T Rx_sf;
						T Ry_sf;
						T Rxy_sc;

						T Ixy_sf;
						T Ixy_sc;

						Region_Rad_2d(): Rx_sf(0), Ry_sf(0), Rxy_sc(1), Ixy_sf(0), Ixy_sc(1)
						{
						}

						void clear()
						{
							Rx.clear();
							Ry.clear();
							R2.clear();
							Ixy.clear();
						}

						void reserve(const size_type& new_size)
						{
							Rx.reserve(new_size);
							Ry.reserve(new_size);
							R2.reserve(new_size);
							Ixy.reserve(new_size);
						}

						void shrink_to_fit()
						{
							Rx.shrink_to_fit();
							Ry.shrink_to_fit();
							R2.shrink_to_fit();
							Ixy.shrink_to_fit();
						}

						TVctr_r sft_Ixy(T bg)
						{
							TVctr_r Ixy_s;
							Ixy_s.reserve(Ixy.size());

							for(auto ixy=0; ixy<Ixy.size(); ixy++)
							{
								Ixy_s.push_back(Ixy[ixy]-bg);
							}
							return Ixy_s;
						}

						TVctr_r sub_region_to_Ixy(Grid_2d<T>& grid_2d, TVctr_r &Im_s, T x, T y)
						{
							TVctr_r v = Ixy;

							T R2_max = pow(R_max, 2);

							R_2d<T> p(x, y);
							auto range = grid_2d.region_ind(p, R_max);
							dt_int32 iv = 0;
							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T r2 = grid_2d.r2(ix, iy, p.x, p.y);
									if (r2 < R2_max)
									{
										v[iv++] -= Im_s[grid_2d.sub_2_ind(ix, iy)]/Ixy_sc;
									}
								}
							}
							return v;
						}

						TVctr_r sft_Ixy(Grid_2d<T>& grid_2d, TVctr_r &Im_s, T x, T y, T a, T s)
						{
							TVctr_r v = sub_region_to_Ixy(grid_2d, Im_s, x*Rxy_sc+Rx_sf, y*Rxy_sc+Ry_sf);

							T alpha = 0.5/pow(s, 2);
							T r2_l = pow(4.0*s, 2);
							for(auto im = 0; im < v.size(); im++)
							{
								T rx = Rx[im]-x;
								T ry = Ry[im]-y;
								T r2 = rx*rx+ry*ry;
								if (r2<r2_l)
								{
									v[im] += a*exp(-alpha*r2);
								}
							}

							return v;
						}

						TVctr_r sft_Ixy(Grid_2d<T>& grid_2d, TVctr_r &Im_s, TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &S)
						{
							TVctr_r v = sub_region_to_Ixy(grid_2d, Im_s, x[0]*Rxy_sc+Rx_sf, y[0]*Rxy_sc+Ry_sf);

							for(auto ip = 0; ip < x.size(); ip++)
							{
								T a = A[ip];
								T b = S[ip];
								T alpha = 0.5/pow(b, 2);
								T r2_l = pow(4.0*b, 2);
								for(auto im = 0; im < v.size(); im++)
								{
									T rx = Rx[im]-x[ip];
									T ry = Ry[im]-y[ip];
									T r2 = rx*rx+ry*ry;
									if (r2<r2_l)
									{
										v[im] += a*exp(-alpha*r2);
									}
								}
							}

							return v;
						}

						T sft_Rx(T x) const 
						{ 
							return (x-Rx_sf)/Rxy_sc;
						}

						T sft_Ry(T y) const 
						{ 
							return (y-Ry_sf)/Rxy_sc;
						}

						R_2d<T> sft_x_y(T x, T y) const 
						{ 
							x = sft_Rx(x);
							y = sft_Ry(y);
							return R_2d<T>(x, y);
						}
				};

				struct Regions
				{
					public:
						using size_type = dt_uint64;

						dt_int32 m_max;
						dt_int32 n_max;
						Regions():m_max(1), n_max(1) {}

						size_type size() const
						{
							return reg.size();
						}

						void set(Stream<edev_cpu>& stream, const Grid_2d<T>& grid_i, TVctr_r &Im_i, 
						Neigh_2d<T>& neigh, TVctr_r &x, TVctr_r &y, dt_int32 iradius=1, T ff=0.5)
						{
							grid_2d = grid_i;
							reg.resize(x.size());
							T R_min = 3.1*grid_2d.dR_min();

							auto sel_radius = [&neigh, R_min, ff](dt_int32 ipk, dt_int32 iradius)->T
							{
								T radius = 0;
								switch (iradius)
								{
									case 1:
										radius = ::fmax(R_min, ff*neigh.d_min(ipk));
										break;
									case 2:
										radius = neigh.d_max(ipk) + ff*neigh.d_min(ipk);
										break;
								}
								return radius;
							};

							auto thr_select_cir_reg = [&](const iThread_Rect_2d& range)
							{
								for(auto ipk = range.ind_0; ipk < range.ind_e; ipk++)
								{
									R_2d<T> p(x[ipk], y[ipk]);
									T radius = sel_radius(ipk, iradius);
									reg[ipk] = select_cir_reg(Im_i, p, radius);
								}
							};

							stream.set_n_stream_act(x.size());
							stream.set_grid(x.size(), 1);
							stream.exec(thr_select_cir_reg);

							T d_max = 0;
							n_max = 1;
							for(auto ipk = 0; ipk < x.size(); ipk++)
							{
								T radius = sel_radius(ipk, iradius);
								d_max = max(d_max, radius);
								n_max = max(n_max, neigh.size(ipk));
							}
							d_max += 3*grid_i.dR_min();
							m_max = static_cast<dt_int32>(std::round(c_pi<T>*pow(d_max, 2)/pow(grid_i.dR_min(), 2)));
						}

						Region_Rad_2d& operator[](const dt_int32 i){ return reg[i]; }

						const Region_Rad_2d& operator[](const dt_int32 i) const { return reg[i]; }
					private:
						vector<Region_Rad_2d> reg;
						Grid_2d<T> grid_2d;

						Region_Rad_2d select_cir_reg(TVctr_r &Im, R_2d<T> p, T radius)
						{
							T R_max = radius;
							T R2_max = pow(R_max, 2);

							auto range = grid_2d.region_ind(p, R_max);

							Region_Rad_2d region;
							region.clear();
							region.reserve(range.ind_e);

							// select circular region
							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T r2_d = grid_2d.r2(ix, iy, p.x, p.y);
									dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
									if (r2_d < R2_max)
									{
										region.Rx.push_back(grid_2d.rx(ix));
										region.Ry.push_back(grid_2d.ry(iy));
										region.R2.push_back(r2_d);
										region.Ixy.push_back(Im[ixy]);
									}
								}
							}

							region.R_max = R_max;
							region.Rx_sf = p.x;
							region.Ry_sf = p.y;
							region.Rxy_sc = R_max;

							// really important that Ixy_mean is zero if not
							// we have to change our reference system to fit the parameters
							region.Ixy_sf = 0;
							region.Ixy_sc = sqrt(fcn_variance(region.Ixy));
							// shift and fcn_scale
							dt_int32 m = region.Ixy.size();
							for(auto ixy = 0; ixy < m; ixy++)
							{
								region.Rx[ixy] = (region.Rx[ixy]-region.Rx_sf)/region.Rxy_sc;
								region.Ry[ixy] = (region.Ry[ixy]-region.Ry_sf)/region.Rxy_sc;
								region.R2[ixy] = region.R2[ixy]/pow(region.Rxy_sc, 2);
								region.Ixy[ixy] = (region.Ixy[ixy]-region.Ixy_sf)/region.Ixy_sc;
							}

							region.shrink_to_fit();
							return region;
						}
				};

				TVctr_r scale_in_image(Stream<edev_cpu>& stream, TVctr_r &image_i)
				{
					fcn_mean_var(stream, image_i, image_mean, image_std);
					image_std = sqrt(image_std);

					auto thr_scale = [=](const iThread_Rect_2d& range, TVctr_r &image_i, TVctr_r &image_o)
					{
						for(auto ixy = range.ind_0; ixy < range.ind_e; ixy++)
						{
							T v = image_i[ixy];
							v = (v-image_mean)/image_std;
							v = v*ref_std + ref_mean;
							image_o[ixy] = ::fmax(v, 4);
						}
					};

					TVctr_r image_o(image_i.size());
					stream.set_n_stream_act(grid_2d.nx);
					stream.set_grid(grid_2d.nx, grid_2d.ny);
					stream.exec(thr_scale, image_i, image_o);

					return image_o;
				}

				void fast_peak_finding(Stream<edev_cpu>& stream, TVctr_r &image, TVctr_r &x, TVctr_r &y)
				{
					auto Im_minmax = fcn_minmax_element(image.begin(), image.end());

					T thres_n = *(Im_minmax.first) + thres*(*(Im_minmax.second)-*(Im_minmax.first));

					image = fcn_threshold_max(stream, image, thres_n);

					// local maximum
					auto krn_maximum = [thres_n](const dt_int32& ix, const dt_int32& iy, Grid_2d<T>& grid_2d, TVctr_r &Im, R_2d<T>& peak)->dt_bool
					{
						auto v = Im[grid_2d.sub_2_ind(ix, iy)];
						peak = R_2d<T>(grid_2d.rx(ix), grid_2d.ry(iy));

						if (v <= thres_n)
						{
								return false;
						}

						T v1 = Im[grid_2d.sub_2_ind(ix-1, iy-1)];
						T v2 = Im[grid_2d.sub_2_ind(ix, iy-1)];
						T v3 = Im[grid_2d.sub_2_ind(ix+1, iy-1)];

						T v4 = Im[grid_2d.sub_2_ind(ix-1, iy)];
						T v6 = Im[grid_2d.sub_2_ind(ix+1, iy)];

						T v7 = Im[grid_2d.sub_2_ind(ix-1, iy+1)];
						T v8 = Im[grid_2d.sub_2_ind(ix, iy+1)];
						T v9 = Im[grid_2d.sub_2_ind(ix+1, iy+1)];

						T v_s = v1+v2+v3+v4+v+v6+v7+v8+v9;

						T x1 = grid_2d.rx(ix-1);
						T x2 = grid_2d.rx(ix);
						T x3 = grid_2d.rx(ix+1);

						T y1 = grid_2d.ry(iy-1);
						T y2 = grid_2d.ry(iy);
						T y3 = grid_2d.ry(iy+1);

						T x = v1*x1 + v2*x2 + v3*x3 + v4*x1 + v*x2 + v6*x3 + v7*x1 + v8*x2 + v9*x3;
						T y = v1*y1 + v2*y2 + v3*y3 + v4*y1 + v*y2 + v6*y3 + v7*y1 + v8*y2 + v9*y3;
						peak = R_2d<T>(x, y)/v_s;

						return (v1<=v) && (v2<=v) && (v3<=v) && (v4<=v) && (v6<=v) && (v7<=v) && (v8<=v) && (v9<=v);
					};

					auto npeaks_m = static_cast<dt_int32>(::ceil(grid_2d.bs_x*grid_2d.bs_y/(c_pi<T>*sigma*sigma)));

					x.reserve(2*npeaks_m);
					y.reserve(2*npeaks_m);

					// get local peaks
					auto thr_peaks = [&](const iThread_Rect_2d& range, Grid_2d<T>& grid_2d, TVctr_r &image, TVctr_r &x_o, TVctr_r &y_o)
					{
						TVctr_r x;
						x.reserve(npeaks_m);

						TVctr_r y;
						y.reserve(npeaks_m);

						auto ix_0 = 1 + range.ix_0;
						auto ix_e = 1 + range.ix_e;
						auto iy_0 = 1 + range.iy_0;
						auto iy_e = 1 + range.iy_e;

						for(auto ix = ix_0; ix < ix_e; ix++)
						{
							for(auto iy = iy_0; iy < iy_e; iy++)
							{
								R_2d<T> peak;
								if (krn_maximum(ix, iy, grid_2d, image, peak))
								{
									x.push_back(peak.x);
									y.push_back(peak.y);
								}
							}
						}

						stream.stream_mutex.lock();
						x_o.insert(x_o.end(), x.begin(), x.end());
						y_o.insert(y_o.end(), y.begin(), y.end());
						stream.stream_mutex.unlock();
					};

					stream.set_n_stream_act(grid_2d.nx-2);
					stream.set_grid(grid_2d.nx-2, grid_2d.ny-2);
					stream.exec(thr_peaks, grid_2d, image, x, y);

					x.shrink_to_fit();
					y.shrink_to_fit();
				}

				void set_init_bg_values(Stream<edev_cpu>& stream, vector<dt_bool>& mask, TVctr_r &image, 
				TVctr_r &x, TVctr_r &y, T& bg, T& bg_min, T& bg_max)
				{
					T R_max = 2.5*grid_2d.dR_min();
					T peak_min = image[grid_2d.rv_2_ir_bfds(x[0], y[0])];

					auto get_mean_peak = [](const R_2d<T>& p, const T& R_max, const Grid_2d<T>& grid_2d, const TVctr_r &image)->T
					{
						T R2_max = pow(R_max, 2);
						auto range = grid_2d.region_ind(p, R_max);
						T Ixy_mean = 0;
						dt_int32 n_Ixy_mean = 0;
						for(auto ix = range.ix_0; ix < range.ix_e; ix++)
						{
							for(auto iy = range.iy_0; iy < range.iy_e; iy++)
							{
								T R2_d = grid_2d.r2(ix, iy, p.x, p.y);
								if (R2_d < R2_max)
								{
									Ixy_mean += image[grid_2d.sub_2_ind(ix, iy)];
									n_Ixy_mean++;
								}
							}
						}
						Ixy_mean /= n_Ixy_mean;
						return Ixy_mean;
					};

					auto thr_min_peak = [&](const iThread_Rect_2d& range, Grid_2d<T>& grid_2d, vector<dt_bool>& mask, TVctr_r &image, 
					TVctr_r &x, TVctr_r &y, T& peak_min)
					{
						T peak_min_t = 100*image[grid_2d.rv_2_ir_bfds(x[0], y[0])];
						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							R_2d<T> p(x[ipk], y[ipk]);
							if (mask[grid_2d.rv_2_ir_bfds(p.x, p.y)])
							{
								peak_min_t = ::fmin(peak_min_t, get_mean_peak(p, R_max, grid_2d, image));
							}
						}
						stream.stream_mutex.lock();
						peak_min = ::fmin(peak_min, peak_min_t);
						stream.stream_mutex.unlock();
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_min_peak, grid_2d, mask, image, x, y, peak_min);

					// fcn_mean bellow thr
					T I_min = 0;
					T I_mean = 0;
					T I_std = 0;
					T I_low = 0;
					for(auto it=0; it<2; it++)
					{
						I_min = peak_min;
						T I_sum = 0;
						T I_sum_ee = 0;
						T I_sum2 = 0;
						T I_sum2_ee = 0;
						dt_int32 I_c = 0;
						for(auto ix = range_ct.ix_0; ix < range_ct.ix_e; ix++)
						{
							for(auto iy = range_ct.iy_0; iy < range_ct.iy_e; iy++)
							{
								T v = image[grid_2d.sub_2_ind(ix, iy)];
								if ((v < peak_min) && (I_low < v))
								{
									fcn_kh_sum(I_sum, v, I_sum_ee);
									fcn_kh_sum(I_sum2, v*v, I_sum2_ee);
									I_min = ::fmin(I_min, v);
									I_c++;
								}
							}
						}
						I_mean = I_sum/I_c;
						I_std = sqrt(fabs(I_sum2/I_c-I_mean*I_mean));
						peak_min = I_mean;
						I_low = I_mean-3*I_std;
					}
					bg = ::fmax(0, I_mean - I_std);
					// bg_o = (0-image_mean)*sf + ref_mean;
					bg_min = ::fmax(I_mean-4*I_std, I_min-0.5*I_std);
					bg_max = I_mean+0.5*I_std;
				}

				void set_init_A_values(Stream<edev_cpu>& stream, TVctr_r &image, Neigh_2d<T>& neigh, 
				TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &A_min, TVctr_r &A_max, T Bg)
				{
					auto thr_A_values = [Bg](const iThread_Rect_2d& range, Grid_2d<T>& grid_2d, TVctr_r &image, 
					Neigh_2d<T>& neigh, TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &A_min, TVctr_r &A_max)
					{
						auto get_A_values = [](const R_2d<T>& p, const T& R_max, const Grid_2d<T>& grid_2d, 
						const TVctr_r &image, T& A, T& A_min, T& A_max)
						{
							T R2_max = pow(R_max, 2);
							T R2_peak_max = pow(1.75*grid_2d.dR_min(), 2);

							auto range = grid_2d.region_ind(p, R_max);

							T Ixy_min = image[grid_2d.rv_2_ir_bfds(p.x, p.y)];

							T Ixy_peak = 0;
							dt_int32 Ixy_peak_c = 0;

							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T R2_d = grid_2d.r2(ix, iy, p.x, p.y);
									if (R2_d < R2_max)
									{
										T v = image[grid_2d.sub_2_ind(ix, iy)];
										Ixy_min = (v<Ixy_min)?v:Ixy_min;
										if (R2_d < R2_peak_max)
										{
											Ixy_peak += v;
											Ixy_peak_c++;
										}
									}
								}
							}
							Ixy_peak /= Ixy_peak_c;

							A = Ixy_peak;
							A_min = Ixy_peak - 0.75*(Ixy_peak-Ixy_min);
							A_max = Ixy_peak + 0.25*(Ixy_peak-Ixy_min);
						};

						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							R_2d<T> p(x[ipk], y[ipk]);
							T R_max = neigh.d_min(ipk);
							get_A_values(p, R_max, grid_2d, image, A[ipk], A_min[ipk], A_max[ipk]);
							A[ipk] -= Bg;
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_A_values, grid_2d, image, neigh, x, y, A, A_min, A_max);
				}

				vector<dt_bool> get_mask(Stream<edev_cpu>& stream, TVctr_r &x, TVctr_r &y, TVctr_r &S)
				{
					vector<dt_bool> image(grid_2d.size(), false);

					auto set_area = [](const R_2d<T>& p, const T& R_max, const Grid_2d<T>& grid_2d, vector<dt_bool>& image)
					{
						T R2_max = pow(R_max, 2);
						auto range = grid_2d.region_ind(p, R_max);
						for(auto ix = range.ix_0; ix < range.ix_e; ix++)
						{
							for(auto iy = range.iy_0; iy < range.iy_e; iy++)
							{
								T R2_d = grid_2d.r2(ix, iy, p.x, p.y);
								if (R2_d < R2_max)
								{
									image[grid_2d.sub_2_ind(ix, iy)] = true;
								}
							}
						}
					};

					auto thr_area = [&](const iThread_Rect_2d& range, Grid_2d<T>& grid_2d, TVctr_r &x, TVctr_r &y, TVctr_r &S, vector<dt_bool>& image)
					{
						for(auto ipk = range.ind_0; ipk < range.ind_e; ipk++)
						{
							R_2d<T> p(x[ipk], y[ipk]);
							T R_max = 3.5*S[ipk];
							set_area(p, R_max, grid_2d, image);
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_area, grid_2d, x, y, S, image);

					return image;
				}

				vector<dt_bool> get_mask(Stream<edev_cpu>& stream, Neigh_2d<T>& neigh, TVctr_r &x, TVctr_r &y)
				{
					vector<dt_bool> image(grid_2d.size(), false);

					auto set_area = [](const R_2d<T>& p, const T& R_max, const Grid_2d<T>& grid_2d, vector<dt_bool>& image)
					{
						T R2_max = pow(R_max, 2);
						auto range = grid_2d.region_ind(p, R_max);
						for(auto ix = range.ix_0; ix < range.ix_e; ix++)
						{
							for(auto iy = range.iy_0; iy < range.iy_e; iy++)
							{
								T R2_d = grid_2d.r2(ix, iy, p.x, p.y);
								if (R2_d < R2_max)
								{
									image[grid_2d.sub_2_ind(ix, iy)] = true;
								}
							}
						}
					};

					auto thr_area = [&](const iThread_Rect_2d& range, Grid_2d<T>& grid_2d, Neigh_2d<T>& neigh, TVctr_r &x, TVctr_r &y, vector<dt_bool>& image)
					{
						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							R_2d<T> p(x[ipk], y[ipk]);
							T R_max = neigh.d_min(ipk);
							set_area(p, R_max, grid_2d, image);
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_area, grid_2d, neigh, x, y, image);

					// left
					for(auto ix = 0; ix < range_ct.ix_0; ix++)
					{
						for(auto iy = 0; iy < grid_2d.ny; iy++)
						{
							image[grid_2d.sub_2_ind(ix, iy)] = false;
						}
					}

					// right
					for(auto ix = range_ct.ix_e; ix < grid_2d.nx; ix++)
					{
						for(auto iy = 0; iy < grid_2d.ny; iy++)
						{
							image[grid_2d.sub_2_ind(ix, iy)] = false;
						}
					}

					// top
					for(auto ix = 0; ix < grid_2d.nx; ix++)
					{
						for(auto iy = 0; iy < range_ct.iy_0; iy++)
						{
							image[grid_2d.sub_2_ind(ix, iy)] = false;
						}
					}

					// bottom
					for(auto ix = 0; ix < grid_2d.nx; ix++)
					{
						for(auto iy = range_ct.iy_e; iy < grid_2d.ny; iy++)
						{
							image[grid_2d.sub_2_ind(ix, iy)] = false;
						}
					}

					return image;
				}

				void set_mask_and_scale(Stream<edev_cpu>& stream, vector<dt_bool>& mask, 
				TVctr_r &image, TVctr_r &image_m, T& image_m_sc)
				{
					image_m.resize(image.size());
					T image_sum = 0;
					T image2_sum = 0;
					dt_int32 c_image_sum = 0;

					// mask image
					auto thr_mask = [&](const iThread_Rect_2d& range, TVctr_r &image, 
					TVctr_r &image_m, dt_int32& c_image_sum, T& image_sum, T& image2_sum)
					{
						// set mask and get fcn_mean and std
						T im_sum = 0;
						T im_sum_ee = 0;
						T im2_sum = 0;
						T im2_sum_ee = 0;
						dt_int32 c_sum = 0;
						for(auto ixy = range.ind_0; ixy<range.ind_e; ixy++)
						{
							if (mask[ixy])
							{
								T v = image[ixy];
								image_m[ixy] = v;
								fcn_kh_sum(im_sum, v, im_sum_ee);
								fcn_kh_sum(im2_sum, v*v, im2_sum_ee);
								c_sum++;
							}
							else
							{
								image_m[ixy] = 0;
							}
						}
						stream.stream_mutex.lock();
						image_sum += im_sum;
						image2_sum += im2_sum;
						c_image_sum += c_sum;
						stream.stream_mutex.unlock();
					};

					stream.set_n_stream_act(image.size());
					stream.set_grid(image.size(), 1);
					stream.exec(thr_mask, image, image_m, c_image_sum, image_sum, image2_sum);

					T image_m_mean = image_sum/c_image_sum;
					T image_m_std = sqrt(fabs(image2_sum/c_image_sum-image_m_mean*image_m_mean));

					// shift fcn_scale image
					image_m_sc = image_m_std;

					auto thr_sft_scale = [image_m_sc](const iThread_Rect_2d& range, vector<dt_bool>& mask, TVctr_r &image)
					{
						for(auto ixy = range.ind_0; ixy<range.ind_e; ixy++)
						{
							if (mask[ixy])
							{
								image[ixy] /= image_m_sc;
							}
						}
					};

					stream.set_n_stream_act(image.size());
					stream.set_grid(image.size(), 1);
					stream.exec(thr_sft_scale, mask, image_m);
				}

				void fcn_scale(Stream<edev_cpu>& stream, TVctr_r &image, TVctr_r &image_m, T& image_m_sc)
				{
					image_m.resize(image.size());
				
					image_m_sc = sqrt(fcn_variance(stream, image));

					auto thr_sft_scale = [image_m_sc](const iThread_Rect_2d& range, TVctr_r &image, TVctr_r &image_m)
					{
						for(auto ixy = range.ind_0; ixy<range.ind_e; ixy++)
						{
							image_m[ixy] = image[ixy]/image_m_sc;
						}
					};

					stream.set_n_stream_act(image.size());
					stream.set_grid(image.size(), 1);
					stream.exec(thr_sft_scale, image, image_m);
				}

				void set_init_sigma_values(T sigma, TVctr_r &S, TVctr_r &S_min, TVctr_r &S_max)
				{
					dt_int32 m = x.size();
					A.resize(m);

					S.resize(m);
					S_min.resize(m);
					S_max.resize(m);

					for(auto ipk=0; ipk<m; ipk++)
					{
						// T radius = neigh.radius_min(ipk, grid_2d, Im);
						S[ipk] = sigma;
						S_min[ipk] = ::fmax(0.5*grid_2d.dR_min(), 0.02*sigma);
						S_max[ipk] = ::fmin(0.9*neigh.d_min(ipk), 3*sigma);
					}
				}

				void set_sigma_limits(T f_sigma, TVctr_r &S, TVctr_r &S_min, TVctr_r &S_max)
				{
					T s_mean = fcn_mean(S);
					T s_std = f_sigma*sqrt(fcn_variance(S));
					for(auto ipk=0; ipk<S.size(); ipk++)
					{
						T sig_l = ::fmax(S_min[ipk], s_mean - s_std);
						T sig_u = ::fmin(S_max[ipk], s_mean + s_std);

						S[ipk] = ((S[ipk]<sig_l)||(S[ipk]>sig_u))?s_mean:S[ipk];
					}
				}

				TVctr_r gauss_sp(Stream<edev_cpu>& stream, TVctr_r &x, TVctr_r &y, 
				TVctr_r &A, TVctr_r &S, T Bg)
				{
					auto thr_gauss_sp = [&stream](const Grid_2d<T>& grid_2d, TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &S, T Bg)
					{
						dt_int32 m = grid_2d.size();
						TVctr_r Ixy(m, Bg);

						auto thr_gauss_sup = [&](const iThread_Rect_2d& range, TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &S)
						{
							auto gaussian = [](R_2d<T> p, T a, T b, T R_max, const Grid_2d<T>& grid_2d, TVctr_r &Ixy)
							{
								T R2_max = pow(R_max, 2);
								T c_alpha = 0.5/pow(b, 2);
								auto range = grid_2d.region_ind(p, R_max);
								for(auto ix = range.ix_0; ix < range.ix_e; ix++)
								{
									for(auto iy = range.iy_0; iy < range.iy_e; iy++)
									{
										T r2 = grid_2d.r2(ix, iy, p.x, p.y);
										if (r2 < R2_max)
										{
											Ixy[grid_2d.sub_2_ind(ix, iy)] += a*exp(-c_alpha*r2);
										}
									}
								}
							};

							TVctr_r Ixy_p(m, T(0));
							for(auto ip = range.ind_0; ip < range.ind_e; ip++)
							{
								R_2d<T> p(x[ip], y[ip]);
								T a = A[ip];
								T b = S[ip];
								T R_max = 4.0*b;
								gaussian(p, a, b, R_max, grid_2d, Ixy_p);
							}

							stream.stream_mutex.lock();
							for(auto ixy = 0; ixy < m; ixy++)
							{
								Ixy[ixy] += Ixy_p[ixy];
							}
							stream.stream_mutex.unlock();
						};

						stream.set_n_stream_act(x.size());
						stream.set_grid(x.size(), 1);
						stream.exec(thr_gauss_sup, x, y, A, S);
						return Ixy;
					};

					return thr_gauss_sp(grid_2d, x, y, A, S, Bg);
				};

				TVctr_r fit_b_0(Stream<edev_cpu>& stream, vector<dt_bool>& mask, TVctr_r &image, 
				TVctr_r x, TVctr_r y, TVctr_r A, T S, T Bg, dt_int32 niter, T d_error)
				{
					// mask and scaling
					TVctr_r Ixy;
					T Ixy_sc;
					set_mask_and_scale(stream, mask, image, Ixy, Ixy_sc);

					const T factor = 2;
					T Rxy_sc = grid_2d.bs_max()/factor;

					Grid_2d<T> grid_s(grid_2d.nx, grid_2d.ny, grid_2d.bs_x/Rxy_sc, grid_2d.bs_y/Rxy_sc);

					// scaling
					for(auto ip = 0; ip<x.size(); ip++)
					{
						A[ip] /= Ixy_sc;
						x[ip] /= Rxy_sc;
						y[ip] /= Rxy_sc;
					}
					T b = S/Rxy_sc;
					T b_min = 0.01*b;
					T b_max = 4*b;

					T c = Bg/Ixy_sc;

					for(auto ixy = 0; ixy < Ixy.size(); ixy++)
					{
						Ixy[ixy] -= c;
					}

					auto get_b = [&stream](const Grid_2d<T>& grid_2d, vector<dt_bool>& mask, TVctr_r &Ixy, 
					TVctr_r &x, TVctr_r &y, TVctr_r &A, T& b, T b_min, T b_max, dt_int32 niter, T d_error)
					{
						auto dgaussian = [&mask](R_2d<T> p, T a, T b, T R_max, const Grid_2d<T>& grid_2d, TVctr_r &J, TVctr_r &d_Ixy)
						{
							T R2_max = pow(R_max, 2);
							T c_alpha = 0.5/pow(b, 2);
							T c_a = a/pow(b, 3);
							auto range = grid_2d.region_ind(p, R_max);

							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T r2 = grid_2d.r2(ix, iy, p.x, p.y);
									dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
									if ((r2 < R2_max) && mask[ixy])
									{
										T f = exp(-c_alpha*r2);
										J[ixy] += c_a*r2*f;
										d_Ixy[ixy] += a*f;
									}
								}
							}
						};

						auto solve_d_b = [&mask](TVctr_r &J, TVctr_r &d_Ixy)->T
						{
							T sxx = 0;
							T sxx_ee = 0;
							T sxy = 0;
							T sxy_ee = 0;
							for(auto ixy = 0; ixy < d_Ixy.size(); ixy++)
							{
								if (mask[ixy])
								{
									T x = J[ixy];
									T y = d_Ixy[ixy];

									fcn_kh_sum(sxx, x*x, sxx_ee);
									fcn_kh_sum(sxy, x*y, sxy_ee);
								}
							}

							return sxy/sxx;
						};

						dt_int32 m = Ixy.size();
						TVctr_r d_Ixy(m);
						TVctr_r J(m);
						for(auto iter = 0; iter<niter; iter++)
						{
							std::fill(J.begin(), J.end(), 0);
							d_Ixy = Ixy;

							auto thr_gauss_sup = [&](const iThread_Rect_2d& range)
							{
								TVctr_r d_Ixy_p(m, T(0));
								TVctr_r J_p(m, T(0));
								for(auto ip = range.ind_0; ip < range.ind_e; ip++)
								{
									R_2d<T> p(x[ip], y[ip]);
									T a = A[ip];
									T R_max = 4.0*b;
									dgaussian(p, a, b, R_max, grid_2d, J_p, d_Ixy_p);
								}
								stream.stream_mutex.lock();
								for(auto ixy = 0; ixy < m; ixy++)
								{
									J[ixy] += J_p[ixy];
									d_Ixy[ixy] -= d_Ixy_p[ixy];
								}
								stream.stream_mutex.unlock();
							};

							stream.set_n_stream_act(x.size());
							stream.set_grid(x.size(), 1);
							stream.exec(thr_gauss_sup);

							T d_b = solve_d_b(J, d_Ixy);

							b += d_b;

							b = min(max(b_min, b), b_max);

							if (fabs(d_b/b)<d_error)
							{
								break;
							}
						}
					};

					get_b(grid_s, mask, Ixy, x, y, A, b, b_min, b_max, niter, d_error);

					TVctr_r S_o(A.size(), b*Rxy_sc);
					return S_o;
				}

				TVctr_r fit_a_0(Stream<edev_cpu>& stream, vector<dt_bool>& mask, TVctr_r &image, 
				TVctr_r x, TVctr_r y, TVctr_r A, TVctr_r &A_min, TVctr_r &A_max, TVctr_r S, T Bg)
				{
					// mask and scaling
					TVctr_r Ixy;
					T Ixy_sc;
					set_mask_and_scale(stream, mask, image, Ixy, Ixy_sc);

					const T factor = 2;
					T Rxy_sc = grid_2d.bs_max()/factor;

					Grid_2d<T> grid_s(grid_2d.nx, grid_2d.ny, grid_2d.bs_x/Rxy_sc, grid_2d.bs_y/Rxy_sc);

					// scaling
					for(auto ip = 0; ip<x.size(); ip++)
					{
						A[ip] /= Ixy_sc;
						S[ip] /= Rxy_sc;
						x[ip] /= Rxy_sc;
						y[ip] /= Rxy_sc;
					}

					T c = Bg/Ixy_sc;

					for(auto ixy = 0; ixy < Ixy.size(); ixy++)
					{
						Ixy[ixy] -= c;
					}

					auto get_a = [&stream](const Grid_2d<T>& grid_2d, vector<dt_bool>& mask, TVctr_r &Ixy, 
						TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &S)
					{
						auto dgaussian = [&mask](R_2d<T> p, T a, T b, T R_max, const Grid_2d<T>& grid_2d, TVctr_r &J)
						{
							T R2_max = pow(R_max, 2);
							T c_alpha = 0.5/pow(b, 2);
							auto range = grid_2d.region_ind(p, R_max);
							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T r2 = grid_2d.r2(ix, iy, p.x, p.y);
									dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
									if ((r2 < R2_max) && mask[ixy])
									{
										J[ixy] += a*exp(-c_alpha*r2);
									}
								}
							}
						};

						dt_int32 m = Ixy.size();
						TVctr_r J(m, T(0));
						auto thr_gauss_sup = [&](const iThread_Rect_2d& range)
						{
							TVctr_r J_p(m, T(0));
							for(auto ip = range.ind_0; ip < range.ind_e; ip++)
							{
								R_2d<T> p(x[ip], y[ip]);
								T a = A[ip];
								T b = S[ip];
								T R_max = 4.0*b;
								dgaussian(p, a, b, R_max, grid_2d, J_p);
							}

							stream.stream_mutex.lock();
							for(auto ixy = 0; ixy < m; ixy++)
							{
								J[ixy] += J_p[ixy];
							}
							stream.stream_mutex.unlock();
						};

						stream.set_n_stream_act(x.size());
						stream.set_grid(x.size(), 1);
						stream.exec(thr_gauss_sup);

						T sxy = 0;
						T sxy_ee = 0;
						T sxx = 0;
						T sxx_ee = 0;
						for(auto ixy = 0; ixy < Ixy.size(); ixy++)
						{
							if (mask[ixy])
							{
								fcn_kh_sum(sxy, Ixy[ixy]*J[ixy], sxy_ee);
								fcn_kh_sum(sxx, J[ixy]*J[ixy], sxx_ee);
							}
						}
						return sxy/sxx;
					};

					T alpha = get_a(grid_s, mask, Ixy, x, y, A, S);
					for(auto ip = 0; ip < x.size(); ip++)
					{
						T a = A[ip]*alpha*Ixy_sc + Bg;
						if ((a<A_min[ip])||(a>A_max[ip]))
						{
							a = A[ip]*Ixy_sc + Bg;
						}
						A[ip] = a-Bg;
					}
					return A;
				};

				T fit_c_0(Stream<edev_cpu>& stream, TVctr_r &image, TVctr_r x, TVctr_r y, 
				TVctr_r A, TVctr_r S, T Bg, T Bg_min, T Bg_max)
				{
					// mask and scaling
					TVctr_r Ixy;
					T Ixy_sc;
					fcn_scale(stream, image, Ixy, Ixy_sc);

					const T factor = 2;
					T Rxy_sc = grid_2d.bs_max()/factor;

					Grid_2d<T> grid_s(grid_2d.nx, grid_2d.ny, grid_2d.bs_x/Rxy_sc, grid_2d.bs_y/Rxy_sc);

					// scaling
					for(auto ip = 0; ip<x.size(); ip++)
					{
						A[ip] /= Ixy_sc;
						S[ip] /= Rxy_sc;
						x[ip] /= Rxy_sc;
						y[ip] /= Rxy_sc;
					}

					T c_0 = Bg/Ixy_sc;
					T c_min = Bg_min/Ixy_sc;
					T c_max = Bg_max/Ixy_sc;

					auto get_c = [&](const Grid_2d<T>& grid_2d, vector<dt_bool>& mask, TVctr_r &Ixy, 
						TVctr_r &x, TVctr_r &y, TVctr_r &A, TVctr_r &S)
					{
						auto dgaussian = [](const R_2d<T>& p, const T& a, const T& b, const T& R_max, const Grid_2d<T>& grid_2d, TVctr_r &J)
						{
							T R2_max = pow(R_max, 2);
							T c_alpha = 0.5/pow(b, 2);
							auto range = grid_2d.region_ind(p, R_max);
							for(auto ix = range.ix_0; ix < range.ix_e; ix++)
							{
								for(auto iy = range.iy_0; iy < range.iy_e; iy++)
								{
									T r2 = grid_2d.r2(ix, iy, p.x, p.y);
									if (r2 < R2_max)
									{
										J[grid_2d.sub_2_ind(ix, iy)] += a*exp(-c_alpha*r2);
									}
								}
							}
						};

						dt_int32 m = Ixy.size();
						TVctr_r d_Ixy = Ixy;
						auto thr_gauss_sup = [&](const iThread_Rect_2d& range)
						{
							TVctr_r d_Ixy_p(m, T(0));
							for(auto ip = range.ind_0; ip < range.ind_e; ip++)
							{
								R_2d<T> p(x[ip], y[ip]);
								T a = A[ip];
								T b = S[ip];
								T R_max = 4.0*b;
								dgaussian(p, a, b, R_max, grid_2d, d_Ixy_p);
							}

							stream.stream_mutex.lock();
							for(auto ixy = 0; ixy < m; ixy++)
							{
								d_Ixy[ixy] -= d_Ixy_p[ixy];
							}
							stream.stream_mutex.unlock();
						};

						stream.set_n_stream_act(x.size());
						stream.set_grid(x.size(), 1);
						stream.exec(thr_gauss_sup);

						T sxy = 0;
						T sxy_ee = 0;
						dt_int32 sxx = 0;
						for(auto ix = range_ct.ix_0; ix < range_ct.ix_e; ix++)
						{
							for(auto iy = range_ct.iy_0; iy < range_ct.iy_e; iy++)
							{
								dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
								fcn_kh_sum(sxy, d_Ixy[ixy], sxy_ee);
								sxx++;
							}
						}
						return sxy/sxx;
					};

					T c = get_c(grid_s, mask, Ixy, x, y, A, S);
					c = ((c<c_min)||(c>c_max))?c_0:c_0+(1.0/2.0)*(c-c_0);
					return c*Ixy_sc;
				};

				void fit_a_b_0(Stream<edev_cpu>& stream, Regions &region, TVctr_r &x, TVctr_r &y, 
				TVctr_r &A, TVctr_r &A_min, TVctr_r &A_max, TVctr_r &S, TVctr_r &S_min, 
				TVctr_r &S_max, T Bg, dt_int32 niter, T d_error)
				{
					auto image_s = gauss_sp(stream, x, y, A, S, Bg);

					auto thr_fit_a_b_0 = [&](const iThread_Rect_2d& range)
					{
						dt_int32 m_max = region.m_max;
						dt_int32 n_max = 2;

						TVctr_r Ja(m_max);
						TVctr_r Jb(m_max);
						TVctr_r d_Ixy(m_max);

						auto solve_a_b = [](dt_int32 m, TVctr_r &Ja, TVctr_r &Jb, 
						TVctr_r &d_Ixy, T lambda, T& d_a, T& d_b, T& g_m, T& rho_f)
						{
							T sx1x1 = 0;
							T sx2x2 = 0;
							T sx1x2 = 0;
							T sx1y = 0;
							T sx2y = 0;
							for(auto ixy = 0; ixy < m; ixy++)
							{
								T x1 = Ja[ixy];
								T x2 = Jb[ixy];
								T y = d_Ixy[ixy];

								sx1x1 += x1*x1;
								sx2x2 += x2*x2;
								sx1x2 += x1*x2;
								sx1y += x1*y;
								sx2y += x2*y;
							}
							T D_a = lambda*(sx1x1 > 1e-10)?sx1x1:1;
							T D_b = lambda*(sx2x2 > 1e-10)?sx2x2:1;

							sx1x1 += D_a;
							sx2x2 += D_b;

							T det = sx1x1*sx2x2 - sx1x2*sx1x2;
							d_a = (sx2x2*sx1y - sx1x2*sx2y)/det;
							d_b = (sx1x1*sx2y - sx1x2*sx1y)/det;

							g_m = ::fmax(fabs(sx1y), fabs(sx2y));
							rho_f = d_a*(D_a*d_a + sx1y) + d_b*(D_b*d_b + sx2y);
						};

						auto get_chi2 = [](TVctr_r &R2, TVctr_r &Ixy, T a, T b)->T
						{
							T c_alpha = 0.5/pow(b, 2);
							T chi2 = 0;
							T chi2_ee = 0;
							for(auto im = 0; im < R2.size(); im++)
							{
								T v = Ixy[im]-a*exp(-c_alpha*R2[im]);
								fcn_kh_sum(chi2, v*v, chi2_ee);
							}
							return chi2;
						};

						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							TVctr_r &R2 = region[ipk].R2;

							T x_0 = region[ipk].sft_Rx(x[ipk]);
							T y_0 = region[ipk].sft_Ry(y[ipk]);
							T a = A[ipk]/region[ipk].Ixy_sc;
							T b = S[ipk]/region[ipk].Rxy_sc;

							TVctr_r Ixy = region[ipk].sft_Ixy(grid_2d, image_s, x_0, y_0, a, b);

							T a_0 = a;
							T b_0 = b;

							T a_min = (A_min[ipk]-Bg)/region[ipk].Ixy_sc;
							T a_max = (A_max[ipk]-Bg)/region[ipk].Ixy_sc;

							T b_min = S_min[ipk]/region[ipk].Rxy_sc;
							T b_max = S_max[ipk]/region[ipk].Rxy_sc;

							dt_int32 m = R2.size();
							dt_int32 n = n_max;

							T chi2 = get_chi2(R2, Ixy, a, b);
							T chi2_f = m-n+1;

							T lambda = 2;
							T lambda_f = 2;

							for(auto iter = 0; iter<niter; iter++)
							{
								T c_alpha = 0.5/pow(b, 2);
								T c_b = a/pow(b, 3);
								for(auto im = 0; im < m; im++)
								{
									T r2 = R2[im];
									T f = exp(-c_alpha*r2);
									Ja[im] = f;
									Jb[im] = c_b*r2*f;
									d_Ixy[im] = Ixy[im]-a*f;
								}

								T d_a, d_b, g_m, f_rho;
								solve_a_b(m, Ja, Jb, d_Ixy, lambda, d_a, d_b, g_m, f_rho);

								T a_t = a + d_a;
								T b_t = b + d_b;
								T chi2_t = get_chi2(R2, Ixy, a_t, b_t);
								T rho = (chi2-chi2_t)/f_rho;

								if ((g_m<1e-6)||(chi2/chi2_f<1e-6)||fabs(chi2-chi2_t)<1e-6)
								{
									break;
								}

								if (rho>0)
								{
									a = a_t;
									b = b_t;
								
									// a = ((a < a_min)||(a>a_max))?a_0:a;
									// b = ((b<b_min)||(b>b_max))?b_0:b;

									a = min(max(a, a_min), a_max);
									b = min(max(b, b_min), b_max);

									chi2 = get_chi2(R2, Ixy, a, b);
								}

								lambda = (rho>1e-3)?(::fmax(lambda/lambda_f, 1e-7)):(::fmin(lambda*lambda_f, 1e+7));

								// a = min(max(a, a_min), a_max);
								// b = min(max(b, b_min), b_max);

								// if (fabs(d_b/b)<d_error)
								// {
								// 	break;
								// }
							}

							A[ipk] = a*region[ipk].Ixy_sc;
							S[ipk] = b*region[ipk].Rxy_sc;
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_fit_a_b_0);

					T f_sigma = 3.0;
					// set_sigma_limits(f_sigma, S, S_min, S_max);
				}

				void fit_x_y_0(Stream<edev_cpu>& stream, Regions &region, TVctr_r &x, TVctr_r &y, 
				TVctr_r &A, TVctr_r &S, T Bg, dt_int32 niter, T d_error)
				{
					auto image_s = gauss_sp(stream, x, y, A, S, Bg);

					auto thr_fit_x_y_0 = [&](const iThread_Rect_2d& range)
					{
						dt_int32 m_max = region.m_max;
						dt_int32 n_max = 2;

						TVctr_r Jx(m_max);
						TVctr_r Jy(m_max);
						TVctr_r d_Ixy(m_max);

						auto solve_x_y = [](dt_int32 m, TVctr_r &Jx, TVctr_r &Jy, TVctr_r &d_Ixy, T& d_x, T& d_y)
						{
							T sx1x1 = 0;
							T sx2x2 = 0;
							T sx1x2 = 0;
							T sx1y = 0;
							T sx2y = 0;
							for(auto ixy = 0; ixy < m; ixy++)
							{
								T x1 = Jx[ixy];
								T x2 = Jy[ixy];
								T y = d_Ixy[ixy];

								sx1x1 += x1*x1;
								sx2x2 += x2*x2;
								sx1x2 += x1*x2;
								sx1y += x1*y;
								sx2y += x2*y;
							}
							sx1x1 *= 1.001;
							sx2x2 *= 1.001;

							T det = sx1x1*sx2x2 - sx1x2*sx1x2;
							d_x = (sx2x2*sx1y - sx1x2*sx2y)/det;
							d_y = (sx1x1*sx2y - sx1x2*sx1y)/det;
						};

						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							TVctr_r &Rx = region[ipk].Rx;
							TVctr_r &Ry = region[ipk].Ry;

							T x_0 = region[ipk].sft_Rx(x[ipk]);
							T y_0 = region[ipk].sft_Ry(y[ipk]);
							T a = A[ipk]/region[ipk].Ixy_sc;
							T b = S[ipk]/region[ipk].Rxy_sc;

							TVctr_r Ixy = region[ipk].sft_Ixy(grid_2d, image_s, x_0, y_0, a, b);

							T a_0 = a;
							T b_0 = b;

							T dxy = ::fmax(1.1*grid_2d.dR_min()/region[ipk].Rxy_sc, 0.25*b);
							R_2d<T> p(x_0, y_0);
							R_2d<T> p_0 = p;
							R_2d<T> p_min(p.x-dxy, p.y-dxy);
							R_2d<T> p_max(p.x+dxy, p.y+dxy);

							dt_int32 m = Ixy.size();
							dt_int32 n = 2;

							TVctr_r d_abxy(n);
							for(auto iter = 0; iter<niter; iter++)
							{
								T c_alpha = 0.5/pow(b, 2);
								T c_xy = a/pow(b, 2);
								for(auto im = 0; im < m; im++)
								{
									T rx = Rx[im]-p.x;
									T ry = Ry[im]-p.y;
									T r2 = rx*rx+ry*ry;
									T f = exp(-c_alpha*r2);

									Jx[im] = c_xy*rx*f;
									Jy[im] = c_xy*ry*f;
									d_Ixy[im] = Ixy[im]-a*f;
								}

								T d_x, d_y;
								solve_x_y(m, Jx, Jy, d_Ixy, d_x, d_y);

								p += R_2d<T>(d_x, d_y);

								T ff = 0.5;
								p.x = ((p.x<p_min.x)||(p.x>p_max.x))?p_0.x:p_0.x+ff*(p.x-p_0.x);
								p.y = ((p.y<p_min.y)||(p.y>p_max.y))?p_0.y:p_0.y+ff*(p.y-p_0.y);

								T d_xy = sqrt(d_x*d_x+d_y*d_y);
								if (d_xy<d_error)
								{
									break;
								}
							}

							x[ipk] = p.x*region[ipk].Rxy_sc+region[ipk].Rx_sf;
							y[ipk] = p.y*region[ipk].Rxy_sc+region[ipk].Ry_sf;
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_fit_x_y_0);

					T f_sigma = 3.0;
					set_sigma_limits(f_sigma, S, S_min, S_max);
				}

				void fit_a_b_x_y_0(Stream<edev_cpu>& stream, Regions &region, TVctr_r &x, TVctr_r &y, 
				TVctr_r &A, TVctr_r &A_min, TVctr_r &A_max, TVctr_r &S, TVctr_r &S_min, 
				TVctr_r &S_max, T Bg, dt_int32 niter, T d_error)
				{
					auto image_s = gauss_sp(stream, x, y, A, S, Bg);

					auto thr_fit_a_b_x_y_0 = [&](const iThread_Rect_2d& range)
					{
						lapack::LSF_1<T> fls_1;
				
						dt_int32 m_max = region.m_max;
						dt_int32 n_max = 4;

						TVctr_r M(m_max*n_max);
						TVctr_r d_Ixy(m_max);

						for(auto ipk=range.ind_0; ipk<range.ind_e; ipk++)
						{
							TVctr_r &Rx = region[ipk].Rx;
							TVctr_r &Ry = region[ipk].Ry;

							T x_0 = region[ipk].sft_Rx(x[ipk]);
							T y_0 = region[ipk].sft_Ry(y[ipk]);
							T a = A[ipk]/region[ipk].Ixy_sc;
							T b = S[ipk]/region[ipk].Rxy_sc;

							TVctr_r Ixy = region[ipk].sft_Ixy(grid_2d, image_s, x_0, y_0, a, b);

							T a_0 = a;
							T b_0 = b;

							T a_min = (A_min[ipk]-Bg)/region[ipk].Ixy_sc;
							T a_max = (A_max[ipk]-Bg)/region[ipk].Ixy_sc;

							T b_min = S_min[ipk]/region[ipk].Rxy_sc;
							T b_max = S_max[ipk]/region[ipk].Rxy_sc;

							T dxy = ::fmax(1.1*grid_2d.dR_min()/region[ipk].Rxy_sc, 0.25*b);
							R_2d<T> p(x_0, y_0);
							R_2d<T> p_0 = p;
							R_2d<T> p_min(p.x-dxy, p.y-dxy);
							R_2d<T> p_max(p.x+dxy, p.y+dxy);

							dt_int32 m = Ixy.size();
							dt_int32 n = 4;

							TVctr_r d_abxy(n);
							for(auto iter = 0; iter<niter; iter++)
							{
								T c_alpha = 0.5/pow(b, 2);
								T c_b = a/pow(b, 3);
								T c_xy = a/pow(b, 2);
								for(auto im = 0; im < m; im++)
								{
									T rx = Rx[im]-p.x;
									T ry = Ry[im]-p.y;
									T r2 = rx*rx+ry*ry;
									T f = exp(-c_alpha*r2);

									M[0*m+im] = f;
									M[1*m+im] = c_b*r2*f;
									M[2*m+im] = c_xy*rx*f;
									M[3*m+im] = c_xy*ry*f;
									d_Ixy[im] = Ixy[im]-a*f;
								}

								fls_1(m, n, M.data(), d_Ixy.data(), d_abxy.data());
								a += d_abxy[0];
								b += d_abxy[1];
								p += R_2d<T>(d_abxy[2], d_abxy[3]);

								T ff = 0.75;
								a = ((a < a_min)||(a>a_max))?a_0:a_0+ff*(a-a_0);
								b = ((b<b_min)||(b>b_max))?b_0:b_0+ff*(b-b_0);
								p.x = ((p.x<p_min.x)||(p.x>p_max.x))?p_0.x:p_0.x+ff*(p.x-p_0.x);
								p.y = ((p.y<p_min.y)||(p.y>p_max.y))?p_0.y:p_0.y+ff*(p.y-p_0.y);

								if (fabs(d_abxy[1]/b)<d_error)
								{
									break;
								}
							}

							A[ipk] = a*region[ipk].Ixy_sc;
							S[ipk] = b*region[ipk].Rxy_sc;
							x[ipk] = p.x*region[ipk].Rxy_sc+region[ipk].Rx_sf;
							y[ipk] = p.y*region[ipk].Rxy_sc+region[ipk].Ry_sf;
						}
					};

					stream.set_n_stream_act(x.size());
					stream.set_grid(x.size(), 1);
					stream.exec(thr_fit_a_b_x_y_0);

					T f_sigma = 3.0;
					set_sigma_limits(f_sigma, S, S_min, S_max);
				}

				void fit_a_b_x_y(Stream<edev_cpu>& stream, Regions &region, TVctr_r &x_i, TVctr_r &y_i, 
				TVctr_r &A_i, TVctr_r &A_min, TVctr_r &A_max, TVctr_r &S_i, TVctr_r &S_min, 
				TVctr_r &S_max, T Bg, dt_int32 niter, T d_error)
				{
					auto image_s = gauss_sp(stream, x, y, A, S, Bg);

					dt_int32 m_max = region.m_max;
					dt_int32 n_max = 4*region.n_max;

					auto thr_fit_a_b_x_y = [&](const iThread_Rect_2d& range)
					{
						TVctr_r d_beta(n_max);
						TVctr_r beta(n_max);
						TVctr_r beta_0(n_max);
						TVctr_r beta_min(n_max);
						TVctr_r beta_max(n_max);

						lapack::LSF_1<T> fls_1;

						auto dgaussian = [](dt_int32 n, TVctr_r &beta, TVctr_r &Rx, TVctr_r &Ry, TVctr_r &J, TVctr_r &d_Ixy)
						{
							dt_int32 m = Rx.size();
							dt_int32 npn = n/4;
							for(auto ip = 0; ip < npn; ip++)
							{
								dt_int32 idx_p = 4*ip;

								R_2d<T> p(beta[idx_p+0], beta[idx_p+1]);
								T a = beta[idx_p+2];
								T b = beta[idx_p+3];
								T r2_lim = pow(4*b, 2);
								T c_alpha = 0.5/pow(b, 2);
								T c_b = a/pow(b, 3);
								T c_xy = a/pow(b, 2);
								for(auto im = 0; im < m; im++)
								{
									T rx = Rx[im]-p.x;
									T ry = Ry[im]-p.y;
									T r2 = rx*rx + ry*ry;
									if (r2 < r2_lim)
									{
										T f = exp(-c_alpha*r2);
										J[(idx_p + 0)*m + im] = c_xy*rx*f;
										J[(idx_p + 1)*m + im] = c_xy*ry*f;
										J[(idx_p + 2)*m + im] = f;
										J[(idx_p + 3)*m + im] = c_b*r2*f;
										d_Ixy[im] -= a*f;
									}
								}
							}
						};

						auto x = x_i;
						auto y = y_i;
						auto A = A_i;
						auto S = S_i;
						for(auto ipk = range.ind_0; ipk < range.ind_e; ipk++)
						{
							TVctr_r &Rx = region[ipk].Rx;
							TVctr_r &Ry = region[ipk].Ry;

							TVctr_r xn = neigh.select(ipk, x, region[ipk].Rx_sf, region[ipk].Rxy_sc);
							TVctr_r yn = neigh.select(ipk, y, region[ipk].Ry_sf, region[ipk].Rxy_sc);
							TVctr_r An = neigh.select(ipk, A, 0, region[ipk].Ixy_sc);
							TVctr_r Sn = neigh.select(ipk, S, 0, region[ipk].Rxy_sc);

							TVctr_r Ixy = region[ipk].sft_Ixy(grid_2d, image_s, xn, yn, An, Sn);

							for(auto ipn = 0; ipn < xn.size(); ipn++)
							{
								dt_int32 idx_p = 4*ipn;

								beta_0[idx_p + 0] = xn[ipn];
								beta_0[idx_p + 1] = yn[ipn];
								beta_0[idx_p + 2] = An[ipn];
								beta_0[idx_p + 3] = Sn[ipn];

								T dxy = ::fmax(1.1*grid_2d.dR_min()/region[ipk].Rxy_sc, 0.25*Sn[ipn]);
								R_2d<T> p(xn[ipn], yn[ipn]);

								beta_min[idx_p + 0] = p.x-dxy;
								beta_min[idx_p + 1] = p.y-dxy;
								beta_min[idx_p + 2] = (A_min[ipn]-Bg)/region[ipk].Ixy_sc;
								beta_min[idx_p + 3] = S_min[ipn]/region[ipk].Rxy_sc;

								beta_max[idx_p + 0] = p.x+dxy;
								beta_max[idx_p + 1] = p.y+dxy;
								beta_max[idx_p + 2] = (A_max[ipn]-Bg)/region[ipk].Ixy_sc;
								beta_max[idx_p + 3] = S_max[ipn]/region[ipk].Rxy_sc;
							}

							beta = beta_0;

							dt_int32 m = Rx.size();
							dt_int32 n = 4*xn.size();

							for(auto iter = 0; iter < niter; iter++)
							{
								TVctr_r J(m*n, 0);
								TVctr_r d_Ixy = Ixy;

								dgaussian(n, beta, Rx, Ry, J, d_Ixy);

								fls_1(m, n, J.data(), d_Ixy.data(), d_beta.data());

								thrust::transform(beta.begin(), beta.end(), d_beta.begin(), beta.begin(), thrust::plus<T>());

								for(auto ip = 0; ip < n; ip++)
								{
									beta[ip] = ((beta[ip]<beta_min[ip])||(beta[ip]>beta_max[ip]))?beta_0[ip]:beta[ip];
								}

								if (fabs(d_beta[3]/beta[3])<d_error)
								{
									break;
								}
							}

							x_i[ipk] = beta[0]*region[ipk].Rxy_sc + region[ipk].Rx_sf;
							y_i[ipk] = beta[1]*region[ipk].Rxy_sc + region[ipk].Ry_sf;
							A_i[ipk] = beta[2]*region[ipk].Ixy_sc;
							S_i[ipk] = beta[3]*region[ipk].Rxy_sc;
						}
					};

					stream.set_n_stream_act(x_i.size());
					stream.set_grid(x_i.size(), 1);
					stream.exec(thr_fit_a_b_x_y);
				}

				Stream<edev_cpu> stream;
				FFT<T, edev_cpu> fft_2d;

				Grid_2d<T> grid_2d;
				TVctr_r image;
				TVctr_r image_den;
				vector<dt_bool> mask;

				T sigma;
				T thres;
				T d_error;
				dt_int32 niter;

				T ref_mean;
				T ref_std;
				T image_mean;
				T image_std;

				T radius_n;
				Neigh_2d<T> neigh;
				Regions region;
				iThread_Rect_2d range_ct;

				TVctr_r x;
				TVctr_r y;
				TVctr_r A;
				TVctr_r S;
				T bg;

				TVctr_r A_min;
				TVctr_r A_max;
				TVctr_r S_min;
				TVctr_r S_max;
		};

	}
#endif