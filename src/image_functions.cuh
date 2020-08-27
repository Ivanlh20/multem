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

#ifndef IMAGE_FUNCTIONS_H
#define IMAGE_FUNCTIONS_H

#include <thread>
#include <algorithm>
#include <deque>
#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "cpu_fcns.hpp"
#include "cgpu_classes.cuh"

namespace mt
{
	// hist
	template <class TVector>
	void hist(TVector &v, int nbins, TVector &y, TVector *x = nullptr)
	{
		using T = Value_type<TVector>;
		auto minmax_temp = std::minmax_element(v.begin(), v.end());
		T v_min = *(minmax_temp.first);
		T v_max = *(minmax_temp.second);
		T dv = (v_max-v_min)/nbins;

		thrust::fill(y.begin(), y.end(), 0);
		for(auto iv = 0; iv<v.size(); iv++)
		{
			int iy = min(static_cast<int>(floor((v[iv]-v_min)/dv)), nbins-1);
			y[iy]++;
		}

		if(x!= nullptr)
		{
			for(auto ix = 0; ix<nbins; ix++)
			{
				(*x)[ix] = v_min + ix*dv;
			}
		}
	}

	template <class TVector>
	Vector<int, e_host> hist(TVector &v, int nbins)
	{
		TVector y(nbins);
		hist(v, nbins, y);
		return y;
	}

	// Thresholding
	template <class TVector>
	Value_type<TVector> otsu_thr(TVector &v, int nbins)
	{
		using T = Value_type<TVector>;

		TVector v_rang(nbins);
		TVector v_hist(nbins);

		hist(v, nbins, v_hist, &v_rang);

		T sum_I = 0;
		for (auto ibins = 0; ibins<nbins; ibins++)
		{
			sum_I += v_rang[ibins]*v_hist[ibins];
		}

		T sum_B = 0;
		int w_B = 0;
		int w_F = 0;

		T var_Max = 0;
		T thr = 0;

		for (int ibins = 0; ibins<nbins; ibins++)
		{
			w_B += v_hist[ibins]; // Weight Background

			if (w_B == 0)
			{
				continue;
			}

			w_F = v.size() - w_B; // Weight Foreground
			if (w_F == 0)
			{
				break;
			}

			sum_B += v_rang[ibins]*v_hist[ibins];

			T m_B = sum_B/w_B; // Mean Background
			T m_F = (sum_I - sum_B)/w_F; // Mean Foreground

			// Calculate Between Class Variance
			auto var_Between = static_cast<T>(w_B)*static_cast<T>(w_F)*(m_B-m_F)*(m_B-m_F);

			// Check if new maximum found
			if (var_Between > var_Max) 
			{
				var_Max = var_Between;
				thr = v_rang[ibins];
			}
		}

		return thr;
	}

	// binarize
	template <class TVector>
	TVector binarize(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> thr)
	{
		using value_type = Value_type<TVector>;
		TVector v_o(v_i.size());

		auto thr_binarization = [&](const Range_2d &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::binarize<value_type>(thr));
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_binarization);

		return v_o;
	}

	// thresholding
	template <class TVector>
	TVector thresholding(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> thr)
	{
		using value_type = Value_type<TVector>;
		TVector v_o(v_i.size());

		auto thr_thring = [&](const Range_2d &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::thresholding<value_type>(thr));
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_thring);

		return v_o;
	}

	/******************************************************************************/
	// gray dilation
	template <class TVector>
	TVector morp_g_dilate(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr)
	{
		using T = Value_type<TVector>;

		TVector Im_o(Im_i.size());

		auto thr_dilate = [nkr](const Range_2d &range, TVector &Im_i, TVector &Im_o)
		{
			std::deque<int> wd;

			const int ny = range.iy_e-range.iy_0;
			const int wkr = 2*nkr+1;

			auto wedge_push = [ny](const int &ix, const int &iy, TVector &Im_i, std::deque<int> &wd)
			{
				while((!wd.empty()) && (Im_i[ix*ny+wd.back()]<= Im_i[ix*ny+iy]))
				{
					wd.pop_back();
				}
				wd.push_back(iy);
			};

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				int iy = 0;
				int iy_tk = 0;
				int iy_o = 0;

				for (; iy < nkr; iy++)
				{
					wedge_push(ix, iy, Im_i, wd);
				}

				for (; iy < ny; iy++, iy_o++)
				{
					wedge_push(ix, iy, Im_i, wd);

					Im_o[ix*ny+iy_o] = Im_i[ix*ny+wd.front()];
					if (iy + 1 >= wkr)
					{
						if (wd.front() == iy_tk)
						{
							wd.pop_front();
						}
						iy_tk++;
					}
				}

				for (; iy_o < ny; iy_o++)
				{
					Im_o[ix*ny+iy_o] = Im_i[ix*ny+wd.front()];
					if (wd.front() == iy_tk)
					{
						wd.pop_front();
					}
					iy_tk++;
				}
				wd.clear();
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_dilate, Im_i, Im_o);

		auto Im_t = Im_o;
		trs(stream, ny_i, nx_i, Im_t);

		stream.set_n_act_stream(ny_i);
		stream.set_grid(ny_i, nx_i);
		stream.exec(thr_dilate, Im_t, Im_o);

		trs(stream, nx_i, ny_i, Im_o);

		return Im_o;
	}

	// gray erosion
	template <class TVector>
	TVector morp_g_erode(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr)
	{
		using T = Value_type<TVector>;

		TVector Im_o(Im_i.size());

		auto thr_erode = [nkr](const Range_2d &range, TVector &Im_i, TVector &Im_o)
		{
			std::deque<int> wd;

			const int ny = range.iy_e-range.iy_0;
			const int wkr = 2*nkr+1;

			auto wedge_push = [ny](const int &ix, const int &iy, TVector &Im_i, std::deque<int> &wd)
			{
				while((!wd.empty()) && (Im_i[ix*ny+wd.back()]>= Im_i[ix*ny+iy]))
				{
					wd.pop_back();
				}
				wd.push_back(iy);
			};

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				int iy = 0;
				int iy_tk = 0;
				int iy_o = 0;

				for (; iy < nkr; iy++)
				{
					wedge_push(ix, iy, Im_i, wd);
				}

				for (; iy < ny; iy++, iy_o++)
				{
					wedge_push(ix, iy, Im_i, wd);

					Im_o[ix*ny+iy_o] = Im_i[ix*ny+wd.front()];
					if (iy + 1 >= wkr)
					{
						if (wd.front() == iy_tk)
						{
							wd.pop_front();
						}
						iy_tk++;
					}
				}

				for (; iy_o < ny; iy_o++)
				{
					Im_o[ix*ny+iy_o] = Im_i[ix*ny+wd.front()];
					if (wd.front() == iy_tk)
					{
						wd.pop_front();
					}
					iy_tk++;
				}
				wd.clear();
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_erode, Im_i, Im_o);

		auto Im_t = Im_o;
		trs(stream, ny_i, nx_i, Im_t);

		stream.set_n_act_stream(ny_i);
		stream.set_grid(ny_i, nx_i);
		stream.exec(thr_erode, Im_t, Im_o);

		trs(stream, nx_i, ny_i, Im_o);

		return Im_o;
	}

	// gray opening
	template <class TVector>
	TVector morp_g_open(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr)
	{
		auto Im = morp_g_erode(stream, ny_i, nx_i, Im_i, nkr);
		return morp_g_dilate(stream, ny_i, nx_i, Im, nkr);
	}

	// gray closing
	template <class TVector>
	TVector morp_g_close(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr)
	{
		auto Im = morp_g_dilate(stream, ny_i, nx_i, Im_i, nkr);
		return morp_g_erode(stream, ny_i, nx_i, Im, nkr);
	}

	// gray tophat
	template <class TVector>
	TVector morp_g_tophat(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr)
	{
		auto Im = morp_g_open(stream, ny_i, nx_i, Im_i, nkr);
		add_scale(stream, 1, Im_i, -1, Im, Im);
		return Im;
	}
	
	/******************************************************************************/
	// forward anscombe transform
	template <class TVector>
	TVector anscombe_forward(Stream<e_host> &stream, TVector &v_i)
	{
		using value_type = Value_type<TVector>;
		TVector v_o(v_i.size());

		auto thr_anscombe_forward = [&](const Range_2d &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::anscombe_forward<value_type>());
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_anscombe_forward);

		return v_o;
	}

	// forward anscombe transform
	template <class TVector>
	TVector anscombe_inverse(Stream<e_host> &stream, TVector &v_i)
	{
		using value_type = Value_type<TVector>;
		TVector v_o(v_i.size());

		auto thr_anscombe_inverse = [&](const Range_2d &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::anscombe_inverse<value_type>());
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_anscombe_inverse);

		return v_o;
	}

	/******************************************************************************/
	// wiener filter 1d
	template <class TVector>
	TVector ftr_wiener_1d(Stream<e_host> &stream, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_o(Im_i.size());

		int nk0 = -nkr;
		int nke = nkr+1;
		int nx_i = Im_i.size();

		TVector Im_mean(nx_i);
		TVector Im_var(nx_i);

		T v2 = 0;
		auto thr_mean_var = [&](const Range_2d &range)
		{
			T v2_partial = 0;
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				int ix_0 = max(ixy+nk0, 0);
				int ixe = min(ixy+nke, nx_i);

				T x_mean = 0;
				T x_var = 0;
				for (auto ix = ix_0; ix < ixe; ix++)
				{
					T x = Im_i[ix];
					x_mean += x;
					x_var += x*x;
				}
				x_mean = x_mean/(ixe-ix_0);
				x_var = x_var/(ixe-ix_0) - x_mean*x_mean;

				Im_mean[ixy] = x_mean;
				Im_var[ixy] = x_var;

				v2_partial += x_var;
			}

			stream.stream_mutex.lock();
			v2 += v2_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, 1);
		stream.exec(thr_mean_var);

		v2 /= nx_i;

		auto thr_ftr_wiener = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_o[ixy] = Im_mean[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_mean[ixy])/::fmax(Im_var[ixy], v2);
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, 1);
		stream.exec(thr_ftr_wiener);

		return Im_o;
	}
	
	// wiener filter by rows
	template <class TGrid, class TVector>
	TVector ftr_wiener_2d_br(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_o(Im_i.size());

		int nk0 = -nkr;
		int nke = nkr+1;

		TVector Im_mean(grid_2d.nxy());
		TVector Im_var(grid_2d.nxy());

		auto krn_mean_var = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_mean, TVector &Im_var)
		{
			int ix_0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, grid_2d.nx);

			T x_mean = 0;
			T x_var = 0;
			for (auto ix = ix_0; ix < ixe; ix++)
			{
				T x = Im_i[grid_2d.ind_col(ix, iy_i)];
				x_mean += x;
				x_var += x*x;
			}
			x_mean = x_mean/(ixe-ix_0);
			x_var = x_var/(ixe-ix_0) - x_mean*x_mean;

			int ixy = grid_2d.ind_col(ix_i, iy_i);
			Im_mean[ixy] = x_mean;
			Im_var[ixy] = x_var;
		};

		auto thr_mean_var = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_mean_var, Im_i, Im_mean, Im_var);
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_mean_var);

		auto thr_ftr_wiener = [&](const Range_2d &range)
		{
			for(auto iy = range.ixy_0; iy < range.ixy_e; iy++)
			{
				T v2 = 0;
				for(auto ix = 0; ix < grid_2d.nx; ix++)
				{
					v2 += Im_var[grid_2d.ind_col(ix, iy)];
				}
				v2 /= grid_2d.nx;

				for(auto ix = 0; ix < grid_2d.nx; ix++)
				{
					auto ixy = grid_2d.ind_col(ix, iy);
					Im_o[ixy] = Im_mean[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_mean[ixy])/::fmax(Im_var[ixy], v2);
				}
			}
		};

		stream.set_n_act_stream(grid_2d.ny);
		stream.set_grid(1, grid_2d.ny);
		stream.exec(thr_ftr_wiener);

		return Im_o;
	}

	// wiener filter 2d
	template <class TGrid, class TVector>
	TVector ftr_wiener_2d(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_mean(Im_i.size());
		TVector Im_var(Im_i.size());
		TVector Im_t(Im_i.size());

		// Kahan summation algorithm
		// https:// en.wikipedia.org/wiki/Kahan_summation_algorithm

		auto thr_sum = [=](const Range_2d &range, TVector &Im_i, TVector &Im_sum)
		{
			const int ny = range.iy_e-range.iy_0;
			const int wkr = 2*nkr+1;

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				int iy = 0;
				int iy_tk = 0;
				int iy_o = 0;

				T sum_v = 0;
				T sum_ee = 0; 
				for (; iy < nkr; iy++)
				{
					T v = Im_i[ix*ny+iy];
					host_device_detail::kh_sum(sum_v, v, sum_ee);
				}

				for (; iy < ny; iy++, iy_o++)
				{
					T v = Im_i[ix*ny+iy];
					host_device_detail::kh_sum(sum_v, v, sum_ee);

					Im_sum[ix*ny+iy_o] = sum_v;
					if (iy + 1 >= wkr)
					{
						T v = Im_i[ix*ny+iy_tk];
						host_device_detail::kh_sum(sum_v, -v, sum_ee);
						iy_tk++;
					}
				}

				for (; iy_o < ny; iy_o++)
				{
					Im_sum[ix*ny+iy_o] = sum_v;

					T v = Im_i[ix*ny+iy_tk];
					host_device_detail::kh_sum(sum_v, -v, sum_ee);
					iy_tk++;
				}
			}
		};

		auto thr_sum2 = [=](const Range_2d &range, TVector &Im_i, TVector &Im_sum2)
		{
			const int ny = range.iy_e-range.iy_0;
			const int wkr = 2*nkr+1;

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				int iy = 0;
				int iy_tk = 0;
				int iy_o = 0;

				T sum_v2 = 0;
				T sum_ee = 0; 
				for (; iy < nkr; iy++)
				{
					T v = Im_i[ix*ny+iy];
					host_device_detail::kh_sum(sum_v2, v*v, sum_ee);
				}

				for (; iy < ny; iy++, iy_o++)
				{
					T v = Im_i[ix*ny+iy];
					host_device_detail::kh_sum(sum_v2, v*v, sum_ee);

					Im_sum2[ix*ny+iy_o] = sum_v2;
					if (iy + 1 >= wkr)
					{
						T v = Im_i[ix*ny+iy_tk];
						host_device_detail::kh_sum(sum_v2, -v*v, sum_ee);
						iy_tk++;
					}
				}

				for (; iy_o < ny; iy_o++)
				{
					Im_sum2[ix*ny+iy_o] = sum_v2;

					T v = Im_i[ix*ny+iy_tk];
					host_device_detail::kh_sum(sum_v2, -v*v, sum_ee);
					iy_tk++;
				}
			}
		};

		// sum
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_sum, Im_i, Im_mean);

		Im_t = Im_mean;
		trs(stream, grid_2d.ny, grid_2d.nx, Im_t);

		stream.set_n_act_stream(grid_2d.ny);
		stream.set_grid(grid_2d.ny, grid_2d.nx);
		stream.exec(thr_sum, Im_t, Im_mean);

		trs(stream, grid_2d.nx, grid_2d.ny, Im_mean);

		// sum^2
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_sum2, Im_i, Im_var);

		Im_t = Im_var;
		trs(stream, grid_2d.ny, grid_2d.nx, Im_t);

		stream.set_n_act_stream(grid_2d.ny);
		stream.set_grid(grid_2d.ny, grid_2d.nx);
		stream.exec(thr_sum, Im_t, Im_var);

		trs(stream, grid_2d.nx, grid_2d.ny, Im_var);

		for(auto ix=0; ix<grid_2d.nx; ix++)
		{
			for(auto iy=0; iy<grid_2d.ny; iy++)
			{
				int nx = min(ix+nkr+1, grid_2d.nx)-max(ix-nkr, 0);
				int ny = min(iy+nkr+1, grid_2d.ny)-max(iy-nkr, 0);

				int ixy = grid_2d.ind_col(ix, iy);

				T x_mean = Im_mean[ixy]/(nx*ny);
				T x_var = Im_var[ixy]/(nx*ny) - x_mean*x_mean;

				Im_mean[ixy] = x_mean;
				Im_var[ixy] = x_var;
			}
		}

		auto v2 = mt::mean(stream, Im_var);

		auto thr_ftr_wiener = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_t[ixy] = Im_mean[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_mean[ixy])/::fmax(Im_var[ixy], v2);
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_ftr_wiener);

		return Im_t;
	}

	/******************************************************************************/	
	// median wiener2
	template <class TVector>
	void ftr_mwiener(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nk0 = -nkr;
		int nke = nkr+1;
		Vector<T, e_host> Im_median(nx_i*ny_i);
		Vector<T, e_host> Im_var(nx_i*ny_i);

		auto krn_median_var = [&](const int &ix_i, const int &iy_i, TVector &Im_i, Vector<T, e_host> &Im_median, Vector<T, e_host> &Im_var)
		{
			int ix_0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nx_i);

			int iy_0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, ny_i);

			T nxy = (ixe-ix_0)*(iye-iy_0);
			Vector<T, e_host> v(nxy);
			int iv = 0;
			T x_mean = 0;
			T x_var = 0;
			for (auto ix = ix_0; ix < ixe; ix++)
			{
				for (auto iy = iy_0; iy < iye; iy++)
				{
					T x = Im_i[ix*ny_i+iy];
					v[iv++] = x;
					x_mean += x;
					x_var += x*x;
				}
			}
			auto median = v.begin() + nxy/2;
			std::nth_element(v.begin(), median, v.end());

			x_mean = x_mean/nxy;
			x_var = x_var/nxy - x_mean*x_mean;

			int ixy = ix_i*ny_i+iy_i;
			Im_median[ixy] = *median;
			Im_var[ixy] = x_var;
		};

		auto thr_median_var = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_median_var, Im_i, Im_median, Im_var);
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_median_var);

		auto v2 = mt::mean(stream, Im_var);

		auto thr_ftr_mwiener = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_o[ixy] = Im_median[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_median[ixy])/::fmax(Im_var[ixy], v2);
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_ftr_mwiener);
	}

	/******************************************************************************/
	// median filter 1d
	template <class TVector>
	TVector ftr_median_1d(Stream<e_host> &stream, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_o(Im_i.size());

		int nk0 = -nkr;
		int nke = nkr+1;
		int nx_i = Im_i.size();

		auto thr_ftr_median = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				int ix_0 = max(ixy+nk0, 0);
				int ixe = min(ixy+nke, nx_i);

				TVector v;
				v.reserve(ixe-ix_0);

				for (auto ix = ix_0; ix < ixe; ix++)
				{
					v.push_back(Im_i[ix]);
				}
				auto median = v.begin() + v.size()/2;
				std::nth_element(v.begin(), median, v.end());
				Im_o[ixy] = *median;
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, 1);
		stream.exec(thr_ftr_median);

		return Im_o;
	}

	// median filter by rows
	template <class TGrid, class TVector>
	TVector ftr_median_2d_br(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_o(grid_2d.nxy());

		int nk0 = -nkr;
		int nke = nkr+1;

		auto krn_ftr_median = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		{
			int ix_0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, grid_2d.nx);

			TVector v;
			v.reserve(ixe-ix_0);

			for (auto ix = ix_0; ix < ixe; ix++)
			{
				v.push_back(Im_i[grid_2d.ind_col(ix, iy_i)]);
			}
			auto median = v.begin() + v.size()/2;
			std::nth_element(v.begin(), median, v.end());
			Im_o[grid_2d.ind_col(ix_i, iy_i)] = *median;
		};

		auto thr_ftr_median = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_ftr_median, Im_i, Im_o);
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_ftr_median);

		return Im_o;
	}

	// median filter 2d
	template <class TGrid, class TVector>
	TVector ftr_median_2d(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr)
	{
		if(nkr<=0)
		{
			return Im_i;
		}

		using T = Value_type<TVector>;
		TVector Im_o(grid_2d.nxy());

		int nk0 = -nkr;
		int nke = nkr+1;
		auto R2_max = pow(nkr*grid_2d.dRx, 2) + Epsilon<T>::abs;

		auto krn_ftr_median = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		{
			auto x_i = grid_2d.Rx(ix_i);
			auto y_i = grid_2d.Ry(iy_i);

			int ix_0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, grid_2d.nx);

			int iy_0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, grid_2d.ny);

			TVector v;
			v.reserve((ixe-ix_0)*(iye-iy_0));

			for (auto ix = ix_0; ix < ixe; ix++)
			{
				for (auto iy = iy_0; iy < iye; iy++)
				{
					auto R2 = grid_2d.R2(ix, iy, x_i, y_i);
					if(R2 < R2_max)
					{
						v.push_back(Im_i[grid_2d.ind_col(ix, iy)]);
					}
				}
			}
			auto median = v.begin() + v.size()/2;
			std::nth_element(v.begin(), median, v.end());
			Im_o[grid_2d.ind_col(ix_i, iy_i)] = *median;
		};

		auto thr_ftr_median = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_ftr_median, Im_i, Im_o);
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_ftr_median);

		return Im_o;
	}

	/******************************************************************************/
	// den poiss 
	template <class TVector>
	TVector ftr_poiss_dnois_1d(Stream<e_host> &stream, TVector &Im_i, int nkr_w, int nkr_m)
	{
		if((nkr_w == 0) && (nkr_m == 0))
		{
		 return Im_i;
		}

		auto Im = anscombe_forward(stream, Im_i);
		if(nkr_w>0)
		{
			Im = ftr_wiener_1d(stream, Im, nkr_w);
		}
		if(nkr_m>0)
		{
			Im = ftr_median_1d(stream, Im, nkr_m);
		}
		Im = anscombe_inverse(stream, Im);

		return Im;
	}
	
	// den poiss 
	template <class TGrid, class TVector>
	TVector ftr_poiss_dnois_2d_br(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr_w, int nkr_m)
	{
		if((nkr_w == 0) && (nkr_m == 0))
		{
			return Im_i;
		}

		auto Im = anscombe_forward(stream, Im_i);
		if(nkr_w>0)
		{
			Im = ftr_wiener_2d_br(stream, grid_2d, Im, nkr_w);
		}
		if(nkr_m>0)
		{
			Im = ftr_median_2d_br(stream, grid_2d, Im, nkr_m);
		}
		Im = anscombe_inverse(stream, Im);

		return Im;
	}

	// den poiss 
	template <class TGrid, class TVector>
	TVector ftr_poiss_dnois_2d(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr_w, int nkr_m)
	{
		if((nkr_w == 0) && (nkr_m == 0))
		{
			return Im_i;
		}

		auto Im = anscombe_forward(stream, Im_i);
		if(nkr_w>0)
		{
			Im = ftr_wiener_2d(stream, grid_2d, Im, nkr_w);
		}
		if(nkr_m>0)
		{
			Im = ftr_median_2d(stream, grid_2d, Im, nkr_m);
		}
		Im = anscombe_inverse(stream, Im);

		return Im;
	}

	/******************************************************************************/

	// get peak signal to noise ratio PSNR 
	template <class TGrid, class TVector>
	Value_type<TVector> get_PSNR(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr_w, int nkr_m)
	{
		auto Im_s = ftr_poiss_dnois_2d(stream, grid_2d, Im_i, nkr_w, nkr_m);

		auto var_signal = variance(stream, Im_s);

		add_scale(stream, 1, Im_i, -1, Im_s, Im_s);
		auto var_no = variance(stream, Im_s);

		// peak signal to noise ratio
		return var_no/var_signal;
	}

	// get peak signal to noise ratio PSNR 
	template <class TVector>
	Value_type<TVector> get_PSNR(Stream<e_host> &stream, TVector &Im_i, TVector &Im_d)
	{
		auto var_signal = variance(stream, Im_d);
		TVector Im_n(Im_d.size());
		add_scale(stream, 1, Im_i, -1, Im_d, Im_n);
		auto var_no = variance(stream, Im_n);

		// peak signal to noise ratio
		return var_no/var_signal;
	}

	// scale_image
	template <class TVector>
	TVector scale_image_mean(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, Value_type<TVector> shrink_factor, int &ny_o, int &nx_o)
	{
		using T = Value_type<TVector>;
		TVector Im_o;

		int nkr = max(int(1.0/(2.0*shrink_factor)), 1); 
		int nk0 = -nkr;
		int nke = nkr+1;

		nx_o = max(int(floor(nx_i*shrink_factor)), 1);
		ny_o = max(int(floor(ny_i*shrink_factor)), 1);

		if((nx_i == nx_o) && (ny_i == ny_o))
		{
			Im_o.assign(Im_i.begin(), Im_i.end());
			return Im_o;
		}
		Im_o.resize(nx_o*ny_o);

		auto krn_scale_image = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		{
			auto ix_t = static_cast<int>(ix_i/shrink_factor);
			auto iy_t = static_cast<int>(iy_i/shrink_factor);

			int ix_0 = max(ix_t+nk0, 0);
			int ixe = min(ix_t+nke, nx_i);

			int iy_0 = max(iy_t+nk0, 0);
			int iye = min(iy_t+nke, ny_i);

			T sum = 0;
			for (auto ix = ix_0; ix < ixe; ix++)
			{
				for (auto iy = iy_0; iy < iye; iy++)
				{
					sum += Im_i[ix*ny_i+iy];
				}
			}

			Im_o[ix_i*ny_o+iy_i] = sum/((ixe-ix_0)*(iye-iy_0));
		};

		auto thr_scale_image = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_scale_image, Im_i, Im_o);
		};

		stream.set_n_act_stream(nx_o);
		stream.set_grid(nx_o, ny_o);
		stream.exec(thr_scale_image);

		return Im_o;
	}

	// copy image
	template <class TVector>
	TVector copy_image(Stream<e_host> &stream, int ny_src, int nx_src, TVector &Im_src, int iy_0, int ix_0, int iye, int ixe, int ny_dst, int nx_dst)
	{
		TVector Im_dst(nx_dst*ny_dst);

		int nx = min(nx_src, nx_dst);
		ixe = (ixe<nx)?ixe+1:nx;
		nx = ixe-ix_0;

		int ny = min(ny_src, ny_dst);
		iye = (iye<ny)?iye+1:ny;
		ny = iye-iy_0;

		auto thr_copy_image = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					Im_dst[ix*ny_dst+iy] = Im_src[(ix+ix_0)*ny_src+(iy+iy_0)];
				}
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, ny);
		stream.exec(thr_copy_image);

		return Im_dst;
	}

	// extract shape
	template <class TGrid, class TVector>
	TVector extract_shape(Stream<e_host> &stream, FFT<Value_type<TGrid>, e_host> &fft_2d, TGrid grid_2d, TVector &Im_i)
	{
		using T = Value_type<TVector>;

		T sigma = 1.0;
		T shrink = 0.5;
		T dR = grid_2d.dR_min();

		int nx_d = 1;
		int ny_d = 1;

		// scale_image image
		auto Im_d = scale_image_mean(stream, grid_2d.ny, grid_2d.nx, Im_i, shrink, ny_d, nx_d);

		// Otsu thr
		int nbins = 256;
		auto thr = mt::otsu_thr(Im_d, nbins);

		// binarize
		Im_d = binarize(stream, Im_d, thr);

		// copy to the first quadrant
		auto Im_o = copy_image(stream, ny_d, nx_d, Im_d, 0, 0, ny_d-1, nx_d-1, grid_2d.ny, grid_2d.nx);

		// dilate binary image
		int nkr = static_cast<int>(0.5*sigma/dR);
		Im_o = morp_g_dilate(stream, grid_2d.ny, grid_2d.nx, Im_o, nkr);

		// gaussian convolution
		Gauss_Cv_2d<T, e_host> gauss_cv_2d(&stream, &fft_2d, grid_2d);
		gauss_cv_2d(sigma, Im_o);

		// binarize
		Im_o = binarize(stream, Im_o, 0.01);

		// copy to the first quadrant
		Im_d = copy_image(stream, grid_2d.ny, grid_2d.nx, Im_o, 0, 0, ny_d-1, nx_d-1, ny_d, nx_d);

		// scale image
		Im_o = scale_image_mean(stream, ny_d, nx_d, Im_d, 1.0/shrink, grid_2d.ny, grid_2d.nx);

		// binarize
		Im_o = binarize(stream, Im_o, 0.5);

		return Im_o;
	}

} // namespace mt

#endif