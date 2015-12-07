/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft2.cuh"
#include "host_device_functions.cuh"
#include "host_functions.hpp"

namespace multem
{
	// histogram
	template<class TVector>
	void histogram(TVector &v, int nbins, Vector<int, e_host> &y, TVector *x =nullptr)
	{
		using T = Value_type<TVector>;
		auto minmax_temp = std::minmax_element(v.begin(), v.end());
		T v_min = *(minmax_temp.first);
		T v_max = *(minmax_temp.second);
		T dv = (v_max-v_min)/nbins;

		thrust::fill(y.begin(), y.end(), 0);
		for(auto iv =0; iv<v.size(); iv++)
		{
			int iy = min(static_cast<int>(floor((v[iv]-v_min)/dv)), nbins-1);
			y[iy]++;
		}

		if(x!=nullptr)
		{
			for(auto ix =0; ix<nbins; ix++)
			{
				(*x)[ix] = v_min + ix*dv;
			}
		}
	}

	template<class TVector>
	Vector<int, e_host> histogram(TVector &v, int nbins)
	{
		Vector<int, e_host> y(nbins);
		histogram(v, nbins, y);
		return y;
	}

	// Thresholding
	template<class TVector>
	Value_type<TVector> otsu_threshold(TVector &v, int nbins)
	{
		using T = Value_type<TVector>;

		TVector v_rang(nbins);
		Vector<int, e_host> v_hist(nbins);

		histogram(v, nbins, v_hist, &v_rang);

		T sum_I = 0;
		for (auto ibins =0; ibins<nbins; ibins++)
		{
			sum_I += v_rang[ibins]*v_hist[ibins];
		}

		T sum_B = 0;
		int w_B = 0;
		int w_F = 0;

		T var_Max = 0;
		T threshold = 0;

		for (int ibins =0; ibins<nbins; ibins++)
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
				threshold = v_rang[ibins];
			}
		}

		return threshold;
	}

	// binarization
	template<class TVector>
	void binary(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> threshold, TVector &v_o)
	{
		using value_type = Value_type<TVector>;
		auto thr_binary = [&](const Range &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::binarization<value_type>(threshold));
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_binary);
	}

	template<class TVector>
	TVector binary(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> threshold)
	{
		TVector v_o(v_i.size());
		binary(stream, v_i, threshold, v_o);
		return v_o;
	}

	// morphology dilation
	template<class TVector>
	void morf_dilation(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nk0 = -nkr;
		int nke = nkr+1;

		auto krn_dilate = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		{
			int ix0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nx_i);

			int iy0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, ny_i);

			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
				{
					if(Im_i[ix*ny_i+iy]>0.5)
					{	
						Im_o[ix_i*ny_i+iy_i] = 1;
						return;
					}
				}
			}
			Im_o[ix_i*ny_i+iy_i] = 0;
		};

		auto thr_dilate = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_dilate, Im_i, Im_o);
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_dilate);
	}

	template<class TVector>
	TVector morf_dilation(Stream<e_host> &stream, int ny, int nx, TVector &Im_i, int nkr)
	{
		TVector Im_o(Im_i.size());
		morf_dilation(stream, ny, nx, Im_i, nkr, Im_o);
		return Im_o;
	}

	// scale_image
	template<class TVector>
	void scale_image_mean(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, Value_type<TVector> shrink_factor, int &ny_o, int &nx_o, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nkr = max(int(1.0/(2.0*shrink_factor)), 1); 
		int nk0 = -nkr;
		int nke = nkr+1;

		nx_o = max(int(floor(nx_i*shrink_factor)), 1);
		ny_o = max(int(floor(ny_i*shrink_factor)), 1);

		if((nx_i == nx_o)&&(ny_i == ny_o))
		{
			Im_o.assign(Im_i.begin(), Im_i.end());
			return;
		}
		Im_o.resize(nx_o*ny_o);

		auto krn_scale_image = [&](const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		{
			auto ix_t = static_cast<int>(ix_i/shrink_factor);
			auto iy_t = static_cast<int>(iy_i/shrink_factor);

			int ix0 = max(ix_t+nk0, 0);
			int ixe = min(ix_t+nke, nx_i);

			int iy0 = max(iy_t+nk0, 0);
			int iye = min(iy_t+nke, ny_i);

			T sum = 0;
			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
				{
					sum += Im_i[ix*ny_i+iy];
				}
			}

			Im_o[ix_i*ny_o+iy_i] = sum/((ixe-ix0)*(iye-iy0));
		};

		auto thr_scale_image = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_scale_image, Im_i, Im_o);
		};

		stream.set_n_act_stream(nx_o);
		stream.set_grid(nx_o, ny_o);
		stream.exec(thr_scale_image);
	}

	// wiener2
	template<class TVector>
	void filter_wiener(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nk0 = -nkr;
		int nke = nkr+1;
		Vector<T, e_host> Im_mean(nx_i*ny_i);
		Vector<T, e_host> Im_var(nx_i*ny_i);

		auto krn_mean_var = [&](const int &ix_i, const int &iy_i, TVector &Im_i, Vector<T, e_host> &Im_mean, Vector<T, e_host> &Im_var)
		{
			int ix0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nx_i);

			int iy0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, ny_i);

			T nxy = (ixe-ix0)*(iye-iy0);
			T x_mean = 0;
			T x_var = 0;
			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
				{
					T x = Im_i[ix*ny_i+iy];
					x_mean += x;
					x_var += x*x;
				}
			}
			x_mean = x_mean/nxy;
			x_var = x_var/nxy - x_mean*x_mean;

			int ixy = ix_i*ny_i+iy_i;
			Im_mean[ixy] = x_mean;
			Im_var[ixy] = x_var;
		};

		auto thr_mean_var = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_mean_var, Im_i, Im_mean, Im_var);
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_mean_var);

		auto v2 = multem::mean(stream, Im_var);

		auto thr_filter_wiener = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_o[ixy] = Im_mean[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_mean[ixy])/::fmax(Im_var[ixy], v2);
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_filter_wiener);
	}
	
	// median wiener2
	template<class TVector>
	void filter_mwiener(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nk0 = -nkr;
		int nke = nkr+1;
		Vector<T, e_host> Im_median(nx_i*ny_i);
		Vector<T, e_host> Im_var(nx_i*ny_i);

		auto krn_median_var = [&](const int &ix_i, const int &iy_i, TVector &Im_i, Vector<T, e_host> &Im_median, Vector<T, e_host> &Im_var)
		{
			int ix0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nx_i);

			int iy0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, ny_i);

			T nxy = (ixe-ix0)*(iye-iy0);
			Vector<T, e_host> v(nxy);
			int iv = 0;
			T x_mean = 0;
			T x_var = 0;
			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
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

		auto thr_median_var = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_median_var, Im_i, Im_median, Im_var);
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_median_var);

		auto v2 = multem::mean(stream, Im_var);

		auto thr_filter_mwiener = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_o[ixy] = Im_median[ixy] + ::fmax(T(0), Im_var[ixy]-v2)*(Im_i[ixy]-Im_median[ixy])/::fmax(Im_var[ixy], v2);
			}
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_filter_mwiener);
	}

	template<class TVector>
	void filter_median(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr, TVector &Im_o)
	{
		using T = Value_type<TVector>;
		int nk0 = -nkr;
		int nke = nkr+1;

		Vector<T, e_host> Im(nx_i*ny_i);
		auto krn_filter_median = [&](const int &ix_i, const int &iy_i, TVector &Im_i, Vector<T, e_host> &Im_o)
		{
			int ix0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nx_i);

			int iy0 = max(iy_i+nk0, 0);
			int iye = min(iy_i+nke, ny_i);

			T nxy = (ixe-ix0)*(iye-iy0);
			Vector<T, e_host> v(nxy);
			int iv = 0;
			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
				{
					v[iv++] = Im_i[ix*ny_i+iy];
				}
			}
			auto median = v.begin() + nxy/2;
			std::nth_element(v.begin(), median, v.end());
			Im_o[ix_i*ny_i+iy_i] = *median;
		};

		auto thr_filter_median = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_filter_median, Im_i, Im);
		};

		stream.set_n_act_stream(nx_i);
		stream.set_grid(nx_i, ny_i);
		stream.exec(thr_filter_median);
		Im_o.assign(Im.begin(), Im.end());
	}

	// forward anscombe transform
	template<class TVector>
	void anscombe_forward(Stream<e_host> &stream, TVector &v_i, TVector &v_o)
	{
		using value_type = Value_type<TVector>;
		auto thr_anscombe_forward = [&](const Range &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::anscombe_forward<value_type>());
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_anscombe_forward);
	}

	// forward anscombe transform
	template<class TVector>
	void anscombe_inverse(Stream<e_host> &stream, TVector &v_i, TVector &v_o)
	{
		using value_type = Value_type<TVector>;
		auto thr_anscombe_inverse = [&](const Range &range)
		{
			thrust::transform(v_i.begin()+range.ixy_0, v_i.begin()+range.ixy_e, 
			v_o.begin()+range.ixy_0, functor::anscombe_inverse<value_type>());
		};

		stream.set_n_act_stream(v_i.size());
		stream.set_grid(1, v_i.size());
		stream.exec(thr_anscombe_inverse);
	}

	// copy image
	template<class TVector>
	void copy_image(Stream<e_host> &stream, int ny_src, int nx_src, TVector &Im_src, int iy0, int ix0, int iye, int ixe, int ny_dst, int nx_dst, TVector &Im_dst)
	{
		int nx = min(nx_src, nx_dst);
		ixe = (ixe<nx)?ixe+1:nx;
		nx = ixe-ix0;

		int ny = min(ny_src, ny_dst);
		iye = (iye<ny)?iye+1:ny;
		ny = iye-iy0;

		auto krn_copy_image = [&](const int &ix, const int &iy, TVector &Im_src, TVector &Im_dst)
		{

			Im_dst[ix*ny_dst+iy] = Im_src[(ix+ix0)*ny_src+(iy+iy0)];
		};

		auto thr_copy_image = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_copy_image, Im_src, Im_dst);
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, ny);
		stream.exec(thr_copy_image);
	}

	// extract shape
	template<class TGrid, class TVector>
	void extract_shape(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid grid, TVector &Im_i, TVector &Im_o)
	{
		using T = Value_type<TVector>;

		T sigma = 1.0;
		T shrink = 0.5;
		T dR = grid.dR_min();

		int nx_d = 1;
		int ny_d = 1;
		TVector Im_d;

		// scale_image image
		scale_image_mean(stream, grid.ny, grid.nx, Im_i, shrink, ny_d, nx_d, Im_d);

		// Otsu threshold
		int nbins = 256;
		T threshold = multem::otsu_threshold(Im_d, nbins);

		// binarization
		binary(stream, Im_d, threshold, Im_d);

		// copy to the first quadrant
		copy_image(stream, ny_d, nx_d, Im_d, 0, 0, ny_d-1, nx_d-1, grid.ny, grid.nx, Im_o);

		// dilate binary image
		int nkr = static_cast<int>(0.5*sigma/dR);
		Im_o = morf_dilation(stream, grid.ny, grid.nx, Im_o, nkr);

		// shift Image
		multem::fft2_shift(stream, grid, Im_o);

		// copy real matrix to complex matrix
		Vector<complex<T>, e_host> Im_c(Im_o.size());
		for(auto ixy =0; ixy<Im_c.size(); ixy++)
		{
			Im_c[ixy] = complex<T>(Im_o[ixy], 0);
		}

		// gaussian convolution
		grid.set_input_data(grid.nx, grid.ny, 2*grid.lx, 2*grid.ly, 0.5, false, true);
		gaussian_convolution(stream, fft2, grid, sigma, Im_c);

		// shift Image
		multem::fft2_shift(stream, grid, Im_c);

		// binarization
		for(auto ixy =0; ixy<Im_o.size(); ixy++)
		{
			Im_o[ixy] = (Im_c[ixy].real()<0.01)?0:1;
		}

		// copy to the first quadrant
		copy_image(stream, grid.ny, grid.nx, Im_o, 0, 0, ny_d-1, nx_d-1, ny_d, nx_d, Im_d);

		// scale image
		scale_image_mean(stream, ny_d, nx_d, Im_d, 1.0/shrink, grid.ny, grid.nx, Im_o);

		// binarization
		binary(stream, Im_o, 0.5, Im_o);
	}

} // namespace multem

#endif