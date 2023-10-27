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

#include "opt_flow.h"

namespace mt
{
	template <class T, eDev Dev>
	Opt_Flow <T, Dev>::Opt_Flow(): stream(nullptr) {}

	template <class T, eDev Dev>
	Opt_Flow <T, Dev>::Opt_Flow(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
	{
		set_in_data(stream_i, fft_2d_i, grid_2d_i);
	}

	template <class T, eDev Dev>
	inline
	void Opt_Flow <T, Dev>::set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
	{
		stream = stream_i;
		grid_2d = grid_2d_i;

		intrpl_2d.set_in_data(stream, grid_2d, grid_2d);

		gauss_cv_2d.set_in_data(stream, fft_2d_i, grid_2d);

		v_x.resize(grid_2d.size());
		v_y.resize(grid_2d.size());

		Rx.resize(grid_2d.size());
		Ry.resize(grid_2d.size());
		M.resize(grid_2d.size());

	}

	template <class T, eDev Dev>
	void Opt_Flow <T, Dev>::operator()(TVctr& M_s, TVctr& M_m, T alpha, T sigma, dt_int32 n_iter, TVctr& v_xt, TVctr& v_yt)
	{
		// create rectangular grid
		set_regular_grid(Rx, Ry);

		// set initial optical flow
		v_x = v_xt;
		v_y = v_yt;

		for(auto iter = 0; iter < n_iter; iter++)
		{
			// create new grid
			mt::add(*stream, v_x, Rx);
			mt::add(*stream, v_y, Ry);

			// resample distored image in a new grid
			intrpl_2d(M_m, Rx, Ry, M);

			// calculate optical flow
			fcn_opt_flow(M_s, M, alpha, v_x, v_y);

			// regularization based on convolution
			if (fcn_is_nzero(sigma))
			{
				gauss_cv_2d(sigma, v_x);
				gauss_cv_2d(sigma, v_y);
			}

			// add optical flow
			mt::add(*stream, v_x, v_xt);
			mt::add(*stream, v_y, v_yt);

		}
	}

	template <class T, eDev Dev>
	void Opt_Flow <T, Dev>::set_fft_plan()
	{
		gauss_cv_2d.set_fft_plan();
	}

	template <class T, eDev Dev>
	void Opt_Flow <T, Dev>::cleanup()
	{
		gauss_cv_2d.cleanup();
	}

	template <class T, eDev Dev>
	void Opt_Flow <T, Dev>::set_regular_grid(TVctr& Rx, TVctr& Ry)
	{
		Vctr<T, edev_cpu> Rx_h;
		Vctr<T, edev_cpu> Ry_h;

		Rx_h.reserve(grid_2d.size());
		Ry_h.reserve(grid_2d.size());

		for(auto ix = 0; ix < grid_2d.nx; ix++)
		{
			for(auto iy = 0; iy < grid_2d.ny; iy++)
			{
				Rx_h.push_back(grid_2d.rx(ix));
				Ry_h.push_back(grid_2d.ry(iy));
			}
		}

		thrust::copy(Rx_h.begin(), Rx_h.end(), Rx.begin());
		thrust::copy(Ry_h.begin(), Ry_h.end(), Ry.begin());
	}

	/***************************************** cpu *****************************************/
	template <class T, eDev Dev>
	template <eDev devn = Dev>
	enable_if_edev_cpu<devn, void>
	Opt_Flow <T, Dev>::fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y)
	{
		stream->set_n_stream_act(grid_2d.nx);
		stream->set_grid(grid_2d.nx, grid_2d.ny);
		stream->exec_2d(detail_cgpu::fcn_opt_flow<Grid_2d<T>, TVctr>, grid_2d, M_s, M_m, alpha, v_x, v_y);
	}

	/**********************Device**********************/
#ifdef __CUDACC__
	template <class T, eDev Dev>
	template <eDev devn = Dev>
	enable_if_edev_gpu<devn, void>
	Opt_Flow <T, Dev>::fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y)
	{
		auto d_grid_blk = grid_2d.d_grid_blk();
		detail_gpu::fcn_opt_flow<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, M_s, M_m, alpha, v_x, v_y);
	}
#endif
}