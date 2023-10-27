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

#pragma once

#include "const_enum.h"
#include "grid_2d.h"
#include "vctr_cpu.h"
#include "vctr_gpu.h"
#include "stream_cpu.h"
#include "stream_gpu.h"
#include "fft_cpu.h"
#include "fft_gpu.h"

namespace mt
{
	template <class T, eDev Dev>
	class Opt_Flow
	{
	public:
		using T_r = T;
		using TVctr = Vctr<T, Dev>;

		static const eDev device = Dev;

		Opt_Flow();

		Opt_Flow(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i);

		inline
		void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i);

		void operator()(TVctr& M_s, TVctr& M_m, T alpha, T sigma, dt_int32 n_iter, TVctr& v_xt, TVctr& v_yt);

		void set_fft_plan();

		void cleanup();

	protected:

		void set_regular_grid(TVctr& Rx, TVctr& Ry);

		/***************************************** cpu *****************************************/
		template <eDev devn = Dev>
		enable_if_edev_cpu<devn, void>
		fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y);

		/**********************Device**********************/
	#ifdef __CUDACC__
		template <eDev devn = Dev>
		enable_if_edev_gpu<devn, void>
		fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y);
	#endif

		Stream<Dev> *stream;
		Grid_2d<T> grid_2d;

		Interp_rn_2d<T, Dev> intrpl_2d;
		Gauss_Cv_2d<T, Dev> gauss_cv_2d;

		TVctr v_x;
		TVctr v_y;

		TVctr Rx;
		TVctr Ry;
		TVctr M;
	};
}

#include "../src/opt_flow.inl"