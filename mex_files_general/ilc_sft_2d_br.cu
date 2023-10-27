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

#define MATLAB_BLAS_LAPACK

#include "math_mt.h"
#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"

#include "cgpu_classes.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;

template <class T, mt::eDev Dev>
void run_sft_2d_br(mt::System_Config &system_config, pMLD &rIm_i, 
T bs_x, T bs_y, pMLD &ry, pMLD &rIm_o)
{
	mt::Grid_2d<T> grid_2d(rIm_i.rows, rIm_i.cols, bs_y, bs_x);
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;

	mt::Vctr<T, Dev> Im(rIm_i.begin(), rIm_i.end());
	mt::Vctr<T, Dev> y(ry.begin(), ry.end());

	mt::fcn_trs_2d(stream, grid_2d.nx, grid_2d.ny, Im);

	mt::Sft_2d_BC<T, Dev> sft_2d_bc(&stream, &fft_2d, grid_2d);
	sft_2d_bc.set_fft_plan();
	sft_2d_bc(1.0, y, Im);
	sft_2d_bc.cleanup();

	mt::fcn_trs_2d(stream, grid_2d.ny, grid_2d.nx, Im);

	thrust::copy(Im.begin(), Im.end(), rIm_o.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	auto rIm_i = mex_get_pvctr<pMLD>(prhs[idx_0+0]);
	auto bs_x = rIm_i.cols*mex_get_num<dt_float64>(prhs[idx_0+1]);
	auto bs_y = rIm_i.rows*mex_get_num<dt_float64>(prhs[idx_0+2]);
	auto rx = mex_get_pvctr<pMLD>(prhs[idx_0+3]);

	/***************************************************************************************/
	auto rIm_o = mex_create_pVctr<pMLD>(rIm_i.rows, rIm_i.cols, plhs[0]);

	if (system_config.is_float32_cpu())
	{
		run_sft_2d_br<dt_float32, mt::edev_cpu>(system_config, rIm_i, bs_x, bs_y, rx, rIm_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_sft_2d_br<dt_float64, mt::edev_cpu>(system_config, rIm_i, bs_x, bs_y, rx, rIm_o);
	}
	else if (system_config.is_float32_gpu())
	{
		run_sft_2d_br<dt_float32, mt::edev_gpu>(system_config, rIm_i, bs_x, bs_y, rx, rIm_o);
	}
	else if (system_config.is_float64_gpu())
	{
		run_sft_2d_br<dt_float64, mt::edev_gpu>(system_config, rIm_i, bs_x, bs_y, rx, rIm_o);
	}
}