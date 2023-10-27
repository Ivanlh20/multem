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

template <class T, mt::eDev Dev>
void run_sd_2d(mt::System_Config &system_config, pMLD &rM_i, 
T bs_x, T bs_y, T ds_x, T phi_x, T ds_y, T phi_y, dt_int32 seed, mt::eFil_Sel_Typ bg_opt, 
dt_float64 &bg, pMLD &rM_o, pMLD &rdx_o, pMLD &rdy_o)
{
	mt::Grid_2d<T> grid_2d(rM_i.cols, rM_i.rows, bs_x, bs_y);
	mt::Stream<Dev> stream(system_config.n_stream);

	mt::Vctr<T, Dev> M(rM_i.begin(), rM_i.end());
	mt::Vctr<T, Dev> mx_o(rM_o.size());
	mt::Vctr<T, Dev> dx_o(grid_2d.ny);
	mt::Vctr<T, Dev> dy_o(grid_2d.ny);

	mt::Sd_2d<T, Dev> sd_2d(&stream, grid_2d, bg_opt, bg);
	sd_2d.generate_dx_dy(ds_x, phi_x, ds_y, phi_y, seed, dx_o, dy_o);
	bg = sd_2d(M, dx_o, dy_o, mx_o);

	thrust::copy(mx_o.begin(), mx_o.end(), rM_o.begin());
	thrust::copy(dx_o.begin(), dx_o.end(), rdx_o.begin());
	thrust::copy(dy_o.begin(), dy_o.end(), rdy_o.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;
	auto rM_i = mex_get_pvctr<pMLD>(prhs[idx_0+0]);
	auto bs_x = rM_i.cols*mex_get_num<dt_float64>(prhs[idx_0+1]);
	auto bs_y = rM_i.rows*mex_get_num<dt_float64>(prhs[idx_0+2]);
	auto ds_x = mex_get_num<dt_float64>(prhs[idx_0+3]);
	auto phi_x = mex_get_num<dt_float64>(prhs[idx_0+4]);
 	auto ds_y = mex_get_num<dt_float64>(prhs[idx_0+5]);
	auto phi_y = mex_get_num<dt_float64>(prhs[idx_0+6]);
	auto seed = (nrhs>idx_0+7)?mex_get_num<dt_int32>(prhs[idx_0+7]):0;
	auto bg_opt = (nrhs > idx_0 + 8)?mex_get_num<mt::eFil_Sel_Typ>(prhs[idx_0 + 8]):mt::efst_min;
	auto bg = (nrhs > idx_0 + 9)?mex_get_num<dt_float64>(prhs[idx_0 + 9]):0.0;

	/***************************************************************************************/
	auto rM_o = mex_create_pVctr<pMLD>(rM_i.rows, rM_i.cols, plhs[0]);
	auto rdx_o = mex_create_pVctr<pMLD>(rM_i.rows, 1, plhs[1]);
	auto rdy_o = mex_create_pVctr<pMLD>(rM_i.rows, 1, plhs[2]);

	if (system_config.is_float32_cpu())
	{
		run_sd_2d<dt_float32, mt::edev_cpu>(system_config, rM_i, bs_x, bs_y, ds_x, phi_x, ds_y, phi_y, seed, bg_opt, bg, rM_o, rdx_o, rdy_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_sd_2d<dt_float64, mt::edev_cpu>(system_config, rM_i, bs_x, bs_y, ds_x, phi_x, ds_y, phi_y, seed, bg_opt, bg, rM_o, rdx_o, rdy_o);
	}
	else if (system_config.is_float32_gpu())
	{
		run_sd_2d<dt_float32, mt::edev_gpu>(system_config, rM_i, bs_x, bs_y, ds_x, phi_x, ds_y, phi_y, seed, bg_opt, bg, rM_o, rdx_o, rdy_o);
	}
	else if (system_config.is_float64_gpu())
	{
		run_sd_2d<dt_float64, mt::edev_gpu>(system_config, rM_i, bs_x, bs_y, ds_x, phi_x, ds_y, phi_y, seed, bg_opt, bg, rM_o, rdx_o, rdy_o);
	}

	if (nrhs == 4)
	{
		mex_create_set_num<pMLD>(plhs[3], bg);
	}
}