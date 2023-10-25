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
void run_opt_flow(mt::System_Config &system_config, pMLD &rM_s, 
pMLD &rM_m, T alpha, T sigma, dt_int32 n_iter, pMLD &rv_x, pMLD &rv_y)
{
	mt::Grid_2d<T> grid_2d(rM_s.cols, rM_s.rows);
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;

	mt::Vctr<T, Dev> M_s(rM_s.begin(), rM_s.end());
	mt::Vctr<T, Dev> M_m(rM_m.begin(), rM_m.end());

	mt::Vctr<T, Dev> v_x(rv_x.begin(), rv_x.end());
	mt::Vctr<T, Dev> v_y(rv_y.begin(), rv_y.end());

	mt::Opt_Flow<T, Dev> fcn_opt_flow(&stream, &fft_2d, grid_2d);

	if (!mt::fcn_is_zero(sigma))
	{
		fcn_opt_flow.set_fft_plan();
	}

	fcn_opt_flow(M_s, M_m, alpha, sigma, n_iter, v_x, v_y);
	fcn_opt_flow.cleanup();

	thrust::copy(v_x.begin(), v_x.end(), rv_x.begin());
	thrust::copy(v_y.begin(), v_y.end(), rv_y.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	auto rM_s = mex_get_pvctr<pMLD>(prhs[idx_0+0]);
	auto rM_m = mex_get_pvctr<pMLD>(prhs[idx_0+1]);
	auto alpha = (nrhs>idx_0+2)?mex_get_num<dt_float64>(prhs[idx_0+2]):1.0;
	auto sigma = (nrhs>idx_0+3)?mex_get_num<dt_float64>(prhs[idx_0+3]):0.0;
	auto n_iter = (nrhs>idx_0+4)?mex_get_num<dt_int32>(prhs[idx_0+4]):1;
	auto rv_xi = (nrhs>idx_0+5)?mex_get_pvctr<pMLD>(prhs[idx_0+5]):pMLD();
	auto rv_yi = (nrhs>idx_0+6)?mex_get_pvctr<pMLD>(prhs[idx_0+6]):pMLD();

	/***************************************************************************************/
	auto rv_x = mex_create_pVctr<pMLD>(rM_s.rows, rM_s.cols, plhs[0]);
	rv_x.assign(rv_xi.begin(), rv_xi.end());

	auto rv_y = mex_create_pVctr<pMLD>(rM_s.rows, rM_s.cols, plhs[1]);
	rv_y.assign(rv_yi.begin(), rv_yi.end());

	if (system_config.is_float32_cpu())
	{
		run_opt_flow<dt_float32, mt::edev_cpu>(system_config, rM_s, rM_m, alpha, sigma, n_iter, rv_x, rv_y);
	}
	else if (system_config.is_float64_cpu())
	{
		run_opt_flow<dt_float64, mt::edev_cpu>(system_config, rM_s, rM_m, alpha, sigma, n_iter, rv_x, rv_y);
	}
	else if (system_config.is_float32_gpu())
	{
		run_opt_flow<dt_float32, mt::edev_gpu>(system_config, rM_s, rM_m, alpha, sigma, n_iter, rv_x, rv_y);
	}
	else if (system_config.is_float64_gpu())
	{
		run_opt_flow<dt_float64, mt::edev_gpu>(system_config, rM_s, rM_m, alpha, sigma, n_iter, rv_x, rv_y);
	}
}