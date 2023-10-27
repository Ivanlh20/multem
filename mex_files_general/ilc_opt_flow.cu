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

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto pmx_s = mex_get_pvctr<pMLD>(prhs[system_config.idx_0+0]);
	auto pmx_m = mex_get_pvctr<pMLD>(prhs[system_config.idx_0+1]);
	auto alpha = (nrhs>system_config.idx_0+2)?mex_get_num<T>(prhs[system_config.idx_0+2]):T(0.5);
	auto sigma = (nrhs>system_config.idx_0+3)?mex_get_num<T>(prhs[system_config.idx_0+3]):T(2.0);
	auto n_iter = (nrhs>system_config.idx_0+4)?mex_get_num<dt_int32>(prhs[system_config.idx_0+4]):1;
	auto prx_i = (nrhs>system_config.idx_0+5)?mex_get_pvctr<T>(prhs[system_config.idx_0+5]):pMLD();
	auto pry_i = (nrhs>system_config.idx_0+6)?mex_get_pvctr<T>(prhs[system_config.idx_0+6]):pMLD();

	/***************************************************************************************/
	mt::Grid_2d<T> grid_2d(pmx_s.s1(), pmx_s.s0());
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;

	mt::Vctr<T, Dev> mx_s(pmx_s);
	mt::Vctr<T, Dev> mx_m(pmx_m);
	mt::Vctr<T, Dev> rx(prx_i);
	mt::Vctr<T, Dev> ry(pry_i);

	mt::Opt_Flow<T, Dev> fcn_opt_flow(&stream, &fft_2d, grid_2d);

	if (!mt::fcn_is_zero(sigma))
	{
		fcn_opt_flow.set_fft_plan();
	}

	fcn_opt_flow(mx_s, mx_m, alpha, sigma, n_iter, rx, ry);
	fcn_opt_flow.cleanup();

	mex_create_set_pVctr<T>(plhs[0], rx.ptr_64());
	mex_create_set_pVctr<T>(plhs[1], ry.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}