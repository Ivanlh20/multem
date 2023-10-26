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
#include "vctr_cpu.h"
#include "vctr_gpu.h"
#include "stream_cpu.h"
#include "stream_gpu.h"
#include "fcns_cpu.h"
#include "fcns_gpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto mx_i_cpu = mex_get_vctr<T>(prhs[system_config.idx_0+0]);

	/***************************************************************************************/
	auto pdmx_x = mex_create_pVctr<T>(mx_i_cpu.shape(), plhs[0]);
	auto pdmx_y = mex_create_pVctr<T>(mx_i_cpu.shape(), plhs[1]);

	mt::Vctr<T, Dev> mx_i = mx_i_cpu;
	mt::Vctr<T, Dev> dmx_x(mx_i_cpu.shape());
	mt::Vctr<T, Dev> dmx_y(mx_i_cpu.shape());

	mt::Stream<Dev> stream(system_config.n_stream);

	mt::fcn_grad(mx_i, dmx_x, dmx_y, &stream);

	dmx_x.cpy_to_cpu_ptr(pdmx_x.begin(), pdmx_x.end());
	dmx_y.cpy_to_cpu_ptr(pdmx_y.begin(), pdmx_y.end());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}