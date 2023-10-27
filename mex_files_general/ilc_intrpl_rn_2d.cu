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

#include "math_mt.h"
#include "grid_2d.h"
#include "vctr_cpu.h"
#include "vctr_gpu.h"
#include "stream_cpu.h"
#include "stream_gpu.h"
#include "intrpl_rn_2d.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto pmx_i = mex_get_pvctr<T>(prhs[system_config.idx_0+0]);
	auto prx_i = mex_get_pvctr<T>(prhs[system_config.idx_0+1]);
	auto pry_i = mex_get_pvctr<T>(prhs[system_config.idx_0+2]);
	auto bg_opt = (nrhs>system_config.idx_0+3)?mex_get_enum<mt::eFil_Sel_Typ>(prhs[system_config.idx_0+3]):mt::efst_min;
	auto bg = (nrhs>system_config.idx_0+4)? mex_get_num<dt_float64>(prhs[system_config.idx_0+4]):0;

	mt::Grid_2d<T> grid_2d(pmx_i.s1(), pmx_i.s0());
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::Vctr<T, Dev> mx_i(pmx_i);
	mt::Vctr<T, Dev> rx_i(prx_i);
	mt::Vctr<T, Dev> ry_i(pry_i);
	mt::Vctr<T, Dev> mx_o(rx_i.shape());

	auto pmx_o = mex_create_pVctr<T>(prx_i.shape(), plhs[0]);

	mt::Interp_rn_2d<T, Dev> intrpl_2d(&stream, grid_2d, bg_opt, bg);
	bg = intrpl_2d(mx_i, rx_i, ry_i, mx_o);

	mx_o.cpy_to_cpu_ptr(pmx_o.begin(), pmx_o.end());

	if (nlhs == 2)
	{
		mex_create_set_num<T>(plhs[1], bg);
	}
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}