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
#include "wd_fermi.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto shape = mex_get_shape(prhs[system_config.idx_0+0]);
	auto bs_i = mex_get_r_3d_bs<T>(prhs[system_config.idx_0+1], shape);
	auto alpha = mex_get_num<T>(prhs[system_config.idx_0+2]);
	auto r_wd = (nrhs>system_config.idx_0+3)?mex_get_num<T>(prhs[system_config.idx_0+3]):alpha;
	auto sft = (nrhs>system_config.idx_0+4)?mex_get_bool(prhs[system_config.idx_0+4]):false;
	auto r_c_i = (nrhs>system_config.idx_0+6)?mex_get_r_3d_r_c<T>(prhs[system_config.idx_0+6], bs_i/T(2)):bs_i/T(2);

	/***************************************************************************************/
	auto pmx_o = mex_create_pVctr<T>(shape, plhs[0]);
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::Vctr<T, Dev> mx_o(shape);

	if (pmx_o.is_1d())
	{
		const auto bs = mt::fcn_max(bs_i.x, bs_i.y);
		const auto r_max = (nrhs>system_config.idx_0+5)?mex_get_num<T>(prhs[system_config.idx_0+5]):bs/T(2);
		const auto r_c = (bs_i.x>bs_i.y)?r_c_i.x:r_c_i.y;
		mt::Wd_Fermi_1d<T> wd_fcn(r_c, alpha, r_wd, r_max);
		mt::Grid_1d<T> grid(bs, shape.shape_size());

		mt::fcn_wd_fermi_1d(grid, wd_fcn, sft, mx_o, &stream);
	}
	else
	{
		const auto bs = mt::R_2d<T>(bs_i.x, bs_i.y);
		const auto r_max = (nrhs>system_config.idx_0+5)?mex_get_num<T>(prhs[system_config.idx_0+5]):bs.min()/T(2);
		const auto r_c = mt::R_2d<T>(r_c_i.x, r_c_i.y);
		mt::Wd_Fermi_2d<T> wd_fcn(r_c, alpha, r_wd, r_max);
		mt::Grid_2d<T> grid(bs.x, bs.y, shape[1], shape[0]);

		mt::fcn_wd_fermi_2d(grid, wd_fcn, sft, mx_o, &stream);
	}

	mx_o.cpy_to_cpu_ptr(pmx_o.begin(), pmx_o.end());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}