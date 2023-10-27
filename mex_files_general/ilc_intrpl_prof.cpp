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
#include "fcns_cpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto mx_i = mex_get_vctr<T>(prhs[0]);
	auto bs = mex_get_r_2d_bs<T>(prhs[1], mx_i.shape());
	auto p_1 = mex_get_r_2d<T>(prhs[2]);
	auto p_2 = mex_get_r_2d<T>(prhs[3]);
	auto n_p = mex_get_num<dt_int32>(prhs[4]);

	/*******************************************************************/
	mt::Grid_2d<T> grid_2d(bs.x, bs.y, mx_i.s1_32(), mx_i.s0_32());

	auto profile = mt::intrpl_prof(grid_2d, mx_i, p_1, p_2, n_p);

	mex_create_set_pVctr<T>(plhs[0], profile.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT(mex_run, 0);
}