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

#include <algorithm>
#include "types.cuh"
#include "quad_data.cuh"

#include <mex.h>
#include "matlab_mex.h"
	
template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto q_type = mex_get_enum<mt::eQuad_Typ>(prhs[0]);
	auto quad_pts = mex_get_num<dt_int32>(prhs[1]);
	auto alpha = (nrhs>2)?mex_get_num<T>(prhs[2]):0;
	auto beta = (nrhs>3)?mex_get_num<T>(prhs[3]):0;
	auto a = (nrhs>4)?mex_get_num<T>(prhs[4]):0;
	auto b = (nrhs>5)?mex_get_num<T>(prhs[5]):1;

	/***************************************************************************************/
	mt::Quad_Coef_1d_cpu<T> quad_coef;
	mt::Quad_Data quad_data(q_type, quad_pts, quad_coef, alpha, beta, a, b);

	mex_create_set_pVctr<T>(plhs[0], quad_coef.x.ptr_64());
	mex_create_set_pVctr<T>(plhs[1], quad_coef.w.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	mex_run<dt_float64>(nlhs, plhs, nrhs, prhs);
}