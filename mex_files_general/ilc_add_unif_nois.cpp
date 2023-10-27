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
#include "vctr_cpu.h"
#include "rnd_1d_cpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto mx_i = mex_get_vctr<T>(prhs[0]);
	auto sc = (nrhs>1)?mex_get_num<T>(prhs[1]):T(1);
	auto seed = (nrhs>2)?mex_get_num<dt_int32>(prhs[2]):0;

	/***************************************************************************************/
	mt::Vctr_cpu<T> mx_o(mx_i.shape());

	mt::fcn_add_unif_nois(mx_i, sc, seed, mx_o);

	mex_create_set_pVctr<T>(plhs[0], mx_o.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT_OUT(mex_run, 0);
}