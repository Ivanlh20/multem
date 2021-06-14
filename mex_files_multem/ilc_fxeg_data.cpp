/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "fxeg_data.hpp"

#include "mex.h"
#include "matlab_mex.cuh"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float64;

	auto Z = mex_get_num<dt_int32>(prhs[0]);
	auto ndg = mex_get_num<dt_int32>(prhs[1]);

	/***************************************************************************************/
	mt::Vctr_cpu<T> g, fxg, feg;
	mt::fxeg_Data<T>(1, Z, ndg, g, fxg, feg);
	
	mex_create_set_pVctr<T>(plhs[0], g.ptr_64());
	mex_create_set_pVctr<T>(plhs[1], fxg.ptr_64());
	mex_create_set_pVctr<T>(plhs[2], feg.ptr_64());
}