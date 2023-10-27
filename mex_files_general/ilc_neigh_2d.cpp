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

#include "types.cuh"
#include "cgpu_stream.cuh"
#include "box_occ.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;
using mt::edev_cpu;
using mt::Value_type;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rx = mex_get_pvctr<pMLD>(prhs[0]);
	auto ry = mex_get_pvctr<pMLD>(prhs[1]);
	auto r_neigh = mex_get_num<dt_float64>(prhs[2]);

	mt::Vctr<dt_float64, edev_cpu> x(rx.begin(), rx.end());
	mt::Vctr<dt_float64, edev_cpu> y(ry.begin(), ry.end());

	/***************************************************************************************/
	mt::Stream<edev_cpu> stream(4);
	mt::Neigh_2d<dt_float64> list(stream, x, y, r_neigh);

	const char *field_names[] = {"index"};
	dt_int32 number_of_fields = 1;
	mwSize dims_output[2] = {list.size(), 1};

	plhs[0] = mxCreateStructArray(2, dims_output, number_of_fields, field_names);

	for(auto idx=0; idx<x.size(); idx++)
	{
		dt_int32 nx = list[idx].size();
		dt_int32 ny = (nx == 0)?0:1;
		mex_create_set_pVctr_field<pMLD>(plhs[0], idx, "index", ny, nx, list[idx]);
	}
}