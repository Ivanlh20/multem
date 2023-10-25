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

#include <mex.h>
#include "matlab_mex.h"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto pbs = mex_get_pvctr<dt_float64>(prhs[0]);
	auto pn = mex_get_pvctr<dt_float64>(prhs[1]);

	/***************************************************************************************/
	auto pg_max = mex_create_num<dt_float64>(plhs[0]);

	dt_float64 v = pn[0]/pbs[0];
	for (auto ik = 1; ik < pbs.size(); ik++)
	{
		v = min(v, pn[0]/pbs[0]);
	}

	pg_max[0] = 0.5*v;
}
