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

#include "const_enum.h"
#include "vctr_cpu.h"
#include "fcns_image_cpu.h"
#include "fcns_cpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto pr = mex_get_pvctr<T>(prhs[0]);
	auto pfr = mex_get_pvctr<T>(prhs[1]);
	auto r_0 = mex_get_num<T>(prhs[2]);
	auto r_e = mex_get_num<T>(prhs[3]);
	auto n_r = mex_get_num<T>(prhs[4]);
	auto typ = (nrhs>5)?mex_get_num<dt_int32>(prhs[5]):1;

	/***************************************************************************************/
	auto prl = mex_create_pVctr<T>({n_r}, plhs[0]);
	auto pfrl = mex_create_pVctr<T>({n_r}, plhs[1]);
	auto pcfrl = mex_create_pVctr<T>({n_r}, plhs[2]);

	const bool reg = true;

	if (reg)
	{
		mt::fcn_rad_dist_by_div(pr, pfr, r_0, r_e, n_r, prl, pfrl, pcfrl, typ);
	}
	else
	{
		mt::fcn_rad_dist_by_srch(pr, pfr, r_0, r_e, n_r, prl, pfrl, pcfrl, typ);
	}
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT_OUT(mex_run, 0);
}