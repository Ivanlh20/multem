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
#include "types.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::Vctr;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIm_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto np_min = mex_get_num<dt_float64>(prhs[1]);
	auto delta = mex_get_num<dt_float64>(prhs[2])*mt::c_pi/180;

	/***************************************************************************************/
	vector<dt_float32> Im(rIm_i.begin(), rIm_i.end());

	mt::Grid_2d<dt_float32> grid_2d(rIm_i.cols, rIm_i.rows);

	auto pIm = mt::projected_intensity(grid_2d, Im, np_min, delta);

	auto rpIm = mex_create_pVctr<pMLD>(pIm.size(), 1, plhs[0]);

	rpIm.assign(pIm.begin(), pIm.end());
}