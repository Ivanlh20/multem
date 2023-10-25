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

#include "types.cuh"
#include "lapack.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIm_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto bs_x = rIm_i.size()*mex_get_num<dt_float64>(prhs[1]);
	auto x = mex_get_num<dt_float64>(prhs[2])-1;
	auto sigma = mex_get_num<dt_float64>(prhs[3]);
	auto radius = (nrhs>4)?mex_get_num<dt_float64>(prhs[4]):1.25*sigma;
	radius = (radius<sigma)?1.25*sigma:radius;
	/***************************************************************************************/
	auto rpar_o = mex_create_pVctr<pMLD>(3, 1, plhs[0]);

	using mt_type = dt_float32;

	vector<mt_type> Im_i(rIm_i.begin(), rIm_i.end());

	mt::Grid_1d<mt_type> grid_1d(rIm_i.size(), bs_x);

	auto parm = mt::fit_gauss_1d(grid_1d, Im_i, x, sigma, radius);

	rpar_o[0] = parm[0] + 1;
	rpar_o[1] = parm[1];
	rpar_o[2] = parm[2];
}