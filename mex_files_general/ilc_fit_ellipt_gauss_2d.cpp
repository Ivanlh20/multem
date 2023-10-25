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
#include "eval_fit_gaussians.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIxy_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto bs_x = rIxy_i.cols*mex_get_num<dt_float64>(prhs[1]);
	auto bs_y = rIxy_i.rows*mex_get_num<dt_float64>(prhs[2]);
	auto pxy = mex_get_r_2d<dt_float64>(prhs[3]) - mt::R_2d<dt_float64>(1, 1);
	auto sigma = mex_get_num<dt_float64>(prhs[4]);
	auto radius = (nrhs>5)?mex_get_num<dt_float64>(prhs[5]):1.25*sigma;
	radius = (radius<sigma)?1.25*sigma:radius;

	/***************************************************************************************/
	auto rpar_o = mex_create_pVctr<pMLD>(1, 6, plhs[0]);

	vector<dt_float32> Ixy_i(rIxy_i.begin(), rIxy_i.end());
	mt::Grid_2d<dt_float32> grid_2d(rIxy_i.cols, rIxy_i.rows, bs_x, bs_y);

	mt::Fit_Ellipt_Gauss_2d<dt_float32, edev_cpu> Fit_Ellipt_Gauss_2d;
	Fit_Ellipt_Gauss_2d.init_variables(grid_2d, 1.0);
	auto parm = Fit_Ellipt_Gauss_2d.fit(Ixy_i, pxy, sigma, radius);

	rpar_o[0] = parm[0] + 1;
	rpar_o[1] = parm[1] + 1;
	rpar_o[2] = parm[2];
	rpar_o[3] = parm[3];
	rpar_o[4] = parm[4];
	rpar_o[5] = parm[5];
}