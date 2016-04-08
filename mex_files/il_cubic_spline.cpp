/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include <vector>
#include "types.cuh"
#include "matlab_types.cuh"
#include "cubic_spline.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::rmatrix_c;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto xr = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto yr = mx_get_matrix<rmatrix_r>(prhs[1]);
	auto xi = mx_get_matrix<rmatrix_r>(prhs[2]);

	std::vector<double> x(xr.begin(), xr.end());
	std::vector<double> y(yr.begin(), yr.end());

	/*******************************************************************/
	auto yi = mx_create_matrix<rmatrix_r>(xi.rows, xi.cols, plhs[0]);

	multem::Cubic_Spline<double> spline;
	spline.set_points(x, y);
	spline.eval_function(xi, yi);
}