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

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::Vector;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto np_min = mx_get_scalar<double>(prhs[1]);
	auto delta = mx_get_scalar<double>(prhs[2])*multem::c_Pi/180;

	/*******************************************************************/
	vector<float> Im(rIm_i.begin(), rIm_i.end());

	multem::Grid<float> grid(rIm_i.cols, rIm_i.rows);

	auto pIm = multem::projected_intensity(grid, Im, np_min, delta);

	auto rpIm = mx_create_matrix<rmatrix_r>(pIm.size(), 1, plhs[0]);

	rpIm.assign(pIm.begin(), pIm.end());
}