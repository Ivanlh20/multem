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

#include "types.cuh"
#include "matlab_types.cuh"
#include "lapack.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	double lx = rIm_i.cols*mx_get_scalar<double>(prhs[1]);
	double ly = rIm_i.rows*mx_get_scalar<double>(prhs[2]);
	auto x = mx_get_scalar<double>(prhs[3]);
	auto y = mx_get_scalar<double>(prhs[4]);
	auto radius = 0.5*min(lx, ly);

	/*******************************************************************/
	auto rx_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[0]);
	auto ry_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[1]);

	vector<float> Im_i(rIm_i.begin(), rIm_i.end());
	mt::r2d<float> R(x, y);

	mt::Grid<float> grid(rIm_i.cols, rIm_i.rows, lx, ly);
	mt::Stream<e_host> stream(4);

	R = mt::max_pos(stream, grid, Im_i, R, radius);

	rx_o[0] = R.x + 1;
	ry_o[0] = R.y + 1;
}