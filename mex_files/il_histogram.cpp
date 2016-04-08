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
#include "image_functions.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::rmatrix_c;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto nbins = mx_get_scalar<int>(prhs[1]);

	/******************************************************************/
	auto rx = mx_create_matrix<rmatrix_r>(nbins, 1, plhs[0]);
	auto ry = mx_create_matrix<rmatrix_r>(nbins, 1, plhs[1]);

	vector<float> Im(rIm_i.begin(), rIm_i.end());
	vector<float> x(nbins);
	vector<float> y(nbins);
	multem::histogram(Im, nbins, y, &x);

	rx.assign(x.begin(), x.end());
	ry.assign(y.begin(), y.end());
}