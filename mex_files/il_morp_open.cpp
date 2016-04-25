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

using mt::rmatrix_r;
using mt::rmatrix_c;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto nkr = mx_get_scalar<int>(prhs[1]);

	/******************************************************************/
	auto rIm_o = mx_create_matrix<rmatrix_r>(rIm_i.rows, rIm_i.cols, plhs[0]);

	vector<double> Im(rIm_i.begin(), rIm_i.end());

	mt::Stream<e_host> stream(4);
	Im = mt::morp_g_open(stream, rIm_i.rows, rIm_i.cols, Im, mt::eSE_Square, nkr);

	rIm_o.assign(Im.begin(), Im.end());
}