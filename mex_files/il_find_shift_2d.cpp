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
#include "traits.cuh"
#include "stream.cuh"
#include "fft2.cuh"
#include "host_device_functions.cuh"
#include "host_functions.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::rmatrix_c;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	auto rM_1i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto rM_2i = mx_get_matrix<rmatrix_r>(prhs[1]);
	double lx = rM_1i.cols*mx_get_scalar<double>(prhs[2]);
	double ly = rM_1i.rows*mx_get_scalar<double>(prhs[3]);
	auto k = mx_get_scalar<double>(prhs[4]);
	auto sigma = mx_get_scalar<double>(prhs[5]);

	/*******************************************************************/
	auto rdx_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[0]);
	auto rdy_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[1]);

	vector<float> M1_i(rM_1i.begin(), rM_1i.end());
	vector<float> M2_i(rM_2i.begin(), rM_2i.end());

	mt::Grid<float> grid(rM_1i.cols, rM_1i.rows, lx, ly);
	mt::Stream<e_host> stream(4);
	mt::FFT2<float, e_host> fft2;

	auto R = mt::find_shift_2d(stream, fft2, grid, M1_i, M2_i, k, sigma);

	fft2.cleanup();

	rdx_o[0] = R.x;
	rdy_o[0] = R.y;
}