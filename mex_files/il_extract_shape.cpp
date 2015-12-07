/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
	// read input data from matlab
	auto mIm = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto mdR = mx_get_scalar<double>(prhs[1]);

	// create output data to matlab
	auto mIm_s = mx_create_matrix<rmatrix_r>(mIm.rows, mIm.cols, plhs[0]);

	// create grid
	multem::Grid<float> grid;
	grid.set_input_data(mIm_s.cols, mIm_s.rows, mdR*mIm_s.cols, mdR*mIm_s.rows, 0.5, false, true);

	// crate stream
	multem::Stream<e_host> stream;
	stream.resize(4);

	// create fourier transform plan
	multem::FFT2<float, e_host> fft2;
	fft2.create_plan(grid.ny, grid.nx, 4);

	// create vectors
	multem::Vector<float, e_host> Im(mIm.begin(), mIm.end());
	multem::Vector<float, e_host> Im_s(mIm.size());

	// extract shape
	multem::extract_shape(stream, fft2, grid, Im, Im_s);

	// clean fft2 plan
	fft2.cleanup();

	// copy to matlab output
	multem::copy_to_host(stream, Im_s, mIm_s);
}