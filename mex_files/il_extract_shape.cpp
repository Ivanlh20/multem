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
	// read input data from matlab
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto dR = mx_get_scalar<double>(prhs[1]);
	double lx = rIm_i.cols*dR;
	double ly = rIm_i.rows*dR;
	bool bwl = false;
	bool pbc_xy = true;
	double dz = 0.5;


	multem::Grid<float> grid;
	grid.set_input_data(rIm_i.cols, rIm_i.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Stream<e_host> stream(4);

	// create fourier transform plan
	multem::FFT2<float, e_host> fft2;
	fft2.create_plan_2d(grid.ny, grid.nx, 4);

	// create vectors
	multem::Vector<float, e_host> Im_i(rIm_i.begin(), rIm_i.end());

	// extract shape
	auto Im_o = multem::extract_shape(stream, fft2, grid, Im_i);

	// clean fft2 plan
	fft2.cleanup();

	// create output data to matlab
	auto rIm_o = mx_create_matrix<rmatrix_r>(rIm_i.rows, rIm_i.cols, plhs[0]);

	// copy to matlab output
	multem::copy_to_host(stream, Im_o, rIm_o);
}