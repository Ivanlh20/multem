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
#include "image_functions.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::rmatrix_c;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	double lx = rIm_i.cols*mx_get_scalar<double>(prhs[1]);
	double ly = rIm_i.rows*mx_get_scalar<double>(prhs[2]);
	auto sigma = mx_get_scalar<double>(prhs[3]);
	auto thres = mx_get_scalar<double>(prhs[4]);
	auto niter = mx_get_scalar<int>(prhs[5]);

	/*******************************************************************/
	vector<float> Im(rIm_i.begin(), rIm_i.end());
	vector<float> x, y, A, S;

	multem::Grid<float> grid(rIm_i.cols, rIm_i.rows, lx, ly);
	multem::Stream<e_host> stream(4);
	multem::FFT2<float, e_host> fft2;

	multem::find_peaks_2d(stream, fft2, grid, Im, sigma, thres, niter, x, y, A, S);

	fft2.cleanup();

	auto rpeaks = mx_create_matrix<rmatrix_r>(x.size(), 4, plhs[0]);

	std::copy(x.begin(), x.end(), rpeaks.real+0*x.size());
	std::copy(y.begin(), y.end(), rpeaks.real+1*x.size());
	std::copy(A.begin(), A.end(), rpeaks.real+2*x.size());
	std::copy(S.begin(), S.end(), rpeaks.real+3*x.size());
}