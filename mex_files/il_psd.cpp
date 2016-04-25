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
#include "image_functions.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::Vector;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto np_min = mx_get_scalar<double>(prhs[1]);
	auto d_delta = mx_get_scalar<double>(prhs[2]);
	double d_delta_0 = -90;
	double d_delta_e = 90;

	/*******************************************************************/
	vector<float> Im(rIm_i.begin(), rIm_i.end());

	mt::Grid<float> grid(rIm_i.cols, rIm_i.rows);
	mt::Stream<e_host> stream(4);

	vector<float> angle;
	vector<float> psd;
	mt::PSD(stream, grid, Im, np_min, d_delta_0, d_delta_e, d_delta, angle, psd);
	auto peaks = mt::PSD_find_peaks(stream, grid, Im, np_min, angle, psd);

	/*******************************************************************/
	auto rangle = mx_create_matrix<rmatrix_r>(angle.size(), 1, plhs[0]);
	auto rpsd = mx_create_matrix<rmatrix_r>(psd.size(), 1, plhs[1]);
	auto rpeaks = mx_create_matrix<rmatrix_r>(peaks.size(), 1, plhs[2]);

	rangle.assign(angle.begin(), angle.end());
	rpsd.assign(psd.begin(), psd.end());
	rpeaks.assign(peaks.begin(), peaks.end());
}