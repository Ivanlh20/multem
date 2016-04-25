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

using mt::rmatrix_r;
using mt::Vector;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rx_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto ry_i = mx_get_matrix<rmatrix_r>(prhs[1]);
	auto y_thr = mx_get_scalar<double>(prhs[2]);

	std::vector<double> x(rx_i.begin(), rx_i.end());
	std::vector<double> y(ry_i.begin(), ry_i.end());
	auto peaks = mt::find_peaks_vector_typ_1(x, y, y_thr);

	/*******************************************************************/
	auto rpeaks = mx_create_matrix<rmatrix_r>(peaks.size(), 1, plhs[0]);
	std::copy(peaks.begin(), peaks.end(), rpeaks.real);
}