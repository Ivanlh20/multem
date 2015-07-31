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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#include "host_functions.hpp"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	int n, nx, ny, shift;
	double dx, dy, Radius;
	rmatrix_r f;
 
	ny = mx_get_scalar<int>(prhs[0]);
	nx = mx_get_scalar<int>(prhs[1]);
	dy = mx_get_scalar<double>(prhs[2]);
	dx = mx_get_scalar<double>(prhs[3]);
	Radius = mx_get_scalar<double>(prhs[4]);
	n = mx_get_scalar<int>(prhs[5]);
	shift = mx_get_scalar<int>(prhs[6]);

	f = mx_create_matrix<rmatrix_r>(ny, nx, plhs[0]);

	multem::filter_Butterworth_2D(ny, nx, dy, dx, Radius, n, 1, shift, f.real);
}