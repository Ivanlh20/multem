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
#include "matlab_mex.cuh"

using multem::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	rmatrix_r R, fR, Rl;
	rmatrix_r rl, frl, cfrl;

	R = mx_get_matrix<rmatrix_r>(prhs[0]);
	fR = mx_get_matrix<rmatrix_r>(prhs[1]);
	Rl = mx_get_matrix<rmatrix_r>(prhs[2]);

	int typ = mx_get_scalar<int>(prhs[3]);

	/**************************************************************/
	rl = mx_create_matrix<rmatrix_r>(Rl.size-1, 1, plhs[0]);
	frl = mx_create_matrix<rmatrix_r>(Rl.size-1, 1, plhs[1]);
	cfrl = mx_create_matrix<rmatrix_r>(Rl.size-1, 1, plhs[2]);

	multem::radial_distribution_2D(R.size, R.real, fR.real, Rl.size, Rl.real, rl.real, frl.real, cfrl.real, true, typ);
}