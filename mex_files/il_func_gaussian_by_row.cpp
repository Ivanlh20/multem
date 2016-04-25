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
#include "stream.cuh"
#include "host_functions.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto nx = mx_get_scalar<int>(prhs[0]);
	auto ny = mx_get_scalar<int>(prhs[1]);
	auto dx = mx_get_scalar<double>(prhs[2]);
	auto dy = mx_get_scalar<double>(prhs[3]);
	auto sigma = mx_get_scalar<double>(prhs[4]);
	auto shift = mx_get_scalar<bool>(prhs[5]);

	/******************************************************************/
	auto rf = mx_create_matrix<rmatrix_r>(ny, nx, plhs[0]);

	using TVector = vector<float>;
	mt::Stream<e_host> stream(4);

	TVector f;
	if(min(nx, ny) == 1)
	{
		auto nr = (nx>ny)?nx:ny;
		auto dr = (nx>ny)?dx:dy;
		f = mt::func_gaussian_1d<TVector>(stream, nr, dr, sigma, shift);
	}
	else
	{
		mt::Grid<float> grid(nx, ny, nx*dx, ny*dy);
		f = mt::func_gaussian_2d_by_row<TVector>(stream, grid, sigma, shift);
	}

	rf.assign(f.begin(), f.end());
}