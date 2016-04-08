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

using multem::rmatrix_r;
using multem::rmatrix_c;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rM_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto border_x = mx_get_scalar<double>(prhs[1]);
	auto border_y = border_x;
	double lx = rM_i.cols;
	double ly = rM_i.rows;
	bool bwl = false;
	bool pbc_xy = true;
	double dz = 0.5;

	multem::Grid<double> grid_i, grid_o;
	grid_i.set_input_data(rM_i.cols, rM_i.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Stream<e_host> stream(4);

	multem::Vector<double, e_host> M_i(grid_i.nxy()), M_o;
	multem::assign(stream, rM_i, M_i);

	/*******************************************************************/
	multem::add_PB_border(grid_i, M_i, border_x, border_y, grid_o, M_o);

	auto rM_o = mx_create_matrix<rmatrix_r>(grid_o.ny, grid_o.nx, plhs[0]);

	multem::copy_to_host(stream, M_o, rM_o);

}