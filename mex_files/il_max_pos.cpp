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
#include "lapack.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto dy = mx_get_scalar<double>(prhs[1]);
	auto dx = mx_get_scalar<double>(prhs[2]);
	auto x = mx_get_scalar<double>(prhs[3]);
	auto y = mx_get_scalar<double>(prhs[4]);
	auto radius = mx_get_scalar<double>(prhs[5]);

	double lx = rIm.cols*dx;
	double ly = rIm.rows*dy;
	bool bwl = false;
	bool pbc_xy = false;
	double dz = 0.5;

	multem::Grid<double> grid;
	grid.set_input_data(rIm.cols, rIm.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Stream<e_host> stream(4);
	multem::Vector<double, e_host> Im(rIm.size());
	multem::assign(stream, rIm, Im);

	/*****************************************************************************/
	auto rx_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[0]);
	auto ry_o = mx_create_matrix<rmatrix_r>(1, 1, plhs[1]);

	multem::max_pos(stream, grid, Im, x, y, radius, rx_o.real[0], ry_o.real[0]);
}