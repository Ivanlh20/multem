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

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "stream.cuh"
#include "image_functions.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::rmatrix_c;
using multem::Vector;
using multem::e_host;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto dy = mx_get_scalar<double>(prhs[1]);
	auto dx = mx_get_scalar<double>(prhs[2]);
	auto x = mx_get_scalar<double>(prhs[3]);
	auto y = mx_get_scalar<double>(prhs[4]);
	double lx = rIm_i.cols*dx;
	double ly = rIm_i.rows*dy;
	bool bwl = false;
	bool pbc_xy = true;
	double dz = 0.5;

	multem::Grid<double> grid;
	grid.set_input_data(rIm_i.cols, rIm_i.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Stream<e_host> stream(4);
	multem::FFT2<double, e_host> fft2;
	fft2.create_plan(grid.ny, grid.nx, 4);

	/*****************************************************************************/
	auto rIm_o = mx_create_matrix<rmatrix_r>(rIm_i.rows, rIm_i.cols, plhs[0]);

	Vector<complex<double>, e_host> M_1(grid.nxy());
	multem::assign(stream, rIm_i, M_1);
	multem::shift_2D(stream, fft2, grid, x, y, M_1);
	multem::copy_to_host(stream, M_1, rIm_o);
	fft2.cleanup();
}