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

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
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
	double lx = rM_i.cols*mx_get_scalar<double>(prhs[1]);
	double ly = rM_i.rows*mx_get_scalar<double>(prhs[2]);
	auto sigma = mx_get_scalar<double>(prhs[3]);
	bool bwl = false;
	bool pbc_xy = true;
	double dz = 0.5;
	multem::Grid<double> grid;
	grid.set_input_data(rM_i.cols, rM_i.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Vector<complex<double>, e_host> M(grid.nxy());
	multem::assign(rM_i, M);
	multem::fft2_shift(grid, M);

	/*****************************************************************************/
	auto rM_o = mx_create_matrix<rmatrix_r>(grid.ny, grid.nx, plhs[0]);

	multem::FFT2<double, e_host> fft2;
	fft2.create_plan(grid.ny, grid.nx, 4);

	multem::gaussian_convolution(grid, fft2, sigma, M);
	for(auto ixy=0; ixy < grid.nxy(); ixy++)
	{
		rM_o[ixy] = std::max(0.0, M[ixy].real());
	}
	multem::fft2_shift(grid, rM_o);

	fft2.cleanup();
}