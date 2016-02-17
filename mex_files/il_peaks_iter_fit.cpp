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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto dy_i = mx_get_scalar<double>(prhs[1]);
	auto dx_i = mx_get_scalar<double>(prhs[2]);
	auto peaks_i = mx_get_matrix<rmatrix_r>(prhs[3]);
	auto radius_i = mx_get_scalar<double>(prhs[4]);
	auto niter_i = mx_get_scalar<double>(prhs[5]);

	double lx = rIm_i.cols*dx_i;
	double ly = rIm_i.rows*dy_i;
	bool bwl = false;
	bool pbc_xy = false;
	double dz = 0.5;

	multem::Grid<double> grid;
	grid.set_input_data(rIm_i.cols, rIm_i.rows, lx, ly, dz, bwl, pbc_xy);

	multem::Stream<e_host> stream(4);
	multem::Vector<double, e_host> Im(rIm_i.size());
	multem::assign(stream, rIm_i, Im);

	/*****************************************************************************/
	auto peaks_o = mx_create_matrix<rmatrix_r>(peaks_i.rows, peaks_i.cols, plhs[0]);
	auto npeaks = peaks_i.rows;

	peaks_o.assign(peaks_i.begin(), peaks_i.end());
	for(auto iter= 0; iter<niter_i; iter++)
	{
		for(auto i= 0; i<npeaks; i++)
		{
			double &x = peaks_o[0*npeaks+i];
			double &y = peaks_o[1*npeaks+i];
			multem::Rx_Ry_fit(grid, Im, x, y, radius_i, x, y);
		}
	}
}