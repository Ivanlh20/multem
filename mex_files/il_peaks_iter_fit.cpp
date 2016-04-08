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
using multem::r2d;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	auto rIm_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	double lx_i = rIm_i.cols*mx_get_scalar<double>(prhs[1]);
	double ly_i = rIm_i.rows*mx_get_scalar<double>(prhs[2]);
	auto peaks_i = mx_get_matrix<rmatrix_r>(prhs[3]);
	auto radius_i = mx_get_scalar<double>(prhs[4]);
	auto niter_i = mx_get_scalar<double>(prhs[5]);

	/*******************************************************************/
	auto peaks_o = mx_create_matrix<rmatrix_r>(peaks_i.rows, peaks_i.cols, plhs[0]);

	vector<float> Im(rIm_i.begin(), rIm_i.end());

	multem::Grid<float> grid(rIm_i.cols, rIm_i.rows, lx_i, ly_i);

	auto npeaks = peaks_i.rows;

	peaks_o.assign(peaks_i.begin(), peaks_i.end());
	for(auto iter= 0; iter<niter_i; iter++)
	{
		for(auto i = 0; i<npeaks; i++)
		{
			double &x = peaks_o[0*npeaks+i];
			double &y = peaks_o[1*npeaks+i];
			auto fit = multem::Rx_Ry_fit(grid, Im, r2d<float>(x, y), radius_i, true);

			x = fit[0];
			y = fit[1];
		}
	}
}