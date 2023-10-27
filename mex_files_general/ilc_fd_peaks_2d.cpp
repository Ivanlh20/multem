/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#define MATLAB_BLAS_LAPACK

#include "types.cuh"
#include "type_traits_gen.h"
#include "peak_finding.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIm_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto bs_x = rIm_i.cols*mex_get_num<dt_float64>(prhs[1]);
	auto bs_y = rIm_i.rows*mex_get_num<dt_float64>(prhs[2]);
	auto sigma_r = mex_get_num<dt_float64>(prhs[3]);
	auto thres = mex_get_num<dt_float64>(prhs[4]);
	auto niter = mex_get_num<dt_int32>(prhs[5]);

	/***************************************************************************************/
	mt::Vctr<dt_float32, mt::edev_cpu> Im(rIm_i.begin(), rIm_i.end());
	mt::Vctr<dt_float32, mt::edev_cpu> x, y, A, S;
	dt_float32 bg;
	mt::Grid_2d<dt_float32> grid_2d(rIm_i.cols, rIm_i.rows, bs_x, bs_y);

	mt::Peak_Finding<dt_float32, mt::edev_cpu> peaks(grid_2d, Im, sigma_r, thres, niter);

	peaks.find();
	// peaks.fit();
	peaks.get_data(x, y, A, S, bg);
	peaks.cleanup();

	auto rpeaks = mex_create_pVctr<pMLD>(x.size(), 4, plhs[0]);
	auto rbg = mex_create_pVctr<pMLD>(1, 1, plhs[1]);

	for(auto ix = 0; ix<x.size(); ix++)
	{
		rpeaks[0*x.size()+ix] = x[ix] + 1;
		rpeaks[1*x.size()+ix] = y[ix] + 1;
		rpeaks[2*x.size()+ix] = A[ix];
		rpeaks[3*x.size()+ix] = S[ix];
	}

	rbg.real[0] = bg;
}