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

#include "math_mt.h"
#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"
#include "eval_fit_gaussians.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIxy_i = mex_get_pvctr<pMLD>(prhs[0]);
	dt_float64 lx_i = rIxy_i.cols*mex_get_num<dt_float64>(prhs[1]);
	dt_float64 ly_i = rIxy_i.rows*mex_get_num<dt_float64>(prhs[2]);
	auto peaks_i = mex_get_pvctr<pMLD>(prhs[3]);
	auto radius_i = mex_get_num<dt_float64>(prhs[4]);
	auto niter_i = mex_get_num<dt_float64>(prhs[5]);
	auto sigma_i = 0.5*radius_i;

	/***************************************************************************************/
	auto peaks_o = mex_create_pVctr<pMLD>(peaks_i.rows, peaks_i.cols, plhs[0]);

	vector<dt_float64> Im(rIxy_i.begin(), rIxy_i.end());

	mt::Grid_2d<dt_float64> grid_2d(rIxy_i.cols, rIxy_i.rows, lx_i, ly_i);
	mt::Fit_Ellipt_Gauss_2d<dt_float64, edev_cpu> Fit_Ellipt_Gauss_2d;
	Fit_Ellipt_Gauss_2d.init_variables(grid_2d, 1.0);

	auto npeaks = peaks_i.rows;

	peaks_o.assign(peaks_i.begin(), peaks_i.end());
	for(auto iter = 0; iter<niter_i; iter++)
	{
		for(auto i = 0; i<npeaks; i++)
		{
			auto pxy = mt::R_2d<dt_float64>(peaks_o[0*npeaks+i], peaks_o[1*npeaks+i]);
			pxy = Fit_Ellipt_Gauss_2d.fit(Im, pxy, sigma_i, radius_i);

			peaks_o[0*npeaks+i] = pxy.x;
			peaks_o[1*npeaks+i] = pxy.y;
		}
	}
}