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

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::Vctr;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rx_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto ry_i = mex_get_pvctr<pMLD>(prhs[1]);
	auto y_thr = mex_get_num<dt_float64>(prhs[2]);
	auto typ = (nrhs>3)?mex_get_num<dt_int32>(prhs[3]):0;

	std::vector<dt_float64> x(rx_i.begin(), rx_i.end());
	std::vector<dt_float64> y(ry_i.begin(), ry_i.end());
	std::vector<dt_float64> peaks;
	
	if (typ == 0)
	{
		peaks = mt::fd_peaks_vector_typ_1(x, y, y_thr);
	}
	else
	{
		peaks = mt::fd_peaks_vector_typ_2(x, y, y_thr);
	}

	/***************************************************************************************/
	auto rpeaks = mex_create_pVctr<pMLD>(peaks.size(), 1, plhs[0]);
	std::copy(peaks.begin(), peaks.end(), rpeaks.real);
}