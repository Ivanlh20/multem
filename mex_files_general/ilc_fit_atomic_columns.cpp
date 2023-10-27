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
#include "eval_fit_atomic_columns.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rIxy_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto rRx_i = mex_get_pvctr<pMLD>(prhs[1]);
	auto rRy_i = mex_get_pvctr<pMLD>(prhs[2]);
	auto rpos_i = mex_get_pvctr<pMLD>(prhs[3]);
	auto rcoef_i = mex_get_pvctr<pMLD>(prhs[4]);
	auto lambda_0 = (nrhs>5)?mex_get_num<dt_float64>(prhs[5]):1e-2;
	auto rcstr_i = (nrhs>6)?mex_get_pvctr<pMLD>(prhs[6]):pMLD();

	std::vector<dt_float64> Ixy(rIxy_i.begin(), rIxy_i.end());
	std::vector<dt_float64> rx(rRx_i.begin(), rRx_i.end());
	std::vector<dt_float64> ry(rRy_i.begin(), rRy_i.end());
	std::vector<dt_float64> pos(rpos_i.begin(), rpos_i.end());
	std::vector<dt_float64> coef(rcoef_i.begin(), rcoef_i.end());
	std::vector<dt_float64> cstr(rcstr_i.begin(), rcstr_i.end());

	/***************************************************************************************/
	auto pcoef = mex_create_pVctr<pMLD>(1, coef.size(), plhs[0]);

	/***************************************************************************************/
	mt::Fit_Atomic_Columns<dt_float64, edev_cpu> fit_atomic_columns;
	auto norm = fit_atomic_columns.fit(Ixy, Rx, Ry, pos, lambda_0, coef, cstr);

	thrust::copy(coef.begin(), coef.end(), pcoef.begin());

	if (nlhs>1)
	{
		auto rnorm = mex_create_num<pMLD>(plhs[1]);
		rnorm[0] = norm;
	}
}