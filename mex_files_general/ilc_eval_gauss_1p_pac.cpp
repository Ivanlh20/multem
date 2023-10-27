/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "eval_fit_gaussians.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rRx_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto rRy_i = mex_get_pvctr<pMLD>(prhs[1]);
	auto rpos_i = mex_get_pvctr<pMLD>(prhs[2]);
	auto rcoef_i = mex_get_pvctr<pMLD>(prhs[3]);

	std::vector<dt_float64> rx(rRx_i.begin(), rRx_i.end());
	std::vector<dt_float64> ry(rRy_i.begin(), rRy_i.end());
	std::vector<dt_float64> pos(rpos_i.begin(), rpos_i.end());
	std::vector<dt_float64> coef(rcoef_i.begin(), rcoef_i.end());

	std::vector<dt_float64> Ixy_o(rRx_i.size());
	/***************************************************************************************/
	auto rIxy_o = mex_create_pVctr<pMLD>(rRx_i.rows, rRx_i.cols, plhs[0]);

	/***************************************************************************************/
	mt::Eval_Gauss_1P_PAC<dt_float64, edev_cpu> Eval_Gauss_1P_PAC;
	Eval_Gauss_1P_PAC(Rx, Ry, pos, coef, Ixy_o);

	thrust::copy(Ixy_o.begin(), Ixy_o.end(), rIxy_o.begin());
}