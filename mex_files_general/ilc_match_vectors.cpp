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

#include "types.cuh"
#include "fcns_cpu.h"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rv_1 = mex_get_pvctr<pMLD>(prhs[0]);
	auto rv_2 = mex_get_pvctr<pMLD>(prhs[1]);

	/***************************************************************************************/
	mt::Vctr_cpu<dt_float32> v1(rv_1.begin(), rv_1.end());
	mt::Vctr_cpu<dt_float32> v2(rv_2.begin(), rv_2.end());

	mt::fcn_uniq_vctr(v1);
	mt::fcn_match_vctr(v1.begin(), v1.end(), v2, 0.1);

	auto rv_o = mex_create_pVctr<pMLD>(v2.size(), 1, plhs[0]);

	rv_o.assign(v2.begin(), v2.end());
}