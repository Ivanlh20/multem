/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "cgpu_fcns_gen.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto pE_0 = mex_get_pvctr<dt_float64>(prhs[0]);
	auto pc_30 = mex_get_pvctr<dt_float64>(prhs[1]);

	/***************************************************************************************/
	auto pscherzer_defocus = mex_create_pVctr<dt_float64>(pE_0.shape(), plhs[0]);

	for (auto ik = 0; ik < pscherzer_defocus.size(); ik++)
	{
		pscherzer_defocus[ik] = mt::fcn_scherzer_defocus(pE_0[ik], pc_30[ik]*mt::c_mm_2_angs<dt_float64>);
	}
}