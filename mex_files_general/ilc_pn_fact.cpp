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

#include "fcns_cgpu_gen.h"
#include "pn_fact.hpp"

#include <mex.h>
#include "matlab_mex.h"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto psize = mex_get_pvctr<dt_float64>(prhs[0]);
	auto opt = (nrhs>1)?mex_get_num<mt::eDat_Sel_Typ>(prhs[1]):mt::edst_closest;

	/***************************************************************************************/
	auto psize_o = mex_create_pVctr<dt_float64>(psize.shape(), plhs[0]);

	mt::PN_Fact pn_fact;

	for (auto ik = 0; ik < psize.size(); ik++)
	{
		psize_o[ik] = pn_fact(psize[ik], opt);
	}
}