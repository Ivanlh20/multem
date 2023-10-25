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
	auto E_0 = mex_get_num<dt_float64>(prhs[0]);
	auto pbs = mex_get_pvctr<dt_float64>(prhs[1]);
	auto ptheta = mex_get_pvctr<dt_float64>(prhs[2]);

	/***************************************************************************************/
	auto psize = mex_create_pVctr<dt_float64>({ptheta.size(), pbs.size()}, plhs[0]);

	mt::PN_Fact pn_fact;
	for (auto ik = 0; ik < ptheta.size(); ik++)
	{
		// from milliradians to g
		auto g_max = mt::fcn_rad_2_rangs(E_0, ptheta[ik]*mt::c_mrad_2_rad<dt_float64>);

		// g_max = nx/(2*lx) -> nx = 2*g_max*lx
		for(auto ix = 0; ix< pbs.size(); ix++)
		{
			psize(ik, ix) = pn_fact(dt_int64(::round(2*pbs[ix]*g_max)), mt::edst_egreater_than);
		}
	}
}