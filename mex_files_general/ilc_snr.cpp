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
#include "fcns_gpu.h"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::edev_cpu;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto rM_i = mex_get_pvctr<pMLD>(prhs[0]);
	auto rM_r = mex_get_pvctr<pMLD>(prhs[1]);
	dt_int32 type = (nrhs>2)?mex_get_num<dt_int32>(prhs[2]):0;
	/***************************************************************************************/
	auto rsnr_o = mex_create_pVctr<pMLD>(1, 1, plhs[0]);
	
	vector<dt_float32> mx_i(rM_i.begin(), rM_i.end());
	vector<dt_float32> M_r(rM_r.begin(), rM_r.end());

	if (type==0)
	{
		rsnr_o[0] = mt::fcn_snr_ma1(mx_i, M_r);
	}
	else if (type==1)
	{
		rsnr_o[0] = mt::fcn_snr_ma2(mx_i, M_r);
	}
	else
	{
		rsnr_o[0] = min(mt::fcn_snr_ma1(mx_i, M_r), mt::fcn_snr_ma2(mx_i, M_r));
	}

}