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

#include "atomic_fcns_mt.cuh"

#include "mex.h"
#include "matlab_mex.h"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float64;

	auto atomic_pot_parm_typ = mex_get_num<mt::eAtomic_Pot_Parm_Typ>(prhs[0]);
	auto Z = mex_get_num<dt_int32>(prhs[1]);
	auto charge = mex_get_num<dt_int32>(prhs[2]);
	auto r = mex_get_pvctr<T>(prhs[3]);

	/***************************************************************************************/
	auto vzp = mex_create_pVctr<T>(r.shape(), plhs[0]);
	auto dvzp = mex_create_pVctr<T>(r.shape(), plhs[1]);

	mt::Atomic_Info_cpu<T> atomic_info = mt::Atomic_Data(Z, atomic_pot_parm_typ, charge);
	mt::Atomic_Fcns<T> atomic_fcns_mt(atomic_info.coef[0]);

	atomic_fcns_mt.vzp_dvzp(r.ptr_32(), vzp.ptr_32(), dvzp.ptr_32());
}