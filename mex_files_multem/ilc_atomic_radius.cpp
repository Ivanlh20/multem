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

#include "atomic_data_mt.cuh"
#include "atomic_fcns_mt.cuh"

#include "mex.h"
#include "matlab_mex.cuh"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	using T = dt_float64;

	auto atomic_pot_parm_typ = mex_get_enum<mt::eAtomic_Pot_Parm_Typ>(prhs[0]);
	auto dim = mex_get_num<dt_int32>(prhs[1]);
	auto vr_lim = mex_get_num<dt_float64>(prhs[2]);

	/***************************************************************************************/
	auto pradius = mex_create_pVctr<T>({mt::c_n_atom_typ, 3}, plhs[0]);

	mt::Atomic_Data atomic_data;
	mt::Atomic_Fcns<T> atomic_fcns;

	for(auto ik = 0; ik< pradius.s0(); ik++)
	{
		auto Z = ik+1;
		auto atomic_info = atomic_data(Z, atomic_pot_parm_typ);
		atomic_fcns.set_atomic_coef(atomic_info.coef[0]);

		pradius(ik, 0) = atomic_fcns.atomic_radius_rms(dim);
		pradius(ik, 1) = atomic_fcns.atomic_radius_cutoff(dim, vr_lim);
		pradius(ik, 2) = atomic_info.atomic_radius();
	}
}