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

#include "const_enum.h"
#include "particles.cuh"
#include "xtl_build_in_parm.hpp"
#include "xtl_build.hpp"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void read_xtl_build_in_parm(const mxArray* mex_xtl_build_in_parm, mt::Xtl_Build_In_Parm<T>& xtl_build_in_parm)
{
	xtl_build_in_parm.a = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "a");
	xtl_build_in_parm.b = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "b");
	xtl_build_in_parm.c = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "c");

	xtl_build_in_parm.alpha = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "alpha")*mt::c_deg_2_rad<T>;
	xtl_build_in_parm.beta = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "beta")*mt::c_deg_2_rad<T>;
	xtl_build_in_parm.gamma = mex_get_num_from_field<T>(mex_xtl_build_in_parm, "gamma")*mt::c_deg_2_rad<T>;

	xtl_build_in_parm.n_a = mex_get_num_from_field<dt_int32>(mex_xtl_build_in_parm, "na");
	xtl_build_in_parm.n_b = mex_get_num_from_field<dt_int32>(mex_xtl_build_in_parm, "nb");
	xtl_build_in_parm.n_c = mex_get_num_from_field<dt_int32>(mex_xtl_build_in_parm, "nc");

	xtl_build_in_parm.sgn = mex_get_num_from_field<dt_int32>(mex_xtl_build_in_parm, "sgn");

	xtl_build_in_parm.pbc = mex_get_bool_from_field(mex_xtl_build_in_parm, "pbc");

	auto pasym_uc = mex_get_pvctr_from_field<T>(mex_xtl_build_in_parm, "asym_uc");
	xtl_build_in_parm.asym_uc.set_ptc(pasym_uc, {0, 0, 0}, false);

	auto pbase = mex_get_pvctr_from_field<T>(mex_xtl_build_in_parm, "base");
	xtl_build_in_parm.base.set_ptc(pbase, {0, 0, 0}, false);
 }

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float64;

	mt::Xtl_Build_In_Parm<T> xtl_build_in_parm;
	read_xtl_build_in_parm(prhs[0], xtl_build_in_parm);

	mt::Xtl_Build<T> xtl_build(xtl_build_in_parm);
	auto atoms = xtl_build();

	auto patoms = mex_create_pVctr<T>({atoms.size(), atoms.cols_used}, plhs[0]);
	atoms.cpy_to_ptr(patoms.data(), atoms.size(), 0, atoms.cols_used);
}