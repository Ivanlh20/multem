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

#include "const_enum.cuh"
#include "particles.cuh"
#include "xtl_build.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

template <class TIn_Xtl_Build>
void read_xtl_build(const mxArray *mex_in_xtl_build, TIn_Xtl_Build &in_xtl_build)
{
	using T = mt::Value_type<TIn_Xtl_Build>;

	in_xtl_build.a = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "a");
	in_xtl_build.b = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "b");
	in_xtl_build.c = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "c");

	in_xtl_build.alpha = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "alpha")*mt::c_deg_2_rad<T>;
	in_xtl_build.beta = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "beta")*mt::c_deg_2_rad<T>;
	in_xtl_build.gamma = mex_get_num_from_field<dt_float64>(mex_in_xtl_build, "gamma")*mt::c_deg_2_rad<T>;

	in_xtl_build.n_a = mex_get_num_from_field<dt_int32>(mex_in_xtl_build, "na");
	in_xtl_build.n_b = mex_get_num_from_field<dt_int32>(mex_in_xtl_build, "nb");
	in_xtl_build.n_c = mex_get_num_from_field<dt_int32>(mex_in_xtl_build, "nc");

	in_xtl_build.sgn = mex_get_num_from_field<dt_int32>(mex_in_xtl_build, "sgn");

	in_xtl_build.pbc = mex_get_num_from_field<dt_bool>(mex_in_xtl_build, "pbc");

	auto pasym_uc = mex_get_pvctr_from_field<T>(mex_in_xtl_build, "asym_uc");
	in_xtl_build.asym_uc.set_ptc(pasym_uc, {0, 0, 0}, false);

	auto pbase = mex_get_pvctr_from_field<T>(mex_in_xtl_build, "base");
	in_xtl_build.base.set_ptc(pbase, {0, 0, 0}, false);
 }

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float64;
	mt::In_Xtl_Build<T> in_xtl_build;
	read_xtl_build(prhs[0], in_xtl_build);

	mt::Xtl_Build<T> xtl_build(in_xtl_build);
	auto atoms = xtl_build();

	auto patoms = mex_create_pVctr<T>({atoms.size(), atoms.cols_used}, plhs[0]);
	atoms.cpy_to_ptr(patoms.data(), atoms.size(), 0, atoms.cols_used);
}