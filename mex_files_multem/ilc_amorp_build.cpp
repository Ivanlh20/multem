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

#include "types_mt.cuh"
#include "particles.cuh"
#include "amorp_build.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float32;

	const auto Z = mex_get_num<dt_int32>(prhs[0]);
	const auto rms_3d = mex_get_num<T>(prhs[1]);
	const auto occ = mex_get_num<T>(prhs[2]);
	const auto tag = mex_get_num<dt_int32>(prhs[3]);
	const auto bs = mex_get_r_3d<T>(prhs[4]);
	const auto d_min = mex_get_num<T>(prhs[5]);
	const auto rho = mex_get_num<T>(prhs[6]);
	const dt_int32 seed = (nrhs>7)?mex_get_num<dt_int32>(prhs[7]):300183;

	/***************************************************************************************/
	mt::R_3d<T> r_0(0, 0, 0);
	mt::Spec_Lay_Info<T> spec_lay_info(bs, r_0, tag);

	mt::Amorp_Build<T> amorp_build;
	auto atoms = amorp_build(Z, rms_3d, occ, d_min, rho, seed, spec_lay_info);

	auto patoms = mex_create_pVctr<dt_float64>({atoms.size(), atoms.cols_used}, plhs[0]);
	atoms.cpy_to_ptr(patoms.data(), atoms.size(), 0, atoms.cols_used);
}