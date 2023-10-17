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

#include "particles.cuh"
#include "space_group.hpp"

#include <mex.h>
#include "matlab_mex.h"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	using T = dt_float64;

	auto pasym_uc = mex_get_pvctr<T>(prhs[0]);
	auto sgn = mex_get_num<dt_int32>(prhs[1]);

	mt::Ptc_Atom<T> asym_uc(pasym_uc, {0, 0, 0}, false);
	mt::Ptc_Atom<T> base = mt::Space_Group<T>(asym_uc, sgn);

	auto pbase = mex_create_pVctr<T>({base.size(), base.cols_used}, plhs[0]);
	base.cpy_to_ptr(pbase.data(), base.size(), 0, base.cols_used);
}