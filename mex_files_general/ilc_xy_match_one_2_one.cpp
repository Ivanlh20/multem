/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with Multem. If not, see <http://www.gnu.org/licenses/>.
 */

#define MATLAB_BLAS_LAPACK

#include "math_mt.h"
#include "vctr_cpu.h"
#include "particle_fcns.hpp"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
    auto pvr_2d_r = mex_get_pvctr<T>(prhs[0]);
    auto pvr_2d = mex_get_pvctr<T>(prhs[1]);

    mt::Vctr_r_2d_cpu<T> vr_2d_r(pvr_2d_r.data(), pvr_2d_r.s0());
    mt::Vctr_r_2d_cpu<T> vr_2d(pvr_2d.data(), pvr_2d.s0());
    mt::Vctr_cpu<dt_uint64> ind(pvr_2d.s0());

    mt::fcn_xy_match_one_2_one(vr_2d_r, vr_2d, ind);

    mex_create_set_pVctr<T>(plhs[0], vr_2d.ptr_64());

    if (nlhs > 1) {
        mex_create_set_pVctr<dt_int64>(plhs[1], ind.ptr_64());
    }
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
    MEX_RUN_FCN_FLOAT(mex_run, 0);
}