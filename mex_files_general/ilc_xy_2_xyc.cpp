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
    auto pvr_2d = mex_get_pvctr<T>(prhs[0]);
    auto radius = (nrhs > 1) ? mex_get_num<T>(prhs[1]) : T(0.25);

    mt::Vctr_r_2d_cpu<T> vr_2d(pvr_2d.data(), pvr_2d.s0());

    auto r_c = mt::fcn_xy_2_xyc(vr_2d, radius);

    mex_create_set_pVctr<T>(plhs[0], r_c.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
    MEX_RUN_FCN_FLOAT(mex_run, 0);
}