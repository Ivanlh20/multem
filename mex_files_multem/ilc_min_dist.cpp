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

#include "types.cuh"
#include "particle_fcns.hpp"

#include <mex.h>
#include "matlab_mex.cuh"


template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto pvr = mex_get_pvctr<T>(prhs[0]);

	/***************************************************************************************/
	std::tuple<T, dt_int64, dt_int64> dist_t;

	if (pvr.m_s1==2)
	{
		mt::Vctr_r_2d_cpu<T> vctr_r_2d(pvr.data(), pvr.s0());
		dist_t = mt::fcn_vctr_r_nd_dist(vctr_r_2d, std::less<T>());
	}
	else
	{
		mt::Vctr_r_3d_cpu<T> vctr_r_3d(pvr.data(), pvr.s0());
		dist_t = mt::fcn_vctr_r_nd_dist(vctr_r_3d, std::less<T>());
	}

	mex_create_set_num<T>(plhs[0], std::get<0>(dist_t));

	if (nlhs>1)
	{
		mex_create_set_num<T>(plhs[1], std::get<1>(dist_t));
	}

	if (nlhs>2)
	{
		mex_create_set_num<T>(plhs[2], std::get<2>(dist_t));
	}

}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT(mex_run, 0);
}