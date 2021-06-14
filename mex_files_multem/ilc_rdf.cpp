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
	auto r_max = mex_get_num<T>(prhs[1]);
	auto nr = (nrhs>2)?mex_get_num<dt_int32>(prhs[2]):10;

	/***************************************************************************************/
	mt::Vctr_cpu<T> r(nr);
	mt::Vctr_cpu<T> rdf(nr);

	if (pvr.m_s1==2)
	{
		mt::Vctr_r_2d_cpu<T> vctr_r_2d(pvr.data(), pvr.s0());
		mt::fcn_rdf(vctr_r_2d, r_max, nr, r, rdf);
	}
	else
	{
		mt::Vctr_r_3d_cpu<T> vctr_r_3d(pvr.data(), pvr.s0());
		mt::fcn_rdf(vctr_r_3d, r_max, nr, r, rdf);
	}

	mex_create_set_pVctr<T>(plhs[0], r.ptr_64());
	mex_create_set_pVctr<T>(plhs[1], rdf.ptr_64());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT(mex_run, 0);
}