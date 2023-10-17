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
#include "types_mt.cuh"
#include "particles.cuh"

#include <mex.h>
#include "matlab_mex.h"
#include "matlab_multem_io.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	const auto patoms_i = mex_get_pvctr<T>(prhs[0]);
	auto theta = mex_get_num<T>(prhs[1])*mt::c_deg_2_rad<T>;
	auto u_0 = mex_get_r_3d<T>(prhs[2], mt::R_3d<T>(0, 0, 1));
	auto ctr_type = (nrhs>3)?mex_get_enum<mt::eRot_Ctr_Typ>(prhs[3]):mt::erct_geometric_ctr;
	auto ctr_p = (nrhs>4)?mex_get_r_3d<T>(prhs[4]):mt::R_3d<T>(0, 0, 0);

	/***************************************************************************************/
	mt::Ptc_Atom<T> atoms(patoms_i, {0, 0, 0}, false, true);

	u_0.normalize();

	if (mt::is_rot_ctr_geometric_ctr(ctr_type))
	{
		ctr_p = mt::R_3d<T>(atoms.r_mean.x, atoms.r_mean.y, atoms.r_mean.z);
	}

	atoms.rotate(theta, u_0, ctr_p);

	auto patoms_o = mex_create_pVctr<T>({atoms.size(), atoms.cols_used}, plhs[0]);
	atoms.cpy_to_ptr(patoms_o.data(), atoms.size(), 0, atoms.cols_used);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT_OUT(mex_run, 0);
}