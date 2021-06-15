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

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	const auto patoms_i = mex_get_pvctr<T>(prhs[0]);
	const auto Z = mex_get_num<dt_int32>(prhs[1]);
	const auto rms_3d = mex_get_num<T>(prhs[2]);
	const auto occ = mex_get_num<T>(prhs[3]);
	const auto region = mex_get_num<dt_int32>(prhs[4]);
	const auto bs = mex_get_r_3d<T>(prhs[5]);
	const auto d_min = mex_get_num<T>(prhs[6]);
	const auto rho = mex_get_num<T>(prhs[7]);
	const auto spec_lay_typ = (nrhs>8)?mex_get_enum<mt::eSpec_Lay_Typ>(prhs[8]):mt::eSLT_Top;
	const dt_int32 seed = (nrhs>9)?mex_get_num<dt_int32>(prhs[9]):300183;

	/***************************************************************************************/
	mt::Ptc_Atom<T> atoms(patoms_i, {0, 0, 0}, false, true);

	mt::R_3d<T> r_0 = (mt::is_spec_lay_top(spec_lay_typ))?mt::R_3d<T>(0, 0, atoms.z_lim.x-bs.z):mt::R_3d<T>(0, 0, atoms.z_lim.y);
	mt::Spec_Lay_Info<T> spec_lay_info(bs, r_0, region, spec_lay_typ);

	mt::Amorp_Build<T> amorp_build;
	amorp_build(atoms, Z, rms_3d, occ, d_min, rho, seed, spec_lay_info);

	auto patoms_o = mex_create_pVctr<T>({atoms.size(), atoms.cols_used}, plhs[0]);
	atoms.cpy_to_ptr(patoms_o.data(), atoms.size(), 0, atoms.cols_used);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT_OUT(mex_run, 0);
}