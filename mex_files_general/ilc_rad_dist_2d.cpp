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
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#define MATLAB_BLAS_LAPACK

#include "const_enum.h"
#include "vctr_cpu.h"
#include "fcns_cpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void mex_run(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto pmx_i = mex_get_pvctr<T>(prhs[0]);
	auto bs_i = mex_get_r_3d_bs<T>(prhs[1], pmx_i.shape());
	auto r_c_i = mex_get_r_3d_r_c<T>(prhs[2], bs_i/T(2));
	auto radius_i = mex_get_num<T>(prhs[3]);
	auto typ = (nrhs>4)?mex_get_num<dt_int32>(prhs[4]):1;

	/***************************************************************************************/
	mt::Grid_2d<T> grid(bs_i.x, bs_i.y, pmx_i.s1(), pmx_i.s0(), T(0), T(0), false, false);
	const auto r_c = mt::R_2d<T>(r_c_i.x, r_c_i.y);
	const auto n_r = mt::fcn_max(grid.rx_2_irx_cd(radius_i), grid.ry_2_iry_cd(radius_i));

	auto prl = mex_create_pVctr<T>({n_r}, plhs[0]);
	auto pfrl = mex_create_pVctr<T>({n_r}, plhs[1]);
	auto pcfrl = mex_create_pVctr<T>({n_r}, plhs[2]);

	const bool reg = true;

	if (reg)
	{
		mt::fcn_rad_dist_2d_by_div(grid, pmx_i, r_c, radius_i, prl, pfrl, pcfrl, typ);
	}
	else
	{
		mt::fcn_rad_dist_2d_by_srch(grid, pmx_i, r_c, radius_i, prl, pfrl, pcfrl, typ);
	}
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	MEX_RUN_FCN_FLOAT_OUT(mex_run, 0);
}