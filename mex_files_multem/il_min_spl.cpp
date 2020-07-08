/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include "matlab_types.cuh"
#include "cgpu_fcns.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto E0 = mx_get_scalar<double>(prhs[0]);
	auto lx = mx_get_scalar<double>(prhs[1]);
	auto ly = (nrhs > 3)?mx_get_scalar<double>(prhs[2]):lx;
 int idx_theta = (nrhs==3)?2:3;
	auto theta = mx_get_scalar<double>(prhs[idx_theta])*mt::c_mrad_2_rad;

 auto g_req = mt::rad_2_rAngs(E0, theta);
 auto nx = static_cast<int>(round(2*lx*g_req));
 auto ny = static_cast<int>(round(2*ly*g_req));

	mt::Prime_Num pm;
	nx = pm(nx, mt::eDST_Closest);
	ny = pm(ny, mt::eDST_Closest);

	/******************************************************************/
	auto rnx = mx_create_scalar<rmatrix_r>(plhs[0]);
 rnx[0] = nx;
 if(nlhs>1)
 {
	 auto rny = mx_create_scalar<rmatrix_r>(plhs[1]);
 rny[0] = ny;
 }
}