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

#include "types.cuh"
#include "matlab_types.cuh"
#include "type_traits_gen.h"
#include "peak_finding.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMx_r;
using mt::pMx_c;
using mt::edev_cpu;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) 
{
	auto rIm = mex_get_mx<pMx_r>(prhs[0]);
	auto lx = rIm.cols*mex_get_num<double>(prhs[1]);
	auto ly = rIm.rows*mex_get_num<double>(prhs[2]);
	auto rp1 = mex_get_mx<pMx_r>(prhs[3]);
	auto rp2 = mex_get_mx<pMx_r>(prhs[4]);
	auto npro = mex_get_num<int>(prhs[5]);

	/*******************************************************************/
	vector<float> Im(rIm.begin(), rIm.end());
	mt::Grid_2d<float> grid_2d(rIm.cols, rIm.rows, lx, ly);
	mt::R_2d<float> p1(rp1[0], rp1[1]);
	mt::R_2d<float> p2(rp2[0], rp2[1]);

	auto profile = mt::intrpl_profile(grid_2d, Im, p1, p2, npro);
	auto rprofile = mex_create_mx<pMx_r>(profile.size(), 1, plhs[0]);
	rprofile.assign(profile.begin(), profile.end());
}