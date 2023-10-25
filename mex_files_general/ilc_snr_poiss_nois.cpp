/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "types.cuh"
#include "matlab_types.cuh"
#include "type_traits_gen.h"
#include "host_device_functions.cuh"
#include "host_functions.hpp"

#include <mex.h>
#include "matlab_mex.h"

using mt::rmatrix_r;
using mt::e_host;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[ ]) 
{
	auto rM_i = mx_get_matrix<rmatrix_r>(prhs[0]);
	double SNR_i = (nrhs>1)?mx_get_scalar<double>(prhs[1]):0;

	/*******************************************************************/
	auto rM_o = mx_create_matrix<rmatrix_r>(rM_i.rows, rM_i.cols, plhs[0]);
	
	vector<float> M(rM_i.begin(), rM_i.end());
	float scf_o = 0;

	mt::Stream<e_host> stream(4);

	if (mt::isZero(SNR_i))
	{
		M = mt::fcn_add_poiss_nois(stream, M, 1.0);
	}
	else
	{
		M = mt::add_poiss_nois_by_SNR(stream, M, SNR_i, scf_o);
		auto rk = mx_create_scalar<rmatrix_r>(plhs[1]);
		rk[0] = scf_o;
	}

	rM_o.assign(M.begin(), M.end());
}