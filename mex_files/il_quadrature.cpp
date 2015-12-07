/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#include <algorithm>
#include "types.cuh"
#include "quadrature.hpp"

#include <mex.h>
#include "matlab_mex.cuh"
	
using multem::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto q_type = mx_get_scalar<int>(prhs[0]);
	auto nq = mx_get_scalar<int>(prhs[1]);

	/*****************************************************************************/
	multem::Q1<double, multem::e_host> q;
	multem::Quadrature quadrature;
	quadrature.get(q_type, nq, q);

	rmatrix_r x = mx_create_matrix<rmatrix_r>(nq, 1, plhs[0]);
	rmatrix_r w = mx_create_matrix<rmatrix_r>(nq, 1, plhs[1]);

	std::copy(q.x.begin(), q.x.end(), x.begin());
	std::copy(q.w.begin(), q.w.end(), w.begin());
}