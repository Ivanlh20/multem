/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hGeneral_CPU.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	int nx, ny, shift;
	double dx, dy, k;
	double *f;
 
	ny = (int)mxGetScalar(prhs[0]);
	nx = (int)mxGetScalar(prhs[1]);
	dy = mxGetScalar(prhs[2]);
	dx = mxGetScalar(prhs[3]);
	k = mxGetScalar(prhs[4]);
	if(k>1) k = 1.0/k;
	shift = (int)mxGetScalar(prhs[5]);

	plhs[0] = mxCreateDoubleMatrix(ny, nx, mxREAL);
	f = mxGetPr(plhs[0]);

	f_getHanning_Filter_2D(ny, nx, dy, dx, k, shift, f);
}