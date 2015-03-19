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

#include "hMT_General_CPU.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	int nx, ny, shift;
	double *fI, *radius;
 
	ny = (int)mxGetM(prhs[0]);
	nx = (int)mxGetN(prhs[0]);
	fI = mxGetPr(prhs[0]);
	shift = (int)mxGetScalar(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	radius = mxGetPr(plhs[0]);

	*radius = f_getFFT_InformationLimit_2D(ny, nx, shift, fI);
}