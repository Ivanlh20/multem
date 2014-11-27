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
	int nR, nRl, typ;
	double *R, *fR, *Rl, *rl, *frl, *cfrl;
 
	nR = (int)mxGetM(prhs[0])*(int)mxGetN(prhs[0]);
	R = mxGetPr(prhs[0]);
	fR = mxGetPr(prhs[1]);

	nRl = (int)mxGetM(prhs[2])*(int)mxGetN(prhs[2]);
	Rl = mxGetPr(prhs[2]);
	typ = (int)mxGetScalar(prhs[3]);

	/**************************************************************/
	plhs[0] = mxCreateDoubleMatrix(nRl-1, 1, mxREAL);
	rl = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(nRl-1, 1, mxREAL);
	frl = mxGetPr(plhs[1]); 

	plhs[2] = mxCreateDoubleMatrix(nRl-1, 1, mxREAL);
	cfrl = mxGetPr(plhs[2]);

	f_get2DRadDist(nR, R, fR, nRl, Rl, rl, frl, cfrl, true, typ);
}