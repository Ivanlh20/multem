/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hfxg_CPU.h"
#include <mex.h>

void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int PotPar, Z;
	double *g, *fxg, *dfxg;
	int m, n, ng;

	PotPar = (int)mxGetScalar(prhs[0]);
	Z = (int)mxGetScalar(prhs[1]);
	g = mxGetPr(prhs[2]);
	m = (int)mxGetM(prhs[2]); 
	n = (int)mxGetN(prhs[2]); 
	ng = m*n;
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	fxg = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
	dfxg = mxGetPr(plhs[1]);

	sAtomTypesCPU AtomTypesCPU;
	f_SetAtomTypes(Z, PotPar, 0, stVrl, AtomTypesCPU);

	cfxg_CPU fxg_CPU;
	fxg_CPU.SetAtomT(PotPar, AtomTypesCPU);
	fxg_CPU.fxg(ng, g, fxg, dfxg);
}